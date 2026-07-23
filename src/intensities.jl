################################################################################
# Time-of-flight intensities calculations 
################################################################################
function apply_observation_nan_mask!(res, observation)
    if isnothing(observation)
        return res
    end

    @assert size(res) == size(observation.ints) "observation intensities must match calculated intensity dimensions"

    for i in eachindex(res)
        if isnan(observation.ints[i])
            res[i] = NaN
        end
    end

    return res
end

function Sunny.intensities(swt::Sunny.AbstractSpinWaveTheory, broadening_spec::StationaryQConvolution;
    params=nothing,
    observation = nothing,
    R::Sunny.Mat3=Sunny.Mat3(I),
    kwargs...
)
    (; qpoints, epoints, qidcs, eidcs, qkernel, ekernel, binning) = broadening_spec
    (; qcenters, Es, binvol, crystvol) = binning

    # Calculate intensities for all points in subsuming grid around bin.
    # `R` rotates only the physical Q-values sampled from the model (e.g. for
    # crystallographic domain averaging); the bin/kernel geometry (qkernel,
    # qidcs, eidcs) is fixed in the lab frame and does not depend on `R`.
    qs = map(q -> R*q, qpoints[:])
    res = Sunny.intensities(swt, qs; energies=epoints, kernel=ekernel, kwargs...)
    data = reshape(res.data, (length(epoints), size(qpoints)...))

    # Convolve along Q-axes only using an FFT. Unfortunately, energy is the fast
    # axis. 
    data_ft = fft(data, (2, 3, 4))
    for i in axes(data_ft, 1)
        data_ft[i,:,:,:] .*= qkernel
    end
    data_conv = real.(ifft(data_ft, (2, 3, 4))) ./ prod(size(qpoints))

    # Sum over samples that lie within each bin and normalize by number of
    # samples.
    res = zeros(length(Es), size(qcenters)...)
    for i in CartesianIndices(qcenters), j in eachindex(Es)
        for (ei, qi) in Iterators.product(eidcs[j], qidcs[i])
            res[j, i] += data_conv[ei,qi]
        end
        res[j, i] /= length(eidcs[j]) * length(qidcs[i])
    end
    res .*= binvol

    apply_observation_nan_mask!(res, observation)

    # spec and params
    ModelCalculation(res, binning, broadening_spec, params)
end


function Sunny.intensities(swt::Sunny.AbstractSpinWaveTheory, broadening_spec::UniformSampling;
    params=nothing,
    unit_intensity=false,
    thresh=1e-12,
    observation = nothing,
    R::Sunny.Mat3=Sunny.Mat3(I),
    kwargs...
)
    (; qpoints, epoints, qidcs, eidcs, ekernel, binning) = broadening_spec
    (; qcenters, Es, binvol) = binning

    qs = map(q -> R*q, qpoints[:])
    dispersion_and_intensities = Sunny.intensities_bands(swt, qs)
    if unit_intensity
        dispersion_and_intensities.data .= map(dispersion_and_intensities.data) do val
            val > thresh ? 1.0 : 0.0
        end
    end
    res = Sunny.broaden(dispersion_and_intensities; energies=epoints, kernel=ekernel, kwargs...)
    data = reshape(res.data, (length(epoints), size(qpoints)...))

    # Sum over samples that lie within each bin and normalize by number of
    # samples.
    res = zeros(length(Es), size(qcenters)...)
    for i in CartesianIndices(qcenters), j in eachindex(Es)
        res[j, i] = accumulate_bin_average(data, eidcs[j], qidcs[i])
    end
    res .*= abs(binvol)

    apply_observation_nan_mask!(res, observation)

    # spec and params
    ModelCalculation(res, binning, broadening_spec, params)
end

"""
    Sunny.domain_average(swt::Sunny.AbstractSpinWaveTheory, spec::AbstractCalculationSpec;
                          rotations, weights, params=nothing, observation=nothing, kwargs...)

Average a TOF/MDE intensities calculation over a discrete set of crystallographic
domains. Each element of `rotations` (an `(axis, angle)` pair or a 3×3 matrix, in
global Cartesian coordinates) is converted via `Sunny.rotation_in_rlu` and passed
as the `R` keyword to `Sunny.intensities(swt, spec; R, kwargs...)`; the results
are combined with the matching entry of `weights` and normalized by
`sum(weights)`. Works for any `AbstractCalculationSpec` that accepts an `R`
keyword (`StationaryQConvolution`, `UniformSampling`, `LatinHyperCube`).

`observation`, if given, is applied once to the final averaged result (not
per-domain).

Note: for `LatinHyperCube` specs, each domain draws independent stochastic
samples from `spec.rng` (unbiased, but not variance-reduced via common random
numbers across domains).
"""
function Sunny.domain_average(swt::Sunny.AbstractSpinWaveTheory, spec::AbstractCalculationSpec;
    rotations, weights,
    params=nothing,
    observation=nothing,
    kwargs...
)
    isempty(rotations) && error("Rotations must be nonempty list")
    length(rotations) == length(weights) || error("Rotations and weights must be same length")
    sum(weights) == 0 && error("Sum of weights must be nonzero")

    (; binning) = spec
    (; qcenters, Es, crystal) = binning
    Rs = Sunny.rotation_in_rlu.(Ref(crystal), rotations)

    data = zeros(length(Es), size(qcenters)...)
    for (R, w) in zip(Rs, weights)
        res = Sunny.intensities(swt, spec; R, kwargs...)
        data .+= w .* res.data
    end
    data ./= sum(weights)

    apply_observation_nan_mask!(data, observation)

    ModelCalculation(data, binning, spec, params)
end


accumulate_bin_average(data, einds, qinds) = sum(data[ei, qi] for (ei, qi) in Iterators.product(einds, qinds)) / (length(einds) * length(qinds))

uniform_bin_samples(Ecenter, ΔE, nepoints) = [Ecenter - ΔE/2 + (i - 0.5) * ΔE / nepoints for i in 1:nepoints]


"""
Note: each Q-bin's Latin hypercube samples are drawn fresh from `spec.rng` on
every call, so repeated calls (including successive domains inside
`Sunny.domain_average`) are stochastically independent rather than correlated.
This is unbiased (each domain's estimate is unbiased regardless of which draw
produced it) but forgoes any variance reduction from reusing common random
numbers across domains/calls.
"""
function Sunny.intensities(swt::Sunny.AbstractSpinWaveTheory, broadening_spec::LatinHyperCube;
    params=nothing,
    unit_intensity=false,
    thresh=1e-12,
    observation = nothing,
    R::Sunny.Mat3=Sunny.Mat3(I),
    kwargs...
)
    (; binning, nqpoints, nepoints, rng, ekernel) = broadening_spec
    (; qcenters, Es, binvol, directions, Δs) = binning

    bounds = [(-Δ/2, Δ/2) for Δ in Δs[1:3]]
    ΔE = Δs[4]

    # Sample Q and E locally within each bin using Latin hypercubes.
    # Q samples and their dispersions are shared across all energy bins at a
    # fixed q-bin to avoid repeating Sunny.intensities_bands work.
    res = zeros(length(Es), size(qcenters)...)
    for i in CartesianIndices(qcenters)
        qcenter = SVector{3, Float64}(qcenters[i]...)
        qsamples = latin_hypercube_points(qcenter, directions, bounds, nqpoints; rng)
        qsamples_rot = map(q -> R*q, qsamples)

        dispersion_and_intensities = Sunny.intensities_bands(swt, qsamples_rot)
        if unit_intensity
            dispersion_and_intensities.data .= map(dispersion_and_intensities.data) do val
                val > thresh ? 1.0 : 0.0
            end
        end

        # Build all per-energy-bin uniform samples for this q-bin,
        # then issue a single broaden call with globally sorted energies.
        nEs = length(Es)
        all_esamples_unsorted = Vector{Float64}(undef, nepoints * nEs)
        for j in eachindex(Es)
            offset = (j - 1) * nepoints
            esamples_j = uniform_bin_samples(Es[j], ΔE, nepoints)
            all_esamples_unsorted[offset+1:offset+nepoints] = esamples_j
        end

        perm = sortperm(all_esamples_unsorted)
        all_esamples = all_esamples_unsorted[perm]
        row_for_flat = invperm(perm)

        broadened = Sunny.broaden(dispersion_and_intensities; energies=all_esamples, kernel=ekernel, kwargs...)
        data = reshape(broadened.data, (length(all_esamples), length(qsamples)))

        for j in eachindex(Es)
            offset = (j - 1) * nepoints
            erows = row_for_flat[offset+1:offset+nepoints]
            res[j, i] = accumulate_bin_average(data, erows, eachindex(qsamples))
        end
    end
    res .*= abs(binvol)

    apply_observation_nan_mask!(res, observation)

    ModelCalculation(res, binning, broadening_spec, params)
end


################################################################################
# TAX Intensities Functions
################################################################################

# For a single HKLE point and resolution kernel, calculate the convolved
# intensity using a Sunny SpinWaveTheory (swt). Sums over a grid of intensities
# at neihboring HKLs about the given point, with the intensities weighted by the
# convolution kernel. No effort here is made to normalize (i.e., there is no
# differential element, ΔHΔKΔLΔE, included in the sum.)
function tax_convolved_intensity_grid(intfunc, qe0, K, directions, bounds, counts)
    qh, qk, ql, _ = qe0

    HKLs = grid_points(SVector{3, Float64}(qh, qk, ql), directions, bounds, counts)
    HKLs = reshape(HKLs, length(HKLs))  # Interpret as linear array so intensities_bands will accept it

    (; data, disp) = intfunc(HKLs)

    cumval = 0.0
    for (iq, q) in enumerate(HKLs), iband in axes(disp, 1)
        qe = SVector{4, Float64}(q..., disp[iband, iq])
        cumval += data[iband, iq] * gaussian_func(qe, qe0, K)
    end

    return cumval # Multiply by differential when considering absolute units
end


function principal_axes_of_gaussian(Σ)
    vals, vecs = eigen(Σ)
    σs = sqrt.(vals)
    [σ*axis for (σ, axis) in zip(σs, eachcol(vecs))]
end

# Get the principal axes of the distribution as well as their relative 
# magnitudes in multiples of the corresponding eigenvalues corresponding 
# the principal axes. This is useful for defining a bounding box for the 
# grid of sampled qs.
function directions_and_bounds(Σ; nsigmas=3)
    vals, directions = eigen(Σ)
    σs = sqrt.(vals)
    bounds = [nsigmas .* (-σ, σ) for σ in σs]
    return (; directions, bounds)
end

function calculate_intensities(swt::Sunny.SpinWaveTheory, taxspec::TripleAxisGrid{2}; kwargs...)
    (; path, Ks, nsigmas, counts) = taxspec
    (; HKLs, Es, projection) = path
    buf = zeros(length(HKLs), length(path.Es))
    intfunc(hkls) = Sunny.intensities_bands(swt, hkls; kwargs...)
    for (n, ((j, HKL), (k, E))) in enumerate(Iterators.product(enumerate(HKLs), enumerate(Es)))
        K = Ks[n]
        q = projection*HKL
        (; directions, bounds) = directions_and_bounds(inv(K); nsigmas)
        directions = directions[1:3, 1:3]
        bounds = bounds[1:3]
        buf[j, k] = tax_convolved_intensity_grid(intfunc, SVector{4, Float64}(q..., E), K, directions, bounds, counts)  
    end
    return buf
end

function calculate_intensities(intfunc::Function, taxspec::TripleAxisGrid{2}; kwargs...)
    (; path, Ks, nsigmas, counts) = taxspec
    (; HKLs, Es, projection) = path
    buf = zeros(length(HKLs), length(path.Es))
    for (n, ((j, HKL), (k, E))) in enumerate(Iterators.product(enumerate(HKLs), enumerate(Es)))
        K = Ks[n]
        q = projection*HKL
        (; directions, bounds) = directions_and_bounds(inv(K); nsigmas)
        directions = directions[1:3, 1:3]
        bounds = bounds[1:3]
        buf[j, k] = tax_convolved_intensity_grid(intfunc, SVector{4, Float64}(q..., E), K, directions, bounds, counts)  
    end
    return buf
end

function calculate_intensities(intfunc::Function, path, Ks, nsigmas, counts; kwargs...)
    # (; path, Ks, nsigmas, counts) = taxspec
    (; HKLs, Es, projection) = path
    buf = zeros(length(HKLs), length(path.Es))
    for (n, ((j, HKL), (k, E))) in enumerate(Iterators.product(enumerate(HKLs), enumerate(Es)))
        K = Ks[n]
        q = projection*HKL
        (; directions, bounds) = directions_and_bounds(inv(K); nsigmas)
        directions = directions[1:3, 1:3]
        bounds = bounds[1:3]
        buf[j, k] = tax_convolved_intensity_grid(intfunc, SVector{4, Float64}(q..., E), K, directions, bounds, counts)  
    end
    return buf
end

function tax_convolved_intensity_mc(intfunc, qe0, K, nsamps)
    Σ = inv(K)
    qes = sample_q(qe0, Σ, nsamps)
    hkls = [Sunny.Vec3(qe[1], qe[2], qe[3]) for qe in eachcol(qes)]
    hkls = reshape(hkls, length(hkls))

    # (; data, disp) = intensities_bands(swt, hkls)
    (; data, disp) = intfunc(hkls)

    cumval = 0.0
    for (iq, q) in enumerate(hkls), iband in axes(disp, 1)
        qe = SVector{4, Float64}(q..., disp[iband, iq])
        cumval += data[iband, iq] * gaussian_func(qe, qe0, K) # Shouldn't have to multiply by gaussian_func...
        # cumval += data[iband, iq] 
    end

    return cumval/nsamps 
end

function calculate_intensities(swt::Sunny.SpinWaveTheory, taxspec::TripleAxisMC{1}; kwargs...)
    (; path, N, Ks) = taxspec 
    (; HKLs, Es, projection) = path

    buf = zero(path.Es)
    intfunc(hkls) = Sunny.intensities_bands(swt, hkls; kwargs...)
    for (n, (HKL, K, E)) in enumerate(zip(HKLs, Ks, Es))
        q = projection*HKL
        buf[n] = tax_convolved_intensity_mc(intfunc, SVector{4, Float64}(q..., E), K, N)
    end
    return buf
end

function calculate_intensities(swt::Sunny.SpinWaveTheory, taxspec::TripleAxisMC{2}; kwargs...)
    (; path, N, Ks) = taxspec 
    (; HKLs, Es, projection) = path

    buf = zeros(length(HKLs), length(path.Es))
    intfunc(hkls) = Sunny.intensities_bands(swt, hkls; kwargs...)
    for (n, ((j, HKL), (k, E))) in enumerate(Iterators.product(enumerate(HKLs), enumerate(Es)))
        K = Ks[n]
        q = projection*HKL
        buf[j,k] = tax_convolved_intensity_mc(intfunc, SVector{4, Float64}(q..., E), K, N)
    end
    return buf
end

# function calculate_intensities(intfunc::Function, path, Ks, nsigmas, counts; kwargs...)
function calculate_intensities(intfunc::Function, path, Ks, N; kwargs...)
    # (; path, N, Ks) = taxspec 
    (; HKLs, Es, projection) = path

    buf = zeros(length(HKLs), length(path.Es))
    for (n, ((j, HKL), (k, E))) in enumerate(Iterators.product(enumerate(HKLs), enumerate(Es)))
        K = Ks[n]
        q = projection*HKL
        buf[j,k] = tax_convolved_intensity_mc(intfunc, SVector{4, Float64}(q..., E), K, N)
    end
    return buf
end