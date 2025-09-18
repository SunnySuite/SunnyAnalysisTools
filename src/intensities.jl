################################################################################
# Time-of-flight intensities calculations 
################################################################################
function calculate_intensities(swtmodel::SWTModel, broadening_spec::StationaryQConvolution; 
    observation = nothing,
    kwargs...
)
    (; swt) = swtmodel
    (; qpoints, epoints, qidcs, eidcs, qkernel, ekernel, binning) = broadening_spec
    (; qcenters, Es, binvol, crystvol) = binning

    # Calculate intensities for all points in subsuming grid around bin.
    res = Sunny.intensities(swt, qpoints[:]; energies=epoints, kernel=ekernel, kwargs...)
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

    if !isnothing(observation)
        for i in eachindex(res)
            if isnan(observation.ints[i])
                res[i] = NaN
            end
        end
    end

    # spec and params
    ModelCalculation(res, binning, broadening_spec, swtmodel.params)
end


function calculate_intensities(swtmodel::SWTModel, broadening_spec::UniformSampling; 
    unit_intensity=false, 
    thresh=1e-12, 
    observation = nothing, 
    kwargs...
)
    (; swt) = swtmodel
    (; qpoints, epoints, qidcs, eidcs, ekernel, binning) = broadening_spec
    (; qcenters, Es, binvol) = binning

    # Calculate intensities for all points in subsuming grid around bin.
    # res = Sunny.intensities(swt, qpoints[:]; energies=epoints, kernel=ekernel, kwargs...)

    dispersion_and_intensities = Sunny.intensities_bands(swt, qpoints[:])
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
        for (ei, qi) in Iterators.product(eidcs[j], qidcs[i])
            res[j, i] += data[ei,qi]
        end
        res[j, i] /= length(eidcs[j]) * length(qidcs[i])
    end
    res .*= abs(binvol)

    if !isnothing(observation)
        for i in eachindex(res)
            if isnan(observation.ints[i])
                res[i] = NaN
            end
        end
    end

    # spec and params
    ModelCalculation(res, binning, broadening_spec, swtmodel.params)
end


################################################################################
# TAX Intensities Functions
################################################################################

# Take some vectors to define a coordinate system, and then bounding values
# along each axis of this coordinate system. Use these to define a
# multidimensional parallelpiped.
function corners_of_parallelepiped(directions, bounds; offset=nothing)
    N = size(directions, 1)
    @assert length(bounds) == N "Number of bounds must equal number of dimensions (direction vectors)."
    offset = @something offset zeros(N)
    points = []
    for idx in Iterators.product([1:2 for _ in 1:N]...) 
        boundsvec = [bounds[i][idx[i]] for i in 1:N] 
        q_corner = offset + directions * boundsvec
        push!(points, q_corner)
    end
    return points
end

# Given some (hopefully linearly independent) vectors defining a coordinate
# system and a set of bounds on each dimension (in terms of the coordinate
# system), define a grid of points that encompasses the parallelpiped defined by
# these vectors and bounds. This is very rudimentary.
function grid_points(q::SVector{N, Float64}, directions, bounds, counts) where N
    corners = corners_of_parallelepiped(directions, bounds)
    extremas = [extrema([corner[i] for corner in corners]) for i in 1:N]
    grid = Iterators.product([range(extrema[1], extrema[2], counts[i]) for (i, extrema) in enumerate(extremas)]...)
    return [q + directions*SVector{N, Float64}(point) for point in grid]
end

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

function calculate_intensities(swtmodel::SWTModel, taxspec::TripleAxisGrid{2}; kwargs...)
    (; path, Ks, nsigmas, counts) = taxspec
    (; HKLs, Es, projection) = path
    buf = zeros(length(HKLs), length(path.Es))
    intfunc(hkls) = Sunny.intensities_bands(swtmodel.swt, hkls; kwargs...)
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


gaussian_func(v, μ, K) = exp(-((v-μ)' * K * (v-μ))/2)
sample_q(μ, Σ, N) = rand(MvNormal(μ, (Σ + Σ')/2), N)
sample_q(rng, μ, Σ, N) = rand(rng, MvNormal(μ, (Σ + Σ')/2), N)

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

function calculate_intensities(swtmodel::SWTModel, taxspec::TripleAxisMC{1}; kwargs...)
    (; path, N, Ks) = taxspec 
    (; HKLs, Es, projection) = path

    buf = zero(path.Es)
    intfunc(hkls) = Sunny.intensities_bands(swtmodel.swt, hkls; kwargs...)
    for (n, (HKL, K, E)) in enumerate(zip(HKLs, Ks, Es))
        q = projection*HKL
        buf[n] = tax_convolved_intensity_mc(intfunc, SVector{4, Float64}(q..., E), K, N)
    end
    return buf
end

function calculate_intensities(swtmodel::SWTModel, taxspec::TripleAxisMC{2}; kwargs...)
    (; path, N, Ks) = taxspec 
    (; HKLs, Es, projection) = path

    buf = zeros(length(HKLs), length(path.Es))
    intfunc(hkls) = Sunny.intensities_bands(swtmodel.swt, hkls; kwargs...)
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