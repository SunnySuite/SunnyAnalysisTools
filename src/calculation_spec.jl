
################################################################################
# Q-broadening only
################################################################################

# An AbstractConvolution should provide full specifications for how to treat
# intensities calculations from Sunny, including procedures for both energy
# and momentum convolution.
abstract type AbstractCalculationSpec end

# Add IFFT plan
struct StationaryQConvolution <: AbstractCalculationSpec
    binning  :: UniformBinning

    # Q-broadening
    qfwhm     :: Union{Float64, Array{Float64, 2}}   # Kernel FWHM. Not needed after construction, but useful as a record. 
    qkernel   :: Array{ComplexF64, 3}                # Fourier trasnformed convolution kernel
    qpoints   :: Array{Sunny.Vec3, 3}                # Sampled points in momentum space
    qidcs                                            # Indices corresponding to interior of q bins. Same dimensions as binning.qcenters.

    # E-broadening
    ekernel   :: Sunny.AbstractBroadening            # Sunny broadening to be passed to `Sunny.intensities`
    epoints   :: Vector{Float64}                     # Sampled points in momentum space
    eidcs                                            # Indices corresponding to interior of energy bins. Same dimensions as binning.Es
end

struct ModelCalculation
    data    :: Array{Float64, 4}
    binning :: AbstractBinning
    spec    :: AbstractCalculationSpec
    params  :: Union{Nothing, Dict{Any, Any}}
end

function Base.show(io::IO, ::ModelCalculation)
    printstyled(io, "Analog Calculation\n"; bold=true, color=:underline)
end

"""
    UniformQBroadening(binning::UniformBinning, qfwhm, ekernel; nperbin, nghosts, nperebin=1)

Generates an intensities calculation specification that performs convolution
over both energy and inverse Anstroms dimensions. The convolution kernel is
Gaussian and separable, i.e., the energy and spatial directions are
uncorrelated.

- `qfwhm` specifies the spatial convolution kernel: either a single number for
  an isotropic (spherical) kernel, or a 3x3 matrix for a general kernel (its
  columns are treated as the FWHM extent along each of that kernel's own
  principal axes). Units are inverse Angstrom.
- `ekernel` is a Sunny energy broadening kernel, e.g. `gaussian(fwhm)`.
- `nperbin` determines how many samples are included per linear dimension of a
  bin. This may be either a single number, or three numbers (one for each linear
  dimension).
- `nghosts` determines a number of padding bins to be included about the given
  binning scheme. These should be included if the size of the q-convolution
  kernel is large relative to the bin size.
- `nperebin` determines the number of samples to be included along the energy
  dimension of each bin. By default just the bin center is used.

"""
function StationaryQConvolution(binning::UniformBinning, qfwhm, ekernel; nperqbin, nghosts, nperebin=1)
    (; crystal, Δs, qcenters, Es, directions) = binning

    # Convert the FWHM specification to a covariance matrix. `S` is a
    # "square-root" factor of the covariance (a standard deviation for the
    # scalar/isotropic case, or a matrix whose S*S' gives the covariance for
    # the general case).
    S = qfwhm ./ (2√(2log(2)))
    Σ = isa(S, Number) ? (S^2)*I(3) : S*S'

    # Generate nperbin uniformly spaced samples for each bin, including in
    # padding bins.
    (; qpoints, epoints) = sample_binning(binning; nperqbin, nghosts, nperebin)

    # Keep track of which sample points are in which q-bins. *Note that the
    # indices need to be FFT shifted.*
    bounds = [(-Δ/2, Δ/2) for Δ in Δs[1:3]]
    points_shifted = fftshift(qpoints)
    qidcs = map(q0 -> find_points_in_bin(q0, directions, bounds, points_shifted), qcenters)

    # Calculate the q convolution kernel. Normalize so the kernel sums to
    # exactly 1 on the sample grid -- `gaussian_md` returns values of the
    # continuum-normalized PDF, which integrate (not sum) to 1, so using them
    # unnormalized as discrete convolution weights would scale the result by
    # the (grid-dependent) sample volume.
    binning_center = sum(qcenters)/length(qcenters)
    qkernel = gaussian_md(map(p -> crystal.recipvecs*(p-binning_center), qpoints), [0., 0, 0], Σ)
    qkernel ./= sum(qkernel)
    qkernel = fft(qkernel, (1, 2, 3))

    # Binning indices for energy axis.
    ΔE = Δs[4]
    eidcs = map(Es) do E0
        findall(E -> abs(E0 - E) < ΔE/2, epoints)
    end

    return StationaryQConvolution(binning, qfwhm, qkernel, qpoints, qidcs, ekernel, epoints, eidcs)
end

function Base.show(io::IO, ::StationaryQConvolution)
    printstyled(io, "Calculation Specification: Uniform Q-broadening\n"; bold=true, color=:underline)
end

################################################################################
# Bin sampling without convolution
################################################################################

struct UniformSampling <: AbstractCalculationSpec
    binning  :: UniformBinning

    # Q-sampling
    qpoints   :: Array{Sunny.Vec3, 3}                # Sampled points in momentum space
    qidcs                                            # Indices corresponding to interior of q bins. Same dimensions as binning.qcenters.

    # E-sampling and broadening
    ekernel   :: Sunny.AbstractBroadening            # Sunny broadening to be passed to `Sunny.intensities`
    epoints   :: Vector{Float64}                     # Sampled points in momentum space
    eidcs                                            # Indices corresponding to interior of energy bins. Same dimensions as binning.Es.
end

function Base.show(io::IO, ::UniformSampling)
    printstyled(io, "Uniform Bin Sampling Specification\n"; bold=true, color=:underline)
end


struct LatinHyperCube <: AbstractCalculationSpec
    binning  :: UniformBinning
    nqpoints :: Int
    nepoints :: Int
    rng      :: AbstractRNG

    # E-sampling and broadening
    ekernel   :: Sunny.AbstractBroadening            # Sunny broadening to be passed to `Sunny.intensities`
end

function Base.show(io::IO, ::LatinHyperCube)
    printstyled(io, "Latin Hypercube Bin Sampling Specification\n"; bold=true, color=:underline)
end

function LatinHyperCube(binning::UniformBinning, ekernel; nqpoints, nepoints=1, rng=Random.default_rng())
    return LatinHyperCube(binning, nqpoints, nepoints, rng, ekernel)
end

"""
    UniformSampling(binning::UniformBinning, ekernel; nperbin, nperebin=1)

"""
function UniformSampling(binning::UniformBinning, ekernel; nperqbin, nperebin=1)
    (; Δs, qcenters, Es, directions) = binning

    # Generate nperbin uniformly spaced samples for each bin, including in
    # padding bins.
    (; qpoints, epoints) = sample_binning(binning; nperqbin, nghosts=0, nperebin)

    # Keep track of which sample points are in which q-bins. 
    bounds = [(-Δ/2, Δ/2) for Δ in Δs[1:3]]
    qidcs = map(q0 -> find_points_in_bin(q0, directions, bounds, qpoints), qcenters)

    # Binning indices for energy axis.
    ΔE = Δs[4]
    eidcs = map(Es) do E0
        findall(E -> abs(E0 - E) < ΔE/2, epoints)
    end

    return UniformSampling(binning, qpoints, qidcs, ekernel, epoints, eidcs)
end


################################################################################
# Triple-axis calculations
################################################################################
struct TripleAxisMC{N}
    path  :: Union{TripleAxisPath, TripleAxis2DContour}
    Ks    # :: Array{SMatrix{4, 4, Float64, 16}, N}
    N     :: Int64
end

# Right now this is a useless type which is only used for dispatch. In the
# future, it can be used to store any precalculated heuristics, e.g. directions
# and bounds for each point or, for example, a complete and minimal set of
# qs that need to be calculate to convolve all points in the path.
# 
struct TripleAxisGrid{N}
    path  :: Union{TripleAxisPath, TripleAxis2DContour}
    Ks    # :: Array{SMatrix{4, 4, Float64, 16}, N}
    nsigmas
    counts
end

# See PythonToolsExt for constructor
# TripleAxisMC(path, instrument::TAVISpec; N)