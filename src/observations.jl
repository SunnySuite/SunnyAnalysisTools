abstract type AbstractObservation end

################################################################################
# Time-of-flight observations
################################################################################
struct TimeOfFlightObservation <: AbstractObservation
    instrument  :: Union{Nothing, AbstractInstrumentSpec}   # Data about the instrument -- stub for now
    binning     :: AbstractBinning                          # Binning used in observation
    ints        :: Array{Float64, 4}                        # Observed intensities
    errs        :: Array{Float64, 4}                        # Errors
    mask        :: Array{Float64, 4}                        # 1.0s and NaNs -- useful for plotting
    mask_idcs   :: Vector{CartesianIndex{4}}                # Indices of non-NaN values in ints and errs
    model       :: Union{Nothing, AbstractDataModel}        # Model of data -- stub for now.
    background  :: Union{Nothing, Function}                 # Model of background, given as a function f(q, E) (not sure worthy of a type yet)
    
    params      :: Dict{String, Any}
end

function TimeOfFlightObservation(binning, ints, errs; instrument=nothing, model=nothing, background=nothing, filtersubzeros=false, params=nothing)

    # Check that intensities and errors are compatible with given binning.
    (; Es, qcenters) = binning
    composite_size = (length(Es), size(qcenters)...)
    @assert composite_size == size(ints) == size(errs) "Size of errors and/or intensities not compatible with given binning scheme"

    # Remove negative values (from background subtraction)
    if filtersubzeros
        for idx in eachindex(ints)
            if ints[idx] < 0
                ints[idx] = errs[idx] = NaN
            end
        end
    end

    # Determine the mask (assuming empty bins are NaNs).
    mask_idcs = findall(val -> !isnan(val), ints)
    mask = NaN * ones(composite_size)
    mask[mask_idcs] .= 1.0

    params = @something params Dict{String, Any}()

    return TimeOfFlightObservation(instrument, binning, ints, errs, mask, mask_idcs, model, background, params)
end

function Base.show(io::IO, ::TimeOfFlightObservation)
    printstyled(io, "Time-of-Flight Observation\n"; bold=true, color=:underline)
end


"""
    read_shiver_ascii(file, binning; instrument=nothing)

Parse a Shiver file and return an Observation. An appropriate UniformBinning
must be provided. Instrument information may optionally be attached to the
resulting Observation.
"""
function read_shiver_ascii(file, binning; instrument=nothing, filtersubzeros=false, permutation=nothing)
    data = readdlm(file)

    # Read metadata from the shiver file
    labels = data[1,4:7]
    shape = parse.(Int64, split(data[2,3], "x"))

    # Set up a permutation that moves the energy axis to the first index.
    eidx = findfirst(label -> label==("DeltaE"), labels)
    permutation = @something permutation [1, 2, 3, 4]
    permutation = [eidx, deleteat!(permutation, eidx)...]  # Shift energy to first position

    # Read the intensities and errors, reshaping (to deal with data that must be
    # interpreted as row major) and permuting the dimensions to move the energy
    # axis first.
    ints, errs = map([data[3:end,1], data[3:end,2]]) do vals
        vals = Float64.(vals)
        vals = PermutedDimsArray(reshape(vals, reverse(shape)...), [4, 3, 2, 1])
        vals = permutedims(vals, permutation)
        vals
    end

    binning.labels = labels[permutation]

    return TimeOfFlightObservation(binning, ints, errs; instrument, filtersubzeros)
end




function StationaryQConvolution(obs::TimeOfFlightObservation; nperqbin, nperebin=1, nghosts=[1,1,1])
    (; binning, instrument) = obs
    ekernel = nonstationary_gaussian(instrument)
    qfwhm = 0.002 # Potemkin village here. Use the chopper spec details to calculate with resolution.jl
    StationaryQConvolution(binning, qfwhm, ekernel; nperqbin, nperebin, nghosts)
end

function UniformSampling(obs::TimeOfFlightObservation; nperqbin, nperebin=1)
    (; binning, instrument) = obs
    ekernel = nonstationary_gaussian(instrument)
    UniformSampling(binning, ekernel; nperqbin, nperebin)
end

function LatinHyperCube(obs::TimeOfFlightObservation; nqpoints, nepoints=1, rng=Random.default_rng())
    (; binning, instrument) = obs
    ekernel = nonstationary_gaussian(instrument)
    LatinHyperCube(binning, ekernel; nqpoints, nepoints, rng)
end


################################################################################
# Triple-axis observations
################################################################################

struct TripleAxisObservation <: AbstractObservation
    instrument  :: Union{Nothing, AbstractInstrumentSpec}   # Data about the instrument -- stub for now
    Qs          :: Vector{Sunny.Vec3}
    Es          :: Vector{Float64}
    ints        :: Array{Float64, 4}                        # Observed intensities
    errs        :: Array{Float64, 4}                        # Errors
    model       :: Union{Nothing, AbstractDataModel}        # Model of data -- stub for now.
    background  :: Union{Nothing, Function}                 # Model of background, given as a function f(q, E) (not sure worthy of a type yet)
    
    params      :: Dict{String, Any}
end