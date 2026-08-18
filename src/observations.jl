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


################################################################################
# Shiver ASCII parsing
################################################################################

const _SHIVER_AXIS_RE = r"^#\s*([^\s:]+):\s*(\d+)\s*bins?\s*from\s*([^\s,]+)\s*to\s*([^\s,]+)(?:,\s*step\s*([^\s,]+))?\s*$"
const _SHIVER_LATTICE_RE = r"^#\s*Lattice parameters:\s*a\s*=\s*(\S+)\s*b\s*=\s*(\S+)\s*c\s*=\s*(\S+)\s*alpha\s*=\s*(\S+)\s*beta\s*=\s*(\S+)\s*gamma\s*=\s*(\S+)\s*$"
const _SHIVER_ORIENTATION_RE = r"^#\s*Orientation\s*u:\[([^\]]+)\]\s*v:\[([^\]]+)\]\s*$"
const _SHIVER_NAME_RE = r"^#\s*Name:\s*(.+?)\s*$"
const _SHIVER_INSTRUMENT_RE = r"^#\s*Instrument:\s*(.+?)\s*$"
const _SHIVER_EI_RE = r"^#\s*Ei:\s*(\S+)\s*$"
const _SHIVER_BINNING_MARKER_RE = r"^#\s*Binning:\s*$"
const _SHIVER_UNITCELL_MARKER_RE = r"^#\s*UnitCell:\s*$"
const _SHIVER_COLUMNS_RE = r"^#\s*Intensity\s+Error\s+(.+?)\s*$"
const _SHIVER_SHAPE_RE = r"^#\s*shape:\s*([\dx]+)\s*$"

# Parses the header lines of a Shiver ASCII export (everything up to the first
# non-'#' line). Section markers ("Binning:"/"UnitCell:") are matched rather
# than relying on fixed line numbers, since the exact set/order of metadata
# lines is not guaranteed to be stable across Shiver versions.
function _parse_shiver_header(lines)
    name = instrument_name = Ei = lattice = orientation_u = orientation_v = nothing
    axis_order = String[]
    axis_info = Dict{String, NamedTuple}()
    data_column_labels = String[]
    shape = Int[]

    for line in lines
        if (m = match(_SHIVER_NAME_RE, line)) !== nothing
            name = m[1]
        elseif (m = match(_SHIVER_INSTRUMENT_RE, line)) !== nothing
            instrument_name = m[1]
        elseif (m = match(_SHIVER_EI_RE, line)) !== nothing
            Ei = parse(Float64, m[1])
        elseif match(_SHIVER_BINNING_MARKER_RE, line) !== nothing || match(_SHIVER_UNITCELL_MARKER_RE, line) !== nothing
            continue
        elseif (m = match(_SHIVER_AXIS_RE, line)) !== nothing
            label = m[1]
            nbins = parse(Int, m[2])
            lo = parse(Float64, m[3])
            hi = parse(Float64, m[4])
            step = isnothing(m[5]) ? nothing : parse(Float64, m[5])
            push!(axis_order, label)
            axis_info[label] = (; nbins, lo, hi, step)
        elseif (m = match(_SHIVER_LATTICE_RE, line)) !== nothing
            lattice = parse.(Float64, (m[1], m[2], m[3], m[4], m[5], m[6]))
        elseif (m = match(_SHIVER_ORIENTATION_RE, line)) !== nothing
            orientation_u = parse.(Float64, split(m[1], ","))
            orientation_v = parse.(Float64, split(m[2], ","))
        elseif (m = match(_SHIVER_COLUMNS_RE, line)) !== nothing
            data_column_labels = split(strip(m[1]))
        elseif (m = match(_SHIVER_SHAPE_RE, line)) !== nothing
            shape = parse.(Int, split(m[1], "x"))
        end
    end

    isempty(axis_order) && error("Could not find any \"Binning:\" axis entries in Shiver header")
    isnothing(lattice) && error("Could not find a \"Lattice parameters:\" line in Shiver header")
    isempty(data_column_labels) && error("Could not find the \"Intensity Error ...\" column header line")
    isempty(shape) && error("Could not find a \"shape:\" line in Shiver header")

    return (; name, instrument_name, Ei, axis_order, axis_info, lattice, orientation_u, orientation_v, data_column_labels, shape)
end

# Parses a bracketed projection-direction label like "[H,H,0]" or "[H,-H,0]"
# into the RLU direction vector it names, e.g. (1,1,0) or (1,-1,0). Each
# comma-separated slot gives that vector's H/K/L component directly: a bare
# letter (H/K/L) contributes coefficient ±1 (its sign), an optional numeric
# prefix scales it (e.g. "2H" -> 2.0), and a pure number is a literal
# component (almost always 0 in practice, for an axis held fixed).
function _parse_shiver_direction(label::AbstractString)
    m = match(r"^\[(.*)\]$", label)
    isnothing(m) && error("Expected a bracketed direction label like \"[H,H,0]\", got \"$label\"")
    tokens = split(m[1], ",")
    length(tokens) == 3 || error("Direction label \"$label\" must have exactly 3 comma-separated components")
    return SVector{3,Float64}(map(_parse_shiver_direction_component, tokens))
end

function _parse_shiver_direction_component(token::AbstractString)
    token = strip(token)
    m = match(r"^([+-]?\d*\.?\d*)([HKL])$", token)
    if !isnothing(m)
        numpart = m[1]
        isempty(numpart) && return 1.0
        numpart in ("+", "-") && return numpart == "-" ? -1.0 : 1.0
        return parse(Float64, numpart)
    end
    return parse(Float64, token)
end

# Given an axis's parsed (nbins, lo, hi, step), returns (centers, Δ). For a
# single-bin axis there is no step in the header (only one bin, so no spacing
# to report) -- its width is just hi - lo, and its one center is the midpoint.
function _shiver_axis_centers(info)
    (; nbins, lo, hi, step) = info
    if nbins == 1
        return ([(lo + hi) / 2], hi - lo)
    else
        centers = [lo + step * (i + 0.5) for i in 0:nbins-1]
        return (centers, step)
    end
end

"""
    read_shiver_ascii(file; instrument=nothing, filtersubzeros=false)

Parse a Shiver ASCII-histogram export and return a `TimeOfFlightObservation`.
The crystal, projection directions, and bin centers/widths are all extracted
directly from the file's header -- no `UniformBinning` needs to be
constructed by hand. Instrument information in the header is currently
ignored (stored unparsed in `params`); it may optionally be attached to the
resulting Observation via the `instrument` keyword.
"""
function read_shiver_ascii(file; instrument=nothing, filtersubzeros=false)
    lines = readlines(file)
    nheader = findfirst(l -> !startswith(l, "#"), lines)
    isnothing(nheader) && error("File appears to contain only header lines (nothing but '#'-prefixed lines)")
    nheader -= 1

    header = _parse_shiver_header(@view lines[1:nheader])
    (; axis_order, axis_info, lattice, data_column_labels, shape,
       orientation_u, orientation_v, name, instrument_name, Ei) = header

    eidx = findfirst(==("DeltaE"), axis_order)
    isnothing(eidx) && error("Could not find a \"DeltaE\" axis among Binning: entries $(axis_order)")
    energy_label = axis_order[eidx]
    spatial_labels = filter(!=(energy_label), axis_order)
    length(spatial_labels) == 3 || error("Expected exactly 3 spatial binning axes plus DeltaE, got $(axis_order)")

    directions = Sunny.Mat3(hcat(map(_parse_shiver_direction, spatial_labels)...))

    a, b, c, alpha, beta, gamma = lattice
    crystal = Sunny.Crystal(Sunny.lattice_vectors(a, b, c, alpha, beta, gamma), [[0.0, 0.0, 0.0]])

    Us, ΔU = _shiver_axis_centers(axis_info[spatial_labels[1]])
    Vs, ΔV = _shiver_axis_centers(axis_info[spatial_labels[2]])
    Ws, ΔW = _shiver_axis_centers(axis_info[spatial_labels[3]])
    Es, ΔE = _shiver_axis_centers(axis_info[energy_label])

    labels = [energy_label, spatial_labels...]
    binning = UniformBinning(crystal, directions, Us, Vs, Ws, Es, [ΔU, ΔV, ΔW, ΔE]; labels)

    # Sanity check the data table's own declared shape against the header's
    # per-axis bin counts, since the two are written independently.
    for (label, n) in zip(data_column_labels, shape)
        axis_info[label].nbins == n || error("shape entry for \"$label\" ($n) does not match its Binning: entry ($(axis_info[label].nbins))")
    end

    # Read the intensities and errors, reshaping (to deal with data that must be
    # interpreted as row major) and permuting the dimensions to move the energy
    # axis first, in the same spatial order as `labels`.
    raw = readdlm(file; skipstart=nheader)
    perm = [findfirst(==(l), data_column_labels) for l in labels]
    ints, errs = map((1, 2)) do col
        vals = Float64.(raw[:, col])
        vals = PermutedDimsArray(reshape(vals, reverse(shape)...), [4, 3, 2, 1])
        permutedims(vals, perm)
    end

    params = Dict{String, Any}()
    isnothing(name) || (params["name"] = name)
    isnothing(instrument_name) || (params["instrument_name"] = instrument_name)
    isnothing(Ei) || (params["Ei"] = Ei)
    isnothing(orientation_u) || (params["orientation_u"] = orientation_u)
    isnothing(orientation_v) || (params["orientation_v"] = orientation_v)

    return TimeOfFlightObservation(binning, ints, errs; instrument, filtersubzeros, params)
end




function StationaryQConvolution(obs::TimeOfFlightObservation; qfwhm, nperqbin, nperebin=1, nghosts=[1,1,1])
    (; binning, instrument) = obs
    ekernel = nonstationary_gaussian(instrument)
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