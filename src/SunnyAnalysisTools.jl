module SunnyAnalysisTools

import Sunny
using LinearAlgebra, Random, LsqFit, FFTW, StaticArrays
import DelimitedFiles: readdlm
import Distributions: MvNormal

################################################################################
# Files and exports 
################################################################################
include("constants.jl")

include("util.jl")

include("instruments.jl")
export ChopperSpec

include("data_modeling.jl")
export gaussian_mixture_model

include("binning.jl")
export UniformBinning

include("paths.jl")
export TripleAxisPath, TripleAxis2DContour

include("calculation_spec.jl")
export Calculation, StationaryQConvolution, UniformSampling, TripleAxisMC, TripleAxisGrid

include("resolution.jl")
export nonstationary_gaussian

include("observations.jl")
export TimeOfFlightObservation, read_shiver_ascii

include("models.jl")
export SWTModel

include("intensities.jl")
export ModelCalculation, calculate_intensities

include("parsing.jl")
# export read_shiver_ascii

include("errors.jl")
export chi_square, weighted_residual, squared_difference, similarity_measure


################################################################################
# Only load plotting functions if GLMakie imported. See Sunny.jl as reference
################################################################################
function is_pkg_loaded(pkg::Symbol)
    return any(k -> Symbol(k.name) == pkg, keys(Base.loaded_modules))
end

extension_fns = [
    :GLMakie => [:draw_boundary!, :visualize_binning, :plot_binned_data!, :plot_binnned_data],
    :PythonCall => [:cncs, :spins, :hmi, :TAVISpec, :TripleAxisMC, :TripleAxisGrid]
]

for (_pkg, fns) in extension_fns
    for fn in fns
        @eval function $fn end
        @eval export $fn
    end
end

function __init__()
    # Notify user if extension function requires package import
    if isdefined(Base.Experimental, :register_error_hint)
        Base.Experimental.register_error_hint(MethodError) do io, exc, _argtypes, _kwargs
            fn = Symbol(exc.f)
            for (pkg, fns) in extension_fns
                if in(fn, fns) && !is_pkg_loaded(pkg)
                    pkgstr = (pkg == :Makie) ? "a variant of Makie" : "package $pkg"
                    printstyled(io, "\nImport $pkgstr to enable `$fn`.\n"; bold=true)
                end
            end
        end
    end
end

# Access to PlottingExt module for developer convenience
PlottingExt() = Base.get_extension(@__MODULE__, :PlottingExt)
PythonToolsExt() = Base.get_extension(@__MODULE__, :PythonToolsExt)

end