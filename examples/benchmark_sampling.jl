using Sunny
using SunnyAnalysisTools
using Statistics
using Printf
using Random

include(joinpath(@__DIR__, "data", "model.jl"))

function make_demo_problem()
    crystal_full = Crystal(joinpath(@__DIR__, "data", "Ba3Mn2O8_OCD_2008132.cif"); symprec=1e-3)
    crystal = subcrystal(crystal_full, "Mn1")

    directions = [
        1  0  1
        1  0 -1
        0  1  0
    ]

    # Keep this moderate so benchmarks complete quickly.
    Us = collect(0.2:0.05:0.5)
    Vs = collect(0.0:0.2:1.0)
    Ws = [-0.05, 0.05]
    Es = collect(0.8:0.1:1.2)

    binning = UniformBinning(crystal, directions, Us, Vs, Ws, Es)

    params = (; J0=1.575, J1=0.12, J2=0.256, J3=0.142, J4=0.037, D=0.03)
    model = make_swt_model(crystal, params)

    return (; model, binning)
end

function summarize_times(name, times)
    @printf("%-24s mean=%8.3f s   min=%8.3f s   max=%8.3f s\n", name, mean(times), minimum(times), maximum(times))
end

function benchmark_sampling_methods(; repeats=3, nperqbin=2, nperebin=1, latin_nqpoints=nothing, latin_nepoints=nperebin, ekernel_fwhm=0.1)
    latin_nqpoints = isnothing(latin_nqpoints) ? nperqbin^3 : latin_nqpoints

    @assert repeats >= 1 "repeats must be >= 1"
    @assert nperqbin >= 1 "nperqbin must be >= 1"
    @assert nperebin >= 1 "nperebin must be >= 1"
    @assert latin_nqpoints >= 1 "latin_nqpoints must be >= 1"
    @assert latin_nepoints >= 1 "latin_nepoints must be >= 1"

    (; model, binning) = make_demo_problem()

    # Use a stationary Gaussian broadening for reproducible benchmark setup.
    ekernel = gaussian(; fwhm=ekernel_fwhm)

    uniform_spec = UniformSampling(binning, ekernel; nperqbin, nperebin)
    latin_spec = LatinHyperCube(binning, ekernel; nqpoints=latin_nqpoints, nepoints=latin_nepoints, rng=MersenneTwister(1234))

    println("Benchmark configuration")
    println("  repeats:       ", repeats)
    println("  nperqbin:      ", nperqbin)
    println("  nperebin:      ", nperebin)
    println("  latin_nqpoints:", latin_nqpoints, " (", nperqbin^3, " matches uniform q count when scalar nperqbin)")
    println("  latin_nepoints:", latin_nepoints)
    println("  ekernel_fwhm:  ", ekernel_fwhm)
    println()

    # Warmup: compile both paths before timing.
    calculate_intensities(model, uniform_spec)
    calculate_intensities(model, latin_spec)

    uniform_times = Float64[]
    latin_times = Float64[]

    for _ in 1:repeats
        push!(uniform_times, @elapsed calculate_intensities(model, uniform_spec))
        push!(latin_times, @elapsed calculate_intensities(model, latin_spec))
    end

    summarize_times("UniformSampling", uniform_times)
    summarize_times("LatinHyperCube", latin_times)

    ratio_mean = mean(latin_times) / mean(uniform_times)
    ratio_min = minimum(latin_times) / minimum(uniform_times)

    println()
    @printf("Relative speed (Latin / Uniform): mean=%6.2fx   min=%6.2fx\n", ratio_mean, ratio_min)

    return (; uniform_times, latin_times, ratio_mean, ratio_min)
end

if abspath(PROGRAM_FILE) == @__FILE__
    benchmark_sampling_methods()
end
