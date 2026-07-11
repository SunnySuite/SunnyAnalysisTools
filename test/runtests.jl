using SunnyAnalysisTools
using Sunny
using StaticArrays
using Random
using LinearAlgebra
using Test

@testset "SunnyAnalysisTools.jl" begin
    @testset "LatinHyperCube sampling" begin
        qcenter = SVector{3, Float64}(0.0, 0.0, 0.0)
        directions = Matrix{Float64}(I, 3, 3)
        bounds = [(-1.0, 1.0), (-2.0, 2.0), (-3.0, 3.0)]
        rng = MersenneTwister(42)
        points = SunnyAnalysisTools.latin_hypercube_points(qcenter, directions, bounds, 5; rng)

        @test length(points) == 5

        strata = [Int[] for _ in 1:3]
        for p in points
            local_coords = inv(directions) * (p - qcenter)
            for d in 1:3
                lo, hi = bounds[d]
                t = (local_coords[d] - lo) / (hi - lo)
                push!(strata[d], floor(Int, t * 5) + 1)
            end
        end

        for d in 1:3
            @test sort(strata[d]) == collect(1:5)
        end
    end

    @testset "LatinHyperCube constructor" begin
        crystal = Sunny.Crystal(Sunny.lattice_vectors(5.0, 5.0, 5.0, 90, 90, 90), [[0.0, 0.0, 0.0]])
        directions = Matrix{Float64}(I, 3, 3)
        binning = UniformBinning(crystal, directions, [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0])
        instrument = ChopperSpec("test", 5.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0)
        ekernel = nonstationary_gaussian(instrument)
        spec = LatinHyperCube(binning, ekernel; nqpoints=4, nepoints=3, rng=MersenneTwister(1))

        @test spec isa LatinHyperCube
        @test spec.nqpoints == 4
        @test spec.nepoints == 3
        @test spec.rng isa AbstractRNG
    end
end
