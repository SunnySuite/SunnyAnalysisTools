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

    # Shared fixture for a small, stable, Q-dependent LSWT model: single-site
    # ferromagnet with weak exchange and a strong field, ground state found via
    # minimize_energy! (a manually-guessed polarization direction is not
    # guaranteed to be the true energy minimum and can make excitations!
    # fail with "Not an energy-minimum").
    function make_test_swt(crystal; seed=1)
        sys = Sunny.System(crystal, [1 => Sunny.Moment(; s=1, g=2)], :dipole)
        Sunny.set_exchange!(sys, -0.1, Sunny.Bond(1, 1, [1, 0, 0]))
        Sunny.set_field!(sys, [0, 0, 5.0])
        Random.seed!(seed)
        Sunny.randomize_spins!(sys)
        Sunny.minimize_energy!(sys)
        return Sunny.SpinWaveTheory(sys; measure=Sunny.ssf_trace(sys))
    end

    @testset "TOF intensities: renamed to `intensities`" begin
        crystal = Sunny.Crystal(Sunny.lattice_vectors(5.0, 5.0, 5.0, 90, 90, 90), [[0.0, 0.0, 0.0]])
        swt = make_test_swt(crystal)

        directions = Matrix{Float64}(I, 3, 3)
        binning = UniformBinning(crystal, directions, [-0.1, 0.1], [-0.1, 0.1], [-0.1, 0.1], [9.9, 10.5])
        ekernel = Sunny.gaussian(; fwhm=0.5)

        sq_spec = StationaryQConvolution(binning, 0.05, ekernel; nperqbin=2, nghosts=1, nperebin=2)
        us_spec = UniformSampling(binning, ekernel; nperqbin=2, nperebin=2)
        lhc_spec = LatinHyperCube(binning, ekernel; nqpoints=4, nepoints=3, rng=MersenneTwister(1))

        res_sq = intensities(swt, sq_spec)
        res_us = intensities(swt, us_spec)
        res_lhc = intensities(swt, lhc_spec)

        # Golden values captured from this exact fixture prior to the
        # calculate_intensities -> Sunny.intensities rename.
        @test isapprox(res_sq.data, [0.05972051613090612;;;;]; atol=1e-10)
        @test isapprox(res_us.data, [0.02371996622230557;;;;]; atol=1e-10)
        @test isapprox(res_lhc.data, [0.023588094826783205;;;;]; atol=1e-10)

        # Explicit identity rotation must reproduce default (no-`R`) behavior.
        @test intensities(swt, sq_spec; R=Sunny.Mat3(I)).data ≈ res_sq.data
        @test intensities(swt, us_spec; R=Sunny.Mat3(I)).data ≈ res_us.data
    end

    @testset "domain_average" begin
        crystal = Sunny.Crystal(Sunny.lattice_vectors(5.0, 5.0, 5.0, 90, 90, 90), [[0.0, 0.0, 0.0]])
        swt = make_test_swt(crystal)

        directions = Matrix{Float64}(I, 3, 3)
        binning = UniformBinning(crystal, directions, [-0.1, 0.1], [-0.1, 0.1], [-0.1, 0.1], [9.9, 10.5])
        ekernel = Sunny.gaussian(; fwhm=0.5)
        us_spec = UniformSampling(binning, ekernel; nperqbin=2, nperebin=2)

        @testset "identity rotation matches plain intensities" begin
            res_plain = intensities(swt, us_spec)
            res_da = domain_average(swt, us_spec; rotations=[([0, 0, 1], 0.0)], weights=[1.0])
            @test res_da.data ≈ res_plain.data
        end

        @testset "guard clauses" begin
            @test_throws ErrorException domain_average(swt, us_spec; rotations=[], weights=[])
            @test_throws ErrorException domain_average(swt, us_spec; rotations=[([0, 0, 1], 0.0)], weights=[1.0, 2.0])
            @test_throws ErrorException domain_average(swt, us_spec; rotations=[([0, 0, 1], 0.0)], weights=[0.0])
        end

        @testset "3-fold symmetry invariance (deterministic specs)" begin
            # A triangular crystal's 6-fold axis is exactly global ẑ when
            # α=β=90°, so a z-field/z-polarized ferromagnet on it is exactly
            # invariant under rotations about [0,0,1] by multiples of 2π/3 --
            # a genuine physical check that domain averaging is implemented
            # correctly, not just self-consistent.
            cryst3 = Sunny.triangular_crystal()
            swt3 = make_test_swt(cryst3; seed=2)

            directions3 = Matrix{Float64}(I, 3, 3)
            binning3 = UniformBinning(cryst3, directions3, [-0.1, 0.1], [-0.1, 0.1], [-0.1, 0.1], [9.9, 10.5])
            us_spec3 = UniformSampling(binning3, ekernel; nperqbin=2, nperebin=2)
            sq_spec3 = StationaryQConvolution(binning3, 0.05, ekernel; nperqbin=2, nghosts=1, nperebin=2)

            rotations3 = [([0, 0, 1], n * (2π / 3)) for n in 0:2]
            weights3 = [1.0, 1.0, 1.0]

            res_us3 = intensities(swt3, us_spec3)
            da_us3 = domain_average(swt3, us_spec3; rotations=rotations3, weights=weights3)
            @test isapprox(da_us3.data, res_us3.data; atol=1e-8)

            res_sq3 = intensities(swt3, sq_spec3)
            da_sq3 = domain_average(swt3, sq_spec3; rotations=rotations3, weights=weights3)
            @test isapprox(da_sq3.data, res_sq3.data; atol=1e-8)
        end
    end

    @testset "Q-resolution covariance (Phase 2)" begin
        crystal = Sunny.Crystal(Sunny.lattice_vectors(5.0, 5.0, 5.0, 90, 90, 90), [[0.0, 0.0, 0.0]])
        instrument = ChopperSpec("test", 5.0, 30.0, 2.0, 3.0, 5e-5, 3e-5, 1e-5, 1.5*π/180)
        directions = Matrix{Float64}(I, 3, 3)
        # Centered near Q=(1,0,0), a kinematically realistic (non-degenerate) point --
        # Q=(0,0,0) with nonzero energy transfer is never kinematically achievable
        # (|Q| has a nonzero minimum for inelastic scattering), so a binning region
        # centered exactly at the origin would leave no valid representative points.
        binning = UniformBinning(crystal, directions, [0.9, 1.1], [-0.1, 0.1], [-0.1, 0.1], [0.5, 1.5])

        @testset "orientation_matrix" begin
            o = SampleOrientation(SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0))
            U0 = SunnyAnalysisTools.orientation_matrix(crystal, o)
            @test isapprox(U0'*U0, I(3); atol=1e-10)
            @test isapprox(det(U0), 1.0; atol=1e-10)
        end

        @testset "solve_phi and rotate_covariance round-trip" begin
            o = SampleOrientation(SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0))
            U0 = SunnyAnalysisTools.orientation_matrix(crystal, o)

            Q_crystal = crystal.recipvecs * (directions * [1.0, 0.05, 0.0])
            θ = SunnyAnalysisTools.theta(5.0, 1.0, norm(Q_crystal))
            vi = SunnyAnalysisTools.energy_to_velocity(5.0)
            vf = SunnyAnalysisTools.energy_to_velocity(5.0 - 1.0)
            Qx = (SunnyAnalysisTools.mₙ/SunnyAnalysisTools.ħ) * (vi - vf*cos(2θ)) / SunnyAnalysisTools.angstrom_per_meter
            Qy = (SunnyAnalysisTools.mₙ/SunnyAnalysisTools.ħ) * (-vf*sin(2θ)) / SunnyAnalysisTools.angstrom_per_meter

            φ = SunnyAnalysisTools.solve_phi(U0, Q_crystal, Qx, Qy)
            w = U0' * Q_crystal
            Rz = [cos(φ) -sin(φ) 0; sin(φ) cos(φ) 0; 0 0 1]
            @test isapprox((Rz*w)[1:2], [Qx, Qy]; atol=1e-10)

            # solve_phi must remain well-defined (no error, no NaN) even for a
            # deliberately out-of-plane Q -- it always solves via the in-plane
            # projection rather than gating on reachability.
            Q_offplane = crystal.recipvecs * (directions * [1.0, 0.05, 0.3])
            φ_off = SunnyAnalysisTools.solve_phi(U0, Q_offplane, Qx, Qy)
            @test isfinite(φ_off)

            Σ_local = diagm([1.0, 2.0, 3.0])
            Σ_rot = SunnyAnalysisTools.rotate_covariance(U0, φ, Σ_local)
            @test isapprox(Σ_rot, Σ_rot'; atol=1e-10)
            @test isapprox(tr(Σ_rot), tr(Σ_local); atol=1e-10)  # rotation preserves trace
        end

        @testset "mosaic_covariance" begin
            Q = [1.5, 0.3, -0.2]
            @test isapprox(SunnyAnalysisTools.mosaic_covariance(Q, 0.0), zeros(3,3); atol=1e-12)
            @test isapprox(SunnyAnalysisTools.mosaic_covariance([0.0,0.0,0.0], 0.5), zeros(3,3); atol=1e-12)
            M = SunnyAnalysisTools.mosaic_covariance(Q, 0.01)
            @test isapprox(Q'*M*Q, 0.0; atol=1e-15*norm(Q)^4)  # purely transverse to Q
        end

        @testset "q_resolution_covariance" begin
            sample_none = Sample()
            Σ_none = SunnyAnalysisTools.q_resolution_covariance(binning, instrument, sample_none)
            @test isapprox(Σ_none, (tr(Σ_none)/3)*I(3); atol=1e-12)  # isotropic when no orientation given

            o = SampleOrientation(SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0))
            sample_o = Sample(; orientation=o, height=0.02)
            Σ_o = SunnyAnalysisTools.q_resolution_covariance(binning, instrument, sample_o)
            @test isapprox(Σ_o, Σ_o'; atol=1e-10)
            @test all(eigvals(Σ_o) .> 0)  # positive-definite

            sample_mosaic = Sample(; mosaicity_deg=0.5)
            Σ_mosaic = SunnyAnalysisTools.q_resolution_covariance(binning, instrument, sample_mosaic)
            @test tr(Σ_mosaic) > tr(Σ_none)  # mosaicity strictly adds broadening

            # A genuinely out-of-plane orientation should warn, not error, and
            # still return a valid (symmetric, positive-definite) covariance --
            # this tool models resolution broadening, not acceptance/coverage.
            bad_o = SampleOrientation(SVector(0.0, 0.0, 1.0), SVector(0.0, 1.0, 0.0))
            sample_bad = Sample(; orientation=bad_o)
            Σ_bad = @test_logs (:warn, r"vertical-resolution-widths") SunnyAnalysisTools.q_resolution_covariance(binning, instrument, sample_bad)
            @test isapprox(Σ_bad, Σ_bad'; atol=1e-10)
            @test all(eigvals(Σ_bad) .> 0)

            # All-invalid guard clause: a binning region collapsed to the
            # origin is kinematically unreachable for any nonzero energy
            # transfer, so no representative point should survive.
            binning_origin = UniformBinning(crystal, directions, [-0.1, 0.1], [-0.1, 0.1], [-0.1, 0.1], [0.5, 1.5])
            @test_throws ErrorException SunnyAnalysisTools.q_resolution_covariance(binning_origin, instrument, sample_none)
        end

        @testset "StationaryQConvolution(binning, instrument, sample, ekernel)" begin
            o = SampleOrientation(SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0))
            sample = Sample(; orientation=o, mosaicity_deg=0.3, height=0.02)
            ekernel = Sunny.gaussian(; fwhm=0.5)
            spec = StationaryQConvolution(binning, instrument, sample, ekernel; nperqbin=2, nghosts=1, nperebin=2)
            @test spec isa StationaryQConvolution

            # qfwhm must round-trip back to the same covariance via the
            # existing σ = qfwhm/2√(2log2) relation used inside StationaryQConvolution's
            # own (unmodified) constructor -- guards against a mistake in the
            # conversion-factor direction.
            Σ = SunnyAnalysisTools.q_resolution_covariance(binning, instrument, sample)
            @test isapprox(spec.qfwhm ./ (2*sqrt(2*log(2))), Σ; atol=1e-12)
        end
    end
end
