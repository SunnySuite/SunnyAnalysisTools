################################################################################
# Energy broadening
################################################################################

function nonstationary_gaussian(instrument::ChopperSpec; intrinsic_width=0)
    (; Ei, L1, L2, L3, Δtp, Δtc, Δtd) = instrument

    FWHM_to_sigma = 1/2sqrt(2log(2))

    # Calculate the resolution kernel at some representative points.
    Es = range(-Ei, Ei, 200)
    dEs = [FWHM_to_sigma*energy_resolution(Ei, E, L1, L2, L3, Δtp, Δtd, Δtc) for E in Es]
    dEs = map(dEs) do dE
        sqrt(dE^2 + intrinsic_width^2)
    end

    # First a third-order polynomial to resulting values to determine a FWHM function.
    @. poly3(x, p) = p[1]*x^3 + p[2]*x^2 + p[3]*x + p[4]
    fit = curve_fit(poly3, Es, dEs, 0.1*ones(4))
    sigma(E) = poly3(E, fit.param)

    return Sunny.NonstationaryBroadening((b, ω) -> exp(-(ω-b)^2/2sigma(b)^2) / √(2π*sigma(b)^2))
end
# σ = FWHM_to_sigma*SunnyAnalysisTools.energy_resolution(Ei, μ > Ei ? Ei : μ, L1, L2, L3, Δtp, Δtd, Δtc)
# σ_tot = sqrt(σ^2 + λ^2)


# This is a stub -- insert real model of triple-axis instrument later
function nonstationary_gaussian(instrument::TripleAxisSpec)
    (; name, params) = instrument

    sigma = if name == "SPINS"
        E -> 0.9*(0.0949 + 0.0182E) / √2
    elseif name == "HMI"
        if params["resolution"] == "high"
            E -> (0.0886 + 0.0015E) / √2
        elseif params["resolution"] == "low"
            E ->  max( (-0.0594 + 0.0981E) / √2, 0.1 / √2 )
        else
            error("Unknown resolution setting for HMI")
        end
    else
        error("Unknown triple-axis instrument")
    end

    return Sunny.NonstationaryBroadening((b, ω) -> exp(-(ω-b)^2/2sigma(b)^2) / √(2π*sigma(b)^2))
end


################################################################################
# Daniel's instrument model
################################################################################

# Takes meV, returns m/s
energy_to_velocity(E_meV) = sqrt(2 * E_meV * J_per_meV / mₙ)

# Takes m/s, returns m (recall h has units J*s = kg*m*s^-1)
velocity_to_wavelength(v_mps) = h / (v_mps * mₙ)

# Input units: meV, meV, inverse angstrom
function theta(Ei_meV, E_meV, Q)
    clip(val, lo, hi) = min(max(val, lo), hi) 

    Ef_meV = Ei_meV - E_meV
    vi = energy_to_velocity(Ei_meV)
    vf = energy_to_velocity(Ef_meV)
    Q = Q * angstrom_per_meter 

    cos_2theta = (vi^2 + vf^2 - (Q * ħ/mₙ)^2) / (2 * vi * vf)
    if !(-1 <= cos_2theta <= 1)
        @warn "Cos(2θ) not between -1 and 1! Resulted in clipping."
        cos_2theta = clip(cos_2theta, -1., 1.)
    end

    return acos(cos_2theta)/2
end

function energy_resolution(Ei_meV, E_meV, L1, L2, L3, Δtp, Δtc, Δtd)
    Ef_meV = Ei_meV - E_meV
    vi = energy_to_velocity(Ei_meV)
    vf = energy_to_velocity(Ef_meV)
    return mₙ * sqrt(
        ( (vi^3) + (vf^3)*L2/L3      )^2 * (Δtp/L1)^2 +
        ( (vi^3) + (vf^3)*(L1+L2)/L3 )^2 * (Δtc/L1)^2 +
        ( (vf^3) )^2 * (Δtd/L3)^2
    ) / J_per_meV
end


# meV (Es), inverse Angstrom (Q), meters (Ls). Return inverse angstrom
function dQx(Ei, E, Q, L1, L2, L3, Δtp, Δtd, Δtc, Δθ)
    Ef = Ei - E
    vi = energy_to_velocity(Ei)
    vf = energy_to_velocity(Ef)
    θ = theta(Ei, E, Q)

    return (mₙ/ħ) * sqrt(
        ( vi^2 + vf^2 * (L2/L3) * cos(2θ)      )^2 * (Δtp/L1)^2 +
        ( vi^2 + vf^2 * ((L1+L2)/L3) * cos(2θ) )^2 * (Δtc/L1)^2 +
        ( vf^2 * cos(2θ) )^2 * (Δtd/L3)^2 +
        ( vf * sin(2θ) )^2 * Δθ^2 
    ) / angstrom_per_meter
end

function dQy(Ei, E, Q, L1, L2, L3, Δtp, Δtd, Δtc, Δθ)
    Ef = Ei - E
    vf = energy_to_velocity(Ef)
    θ = theta(Ei, E, Q)

    return (mₙ/ħ) * sqrt(
        ( vf^2 * (L2/L3) * sin(2θ)        )^2 * (Δtp/L1)^2 + 
        ( vf^2 * ((L1+L2)/L3) * sin(2θ)   )^2 * (Δtc/L1)^2 +
        ( vf^2 * sin(2θ) )^2 * (Δtd/L3)^2 + 
        ( vf * cos(2θ) * Δθ )^2
    ) / angstrom_per_meter
end

function dQ(Ei, E, Q, L1, L2, L3, Δtp, Δtd, Δtc, Δθ)
    Ef = Ei - E
    vi = energy_to_velocity(Ei)
    vf = energy_to_velocity(Ef)
    θ = theta(Ei, E, Q)

    dQx_val = dQx(Ei, E, Q, L1, L2, L3, Δtp, Δtd, Δtc, Δθ)
    dQy_val = dQy(Ei, E, Q, L1, L2, L3, Δtp, Δtd, Δtc, Δθ)

    Qx = (mₙ/ħ) * (vi - vf * cos(2θ)) / angstrom_per_meter
    Qy = (mₙ/ħ) * (-vf * sin(2θ)) / angstrom_per_meter

    return (1/Q) * sqrt((Qx*dQx_val)^2 + (Qy*dQy_val)^2)
end

# Vertical (out-of-plane) Q-resolution estimate. This is a deliberately simple
# geometric approximation (not the full Violini et al. detector-geometry
# treatment -- see project memory phase2_violini_deferred.md), valid near the
# horizontal/equatorial detector bank. If `height` (sample height, meters) is
# given, use it to estimate the vertical angular divergence; otherwise fall
# back to reusing the instrument's horizontal divergence Δθ.
function dQz(Ei, E, height, L3, Δθ)
    Δφ_vert = height > 0 ? atan(height / (2L3)) : Δθ
    vf = energy_to_velocity(Ei - E)
    return (mₙ/ħ) * vf * Δφ_vert / angstrom_per_meter
end

# Returns whether (Ei, E, Q) corresponds to a kinematically valid scattering
# triangle (i.e., theta()'s cos(2θ) computation would not need to clip). Used
# to exclude degenerate representative points from q_resolution_covariance,
# rather than duplicating theta()'s own clipping/warning behavior.
function kinematically_valid(Ei_meV, E_meV, Q)
    Ef_meV = Ei_meV - E_meV
    Ef_meV <= 0 && return false
    vi = energy_to_velocity(Ei_meV)
    vf = energy_to_velocity(Ef_meV)
    Qm = Q * angstrom_per_meter
    cos_2theta = (vi^2 + vf^2 - (Qm * ħ/mₙ)^2) / (2 * vi * vf)
    return -1 <= cos_2theta <= 1
end


################################################################################
# Q-resolution covariance construction
################################################################################

"""
    q_local_covariance(instrument::ChopperSpec, sample::Sample, Q, E)

The local (beam-frame: x=ki, y=in-plane-⊥, z=vertical) Q-resolution covariance
matrix (units Å⁻², i.e. already variance not FWHM) at a single (Q,E) point,
from the diagonal `dQx`/`dQy`/`dQz` treatment. No off-diagonal (Qx-Qy tilt)
term -- deferred, see project memory phase2_violini_deferred.md.
"""
function q_local_covariance(instrument::ChopperSpec, sample::Sample, Q, E)
    (; Ei, L1, L2, L3, Δtp, Δtc, Δtd, Δθ) = instrument
    FWHM_to_sigma = 1/2√(2log(2))

    σx = FWHM_to_sigma * dQx(Ei, E, Q, L1, L2, L3, Δtp, Δtd, Δtc, Δθ)
    σy = FWHM_to_sigma * dQy(Ei, E, Q, L1, L2, L3, Δtp, Δtd, Δtc, Δθ)
    σz = FWHM_to_sigma * dQz(Ei, E, sample.height, L3, Δθ)

    return diagm([σx^2, σy^2, σz^2])
end

"""
    mosaic_covariance(Q_crystal, mosaicity_rad)

Sample-mosaicity contribution to the Q-resolution covariance, in the crystal's
Cartesian frame: a broadening transverse to `Q_crystal`, scaling with
`|Q_crystal|². Independent of sample orientation -- always computable from
`Q_crystal` alone, since mosaic (rotational) smearing is by definition
transverse to `Q`.
"""
function mosaic_covariance(Q_crystal, mosaicity_rad)
    Qnorm = norm(Q_crystal)
    Qnorm < 1e-12 && return zeros(3, 3)
    Q̂ = Q_crystal / Qnorm
    P = Matrix{Float64}(I, 3, 3) - Q̂ * Q̂'
    return (mosaicity_rad * Qnorm)^2 .* P
end

"""
    q_resolution_covariance(binning::UniformBinning, instrument::ChopperSpec, sample::Sample)

Builds the single, fixed 3×3 Q-resolution covariance matrix (crystal Cartesian
frame, units Å⁻²) used by `StationaryQConvolution`, by evaluating the local
covariance (plus mosaicity, plus -- if `sample.orientation` is given -- proper
orientation into the crystal frame) at a handful of representative points
spanning the requested binning region (bin-set corners + center), then
averaging.

If `sample.orientation === nothing`, each point's covariance is collapsed to
an isotropic matrix (no orientation information available). Otherwise, each
point's local covariance is rotated into the crystal frame via a goniometer
angle solved exactly from that point's own Q; points outside the sample's
`(u,v)` plane are handled via their in-plane projection (always well-defined,
no failure mode), with the dropped out-of-plane offset retained only as a
diagnostic (see the warnings below) -- this function models resolution
broadening only, not instrument acceptance/coverage (whether a given Q was
actually reachable by this rotation scan).
"""
function q_resolution_covariance(binning::UniformBinning, instrument::ChopperSpec, sample::Sample)
    (; crystal, directions, Us, Vs, Ws, Es) = binning
    (; Ei) = instrument

    reps(xs) = length(xs) > 1 ? [xs[1], xs[end], (xs[1]+xs[end])/2] : [xs[1]]
    Ucands, Vcands, Wcands, Ecands = reps(Us), reps(Vs), reps(Ws), reps(Es)

    U0 = isnothing(sample.orientation) ? nothing : orientation_matrix(crystal, sample.orientation)

    Σs = Matrix{Float64}[]
    traces = Float64[]
    offset_ratios = Float64[]

    for (U, V, W, E) in Iterators.product(Ucands, Vcands, Wcands, Ecands)
        Q_hkl = directions * [U, V, W]
        Q_crystal = crystal.recipvecs * Q_hkl
        Qmag = norm(Q_crystal)

        kinematically_valid(Ei, E, Qmag) || continue

        Σ_local = q_local_covariance(instrument, sample, Qmag, E)

        if isnothing(U0)
            Σ_c = diagm(fill(tr(Σ_local)/3, 3))
        else
            θ = theta(Ei, E, Qmag)
            vi = energy_to_velocity(Ei)
            vf = energy_to_velocity(Ei - E)
            Qx = (mₙ/ħ) * (vi - vf*cos(2θ)) / angstrom_per_meter
            Qy = (mₙ/ħ) * (-vf*sin(2θ)) / angstrom_per_meter

            φ = solve_phi(U0, Q_crystal, Qx, Qy)
            Σ_c = rotate_covariance(U0, φ, Σ_local)

            w = U0' * Q_crystal
            push!(offset_ratios, abs(w[3]) / sqrt(Σ_local[3,3]))
        end

        Σ_c = Σ_c + mosaic_covariance(Q_crystal, sample.mosaicity)

        push!(Σs, Σ_c)
        push!(traces, tr(Σ_c))
    end

    isempty(Σs) && error("No kinematically valid representative Q,E points found for this binning region -- cannot build a StationaryQConvolution kernel.")

    if maximum(traces) / minimum(traces) > 5
        @warn "Q-resolution varies substantially (>5x in trace) across this binning region; StationaryQConvolution's single fixed kernel may be a poor approximation here."
    end

    if !isempty(offset_ratios) && maximum(offset_ratios) > 3
        @warn "Some representative Q,E points lie several vertical-resolution-widths outside the sample's (u,v) scattering plane; this binning region may extend beyond what this rotation scan geometrically covers. Note: this tool models resolution broadening only, not instrument acceptance/coverage."
    end

    return sum(Σs) / length(Σs)
end

# Probably eliminate
function reciprocal_lattice(a, b, c, α, β, γ)
    α, β, γ = [α, β, γ] .* (π/180)

    # Volume of the direct unit cell
    V = a * b * c * sqrt(
        1 - cos(α)^2 - cos(β)^2 - cos(γ)^2
        + 2 * cos(α) * cos(β) * cos(γ)
    )

    # Reciprocal lattice parameters
    a_star = (b * c * sin(α)) / V
    b_star = (a * c * sin(β)) / V
    c_star = (a * b * sin(γ)) / V

    return 2π * a_star, 2π * b_star, 2π * c_star  # Include 2π factor
end
