################################################################################
# Sample orientation and mosaicity
################################################################################

"""
    SampleOrientation(u, v)

Defines the crystal's mounting orientation relative to the fixed incident beam
direction of a direct-geometry (rotation-scan) instrument, using the
Busing-Levy convention: `u` is the reciprocal-lattice direction (HKL) that
points along the incident beam at the reference goniometer angle (φ=0), and
`v` is a second HKL direction lying in the (nominally horizontal) scattering
plane at φ=0, chosen so that `u × v` points "up." A single fixed vertical lab
rotation axis is assumed, matching standard direct-geometry sample-rotation
stages.
"""
struct SampleOrientation
    u :: SVector{3, Float64}
    v :: SVector{3, Float64}
end

"""
    Sample(; mosaicity_deg=0.0, height=0.0, orientation=nothing)

Physical properties of the sample used in a neutron scattering measurement:

- `mosaicity_deg`: FWHM rotational spread of the sample mosaic, in degrees.
  Contributes a Q-space broadening transverse to `Q`, scaling with `|Q|`.
- `height`: sample height in meters, used to estimate the vertical
  (out-of-plane) instrumental Q-resolution via a simple geometric estimate. If
  `0.0` (not supplied), the vertical resolution falls back to reusing the
  instrument's horizontal divergence `Δθ` (a cruder placeholder).
- `orientation`: an optional `SampleOrientation`. If `nothing` (default), the
  Q-resolution covariance falls back to an isotropic collapse rather than
  being correctly oriented in the crystal's Cartesian frame.
"""
struct Sample
    mosaicity   :: Float64
    height      :: Float64
    orientation :: Union{Nothing, SampleOrientation}
end
function Sample(; mosaicity_deg=0.0, height=0.0, orientation=nothing)
    Sample(mosaicity_deg*(π/180), height, orientation)
end


"""
    orientation_matrix(crystal, o::SampleOrientation)

Builds `U0`, the 3×3 matrix whose columns are the crystal-frame Cartesian
representations of the beam-frame basis vectors (x=ki, y=in-plane-⊥, z=up) at
the reference goniometer angle (φ=0), via the Busing-Levy construction. Since
`crystal.recipvecs` already supplies the B-matrix (HKL → crystal Cartesian
Å⁻¹), only the orientation matrix `U0` itself (not a full UB matrix) is needed
here.
"""
function orientation_matrix(crystal::Sunny.Crystal, o::SampleOrientation)
    t1 = crystal.recipvecs * o.u
    t2p = crystal.recipvecs * o.v
    t3 = cross(t1, t2p)
    t2 = cross(t3, t1)

    û1 = t1 / norm(t1)
    û2 = t2 / norm(t2)
    û3 = t3 / norm(t3)

    return hcat(û1, û2, û3)
end


"""
    solve_phi(U0, Q_crystal, Qx, Qy)

Solves for the goniometer rotation angle `φ` (about the fixed vertical lab
axis) that places `Q_crystal` at beam-frame coordinates `(Qx, Qy)`. Exact and
single-valued, since only the in-plane components of `U0'*Q_crystal` enter the
formula: points outside the plane spanned by the sample's `(u,v)` orientation
are implicitly handled via their in-plane projection. This is not an
approximation of convenience -- it exactly reflects that a point's
out-of-plane offset is a physically separate effect (instrumental
acceptance/coverage, not orientation) that this function does not, and should
not, attempt to resolve. See `q_resolution_covariance` for how the dropped
out-of-plane component is retained as a diagnostic.
"""
function solve_phi(U0, Q_crystal, Qx, Qy)
    w = U0' * Q_crystal
    return atan(Qy, Qx) - atan(w[2], w[1])
end


"""
    rotate_covariance(U0, φ, Σ_local)

Rotates a beam-frame covariance matrix `Σ_local` into the crystal's Cartesian
frame, given the orientation matrix `U0` and the goniometer angle `φ` at which
`Σ_local` was evaluated.
"""
function rotate_covariance(U0, φ, Σ_local)
    Rz = [cos(-φ) -sin(-φ) 0; sin(-φ) cos(-φ) 0; 0 0 1]
    R = U0 * Rz
    return R * Σ_local * R'
end
