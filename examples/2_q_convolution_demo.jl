# # Q-space convolution with a Sample (orientation + mosaicity) on CNCS
#
# This script builds a simple Sunny model, sets up a CNCS instrument via
# PyChop, defines a `Sample` (crystal orientation + mosaicity + height), and
# computes a Q-convolved, binned calculation with `StationaryQConvolution`.
#
# It then compares this against a "reference" calculation: Sunny's own native
# `intensities`, with energy broadening only (no Q-space convolution), along a
# much more densely sampled Q-path. The comparison shows what the Q-space
# convolution smooths out relative to the underlying, infinitely-sharp-in-Q
# theoretical dispersion.

ENV["JULIA_CONDAPKG_BACKEND"] = "Current"
ENV["JULIA_PYTHONCALL_EXE"] = "path/to/your/python"    # <- adjust for your environment

using Sunny, GLMakie, LinearAlgebra, StaticArrays
using PythonCall
using SunnyAnalysisTools

## 1. A simple Sunny model with a nontrivial dispersion.
#
# Single-site ferromagnet with weak exchange and a strong field. The ground
# state is found via `minimize_energy!` rather than a hand-guessed
# polarization direction, since the latter is not guaranteed to be the true
# energy minimum (and `SpinWaveTheory` will error on an unstable ground
# state).

crystal = Crystal(lattice_vectors(5.0, 5.0, 5.0, 90, 90, 90), [[0.0, 0.0, 0.0]])
sys = System(crystal, [1 => Moment(; s=1, g=2)], :dipole)
set_exchange!(sys, -0.3, Bond(1, 1, [1, 0, 0]))
set_field!(sys, [0, 0, 3.0])
randomize_spins!(sys)
minimize_energy!(sys)
swt = SpinWaveTheory(sys; measure=ssf_trace(sys))

# Check the dispersion range along the path we'll examine, to sanity check the
# incident energy and energy-axis binning chosen below.
disp = dispersion(swt, [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
println("Dispersion at (0,0,0) and (0.5,0,0): ", disp, " meV")

## 2. CNCS instrument, via PyChop.
#
# Ei is chosen comfortably above the dispersion's maximum (~7.2 meV, see the
# printed range above) so the whole feature is kinematically accessible
# (energy transfer cannot exceed Ei) -- otherwise part of the Q-range would
# show zero intensity for kinematic reasons unrelated to the resolution
# convolution this script is meant to illustrate.

instrument = cncs(; Ei=15.0)

## 3. Sample: crystal mounted with a* along the incident beam and b*
## in-plane, some mosaic spread, and a modest sample height (used to estimate
## the vertical/out-of-plane Q-resolution).

orientation = SampleOrientation(SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0))
sample = Sample(; orientation, mosaicity_deg=0.5, height=0.01)

## 4. Binning: a Q-path along H (K, L held to a thin slice around 0) and an
## energy axis spanning the dispersion.

directions = Matrix{Float64}(I, 3, 3)
Us = collect(0.0:0.02:1.0)
Vs = [-0.025, 0.025]
Ws = [-0.025, 0.025]
Es = collect(5.0:0.1:8.0)
binning = UniformBinning(crystal, directions, Us, Vs, Ws, Es)

## 5. Q-convolved, binned calculation.

ekernel = nonstationary_gaussian(instrument)
calc_spec = StationaryQConvolution(binning, instrument, sample, ekernel; nperqbin=3, nghosts=2, nperebin=2)
calc_convolved = intensities(swt, calc_spec)

println("Q-resolution covariance used (crystal frame, Å⁻²):")
println(SunnyAnalysisTools.q_resolution_covariance(binning, instrument, sample))

## 6. Reference: Sunny's own `intensities`, energy-broadening only, along a
## much more finely sampled Q-path (no Q-space convolution at all).

Us_fine = collect(0.0:0.002:1.0)
qpath = [Sunny.Vec3(u, 0.0, 0.0) for u in Us_fine]
Es_fine = collect(5.0:0.02:8.0)
res_fine = Sunny.intensities(swt, qpath; energies=Es_fine, kernel=ekernel)

## 7. Compare visually: the Q-convolved/binned calculation should look like a
## blurred-out version of the dense reference path.

fig = Figure(; size=(1000, 450))
plot_binned_data!(fig[1, 1], calc_convolved; title="StationaryQConvolution\n(CNCS + Sample)")
ax = Axis(fig[1, 2]; xlabel="H (r.l.u.)", ylabel="E (meV)", title="Dense Q-path\n(energy broadening only)")
heatmap!(ax, Us_fine, Es_fine, res_fine.data'; colormap=:viridis)
fig
