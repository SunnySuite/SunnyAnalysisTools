# # Going from data to optimization

# Initialize the relevant projects.
using Pkg

using Sunny, GLMakie, LinearAlgebra
using SunnyAnalysisTools

# include("examples/data/model.jl")
include(joinpath(@__DIR__, "data", "model.jl"))

# We'll be using Ba₃Mn₂O₈ as an example, a system consisting of ABC-stacked
# triangular lattice bilayers. Begin by loading the crystal. 

crystal_full = Crystal("examples/data/Ba3Mn2O8_OCD_2008132.cif"; symprec=1e-3)
crystal = subcrystal(crystal_full, "Mn1")

# Load the experimental data. (This will become a PythonCall to shiver, but
# for now we load an ASCII export from Shiver.) The binning (projection
# directions, bin centers and widths) is parsed directly from the file's own
# header -- no `UniformBinning` needs to be constructed by hand. Attach the
# instrument during this process; its parameters are computed locally from
# bundled PyChop instrument data (no PyChop/PythonCall install needed).

instrument = cncs(; Ei=2.7, variant="High Flux")
observation = read_shiver_ascii("examples/data/ascii_slice.dat"; instrument)

# Visualize the loaded data.

fig = Figure()
plot_binned_data!(fig[1,1], observation)
fig

# Now set up our Sunny model. 

units = Units(:meV, :angstrom)
fielddir = [1/sqrt(2), -1/sqrt(2), 0]
field = fielddir*2units.T
params = (; J0=1.575, J1=0.12, J2=0.256, J3=0.142, J4=0.037, D=0.03)
model = make_swt_model(crystal, params; field)

# Perform a Sunny calculation using this model using the binning and instrument
# information above.

calc_spec = UniformSampling(observation; nperqbin=5, nperebin=5)
calc_binned = calculate_intensities(model, calc_spec)

# Visualize the result.

fig = Figure()
plot_binned_data!(fig[1,1], calc_binned)
fig

# Set up a corresponding Sunny calculation using binning and q-convolution.
# `qfwhm` is a placeholder pending a real instrument-based Q-resolution
# calculation (e.g. from resolution.jl's dQx/dQy).

calc_spec = StationaryQConvolution(observation; qfwhm=0.002, nperqbin=8, nperebin=8)
@time calc_conv = calculate_intensities(model, calc_spec)

# Visualize the result.

fig = Figure()
plot_binned_data!(fig[1,1], calc_conv)
fig

# Calculate the χ² associated with each calculation, using Sunny's own
# fitting utilities. `scale=true` lets Sunny fit the best-fit overall
# intensity scale between calculation and data, rather than guessing one by
# hand; weights come from the experimental errors (zero-error bins excluded).

weights = ifelse.(iszero.(observation.errs), 0.0, 1 ./ observation.errs .^ 2)
err_binned = squared_error_fitted(observation.ints, calc_binned.data; scale=true, weights).err
err_conv = squared_error_fitted(observation.ints, calc_conv.data; scale=true, weights).err
println("χ² for binned data: $err_binned")
println("χ² for convolved data: $err_conv")

# Compare the calculation results with the experimental data.

fig = Figure(; size=(700, 350))
plotopts1 = (; colorrange=(0, 0.006))
plotopts2 = (; colorrange=(0, 0.006))
plotopts3 = (; colorrange=(0, 0.004))
plot_binned_data!(fig[1,1], observation; plotopts=plotopts1, title="CNCS")
plot_binned_data!(fig[1,2], calc_binned; plotopts=plotopts2, title="Binned")
plot_binned_data!(fig[1,3], calc_conv;   plotopts=plotopts3, title="Convolved and binned")
fig