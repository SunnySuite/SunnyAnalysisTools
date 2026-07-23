# # Going from data to optimization

# Initialize the relevant projects.

using Sunny, GLMakie, LinearAlgebra
ENV["JULIA_CONDAPKG_BACKEND"] = "Current"
ENV["JULIA_CONDAPKG_EXE"] = "C:\\Users\\uud\\AppData\\Local\\miniforge3\\condabin\\mamba.bat"
ENV["JULIA_PYTHONCALL_EXE"] = "C:\\Users\\uud\\AppData\\Local\\miniforge3\\envs\\neutron\\python.exe"
using PythonCall

using SunnyAnalysisTools

include("examples/data/model.jl")

# We'll be using Ba₃Mn₂O₈ as an example, a system consisting of ABC-stacked
# triangular lattice bilayers. Begin by loading the crystal. 

crystal_full = Crystal("examples/data/Ba3Mn2O8_OCD_2008132.cif"; symprec=1e-3)
crystal = subcrystal(crystal_full, "Mn1")

# Next establish the binning that will be used for the experimental data.

directions = [
    1  0  1 
    1  0 -1 
    0  1  0 
]
Us = 0.2:0.025:0.5 |> collect
Vs = 0:0.1:5.0 |> collect
Ws = [-0.05, 0.05]
Es = [1.05, 1.1]

binning = UniformBinning(crystal, directions, Us, Vs, Ws, Es)

# Load the experimental data. (This will become a PythonCall to shiver, but
# for now we load an ASCII export from Shiver using the above binning scheme.)
# Attach the instrument during this process, the data for which is loaded from PyChop.

instrument = cncs(; Ei=2.5)
observation = read_shiver_ascii("examples/data/CNCS_data.dat", binning; instrument)

# Visualize the loaded data.

fig = Figure()
plot_binned_data!(fig[1,1], observation)
fig

# Now set up our Sunny model. 

params = (; J0=1.575, J1=0.12, J2=0.256, J3=0.142, J4=0.037, D=0.03)
model = make_swt_model(crystal, params)

# Perform a Sunny calculation using this model using the binning and instrument
# information above.

calc_spec = UniformSampling(observation; nperqbin=2, nperebin=3)
calc_binned = intensities(model, calc_spec)

# Visualize the result.

fig = Figure()
plot_binned_data!(fig[1,1], calc_binned)
fig

# Set up a corresponding Sunny calculation using binning and q-convolution.

calc_spec = StationaryQConvolution(observation; nperqbin=2, nperebin=3)
calc_conv = intensities(model, calc_spec)

# Visualize the result.

fig = Figure()
plot_binned_data!(fig[1,1], calc_conv)
fig

# Calculate the χ² associated with each calculation.

err_binned = chi_square(observation, calc_binned; scale=4.0, intensity_normalization=true)
err_conv = chi_square(observation, calc_conv; scale=2.0, intensity_normalization=true)
println("χ² for binned data: $err_binned")
println("χ² for convolved data: $err_conv")

# Compare the calculation results with the experimental data.

fig = Figure(; size=(700, 350))
plot_binned_data!(fig[1,1], observation; title="CNCS")
plot_binned_data!(fig[1,2], calc_binned; title="Binned")
plot_binned_data!(fig[1,3], calc_conv; title="Convolved and binned")
fig