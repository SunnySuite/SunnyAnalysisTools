# SunnyAnalysisTools 

This will be a collection of tools that fill in the gap between the [Sunny.jl](https://github.com/SunnySuite/Sunny.jl)
calculator and the fitting of models to neutron scattering data. Currently
contains simples tools for calculating binned intensities, including convolution
across inverse spatial dimensions in momentum space.

[![Build Status](https://github.com/ddahlbom/SunnyHelpersORNL.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ddahlbom/SunnyHelpersORNL.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Release status

Version 0.1.0 is a compatibility baseline release intended to provide a stable,
retrievable package version for existing scripts.

The next planned development line will focus the public API around
`intensities_binned`, dimensional binning/intensity types, and direct use of
Sunny spin-wave objects.
