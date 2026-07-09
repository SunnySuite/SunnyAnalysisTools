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

## LatinHyperCube Sampling

Time-of-flight calculations now include a `LatinHyperCube` calculation spec for
bin integration with configurable Q and E sampling counts:

```julia
spec = LatinHyperCube(binning, ekernel; nqpoints=8, nepoints=2)
calc = calculate_intensities(model, spec)
```

Behavior:

- Q sampling uses a Latin hypercube design with `nqpoints` samples per q-bin.
- E sampling uses `nepoints` uniformly spaced points per energy bin,
	centered on the bin center and spanning the bin width.
- Broadening energies are sorted before calling Sunny, then mapped back to the
	requested bins for accumulation.

Notes:

- `nqpoints` is required.
- `nepoints` defaults to `1`.

## Release Automation

This repository is configured to use Julia's standard registration/release flow:

1. Bump `version` in `Project.toml` and merge to the default branch.
2. Comment `@JuliaRegistrator register` on the commit (or PR) to trigger Registrator.
3. After the General registry PR is merged, TagBot creates the git tag/release.

Configured workflows:

- `.github/workflows/TagBot.yml`
- `.github/workflows/CompatHelper.yml`

Required GitHub secrets:

- `GITHUB_TOKEN` (provided automatically by GitHub Actions)
- One SSH private key usable by automation workflows:
	- preferred: `TAGBOT_SSH_KEY` (TagBot) and `COMPATHELPER_PRIV_KEY` (CompatHelper)
	- fallback: `DOCUMENTER_KEY` (used by both workflows if dedicated keys are absent)

Notes:

- CompatHelper opens dependency-compatibility PRs on a daily schedule.
- TagBot listens for `JuliaTagBot` comments and can also be run manually via workflow dispatch.
