# MDE export schema: Mantid → Julia hand-off for MDNorm-style rehistogramming

## Purpose

This document specifies a portable file format for exporting everything needed
to rehistogram (normalize + bin) direct-geometry inelastic neutron data
*outside* Mantid, in Julia, using an algorithm equivalent to Mantid's `MDNorm`
(specifically its `MDNormDirectSC` branch). The goal is to hand off from
Mantid exactly once per dataset, then do all subsequent rehistogramming
(different binning schemes, symmetrizations, etc.) natively in Julia/Sunny
without needing a live Mantid session.

This is scoped to **direct-geometry, single-crystal, inelastic** data (CNCS,
SEQUOIA, HYSPEC, ARCS) — the case relevant to SunnyAnalysisTools. It is *not*
a general MDNorm replacement: powder/diffraction normalization uses a
`FluxWorkspace` (incident-flux spectrum integral) that this schema omits
entirely, because `MDNormDirectSC`'s inelastic branch never consults it — the
per-segment normalization weight there is just
`solid_angle_efficiency × proton_charge × energy_segment_width`.

Background research (Mantid source citations, algorithm walkthrough) that
motivated this schema lives in this session's history; the key source files
referenced are `Framework/MDAlgorithms/src/MDNorm.cpp`,
`MDNormDirectSC.cpp`, and `MDNormSCDPreprocessIncoherent.py`.

## Why the schema is split the way it is

MDNorm's efficiency trick is that the **normalization** (denominator) is
computed once per (detector, run) pair — a trajectory line through Q-E space,
intersected against the output bin grid — never per event. The **data**
(numerator) is an ordinary weighted histogram of the actual recorded events.
Keeping these two data sources separate in the export preserves that
efficiency: the event table can be huge, but the geometry/normalization
tables are `Ndetectors × Nruns` in size, typically orders of magnitude
smaller.

## Conventions (must hold across both files)

- Angles in **radians**.
- Energies in **meV**.
- `Q` (or `H,K,L`) and all rotation matrices (`goniometer_matrix`,
  `UB_matrix`, `symmetry_operations`) expressed in a **single, consistent
  basis** — recommended: reciprocal lattice units (RLU), matching
  SunnyAnalysisTools's existing convention (`crystal.recipvecs * q`)
  elsewhere in the package.
- Every file carries top-level attributes: `schema_version`, `instrument`,
  `q_convention` (`"RLU"` or `"Q_sample_cartesian"`), and, if RLU, the basis
  vectors used (equivalent to Mantid's `QDimension0/1/2`).
- Format: plain HDF5, written by a small Mantid/Python export script
  (`h5py`), read in Julia via `HDF5.jl`. Deliberately **not** Mantid's native
  `SaveMD` NeXus layout — that schema embeds `ExperimentInfo` blocks
  (instrument/logs/sample) in a form that isn't meant for outside
  consumption. This schema is self-designed and fully documented instead.

## File A — instrument/detector geometry (reusable across experiments)

One row per detector/spectrum, as Mantid's `SpectrumInfo` sees it. Reusable
for every dataset on the same instrument as long as the vanadium calibration
is unchanged — generate once per calibration cycle, not per sample run.

```
/detectors
    detector_id              int32   [Ndet]   # original Mantid detector ID (provenance/debugging)
    two_theta                float64 [Ndet]   # scattering angle, radians
    azimuthal                float64 [Ndet]   # "phi", radians
    solid_angle_efficiency   float64 [Ndet]   # scalar, from vanadium SolidAngleWorkspace
```

## File B — per-dataset (one per sample/experiment)

```
/runs                                          # one row per goniometer setting
    run_index                 int32   [Nruns]  # matches events/run_index
    goniometer_matrix         float64 [Nruns, 3, 3]
    Ei                        float64 [Nruns]  # meV; may be constant across runs
    proton_charge             float64 [Nruns]  # monitor/live-time normalization

/sample
    UB_matrix                 float64 [3, 3]   # shared across all runs (Mantid itself only reads run 0's)
    symmetry_operations       float64 [Nsym, 3, 3]  # resolved rotation matrices, NOT space-group names/symbols

/detector_kinematics                           # valid trajectory range per detector (Mantid's MDNorm_low/high)
    mdnorm_low                float64 [Ndet]  (or [Nruns, Ndet] if Ei varies by run)
    mdnorm_high                float64 [Ndet]  (or [Nruns, Ndet])

/events                                        # the large table
    Q                          float32 [Nevents, 3]   # or H,K,L — pick one, must match q_convention attribute
    deltaE                     float32 [Nevents]       # meV
    weight                     float32 [Nevents]       # typically 1.0; carries any prior corrections
    detector_index             int32   [Nevents]  # OPTIONAL — index into File A /detectors; only needed for masking/debugging
    run_index                  int32   [Nevents]  # OPTIONAL — index into /runs; only needed for masking/debugging
```

Note `/events` deliberately does not require `detector_index`/`run_index` for
the core rehistogramming path: building the data histogram is an ordinary
weighted bin of `(Q, deltaE, weight)`. The normalization histogram is built
entirely from `/detectors` × `/runs` × `/detector_kinematics`, independent of
`Nevents`.

## Open questions / known gaps (flag before implementing)

- **`mdnorm_low`/`mdnorm_high` are exported rather than re-derived.** These
  are plausibly a deterministic function of `Ei` and detector geometry
  (kinematically accessible `kf` range), but Mantid's exact formula — and
  whether it also encodes dead-pixel/masking information — has not been
  verified against source. Exporting them sidesteps a subtle correctness bug;
  revisit only if re-deriving them in Julia becomes worthwhile.
  See [[phase2_violini_deferred]] for the related precedent of deferring an
  unverified formula rather than guessing at it.
- **Detector grouping**: instruments that group multiple physical pixels into
  one spectrum (tube-summed, etc.) are assumed already resolved at the
  `/detectors` row level, matching Mantid's own `SpectrumInfo` — i.e. one row
  per *spectrum*, not necessarily one row per raw hardware pixel.
  Confirm this holds for all four target instruments before relying on it.
- **Multi-Ei / event-mode chopper data**: if a run mixes multiple incident
  energies (multi-rep chopper operation), `Ei` may need to vary event-by-event
  rather than per-run — the schema above assumes one `Ei` per run, which may
  need revisiting for those instruments/configurations.

## Status

Design sketch only — not yet implemented. No export script or Julia reader
exists yet.
