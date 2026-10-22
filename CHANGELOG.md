# Changelog

## 0.1.0

First public release.

- GRAVITY calibration of V², T3AMP and T3PHI (`virgil-vlti-calibrate`): transfer functions from calibrator or in-field calibrator exposures, interpolated in time, with leave-one-out errors; a port of ESO's `gravity_viscal` (`--match-gravity-viscal`).
- Nuisance modes as `.npz` sidecars (calibrator diameters, transfer function, calibration templates), carrying a schema string (`virgil_vlti.modes/1`) and per-array units; files without them are still read.
- `to_oidata`: calibrated files and sidecars to a virgil `OIData` with gains and closure offsets (needs `virgil-astro>=0.4`).
- Closure-phase-only output for nights without calibrators.
- Calibrator PCA (`virgil-vlti-calpca`) and calibration templates (`virgil-vlti-template`).
- PIONIER pair test (`virgil-vlti-pairs`) and MATISSE header handling; neither is calibrated here.
- Packaging: the `virgil` and `dev` extras require `virgil-astro>=0.4`; releases are published to PyPI by a trusted-publishing workflow gated on the tests passing for the tagged commit.
