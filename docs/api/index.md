# API Reference

The calibration needs only numpy, scipy and astropy. virgil never imports this package: the two meet in files (calibrated OIFITS, and `.npz` nuisance modes).

- [Calibration](calibrate.md): transfer function, interpolation, errors, diameter modes, flux ratio
- [Products](products.md): reading VLTI OIFITS products; labelling science and calibrator exposures
- [Writing](write.md): calibrated OIFITS files and `.npz` nuisance modes
- [Loading into virgil](load.md): calibrated files and their modes as a virgil `OIData`
- [gravity_viscal port](viscal.md): calibration exactly as ESO's `gravity_viscal`
- [GRAVITY](gravity.md): header facts, dual-field geometry, anisoplanatic loss
- [MATISSE](matisse.md): header facts, refusals and BCD tables (Stage 1 of the MATISSE plan)
- [Pair test](pairs.md): error structure from consecutive frames
- [Calibrator PCA](calpca.md): residuals of calibrators and their principal components
- [Null datasets](nulls.md): calibrators calibrated as science targets, for false-alarm tests
- [Templates](templates.md): calibration templates and their nuisance modes
- [Simulated products](synthetic.md): GRAVITY products with a known truth, for examples and tests
- [Command-line entry points](cli.md): the five `virgil-vlti-*` commands
