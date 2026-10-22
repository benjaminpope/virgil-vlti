# `virgil_vlti.templates`

Calibration templates: the calibrator PCA's components, packaged as nuisance modes per frame for virgil's gains (`OIData.with_gains`, `OIData.with_closure_offsets`). `virgil-vlti-template` builds one from the PCA's residuals, and `virgil-vlti-calibrate --template` writes its modes beside each calibrated file.

::: virgil_vlti.templates
