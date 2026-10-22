# virgil-vlti: plan

Status: 2026-10-04. Created from the "separate GRAVITY project" that virgil's design notes deferred until the core worked (2026-10-03), and widened to PIONIER and MATISSE. Aperture masking stays out of scope: it may get its own repository.

## The boundary with virgil
**virgil (core)** keeps the generic machinery, and must never import this package:
- `read_oifits(insname=...)` and its refusal to merge GRAVITY FT and SC tables;
- the refusal of differential `VISPHI` as an absolute phase;
- whitened, independent closure phases (virgil#120);
- times and frames in `OIData` (Stage 6a.0);
- differential phases, smearing and a generic primary beam (6a);
- low-rank calibration nuisances that accept externally supplied modes and widths (6d).

**virgil-vlti** holds what depends on one instrument: calibration from pipeline products, empirical error models, templates and priors, and the instrument facts behind them (layouts, fibre fields, wavelength accuracy). It depends on numpy, scipy and astropy, and on virgil only where it hands data or nuisance modes to it.

## Done
- **Calibration from calibrator exposures** (`virgil_vlti.calibrate`, `virgil-vlti-calibrate`), generalised from an earlier target-specific GRAVITY script:
  - calibrators found by the GRAVITY dual-field fibre offset (in-field calibrator), by name, or by `ESO PRO CATG`;
  - optional uniform-disk calibrator diameters;
  - linear or Gaussian-weighted (the pipeline manual's formula) interpolation in time; closure phases as unit vectors;
  - the transfer-function uncertainty as a leave-one-out prediction error (review §7c.1) or the night's scatter;
  - split polarisations calibrated separately; FT tables never read;
  - one MJD per frame across baselines and triangles;
  - a science-to-calibrator `OI_FLUX` ratio;
  - calibrator diameter errors as nuisance modes, one per calibrator and night (review §7c.3), written beside each file for virgil's Stage 6d;
  - in dual field, the predicted anisoplanatic V² loss of each science exposure (GRAVITY Wide's model, from the header's seeing and θ₀), reported in `index.csv` and warned about above 2% (review §7c.2);
  - the pair test (`virgil_vlti.pairs`, `virgil-vlti-pairs`): white, per-baseline and per-telescope V² errors and closure-phase offsets, fitted from consecutive frames; rerun on IW Car (design/pionier_systematics.md);
  - **empirical null datasets** (`virgil_vlti.nulls`): each calibrator exposure calibrated as a science target by the night's other calibrators and written as OIFITS (uniform-disk truth of the catalogue diameter, recorded in the manifest), with possible binaries flagged and kept out of the references only, and block resampling of the residuals; the false-alarm reference of [`design/null_false_alarm.md`](design/null_false_alarm.md);
  - calibrator diameters from GRAVITY `*_CAL_VIS` products, and optional spectral smoothing of the transfer function (`--smooth-tf 5`) for parity with `gravity_viscal`, against whose calibrated products the calibration has been validated (V² level 0.12%, closure phases 0.06°);
  - **calibration templates** (`virgil_vlti.templates`, `virgil-vlti-template`, `--template`): the calibrator PCA's components as per-telescope and per-baseline log|V| gain modes and optional closure-phase offsets, for virgil's Stage 6d; the pilot template (GRAVITY single field, MEDIUM, combined, 2023–25) ships with the package;
  - T3AMP, and optionally `OI_VIS` (VISAMP and VISPHI, `--vis`), as well as V² and T3PHI;
  - calibrators in another instrument mode (spectral resolution, polarisation) or reduced with another `vis-correction-sc` are refused (review §7c.5); differences in DIT or pipeline version are warned about.

  Checked against that script on real data: with the original choices (`--tf-error scatter`), the results agree except for two deliberate changes: each baseline's calibrator estimate takes its own mean time, and calibrator closure phases are interpolated as unit vectors (checked privately against the original script).
- **Design notes moved from virgil:** the GRAVITY calibration review, the calibrator-PCA plan and the questions for the GRAVITY team. New: PIONIER systematics notes and a MATISSE stub.

## Next, roughly in order
### Calibration (no core dependency)

### Error models (some need virgil 6d)
1. **Per-DIT covariances** for GRAVITY from `P2VMRED`/`ASTROREDUCED` (review §12.1; Nowak et al. 2020, Kammerer et al. 2020).
2. **TF error as a per-(night, baseline) gain mode**: done (`--tf-modes`, #25).
3. **Dual-field GRAVITY templates**: the code is in (#26); the template awaits a cluster job (`gravity_cal_pca_dual`).
4. **Check that 2023 CAL V² > 1 (FT coherence correction) cancels between SCI and CAL.**

### Instrument facts
5. **GRAVITY fibre injection** (UT and AT fields, tip-tilt jitter) for virgil's generic primary beam (review §1c.2): done (`gravity.fibre_fwhm_mas`, `gravity.fibre_coupling`); the beam itself is virgil's.
6. **Wavelength priors** per instrument mode (review §6): provisionally N(0, 2×10⁻⁴) in scale plus a 0.5 nm offset for GRAVITY MEDIUM/HIGH, coded as `gravity.wavelength_prior`; final values await the GRAVITY team.
7. **Tellurics and flux calibration** for `OI_FLUX` (review §1c.3).

### MATISSE
The literature review is done ([`design/matisse_calibration_review.md`](design/matisse_calibration_review.md)). Its first work item, explicit product categories in `classify`, is done. The staged implementation plan is [`design/matisse_plan.md`](design/matisse_plan.md). The rest, in the review's order:
- a MATISSE header reader (band, resolution, BCD, chopping, `corrFlux`, spectral binning, pipeline versions, merges), refusing `corrFlux` products as V² and flagging N-band phases from pipelines 1.5.2–1.x (180° convention) (done: plan Stage 1, `virgil_vlti.matisse` and `virgil-vlti-inspect`);
- calibration per BCD position and chopping status, with orientation-normalised matching, then BCD combination with propagated errors;
- N-band correlated-flux calibration with calibrator spectra;
- MATISSE nuisance widths, instrument facts (R per mode after spectral averaging, pinhole sizes, wavelength windows), and validation against `mat_cal_oifits`.

For virgil core it recommends `AMPTYP` checks in `read_oifits` (done: virgil#150), smearing with an explicit spectral response (MATISSE's `EFF_BAND` is the channel spacing), per-telescope low-rank modes, and a pinhole primary beam.

## Validation
Independent checks of this code are made in a private validation suite, built on separate code; their status is in [PROVENANCE.md](PROVENANCE.md). The most important one: the same products calibrated here and by the instrument pipeline (`gravity_viscal`), which should agree when the pipeline's interpolation is used.
