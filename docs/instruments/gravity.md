# GRAVITY

GRAVITY is the instrument virgil-vlti supports best. Give it the *uncalibrated* products of ESO's GRAVITY pipeline (`*_SCI_VIS` and `*_CAL_VIS`, or the in-field calibrator frames of dual-field data) and it writes calibrated OIFITS files (V², T3AMP, T3PHI, and optionally `OI_VIS`), with the calibration's shared errors as nuisance modes for virgil to marginalise. It can also reproduce ESO's own calibration (`gravity_viscal`) exactly, learn the systematics calibration leaves behind from many calibrator stars, and ship them as templates.

It does not reduce raw frames (use ESO's pipeline for that), never reads the fringe-tracker (FT) tables, and does not fit models: that is virgil's job.

Start with the tutorial [Calibrating GRAVITY products](../calibration.md); the command-line options are on the [home page](../index.md) and in [Command-line tools](../command_line.md).

## Implemented

**Calibration** (`virgil-vlti-calibrate`, [`virgil_vlti.calibrate`](../api/calibrate.md))

- Transfer function from calibrator stars, interpolated in time to each science frame: linear (default) or the pipeline manual's Gaussian weights (`--interp gaussian`) ([tutorial](../calibration.md)).
- Calibrators found by `ESO PRO CATG`, by name (`--cal NAME[=DIAM±ERR]`), or in dual field by the fibre offset (`--dual-field`, the fringe-tracking star as an in-field calibrator) ([products](../api/products.md)).
- V², T3AMP and T3PHI calibrated; `OI_VIS` (VISAMP, VISPHI) with `--vis` ([calibration](../api/calibrate.md)).
- Closure phases averaged and interpolated as unit vectors, so nothing wraps at ±180° ([calibration](../api/calibrate.md)).
- Each baseline's calibrator estimate at the mean time of its own frames, after GRAVITY's per-baseline frame rejection ([calibration](../api/calibrate.md)).
- Transfer-function error as a leave-one-out prediction error (default), the night's scatter, the calibrators' statistical error, or none (`--tf-error`) ([calibration](../api/calibrate.md)).
- Transfer-function error as one nuisance mode per night and baseline or triangle (`--tf-modes`), with a smoothed channel shape (`--tf-mode-shape`, default `poly1`) ([tutorial](../calibration.md)).
- Calibrator uniform-disk diameters, from `--cal`, `--diameter`, or the pipeline's `ESO QC TF CAL_SC DIAM` cards ([calibration](../api/calibrate.md)).
- Calibrator diameter errors as nuisance modes, absolute and model-relative (`.modes.npz`) ([writing](../api/write.md)).
- Calibrators in another resolution, polarisation, wavelength grid or `vis-correction-sc` refused; unrecorded keys refused unless `--allow-unrecorded` ([GRAVITY facts](../api/gravity.md)).
- Optional spectral smoothing of the amplitude transfer function, as `gravity_viscal` does it (`--smooth-tf DEG`) ([calibration](../api/calibrate.md)).
- Closure-phase-only output for nights with no calibrator (`--cp-only-when-uncalibrated`), with **uncalibrated** T3PHI, and for named files (`--cp-only`) ([home page](../index.md)).
- Science-to-calibrator `OI_FLUX` ratio and per-telescope injection scatter (`flux_ratio`, `injection_scatter`) ([calibration](../api/calibrate.md)).
- Predicted anisoplanatic V² loss of each off-axis dual-field exposure, reported in `index.csv` and warned about above 2% ([GRAVITY facts](../api/gravity.md)).
- Calibrated files with provenance cards (`HIERARCH VIRGIL …`) and OIFITS 2 columns only ([writing](../api/write.md)).

**ESO parity** ([`virgil_vlti.viscal`](../api/viscal.md))

- `--match-gravity-viscal`: a port of `gravity_viscal` (GRAVITY pipeline 1.9.7), to reproduce ESO's `*_SCI_VIS_CALIBRATED` products; see [Matching ESO's gravity_viscal](#matching-esos-gravity_viscal) below.

**Error models and templates**

- Calibrator PCA (`virgil-vlti-calpca`): leave-one-out calibrations of calibrators, screened, clipped and decomposed by weighted PCA, with held-out widths ([calpca](../api/calpca.md)).
- Calibration templates (`virgil-vlti-template`, `--template NAME`): per-telescope and per-baseline log |V| gains and closure-phase offsets, refused for data in another configuration ([tutorial](../templates.md), [templates](../api/templates.md)).
- Three shipped single-field templates: `gravity_single_medium_combined_2023-2025` (the MEDIUM pilot), `gravity_single_high_combined_2023-2025` (HIGH, ATs only) and `gravity_single_medium_combined_2023-2025_nowr` (a MEDIUM rebuild on a smaller calibrator set) ([calibrator PCA note](../design/gravity_calibrator_pca.md)).
- Flat gain templates from a target's own pair test (`virgil-vlti-template --from-pairs`) ([templates](../api/templates.md)).
- The pair test on GRAVITY frames (`virgil-vlti-pairs`) ([tutorial](../pair_test.md)).
- Empirical null datasets: each calibrator exposure calibrated as a science target by the night's other calibrators, written as OIFITS, with block-resampled versions ([nulls](../api/nulls.md)).
- Loading calibrated files and every mode into a virgil `OIData` (`to_oidata`; needs virgil-astro 0.4.0 or virgil's main branch) ([loading](../api/load.md)).

**Instrument facts** ([`virgil_vlti.gravity`](../api/gravity.md))

- Field configuration from the headers: `SINGLE`, `DUAL_ON`, `DUAL_OFF`, `DUAL_WIDE` (`field_mode`).
- Fibre beam and tip-tilt jitter: `fibre_fwhm_mas`, `fibre_coupling`, `fibre_response`.
- Smearing diagnostics: `spectral_fwhm`, `bandwidth_smearing_v2_factor`, `time_smearing_v2_factor`.
- Calibrator K-band SED slopes and Brγ and telluric masks: `calibrator_sed_slope`, `brgamma_mask`, `telluric_mask`.
- Provisional wavelength priors for MEDIUM and HIGH: `wavelength_prior`.

## Not implemented / planned

- Reducing raw frames: stays with ESO's pipeline.
- Per-DIT covariances from `P2VMRED` or `ASTROREDUCED` products ([plan](../PLAN.md), error models 1).
- Dual-field calibration templates: the code is in (`virgil-vlti-calpca --field`), the template awaits a run ([plan](../PLAN.md), error models 3).
- A check that 2023 calibrator V² above 1 (the FT coherence correction) cancels between science and calibrator ([plan](../PLAN.md), error models 4).
- Tellurics and flux calibration of `OI_FLUX` beyond the ratio and the masks above ([plan](../PLAN.md), instrument facts 7).
- Final wavelength priors: the current ones are provisional, pending the GRAVITY team; none for LOW ([questions for the GRAVITY team](../design/questions_for_gravity_team.md)).
- A HIGH-resolution template for the UTs, and templates for split polarisation or LOW resolution ([calibrator PCA note](../design/gravity_calibrator_pca.md)).
- A correction for the coherence loss of off-axis dual-field targets: only its predicted size is reported ([review](../design/gravity_calibration_review.md) §7).
- Channel correlations: GRAVITY channels are about twofold oversampled, so bin by 2 or more ([home page](../index.md)).

## Checked against

Full rows, with sources and tests, are in [Provenance](../PROVENANCE.md). "Independent" means the private validation suite, which shares no code with this package.

| What | Compared with | Result | Status |
|---|---|---|---|
| `--match-gravity-viscal` | ESO's calibrated products from `gravity_viscal` (pipeline 1.9.7; source `gravi/gravi_tf.c`, `gravi/gravi_vis.c`), 4 public cases with public inputs | Every observable within about 1e-11 relative (amplitudes) and 4e-8° (phases), errors identical; independent check V² ≤ 1.2e-11 | **Passed** |
| Default calibration (`virgil-vlti-calibrate`) | The same 4 ESO products (independent check 1, reported 2026-10-04) | Broadband V² within 0.12% (median), every row within 1.1%; closure phases within 0.06° | **Passed**, for one calibrator per file |
| Default calibration, V² errors | The same products (independent check 1) | Equal to ESO's when checked, but that was before the one-calibrator floor: with one calibrator the default errors now also include the calibrator's propagated error, which `gravity_viscal` does not add, so they exceed ESO's. Not rerun since | **Not yet** (for today's default) |
| Default calibration, other outputs | The same products | Diameter errors, `OI_VIS`, T3AMP, flux ratio and several calibrators per night not yet compared | **Not yet** |
| Gaussian time weights (`--interp gaussian`) | GRAVITY pipeline manual 1.11.0, §9.20, eqs. 78–80 | Follows the manual's kernel, which differs from the 1.9.7 source (an exponential in \|ΔT\|); interpolation between several calibrators not yet compared with ESO | **Partial** |
| `--smooth-tf` | `gravi_vis_fit_amp` in `gravi/gravi_vis.c` (1.9.7) | Same normal equations (unit test); on one public star the channel scatter of the ratio to ESO falls from 19.7% to 0.2% | **Passed** (development check, not independent) |
| Calibrator compatibility rules | Pipeline manual §6 and §9.10, via the GRAVITY review §7c.5 | Unit tests on headers encode the review's reading; the manual was not re-read | **Partial** (internal check only) |
| Which `ESO PRO RECn` holds `gravity_vis` | A real GRAVITY product header | Only synthetic multi-REC headers tested | **Not yet** (unverified) |
| Uniform-disk V², \|V\| and T3AMP | A first-principles closed form (independent check 3) | Agree to 1e-13 at 0.3, 1.2 and 3 mas | **Passed** (the validator saw these function bodies; see Provenance) |
| Anisoplanatic V² loss | GRAVITY+ Collaboration 2022, eqs. 1, 2, 4, 5 | Reproduces the paper's HD 48543B example | **Partial**: Elhalkouj et al. 2008's full expression and the paper's Figs. 9–11 not yet (independent check 4) |
| Diameter modes | Monte Carlo over calibrator diameters (independent check 7) | Finite-difference unit tests only | **Not yet** |
| `--tf-mode-shape poly1` default | Our own Monte Carlo (GLS pull rms about 1.3 for transfer-function colours 0–0.5) | Not compared with ESO/VLTI documentation of transfer-function chromaticity | **Not yet** (unverified) |
| Fibre beam (Gaussian of FWHM λ/D) | Lacour et al. 2019 §2 (about 60 mas at K on the UTs), via the GRAVITY review §1c.2 | 55 mas at 2.2 µm on the UTs; not compared with a measured injection profile | **Partial** (unverified) |
| Tip-tilt jitter in the beam | GRAVITY Collaboration 2020 App. A.3 (≈15 mas r.m.s. per axis) | Monte Carlo over pointing errors | **Passed** (internal check only) |
| Spectral kernel (1.4 channels) | PMOIRED's GRAVITY `wl kernel` (Mérand 2022, code), via the review §8a | Unit tests against numerical integration | **Partial**: no fit on real data near the smearing limit yet |
| Wavelength priors | Pipeline manual §5.1 (0.5 nm absolute, 0.1 nm between baselines) | Lookup and refusal unit-tested | **Partial**: provisional, awaiting the GRAVITY team |
| Calibrator SED slopes | Pecaut & Mamajek 2013, Table 5 (dwarfs) | Giant and supergiant temperatures entered from memory | **Not yet** (unverified) |
| Weighted PCA | Bailey 2012 (PASP 124, 1015), as recalled | Recovers an injected shape and width on synthetic data | **Not yet**: no comparison with an established implementation |
| Calibrated files read by other OIFITS readers | PMOIRED, CANDID (independent check 8) | — | **Not yet** |
| Earlier target-specific calibration script | The script's own outputs on real data | Agree except for two deliberate changes (per-baseline mean times, unit-vector closure phases) | **Passed** (checked privately) |

## Matching ESO's gravity_viscal

`virgil-vlti-calibrate --match-gravity-viscal` calibrates exactly as ESO's `gravity_viscal` recipe (GRAVITY pipeline 1.9.7) does, to reproduce ESO's calibrated products (`*_SCI_VIS_CALIBRATED`). The module is `virgil_vlti.viscal` ([API](../api/viscal.md)). It was written from the pipeline's source (`gravi/gravi_tf.c` and `gravi/gravi_vis.c`). Where the source and the pipeline manual differ, it follows the source. It is a separate design from the default calibration in `virgil_vlti.calibrate`.

```bash
virgil-vlti-calibrate reduced/ --out calibrated --match-gravity-viscal
```

Its own defaults are used: `--window`, `--bin`, `--interp`, `--tf-error` and `--smooth-tf` do not apply. The parity options mirror `gravity_viscal`'s own:

| Option | Default | `gravity_viscal` parameter |
|---|---|---|
| `--delta-time-calib S` | 3600 | `delta-time-calib` (s) |
| `--maxdeg-tfvis DEG` | 5 | `maxdeg-tfvis-sc` |
| `--nsmooth-tfvis N` | 0 | `nsmooth-tfvis-sc` |
| `--no-smoothing` | off | `smoothing=FALSE` |
| `--force-calib` | off | `force-calib`: keep calibrators with another setup |

Per instrument and for every calibrator frame (row) separately:

1. **Transfer function.** VIS2 and its error are divided by the uniform-disk V², VISAMP and its error by |V|. T3AMP is not corrected for the calibrator's size, and phases are used as they are.
2. **Spectral smoothing** (not in LOW resolution): flags are extended and shared by a table's amplitude and phase; an optional running average (`nsmooth`); then a polynomial of degree `maxdeg` in the channel index, weighted by σ⁻², fitted to VIS2, VISAMP, T3AMP, VISPHI and T3PHI.
3. **Weights.** Each calibrator frame of the same baseline or triangle gets `W = exp(−2|ΔT|/Δ)/σ²`, with σ floored at 2% of its median transfer function for amplitudes and at 0.2° for phases. The manual (§9.20, eq. 79) gives a Gaussian in ΔT; the source uses the exponential in |ΔT|, and so does the port. Across UT midnight `TIME` jumps, so calibrators on the other side of midnight are effectively ignored, which is also reproduced.
4. **Calibration.** Amplitude over the transfer function, error over the transfer function (no transfer-function error: the source marks it `FIXME`), and phase minus the weighted mean phase.
5. **Setups.** A calibrator is used only if its spectral resolution, both polarisation modes and both DITs match the science, unless `--force-calib`.

The calibrator diameter is the one `gravity_vis` recorded (`ESO QC TF CAL_SC DIAM`), unless given; its error is not propagated.

## Background

- [GRAVITY calibration review](../design/gravity_calibration_review.md): the literature and pipeline review behind the design.
- [GRAVITY calibrator PCA](../design/gravity_calibrator_pca.md): how the templates were built, and their caveats.
- [Questions for the GRAVITY team](../design/questions_for_gravity_team.md): what is still open.
