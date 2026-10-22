# PIONIER: what we know about its calibration errors

Status: **notes**, 2026-10-04. Collected from analyses done elsewhere; nothing here is implemented in this repository yet.

## Sources
- **The systematics study:** `pionier_iwcar_systematics.ipynb` (a private IW Car notebook, not in this repository). It used IW Car, PIONIER, the 23 files from 2019-11-25 to 2020-02-14 (ESO 1103.D-0797), with virgil at commit `823ccab`. **That commit predates virgil#120**, which keeps only the independent closure phases, so the notebook's own closure-phase numbers ("What the data show" and the "Notebook" column below) over-count by about 4/3. The recomputed columns in "Recomputed with `virgil_vlti.pairs`" are after #120 (independent triangles, ±1/3 correlations), and supersede them.
- **Real-data analyses** in a private repository (IRAS 08544-4431, IW Car, HR 4049).
- **The GRAVITY review** ([`gravity_calibration_review.md`](gravity_calibration_review.md)), whose footgun 9 (closure-phase redundancy) and §4 (correlated gains) apply equally to PIONIER.

## What the data show
The test compares consecutive exposures of one configuration, about 20 min apart (13 pairs). Over that interval the uv points barely move, so the difference between the two exposures is almost pure noise. Its covariance is fitted with independent noise, a gain per (file, baseline) and a gain per (file, telescope), all shared by the six channels.

- **V²: about the right total size, but the wrong structure.**
  - Only about 0.32σ is independent between channels.
  - The rest is an achromatic offset per file: about 1.3% per telescope plus 0.9% per baseline, ΔlnL = +352 over independent errors. In the joint fit the telescope term carries about 80% of the correlated variance.
  - Together these reproduce the pipeline's total per-channel error of about 2.4%.
- **Closure phases.** The errors are about 1.6× too large (`s_white` ≈ 0.6), with a common offset per (file, triangle) of about 0.2°. As expected, closure phases are almost free of telescope-based errors.
- **In the fits.** In every image reconstruction, 64–77% of the V² residual χ² sits in the six-channel mean of each (file, baseline). Independent errors would put 17% there.
- **A lower bound.** Consecutive exposures share calibrators, so any calibration error they share cancels in the difference. The correlated amplitudes are lower bounds.

The overall scale is not consistent between targets. virgil's [`regulariser_weight_selection.md`](https://github.com/benjaminpope/virgil/blob/main/design/regulariser_weight_selection.md) reports χ² per point of 0.5–0.8 for smooth images on real PIONIER data, so the errors there look conservative. HR 4049 instead gave `error_scale` ≈ 1.13 (errors about 13% too small). Both predate #120.

## Recomputed with `virgil_vlti.pairs` (2026-10-04)
The pair test is now a tool (`virgil_vlti.pairs`, `virgil-vlti-pairs`). Three things changed from the notebook:
- **The gain terms are built from a reference V²** (the pair's mean, or a model), not from each frame's noisy V². With the noisy values the gain terms are nearly full rank and can absorb the white noise; on synthetic data with many channels the fit then collapsed to no white noise at all. With six PIONIER channels the effect was small.
- **Closure-phase white noise has the ±1/3 structure** between triangles that share a baseline, and only independent triangles are kept (per channel, from the valid ones), as in virgil since #120.
- **Pairs are disjoint.** (A, B) and (B, C) share B's noise, so their likelihoods cannot simply be added.

On IW Car (23 files, up to 2020-02-14). The tool uses disjoint pairs, so frames shared between overlapping pairs (2 of the notebook's 13) are not counted twice. That leaves 11 pairs, and the log-likelihoods add up to a proper joint likelihood:

| Term | Notebook (13 overlapping pairs, model subtracted) | 11 disjoint pairs, no model | 11 disjoint pairs, model subtracted |
|---|---|---|---|
| V² white (× σ) | 0.32 | 0.32 | 0.31 |
| V² per baseline | 0.9% | 1.2% | 0.9% |
| V² per telescope | 1.3% | 1.1% | 1.1% |
| V² ΔlnL over independent | +352 | +275 | +250 |
| CP white (× σ) | ≈ 0.6 | 0.68 | 0.64 |
| CP offset per triangle | ≈ 0.2° | 0.40° | **0.21°** |
| CP ΔlnL for the offset | — | +48.9 | +19.5 |

**The model** is the notebook's single-ring parametric model (a private script). It was refitted on a cluster on 2026-10-04 (virgil eabcbcc, virgil-vlti 9d70f1a) and agrees with the PIONIER session's own run: ring 15.7 mas across, i = 53°, PA 164°; secondary at (1.05, −1.72) mas. Its predicted change between paired frames is subtracted (`pair_test(..., model=...)`).

**Conclusions:**
- **V²:** about a third of the reported error is independent between channels. The rest is an achromatic gain of about 1% per telescope plus about 1% per baseline, common to the six channels. With or without the model, and with or without overlapping pairs, the numbers hardly move.
- **Closure phases:** the reported errors are about 1.5× too large. The offset per (frame, triangle) is **0.21°** once the target's own change is removed. Without the model it doubled to 0.40°, because IW Car's closure phases change in the ~20 min between frames. The notebook's 0.2° stands.
- As before, these are lower bounds: consecutive frames share calibrators.

## What this means for virgil-vlti
1. **The pair test as a tool:** done (above).
2. **Default priors for virgil's Stage 6d gains.** Gains on |V| per telescope of about 0.55% (half of 1.1% on V²) and per baseline of about 0.45% (half of 0.9%), as fitted widths with these as starting values. Once 6d exists these become a small, versioned table here, like the GRAVITY templates of [`gravity_calibrator_pca.md`](gravity_calibrator_pca.md).
3. **Until 6d.** Adding the calibration error in quadrature (about 2% of V², 0.21° on closure phase) is the honest stop-gap. It treats correlated errors as independent, though, so it still overweights common offsets.
4. **Recompute after #120:** done, with and without the target's model (above).

## Open
- Whether the same structure holds on other PIONIER programmes and epochs (GRAVITY's calibrator PCA could be run on PIONIER calibrators too).
- PIONIER products are usually calibrated by `pndrs` before we see them. Re-calibrating from the raw transfer function with `virgil-vlti-calibrate` would allow a leave-one-out error, but has not been tried.
