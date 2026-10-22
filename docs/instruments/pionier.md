# PIONIER

For PIONIER, virgil-vlti measures the structure of the errors in products that are already calibrated, and turns it into nuisance modes for virgil. It does **not** calibrate PIONIER data: `virgil-vlti-calibrate` refuses any `INSTRUME` other than GRAVITY and MATISSE. PIONIER products are normally calibrated by ESO's PIONIER pipeline, `pndrs`, before they reach us.

The tool is the pair test: two consecutive exposures of one target, about 20 minutes apart, see almost the same uv points, so their difference is almost pure noise, and fitting its covariance shows which part of the reported error is independent between channels and which part is a gain shared by a telescope or a baseline. The tutorial is [The pair test](../pair_test.md); the full analysis is in [PIONIER systematics](../design/pionier_systematics.md).

## Implemented

- The pair test on calibrated PIONIER products (`virgil-vlti-pairs`, [`virgil_vlti.pairs`](../api/pairs.md)): white noise per channel, gains per (frame, baseline) and (frame, telescope) on V², and closure-phase offsets per (frame, triangle) ([tutorial](../pair_test.md)).
- PIONIER's one-frame products and per-file `INSNAME` (which carries the wavelength solution) handled: by default each file's only instrument is used ([pairs](../api/pairs.md)).
- Gain covariances built from a reference V² (the target's model, or the pair's mean), not the noisy data ([pairs](../api/pairs.md)).
- Closure-phase white noise with the ±1/3 correlation between triangles sharing a baseline; only independent triangles kept ([pairs](../api/pairs.md)).
- Disjoint pairs only, so log-likelihoods add up to a joint likelihood ([pairs](../api/pairs.md)).
- Subtracting a target model's predicted change between paired frames (`pair_test(..., model=...)`) ([pairs](../api/pairs.md)).
- Flat gain templates from the pair test's widths (`virgil-vlti-template --from-pairs`), loaded into virgil by `to_oidata(..., template=...)` ([templates](../api/templates.md), [loading](../api/load.md)).

## Not implemented / planned

- Calibration of PIONIER products from their raw transfer functions with `virgil-vlti-calibrate`: not tried ([PIONIER systematics](../design/pionier_systematics.md), Open).
- A shipped PIONIER calibration template (a small versioned table of gain widths, once virgil's correlated gains are in use) ([PIONIER systematics](../design/pionier_systematics.md), item 2).
- The calibrator PCA on PIONIER calibrators ([PIONIER systematics](../design/pionier_systematics.md), Open).
- Whether the IW Car error structure holds on other programmes and epochs ([PIONIER systematics](../design/pionier_systematics.md), Open).

## What the pair test found

On PIONIER data of IW Car (23 files, ESO programme 1103.D-0797, six spectral channels), on 11 disjoint pairs, with the target's own model subtracted:

| Term | Value |
|---|---|
| V² white (× σ) | 0.31 |
| V² per baseline | 0.9% |
| V² per telescope | 1.1% |
| CP white (× σ) | 0.64 |
| CP offset per triangle | 0.21° |

In words: about a third of the reported V² error is independent between channels; the rest is an achromatic gain of about 1% per telescope plus about 1% per baseline. The closure-phase errors are about 1.5 times too large, with a common offset per (frame, triangle). These are lower bounds, because consecutive frames share calibrators. The numbers hardly move with or without the model, but the closure-phase offset doubles to 0.40° without it, because IW Car's closure phases change between frames.

**Using this in virgil.** Until correlated gains are in use, the honest stop-gap is to add the calibration error in quadrature (about 2% of V², 0.21° on closure phase), which still overweights common offsets. With gains, about 0.55% on |V| per telescope and 0.45% per baseline are the suggested starting widths ([the note](../design/pionier_systematics.md); [Calibration templates for virgil](../templates.md) shows how modes and widths go in).

## Checked against

| What | Compared with | Result | Status |
|---|---|---|---|
| Calibration | `pndrs` (ESO's PIONIER pipeline) | No PIONIER calibration here, so nothing to compare | **Not yet** |
| Closure-phase correlation ±1/3 and `independent_triangles` | Kammerer et al. 2020, §2.2; Monte Carlo of 200,000 draws of baseline-phase noise for 4, 6 and 7 telescopes (independent check 6) | Agree | **Passed** |
| Telescope-shared gain terms | Kammerer et al. 2020, §2.2 (telescope-based correlations) | Injected widths recovered on synthetic data | **Passed** (internal check only) |
| Pair-test widths | Bias and coverage on independently simulated products (independent check 5) | — | **Not yet** |
| IW Car numbers | The earlier private notebook on the same 23 files | V² white 0.32 against 0.31; V² per telescope 1.3% against 1.1%; per baseline 0.9% both; CP offset about 0.2° against 0.21° | **Passed** (our own analysis, not a pipeline comparison) |

Sources and tests for each row are in [Provenance](../PROVENANCE.md), under "Error structure".
