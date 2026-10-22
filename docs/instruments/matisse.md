# MATISSE

For MATISSE, virgil-vlti reads product headers and refuses what must not be calibrated: already-calibrated or merged products, correlated fluxes posing as visibilities, and N-band phases whose sign is uncertain. `virgil-vlti-inspect` reports all of this per file. It does **not** yet calibrate MATISSE data in a way that has been validated: there is no calibration per BCD position, no BCD combination, no N-band correlated-flux calibration, and no comparison with ESO's `mat_cal_oifits`.

MATISSE is outside the maintainer's core expertise, so every convention below cites its source, ambiguous conventions refuse by default, and all of it awaits review by a MATISSE expert. The staged plan is [MATISSE implementation plan](../design/matisse_plan.md) (its §7.1 lists what is unverified); the literature review is [MATISSE calibration review](../design/matisse_calibration_review.md).

Sources, as in [Provenance](../PROVENANCE.md): [MPM] the MATISSE pipeline manual 2.5.0; [UM] the MATISSE user manual 1.0; [MT] the MATISSE consortium tools; [Lo22] Lopez et al. 2022; [GR22] Gámez Rosas et al. 2022; [ESO-news] the ESO MATISSE news page (read 2026-10-07).

## Implemented

All in [`virgil_vlti.matisse`](../api/matisse.md) and [`virgil_vlti.products`](../api/products.md) (Stages 0 and 1 of the plan).

- Product categories: `TARGET_RAW_INT` (science) and `CALIB_RAW_INT` (calibrator) accepted; `*_CAL_INT` (already calibrated) refused ([MPM] §7.5, §9.1.2).
- Setup from headers: band from `ESO DET NAME`, resolution from `ESO INS DIL NAME` (LM) or `DIN NAME` (N), BCD state from `ESO INS BCD1/BCD2 NAME`, telescopes, DIT, template (keys read from real headers; band and resolution as plan §4 reads [MPM]).
- Chopping from `ESO ISS CHOP ST`, never `ESO DET CHOP ST` (from real headers only; neither manual names either key: **unverified**).
- Spatial filter from `ESO INS SFL NAME` (LM) or `SFN NAME` (N), compared only (from real headers only: **unverified**).
- A product whose `ESO CFG BCD MODE` disagrees with its BCD keys is refused (new here; which key is which device is undocumented).
- BCD positions in device order; naming, combining or calibrating across positions refused until the `bcd_key_order` fact is recorded (swaps from [MPM] §4.1 and [UM] §2.2.1; the conflict between [MPM] §11.9 and [MT] `mat_removeBCD` is in plan §3.1: **unverified**).
- BCD peak tables: which `OI_VIS2` and `OI_T3` row holds which baseline or triangle at each position, with closure-phase signs (derived in plan §3.1 from the [MPM] §4.1 swaps).
- N-band phase convention by `mat_raw_estimates` version: `flipped` in [1.5.2, 2.0.0), `standard` from 2.0.0 and in L/M, `unknown` otherwise; refused by default, `--n-phase` and `--n-phase-unknown` to negate or keep ([MPM] §2.5; [GR22] Methods §4.3, §4.8; [ESO-news] December 2023; the 1.5.2 lower bound and VISPHI are **unverified**).
- `corrFlux=true`, or unrecorded, refused for V², T3AMP and VISAMP ([MPM] §7: `VIS2DATA` may hold squared correlated fluxes; the review, footgun 1).
- Merged products refused: `INT_TIME` above 1.5 × `ESO DET SEQ1 EXPTIME`, consortium or `_CAL_INT` names, `mat_mergeOifits` in `HISTORY`; one baseline at several MJDs is reported only (plan §4 Stage 1; names from [MT] `mat_mergeAllOiFits.py`).
- Calibrator compatibility: band, resolution, BCD keys, telescopes, spatial filter, chopping and DIT must agree, and in N the phase convention; `mat_raw_estimates` parameters `corrFlux`, `coherentAlgo`, `cpAlgo`, `spectralBinning`, `spectralAverage`, `replaceTel` must agree; the pipeline version is warned about ([MPM] §2.1 for the 2.5.0 renames, which are never equated; [MPM] §9.1.4 for `replaceTel`).
- `replaceTel` ≠ 0 warned about ([MPM] §9.1.4: one telescope's photometry replaced by the others' mean).
- Group key (night, band, resolution, DIT, BCD keys, chopped), stricter than [MT] `mat_autoCalib.py` `findClosestCal` (`group_key`; plan §3.3). Not yet used by `virgil-vlti-calibrate`.
- `virgil-vlti-inspect`: one row per file with category, band, resolution, BCD, chopping, pipeline version, phase convention, merge evidence, `replaceTel` and a verdict ([command-line tools](../command_line.md)).

**What `virgil-vlti-calibrate` does with MATISSE today.** It screens every MATISSE product with the refusals above, then groups by night only. On the default path a night with several BCD positions fails compatibility and nothing is calibrated, and a night with one position goes through the GRAVITY method unchanged, which has not been validated for MATISSE. Two options bypass the MATISSE compatibility rules altogether: `--match-gravity-viscal` (its setup check reads GRAVITY-only keys, which a MATISSE header lacks, so calibrators at different BCD positions are pooled into one wrong transfer function, with no error raised) and `--cp-only-when-uncalibrated` (closure phases written with no BCD handling). Do not use any of these on MATISSE data for science.

## Not implemented / planned

In the order of the [plan](../design/matisse_plan.md):

- Stage 2: validation data and pipeline reductions, which record the instrument facts `bcd_key_order` and `sta_index_frame`.
- Stage 3: synthetic MATISSE products and an orientation-neutral core.
- Stage 4: L/M V² and closure phases calibrated per BCD position and chopping status, then combined across positions with propagated errors; a `--match-mat-cal-oifits` parity mode.
- Stage 5: the MATISSE error model (pair test, calibrator PCA, templates, photometric modes after [Lo22] eq. A.9). The pair test and calibrator PCA skip MATISSE products until then.
- Stage 6: instrument facts and default widths (R per mode after spectral averaging, pinhole sizes, wavelength windows).
- Stage 7: N-band correlated-flux calibration with calibrator spectra.
- Stage 8: N-band V², differential phases, MED and HIGH resolution.

For virgil core the review also recommends smearing with an explicit spectral response (MATISSE's `EFF_BAND` is the channel spacing), per-telescope low-rank modes and a pinhole primary beam; `AMPTYP` checks in `read_oifits` are done (virgil#150).

## Checked against

| What | Compared with | Result | Status |
|---|---|---|---|
| Calibration | `mat_cal_oifits` (`tfInterp=0`, per BCD position; independent check 2) | No validated MATISSE calibration here yet | **Not yet** |
| Header reader | Three real L/M LOW chopped `OUT-OUT` products of a public target (`matisse/2.0.0`) | Read as expected; refused as merged (`INT_TIME` 120.46 s against 61.8 s) and as `*_CAL_INT` | **Passed** (unit tests on header literals) |
| BCD peak tables | [MPM] §11.9's closure-phase sign row | Consistent; a simulated binary's closure phases recovered to 1e-6° through the tables. This catches transcription errors only: the simulation uses the same input-swap model | **Partial** |
| BCD position order | [MT] `mat_removeBCD` (`libPostTools.py`) | Contradicts the order inferred from [MPM] §11.9 if read with the same row order; unresolved, so cross-position work is refused | **Not yet** (unverified) |
| N-band phase convention | [MPM] §2.5, [GR22] Methods §4.3, §4.8, [ESO-news] | Version boundaries as cited; the 1.5.2 lower bound is not confirmed by [MPM] §2.11 | **Partial** (unverified before 1.5.2 and for VISPHI) |
| Every convention above | A MATISSE expert | Review pending | **Not yet** |

## Background

- [MATISSE calibration review](../design/matisse_calibration_review.md): the literature, pipeline and consortium-tool review.
- [MATISSE implementation plan](../design/matisse_plan.md): the stages, conventions and open questions.
- [MATISSE status](../design/matisse.md): the note that existed before the review.
