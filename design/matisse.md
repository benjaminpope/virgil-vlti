# MATISSE: status and what to find out

Status: **answered**, 2026-10-04: the questions below are answered in [`matisse_calibration_review.md`](matisse_calibration_review.md), which also lists the work items. The staged implementation plan is [`matisse_plan.md`](matisse_plan.md). This note records what existed before that review.

## What exists
- **Data.** Three ν Hor files, L/M band, LOW resolution, chopped, calibrated according to their names (`2023-12-*_nu.Hor_K0G2D0J3_IR-LM_LOW_Chop_cal_oifits_0.fits` local only, not in this repository).
- **Code.** An early scratch notebook, a private scratch notebook, which only loads the files with AMICAL's reader. virgil's `coverage.vlti_oidata` makes synthetic coverage modelled on MATISSE LOW in L.
- **Known limitation in virgil.** Bandwidth smearing is not modelled. At LOW resolution (R ≈ 30) this limits the field of view to about Rλ/B (virgil `design/imaging_plan.md`, Stage 4).

## A review should answer, as the GRAVITY one did
1. **Products.** Which observables MATISSE writes in each band (L/M and N), and which are calibrated: V², the correlated flux, closure phases and differential phases. Its OIFITS layout: `INSNAME`s, `PHITYP` (differential phases are likely), and whether chopped and non-chopped data are separate tables.
2. **Incoherent and coherent processing.** What each gives, and the error model of each.
3. **Calibration.** How the pipeline (and the consortium) builds the transfer function, how photometry and chopping enter, and how errors are estimated. Is there a calibration term?
4. **The N band.** Thermal background, the strong variation of the transfer function, and what that means for calibration nuisances.
5. **Correlations.** Between channels, between the BCD configurations, and between baselines.
6. **Closure phases.** Four telescopes give four triangles, of which three are independent. virgil now handles this (#120), but check how MATISSE's errors are built.
7. **Smearing.** Channel widths in LOW and MED, and when bandwidth smearing must be modelled.

## Prerequisites in virgil
The same generic interfaces as for GRAVITY: differential phases (Stage 6a), smearing (6a), and low-rank calibration nuisances (6d).
