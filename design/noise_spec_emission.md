# Emitting a noise spec from the pair test

Status (2026-10-18): **design, nothing implemented.** Counterpart of virgil's
`design/vlti_noise_in_pipelines.md`, which defines the spec and how pipelines
use it; read that first. Target virgil 0.4.1.

## What virgil-vlti adds

`virgil_vlti.pairs.fit_noise` returns `Fit(widths, delta_lnl)` per nested model,
with widths keyed by term name. `pair_test` returns `PairTest(n_pairs, v2, cp)`.
Neither is a record a pipeline can read. This note adds one function and one
CLI option; virgil-vlti stays independent of virgil's pipeline code (it already
uses virgil only for the data classes).

## `noise_spec(test, *, v2_model="+ baseline + telescope", cp_model="+ triangle offset", **provenance) -> dict`

Builds the spec of virgil's design (consumed by the pipeline setting
`noise_spec=`, CLI `--noise-spec`) from a `PairTest`. The V² and closure-phase
fits are separate (`pairs.V2_MODELS`, `pairs.CP_MODELS`; there is no joint
model), so the two labels are named separately, with the defaults of
`templates.gain_template_from_pairs`.

**It reuses that code, not a second implementation.** `templates.gain_template`
already halves the V² widths onto log |V| and keeps offsets in degrees;
`gain_template_from_pairs` already reads one model's widths and records the
white factors in provenance. `noise_spec` calls `gain_template_from_pairs` (or
shares the width-reading helper) and only adds what the template does not
carry: the white scales as spec fields, the degrees-to-radians conversion of
the offset, and the spec's provenance block. Each conversion gets a test.

| Pair-test term | Spec field | Conversion |
| --- | --- | --- |
| V² white scale (× σ) | `white.v2` | none |
| closure-phase white scale | `white.t3` | none |
| V² gain per telescope, per baseline (`Fit.widths`, on V²: IW Car 1.1 %, 0.9 %) | `gains.telescope`, `gains.baseline` | halve (done in `gain_template`; once only). IW Car: 0.0055, 0.0045 on log \|V\| |
| closure offset per triangle (degrees) | `closure_offsets.triangle` | degrees to radians |
| (not fitted) | `gains.chromatic` | `null` |

**The white scale is relative to the reported σ.** The pair test divides by the
σ in the product, and for calibrated products that σ may include a
transfer-function error (`--tf-error loo/scatter`). Consecutive frames share
calibrators, so such a σ already holds part of the correlated error and the
white scale is then too small to use on data whose σ does not include it. The
provenance therefore records `sigma_convention` (the `tf_error` mode the
products were calibrated with, or `"as_reported"`), and the pipeline side warns
when the data's convention is known and differs.

Provenance filled from the test and its inputs: tool and version, the run
description (number of disjoint pairs, whether a model was subtracted), date,
target, programme ID if given, `night` if a single night, and
`bound: "lower"` (consecutive frames share calibrators). `bandwidth` is a
`read_oifits` option of virgil; `products.read_exposure` has no such argument, so
the emitter cannot know it and leaves the field out unless the caller passes it
as provenance.

The function writes nothing; `json.dump(..., sort_keys=True)` is the caller's.
It does not import virgil's validator; its test validates the output against
`virgil.noise_spec.validate_noise_spec` when virgil provides it (skipped
otherwise), and against a pinned copy of the 0.4.1 field list always.

## CLI

`virgil-vlti-pairs ... --emit-spec spec.json` runs the pair test as now and
writes the spec. `pair_test` passes the blocks of every target to one
`fit_noise`, so a multi-target run yields **one pooled set of widths**, not one
per target; pairs never mix targets, but the fit does. A pooled spec carries
`target: null` and the list of targets in provenance. For a per-target spec
(what the post-AGB work used) run the tool once per target (`--target`, or one
input directory per target) and write one file each.

## Checks on emission

- Refuse widths that look like a unit slip, matching virgil's validator: a
  closure offset above 0.05 rad (0.21 written in degrees would pass a looser
  guard as 12°) and a gain width above 0.1 on log |V|. Warn above 0.02 rad.
- Refuse a V² white scale below 0.05, with the numbers. `fit_noise` optimises
  log-widths without bounds (only |log w| > 30 is penalised), so no width sits
  at a bound and none is exactly zero; the check needs an explicit, documented
  tolerance. The collapse seen in `pionier_systematics.md` (gain terms built
  from noisy V², absorbing the white noise) shows as white → 0 with large
  gains, so the white scale is the quantity to test. Warn if both gain widths
  are below 1e-4 (no gain structure found).
- Warn when `n_pairs` is below 6; the widths are then poorly constrained
  (HR 4049 had 9, EN TrA 7).
- Warn that the closure offset is an upper bound when no model was subtracted
  (IW Car: 0.40° without, 0.21° with).

## Open

- Whether the stored spec should carry the ΔlnL and `n_pairs` as `fit` metadata
  (useful for judging a spec; not interpreted by virgil).
- A per-programme pooled spec, once enough targets are measured.

## Order

1. This PR: the design.
2. `noise_spec` and `--emit-spec` with the conversion tests (before virgil PR 2
   of the sequence in virgil's note (virgil#420), so pipelines have something to read).
