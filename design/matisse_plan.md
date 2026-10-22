# MATISSE support: implementation plan

Status: **plan**, 2026-10-07. **Mostly unimplemented** (a background working note, not user documentation). Stages 0 (product classification) and 1 (header reader, refusals and `virgil-vlti-inspect`) are implemented; Stage 1's deviations from this text are listed in its PR (`matisse-stage1`). Until Stage 4 wires `group_key` into `virgil-vlti-calibrate`, the CLI groups by night only, so `compatibility` also refuses a calibrator of another DIT, and a night with several BCD positions fails compatibility on `bcd_keys` (safely: nothing is calibrated). `calpca` and the pair test skip MATISSE products until they take `group_key` (Stages 4 and 5). This note turns the work items of the [MATISSE calibration review](matisse_calibration_review.md) ("the review", cited by § and footgun number) into stages that can each be implemented as one PR without further design. It respects the boundary in [`PLAN.md`](../PLAN.md): virgil never imports virgil-vlti, virgil-vlti imports virgil only at the hand-off, and the hand-off is OIFITS plus `.npz` nuisance modes. Generic machinery (smearing, primary beams, nuisance modes) belongs in virgil.

Sources (full list in §9): the MATISSE pipeline manual, issue 2.5.0 ([PM]); the MATISSE user manual, issue 1.0 for P118 ([UM]); Lopez et al. 2022, A&A 659, A192 ([Lo22]); Varga et al. 2021, arXiv:2012.05697 ([Va21]); Gámez Rosas et al. 2022, arXiv:2112.13694 ([GR22]); Isbell et al. 2022, A&A 663, A35 ([Is22]); Millour et al. 2016, arXiv:1608.01913 ([Mi16]); the consortium tools ([MT]). The review's §11 says how far each was read. §7 and the facts below were re-checked against [PM], [UM] and [Lo22] in a conservative review on 2026-10-07; what could not be verified is listed in §7.1 for the MATISSE expert. Code references are to `origin/main` at 11b4ada. Archive numbers are from read-only ESO TAP and CalSelector queries of 2026-10-07; no data were downloaded.

## Contents

1. Scope
2. What is GRAVITY-specific in the code today
3. Conventions used by every stage: BCD tables, station keys, groups and visits, names
4. Stages 0–8, each with the same subheadings: goal; interfaces; algorithm; synthetic products and tests; refactors and GRAVITY guards; validation; acceptance; risks and fallbacks; effort and dependencies
5. Dependencies, parallel work and effort
6. The error model: what differs from GRAVITY
7. Decisions (2026-10-07), with the conservative review and §7.1, unverified items
8. virgil issues to open (drafts)
9. References

## 1. Scope

**In scope, in this order:**

- **Bands.** L and M (one detector, `ESO DET NAME = MATISSE-LM`) first; N (`MATISSE-N`) second. The two bands are always separate datasets, even when their `INSNAME`s agree (review §1c.3).
- **Resolutions.** LOW in L/M and N first, since most archival data use it. MED and HIGH (L/M) and HIGH (N) follow with new instrument facts and tests, not new code. VHIGH is accepted but not validated.
- **Telescopes.** ATs and UTs. They differ only in facts (pinhole size, photometric error, widths), never in code paths.
- **Inputs.** Per-exposure `TARGET_RAW_INT` and `CALIB_RAW_INT` products of `mat_raw_estimates`. Calibrated `*_CAL_INT` products stay refused (`products.is_calibrated`).
- **Observables.** V² and closure phases in L/M; then correlated flux, closure phases and (optionally) V² in N; then the two BCD-correctable differential phases.

**Out of scope at first:**

- **GRA4MAT-specific processing.** Data *taken* with GRAVITY as fringe tracker are in scope (ν Hor used it: `ESO TPL ID = MATISSE_hyb_obs_ft`, `ESO DEL FT SENSOR = GRAVITY`), since their `RAW_INT` products are calibrated the same way. Out of scope: a posteriori co-phasing of N from the K-band delays, and the wide-band MED/HIGH modes it enables ([PM] §11.7–11.8).
- **Reducing raw frames** in the package. That stays with the ESO pipeline, run in the container of `containers/` (Stage 2).
- **Total-flux calibration** of `OI_FLUX` with spectrophotometric standards ([GR22] §4.6). Only flux *ratios* are kept, as `calibrate.flux_ratio` does now.
- **Merged products** of the consortium tools (`mat_mergeAllOiFits.py`). They are detected and refused as calibration inputs (review footgun 7), not re-split.
- **Channel covariances** (`OI_CORR`). Binning to one resolution element remains the rule (review §5c.1) until virgil has a banded covariance (§8, issue 4).

## 2. What is GRAVITY-specific in the code today

Read from the code at 11b4ada. "Generic" means the logic needs no change for MATISSE, only its inputs. The last column names the stage that changes it.

| Where | What | GRAVITY-specific? | Stage |
|---|---|---|---|
| `products.read_exposure` | reads `OI_WAVELENGTH`, `OI_VIS2`, `OI_T3`, optional `OI_VIS` per `INSNAME` | generic | — |
| `products.default_insnames` | a single `INSNAME` is used directly; else `gravity.science_insnames` | generic for MATISSE (one `INSNAME` per file) | — |
| `products.Exposure.pipeline_diameter` | `ESO QC TF CAL_SC DIAM` from `*_CAL_VIS` | yes | 4 |
| `products.Exposure.night` | noon-UT boundary | generic | — |
| `products.classify` | dual field from `ESO INS SOBJ X/Y`; MATISSE categories handled (Stage 0) | dual field only for GRAVITY | — |
| `calibrate.check_compatible` | calls `gravity.compatibility` unconditionally | yes: `MODE_KEYS`, `REDUCTION_PARAMS = ("vis-correction-sc",)`, `WARN_KEYS` with `ESO DET2 SEQ1 DIT` | 1 |
| `gravity.pipeline_params` | reads `ESO PRO REC1 PARAMn` only | yes: MATISSE products carry several recipes | 1 |
| `calibrate.calibrator_estimates` | keyed on the raw `STA_INDEX` tuple | generic, but wrong when a BCD position reverses a baseline (review §8b.2) | 3 |
| `calibrate.predict`, `prediction_error`, `prediction_residuals`, `binned_estimates`, `_apply` | interpolation, LOO error, binned TF error (#21) | generic; LOO is per file | 4 |
| `calibrate.channel_correlation` | rank-one split of channel residuals | generic | 5 |
| `calibrate.transfer_function_modes` | names `tf-vis2-base<a>-<b>@<night>`, `tf-t3phi-tri<a>-<b>-<c>@<night>` | generic, but one group per night | 4 |
| `calibrate.diameter_modes`, `frame_times`, `group_by_night` | finite-difference diameter modes; one MJD per frame | generic | — |
| `calibrate` module docstring, `nbin` advice | "twofold oversampled" | GRAVITY wording; MATISSE LOW is about fivefold (review footgun 4) | 3 |
| `cli.py` | `--dual-field`, `--no-pipeline-diameters`, `--onaxis-mas`, `--aniso-warn`, the `gravity_viscal` parity group, `--insname` help | GRAVITY options | 4 |
| `pairs.frames_of`, `find_pairs` | frames split by 60 s gaps; pairs by time and grid | generic, but pairs across BCD positions are not like pairs | 5 |
| `calpca._group_key` | night, `viscal.setup_string`, stations, `gravity.field_mode` (#25) | yes: `setup_string` reads GRAVITY keys | 5 |
| `calpca.calibrator_residuals` | requires `pipeline_diameter`; `insname="GRAVITY_SC"`; `references` (#26) | yes | 5 |
| `templates.Template.mismatches` | resolution, polarization, telescopes, field (GRAVITY keys) | yes | 5 |
| `templates.template_name` | `instrument="gravity"` default | generic | 5 |
| `viscal.setup_string`, `calibrate_like_viscal` | GRAVITY keys; `gravity_viscal` rules | yes, and stays GRAVITY-only | — |
| `synthetic.write_product` | GRAVITY layout, FT table, GRAVITY header cards | yes, and stays so; MATISSE gets its own writer | 3 |
| `write.calibrated_hdus`, `write_calibrated`, `write_modes` | OIFITS and `.npz` writing | generic | — |

The refactor is small and local: one instrument dispatch for headers (compatibility, setup, pipeline parameters), one orientation normalization in the calibrator estimates, and an observable-aware group key that includes BCD, and chopping for V² only.

## 3. Conventions used by every stage

### 3.1 BCD positions, inputs and peaks

Each telescope k feeds an interferometric input `ESO ISS CONF INPUTk` (1, 3, 5 or 7; ν Hor: K0→1, G2→3, D0→5, J3→7). A detector fringe peak is labelled by the pair of inputs it combines. The BCDs swap inputs before the beam combiner: one swaps IP1↔IP3, the other IP5↔IP7 ([PM] §4.1; [UM] §2.2.1 gives the same swap as beam order 1-2-3-4 → 2-1-4-3). [Lo22] §2.4–2.5 describes it in beam numbers (1↔2 and 3↔4) and uses no IP labels. *Corrected 2026-10-07: this sentence cited [PM] §11.9 and [Lo22] §2.3; §2.3 of [Lo22] is about signal encoding.* Neither manual says which device is "BCD1", nor gives the header keys. This note writes a position as `<IP1/IP3 BCD>-<IP5/IP7 BCD>`, so `IN-OUT` means the IP1↔IP3 BCD is in. That order is *inferred*: the [PM] §11.9 closure-phase formula for triangle IP1-IP3-IP5, (φ_OUT-OUT + φ_OUT-IN − φ_IN-OUT − φ_IN-IN)/4, negates exactly the positions in which the first word is `IN`, and only the IP1↔IP3 swap reverses that triangle. **The consortium tools contradict this order** if read with the row order below: `mat_removeBCD` ([MT] `libPostTools.py` at b44f965, the `BCDsign`, `BCD`, `BCDcp` and `BCDcpsign` tables) indexes positions by `2·(BCD1 = IN) + (BCD2 = IN)` from `ESO INS BCD1/BCD2 NAME`, and for BCD1 in, BCD2 out it reverses row 0, (5,7), swaps rows (3,5)↔(3,7) and (1,5)↔(1,7), and swaps triangles (1,3,5)↔(1,3,7). That is the IP5↔IP7 swap. So either `ESO INS BCD1` is the IP5↔IP7 device and the manual's position names are not in BCD1-BCD2 order, or the tools' row order differs (a comment in those tables asks why one sign is what it is). **Therefore the tables of this section are labelled by device, not by header key:** the first word of a position is the state of the IP1↔IP3 device. Which header key is which device is the Stage 2 fact `bcd_key_order` (`"BCD1=IP13"` or `"BCD1=IP57"`), decided on a raw product by the [PM] Table 4.1 shutter sequence and the row moves, and until it is recorded every MATISSE product is refused for BCD combination and for calibration across positions. *Revised 2026-10-07 (conservative review): this paragraph took `IN-OUT` to mean "BCD1 in", with only `ESO CFG BCD MODE` to verify.* Also to verify: that `ESO CFG BCD MODE` uses the same order as the INS keys. The ν Hor products give `OUT-OUT` for both, which does not decide it. The `IN-IN` photometric shutter sequence of [PM] Table 4.1 (IP3, IP1, IP7, IP5 against IP1, IP3, IP5, IP7 at `OUT-OUT`) is a cross-check on a raw product.

The `OI_VIS2` row order of an `OUT-OUT` product, in ν Hor (`STA_INDEX` K0=28, G2=19, D0=13, J3=24), is (13,24), (28,19), (19,13), (19,24), (28,13), (28,24). In inputs that is (5,7), (1,3), (3,5), (3,7), (1,5), (1,7). The baseline seen at each peak, by position, labelled by the `OUT-OUT` input of each telescope, is:

| Peak | OUT-OUT | IN-OUT (BCD1 in) | OUT-IN (BCD2 in) | IN-IN |
|---|---|---|---|---|
| (5,7) | (5,7) | (5,7) | **(7,5)** | **(7,5)** |
| (1,3) | (1,3) | **(3,1)** | (1,3) | **(3,1)** |
| (3,5) | (3,5) | (1,5) | (3,7) | (1,7) |
| (3,7) | (3,7) | (1,7) | (3,5) | (1,5) |
| (1,5) | (1,5) | (3,5) | (1,7) | (3,7) |
| (1,7) | (1,7) | (3,7) | (1,5) | (3,5) |

Bold marks a reversed baseline. Only (1,3) and (5,7) stay on their peak; the other four move between peaks. These two are the BCD-correctable baselines of [PM] §11.9 for differential phases; the other four are only averaged over positions there. [PM] §11.9 writes one sign pattern, (OO + OI − IO − II)/4, for both correctable differential phases. By the swap above, (5,7) should instead reverse under the IP5↔IP7 device (OO +, IO +, OI −, II −); this is listed in §7.1.

The `OI_T3` rows are, in inputs, (3,5,7), (1,3,5), (1,3,7), (1,5,7). Every triangle contains exactly one of (1,3) and (5,7). The physical triangle measured at each peak, with the sign s of the permutation that brings it to the peak's order:

| Peak | OUT-OUT | IN-OUT | OUT-IN | IN-IN |
|---|---|---|---|---|
| (3,5,7) | (3,5,7) + | (1,5,7) + | (3,5,7) − | (1,5,7) − |
| (1,3,5) | (1,3,5) + | (1,3,5) − | (1,3,7) + | (1,3,7) − |
| (1,3,7) | (1,3,7) + | (1,3,7) − | (1,3,5) + | (1,3,5) − |
| (1,5,7) | (1,5,7) + | (3,5,7) + | (1,5,7) − | (3,5,7) − |

Triangles (1,3,5) and (1,3,7) flip sign under the IP1↔IP3 device; (3,5,7) and (1,5,7) flip under the IP5↔IP7 device. Position codes are in device order (first letter: the IP1↔IP3 device). For a triangle that flips under the IP1↔IP3 device the signs are OO +, OI +, IO −, II −, which is the manual's (φ_OO + φ_OI − φ_IO − φ_II)/4 ([PM] §11.9, given there for the first triangle, IP1-IP3-IP5, only; the other three are not written out). For a triangle that flips under the IP5↔IP7 device they are OO +, IO +, OI −, II −.

**Why the BCD combination needs equal weights.** After the sign is undone, the phase of physical triangle τ at position p is φ_p = φ_true(τ) + s_p φ_inst(peak_p), where φ_inst is what calibration left of the instrumental closure phase. For τ = (1,3,5): OO gives +φ_i(135), IO −φ_i(135), OI +φ_i(137), II −φ_i(137). Cancelling both residuals needs only that each triangle's weights obey w_OO = w_IO and w_OI = w_II for a triangle flipping under the IP1↔IP3 device, and w_OO = w_OI and w_IO = w_II for one flipping under the IP5↔IP7 device. The **equal-weight** mean over all four positions is the smallest rule that satisfies all four triangles at once. Inverse-variance weights break the cancellation (they differ between positions, as in the 5/1/1/5 exposure balance of the LM template), and OO plus II alone leaves (φ_i(135) − φ_i(137))/2. So exposures are averaged within a position first (any weights), then across positions with weight 1/4 each.

**Whether the pipeline relabels.** If `mat_raw_estimates` writes `STA_INDEX` in physical order (the telescopes actually combined, in their actual orientation), the station normalization of §3.2 does all the bookkeeping. If it writes the peak's `OUT-OUT` labels whatever the position, virgil-vlti must apply the tables above. Stage 2's explore job decides this from one template: in an `IN-IN` exposure, the row on peak (1,3) has `STA_INDEX` (19,28) if physical, (28,19) if not. The result is recorded as the instrument fact `sta_index_frame` (`"physical"` or `"peak"`), and `matisse.relabel_bcd` (Stage 4) applies the tables only for `"peak"`.

### 3.2 Orientation-normalized station keys

- **Key.** The sorted `STA_INDEX` tuple. Its parity π = ±1 is the sign of the permutation from the file's order to the sorted one.
- **Baseline reversed** (π = −1): VISPHI → −VISPHI, (u, v) → (−u, −v); V², VISAMP and fluxes unchanged.
- **Triangle permuted:** T3PHI → π · T3PHI. The legs are recomputed from relative station positions: with p_a = 0, p_b = uv₁, p_c = uv₁ + uv₂ in the file's order (a, b, c), the new order (a′, b′, c′) has uv₁′ = p_b′ − p_a′ and uv₂′ = p_c′ − p_b′. T3AMP is unchanged.
- **When to flip.** A calibrator row is re-oriented only when its orientation differs from the science file's. GRAVITY and PIONIER files never differ, so their outputs stay byte-identical.
- **Names.** Mode names and written rows use the instrument fact `canonical_orientation`: `"file"` for GRAVITY and PIONIER (today's names) and `"sorted"` for MATISSE. Combined MATISSE files are therefore written in sorted order, so virgil's triangle legs (a,b), (b,c), (a,c) match.

### 3.3 Groups, visits and cycles

- **Group** (what may calibrate what): `matisse.group_key(exposure, observable)` = (night, band, resolution, DIT, BCD position, chopped), stricter than `mat_autoCalib` ([MT] `mat_autoCalib.py` at b44f965, `findClosestCal`), which requires equal `ESO INS BCD1 NAME`, `BCD2 NAME`, `ESO DET CHIP NAME` (the band), `ESO DET SEQ1 DIT` and `ESO ISS CHOP ST`, for every observable, and takes calibrators within 3 h rather than by night; it does not compare resolution. *Corrected 2026-10-07: this said the key is "as `mat_autoCalib` groups".* `chopped` is part of the key for every observable by default. *Revised 2026-10-07 (conservative review): the key was observable-aware, with `chopped` dropped for `"t3phi"`. [PM] §4.1 says chopping is unnecessary to measure closure phases, which is about the estimator, not a statement that chopped and non-chopped exposures share an instrumental closure phase; chopped L exposures exist only at `OUT-OUT` and `IN-IN` ([PM] Table 4.1), so pooling them for closure phases also changes the balance of positions. Splitting costs little, since calibrator templates follow the same sequence. A `--pool-chopping-t3` option restores the pooled key once D1's calibrators show no closure-phase difference between chopped and non-chopped exposures of one position (Stage 4 validation), and §7.1 asks the MATISSE expert.* GRAVITY keeps (night, setup). DIT is a hard split for both observables.
- **Visit**: (target, `ESO TPL START`). All exposures of one calibrator template share one transfer-function sample in time, so they are not independent for a leave-one-out error.
- **Cycle**: the exposures of one science template; the unit the BCD combination works on.

### 3.4 File and mode naming

- Per-position outputs: `<out>/<night>/<slug>_<stamp>_<bcd>[_chop].fits` (e.g. `_IN-IN_chop`), with the `.modes.npz`, `.tf_modes.npz` and `.template_modes.npz` beside them, as for GRAVITY.
- Combined output: `<out>/<night>/<slug>_<tplstart>_bcd.fits` and its `.npz` files. One frame per cycle, at the mean MJD of the cycle's exposures. The scatter across positions goes to a sidecar `<...>_bcd.scatter.npz`, never into a modes file.
- Mode names follow GRAVITY's, with the position code (OO, IO, OI, II) before the `@`: `tf-vis2-base13-24-OO@2023-12-11`, `tf-t3phi-tri13-19-24-II@2023-12-11`. Names shared across positions (§7, decision 10) drop the code. New MATISSE names: `spec-<calibrator>@<night>-k<order>` (Stage 7) and `phot-tel<sta>@<file>` (Stage 5).
- The `.npz` layout is unchanged: `write_modes(path, {insname: {mode: {observable: array}}})`, and `stack_modes` keys stay `"<insname>/<observable>/<mode>"`, aligned by file, with `u` for vis2 and `phi` for t3phi (`docs/templates.md`).
- `index.csv` keeps its eight GRAVITY columns, byte for byte. MATISSE runs append `band`, `resolution`, `bcd`, `chopped`, `tf_error_method`, `n_cal_visits` and `status`; GRAVITY runs never write them. `calibrated` keeps its meaning, the output path relative to `--out` (`cli.py` writes `os.path.relpath(path, args.out)`), and is empty when the exposure was not calibrated; `status` is `ok` or `no-calibrator`. MATISSE runs also write a row for a science exposure with no calibrator in its group, which GRAVITY runs do not.

## 4. Stages

Each stage is one PR (Stage 3 is two) and useful on its own. Every stage that adds a method also adds its rows to [`PROVENANCE.md`](../PROVENANCE.md): the source, how it is checked, and whether an independent check is possible (AGENTS.md). Sizes are rough lines of new or changed code plus tests: S < 300, M 300–800, L > 800.

**GRAVITY guards, used by every stage from 1 on.**

- **G1, the existing suite.** `tests/test_calibrate.py` (including `test_calibrators_in_another_mode_are_always_refused`, `test_other_pipeline_settings_refused_unless_not_strict`, `test_tf_modes_are_named_by_night_and_baseline`, `test_missing_stations_are_flagged`, `test_cli_writes_tf_modes_and_drops_the_diagonal_term`, `test_point_source_calibration_is_exact`, `test_split_polarisation_calibrated_separately`), `test_viscal.py` (`test_cli_parity_mode`, `test_setup_mismatch_dropped_unless_forced`), `test_dual_field.py`, `test_templates.py` (`test_mismatch_refused_and_range_warnings`, `test_shipped_pilot_template`), `test_pairs.py` and `test_calpca.py`. All pass unchanged; no assertion is edited.
- **G2, a golden synthetic night.** Stage 1 adds `tests/test_gravity_golden.py::test_gravity_outputs_unchanged`. It runs `virgil-vlti-calibrate` on a synthetic GRAVITY night (two calibrators, one science file, `--tf-modes`, `--tf-error loo`), and compares every array of the `.fits` and `.npz` outputs, and `index.csv`, with `tests/golden/gravity_golden.npz` by `np.allclose(rtol=1e-9, atol=1e-12)`. The file is under `tests/golden/`, not `tests/data/`, because `.gitignore` has `data/`, which ignores a `tests/data/` directory (`git check-ignore tests/golden/gravity_golden.npz` reports it as not ignored). The golden file (a few kB) is generated from 11b4ada in Stage 1 by the named script `tests/golden/make_gravity_golden.py`, and is regenerated only deliberately, by running that script, with the reason in the PR. Exact equality is not used, since interpolation and `j1` may differ in the last bits across numpy and scipy versions and platforms. The primary check is the before/after comparison within one environment: each refactor PR runs the script at 11b4ada and at its own head in the same environment and compares the two outputs (exactly, `rtol=0`); the golden file is the portable check for CI.
- **G3, real data, by hand in each refactoring PR (1, 3b, 4, 5).** a private target-specific script on real GRAVITY data before and after; outputs compared with `astropy.io.fits.FITSDiff(ignore_keywords=["DATE"])` and `np.array_equal` on every `.npz` array. The `gravity_viscal` parity numbers (V² 0.12%, closure phases 0.06°; private validation suite) are re-run and must agree to every printed digit. This needs numpy only, so it runs on the laptop.

### Stage 0 (done). Product categories

`products.classify` recognizes `TARGET_RAW_INT` and `CALIB_RAW_INT`, and refuses `*_CAL_INT` (review footgun 13).

### Stage 1. Header reader, refusals and `virgil-vlti-inspect`

#### Goal

Read every MATISSE fact the calibration needs from a header, and refuse what must never reach virgil. Useful alone: a user can ask what a MATISSE file is and whether virgil may read it.

#### Interfaces

New `src/virgil_vlti/matisse.py`, the analogue of `gravity.py`:

```python
BANDS = {"MATISSE-LM": "LM", "MATISSE-N": "N"}
POSITIONS = ("OUT-OUT", "IN-OUT", "OUT-IN", "IN-IN")
CODES = {"OUT-OUT": "OO", "IN-OUT": "IO", "OUT-IN": "OI", "IN-IN": "II"}

@dataclass(frozen=True)
class Setup:
    band: str               # "LM" or "N", from ESO DET NAME
    resolution: str         # ESO INS DIL NAME (LM) or ESO INS DIN NAME (N)
    bcd: str                # position in device order, "<IP1-IP3 device>-<IP5-IP7 device>", reordered from
                            # the INS keys by the bcd_key_order fact (§3.1); setup() refuses if the fact is absent
    bcd_keys: tuple[str, str]   # raw (ESO INS BCD1 NAME, ESO INS BCD2 NAME), never used for signs
    chopped: bool           # ESO ISS CHOP ST == "T"; never ESO DET CHOP ST (keys from the ν Hor headers; neither manual names them, §7.1)
    dit: float | None       # ESO DET SEQ1 DIT [s]
    telescopes: str | None  # "AT" or "UT", from ESO ISS CONF T1NAME
    spatial_filter: str | None   # pinhole or slit name, from the INS keys found in Stage 2 (§7.1); never assumed
    fringe_tracker: str | None   # ESO DEL FT SENSOR
    template: str | None    # ESO TPL ID
    tpl_start: str | None   # ESO TPL START
    inputs: tuple[int, ...] # ESO ISS CONF INPUT1..4

def setup(header) -> Setup
def position_code(setup: Setup) -> str                 # "OO", "IO", "OI", "II"

@dataclass(frozen=True)
class Recipe:
    index: int                    # i in ESO PRO RECi
    name: str                     # ESO PRO RECi ID
    pipe_id: str | None           # ESO PRO RECi PIPE ID, e.g. "matisse/2.0.0"
    version: tuple[int, int, int] | None
    params: dict[str, str]        # PARAMn NAME -> PARAMn VALUE

def parse_version(pipe_id: str | None) -> tuple[int, int, int] | None
def recipes(header) -> dict[str, Recipe]               # by name; lowest i wins
def raw_estimates(header) -> Recipe | None             # the mat_raw_estimates entry
def merge_evidence(header, tables=None) -> list[str]   # empty if none
def is_merged(header, tables=None) -> bool
def phase_convention(header, wave=None) -> Literal["standard", "flipped", "unknown"]
MODE_KEYS: tuple[str, ...]          # band, resolution, BCD, telescopes, spatial filter, chopping
REDUCTION_PARAMS = ("corrFlux", "coherentAlgo", "cpAlgo", "spectralBinning", "spectralAverage", "replaceTel")
# 2.5.0 renamed spectralBinning -> spectralAverage and coherentAlgo -> cpAlgo ([PM] §2.1), and
# cpAlgo=2 means coherent closure phases, not the bispectrum. A product with the old names and one
# with the new are a reduction mismatch (refused unless not strict), never silently equated.
WARN_KEYS = ("ESO PRO REC1 PIPE ID",)   # the version of mat_raw_estimates, as gravity.WARN_KEYS reads it
def compatibility(sci_header, cal_header, observable="vis2") -> tuple[list[str], list[str], list[str]]
def group_key(exposure, observable="vis2") -> tuple
```

In `products.py`:

- `instrument_of(header) -> str`: stripped, upper-case `INSTRUME`.
- `INSTRUMENT_MODULES = {"GRAVITY": gravity, "MATISSE": matisse}`. Any other instrument (PIONIER today) falls back to `gravity`, so its behaviour is unchanged.
- `pipeline_params(header, recipe=None)`, generic. `gravity.pipeline_params(header)` becomes `products.pipeline_params(header, recipe=1)`.

New CLI `virgil-vlti-inspect PATH... [--json]` (`virgil_vlti.inspect_cli:main`, an entry point in `pyproject.toml`). It prints one row per file: instrument, category, band, resolution, BCD, chopped, pipeline version, phase convention, merged (with evidence), calibrated, and a verdict (`ok`, or the refusal). Exit status 0 if every file is `ok`, else 1. It also prints `replaceTel` and warns when it is not 0: [PM] §9.1.4 defines it as a two-digit code (UT digit, AT digit), each digit k ≠ 0 replacing telescope k's photometry by the mean of the other three, and its tabulated default is 3, which replaces AT3's photometry.

`virgil-vlti-calibrate` gains `--n-phase {refuse,negate,keep}` (default `refuse`), which applies only to N-band products whose convention is `flipped`, and `--n-phase-unknown {refuse,keep,negate}` (default `refuse`), which applies only to `unknown` ones (no version, or before 1.5.2). The two are separate so that negating flipped data never rotates older, probably standard, data by 180° without being asked. `negate` multiplies T3PHI and VISPHI by −1 on read and writes `HISTORY virgil-vlti: N-band phases negated (mat_raw_estimates <version>)`; `keep` trusts the file.

#### Algorithm

1. **`setup`.** Read the keys above; strip and upper-case strings. `ESO ISS CHOP ST` is `'T'`/`'F'`. In the ν Hor `_Chop_` products `ESO DET CHOP ST` is `False` although the ISS key is `T`, so the detector key is never read. If `ESO CFG BCD MODE` exists and disagrees with `BCD1-BCD2`, the exposure is refused, naming both values. *Revised 2026-10-07 (conservative review): this was a warning with the INS keys winning. Which key gives which device is not documented (§3.1, §7.1), so a disagreement is an unresolved convention.*
2. **`parse_version`.** Match `r"(\d+)\.(\d+)\.(\d+)"` in the `PIPE ID` (`matisse/1.7.9` → (1, 7, 9)); anything else is `None`.
3. **`phase_convention`.** N band is `setup.band == "N"`, or, without `ESO DET NAME`, a median wavelength above 7 µm. In N, a `mat_raw_estimates` version in [(1,5,2), (2,0,0)) gives `"flipped"`, a version at or after (2,0,0) gives `"standard"`, and anything else (`None`, or before (1,5,2)) gives `"unknown"`. L and M are always `"standard"` (review footgun 3; [PM] §2.5; [GR22] §4.8; [ESO-news], December 2023). *Revised 2026-10-07 (conservative review): versions before 1.5.2 were `"standard"`. [PM] §2.5 (the 2.0.0 release notes) says the N-band phases and closure phases were inconsistent with L and M "as of" 1.5.2, giving a 180° rotation, and corrected in 2.0.0; but the 1.5.2 release notes ([PM] §2.11) do not mention a sign change, and the ESO news item of December 2023 says only that 2.0.0 fixed it. The lower bound is therefore not established, and pre-1.5.2 N-band phases are refused by default like the flipped ones (§7, decision 2; §7.1).* [PM] §2.5 names phases and closure phases; that `VISPHI` flips with `T3PHI` follows from the rotation and is listed in §7.1. Negating T3PHI and VISPHI is the same as (u, v) → (−u, −v): a rotation of the image by 180°.
4. **`merge_evidence`** returns every reason it finds; any one marks a merge:
    - the table `INT_TIME` exceeds 1.5 × `ESO DET SEQ1 EXPTIME` (ν Hor: 120.46 s against 61.8 s);
    - one baseline appears in `OI_VIS2` at more than one MJD. This is a warning only until Stage 2's explore job has shown how many rows and MJDs per baseline a single `RAW_INT` has (if `mat_raw_estimates` writes one frame per chop cycle or block, every product would be refused); it then becomes fatal if the check allows;
    - `PIPEFILE` or `ARCFILE` matches `_cal_oifits_`, `_Chop`, `_noChop` or `_CAL_INT`, or a `HISTORY` card names `mat_mergeOifits`.

    A count of `ESO PRO REC1 RAWn` cannot detect a merge: a single ν Hor product lists all 12 template exposures and 2 skies.
5. **`compatibility`.** Mode keys: band, resolution, BCD, telescopes, spatial filter, and chopped, for every observable; with `--pool-chopping-t3`, a chopping difference is not a mismatch for `"t3phi"` (§3.3, revised). DIT (`ESO DET SEQ1 DIT`) is a hard split in `group_key`, so it never reaches `compatibility` and is not a warning key. Reduction parameters are compared as stripped strings, from the `mat_raw_estimates` recipe. Warning: the pipeline version. It is read as GRAVITY's is, from the header key `ESO PRO REC1 PIPE ID` (`gravity.WARN_KEYS`); because MATISSE products can carry several recipes, `Recipe.pipe_id` of the `mat_raw_estimates` entry (`ESO PRO RECi PIPE ID`, found with `raw_estimates(header)`) is used when it is not REC1. In N, a different phase convention between science and calibrator is a mode mismatch, since calibration subtracts phases.
6. **Dispatch.** `calibrate.check_compatible` calls `products.INSTRUMENT_MODULES.get(instrument_of(h), gravity).compatibility`.
7. **Refusals** in the callers of `read_exposure` (`cli.py`, `calpca`, `pairs`):
    - `corrFlux=true` products are refused for V², T3AMP and VISAMP as a visibility, naming the keyword (review footgun 1); only Stage 7 reads them;
    - merged products are refused as calibration inputs;
    - N-band phases follow `--n-phase`.

#### Synthetic products and tests

Headers only (full products come in Stage 3a). New `tests/test_matisse_headers.py`. The ν Hor header is a dict literal of the cards in this note; no data file is committed.

| Test | Asserts |
|---|---|
| `test_setup_reads_lm_low_from_dil` | band LM; resolution from `DIL NAME` when `DIL` and `DIN` differ |
| `test_setup_reads_n_resolution_from_din` | band N uses `DIN NAME` |
| `test_chopping_read_from_iss_not_det` | `ISS CHOP ST = 'T'`, `DET CHOP ST = False` gives `chopped=True` |
| `test_bcd_position_from_ins_keys` | with `bcd_key_order="BCD1=IP13"`, `BCD1=IN`, `BCD2=OUT` gives `IN-OUT`, code `IO` |
| `test_cfg_bcd_mode_disagreement_refused` | refused, naming both values |
| `test_bcd_position_in_device_order` | with `bcd_key_order="BCD1=IP57"`, keys BCD1=IN, BCD2=OUT give `bcd="OUT-IN"` and `bcd_keys=("IN","OUT")`; with `"BCD1=IP13"`, `"IN-OUT"` |
| `test_recipes_parsed_per_recipe` | REC1 `mat_raw_estimates` and REC2 `mat_cal_oifits`, each with its own params |
| `test_parse_version` | `matisse/1.7.9` → (1,7,9); `""`, `None`, `"matisse"` → `None` |
| `test_phase_convention_by_version` | over 1.5.1, 1.5.2, 1.7.9, 2.0.0, 2.5.0 × LM, N: N 1.5.2 and 1.7.9 are flipped, N 1.5.1 is unknown, all LM standard |
| `test_unknown_version_in_n_band` | no `PIPE ID`: `unknown` in N, `standard` in LM |
| `test_negate_does_not_touch_unknown` | a 1.5.1 N product with `--n-phase negate` is still refused; with `--n-phase-unknown keep` it is read unchanged; only `--n-phase-unknown negate` negates it |
| `test_corrflux_products_refused_for_vis2` | `corrFlux=true` raises, naming `corrFlux` |
| `test_merge_from_int_time` | the ν Hor values are a merge |
| `test_merge_from_repeated_baselines` | two MJDs on one baseline are a merge |
| `test_merge_from_consortium_names` | `_noChop` in `PIPEFILE` is a merge |
| `test_single_exposure_not_merged` | `INT_TIME` = EXPTIME, one MJD per baseline: no evidence |
| `test_compatibility_refuses_other_bcd`, `_other_resolution`, `_other_band`, `_other_spectral_binning` | each lands in the mode or reduction list |
| `test_compatibility_refuses_other_chopping` | a chopped science file and a non-chopped calibrator: a mode mismatch for `"vis2"` and `"t3phi"`; for `"t3phi"` none with `pool_chopping_t3=True` |
| `test_group_key_dit_and_chopping_split_groups` | a different DIT or `chopped` gives a different key for both observables; `chopped` drops from the `"t3phi"` key only with `pool_chopping_t3=True` |
| `test_bcd_key_order_required` | without a recorded `bcd_key_order`, calibration across positions and `combine_bcd` refuse, naming the fact |
| `test_compatibility_warns_on_pipeline_version` | `ESO PRO REC1 PIPE ID` differing: in the warning list only |
| `test_n_band_mixed_phase_convention_refused` | 1.7.9 science, 2.0.0 calibrator: refused |
| `test_gravity_dispatch_unchanged` | a GRAVITY polarization mismatch still raises |
| `test_pionier_falls_back_to_gravity_rules` | PIONIER headers give today's result |
| `test_inspect_cli_verdicts` | ν Hor-like header: merged, calibrated, exit status 1 |

Plus G2's golden file and test.

#### Refactors and GRAVITY guards

`gravity.pipeline_params` becomes a thin wrapper of `products.pipeline_params(header, recipe=1)`, and `check_compatible` dispatches. Guards G1, G2 and G3.

#### Validation

`virgil-vlti-inspect` on the three ν Hor files (local only).

#### Acceptance

The three ν Hor files are reported as LM, LOW, `OUT-OUT`, chopped, `matisse/2.0.0`, merged (`INT_TIME`), calibrated, and refused. G1–G3 pass.

#### Risks and fallbacks

- Some archival products may lack `ESO PRO RECi PIPE ID` (very old reductions). They are `unknown` in N and need `--n-phase`.
- A header key renamed between pipeline versions is caught by Stage 2's header survey, which adds a test per version.

#### Effort and dependencies

S–M. Depends on nothing.

### Stage 2. Validation data and pipeline reductions

#### Goal

Per-exposure `RAW_INT` products for real nights, and the ESO pipeline's own calibration of them, so that Stages 4, 5, 7 and 8 have something to be checked against. No package code. It runs while Stages 1 and 3 are written. Plan only until the maintainer approves the downloads and the container build; decision 1 (§7) fixes the plan, not the go-ahead.

#### Interfaces

Two directories of cluster scripts (private, not in this repository):

- `matisse_dl/`: `build_manifests.py` (runs on the Mac), `download.py`, `download.sbatch`, `submit.sh`, `README.md`, `VERSION`, `manifests/<dataset>.csv`.
- `matisse_reduce/`: `explore.sbatch`, `reduce.sbatch`, `make_sof.py`, `check_products.py`, `submit.sh`, `README.md`, `VERSION`.

Data root: `<cluster data root>/matisse/`, with `data/<dataset>/{raw,master,eso_products}/`, `products/<dataset>/<pipeline version>/{raw_estimates,cal_oifits_block,cal_oifits_bcd}/`, `logs/` and `manifests/`.

#### Datasets

Sizes are the archive's `access_estsize`, summed.

| Id | Dataset | Programme | Files | Size | Use |
|---|---|---|---|---|---|
| D1 | ν Hor L/M, nights 2023-12-11 and 12, ATs K0 G2 D0 J3, GRA4MAT | 112.26WM.004, .005 | 208 | 5.34 GB | **Primary.** Stage 4 parity, Stage 5 PCA |
| D1N | ν Hor N, same templates | same | 112 | 45.9 GB | Stage 7 on the same nights |
| D2 | HD 163296, 2019-03-22, ATs A1234 | 0103.D-0294(A) | 88 | 13.7 GB (both bands) | Stage 7 against [Va21]; a pre-2.0 N-band phase case |
| D3 | HD 163296, 2022-03-21, UTs | 108.225V.011 | 120 | 19.6 GB (both bands) | UT case; calibrators only before the science (extrapolation) |
| D4 | NGC 1068, 2018-09-21 to 24, UTs | 60.A-9257(D) | 92 science; calibrators to count | 18.9 GB science | N-band phase convention ([GR22]); lower priority |

- **D1.** Science: 36 LM `OBJECT` (0.70 GB), 36 `OBJECT,RMNREC` (1.19 GB), 16 LM `SKY` (0.28 GB; 8 templates × 2, so the calibrators' skies are included). Calibrators: 60 LM `STD` (1.18 GB), 60 `STD,RMNREC` (1.99 GB). D1 plus D1N is 320 files and 51.3 GB. Three science templates (2023-12-12T00:47:05, 02:04:14, 2023-12-13T01:09:31) and five calibrator templates: HD 12477 (2023-12-12T00:17:55), HD 16419 (01:13:13), HD 14682 (01:38:26, and 2023-12-13T00:44:18) and κ Hyi (01:35:34). Each template has 12 LM, 12 N, 12 RMNREC and 2+2 SKY files. LM files are 17–20 MB; N files about 226 MB.
- **D2.** The science template 2019-03-23T08:41:19 and the three bracketing calibrator templates 08:15:36 (ζ Ara), 08:53:41 (δ Sgr) and 09:19:37 (ε Sco). In 2019 the `object` column reads `STD` and `SKY`, so names come from the headers.
- **D3.** The science template 2022-03-22T09:40:29 (24 `OBJECT`, 12 RMNREC, 4 SKY; 6.47 GB) and the calibrator templates 07:38:11 (HD 147084) and 08:39:05 (HD 151051), 6.5 GB each. No calibrator follows the science.

#### Queries

ESO TAP at `https://archive.eso.org/tap_obs`, synchronous, run by `build_manifests.py`:

```sql
-- D1 and D1N (add  AND origfile NOT LIKE '%HISENS_N%'  for D1 alone)
SELECT dp_id, origfile, dp_cat, dp_type, object, tpl_start, tpl_expno, tpl_nexp,
       exp_start, access_estsize, release_date, prog_id
FROM dbo.raw
WHERE instrument = 'MATISSE'
  AND prog_id IN ('112.26WM.004', '112.26WM.005')
  AND exp_start BETWEEN '2023-12-11T12:00:00' AND '2023-12-13T12:00:00'
ORDER BY exp_start

-- D2 (same columns)
SELECT ... FROM dbo.raw
WHERE instrument = 'MATISSE'
  AND tpl_start IN ('2019-03-23T08:15:36', '2019-03-23T08:41:19',
                    '2019-03-23T08:53:41', '2019-03-23T09:19:37')

-- D3 (same columns)
SELECT ... FROM dbo.raw
WHERE instrument = 'MATISSE'
  AND tpl_start IN ('2022-03-22T07:38:11', '2022-03-22T08:39:05',
                    '2022-03-22T09:40:29')
```

The counting check is a `GROUP BY dp_type` of the same `WHERE` clause, with `COUNT(*)` and `SUM(access_estsize)`. ESO TAP has no `CASE`, so band splits are made on `origfile` in Python. `build_manifests.py` stops if a count differs from the numbers above.

#### Calibration files

For each science and calibrator frame, `build_manifests.py` asks CalSelector, `https://archive.eso.org/calselector/v1/associations?dp_id=<dp_id>&mode=Raw2Master&responseformat=votable`, and parses the reply with `astropy.io.votable.parse_single_table`.

- For an LM science frame of D1 it returns ESO master `BADPIX`, `NONLINEARITY`, `OBS_FLATFIELD`, `SHIFT_MAP` and `KAPPA_MATRIX` (`M.MATISSE.*`; about 0.13 GB per set, mostly `NONLINEARITY` at 118 MB), the template's `SKY_RAW` and `TARGET_RAW_RMNREC`, and 12 ESO-processed `CALIB_RAW_INT` products (219 kB each). An N frame gets 8; a calibrator frame also gets `JSDC_CAT`.
- Masters are de-duplicated by `dp_id`.
- ESO's `CALIB_RAW_INT` products go to `eso_products/`. They give a free parity check on our container's reduction, and an input for the Stage 5 PCA.
- **Fallback** if a master is missing or the recipe rejects it: `mode=Raw2Raw`. It returns the raw `DISTOR_*`, `KAPPA_*`, `SPECTRA_*`, `OBSDARK` and `FLATFIELD` frames (about 1.3 GB per association), reduced by the calibration cascade in `reduce.sbatch --masters`.

#### Download job

`download.sbatch`, after `eso_binaries_dl`:

- array over datasets, `%3`; 12 h, 2 CPUs, 8 GB; submitted from a trevor node;
- four parallel streamed downloads from `https://dataportal.eso.org/dataPortal/file/{dp_id}`, with `.part` renames;
- each size checked against the manifest, then a `sha256sum` list written;
- all data are public; no credentials.

Order: nothing is downloaded before the maintainer approves, since D1 (5.34 GB plus masters, under 7 GB) is above the 5 GB threshold in AGENTS.md; once approved, D1 first, and D1N, D2 and D3 only after their own approval.

#### Container

1. `gh workflow run eso-pipelines-sif.yml -R benjaminpope/virgil-vlti -f pipelines=matisse`.
2. `gh run download <run-id> -R benjaminpope/virgil-vlti -n eso-pipelines-matisse -D eso-pipelines-sif`.
3. Check the SHA-256 and the `esorex --recipes` list in the artifact.
4. Copy the `.sif` to the cluster's container directory as `eso-pipelines-matisse-<date>.sif` (`containers/README.md`).

#### Explore job

`explore.sbatch`: 30 min, 1 CPU, 4 GB. In `apptainer exec --cleanenv --home $EDPS_HOME -B <cluster data root> $SIF`, it runs:

- `esorex --version` and `esorex --recipes | grep mat_`;
- `esorex --man-page` for `mat_raw_estimates`, `mat_est_tf`, `mat_cal_oifits` and any flux-calibration recipe the image has ([PM] §9.5 says correlated-flux calibration is not supported by the pipeline, so there may be none);
- `edps --help`.

Then `check_products.py` reads ESO's `CALIB_RAW_INT` products and records:

- pipeline version, `INSNAME` and `OI_REVN` of each table;
- whether a calibrator diameter is recorded, and where;
- the row order and `STA_INDEX` per BCD position (§3.1; this decides `sta_index_frame`);
- `ESO CFG BCD MODE` against `BCD1`/`BCD2`, and which device each INS key is (`bcd_key_order`, §3.1): in a raw N template, the shutter order of the chopped photometric exposures at `IN-IN` ([PM] Table 4.1), and in `RAW_INT` products the rows that move at each single-device position;
- the number of N-band products per template, and of files per N template (D1N lists 12 N files per template; `write_matisse_template(band="N")` writes 4);
- rows and MJDs per baseline in a single `RAW_INT` (decides whether the multi-MJD rule of `merge_evidence` can be fatal).

Output: `explore/<date>.txt`, pulled to the Mac.

#### Reduce job

`reduce.sbatch`: an array, one task per (dataset, night, band); 6 h, 4 CPUs, 32 GB (to be measured on D1). `make_sof.py` writes one SOF per template from the manifests: raw `TARGET_RAW`/`CALIB_RAW`, `SKY_RAW`, `*_RAW_RMNREC`, the five masters, and `JSDC_CAT` for calibrators (a required input of `mat_est_tf`, which writes the transfer function `RAW_TF2`, [PM] §9.21; `mat_raw_estimates` lost its own diameter options in 1.4.4, [PM] §2.15). Recipes:

1. `esorex mat_raw_estimates --spectralBinning=5 --corrFlux=FALSE --replaceTel=3 --hampelFilterKernel=10 --coherentAlgo=2 --cumulBlock=FALSE <sof>` for LM (the ν Hor `ESO PRO REC1` values, from a 2.0.0 reduction). The recipe's own default is `spectralBinning=1`; 5 in LM and 7 in N are the release-note recommendation ([PM] §2.14). In 2.5.0 the options are `spectralAverage` and `cpAlgo` ([PM] §2.1), so the command is written per container version, and `hampelFilterKernel=10` is checked against [PM] §9.1.4, which recommends 0. For N, `--spectralBinning=7` and two runs, `--corrFlux=FALSE` and `--corrFlux=TRUE`, into separate directories. Each LM and N run is made twice, with `--replaceTel=3` (the ν Hor value, into `raw_estimates_rt3/`, for parity with the ν Hor products only) and `--replaceTel=0` (into `raw_estimates/`). The `=0` products feed calibration, the PCA and Stage 4–5 acceptance; the `=3` products are never calibrated (Stage 5 step 4 would refuse them for V²).
2. `esorex mat_cal_oifits --tfInterp=0 --tfKeep=0 --cumulBlock=FALSE` (`tfInterp=0` is an average of the night's transfer function, 2 a linear function of time; `cumulBlock` defaults to TRUE in this recipe; [PM] §9.5) on each night's `TARGET_RAW_INT` and `CALIB_RAW_INT`, into `cal_oifits_block/` (per position).
3. The same with `--cumulBlock=TRUE`, which triggers the BCD combination of [PM] §11.9 (combined files carry `DRS.BCD.COMBINED = T`, [PM] §2.8), into `cal_oifits_bcd/`. *Corrected 2026-10-07: this step was described as "the ν Hor `bcdMode=5` value". `bcdMode` is a `mat_merge_results` option taking 0–3 ([PM] §9.13); the value 5 in the ν Hor headers is not documented there (§7.1).*

Option names and defaults are checked against the explore job's man pages before the first run; any change goes into the job's README. Outputs: the products, an `esorex.log` per recipe, and `products.csv` (file, category, band, BCD, chopped, version, rows). `check_products.py` then checks that every product opens with `products.read_exposure`, and that every science exposure has calibrators in its group (§3.3).

#### Synthetic products and tests

None in the package. The jobs check themselves, as above.

#### Refactors and GRAVITY guards

None.

#### Acceptance

For D1: 36 LM science and 60 LM calibrator `RAW_INT` products, the per-position and combined `mat_cal_oifits` products, the explore report, and checksums, all under the cluster data root.

#### Risks and fallbacks

- **The container will not build** (ESO's RPM repository lacks `esopipe-matisse` for the Fedora release, or the build runs out of disk):
    1. build on Fedora 42 (the build argument);
    2. otherwise, ESO's own Docker image, or a `conda` install of the pipeline on the cluster;
    3. otherwise, skip our reduction and use ESO's archived `CALIB_RAW_INT` products (from Raw2Master). The Stage 4 parity gate then becomes a comparison of calibrator transfer functions only.
- **The masters come from another pipeline version** and the recipe refuses them: use the Raw2Raw cascade.
- **Our `CALIB_RAW_INT` differs from ESO's**: a reduction-setting difference, found by diffing the `ESO PRO REC1 PARAMn` cards. Reported, not a blocker.

#### Effort and dependencies

S in code; about a day of queue time per dataset. Depends on the maintainer's approval of the downloads (§7, decision 1 fixes the plan, not the go-ahead).

### Stage 3. Synthetic MATISSE products and an orientation-neutral core

Two PRs: 3a (the synthetic writer) and 3b (the core).

#### Goal

- **3a:** synthetic MATISSE products that carry every feature the later stages must handle.
- **3b:** make `calibrate.py` match calibrators by orientation-normalized stations, without changing any GRAVITY or PIONIER output. Useful alone: GRAVITY and PIONIER products whose baselines appear reversed are matched rather than flagged.

#### Interfaces

**3a.** New functions in `virgil_vlti.synthetic`; `write_product` is not touched.

```python
MATISSE_STATIONS = ((28, "K0", "AT1"), (19, "G2", "AT2"), (13, "D0", "AT3"), (24, "J3", "AT4"))
MATISSE_INPUTS = (1, 3, 5, 7)

def matisse_wave(band="LM", resolution="LOW") -> tuple[np.ndarray, np.ndarray]
    # wave, band [m]. LM LOW: 64 channels, 4.226 -> 2.825 um, decreasing, EFF_BAND 2.19e-8 m.
    # N LOW: 119 channels, 7.937 -> 13.021 um (as in the 2024 imaging-contest files).

def write_matisse_product(
    path, *, mjd, vis2, t3phi,                    # true values, physical order: (6, nw), (4, nw)
    target="SCI", catg="TARGET_RAW_INT", band="LM", resolution="LOW",
    bcd="OUT-OUT", chopped=True, corr_flux=False, pipeline="2.0.0",
    tpl_start="2023-12-12T00:47:05", expno=1, nexp=12,
    transfer=None,        # {(position, peak): T(lambda)}; default 1
    instrumental_t3=None, # {peak: phi_inst(lambda) [deg]}; default 0
    corr_flux_jy=None,    # (6, nw) correlated flux of the source [Jy]
    total_flux_jy=None,   # (4, nw) or (nw,) total flux [Jy]
    visphi=None,          # (6, nw) true differential phase [deg]
    flux_adu=None, vis2err=None, t3phierr=None,
    sta_index_frame="physical",   # or "peak"
    nonchop_bias=0.0,     # fractional V2 bias added when chopped is False
    int_time=None, header=None, insname=None, seed=0,
) -> None

def write_matisse_template(dirpath, *, kind="science", band="LM", **kw) -> list[Path]
```

**HDU layout** written by `write_matisse_product`, as in a 2.0.0 `RAW_INT` product (checked against Stage 2's explore report when it lands):

| HDU | Contents |
|---|---|
| Primary | No data. Cards: `INSTRUME='MATISSE'`; `ESO DET NAME` (`MATISSE-LM` or `MATISSE-N`); `ESO INS DIL NAME` and `ESO INS DIN NAME` (both the resolution); `ESO INS BCD1 NAME`, `ESO INS BCD2 NAME`, `ESO CFG BCD MODE`; `ESO ISS CHOP ST` (`T`/`F`); `ESO DET CHOP ST` (always `F`, to catch misuse); `ESO ISS CHOP FREQ=0.469`; `ESO DET SEQ1 DIT=0.0752`; `ESO DET NDIT=116`; `ESO DET SEQ1 EXPTIME=61.8`; `ESO TPL ID='MATISSE_hyb_obs_ft'`, `ESO TPL START`, `ESO TPL EXPNO`, `ESO TPL NEXP`; `ESO DEL FT SENSOR='GRAVITY'`; `ESO ISS CONF STATION1..4`, `ESO ISS CONF INPUT1..4`, `ESO ISS CONF T1NAME..T4NAME`; `ESO ISS AMBI TAU0 START/END`; `ESO ISS AIRM START/END`; `ESO PRO CATG`; `ESO PRO REC1 ID='mat_raw_estimates'`, `ESO PRO REC1 PIPE ID='matisse/<pipeline>'`, `ESO PRO REC1 PARAM1..16 NAME/VALUE` (the ν Hor list, with `corrFlux` and `spectralBinning` from the arguments); `OBJECT`, `MJD-OBS`, `DATE-OBS` |
| `OI_TARGET` | one row |
| `OI_ARRAY` | `ARRNAME='VLTI'`; four rows of `TEL_NAME`, `STA_NAME`, `STA_INDEX` (28, 19, 13, 24), `DIAMETER` (1.8 or 8.2) and `STAXYZ` |
| `OI_WAVELENGTH` | `INSNAME` (default `'MATISSE'`, replaced by the real name from Stage 2), `EFF_WAVE`, `EFF_BAND` |
| `OI_VIS2` | six rows in the peak order of §3.1, `STA_INDEX` per `sta_index_frame`; `VIS2DATA` = V² × T, or the squared correlated flux when `corr_flux` (review footgun 1; [PM] §7 says `VIS2DATA` "may contain" squared correlated fluxes but gives no unit, so ADU² is an assumption checked in Stage 2) |
| `OI_T3` | four rows in peak order; `T3PHI` = s · φ_true + φ_inst(peak), negated in N for pipelines in [1.5.2, 2.0.0); `T3AMP` zero before 2.5.0 |
| `OI_VIS` | six rows; `AMPTYP='absolute'` or `'correlated flux'` (the [PM] §7 table spells the key `AMPTYPE`, values `'differential'` or `'correlated flux'`; the reader accepts both spellings and refuses a product with neither, until Stage 2 shows what real products carry); `PHITYP='differential'`; `VISAMP` and `VISPHI` with the BCD signs; the N-band flip as for T3PHI |
| `OI_FLUX` | four rows; `FLUXDATA` in ADU ([PM] §7); `CALSTAT='U'`; `FOV=1.0`, `FOVTYPE='RADIUS'` |

`INT_TIME` is EXPTIME unless `int_time` is given (for merge tests). Noise is added only when errors are given, from a seeded generator.

`write_matisse_template(kind="science", band="LM")` writes the L sequence under one `TPL START`: 4 non-chopped exposures (one per position) and 8 chopped (4 `IN-IN`, 4 `OUT-OUT`). `band="N"` writes 4 fringe exposures, one per position, non-chopped. The N photometric exposures are not `RAW_INT` products of their own; Stage 2 confirms the count.

**3b.** In `calibrate.py`:

```python
def station_key(sta: tuple[int, ...]) -> tuple[tuple[int, ...], int]   # sorted key, parity
def reorient(row: dict, sta_from: tuple, sta_to: tuple, kind: str) -> dict   # one row, §3.2
def calibrator_estimates(cals, insname, observable, diameters=None,
                         smooth_degree=None, orientation=None)          # new keyword
```

`orientation` maps each sorted key to the science file's tuple; `None` keeps today's raw-tuple keying exactly. A per-instrument fact `nbin_advice` (GRAVITY 2; MATISSE LOW 5 in LM, 7 in N) replaces the docstring's "twofold oversampled".

#### Algorithm

3b follows §3.2: key by sorted stations; re-orient a calibrator row only when its tuple differs from the science file's; recompute T3 legs from relative positions. `_apply` writes rows in the science file's orientation, so nothing changes for files that already agree.

#### Synthetic products and tests

`tests/test_synthetic_matisse.py` (3a) and `tests/test_orientation.py` (3b):

| Test | Asserts |
|---|---|
| `test_matisse_product_reads_back` | `read_exposure` gives 64 channels, 6 V² and 4 T3 rows, the `INSNAME` and the category |
| `test_matisse_header_round_trips_through_setup` | `matisse.setup` returns what was written, for all four positions |
| `test_peak_tables_match_section_3_1` | the written `STA_INDEX` and signs per position match both tables |
| `test_bcd_signs_cancel_instrumental_phase` | with φ_inst ≠ 0, the equal-weight mean of sign-normalized phases is φ_true; OO plus II alone is not |
| `test_n_band_flip_written_by_version` | pipeline 1.7.9 N products have negated T3PHI; 2.0.0 do not |
| `test_corrflux_written_to_vis2` | `corr_flux=True` puts ADU² in `VIS2DATA` and `AMPTYP='correlated flux'` |
| `test_template_sequence_lm` | 12 files: 4 non-chopped, 4 `IN-IN` and 4 `OUT-OUT` chopped |
| `test_station_key_parity` | (19,28) → ((19,28), +1); (28,19) → ((19,28), −1); triangles in all 6 orders |
| `test_reoriented_triangle_legs` | recomputed (u₁, v₁, u₂, v₂) close the triangle and match a direct computation |
| `test_reversed_calibrator_baseline_is_matched` | science (28,19), calibrator (19,28): V² and closure phase recovered exactly, no flag |
| `test_reversed_triangle_sign` | a calibrator triangle in odd order calibrates with the right sign |
| `test_same_orientation_is_byte_identical` | `orientation` given or `None` give identical arrays when no row differs |

#### Refactors and GRAVITY guards

The keying change is the only change to GRAVITY code paths. Guards G1, G2 and G3 (real-data `FITSDiff` and `.npz` identity; `gravity_viscal` digits).

#### Validation

G3 only. MATISSE validation starts in Stage 4.

#### Acceptance

All tests above pass, and G1–G3 are unchanged.

#### Risks and fallbacks

- **A GRAVITY or PIONIER output changes.** The cause is a calibrator whose tuple differs from the science file's. The change is then a fix (rows once flagged are now calibrated); it is shown on the golden file with old and new values, and Ben decides.
- **The real `INSNAME` and `OI_REVN` differ from the defaults.** The writer takes them as arguments, so Stage 2's report changes only defaults.

#### Effort and dependencies

3a: M; depends on Stage 1 (its tests use `matisse.setup`). 3b: S–M; depends on nothing, and may run in parallel with Stage 1.

### Stage 4. L/M V² and closure phases per BCD position, then combined

#### Goal

Calibrate L/M V² and closure phases from `RAW_INT` products, per group, then combine the positions (review §3c, §8c.2–4). This is the first stage a user can run end to end on MATISSE data.

#### Interfaces

In `matisse.py`:

```python
def relabel_bcd(exposure) -> Exposure      # applies §3.1 when sta_index_frame == "peak"
def visits(exposures) -> dict[tuple[str, str], list[Exposure]]  # (target, TPL START)
def jsdc_diameter(exposure, band) -> tuple[float, float, str] | None  # (θ, σ_θ, column) from JSDC_CAT, band-matched (§7, decision 13)

@dataclass
class Combined:
    insname: str
    wave: np.ndarray
    vis2: Observed                    # calibrate.Observed
    t3phi: Observed
    positions: tuple[str, ...]        # positions present
    scatter: dict[str, np.ndarray]    # per observable: std across positions
    n_exposures: dict[str, int]       # per position
    weights: dict[Path, float]        # per input file, for mode mapping

def combine_bcd(calibrated: dict[str, list[Calibrated]], *, vis2_chopped_only=True,
                require_all=True) -> Combined
def combine_modes(modes_by_file: dict[Path, dict], combined: Combined) -> dict
```

In `calibrate.py`: `prediction_error(est, mode="loo", ..., groups=None)` and `prediction_residuals(..., groups=None)`, where `groups` assigns each estimate to a visit; `None` keeps today's per-file LOO. `calibrate_exposure(..., groups=None)` passes it through.

In `cli.py`:

- `--loo-by {visit,file}`: default `visit` for MATISSE; `file` for GRAVITY and PIONIER.
- `--bcd {per-position,combined,both}`: default `both` (MATISSE only).
- `--vis2-from {chopped,all}`: default `chopped` (MATISSE LM).
- `--require-all-bcd / --no-require-all-bcd`: default require.
- `--pool-chopping-t3`: off by default (§3.3).
- `--match-mat-cal-oifits`: parity mode (the pipeline's averaging, its `tfInterp=0` rule, statistical errors only).
- `--split-t3-modes`: see §6.
- `--window LO HI` becomes repeatable (`action="append"`, `nargs=2`); the intervals are the union. Today it takes one pair (`cli.py`).
- Instrument-aware defaults: window 2.8–4.2 and 4.48–5.2 µm in LM ([PM] §2.12); `--bin 5` advice. `--dual-field`, `--onaxis-mas`, `--aniso-warn` and the `gravity_viscal` group are refused for MATISSE.
- Outputs per §3.4.

In `products.py`: `Exposure.pipeline_diameter` dispatches on instrument (GRAVITY unchanged; MATISSE from `jsdc_diameter`, else `None`, so `--diameter` is needed). `read_tf2(path)` reads the pipeline's `TF2` table, for comparison only (review §8c.6).

#### Algorithm

1. **Read.** Classify; refuse merged, calibrated and `corrFlux=true` products; apply `--n-phase`; apply `relabel_bcd` if needed.
2. **Group** science and calibrators by `group_key` (§3.3). A science exposure with no calibrator in its group gets an `index.csv` row with `calibrated` empty and `status = "no-calibrator"` (§3.4), and is not calibrated across groups.
3. **Calibrate within each group** as for GRAVITY, with the orientation map of Stage 3b:
    - V² from chopped exposures only (`--vis2-from chopped`; [PM] §11.9, §4.1; review §1c.4);
    - closure phases from all exposures, each calibrated within its own (chopped or non-chopped) group by default (§3.3, revised); `--pool-chopping-t3` merges the two for closure phases once Stage 4's D1 check shows no difference ([PM] §4.1 says chopping is unnecessary to measure closure phases, which does not establish that the instrumental closure phase is the same).
4. **TF error.**
    - `--tf-error loo` with `--loo-by visit`: for each calibrator visit v, predict its exposures from the other visits of the group. The residuals of all of v's exposures are the LOO sample for v.
    - `scatter`: the scatter of the visit means about the prediction. When scatter is pooled over positions, each position's prediction is subtracted first, so BCD-dependent transfer-function differences do not enter the pooled scatter.
    - Fewer than three visits in a group: `loo` falls back to `scatter` over the pooled positions of the night, with a warning and `tf_error_method = "scatter-pooled"` in `index.csv`, and the error is never below the literature width of the facts table (Stage 6: [Lo22] Table 5 at the night's τ₀, and the [Va21] App. B floors), recorded as `"scatter-pooled+floor"` when the floor wins. With one visit there is no scatter, so the floor alone is used and the method is `"floor"`; with no calibrator the exposure is not calibrated (§3.4). *Revised 2026-10-07 (conservative review), §7 decision 9.*
    - Binned as in #21.
5. **Per-position output.** One file per science exposure, as for GRAVITY, written in sorted orientation (`canonical_orientation = "sorted"`).
6. **Combination** per cycle (science template), in sorted orientation:
    - **Within a position p**, average the exposures with inverse-variance weights w_pe = σ_pe⁻² / Σ_e σ_pe⁻²: x_p = Σ_e w_pe x_pe and σ_p² = Σ_e w_pe² σ_pe². Closure phases are averaged as unit vectors: x_p = arg Σ_e w_pe exp(i x_pe).
    - **Across positions**, equal weights. V² = (1/n_p) Σ_p x_p over the positions present in the chopped set (`OUT-OUT` and `IN-IN` in L). Closure phase φ = arg Σ_p exp(i φ_p) over all four positions, after the sign of §3.1 (or of §3.2 when the pipeline relabels).
    - **Errors.** The combination is linear in the inputs once the branch is fixed, so σ² = (1/n_p²) Σ_p σ_p². For phases this holds while the spread across positions is small (< 30°). A channel above that is flagged, and its error set to the circular standard deviation.
    - **Scatter.** The scatter across positions is s, per triangle (or baseline) and channel, on the sign-normalized per-position values x_p (after the §3.1 or §3.2 signs are undone). For closure phases it is the circular standard deviation s = √(−2 ln R̄), with R̄ = |(1/n_p) Σ_p exp(i x_p)| (the population form, with no ddof correction, since n_p ≤ 4 and the circular form has none); for V² it is the sample standard deviation with ddof = 1. A linear standard deviation of phases is never used: positions at ±179° would give about 180°. s goes to the `.scatter.npz` sidecar (review §6c.2).
    - **Floor.** Until a Stage 5 BCD-floor mode exists for the setup (band, resolution, telescopes), the combined closure-phase error is the larger of the propagated σ and s/√n_p, and the channels where s/√n_p wins are counted in `index.csv`. Once the mode exists, it replaces this floor; the two are never both applied, so the cross-position scatter is not counted twice. Where the combination works, s is dominated by the ±φ_inst it cancels, so s/√n_p is an upper bound, which is the conservative side. *Revised 2026-10-07 (conservative review): the scatter was a diagnostic only. [Is22] §2.2 takes the standard deviation across the four BCD positions as the closure-phase uncertainty; a residual the combination does not cancel (a sign or order error, or BCD-dependent systematics) shows there first, and the propagated error cannot see it.*
    - **Missing position.** With `--require-all-bcd`, closure phases of that cycle are not combined (per-position files only), with a warning; V² is combined from the positions present.
7. **Modes through the combination.** Each per-file mode vector m_f (from `.tf_modes.npz`, `.modes.npz` or `.template_modes.npz`) enters the combined file with the file's data weight: m = Σ_f a_f m_f, with a_f = w_pe / n_p for the file's position and exposure, and the §3.1 sign for closure phases. Modes with the same name are summed; different names stay separate. `--tf-modes` and `--template` stay mutually exclusive.
8. **Parity mode** (`--match-mat-cal-oifits`): no LOO, no binning, statistical errors only, the `tfInterp=0` rule as `esorex --man-page mat_cal_oifits` describes it (read in Stage 2), and the pipeline's BCD weights. For acceptance only.

#### Synthetic products and tests

`tests/test_matisse_calibrate.py`, all built from `write_matisse_template`:

| Test | Asserts |
|---|---|
| `test_bcd_dependent_tf_removed_within_group` | a different T per position and peak: V² recovered exactly per position |
| `test_groups_never_mix` | a calibrator only in `OUT-IN` never calibrates `IN-OUT` science |
| `test_combination_cancels_instrumental_closure_phase` | a residual φ_inst of 5° per peak cancels in the combined file to 1e-10 |
| `test_inverse_variance_combination_would_fail` | the same data combined with inverse-variance weights leave a residual (guards against a "fix") |
| `test_combined_errors_match_monte_carlo` | 2000 noise draws: combined σ within 2% of the empirical standard deviation |
| `test_missing_position_warns_and_skips_closure_phases` | chopped-only `OUT-OUT` and `IN-IN`: a warning, no combined T3 |
| `test_nonchopped_v2_not_used_by_default` | `nonchop_bias=0.2` has no effect unless `--vis2-from all` |
| `test_closure_phase_wrap` | +179° and −179° across positions combine near 180° |
| `test_scatter_floor_is_circular` | positions at +179°, −179°, +178°, −178°: s ≈ 1.6° (√(−2 ln R̄)), not ≈ 180°, and the floor is s/2 |
| `test_bcd_floor_mode_replaces_scatter_floor` | with a Stage 5 BCD-floor mode for the setup, the combined error is the propagated σ and the mode is written; without it, the s/√n_p floor; never both |
| `test_loo_by_visit_larger_than_by_file` | correlated exposures within a visit: visit-LOO error exceeds file-LOO |
| `test_few_visits_fall_back_to_pooled_scatter` | two visits: `scatter-pooled` recorded, and the error is at least the facts-table floor |
| `test_single_visit_uses_floor` | one visit: `floor` recorded, error equal to the floor |
| `test_closure_phase_sign_against_model` | an injected asymmetric model (a binary of flux ratio 0.3), with φ_inst ≠ 0 and different on each peak: on all four triangles, the combined closure phase equals, in sign and value to 1e-6 deg, the model's closure phase computed directly from the (u, v) of the physical triangle in sorted orientation, without the §3.1 tables (decision 3). A swapped table row or a reversed `bcd_key_order` fails it. The sign gate against the real pipeline is the D1 check |
| `test_modes_carried_through_combination` | a known mode on one position appears with weight w_pe / 4 and the right sign |
| `test_mode_names_have_position_codes` | `tf-vis2-base13-24-OO@<night>`, in sorted stations |
| `test_index_csv_extra_columns_only_for_matisse` | the GRAVITY `index.csv` header is unchanged |
| `test_uncalibrated_science_row_has_empty_path_and_status` | no calibrator in the group: `calibrated` empty, `status = "no-calibrator"`; calibrated rows keep the relative path |
| `test_peak_relabelling` | `sta_index_frame="peak"` products calibrate to the same result as `"physical"` |
| `test_cli_refuses_gravity_only_options` | `--dual-field` with MATISSE inputs is an error |
| `test_multi_interval_window_and_bins_never_span_gap` | two `--window` intervals keep only their channels, a single pair gives today's mask, and with `--bin 5` no bin contains channels from both sides of the gap |

#### Refactors and GRAVITY guards

`pipeline_diameter` dispatch; `groups=None` default in `prediction_error`; an instrument-aware CLI; multi-interval windows: `calibrate.channel_window(wave, window_um)` accepts either one `(lo, hi)` or a list of `(lo, hi)` intervals (a single pair behaves as today), and binning (`bin_channels`) works within each interval, so a bin never spans a gap such as the ozone band at 4.2–4.48 µm (LM) or 9.18–9.56 µm (N). Guards G1–G3.

#### Validation

On D1 (Stage 2):

1. Run `virgil-vlti-calibrate --match-mat-cal-oifits` on the container's `RAW_INT` products, with the calibrators `mat_cal_oifits` used (from the products' JSDC cards; [PM] App. lists `PRO.JSDC.NMAG/LMAG/NFLUX/LFLUX`, and the calibrator-name card is found in Stage 2), and compare with `cal_oifits_block/` and `cal_oifits_bcd/`.
2. Compare calibrator closure phases of chopped and non-chopped exposures at `OUT-OUT` and `IN-IN`, per triangle; `--pool-chopping-t3` is offered as a default only if they agree within their errors (§7.1).
3. Run the normal mode and compare with the consortium-merged files (kept locally; other calibrators, merged). A report only. The consortium merge weights V² and VISAMP by inverse variance across files and positions, with error √(1/Σw) and no scatter term, and averages phases with equal weights ([MT] `libPostTools.py` at b44f965, `mat_mergeOifits`), so its V² errors are expected to be smaller than ours.

The comparisons live in a private validation suite, as the GRAVITY check does (PROVENANCE.md, check 1).

#### Acceptance

z = (x_ours − x_pipe)/σ_pipe, per (exposure, baseline or triangle, channel) inside the window, where σ_pipe is the pipeline's error. Closure-phase differences are wrapped to (−180°, 180°] before dividing.

- Per position, against `cumulBlock=FALSE`: median |z| ≤ 0.1 for V² and closure phase; max |z| ≤ 0.5.
- Combined, against `cumulBlock=TRUE`: median |z| ≤ 0.2; max |z| ≤ 1.0.

These thresholds are decision 3 (§7).

#### Risks and fallbacks

- **Parity fails because the pipeline averages differently** (inverse-variance weights across positions, or another TF interpolation; [PM] §9.5 documents `tfInterp=0` as an average and 2 as linear, but not how the average is weighted). First compare the calibrator transfer functions alone: ours (calibrator raw over model) against the pipeline's `TF2` table, which needs no averaging rule. If they agree to 0.1σ, the remaining difference is the averaging: it is documented and the gate moves to the transfer functions. If they differ, the cause is the diameter or the input selection, checked against the `JSDC_CAT` row.
- **The pipeline writes peak-ordered `STA_INDEX`.** `relabel_bcd` is already planned.
- **Too few calibrator visits per position** on a night: the pooled scatter of step 4, recorded.

#### Effort and dependencies

L. Depends on Stages 1, 3a and 3b. Acceptance needs Stage 2 (D1).

### Stage 5. The MATISSE error model: pairs, calibrator PCA, templates

#### Goal

Empirical nuisance modes and widths for MATISSE L/M, handed to virgil as `.npz` (§6).

#### Interfaces

- `pairs.find_pairs(exposures, insname, max_gap_min, *, same=None)`: `same` is a key function (default `None`, today's behaviour); MATISSE passes `matisse.group_key`. New `pairs.cross_bcd_pairs(exposures, insname, max_gap_min)` pairs exposures of one cycle in different positions.
- `calpca.calibrator_residuals(exposures, insname=None, grid=None, *, exclude=(), window_um=None, field=None, references="other-star", group_key=None, diameters=None, **kw)`: `insname=None` takes `default_insnames`; `group_key=None` keeps `_group_key`; `diameters` overrides `pipeline_diameter`. `calpca_cli` gains `--diameter` and `--instrument` (inferred by default).
- `templates.Template.config` gains `instrument`, `band` and `bcd` (`"combined"` or a position code). `Template.mismatches(header)` compares them through `INSTRUMENT_MODULES`. `FORMAT_VERSION = 2`; version-1 files load with `instrument="gravity"`.
- `template_name(config, epoch, instrument="matisse")` gives, for example, `matisse_lm_low_at_bcd-combined_2023-2025.npz`.
- `matisse.photometric_modes(exposure, eps_p, flux)`: modes `phot-tel<sta>@<file>`, one per telescope, on its three baselines.

#### Algorithm

1. **Pair test** within a group, as for GRAVITY, and across positions within a cycle. The cross-position residual, after the §3.1 signs, estimates the per-position instrumental offset; its rms (circular for closure phases) is the BCD floor, written as a mode per (band, resolution, telescopes). Once it exists for a setup, it replaces the Stage 4 s/√n_p floor for that setup; the two are never both applied.
2. **Calibrator PCA** per group: the residuals of each calibrator exposure against the other-star prediction (`references`), stacked over (positions × channels) per baseline, so that a component may be common to all positions. Held-out widths as for GRAVITY (`heldout_widths`).
3. **Generalized channel split.** `channel_correlation` on residuals concatenated over positions (estimates × 4·n_chan). The leading mode is the part common to all positions and channels; the per-position remainder is a second, position-local mode. `--tf-modes` uses this split only if the "shared" option of §7 decision 10 is switched on.
4. **Photometric modes** ([Lo22] eq. A.9). A photometric error ε_p on telescope i changes the V² of its three baselines by m_b = V²_b ε_p / F_i (b ∋ i), since V² ∝ |C|²/(F_i F_j) gives δV²/V² = −δF_i/F_i (the sign is absorbed into the mode). *Checked 2026-10-07:* [Lo22] eq. A.9 writes σ_V²/V² as the quadrature sum of the coherent-flux and photon terms and ε_pi/F_i and ε_pj/F_j, one independent term per telescope, with F_i the source flux collected by telescope i; so one mode per telescope with this normalization reproduces it, with no extra factor. The width needs F_i and ε_p in the same units. `OI_FLUX` is in ADU, so ε_p is estimated in ADU from the pair test or the PCA. The [Lo22] Table 6 values (0.11 / 0.19 / 2.3 Jy in L/M/N on the ATs; 0.008 / 0.016 / 0.08 Jy on the UTs; broadband, from chopped sky-frame statistics, [Lo22] §4.1.4; the paper does not say whether they are 1σ per exposure) are used only once Stage 7 gives a flux scale, and then as a lower bound on the width, never an upper one. A product with `replaceTel` ≠ 0 ([PM] §9.1.4: the telescope's photometry is replaced by the mean of the other three) gets one extra mode common to the three baselines of the replaced telescope, with a width of at least the rms difference between that telescope's flux and the mean of the others, measured on unreplaced products of the same night; without such products the exposure is refused for V² unless forced. Which telescope numbering the digit uses (AT position in `ISS CONF`, or input) is checked in Stage 2 (§7.1).
5. **Template** from the PCA components, as `build_template`, per (band, resolution, telescopes, bcd).

#### Synthetic products and tests

`tests/test_matisse_error_model.py`:

| Test | Asserts |
|---|---|
| `test_pairs_stay_within_group` | `same=group_key` never pairs two positions |
| `test_cross_bcd_pairs_measure_instrumental_offsets` | injected per-position offsets recovered to 10% |
| `test_photometric_mode_recovered_by_pca` | a mode ∝ V²·ε_p/F_i on three baselines is the first component |
| `test_shared_mode_across_positions` | a mode common to all positions is found once in the concatenated split |
| `test_replace_tel_mode` | `replaceTel=3` adds the common mode on AT3's baselines |
| `test_template_v1_still_loads` | the shipped GRAVITY template loads and matches as before |
| `test_matisse_template_mismatch` | an `OUT-OUT` product against a `combined` template is a mismatch |

#### Refactors and GRAVITY guards

The `calpca` and `templates` defaults keep today's behaviour. The shipped GRAVITY template is unchanged (`test_shipped_pilot_template`). G1–G3.

#### Validation

D1's calibrators (five visits over two nights; HD 14682 twice), and ESO's archived `CALIB_RAW_INT` products for the same templates (Stage 2), which give a second, independent reduction.

#### Acceptance

On D1, the held-out χ²/N of calibrator residuals falls from its diagonal value towards 1 with k ≤ 3 modes, and the per-telescope fraction (`telescope_fraction`) is reported. The widths agree in order of magnitude with [Lo22] Table 5 at the night's τ₀ (Table 5 gives errors on V, not V², so the fractional width 2σ_V/V_cal is compared) (`ESO ISS AMBI TAU0`: 7.3, 8.6 and 3.6 ms in the three ν Hor files).

#### Risks and fallbacks

Five calibrator visits may be too few for a stable PCA. Then no template ships, the widths are reported, and an archive-wide PCA (OiDB "MATISSE calibrators", d259d56c-7304-4131-99d5-6c7304c13119, if its products are `RAW_INT`) becomes a later stage.

#### Effort and dependencies

M. Depends on Stage 4; acceptance on Stage 2.

### Stage 6. Instrument facts and default widths

#### Goal

The facts virgil needs from virgil-vlti, as a versioned table with a source for each number.

#### Interfaces

`src/virgil_vlti/facts/matisse.json` (package data), loaded by `virgil_vlti.facts.load(instrument) -> dict`:

| Key | Content | Source |
|---|---|---|
| `resolving_power` | per band and mode, R of the line-spread function before binning. [Lo22] Table 2: LM LOW 31.5, N LOW 31.5, LM MED 499, L HIGH 979, L and M VHIGH 3370, N HIGH 218, with no reference wavelength. [UM] Table 3: LM 34 / 506 / 959 / 3666 (HIGH+), N 30 / 218; [UM] §1 quotes 3300 for HIGH+. Both are stored, with the source; the default is [Lo22], and the spread is a reason to fit R (§8 issue 1) rather than fix it | [Lo22] Table 2; [UM] Table 3, §1 |
| `spectral_binning` | recommended sliding-bin width in pixels (LM 5, N 7); the recipe default is 1. Read from the product's `REC` params, never from this table, when present | [PM] §2.14 (1.4.9 notes), §9.1.4; ν Hor `REC1` |
| `pinhole_lambda_over_d` | LM 1.5; N 2.0, the standard pinholes. [UM] Table 2 also offers LM 1.0 and 2.0, N 1.5 and 2.5, and N slits; the product's own spatial-filter key wins, and a product without one is refused for primary-beam use | [Lo22] Table 2 notes, §2.4; [UM] §2.2.2, Table 2 |
| `telescope_diameter_m` | AT 1.8; UT 8.2 | ESO |
| `windows_um` | LM [[2.8, 4.2], [4.48, 5.2]]; N [[8.0, 9.18], [9.56, 13.0]]: the unflagged ranges | [PM] §2.12 (1.5.1 notes) |
| `ozone_um` | [9.18, 9.56], the gap between the N windows; [PM] §2.12 does not name it as ozone | [PM] §2.12 |
| `widths` | per band and τ₀ bin: V², closure phase, correlated flux. [Lo22] Table 5 gives σ_V (visibility, not V²) and σ_CP at τ₀ = 3.2, 6.8 and 7.5 ms, as the rms of the night's transfer-function variation; since a transfer-function error is multiplicative on the science target, the table stores the fractional width f_V² = 2σ_V/V_cal (V_cal ≈ 1 for Table 5's bright calibrators, so f_V² = 2σ_V) and it is applied as σ(V²_sci) = f_V² · V²_sci. The [Va21] floors are converted into the same terms: 1° in L closure phase is absolute; 8% in N correlated flux is relative (f = 0.08); 0.03 in L visibility is read as absolute in V, so σ(V²_sci) = 2V_sci · 0.03 (whether [Va21] means absolute or relative is in §7.1). The width used is the larger of the fractional and the floor, each first converted to a σ on V²_sci. For a night's τ₀, the bin used is the nearest Table 5 τ₀ at or below it (3.2 ms below 3.2 ms), with no interpolation or extrapolation. [Lo22] Table 4 gives high-flux minimum errors (V 3e-3 / 5e-3 / 3e-3, CP 0.1° / 0.2° / 0.1° in L/M/N) | [Lo22] Tables 4–6, §4.1.2; [Va21] App. B floors (0.03 in L visibility, 1° in L closure phase, 8% in N correlated flux) |
| `photometric_error_jy` | AT/UT × L/M/N | [Lo22] Table 6 |
| `sta_index_frame` | `"physical"` or `"peak"` | Stage 2 explore report |
| `bcd_key_order` | `"BCD1=IP13"` or `"BCD1=IP57"`; absent until Stage 2, and its absence is a refusal (§3.1) | Stage 2 explore report; [PM] Table 4.1; [MT] `mat_removeBCD` |
| `canonical_orientation` | `"sorted"` | §3.2 |
| `nbin_advice` | LM LOW 5; N LOW 7 | §3 |

*Corrected 2026-10-07: the resolving powers were "LM LOW 31.5 at 3.5 µm; N LOW about 30; LM MED about 500; LM HIGH about 960"; [Lo22] Table 2 gives no reference wavelength, N LOW 31.5 and L HIGH 979, and [UM] Table 3 differs (above).*

Functions:

- `matisse.response(setup, wave) -> (kernels, weights)`: the two-kernel response of §8 issue 1, for each channel.
- `write.write_widths(path, widths)`: writes widths beside the modes in the `.npz`, under `"<insname>/widths/<mode>"`, a new key family that `stack_modes` ignores.

#### Algorithm

A LOW channel's effective response is a Gaussian line-spread function of FWHM λ/R, acting on the complex visibility, followed by a top-hat average of V² over `spectral_binning` pixels, since the pipeline averages estimates, not fields. `EFF_BAND` is not used (review footgun 4).

#### Synthetic products and tests

`tests/test_facts.py`:

| Test | Asserts |
|---|---|
| `test_facts_have_sources` | every leaf has a `source` |
| `test_units` | µm, m, Jy and degrees as declared |
| `test_effective_response_width` | the 5-pixel average of an R = 31.5 Gaussian has the computed second moment |
| `test_widths_are_fractional_vis2` | a Table 5 σ_V is stored as f_V² = 2σ_V/V_cal, and on a science V² = 0.01 gives σ = 0.01 · f_V², not an absolute 2σ_V |
| `test_va21_floors_converted` | the 0.03 visibility floor gives 2V_sci · 0.03 on V², the 8% floor is relative, 1° absolute; the max is taken in V² units |
| `test_tau0_bin_lookup` | τ₀ = 5 ms uses the 3.2 ms bin; 2 ms uses 3.2 ms; 10 ms uses 7.5 ms |
| `test_widths_written_beside_modes` | `stack_modes` output is unchanged |

#### Refactors and GRAVITY guards

The same table format later serves GRAVITY's PLAN item 4 (wavelength priors) and PIONIER. No GRAVITY file is written now.

#### Validation

None beyond the tests: the facts come from the literature, and the Stage 2 explore report.

#### Acceptance

Every number has a source line, and the units are checked by test.

#### Risks and fallbacks

R from [Lo22] is nominal. Stage 2's reduction of a lamp or a spectral standard could measure it, but that waits until virgil's smearing exists.

#### Effort and dependencies

S. Depends on Stage 1, and on Stage 2 for `sta_index_frame`. Runs in parallel with Stages 3–5.

### Stage 7. N band: correlated-flux calibration

#### Goal

Calibrated N-band correlated fluxes, the consortium's preferred N observable ([Va21] §3.1; [Is22] §2.2), and N-band closure phases.

#### Interfaces

```python
# matisse.py
@dataclass(frozen=True)
class Spectrum:
    wave: np.ndarray          # [m]
    flux_jy: np.ndarray
    rel_error: float          # 1-sigma relative error of the spectrum's scale
    source: str               # "user", "rayleigh-jeans" or "mdfc"

def calibrator_spectrum(name, *, csv=None, n_flux_jy=None, mdfc=None) -> Spectrum
def calibrate_corrflux(sci, cals, insname, spectra: dict[str, Spectrum], *,
                       diameters, transmission=None, **calibrate_kw) -> Calibrated
def spectrum_modes(sci, calibrated, spectra, *, order=1) -> dict
```

CLI flags on `virgil-vlti-calibrate`:

- `--corrflux` (N only);
- `--cal-spectrum NAME=FILE.csv`, repeatable; CSV columns `wave_um,flux_jy[,err_jy]`;
- `--cal-nflux NAME=JY`: the Rayleigh–Jeans fallback, from the flux at 10 µm;
- `--transmission FILE.csv`: columns `wave_um,t`, the transmission at airmass 1;
- `--spec-order 1`.

#### Algorithm

1. **Inputs.** `corrFlux=true` `RAW_INT` products (Stage 2's second N run). The raw correlated flux C(λ) is `VISAMP` (ADU), per position. Phase convention per Stage 1.
2. **Transfer function** per calibrator exposure ([Va21] eqs. 2–3): T(λ) = C_cal(λ) / (F_cal(λ) · |V_UD(λ; θ_cal, B)|), with V_UD = 2 J₁(x)/x, x = π θ B/λ (`calibrate.uniform_disk_visamp`), and B from the row's u, v.
3. **Airmass correction** (optional): T(λ) ← T(λ) · t(λ)^(X_sci − X_cal), with X the mean of `ESO ISS AIRM START/END` and t the user's transmission.
4. **Calibrate.** Interpolate T to the science time with `predict`, per group, as for V²; then F_corr(λ) = C_sci(λ) / T(λ).
5. **Errors.** (σ_F/F)² = (σ_C/C)² + (σ_T/T)², with σ_T/T from `--tf-error` (visit LOO or scatter, in log space). The calibrator spectrum and diameter are left out of these errors and become modes.
6. **`spectrum_modes`.** For each calibrator c and night, and k = 0..order: m_k(λ) = F_corr(λ) · ε_c · w_c · P_k(x(λ)), with P_k a Legendre polynomial on the window, ε_c the spectrum's `rel_error`, and w_c the interpolation weight of c. Names: `spec-<c>@<night>-k<k>`. The diameter mode is `diameter_modes` with V_UD.
7. **Spectrum sources.**
    - A user CSV, with `rel_error` from its error column, or 5%.
    - Rayleigh–Jeans, F ∝ λ⁻² scaled to `--cal-nflux`, with `rel_error` 10% and a warning for K and M giants (SiO absorption near 8 µm).
    - A catalogue, once the reader of §7 decision 5 exists. The MDFC (Cruzalèbes et al. 2019) gives diameters and broadband fluxes, not N-band spectra; spectra come from van Boekel et al. 2004 or the consortium's STARSFLUX and `calib_spec_db` files, as [MT] `libFluxCal_STARSFLUX_new.py` and [Va21] §3.1 (van Boekel spectra, MDFC diameters) use them.
8. **Flags.** Ozone channels (9.18–9.56 µm) are flagged. Channels whose calibrated correlated flux is below the coherent-flux bias limit (about 5 Jy at 8.5 µm and 8 Jy at 11 µm; [Lo22] §4.2, Fig. 18, from AT data; [ESO-news], August 2020, gives 6–8 Jy for the bias from the piston correction) are flagged, not given larger errors. The limit is applied to UT data too until a UT value is cited (§7.1); [Lo22] §6 says GRA4MAT mitigates the bias, but no GRA4MAT limit is quoted, so it is not relaxed for GRA4MAT data.
9. **BCD combination** as in Stage 4, with equal weights.
10. **Output.** `OI_VIS` with `AMPTYP='correlated flux'`, `VISAMP` in Jy, and `VISPHI` empty until Stage 8; `OI_T3` closure phases calibrated as in L. virgil reads it as `"corrflux"` (virgil#150, #198).

#### Synthetic products and tests

`tests/test_matisse_corrflux.py`:

| Test | Asserts |
|---|---|
| `test_known_spectrum_and_diameter_give_exact_corrflux` | noise-free: F_corr equals the input to 1e-10 |
| `test_spectrum_error_is_a_mode_not_diagonal` | a 5% scale error in F_cal shows in `spec-*-k0`, not in σ |
| `test_spectrum_slope_mode` | a tilt shows in `k1` |
| `test_diameter_mode_for_resolved_calibrator` | the finite difference matches the analytic dV/dθ |
| `test_airmass_correction` | with t(λ) and X_sci ≠ X_cal, the correction recovers the input |
| `test_ozone_and_bias_flags` | channels flagged, errors untouched |
| `test_flipped_n_phases_refused_then_negated` | 1.7.9 is refused; with `--n-phase negate` a synthetic binary's position angle is right |
| `test_corrflux_written_as_amptyp_correlated_flux` | `virgil.read_oifits` reads it as `corrflux` (skipped if virgil is absent) |

#### Refactors and GRAVITY guards

`Observable.model` gains an optional spectrum argument; GRAVITY passes none. The per-target `Spectrum` reaches it through a `spectra` dict (calibrator name to `Spectrum`) passed to `calibrator_estimates` like `diameters`, which calls `spec.model(table, wave, diameter, spectrum)` per calibrator. G1–G3.

#### Validation

- D2 (HD 163296, 2019): our F_corr against [Va21]'s published values, and against an independent flux calibration of the same inputs: the consortium's `mat_fluxCal.py` ([MT]), since the ESO pipeline does not calibrate correlated flux ([PM] §9.5).
- D1N: ν Hor's N-band closure phases against its L-band ones. A consistency check, not a gate.

#### Acceptance

- Against [MT] `mat_fluxCal.py` on D2, with the same calibrator spectrum: |ΔF|/σ_ours ≤ 1 in 95% of channels (§7, decision 6, revised).
- Against [Va21]: within their 8% floor, if their products can be had (§7, decision 6).

#### Risks and fallbacks

- **No calibrator spectrum:** the Rayleigh–Jeans fallback, with its larger mode.
- **[Va21] products not on OiDB:** compare with their figures, read by hand at a few wavelengths (decision 6).
- **The consortium tool cannot be run on D2:** no gate is passed; Stage 7 stays unvalidated and says so in PROVENANCE, and the figure comparison is reported. (The ESO pipeline has no correlated-flux calibration, [PM] §9.5, so it is not a fallback.)

#### Effort and dependencies

M–L. Depends on Stages 1, 3a and 4, and on Stage 2 for D2 and D1N.

### Stage 8. N-band V², differential phases, MED and HIGH

#### Goal

The remaining observables: N-band V² from non-simultaneous photometry; the two BCD-correctable differential phases ([PM] §11.9); MED and HIGH.

#### Interfaces

- `cli.py`: `--vis` for MATISSE writes `OI_VIS` with `PHITYP='differential'` for baselines (1,3) and (5,7) only (§3.1); `--n-vis2` enables N-band V².
- `matisse.combine_bcd(..., observables=("vis2", "t3phi", "visphi"))`.
- Facts for MED and HIGH (R, windows).

#### Algorithm

1. **N-band V²**: calibrated as in L, with a warning that photometry is taken after the fringes, so a change of transfer function between them is in the data ([Va21] §3.1). The TF error includes it by construction, since it is in the calibrator residuals.
2. **Differential phases**: calibrated per position (φ_sci minus the predicted φ_cal), then combined over positions with the §3.1 signs and equal weights. Only (1,3) and (5,7) cancel their instrumental term, so the other four are written per position only. virgil's `DifferentialPhase` removes an offset and a slope in wavenumber (#198, #232).
3. **MED and HIGH**: the same code with new facts.

#### Synthetic products and tests

| Test | Asserts |
|---|---|
| `test_differential_phase_bcd_cancellation` | an instrumental differential phase on (1,3) cancels |
| `test_other_baselines_not_combined` | (3,5), (3,7), (1,5), (1,7) appear per position only |
| `test_n_vis2_warns` | the non-simultaneous photometry warning is raised |
| `test_med_lm_reads_and_calibrates` | a synthetic MED grid calibrates exactly |

#### Refactors and GRAVITY guards

None beyond G1–G3.

#### Validation

One MED or UT dataset, chosen from Stage 2's headers (D3 is UT; the resolution comes from `ESO INS DIL/DIN NAME`).

#### Acceptance

Parity with `mat_cal_oifits`, with the Stage 4 thresholds.

#### Risks and fallbacks

There may be no public MED or HIGH dataset with calibrators in all positions. Then test on synthetic data only, and say so in PROVENANCE.

#### Effort and dependencies

M. Depends on Stages 4 and 7.

## 5. Dependencies, parallel work and effort

```text
Stage 1 ──┬──> 3a ──┐
          │         ├──> Stage 4 ──┬──> Stage 5
Stage 3b ─┼─────────┘              └──> Stage 7 ──> Stage 8
          └──> Stage 6
Stage 2 (cluster, needs approval) ··> acceptance of 4, 5, 7, 8
virgil issues 1-5 (§8) ··> use of the Stage 6 facts and of closure-phase TF modes
```

| Stage | Size | Depends on | Can run in parallel with |
|---|---|---|---|
| 1 | S–M | — | 2, 3b |
| 2 | S, plus queue time | approval | everything |
| 3a | M | 1 | 3b, 6 |
| 3b | S–M | — | 1, 3a, 6 |
| 4 | L | 1, 3a, 3b; acceptance needs 2 (D1) | 6 |
| 5 | M | 4; acceptance needs 2 | 7 |
| 6 | S | 1 (2 for `sta_index_frame`) | 3–5 |
| 7 | M–L | 4; acceptance needs 2 (D2, D1N) | 5 |
| 8 | M | 4, 7 | — |

The critical path is 1 → 3a → 4 → 7 → 8. Stage 2 should start first, since its queue time is the longest single delay.

## 6. The error model: what differs from GRAVITY

The GRAVITY machinery carries over: a leave-one-out transfer-function error per group (`calibrate.prediction_error`), binned so that it is not averaged down (#21); the pair test (`pairs.py`); the calibrator PCA (`calpca.py`); and gain templates for virgil's `with_gains(modes=...)` (#22; virgil#180). Five things differ.

1. **BCD-dependent systematics.** Each position has its own transfer function and instrumental closure phase ([Lo22] §3.3, §4.1.3). Within a group they are calibrated like GRAVITY's. Across groups, the equal-weight combination cancels what calibration leaves (§3.1). The cross-position scatter s (circular for closure phases) is a floor on the combined closure-phase error, s/√n_p, until the pair test of Stage 5 gives a BCD-floor mode for the setup; the mode then replaces the floor, and the two are never both applied ([Va21] App. B; [Is22] §2.2).
2. **Visits, not files.** The exposures of one calibrator template sample the transfer function once. Leave-one-visit-out (§3.3, Stage 4) gives an honest LOO error; per-file LOO underestimates it.
3. **Chopping and photometry.** Chop-subtraction residuals are the most important contribution to the photometric error, and so to the visibility error ([PM] §4.1). In L, tip-tilt and Strehl fluctuations are also a dominant cause, which is why simultaneous photometry (SiPhot) is used there ([PM] §4.2). *Corrected 2026-10-07: this said "In L, chop residuals are the main photometric error ([PM] §4.1)"; §4.1 is not specific to L.* Both enter V² as V²_b ε_p/F_i on the three baselines of telescope i ([Lo22] eq. A.9, checked in Stage 5). The mode's shape is per telescope, as for GRAVITY; its width scales with 1/F_i and is written per exposure (Stage 5).
4. **The thermal background in N.** It sets the noise, and errors of the per-frame OPD estimate for faint sources bias the coherent flux, limiting reliable coherent fluxes to about 5 Jy at 8.5 µm and 8 Jy at 11 µm ([Lo22] §4.2, Fig. 18, on AT data). *Corrected 2026-10-07: this said "biases the coherent flux low ... on the ATs"; [Lo22] gives neither the sign of the bias nor the telescopes in that sentence.* That bias is flagged, not given an error. N-band widths depend strongly on τ₀ ([Lo22] Table 5: closure phase 1.75° at τ₀ = 3.2 ms against 0.29° at 7.5 ms; visibility 0.045 against 0.015, an error on V, not V²), so the widths are a function of τ₀.
5. **Absolute correlated-flux calibration.** The calibrator's spectrum and diameter multiply everything it calibrates, smoothly in λ, so they are modes, not diagonal errors ([MT] `libFluxCal_STARSFLUX_new.py` leaves both out of its errors). A per-exposure achromatic gain of several per cent ([Va21] App. B: 8% in good, about 20% in poor conditions) is the transfer-function mode.

The pipeline errors are statistical only (review footgun 2). They are kept as the diagonal, and everything above is added as modes or as the binned LOO term.

**What virgil can take today.** `OIData.with_gains(..., modes=...)` accepts modes of shape (n_mode, n_sample) that span frames, so V² and correlated-flux modes pass. `OIData.with_closure_offsets(..., modes=...)` requires each closure-phase mode to lie within one frame. It therefore refuses the per-night `tf-t3phi` modes GRAVITY already writes, and would refuse MATISSE's. Until §8 issue 5 is done, `--split-t3-modes` splits each closure-phase mode into per-frame pieces (losing the correlation between frames); without the flag, the modes are written for virgil-vlti's own use only.

## 7. Decisions (2026-10-07)

The maintainer accepted every recommendation below on 2026-10-07. MATISSE is outside their core expertise, so all 13 decisions await review by a MATISSE expert. Each entry gives the question for context, then the decision. Accepting decision 1 settles the plan; it is not a go-ahead for any download, which still needs the maintainer's explicit approval (Stage 2).

A conservative literature review on 2026-10-07 checked the MATISSE facts against [PM] issue 2.5.0, [UM] issue 1.0, [Lo22] and [ESO-news]. Decisions it changed carry a *Revised 2026-10-07 (conservative review)* line under the maintainer's accepted decision, with the reason; the accepted text is kept. It also changed two defaults outside the list: chopping is a group key for closure phases too (§3.3), and a disagreement between `ESO CFG BCD MODE` and the INS keys is a refusal (Stage 1).

**1. Downloads and container: may Stage 2 fetch data, and may a `pipelines=matisse` container build be triggered?**
- Options considered: (a) D1 only, L/M (5.3 GB plus under 1 GB of masters), then decide; (b) D1 with N (51 GB) plus D2 (14 GB); (c) everything listed (more than 105 GB with D3 and D4, whose calibrators are still to be counted).
- (a) unblocks Stages 4 and 5 within a day and tests the container cheaply; (b) and (c) are needed only for Stages 7 and 8.
- **Decision:** (a) now, with the container build; (b) when Stage 7 starts.

**2. Pre-2.0 N-band phases: refuse by default, or negate by default?**
- Refuse: no silent 180° rotation, but every user of 2018–2023 N-band data must pass `--n-phase negate`.
- Negate: works out of the box, but trusts a version range from two sources that no one here has tested on a product.
- **Decision:** Refuse by default until D2 (pipeline before 2.0) or D4 shows a known orientation coming out right after negation; then make negation the default.
- *Revised 2026-10-07 (conservative review):* refusal stays the default even after D2 or D4. The check is then recorded in PROVENANCE, and negation becomes the default only with the MATISSE expert's sign-off. N-band products from before 1.5.2 are also refused by default (`"unknown"`, not `"standard"`); they have their own option, `--n-phase-unknown`, so `--n-phase negate` never rotates them. Reason: [PM] §2.5 gives 1.5.2 as the start of the inconsistency, but the 1.5.2 release notes ([PM] §2.11) and the ESO news item ([ESO-news], December 2023) do not confirm the lower bound, and a single dataset checks one version, not the range.

**3. Parity thresholds: are median |z| ≤ 0.1 per position and ≤ 0.2 combined (max 0.5 and 1.0) the right gates for Stage 4?**
- Tighter (0.05): catches subtle averaging differences, but may fail on details of the TF interpolation we cannot reproduce.
- Looser (0.3): always passes, and hides a small sign error in one triangle.
- **Decision:** Keep 0.1 and 0.2, with the transfer-function comparison as the fallback gate if the pipeline's averaging cannot be reproduced (Stage 4 risks).
- *Revised 2026-10-07 (conservative review):* two additions. (i) A sign gate separate from |z|: on every triangle where the pipeline's combined closure phase differs from zero by more than 3σ, ours has the same sign (`test_closure_phase_sign_against_model` on synthetic data, against the model's closure phase computed independently of the §3.1 tables; against the real pipeline on D1 by hand). A sign or BCD-order error in one triangle can hide inside a median |z|. (ii) Moving to the fallback gate needs a written reason in the PR, naming the averaging difference found, and the maintainer's approval; it is not automatic. Reason: [PM] §11.9 gives the combination formula for one triangle only, and documents neither the averaging weights nor the `tfInterp=0` average in detail ([PM] §9.5), so a parity failure is as likely to be a convention error as an averaging difference.

**4. N-band order: correlated flux first, or V² first?**
- Correlated flux first: what the consortium uses and N-band science needs, but it needs calibrator spectra (decision 5).
- V² first: reuses Stage 4 unchanged, but N-band V² is biased by non-simultaneous photometry, and many targets are below the photometric limit.
- **Decision:** Correlated flux first (Stage 7); V² as an option in Stage 8.

**5. Calibrator spectra: which source first for Stage 7?**
- A user CSV only: no catalogue code and honest errors, but users must find spectra.
- MDFC (Cruzalèbes et al. 2019): covers most MATISSE calibrators, but is a catalogue to ship or fetch, and has not been read here. Cohen templates: high quality, but few stars. *Corrected 2026-10-07: the MDFC gives diameters and broadband mid-infrared fluxes, not spectra; the consortium's spectra are van Boekel et al. 2004 and STARSFLUX ([MT] `libFluxCal_STARSFLUX_new.py`; [Va21] §3.1).*
- **Decision:** A user CSV plus the Rayleigh–Jeans fallback with a 10% mode; an MDFC reader later, once its format has been read.
- Not revised. The Rayleigh–Jeans fallback is already opt-in (`--cal-nflux`), and its 10% width is listed in §7.1 for a check against K and M giants' SiO feature.

**6. Published products: if OiDB lacks the Varga et al. HD 163296 products, is a comparison with the paper's figures enough for Stage 7?**
- Figures: about 5% reading precision; enough to catch a wrong factor or sign, not a 1σ difference.
- The container's flux calibration on the same inputs: a tight check of our arithmetic, but not of the method.
- **Decision:** Gate on the container comparison and report the figure comparison; the authors are not emailed.
- *Revised 2026-10-07 (conservative review):* gate on the consortium's `mat_fluxCal.py` ([MT]) run on the same inputs, not on the container. Reason: the ESO pipeline does not calibrate correlated flux ([PM] §9.5), so the container has nothing to compare against. If the consortium tool cannot be run, Stage 7 has no gate and is recorded as unvalidated rather than passed on the figures.

**7. L/M first: should Stages 2 and 4 treat L/M as a complete deliverable before any N-band download?**
- Yes: 5.3 GB instead of 51 GB, and Stage 4 needs only L/M.
- No: one download now, and the N-band phases of the same exposures cross-check Stage 4's BCD handling.
- **Decision:** Yes; Stage 4 does not need the N-band cross-check to pass.

**8. ESO's archived `CALIB_RAW_INT` products: should we use them?**
- As a parity check of our container's `mat_raw_estimates`: free, since CalSelector returns them with each frame.
- As inputs to the Stage 5 PCA: more calibrators without our own reductions, but their pipeline version may differ from ours, so they must be grouped separately.
- **Decision:** Use them for both, grouped by pipeline version, with the version recorded in `index.csv`.
- *Revised 2026-10-07 (conservative review):* group by pipeline version **and** every `mat_raw_estimates` parameter in `REDUCTION_PARAMS` (with `replaceTel` and the 2.5.0 renames), and never mix ESO's products with ours in one group. Reason: the parameters change the estimator itself (`cpAlgo=2` gives coherent rather than bispectrum closure phases; `spectralBinning`/`spectralAverage` change the channel response; `replaceTel` replaces one telescope's photometry; [PM] §2.1, §9.1.4), and archived products need not share our container's settings.

**9. Should leave-one-visit-out be the MATISSE default?**
- Visit: honest when the exposures of a template share their transfer function; larger errors, and only 2–3 LOO samples per group.
- File: matches GRAVITY, but can underestimate the error by up to √12 for a 12-exposure template.
- **Decision:** `--loo-by visit` for MATISSE; keep `file` for GRAVITY (byte identity) and revisit GRAVITY separately.
- *Revised 2026-10-07 (conservative review):* with fewer than three visits in a group, the pooled scatter is floored by the literature widths (Stage 6 facts: [Lo22] Table 5 at the night's τ₀, as the fractional width 2σ_V/V_cal applied to the science V², and the [Va21] App. B floors), and one visit uses the floor alone (Stage 4, step 4). Reason: two or three LOO samples give an error estimate that can be near zero by chance, and [Lo22] §4.1.2 measured the transfer-function rms over a night at 0.015–0.08 in visibility.

**10. Transfer-function modes: one per BCD position, or shared across positions?**
- Per position (`tf-vis2-base13-24-OO@<night>`): simple, but four times as many modes, and it misses the atmospheric part common to all positions.
- Shared (one name, from the concatenated split of Stage 5): fewer modes, correlated across positions, but it relies on a split from few calibrators.
- **Decision:** Per position in Stage 4; shared as an option once Stage 5 shows a common component.

**11. Should MATISSE mode and row names use sorted station order?**
- Sorted: stable across BCD positions and calibrators, but differs from GRAVITY's file order.
- File order: the same as GRAVITY, but one baseline gets two names across positions.
- **Decision:** Sorted for MATISSE, through the `canonical_orientation` fact; GRAVITY unchanged.

**12. Closure-phase TF modes that span frames: wait for virgil, or split them per frame?**
- Wait for virgil (§8, issue 5): correct, and it also fixes GRAVITY's `tf-t3phi` modes.
- Split per frame (`--split-t3-modes`): usable now, but it drops the correlation between frames and so understates the error of an average over frames.
- **Decision:** Open issue 5 now, and offer the split as an explicit option until it lands.
- Not revised, but the split writes a `HISTORY` card and a warning that errors of averages over frames are understated.

**13. Where do calibrator diameters come from if the `CALIB_RAW_INT` products do not record them?**
- `--diameter` only: explicit, but tedious for a night of calibrators.
- The `JSDC_CAT` file that CalSelector returns and the pipeline uses: automatic, and the same as the pipeline, which parity needs. *Corrected 2026-10-07: this said `mat_raw_estimates` uses it; `JSDC_CAT` is a required input of `mat_est_tf` ([PM] §9.21), and `mat_raw_estimates` lost its diameter options in 1.4.4 ([PM] §2.15).*
- **Decision:** Read `JSDC_CAT` when present, let `--diameter` override it, and refuse otherwise.
- *Revised 2026-10-07 (conservative review):* the diameter read is the band-matched uniform-disc diameter and its error, and the column name, star match (name and coordinates) and value are written to `index.csv` and a `HISTORY` card. A star matched ambiguously, or without a band-matched column, is refused. The diameter error always enters as a diameter mode. Reason: [PM] documents `JSDC_CAT` as an input but not which column the pipeline uses, and no diameter card was found in the product description ([PM] §7; the QC cards carry JSDC magnitudes and fluxes only), so a silent column choice cannot be checked.

### 7.1 Unverified, for a MATISSE expert

Claims the plan relies on that the review could not confirm from a primary source. Each says what to check; until checked, the plan's default is the conservative one in brackets.

1. **Which header key is which BCD.** The [PM] §11.9 formula implies the first word of a position name is the IP1↔IP3 device; [MT] `mat_removeBCD` treats `ESO INS BCD1 = IN` as the IP5↔IP7 swap (§3.1). Neither manual names the keys. Which is right, and does `ESO CFG BCD MODE` read in the same order? [Refused for BCD combination until the Stage 2 fact `bcd_key_order` is recorded.]
2. **The sign pattern for the IP5–IP7 differential phase.** [PM] §11.9 gives (OO + OI − IO − II)/4 for both correctable differential phases; by the swap, IP5–IP7 should reverse under the other device (OO + IO − OI − II). Is the manual's text generic, or does the pipeline reorder first? [Stage 8 only; follow §3.1.]
3. **Closure-phase formulas for the other three triangles.** [PM] §11.9 writes only IP1-IP3-IP5. The plan derives the others (§3.1). [Synthetic test `test_peak_tables_match_section_3_1`.]
4. **Whether `STA_INDEX` in `RAW_INT` is in physical or peak order per BCD position** (`sta_index_frame`, §3.1). No manual states it. [Stage 2 explore job decides.]
5. **Chopping keywords.** That `ESO ISS CHOP ST` marks a chopped exposure and `ESO DET CHOP ST` does not. This rests on the ν Hor headers only; neither manual names either key. [ESO-news] (October 2024) reports N-band data of 2024-05-16 to 2024-08-12 whose frames were tagged with the wrong chopping state; the affected key is not named. [Refuse an exposure where the two keys disagree unless forced? Currently the ISS key is read; to decide with the MATISSE expert.]
6. **Chopped and non-chopped closure phases in one group.** [PM] §4.1 says chopping is unnecessary to measure closure phases; is the instrumental closure phase the same with and without chopping, at one BCD position? [Separate groups (§3.3) until D1's calibrators show no difference.]
7. **The N-band phase-sign version range.** [PM] §2.5 and [GR22] Methods §4.3 say "as of 1.5.2" (GR22 §4.8 shows the 180° flip on θ¹ Sco), 2.0.0 fixed it; the 1.5.2 notes do not mention it, and Varga et al. 2024 (arXiv:2401.03437, Fig. A.3) names 1.5.8, probably a typo. [MT] applies no version-dependent correction. Were N-band phases from versions before 1.5.2 (2018–2020 data) standard? And does `VISPHI` flip with `T3PHI` ([PM] §2.5 names phases and closure phases)? [Before 1.5.2: refused as `"unknown"`.]
8. **`replaceTel`.** [PM] §9.1.4 defines a two-digit code (UT, AT) with tabulated default 3, which would replace AT3's photometry by default. Is that the real default in the ESO reductions, why AT3, and is "telescope k" the `ISS CONF` telescope or the input? [Warn; extra mode or refusal, Stage 5.]
9. **`bcdMode=5`** in the ν Hor headers: [PM] §9.13 documents 0–3 for `mat_merge_results`. What does 5 mean, and which tool wrote it? [Merged products refused anyway.]
10. **Units and keys of `corrFlux=TRUE` products.** `VIS2DATA` in ADU² (not stated in [PM] §7), and `AMPTYP` against the manual's `AMPTYPE`. [Reader accepts both spellings and refuses neither.]
11. **Spatial-filter header keys**, and whether 1.5λ/D (LM) and 2λ/D (N) are the service-mode defaults ([UM] Table 2 marks the default by highlighting, which the text extraction lost). [Read from the header; refuse for primary-beam use if absent.]
12. **Resolving powers.** [Lo22] Table 2 and [UM] Table 3 disagree (LM LOW 31.5 against 34; L HIGH 979 against 959; HIGH+ 3370, 3666 or 3300). Which applies to the line-spread function, at which wavelength? [Both stored; R fitted where it matters.]
13. **Table 6 photometric errors** ([Lo22]): 1σ per 1-minute exposure, per telescope? Are they valid for GRA4MAT and post-2022 data? [Used as a lower bound on the width only.]
14. **The coherent-flux bias limit on UTs and with GRA4MAT** ([Lo22] §4.2, §6; [ESO-news] August 2020 gives 6–8 Jy). [AT limit applied everywhere.]
15. **Which JSDC column the pipeline uses** for calibrator diameters in `mat_est_tf`, and where the diameter is recorded in `CALIB_RAW_INT` products. [Band-matched UD diameter, logged; refuse if ambiguous.]
16. **The pipeline's averaging in `mat_cal_oifits`**: equal weights across BCD positions for V² and fluxes ([PM] §11.9 says "averaged" without weights), and the weighting of `tfInterp=0`. [Parity gate with sign check, decision 3.]
17. **The Rayleigh–Jeans fallback's 10% width** for K and M giant calibrators near the SiO feature at 8 µm. [Opt-in only, with a warning.]
18. **MATISSE-py** (Soulain et al., arXiv:2609.31283) reports BCD and GRA4MAT corrections and an error budget; read only as a summary here, so none of its numbers are used. Should the plan follow it?
19. **The sequence of chopped L exposures.** [PM] Table 4.1 has 4 chopped at `IN-IN` and 4 at `OUT-OUT`; [Lo22] §3.1 describes "2 × 4" chopped exposures "for each BCD configuration". The synthetic template follows [PM]. [Stage 2 counts files per position.]
20. **Whether the [Va21] App. B floors are absolute or relative** (0.03 in L visibility; 8% in N correlated flux is stated as a percentage). [The visibility floor is treated as absolute in V.]

## 8. virgil issues to open (drafts)

These are not opened yet. Each goes to `benjaminpope/virgil` when the stage that needs it starts, except issue 5, opened as benjaminpope/virgil#297 (§7, decision 12).

**Issue 1. "Bandwidth smearing with an explicit per-dataset spectral response"**

> Spectro-interferometric data with oversampled spectra (MATISSE LOW: 64 channels at R ≈ 31.5, averaged by the pipeline over a 5- or 7-pixel sliding bin) need smearing computed from the real response. `EFF_BAND` describes the pixel spacing, not the response.
>
> Proposal: `OIData.with_spectral_response(response)`. For each `INSNAME` and channel, `response` gives (a) a coherent line-spread kernel (a Gaussian of FWHM λ/R, or a tabulated kernel), applied to the complex visibility, and (b) an incoherent kernel (a top-hat over the sliding bin), applied to V² and correlated flux after the modulus is taken. The model visibility is integrated over (a), squared, then averaged over (b). virgil-vlti supplies the response from its instrument facts.
>
> Acceptance: a binary at 30 mas on a 130 m baseline in L loses about 5% of its contrast (the MATISSE review's estimate), matching direct numerical integration. See `design/spectro_interferometry_workflow.md` §2.8 ("6a PR C").

**Issue 2. "Time smearing for long or merged integrations"**

> `read_oifits` uses `INT_TIME` only for matching. Merged products and long exposures (MATISSE: 60–120 s) move u and v by Earth rotation during the integration.
>
> Proposal: read `INT_TIME`, `OI_ARRAY.STAXYZ` and the target's RA and Dec (`OI_TARGET`), and sub-sample u, v over the integration, with the number of sub-samples set by the rotation rate and the model's extent. Average the model's V² or correlated flux incoherently. Off by default; on when the u, v motion over `INT_TIME` exceeds a tolerance.
>
> Acceptance: a binary's V² matches an analytic average over its track.

**Issue 3. "Pinhole primary beam beside the fibre beam"**

> MATISSE uses spatial filters (standard pinholes of 1.5λ/D in L/M and 2λ/D in N, [Lo22] Table 2; other pinholes and N-band slits exist, [UM] Table 2), not single-mode fibres. Its field transmission P(θ) is the fraction of the telescope PSF, offset by θ, that passes a circular hole, not a Gaussian fibre mode.
>
> Proposal: a primary-beam option `kind="pinhole"` with a diameter in λ/D, giving P(θ) as the energy of an offset Airy pattern inside the hole, and the coupled flux it implies for `OI_FLUX`. virgil-vlti supplies the facts per band and telescope.
>
> Acceptance: P(0) equals the Airy encircled energy at the hole's radius, and the half-power offset matches a direct numerical integration. See `design/spectro_interferometry_workflow.md` §2.2.

**Issue 4. "Banded channel covariance for sliding-average spectra"**

> A sliding average over w pixels correlates channel j with j ± k by (w − k)/w for k < w. Today the rule is to bin to one resolution element.
>
> Proposal: a covariance option per `INSNAME`, C = C_obs ⊗ C_chan, with C_chan banded (Toeplitz, correlation (w − k)/w at lag k) and C_obs the existing per-observable covariance (the whitened closure phases for T3). It is used in the likelihood without binning. `OI_CORR` is read when present.
>
> Acceptance: on synthetic data, the unbinned likelihood with the banded covariance gives the same parameter errors as the binned likelihood, to within the information that binning loses.

**Issue 5. "Closure-phase modes spanning several frames"**

> `OIData.with_closure_offsets(..., modes=...)` raises "A closure-phase mode spans several frames", while `with_gains(modes=...)` accepts modes across frames. virgil-vlti writes per-night closure-phase transfer-function modes (`tf-t3phi-tri<a>-<b>-<c>@<night>`) for GRAVITY now, and will write per-position ones for MATISSE; both span several frames.
>
> Proposal: allow modes of shape (n_mode, n_phase_total) across frames, applied to the whitened closure phases through the same projection as single-frame modes.
>
> Acceptance: a mode shared by two frames reproduces the covariance of an injected common offset, and the single-frame case is unchanged.

## 9. References

- **[PM]** ESO, *MATISSE Pipeline User Manual*, VLT-MAN-MAT-15860-9310, issue 2.5.0, 2026-08-06 (ftp.eso.org/pub/dfs/pipelines/instruments/matisse/). Sections cited: §2.1 (2.5.0 notes: `spectralAverage`, `cpAlgo`, T3AMP from VIS2, chromatic OPD), §2.5 (2.0.0 notes: N-band phase sign), §2.8 (1.7.6), §2.11 (1.5.2), §2.12 (1.5.1: windows), §2.14 (1.4.9: binning advice), §2.15 (1.4.4), §4.1–4.2 and Table 4.1 (BCD, chopping, sequence), §7 (product columns and units), §9.1.4 (`mat_raw_estimates` options, `replaceTel`), §9.5 (`mat_cal_oifits` options; no correlated-flux calibration), §9.13 (`mat_merge_results`, `bcdMode`), §9.21 (`mat_est_tf`, `JSDC_CAT`), §11.7–11.9 (GRA4MAT; BCD calibration). The section numbers of the release notes shift between issues.
- **[UM]** ESO, *MATISSE User Manual*, VLT-MAN-ESO-289268, issue 1.0, 2026-07-03 (P118). §1, §2.2.1 (BCD), §2.2.2 and Table 2 (spatial filters), Table 3 (resolutions), §2.3–2.4, §2.7, §3.2.
- **[Lo22]** Lopez, B. et al. 2022, A&A 659, A192, *MATISSE, the VLTI mid-infrared imaging spectro-interferometer*, arXiv:2110.15556 (read in the arXiv version). §2.4–2.5 (BCD, spatial filtering), §3.1, §3.3, §4.1.2–4.1.4, §4.2, Fig. 14, Fig. 18, §6, Tables 2 and 4–6, App. A eqs. A.8–A.9.
- **[ESO-news]** ESO, MATISSE news page, https://www.eso.org/sci/facilities/paranal/instruments/matisse/news.html, read 2026-10-07: pipeline 2.0.0 and the N-band phase sign (December 2023); N-band chopping-state tagging (October 2024); N-band piston-correction bias (August 2020).
- **[Va21]** Varga, J. et al. 2021, A&A 647, A56, arXiv:2012.05697 (HD 163296). §3.1 (DRS 1.5.0; correlated flux in N; eqs. 2–3), App. B and Table B.1 (floors as conservative lower limits; diameter and spectrum errors excluded).
- **[GR22]** Gámez Rosas, V. et al. 2022, Nature 602, 403, arXiv:2112.13694 (NGC 1068). Methods §4.3 (pipeline versions, N-band flip), §4.5 (BCD positions calibrated separately), §4.8 (image orientation).
- **[Is22]** Isbell, J. W. et al. 2022, A&A 663, A35, arXiv:2205.01575 (Circinus). §2.2 (DRS 1.5.1; calibrator closure phase subtracted, then the four positions averaged, their standard deviation the uncertainty).
- **[Mi16]** Millour, F. et al. 2016, Proc. SPIE 9907, arXiv:1608.01913.
- **[Wo24]** Woillez, J. et al. 2024, *GRAVITY for MATISSE*, arXiv:2405.20730 (GRA4MAT; no calibration floors quoted in the abstract).
- **[MT]** MATISSE consortium tools, https://github.com/Matisse-Consortium/tools, read at b44f965889ab19a810260d1670b19028c9d49f62 (2026-05-22): `mat_autoCalib.py` (`findClosestCal`), `libPostTools.py` (`mat_removeBCD`, `mat_mergeOifits`), `mat_mergeAllOiFits.py`, `libFluxCal_STARSFLUX_new.py`. Code, not a publication.
- Lykou, F. et al. 2023, arXiv:2311.07250, and Varga, J. et al. 2024, arXiv:2401.03437: independent reports of the N-band closure-phase flip.
