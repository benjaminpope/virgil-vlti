# virgil-vlti: critical review and a general instrument-calibration contract

Status: 2026-10-19. A read-only review of virgil-vlti at `origin/main` 448625a and virgil at `origin/main` e066d7b. Nothing was run on real data. Tests were read, not re-run at these heads. Note that the local `main` checkout was 24 merges behind; every item below was re-checked at the heads above.

Three questions:
1. Is the calibration correct?
2. Can its output flow into virgil workflows in general, with the VLTI knowledge kept here?
3. What should the pattern be for AMI, NRM and CHARA?

## Summary

- **The calibration core is sound.** The maths is careful and well tested by injection and recovery. It agrees with `gravity_viscal` to 1e-11 on four public cases. The leave-one-out transfer-function error, the diameter, TF and template modes, and the refusal of mismatched modes are all good design.
- **There are about six real correctness gaps** (§1). Most matter only in edge cases (one calibrator, reversed stations, NaN errors). One that matters routinely: missing header keys count as a match, so non-GRAVITY data, and any product lacking `vis-correction-sc`, pass the compatibility check silently.
- **The interface with virgil is the weakest part (§2).**
  - virgil-vlti writes three kinds of sidecar `.npz`. None has a schema, units, widths or a version, and **virgil has no reader for any of them**.
  - Users stitch modes into `with_gains(modes=)` by hand, using a docs recipe that assumes virgil's sample order equals the row-major table order. That assumption is untested.
  - virgil's planned `virgil.noise_spec/1` (design/instrument_noise_in_pipelines.md) is a good start. It carries only global widths, though, so it cannot carry the per-datum modes that are this package's main product.
- **Pipeline provenance is thin (§3).** The output does not record which calibrators, diameters, code version or pipeline settings made it. The container is unpinned.
- **The recommendation (§4)** is that calibrated data and its nuisance structure travel together, in one OIFITS file:
  - standard `OI_CORR` for fixed correlations;
  - a small, versioned virgil extension for low-rank modes with named, shared parameters;
  - the same file read by virgil's `read_oifits`, with no instrument code in virgil.
  
  An instrument package becomes "anything that writes a conforming file". A conformance test suite lives in virgil, and every instrument package runs it.

## 1. Correctness: what to fix

These are ranked. "Present" means verified at the current head.

| # | Issue | Where | Fix |
|---|---|---|---|
| 1 | **Missing keys count as a match.** `p.get(k) != q.get(k)` passes when both sides lack the key. So a product with no `vis-correction-sc` recorded, or any non-GRAVITY file under the GRAVITY rules, is "compatible" with no check made. Only `PRO REC1` parameters are read. | gravity.py:137-147 | Treat absent-on-either-side as a named "not recorded" outcome (refuse by default, `--allow-unrecorded`). Locate the `gravity_vis` recipe among all `PRO RECn ID` and confirm on a real header which RECn holds it. |
| 2 | **One calibrator gives no TF error.** `prediction_error` returns 0 for n<2, and the calibrator's own statistical error is not added in `loo`/`scatter`. The TF mode vanishes too. One calibrator per night is common. | calibrate.py:347-375 | Floor at `propagated_error` (quadrature with LOO). For n=1, warn and either refuse or use a per-instrument prior TF scatter (e.g. GRAVITY ~2% V²) stated in the output. |
| 3 | **NaN and ≤0 errors are not flagged** in the generic path (the viscal port flags them). Zero errors become infinite weights downstream. | calibrate.py:131-141 | Add `~isfinite(err) \| (err <= 0)` to `_bad`. |
| 4 | **Baselines and triangles are matched as exact `STA_INDEX` tuples.** Reversed or rotated orderings are dropped as "missing" (safe, but loses data). The same index can mean different stations in two files, since `OI_ARRAY` names are never compared. Configuration (`ISS CONF STATION*`) is not compared between science and calibrator either. | calibrate.py:223, 635 | Canonicalise on station *names*: sort pairs (flip VISPHI sign), cyclic-order triangles (flip T3PHI on odd permutations). This belongs in virgil's reader as a generic utility (§4). |
| 5 | **Bins with half their channels bad pass.** Their error is `√Σe²/n` with no account of the missing channels or of GRAVITY's ~2× spectral oversampling. The calibrators' statistical errors within a bin are also combined as independent. | calibrate.py:593, 710-722 | Require ≥ 75% good channels, and state the oversampling factor per instrument mode in the profile (§4). |
| 6 | **The Gaussian TF weights have no error floor.** They use 1/median σ² with no floor (the port has one), so a tiny-σ calibrator dominates and the weights underflow at long gaps. | calibrate.py:306-320 | Use the pipeline's floors (2%, 0.2°) and warn on underflow. |
| 7 | **The calibrator's absolute VISPHI is assumed zero.** This is wrong for a resolved or offset in-field calibrator. | calibrate.py:229 | Subtract the calibrator model phase (zero for a centred UD, but say so), or refuse VISPHI from off-axis calibrators. |
| 8 | **Nulls are cleaner than science.** `nulls.py` defaults to `error="propagated"` while science uses `loo`, and it ignores reference diameter errors. The null Δχ² distribution is therefore narrower than the science one, and false-alarm thresholds come out too optimistic. | nulls.py:127, 215, 339 | Use the science run's exact error settings and modes, ideally by reading them from the science output's provenance (§3). |
| 9 | **Three frame definitions** exist: `frame_times`, `frames_of` (60 s gap, silently drops a frame holding a baseline twice), and `calpca._rows_meta`. virgil has a fourth. | calibrate.py:1114, pairs.py:85-111, calpca.py:258 | One rule, in virgil (frames are generic). virgil-vlti uses it and warns on any drop. |

The statistical items are lower priority but worth an issue each:
- **calpca widths.** Residuals of calibrator *j* calibrated by *k* carry *k*'s TF error, so widths are inflated by up to √2, and rows within a file are correlated (calpca.py, `error="propagated"`).
- **Telescope fraction.** `telescope_fraction` is not noise-debiased and is pulled toward 2/3 (calpca.py:540-566).
- **log-V² residuals** are biased at low SNR.
- **Pair-test white factor.** It is biased low, because calibrated errors include τ, which cancels in pair differences (pairs.py). The noise-spec design's `sigma_convention` is the right fix, but it should land before the spec is emitted.
- **Templates are per frame,** while the residuals have a per-exposure common part (templates.py:287).
- **Template mode names** use the file *basename*, so two nights with equal basenames silently share a parameter (templates.py:~287). Use a hash of the path or of `ARCFILE`.
- **Paranal-only assumptions.** `Exposure.night` uses a fixed 12 UT boundary, and `read_exposure` requires `OI_T3`. Both break for 2-telescope data and for other sites.

Test gaps match these items:
- reversed stations;
- n=1 calibrator;
- NaN/zero errors;
- frames < 60 s apart;
- multi-calibrator, LOW-resolution, T3AMP and `OI_VIS` viscal parity;
- Monte Carlo coverage of the error bars (listed as "requested" in PROVENANCE).

## 2. The interface with virgil

**Today.** virgil-vlti writes four kinds of output:
- calibrated OIFITS;
- `.modes.npz` (diameter);
- `.tf_modes.npz`;
- `.template_modes.npz`.

The modes are keyed `"<insname>/<obs>/<mode>"` and shaped like the table (rows × channels, dense). V² modes are d log|V| per 1σ, which is now consistent. T3PHI modes are in **degrees**, while virgil's `with_closure_offsets(modes=)` wants **radians** (gains.py:567). Diameter modes are absolute unless `relative=True`. None of this is stamped in the file. virgil has:
- no reader;
- no handling of shared names across files (`<cal>@<night>`);
- no consumer for T3AMP or VISAMP modes;
- only one shared width for all external modes, so template modes cannot be kept apart from virgil's built-in telescope and baseline groups (double-counting).

No test in either repo crosses the boundary.

**Planned.** virgil's `virgil.noise_spec/1` is the right shape for *global* noise properties, and its rules are good:
- declarative;
- unit conversion done once, in the emitter;
- unknown keys are errors;
- major-version bumps for new groups;
- `apply_noise_spec` as the single place a spec acts.

But `gains.modes` is explicitly deferred to a v2. The per-datum modes (diameters, TF, PCA templates) are what this package exists to produce, and they are also what AMI (covariance blocks) and CHARA (TF modes) need. They cannot go in a JSON spec, because they have to be aligned to data.

**The alignment problem is the real risk.** Any sidecar must be matched to virgil's flattened sample order. That order depends on:
- `insname` filtering;
- `vis_index`;
- flags;
- windowing and binning;
- table order.

The docs recipe uses `reshape(-1)` and trusts the order. A misalignment is silent: the fit still runs, with nuisance modes attached to the wrong data.

## 3. ESO pipeline interaction

The input level is right: `*_SCI_VIS`/`*_CAL_VIS`, SC only, with already-calibrated products refused. The MATISSE stage 1 refusals (merges, corrFlux as V², N-band sign) are good. The gaps:

- **Output provenance.** The written file keeps the science file's `ESO PRO CATG` and `PRO REC` cards, so it looks like an ESO product. It has no record of:
  - the calibrator files (`ARCFILE`/MD5);
  - the diameters and their errors;
  - the error mode and interpolation;
  - the virgil-vlti version and git SHA.
  
  The HISTORY line hard-codes "pipeline 1.9.7" (cli.py:302-307). Fix: set `PRO CATG` to a virgil-vlti category (or remove it). Add a `VIRGIL_PROVENANCE` binary table, or `HIERARCH VIRGIL` cards, holding all of the above. The nulls (§1 #8) and any re-run can then read it back.
- **Only `vis-correction-sc` is recorded.** Record the whole `gravity_vis` parameter set per exposure, and the "applied by pipeline" state (vFactor, phase reference, `output-phase-sc`). Then a fitted coherence-loss mode can refuse to double-correct.
- **The pipeline version is only warned about.** Make the comparison strict by default, with an explicit allow-list of compatible versions per instrument.
- **The container is unpinned.** `dnf install esopipe-*` takes the latest, and the EDPS failure is swallowed by `|| echo`. The versions manifest is a 30-day CI artifact. Fix: pin the RPM versions, commit `containers/versions.lock` and the SIF digest, and add a CI check. The viscal port follows 1.9.7 source while the container installs 1.11.x; either pin to 1.9.7 for parity runs or re-port.
- **No reduction wrapper.** There is no in-repo `reduce` wrapper that records the sof, parameters and container digest, so the `PRO REC` cards are the only record. (This matters less if we only ever start from archive products, but say so.)
- **viscal parity mode** (`calibrate_like_viscal`) skips `check_compatible`.
- **Remaining MATISSE and PIONIER gaps.** BCD- and chopping-aware calibration is absent, and `TARGET_RAW_INT` is still accepted as science. PIONIER (`pndrsCalibrate` parity) has nothing yet.

## 4. Recommendation: a general instrument-calibration contract

The aim is that virgil knows nothing about instruments, an instrument package knows nothing about virgil's internals, and they meet at a file format with a conformance suite. Four parts follow.

### 4a. The calibrated OIFITS file carries its own noise structure (virgil owns the spec)

One file, so alignment can never drift:

- **Data**: standard OIFITS2 tables, with `VISCAL`/`CALSTAT = 'C'` and a non-ESO `PRO CATG`.
- **Fixed correlations**: standard **`OI_CORR`** plus `CORRINDX` columns. OIFITS2 already defines these, but neither repo reads or writes them today. They are the natural home for:
  - GRAVITY's per-DIT covariances (PLAN item 1);
  - AMI's closure-phase covariance;
  - spectral-channel oversampling correlations.
  
  Other tools (PMOIRED, CANDID) can read them too.
- **Low-rank nuisances**: two virgil extension tables, versioned by `EXTVER` and a `VIRGIL_SCHEMA = 'virgil.modes/1'` card:
  - `VIRGIL_MODE_DEF`: one row per parameter, with columns `NAME`, `GROUP`, `SPACE` (`logV`, `phase_rad`, `t3amp_rel`, ...), `WIDTH`, `SHARED_KEY` and `DESCRIPTION`. `SHARED_KEY` is globally unique, for example `gravity:diam:HD1234:2024-03-05`, or a hash, never a basename. Equal keys across files are one parameter.
  - `VIRGIL_MODES`: a sparse table with columns `MODE` (row in `MODE_DEF`), `EXTNAME`, `INSNAME`, `ROW` and `VALUE[nchan]`. Rows are addressed by the OIFITS table and row index of *this file*, which removes the alignment problem.
- **Global widths** (white scales, untemplated gain widths): the existing `virgil.noise_spec/1`, either as JSON or as a `VIRGIL_NOISE_SPEC` header card holding the JSON.
- **Provenance**: a `VIRGIL_PROVENANCE` table with free-form key and value pairs that virgil never interprets but carries into `run.json`.

virgil's side is:
- `read_oifits` reads `OI_CORR` and the `VIRGIL_*` tables when present. It maps (EXTNAME, ROW, channel) to its flat sample order *inside the reader*, where the order is known, and builds the gains and closure groups with per-group widths.
- A `combine` merges shared keys across files.

Unknown `SPACE` values or schema versions are refused, following the spec's own rule.

Why not more `.npz` sidecars? A FITS extension:
- cannot be separated from its data by a copy or rename;
- survives OIFITS tools that ignore unknown extensions;
- is addressable by row;
- is readable without Python.

`noise_spec/1` stays for the global case, and the `gains.modes` it defers becomes `virgil.modes/1`.

### 4b. Generic machinery moves to virgil, or to a small shared core

Some things in this repo are not VLTI-specific, and AMI, NRM and CHARA would each re-implement them:

| Generic (should be shared) | Stays instrument-specific |
|---|---|
| Station-name canonicalisation of baselines and triangles with sign flips; frame grouping; `OI_CORR` I/O; modes I/O; unit guards | GRAVITY FT/SC, `PRO REC*`, vis-correction, dual-field rules, fibre response, the viscal port |
| TF estimate, interpolation (linear/Gaussian), leave-one-out error, UD calibrator models, diameter modes, TF modes, binning with correlated errors | MATISSE BCD/chopping, corrFlux, N-band sign |
| Pair test core (`closure_matrix`, `fit_noise`), weighted PCA with held-out widths, nulls by block resampling | Which header keys define a mode; product classification; site; oversampling factor; known pipeline-version hazards |

The I/O row (canonicalisation, frames, `OI_CORR`, modes) belongs in **virgil**, since it is reading and likelihood machinery, as PLAN.md already places it. The calibration row is generic only for *calibrator-star* instruments (VLTI, CHARA, ground SAM via AMICAL's own calibration). AMI calibrates by forward modelling (AMIGO), so it would only use the I/O row. **Do not extract a separate `virgil-calib` package yet.** Instead, make the generic modules in this repo free of GRAVITY imports, behind an instrument profile (4c). Extract them when the second calibrator-star instrument (CHARA/MIRC-X) actually arrives. That is the rule of three, at the cost of one refactor.

### 4c. An `InstrumentProfile` in each instrument package

This continues what `INSTRUMENT_MODULES` (calibrate.py:1054) has started. Each instrument is a small module that declares data, with no `if gravity:` logic in the generic paths:

```python
@dataclass(frozen=True)
class InstrumentProfile:
    name: str  # "GRAVITY"
    site_longitude_deg: float  # night boundary, not "12 UT"
    science_insnames: Callable  # SC only; FT refused
    mode_keys: tuple[str, ...]  # must match between sci and cal; absent = refuse
    reduction_params: Callable  # header -> {name: value}, finds its own recipe
    version_key: str
    compatible_versions: Callable
    classify: Callable  # header -> SCI/CAL/CALIBRATED/REFUSED + reason
    channel_oversampling: Callable  # mode -> factor, for binning
    screen: Callable  # instrument refusals (corrFlux, merges, sign)
    provenance: Callable  # header -> dict for VIRGIL_PROVENANCE
```

The GRAVITY, MATISSE and PIONIER profiles stay here. A CHARA package (MIRC-X/MYSTIC pipeline products) would add one profile and reuse the generic calibration. AMICAL needs no profile: it writes the conforming file itself, as a small AMICAL PR adding `OI_CORR` and `VIRGIL_MODES` for its PSF-calibrator errors. AMI/AMIGO writes `OI_CORR` blocks.

### 4d. Conformance and validation

- **virgil ships `virgil.testing.conformance`**:
  - a synthetic file, written by a reference writer;
  - checks that `read_oifits` aligns every mode to the right sample;
  - checks for the units (deliberate degree-vs-radian traps), for the sign of reversed baselines and triangles, and for shared keys across two files;
  - a check that `apply` equals the manual `with_gains` chain.
- **Every instrument package runs the conformance suite in CI** against its own writer. This is the round-trip test that does not exist today.
- **Validation stays as now.** The pipeline-parity oracle (the viscal port, later `mat_cal_oifits` and `pndrsCalibrate`) remains separate from the default calibration. Independent checks go to virgil-vlti-validation, plus a Monte Carlo coverage test of the calibrated error bars.

## 5. Suggested order

1. **Correctness PR (virgil-vlti):** §1 #1-#3 and #8, with tests. Small; it can happen now.
2. **Provenance PR (virgil-vlti):** the output category, a `VIRGIL_PROVENANCE` table, full recipe parameters, and container pinning with `versions.lock`.
3. **virgil design note and PR:** `virgil.modes/1` (§4a), `OI_CORR` reading, station-name canonicalisation, one frame rule, and the conformance suite. This should be decided alongside `noise_spec/1`'s PR 1 so that the two schemas share their rules.
4. **virgil-vlti writer:** emit `VIRGIL_MODE_DEF`/`VIRGIL_MODES` into the calibrated file, deprecate the `.npz` sidecars, and run the conformance suite.
5. **Profile refactor (virgil-vlti):** move the GRAVITY keys out of the generic paths behind `InstrumentProfile`. Fix the night boundary and the `OI_T3` requirement.
6. Statistical items (§1 tail) as separate issues; the calpca width bias and the pair-test τ come first, because the noise spec will publish their numbers.

## 6. GRAVITY observing modes and the fibre beam

Added after the first pass, at the same heads. Facts about GRAVITY taken from memory, not from the manual or a product, are marked [unverified].

**Headline.** The fibre beam exists only as tested helper functions (`gravity.fibre_coupling` and `fibre_response`, gravity.py:238-317). Nothing in the calibration or the CLI calls them, and virgil has no fibre or primary beam for `OIData` (`virgil.radio.primary_beam` is a radio-dish beam for imaging, and the generic beam is "PR C" in virgil's imaging_plan.md, not done). Mode handling is an offset classifier: single field vs dual field, on-axis vs off-axis vs Wide. Phase referencing and astrometry are absent.

| Mode / effect | Status | Evidence |
|---|---|---|
| Single field (SF): SC only, FT never read | handled | gravity.py:5-6, 25 |
| PHITYP / PHIORDER / VISREFMAP consistent between science and calibrator | handled (refuses mismatch and VISREFMAP) | calibrate.py:736-762; virgil oifits.py:825 |
| Differential VISPHI | partly: virgil's `DifferentialPhase` removes a polynomial continuum | review footgun 6 |
| Dual field (DF) on-axis: FT star as in-field calibrator | handled, but see gap 3 | gravity.py:40, products.py:203-205 |
| DF off-axis and GRAVITY Wide | classified only; used for grouping in PCA, pairs, nulls and templates, not calibrated differently | gravity.py:54-84 |
| Anisoplanatic V² loss | warning only (2 %) | gravity.py:176-216, cli.py:442-451 |
| Phase-referenced (absolute) DF VISPHI relative to the FT star | not modelled | no code |
| Metrology, swaps, dOPD astrometry | absent | — |
| Calibrator matching in DF off-axis (separate star) | partly: no check of field mode, FT magnitude or FT mode, or fibre offset | gravity.py:107-137 |
| Fibre profile | Gaussian, FWHM ∝ λ/D (chromatic), jitter-averaged, power or amplitude; not checked against the true coupled mode | gravity.py:238-317 |
| Per-telescope pointing (differential flexure, field and pupil rotation) | absent: one scalar offset per call, isotropic jitter | gravity.py:302-317 |
| Fibre offset that changes within an exposure (swaps, re-pointing) | first header value only | `fibre_offset_mas` |
| Polarisation | compared as a mode key; P1 and P2 calibrated separately | gravity.py:3-5 |

**The physics the beam has to carry.** The beam belongs in the model, not in the calibration. The transfer function from an on-axis point-source calibrator is the right correction for the on-axis visibility, whatever the science field looks like. A companion at offset **s** is then seen by telescope *k* through an amplitude coupling A_k(**s**, λ). The correctly normalised visibility on baseline *kl* is

V_kl = Σ_i f_i A_k(s_i) A_l*(s_i) e^{−2πi **b**·**s**_i/λ} / √(Σ_i f_i |A_k(s_i)|² · Σ_i f_i |A_l(s_i)|²).

For identical, perfectly pointed telescopes this reduces to a beam-weighted image: each source's flux is attenuated by |A(s_i)|², so a fitted contrast comes out low by that factor, and the closure phases are those of the attenuated image. With per-telescope pointing errors, the attenuation becomes baseline-dependent and the closure phases shift. That is generic: CHARA/MIRC-X (single-mode fibres), AMI (none) and an NRM pupil mask all fit the same per-telescope A_k form. So:

- **virgil (generic):** an optional per-telescope, per-wavelength amplitude pattern hook in the `OIData` forward model, using the formula above. It should also carry a phase-centre offset term for absolute, phase-referenced visibilities: S(λ)·exp(−2πi(Δα u + Δδ v)/λ), fitted, as in review §273. Keep refusing differential phases as absolute.
- **virgil-vlti (GRAVITY facts):** the FWHM per telescope type and band, the profile (Gaussian vs Airy overlap, [unverified] which is closer), jitter per telescope, the header keys for fibre offsets (`ESO INS SOBJ X/Y` now, plus the FT and SC pointing keys still to identify), and a per-exposure (or per-frame) offset array. It should hand these over as a beam description in the `VIRGIL_PROVENANCE`/modes file (§4a), never as code virgil imports.

**Ranked gaps**
1. **No beam reaches any fit.** Off-centre companions are fitted as if the coupling were 1, which biases contrast low by |A|². At FWHM/2 that is a factor of 2 in flux [the Gaussian model]. Fix: the virgil hook above, and virgil-vlti writing the offset and beam parameters per exposure. Until then, the CLI should warn whenever a science target is in DF with a separation above about 0.1 FWHM, and docs should say that contrasts are not beam-corrected.
2. **Absolute DF phases are not modelled.** The phase-referenced VISPHI carries the astrometry, and the code neither uses it nor says so. Fix: the virgil phase-centre term, and a stated refusal in virgil-vlti (with a docs line) until metrology and swaps are handled. Astrometry stays out of scope for now, and the docs should say so.
3. **DF classification by offset alone.** With `--dual-field`, *any* exposure under 50 mas counts as a calibrator (products.py:203-205). An on-axis science pointing, or a science target with a companion inside 50 mas, would be mislabelled. Fix: require agreement with `PRO CATG` or a name list, and warn when a named science target falls under the threshold.
4. **Calibrator compatibility ignores field geometry.** Add to the GRAVITY profile's `mode_keys`: field mode (SF/DF), FT mode, FT magnitude within a tolerance, and fibre separation within a tolerance. Refuse by default, following the §1 #1 rule.
5. **Anisoplanatism is warned about, not modelled.** Emit it as a chromatic gain mode with a width (from the GRAVITY Wide model's uncertainty), measured from the phase centre (review §7c.2), so it enters the fit rather than a log.
6. **Per-telescope pointing and jitter.** Make the fibre offsets and jitter per-telescope arrays, read from the header where they are recorded [unverified which keys], and say which effects (field rotation and differential pointing across the night) are assumed negligible.
7. **The beam profile is unchecked.** Validate the Gaussian against the coupled-mode overlap (open Q12, gravity.py:309), or against on-sky offset scans if the GRAVITY team has them (add to design/questions_for_gravity_team.md). Refuse beyond `FIBRE_OFFSET_LIMIT`.

This adds to the §5 order: the virgil beam and phase-centre hook join step 3 (the virgil design note), and items 3–5 join step 5 (the GRAVITY profile).
