# Empirical nulls as the false-alarm reference

Status: **design, nothing implemented** (a background working note).

A companion search reports a maximum Δχ² over a position grid. Its significance
depends on what the data's errors leave out: correlated chromatic and
telescope-based calibration systematics (`gravity_calibrator_pca.md`), which
make χ²_red 2–10 on real GRAVITY data and break any threshold taken from the
reported σ. A calibrator has no companion to find, so a calibrator calibrated as
a science target is data of the science kind whose truth is "nothing there"
beyond the star itself: a uniform disk of the catalogue diameter (see below).

## What `virgil_vlti.nulls` provides

* `null_datasets(exposures, out_dir, ...)`: each calibrator exposure calibrated
  by the other stars of its night (same setup, stations and field), with the
  science calibration defaults (linear time weights, leave-one-out transfer-function
  error: `calibrate_exposure`'s, so the error bars are as wide as a science target's;
  pass `error=`/`method=` to match a science run), while `calibrator_residuals`
  (and so `resampled_null_datasets`) still uses Gaussian weights and propagated
  errors: `null_datasets` is the reference for false-alarm thresholds,
  written as OIFITS by `write_calibrated`: `OI_VIS2`, `OI_T3`, `OI_VIS`
  (amplitude and phase) and the observed `OI_FLUX`, at full resolution. The real
  uv, wavelength grid, flags and error bars are kept, so the files read like
  science files. Options: `exclude`, `resolution`, `polarisation`, `telescope`,
  `field`, `min_calibrators_per_night` (other stars; default 2). Possible
  binaries are flagged, first by `screen_calibrators`, then by the pair test on
  each star's calibrated exposures (`white` width above `pair_factor` times the
  *statistical* σ: the test runs on a second calibration with `error="none"`,
  because consecutive exposures share their references and the propagated
  transfer-function error cancels in their difference). A flagged star is kept
  out of the reference sets only: its own null files are still written, with
  `flagged: true` and `flag_reasons` in `null_manifest.json` and in
  `NullDatasets.dropped`. **Why:** a science target is never screened on its
  own scatter, so removing flagged stars as targets would measure the max-Δχ²
  distribution on the quiet subset and bias the tail, and the false-alarm
  rate, low. Compute the false-alarm probability with and without the flagged
  nulls and report both. Each manifest entry records the catalogue diameter
  (`diameter_mas`, `diameter_err_mas`, from `ESO QC TF CAL_SC DIAM`, the
  diameter the pipeline calibrated with).
* `block_resample(residuals, n, block="exposure"|"chunk", chunk_channels=...)`:
  more realisations of the residuals, from whole exposures or channel chunks of
  exposures redrawn among those of the same row layout (row counts, station
  pair or triple of every row, and field); it warns when an exposure has no
  other with its layout. It works on residual
  arrays (`log V²`, closure phases), not on files.
* `resampled_null_datasets(residuals, n, out_dir, ...)`: writes `block_resample`
  realisations as OIFITS files of the layout `null_datasets` writes (the slot
  exposure's uv, times and flux; V² its uniform disk × exp(r); closure phases
  the donor's; no `T3AMP` or `OI_VIS`), for nights with too few calibrators for
  `null_datasets`. The residuals must come from `calibrator_residuals` with the
  flagged stars passed as `exclude=` (it does no screening itself); the same
  stars go to `exclude=` here, which refuses residuals that still hold them. The
  manifest records block, seed, donors and `donor_pool` (how many exposures
  each slot can be drawn from).
* `null_datasets(reference_subsets=k, subset_size=m)`: up to `k` nulls per
  exposure from different `m`-star subsets of its night's other stars.

## How to use them

0. **The null's truth is a uniform disk, not a point.** Each null is the
   calibrator minus the other calibrators' transfer function, so what it
   contains is the calibrator itself: a uniform disk of the catalogue diameter
   (`diameter_mas` in the manifest). The search must handle that primary the
   way the science fit does (a fitted or fixed primary diameter, with its
   uncertainty), or the null's long-baseline V² deficit will be read as
   signal. Injections go on top of that disk.
1. **One null class per configuration** (MED or HIGH, AT or UT, single or dual
   field): run `null_datasets` with the filters, about 50–100 nulls per class.
2. **Run the full search**, over the same grid and with the same model and
   nuisances as the science fit, on every null. Record the maximum Δχ² of each.
3. **False-alarm probability**: the distribution of that maximum, including the
   look-elsewhere effect. Fit a GEV to the tail for the 3σ and 5σ thresholds,
   rather than extrapolating the few largest values.
4. **Injection–recovery and limits**: add a model companion to each null's
   observables, recover it, and compare with the threshold. Injections must use
   the same `fibre_response` as the science model.

## Cautions

* The nulls are calibrators of *their* nights, so their systematics are those of
  the programme's calibrators; the science star's own calibrators are included
  among them. Calibrating with the same stars as a science target shares
  calibrator noise: keep the number of references similar to the science case
  (`min_calibrators_per_night`).
* Distinct nulls from one night share references, so they are not independent.
  Count effective independent nulls by night; a large sample from few nights
  understates the tail.
* A calibrator's diameter error is in its null as it is in a science target's;
  a bad catalogue diameter would show as excess null signal on long baselines.
* `resampled_null_datasets` files and `reference_subsets` nulls are not
  independent nulls: resampled slots reuse the residuals of their layout pool
  (about as many independent nulls as distinct donors, `donor_pool` in the
  manifest, not `n`), and subsets of one exposure share its data. A GEV fitted
  to 100 of them has a tail far better constrained than it is: take errors
  from the number of independent donors, exposures and nights, not files.
* `block_resample` redraws residuals, not data: use it for bootstrap errors on
  the tail fit, and chunk resampling as a check that the false-alarm rate does
  not hinge on long-range chromatic structure. It does not create new
  systematics; for those, use more nights.
* The pair-test binary screen is a heuristic (threshold 3, untuned). Check the
  dropped list by eye.
