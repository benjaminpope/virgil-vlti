<!-- AUTO-GENERATED FROM README.md by scripts/sync_docs.py. Edit README.md, not this file. -->
# virgil-vlti

Calibration and error models for VLTI interferometry, for use with [virgil](https://github.com/benjaminpope/virgil).

What it does for each VLTI instrument:

| Instrument | What virgil-vlti does | Status | Details |
|---|---|---|---|
| GRAVITY | Calibrates V², T3AMP and T3PHI (and `OI_VIS`) from calibrator exposures, with nuisance modes; ports ESO's `gravity_viscal`; calibrator-PCA templates | Calibration checked against ESO's calibrated products | [GRAVITY](instruments/gravity.md) |
| PIONIER | Measures the error structure of calibrated products (the pair test) and turns it into gain templates | Error model only; no calibration | [PIONIER](instruments/pionier.md) |
| MATISSE | Reads headers, refuses unsafe products, reports them (`virgil-vlti-inspect`) | Headers only; calibration planned | [MATISSE](instruments/matisse.md) |

Each page lists what is implemented, what is not, and what each method has been checked against (ESO pipeline code or publications).

virgil fits and images interferometric data from any instrument. This repository holds what depends on one instrument: calibrating pipeline products, empirical error models, and the instrument facts behind them. See [PLAN.md](PLAN.md) for the boundary with virgil and what comes next, [PROVENANCE.md](PROVENANCE.md) for where each method comes from and how it is checked, and [`design/`](design/gravity_calibration_review.md) for the GRAVITY calibration review and notes on PIONIER and MATISSE.

## Install

```bash
pip install virgil-vlti             # numpy, scipy, astropy: calibration only
pip install "virgil-vlti[virgil]"   # and virgil-astro, to fit the output with virgil
```

For development, clone the repository and install it editable:

```bash
git clone https://github.com/benjaminpope/virgil-vlti
cd virgil-vlti
pip install -e ".[dev]"
```

## Quickstart on your own data

**Input.** The *uncalibrated* GRAVITY products from ESO's pipeline (the `*_VIS` files), science and calibrator exposures together in one folder. Download them from the ESO archive; nothing else is needed. Nights are found from the headers.

**One command**, here with separate calibrator stars (see below for the dual-field and other options):

```bash
virgil-vlti-calibrate reduced/ --out calibrated --cal "HD 12345=0.45" --tf-modes
```

**Output**, in `calibrated/`:

```text
calibrated/
  index.csv                         one row per exposure: what was done, warnings, diagnostics
  <night>/<target>_<time>.fits      calibrated OIFITS (V², T3AMP, T3PHI)
  <night>/<target>_<time>.modes.npz           calibrator-diameter nuisance modes (with --cal NAME=DIAM±ERR)
  <night>/<target>_<time>.tf_modes.npz        transfer-function nuisance modes (with --tf-modes)
  <night>/<target>_<time>.template_modes.npz  calibration-template modes (with --template)
```

The transfer function (TF) is the instrument's response measured on calibrator stars, which are close to point sources, and divided out of the science data. Its uncertainty is estimated by leave-one-out: each calibrator is predicted from the others, and the prediction errors give the error of the TF at the science exposure's time. The `.npz` files hold the part of that error (and of the calibrators' diameter errors) that is shared between exposures, as nuisance modes for virgil to marginalise rather than as independent errors.

**Fit with virgil**, which attaches the modes for you:

```python
import glob
from virgil_vlti import to_oidata

data = to_oidata(
    sorted(glob.glob("calibrated/*/*.fits")), "GRAVITY_SC"
)  # a virgil OIData with its modes
```

`to_oidata` needs `virgil-astro>=0.3` and, until 0.4.0 is released, virgil's main branch (`pip install git+https://github.com/benjaminpope/virgil`). Its modes' units and conventions are in its docstring ([`load.py`](https://github.com/benjaminpope/virgil-vlti/blob/main/src/virgil_vlti/load.py)).

## Calibrating products

Checked against ESO's own calibrated GRAVITY products (four public cases): the V² level agrees to 0.12% (median), the errors exactly, and the closure phases to 0.06°. `--match-gravity-viscal` is a port of `gravity_viscal` itself, written from the pipeline's source, and reproduces ESO's products to about 1e-12. It follows the source where the manual differs (e.g. the time weights). See [PROVENANCE.md](PROVENANCE.md).

`virgil-vlti-calibrate` calibrates a folder of uncalibrated products night by night. Calibrator exposures give a transfer function per baseline and channel. It is interpolated in time to each science frame and divided out of V², and the calibrator closure phases (averaged as unit vectors) are subtracted. The interpolation's uncertainty is added to the errors.

GRAVITY in dual-field mode, using the fringe-tracking star as an in-field calibrator:

```bash
virgil-vlti-calibrate reduced/ --out calibrated --dual-field --window 2.21 2.40 --bin 4
```

With separate calibrator stars, optionally with uniform-disk diameters in mas:

```bash
virgil-vlti-calibrate reduced/ --out calibrated --cal "HD 12345=0.45" --cal "HD 678"
```

Then in virgil:

```python
import glob
from virgil import OIData
from virgil.oifits import read_oifits

data = OIData(read_oifits(sorted(glob.glob("calibrated/*/*.fits")), insname="GRAVITY_SC"))
```

Main options (`--help` for all):

| Option | Default | Meaning |
|---|---|---|
| `--dual-field`, `--cal NAME[=DIAM]` | `ESO PRO CATG` | how calibrator exposures are recognised |
| `--insname` | GRAVITY science channel, both polarisations | instruments to calibrate; the fringe tracker is never used |
| `--window LO HI`, `--bin N` | all channels, 1 | channel window (µm) and binning |
| `--interp linear\|gaussian` | `linear` | linear in time, or Gaussian weights (the pipeline manual's formula, with error floors adapted from the source) with `--timescale` (1 h) |
| `--tf-error loo\|scatter\|propagated\|none` | `loo` | leave-one-out prediction error of the transfer function, the night's scatter, the calibrators' statistical errors only (with `--template`, whose modes carry the systematic), or none |
| `--tf-modes` | off | write the `loo` or `scatter` transfer-function error on V² and closure phases as nuisance modes (`<file>.tf_modes.npz`), one per night and baseline or triangle, shared by that night's files, instead of adding it to each datum's error |
| `--tf-mode-shape` | `poly1` | with `--tf-modes`, how each mode's channel shape is smoothed: `poly1` or `poly2` (a polynomial in wavelength, for a transfer-function error with a colour), `flat` (its channel mean; underestimates errors at realistic chromatic transfer functions, pull rms 3.1 at colour 0.36 against about 1.3 for `poly1`; our own Monte Carlo, not yet checked against ESO/VLTI documentation), or `raw` (the old channel-by-channel estimate, for comparison only) |
| `--cp-only-when-uncalibrated` | off | for a science file on a night with no calibrator, write a closure-phase-only product instead of skipping it. **Its T3PHI is uncalibrated**: it keeps the instrumental closure phase and has statistical errors only, so fit it with the template's closure-phase offset modes (`--template`, t3phi modes only) and treat the errors as a lower bound. V² and T3AMP are flagged, header `CPONLY` |
| `--cp-only NAME` | none | the same flagging for science files whose input name, output stem (`<target>_<YYYYMMDDTHHMMSS>`) or DATE-OBS contains `NAME` (repeatable; a name that marks nothing is an error). With calibrators the closure phases are calibrated as usual and only V² and T3AMP are set aside |
| `--vis` | off | also calibrate `OI_VIS` (VISAMP, VISPHI) |
| `--smooth-tf DEG` | off | smooth the amplitude transfer function across channels with `gravity_viscal`'s polynomial fit (it uses 5) |
| `--match-gravity-viscal` | off | calibrate exactly as ESO's `gravity_viscal` (pipeline 1.9.7) does, to reproduce ESO's calibrated products; with `--delta-time-calib`, `--maxdeg-tfvis`, `--nsmooth-tfvis`, `--no-smoothing`, `--force-calib` |
| `--no-pipeline-diameters` | off | ignore the calibrator diameters GRAVITY's pipeline recorded in `*_CAL_VIS` products (used unless `--cal` names the calibrator) |
| `--diameter DIAM[±ERR]` | point source | default diameter (and error, mas) of every calibrator; a value given with `--cal` overrides it for that calibrator |

### Caveats
- **The in-field calibrator must be unresolved.** In dual-field mode the fringe-tracking star calibrates the target only if it is a point source, or if its diameter is given with `--cal`.
- **Off-axis targets lose coherence that the on-axis star does not.** Anisoplanatism, off-centre injection and field-dependent aberrations differ between the fringe-tracking star and a target arcseconds away, so calibrated V² carry a loss that this does not remove ([review](design/gravity_calibration_review.md) §7). Closure phases are much less affected. With `--dual-field`, `index.csv` gives each exposure's predicted anisoplanatic V² loss (GRAVITY+ Collaboration 2022, eqs. 4–5, from the seeing and isoplanatic angle in the header), and a warning flags losses above `--aniso-warn` (2%). It is a size estimate, not a correction: single observations scatter about the model.
- **By default the added calibration error is treated as independent between channels and exposures**, although it is common to all channels of a baseline and to every science exposure of the night interpolated from the same calibrators: N exposures would wrongly average it down by √N. `--tf-modes` writes it instead as one nuisance mode per night and baseline (V², on log |V|) or triangle (closure phases, degrees) in `<file>.tf_modes.npz`, named `tf-vis2-base28-23@<night>` and `tf-t3phi-tri28-23-18@<night>`, so that the same name in several files is one parameter; only the part independent between channels stays in the errors. The mode's shape across channels is estimated from a handful of calibrators and so is noisy (a coefficient of variation of about 0.16 for V² and 0.25 for closure phases with seven calibrators); a mode spanning hundreds of data would then pin its amplitude to that noise and leave the shape's scatter times the realized error unmodelled, so fitted parameters came out about twofold too precise. `--tf-mode-shape` therefore smooths the shape (default `poly1`) and returns the variance this removes to the per-channel remainder, so each channel's total transfer-function error is unchanged (where the smoothed mode exceeds it, the remainder is clipped at zero, with a warning). Use either `--tf-modes` or a `--template` (whose modes also carry the calibration residual), not both.
- **GRAVITY channels are about twofold oversampled.** Bin by 2 or more unless channel correlations are modelled.
- **Calibrator diameter errors are nuisance modes, not errors.** Give them with `--cal "NAME=0.45±0.03"` or `--diameter 0.2±0.1` (a default for every calibrator, e.g. the fringe-tracking star; `--cal` values override it field by field). Each calibrated file then gets a `.modes.npz` beside it: the change in V² and T3AMP per 1σ of each calibrator's diameter, one mode per calibrator and night, laid out like the file's tables. `index.csv` reports the largest, in units of the V² error. virgil's correlated-nuisance machinery (`OIData.with_gains` and `with_closure_offsets`) is meant to marginalise these. Each mode is written twice. `<insname>/<observable>/<calibrator>@<night>` is absolute, the change in the calibrated values, built from the data. `<insname>/<observable>/rel-<calibrator>@<night>` is the same parameter relative to the model: d log\|V\| per 1σ for V² (and VISAMP) and d ln T3AMP for T3AMP, independent of the science values (`diameter_modes(..., relative=True)`). Pass the `rel-` V² modes to virgil's `OIData.with_gains(modes=...)`, which avoids the bias of data-built covariances (design/gravity_calibration_review.md, footgun 11). Use one form or the other, never both.
- **V², T3AMP and T3PHI are calibrated; `OI_VIS` only with `--vis`.** VISAMP must be a visibility (`AMPTYP` absolute), not a correlated flux. In `OI_T3` and `OI_VIS` the single `FLAG` column follows the phase, and an amplitude flagged on its own is written as NaN. virgil fits V² and T3PHI when they are present.

### Calibration templates
The calibrator PCA (`virgil-vlti-calpca`; `design/gravity_calibrator_pca.md`) measures the systematics that calibration leaves behind. A **template** packages them for one instrument configuration, as nuisance modes in the form virgil's correlated calibration nuisances (`with_gains`, `with_closure_offsets`) are designed to take:
- **V²:** per-telescope and per-baseline gains on log |V|, with chromatic shapes;
- **closure phases:** optional per-triangle offsets.

`virgil-vlti-calibrate --template NAME` writes the modes beside each calibrated file (`<file>.template_modes.npz`, laid out like its tables, one mode per frame and telescope, baseline or triangle). It refuses data in another configuration (field, resolution, polarisation; `--any-template-config` only warns) and warns about channels outside the template's wavelengths.

Shipped:

| Template | Configuration | V² components | Closure-phase components |
|---|---|---|---|
| `gravity_single_medium_combined_2023-2025` | GRAVITY single field, MEDIUM, combined polarisation, 2023–25, mostly ATs, 2.00–2.45 µm | 4: 1.9%, 0.52%, 0.43%, 0.29% per channel in V², 41–81% telescope-based | 2: 0.33°, 0.25° per channel |
| `gravity_single_high_combined_2023-2025` | GRAVITY single field, HIGH, combined polarisation, 2023–25, ATs only | 3 | 2: 0.37°, 0.17° per channel |
| `gravity_single_medium_combined_2023-2025_nowr` | as the first, rebuilt on a smaller calibrator set (71 files, 25 nights); prefer the first for general use | 4 | 2: 0.29°, 0.22° per channel |

Build others from a `virgil-vlti-calpca` run with `virgil-vlti-template residuals.npz --name … --k-vis2 … --k-t3phi … --config …`.

### Dual-field targets
The dual-field recipe was first written for GRAVITY observations where the science fibre alternates between the target and a companion that carries the fringe tracker.

## Development

```bash
uv venv --python 3.11 .venv
uv pip install --python .venv/bin/python -e ".[dev]"
.venv/bin/python -m pytest
```

The tests build synthetic GRAVITY-layout products with a known truth (`tests/synthetic.py`).

The documentation is built with Zensical, as virgil's is. Its home page is this README, and its tutorials are the executed notebooks in `notebooks/` (simulated GRAVITY data from `virgil_vlti.synthetic`), executed on a cluster rather than a laptop. `scripts/sync_docs.py` copies both into `docs/` and links the design notes there:

```bash
uv pip install --python .venv/bin/python -e ".[docs]"
.venv/bin/python scripts/sync_docs.py
.venv/bin/zensical build --clean --strict
```
