# Command-line tools

Five commands come with the package. The [tutorials](calibration.md) show what each one is for; this page lists their options. `--help` on each gives the full list.

## `virgil-vlti-calibrate`

Calibrates a folder of uncalibrated products, night by night, and writes one OIFITS file per science exposure, plus an `index.csv`. Its options are in the table on the [home page](index.md); [tutorial 1](calibration.md) runs it.

With `--tf-modes` (and `--tf-error loo` or `scatter`), the transfer function's error on V² and closure phases is not added to each datum's error but written beside each file as `<file>.tf_modes.npz`, in the same layout as template modes (`"<insname>/<observable>/<mode>"`, rows × channels):

| Mode name | Meaning |
|---|---|
| `tf-vis2-base28-23@<night>` | the transfer-function error of baseline 28–23 on that night: d log\|V\| per 1σ (half the relative V² error), on every channel and every science file of the night |
| `tf-t3phi-tri28-23-18@<night>` | the same for a triangle's closure phase, degrees per 1σ |

For MATISSE products, `--n-phase {refuse,negate,keep}` (default `refuse`) decides what happens to N-band products whose phases are flipped by 180° (`mat_raw_estimates` 1.5.2 to before 2.0.0), and `--n-phase-unknown` (default `refuse`) to N-band products of an earlier or unrecorded version. `negate` multiplies T3PHI and VISPHI by −1 on read and records it in `HISTORY`. Merged products and `corrFlux=true` products are always refused; see `virgil-vlti-inspect`.

Each mode's per-channel size is the measured transfer-function error, so it is fully correlated across channels. Its shape across channels is smoothed by `--tf-mode-shape {flat,poly1,poly2,raw}` (default `poly1`; `flat` is the channel mean; `poly1` and `poly2` fit a polynomial in wavelength, for a transfer-function error with a colour; `raw` is the unsmoothed channel-by-channel estimate, kept for comparison, whose scatter leaves part of the error unmodelled). The variance the smoothing takes from the mode goes back to the independent per-channel part, so each channel's total is unchanged; where the smoothed mode exceeds a channel's total error, that channel's remainder is clipped at zero and its total grows to the mode's size (a warning gives the count and the added variance). The `.tf_modes.npz` layout and names are unchanged; only the values of the modes (and the errors) differ from earlier versions, which wrote `raw`. `poly1` is the default because `flat` underestimates errors at realistic chromatic transfer functions (pull rms 3.1 at colour 0.36, which is λ⁻² across K, against about 1.3 for `poly1` at every colour; `flat` is calibrated only for achromatic errors). This is from our own Monte Carlo and has not yet been checked against ESO/VLTI documentation of transfer-function chromaticity. The part of the error that the calibrators show to be independent between channels (a rank-one fit to their residuals' channel correlations) stays in the errors. Every science exposure of a night shares one mode per baseline: exact when all are interpolated between the same calibrators, an overstatement of the correlation when they lean on different ones. T3AMP and `OI_VIS` keep the diagonal term.

## `virgil-vlti-pairs`

The [pair test](pair_test.md) on a set of calibrated products. It prints each nested model's ΔlnL and fitted widths.

```bash
virgil-vlti-pairs calibrated/*.fits --json pairs.json
```

| Option | Default | Meaning |
|---|---|---|
| `inputs` | | product files, or folders of `*.fits` |
| `--insname` | each file's only one | the instrument (PIONIER puts its wavelength solution in `INSNAME`, so it changes from file to file) |
| `--max-gap` | 30 | longest gap between paired frames, in minutes |
| `--field` | all | keep only products in this field configuration (as `virgil-vlti-calpca --field`) |
| `--by-star-night` | off | one test per star and night instead of one pooled test; `--json` then gets a list |
| `--json` | off | also write the fitted widths here |

## `virgil-vlti-calpca`

The calibrator PCA ([tutorial 3](templates.md)) on archival calibrator products. Each calibrator is calibrated by the others of its night and compared with its own uniform-disk model; the residuals' principal components are cross-validated by night.

```bash
virgil-vlti-calpca calibrators/ --out pca/ --k-max 6
```

| Option | Default | Meaning |
|---|---|---|
| `inputs` | | calibrator products, folders of `*.fits`, or a `residuals.npz` to re-analyse |
| `--out` | required | output folder |
| `--insname` | `GRAVITY_SC` | the instrument |
| `--k-max` | 6 | most components tried |
| `--folds` | 5 | cross-validation folds |
| `--window LO HI` | 2.0 2.45 | channels kept, in µm (GRAVITY's band edges are noisy) |
| `--min-gain` | 0.02 | keep a component if it lowers the held-out χ² by this fraction beyond noise |
| `--clip` | 5 | clip outlier rows beyond this many robust deviations (0: off) |
| `--screen-factor` | 5 | exclude calibrators scoring above this times the median star (0: off) |
| `--field` | all | keep only `SINGLE`, `DUAL` (any dual field), `DUAL_ON` (science fibre on the fringe-tracking star), `DUAL_OFF` (on another star in the field) or `DUAL_WIDE` products; configurations never calibrate each other |
| `--references` | `other-star` | what calibrates each calibrator: the night's other stars; `same-star`, the star's own other exposures (exposure-to-exposure changes only; dual-field programmes rarely observe two calibrators a night); or `any` |

It writes `calpca.npz` (wavelengths, components, widths), `summary.json` (inputs, cross-validation, calibrators excluded), `residuals.npz` (the input to `virgil-vlti-template`) and `components.png`. On archival data it reads hundreds of products: run it on a cluster.

## `virgil-vlti-template`

Builds a calibration template from `virgil-vlti-calpca`'s residuals, with the numbers of components chosen from its cross-validation.

```bash
virgil-vlti-template pca/residuals.npz --epoch 2023-2025 --k-vis2 4 --k-t3phi 2 \
    --config field=SINGLE resolution=MEDIUM polarisation=COMBINED
```

writes `gravity_single_medium_combined_2023-2025.npz`: the name is `gravity_<field>_<resolution>_<polarisation>_<epoch>`, with the field as `single`, `dual`, `dual-on`, `dual-off` or `dual-wide`.

| Option | Default | Meaning |
|---|---|---|
| `residuals` | | `residuals.npz` from `virgil-vlti-calpca` |
| `--name` | from `--config` and `--epoch` | the template's name (the file is `<name>.npz`) |
| `--epoch` | | the years the data span, for the default name |
| `--k-vis2`, `--k-t3phi` | required, 0 | components to keep |
| `--config KEY=VALUE ...` | | which data it describes: `field`, `resolution`, `polarisation`, optionally `telescopes` (`AT` or `UT`); data in another configuration are refused |
| `--provenance KEY=VALUE ...` | | where the residuals came from |
| `--folds`, `--clip` | 5, 5.0 | as in the PCA |
| `--out` | `.` | output folder |

`virgil-vlti-calibrate --template NAME` then writes `<file>.template_modes.npz` beside each calibrated file: arrays named `"<insname>/<observable>/<mode>"`, each shaped like that observable's table (rows × channels). The mode names say what each one is:

| Mode name | Meaning |
|---|---|
| `vis2-pc1-tel28@<file>#frame0` | a telescope gain: component 1, station 28, frame 0 of that file |
| `vis2-pc1-base28-23@…` | a baseline gain |
| `t3phi-pc1-tri28-23-18@…` | an optional closure-phase offset |

Modes with the same name are one parameter. V² modes are on log |V| per 1σ; closure-phase modes are in degrees per 1σ.

## `virgil-vlti-inspect`

Says what each product is and whether `virgil-vlti-calibrate` may take it as an input. For MATISSE it prints the band, resolution, BCD position, chopping, `mat_raw_estimates` version, N-band phase convention, merge evidence and `replaceTel`, and a verdict: `ok`, or why the product is refused (calibrated, merged, `corrFlux`, N-band phase convention). Exit status 0 if every file is `ok`, else 1.

```bash
virgil-vlti-inspect reduced/*.fits
virgil-vlti-inspect reduced/*.fits --json
```

A position with only one BCD in is printed as its raw keys (`BCD1=IN BCD2=OUT (device order unverified)`): which header key is which BCD device is not yet established (design/matisse_plan.md §3.1, §7.1).
