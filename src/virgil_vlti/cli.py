"""``virgil-vlti-calibrate``: calibrate a folder of VLTI products, night by night.

Examples::

    # GRAVITY dual field, the fringe-tracking star as in-field calibrator
    virgil-vlti-calibrate reduced/ --out calibrated --dual-field --window 2.21 2.40 --bin 4

    # separate calibrator stars, with uniform-disk diameters in mas
    virgil-vlti-calibrate reduced/ --out calibrated --cal "HD 12345=0.45±0.03" --cal "HD 678"

Writes one OIFITS file per science exposure to ``<out>/<night>/``, an index
``<out>/index.csv``, and the science-to-calibrator flux ratio
``<out>/flux_ratio[_<INSNAME>].csv`` when the products carry ``OI_FLUX``.
When calibrator diameters have errors, each calibrated file gets a
``.modes.npz`` beside it: the diameter nuisance modes
(:func:`virgil_vlti.write.write_modes`): ``<insname>/<observable>/<cal>@<night>``
the absolute change per 1σ, and ``<insname>/<observable>/rel-<cal>@<night>``
the same parameter as a model-relative mode (d log|V| for V², the one for
``OIData.with_gains(modes=...)``). Use one form or the other, never both.
With ``--tf-modes``, each also gets a ``.tf_modes.npz``: the transfer
function's error as one mode per night and baseline (V²) or triangle
(closure phases), shared by the night's files
(:func:`virgil_vlti.calibrate.transfer_function_modes`), each with a
smoothed channel shape (``--tf-mode-shape``,
:func:`virgil_vlti.calibrate.regularize_shape`).
"""

import argparse
import csv
import os
import re
import warnings

import numpy as np

from . import gravity, matisse
from .calibrate import (
    TF_MODE_SHAPES,
    calibrate_exposure,
    closure_phase_only,
    diameter_modes,
    diameter_of,
    flux_ratio,
    group_by_night,
    mark_closure_phase_only,
    transfer_function_modes,
)
from .products import classify, expand_paths, instrument_of, normalise_name, read_exposure
from .templates import Template, template_modes
from .viscal import calibrate_like_viscal
from .write import write_calibrated, write_modes

# Mode-name prefix of the model-relative diameter modes in ``.modes.npz``.
RELATIVE_PREFIX = "rel-"


def _diameter(text):
    """``"0.45"`` or ``"0.45±0.03"`` (or ``+-``) -> (diameter, error) in mas."""
    value, _, error = text.replace("+-", "±").partition("±")
    try:
        out = float(value), float(error) if error else None
    except ValueError as err:
        raise argparse.ArgumentTypeError(f"not a diameter: {text!r}") from err
    if not all(np.isfinite(x) and x >= 0 for x in out if x is not None):
        raise argparse.ArgumentTypeError(
            f"diameters and errors must be finite and non-negative: {text!r}"
        )
    return out


def _degree(text):
    try:
        value = int(text)
    except ValueError as err:
        raise argparse.ArgumentTypeError(f"not an integer: {text!r}") from err
    if value < 0:
        raise argparse.ArgumentTypeError(f"must be non-negative: {text!r}")
    return value


def _calibrator(text):
    name, _, diam = text.partition("=")
    return (name.strip(), *_diameter(diam)) if diam else (name.strip(), None, None)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("inputs", nargs="+", help="product files, or folders of *.fits")
    p.add_argument("--out", required=True, help="output folder")
    how = p.add_mutually_exclusive_group()
    how.add_argument(
        "--dual-field",
        action="store_true",
        help="GRAVITY: frames on the fringe-tracking star are the calibrator",
    )
    how.add_argument(
        "--cal",
        action="append",
        type=_calibrator,
        default=[],
        metavar="NAME[=DIAM[±ERR]]",
        help="a calibrator target, with its uniform-disk diameter and error in mas (repeatable)",
    )
    p.add_argument(
        "--diameter",
        type=_diameter,
        default=None,
        metavar="DIAM[±ERR]",
        help="the default uniform-disk diameter (and error) in mas of every calibrator; a "
        "diameter or error given with --cal overrides it for that calibrator. E.g. the "
        "fringe-tracking star with --dual-field",
    )
    p.add_argument(
        "--no-pipeline-diameters",
        action="store_true",
        help="ignore the calibrator diameters that GRAVITY's pipeline recorded in *_CAL_VIS "
        "products (ESO QC TF CAL_SC DIAM), used by default unless --cal names the calibrator",
    )
    p.add_argument(
        "--smooth-tf",
        type=_degree,
        default=None,
        metavar="DEG",
        help="smooth the amplitude transfer function across channels with gravity_viscal's "
        "polynomial fit of this degree (it uses 5)",
    )
    p.add_argument(
        "--template",
        default=None,
        metavar="NAME|PATH",
        help="write a calibration template's nuisance modes beside each file "
        "(<file>.template_modes.npz); see virgil_vlti.templates",
    )
    p.add_argument(
        "--cp-only-when-uncalibrated",
        action="store_true",
        help="for a science file on a night with no calibrator, write a closure-phase-only "
        "product instead of skipping it: the science T3PHI windowed and binned as usual but "
        "UNCALIBRATED (instrumental closure phase included, statistical errors only; use "
        "--template for its offset modes); V² and T3AMP flagged; header CPONLY",
    )
    p.add_argument(
        "--cp-only",
        action="append",
        default=[],
        metavar="NAME",
        help="write a closure-phase-only product (as --cp-only-when-uncalibrated) for science "
        "files whose input file name, output stem (e.g. ll_aqr_20211108T011826) or DATE-OBS "
        "contains NAME, for files whose V² the caller marks as suspect (repeatable; a name "
        "that marks nothing is an error). With calibrators the closure phases are calibrated "
        "as usual and V² and T3AMP are flagged",
    )
    p.add_argument(
        "--any-template-config",
        action="store_true",
        help="warn, rather than stop, when --template was measured in another configuration "
        "(field, resolution, polarisation, telescopes)",
    )
    parity = p.add_argument_group(
        "gravity_viscal parity",
        "Calibrate as ESO's gravity_viscal does, following the 1.9.7 source (checked against "
        "1.9.7 products only; other versions unverified), to reproduce its "
        "*_SCI_VIS_CALIBRATED products: virgil_vlti.viscal. Its own defaults are used; "
        "--window, --bin, --interp, --tf-error and --smooth-tf do not apply.",
    )
    parity.add_argument(
        "--match-gravity-viscal", action="store_true", help="use the gravity_viscal port"
    )
    parity.add_argument(
        "--delta-time-calib",
        type=float,
        default=3600.0,
        metavar="S",
        help="gravity_viscal delta-time-calib (s) [3600]",
    )
    parity.add_argument(
        "--maxdeg-tfvis",
        type=_degree,
        default=5,
        metavar="DEG",
        help="gravity_viscal maxdeg-tfvis-sc [5]",
    )
    parity.add_argument(
        "--nsmooth-tfvis",
        type=_degree,
        default=0,
        metavar="N",
        help="gravity_viscal nsmooth-tfvis-sc [0]",
    )
    parity.add_argument(
        "--no-smoothing", action="store_true", help="gravity_viscal smoothing=FALSE"
    )
    parity.add_argument(
        "--force-calib",
        action="store_true",
        help="gravity_viscal force-calib: keep calibrators with another setup",
    )
    p.add_argument(
        "--aniso-warn",
        type=float,
        default=0.02,
        metavar="LOSS",
        help="dual field: warn when the predicted anisoplanatic V² loss exceeds this "
        "(default 0.02)",
    )
    p.add_argument(
        "--onaxis-mas",
        type=float,
        default=None,
        help="dual field: largest fibre offset counted as on axis",
    )
    p.add_argument(
        "--insname",
        action="append",
        default=None,
        help="instrument(s) to calibrate (default: GRAVITY science channel)",
    )
    p.add_argument("--window", nargs=2, type=float, metavar=("LO", "HI"), help="µm")
    p.add_argument("--bin", type=int, default=1, help="channels per bin")
    p.add_argument("--interp", choices=("linear", "gaussian"), default="linear")
    p.add_argument("--timescale", type=float, default=3600.0, help="seconds, for --interp gaussian")
    p.add_argument(
        "--tf-error",
        choices=("loo", "scatter", "propagated", "none"),
        default="loo",
        help="the transfer function's error added to the data: a leave-one-out prediction "
        "error (default), the night's scatter, the calibrators' statistical errors carried "
        "through the interpolation (use with --template, whose modes carry the systematic), "
        "or none",
    )
    p.add_argument(
        "--tf-modes",
        action="store_true",
        help="write the --tf-error (loo or scatter) on V² and closure phases as nuisance modes "
        "beside each file (<file>.tf_modes.npz), one per night and baseline or triangle, "
        "instead of adding it to each datum's error; only its part independent between "
        "channels is still added",
    )
    p.add_argument(
        "--tf-mode-shape",
        choices=TF_MODE_SHAPES,
        default=None,
        help="with --tf-modes, how each mode's channel shape is smoothed: flat, its "
        "mean over the channels (calibrated only for achromatic errors); poly1 (default) or "
        "poly2, a least-squares polynomial in wavelength, for a transfer-function error with a "
        "colour; raw, the channel-by-channel estimate, whose "
        "scatter leaves part of the error unmodelled (for comparison only). The variance taken "
        "from the shape goes back to the independent per-channel term, so each channel's total "
        "is unchanged. poly1 is the default because flat underestimates errors at realistic "
        "chromatic transfer functions (pull rms 3.1 at colour 0.36, λ^-2 across K; poly1 about "
        "1.3 at every colour); this is our own Monte Carlo, not yet checked against ESO/VLTI "
        "documentation",
    )
    p.add_argument(
        "--flux-norm",
        nargs=2,
        type=float,
        metavar=("LO", "HI"),
        help="µm window to normalise the flux ratio (default: --window)",
    )
    p.add_argument(
        "--vis",
        action="store_true",
        help="also calibrate OI_VIS (VISAMP, VISPHI); virgil reads V² and closure phases first",
    )
    p.add_argument(
        "--allow-mismatch",
        action="store_true",
        help="warn, rather than stop, when calibrators were reduced with other pipeline "
        "settings (another instrument mode is always refused)",
    )
    p.add_argument(
        "--allow-unrecorded",
        action="store_true",
        help="warn, rather than stop, when a header key or pipeline parameter that must match "
        "(resolution, polarisation mode, vis-correction-sc) is missing on either side, so "
        "that no check could be made",
    )
    p.add_argument("--prefix", default=None, help="output file prefix (default: target)")
    mat = p.add_argument_group(
        "MATISSE",
        "N-band phases of mat_raw_estimates 1.5.2 to 1.9.x are rotated by 180 degrees "
        "([PM] §2.5); those of earlier or unrecorded versions are unverified. Both are refused "
        "by default (design/matisse_plan.md, §7, decision 2).",
    )
    mat.add_argument(
        "--n-phase",
        choices=matisse.N_PHASE_CHOICES,
        default="refuse",
        help="N-band products whose phase convention is 'flipped' (1.5.2 <= version < 2.0.0): "
        "refuse (default), negate T3PHI and VISPHI on read, or keep them as they are",
    )
    mat.add_argument(
        "--n-phase-unknown",
        choices=matisse.N_PHASE_CHOICES,
        default="refuse",
        help="N-band products whose phase convention is 'unknown' (no version, or before "
        "1.5.2): refuse (default), keep, or negate",
    )
    return p.parse_args(argv)


CP_ONLY_UNCALIBRATED = (
    "Closure phases only (virgil-vlti), UNCALIBRATED (no calibrator): science T3PHI as observed, "
    "instrumental closure phase included, statistical errors only; window {window} um, binned "
    "by {bin}. V2 and T3AMP flagged and not to be used."
)
CP_ONLY_CALIBRATED = (
    " Closure phases only: V2 and T3AMP flagged and not to be used (marked suspect)."
)


def _marked(names, sci, prefix=None):
    """The ``--cp-only`` names that mark this science file, matched against its input file
    name, its output stem (as written, ``<prefix>_<stamp>``) and its DATE-OBS."""
    stem = f"{prefix or _slug(sci.target)}_{_stamp(sci)}".lower()
    texts = (os.path.basename(str(sci.path)).lower(), stem, str(sci.date).lower())
    return [n for n in names if n.strip() and any(n.strip().lower() in t for t in texts)]


def _provenance(args, cals, diameters, errors):
    """The ``HIERARCH VIRGIL`` cards for one output: calibrators, diameters (mas), error mode."""
    names = sorted({c.target for c in cals})
    out = {"NCAL": len(names)}
    for k, n in enumerate(names, 1):  # one short card each: a long list would overflow a card
        out[f"CAL{k} NAME"] = n
        out[f"CAL{k} DIAM"] = diameter_of(diameters, n)
        out[f"CAL{k} DIAMERR"] = diameter_of(errors, n)
    out["MATCH VISCAL"] = bool(args.match_gravity_viscal)
    if not args.match_gravity_viscal:
        out["INTERP"] = args.interp
        out["TF ERROR"] = args.tf_error if cals else "none"
    return out


def _stamp(sci):
    return sci.date[:19].replace(":", "").replace("-", "")


def _closure_only(modes, cp_only):
    """The template's modes, without the V² ones for a closure-phase-only file."""
    if not cp_only:
        return modes
    return {m: by_obs for m, by_obs in modes.items() if "t3phi" in by_obs}


def _slug(text):
    return re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_").lower()


def main(argv=None):
    args = parse_args(argv)
    exposures = [read_exposure(p, args.insname) for p in expand_paths(args.inputs)]
    diameters = {name: d for name, d, _ in args.cal if d}
    errors = {name: e for name, _, e in args.cal if e}
    if args.diameter:
        diameters[None], errors[None] = args.diameter[0], args.diameter[1] or 0.0
    if not args.no_pipeline_diameters:
        named = {normalise_name(name) for name, _, _ in args.cal}
        for e in exposures:
            if normalise_name(e.target) not in named and e.pipeline_diameter:
                diameters[e.target], err = e.pipeline_diameter
                if err:
                    errors[e.target] = err
    classify(
        exposures,
        dual_field=args.dual_field,
        calibrators=[name for name, _, _ in args.cal],
        onaxis_mas=args.onaxis_mas,
    )
    for e in exposures:  # merged, correlated-flux and N-band phase refusals
        if instrument_of(e.header) == "MATISSE":
            matisse.screen(e, n_phase=args.n_phase, n_phase_unknown=args.n_phase_unknown)
    window = tuple(args.window) if args.window else None
    history = (
        f"Calibrated by virgil-vlti ({'in-field' if args.dual_field else 'calibrator'} "
        f"transfer function, {args.interp} in time, {args.tf_error} error"
        + ("" if args.smooth_tf is None else f", smoothed in wavelength (degree {args.smooth_tf})")
        + "); "
        f"window {window or 'all'} um, binned by {args.bin}."
    )
    if args.match_gravity_viscal:
        if window or args.bin != 1 or args.smooth_tf is not None:
            raise SystemExit("--match-gravity-viscal calibrates every channel: drop --window, "
                             "--bin and --smooth-tf")  # fmt: skip
        history = (
            "Calibrated by virgil-vlti's port of gravity_viscal as in the 1.9.7 source "
            "(verified against 1.9.7 products only; the input's own pipeline version is "
            "recorded in HIERARCH VIRGIL INPUT PIPE ID): "
            f"delta-time-calib {args.delta_time_calib} s, smoothing {not args.no_smoothing}, "
            f"maxdeg-tfvis-sc {args.maxdeg_tfvis}, nsmooth-tfvis-sc {args.nsmooth_tfvis}, "
            f"force-calib {args.force_calib}."
        )
    if args.tf_modes:
        if args.match_gravity_viscal:
            raise SystemExit("--tf-modes does not apply with --match-gravity-viscal")
        if args.tf_error not in ("loo", "scatter"):
            raise SystemExit("--tf-modes needs --tf-error loo or scatter: the error it models")
        if args.template:
            raise SystemExit("use --tf-modes or --template, not both: both carry the calibration "
                             "residual, which would be counted twice")  # fmt: skip
        args.tf_mode_shape = args.tf_mode_shape or "poly1"
        history += (
            " Transfer-function error on V2 and T3PHI as modes (.tf_modes.npz), "
            f"shape {args.tf_mode_shape}."
        )
    elif args.tf_mode_shape is not None:
        raise SystemExit("--tf-mode-shape applies only with --tf-modes")
    template = Template.load(args.template) if args.template else None
    if template is not None and not args.any_template_config:
        for sci in (e for e in exposures if e.kind == "SCI"):  # before anything is written
            problems = template.mismatches(sci.header)
            if problems:
                raise SystemExit(
                    f"{os.path.basename(str(sci.path))}: template {template.name!r} is for "
                    f"another configuration ({'; '.join(problems)}); "
                    "--any-template-config to use it anyway"
                )
    rows = []
    used = set()
    for night, group in group_by_night(exposures).items():
        cals = [e for e in group if e.kind == "CAL"]
        for sci in (e for e in group if e.kind == "SCI"):
            hit = _marked(args.cp_only, sci, args.prefix)
            used.update(hit)
            if not cals and not args.cp_only_when_uncalibrated and not hit:
                print(f"{sci.date}  no calibrator that night; skipped")
                continue
            cp_only = not cals or bool(hit)
            tf_modes = args.tf_modes and bool(cals)  # needs calibrators
            options = {
                "window_um": window,
                "nbin": args.bin,
                "method": args.interp,
                "error": args.tf_error,
                "timescale_s": args.timescale,
                "diameters": diameters,
                "vis": args.vis and not cp_only,
                "smooth_degree": args.smooth_tf,
            }
            if not cals:
                results = [
                    closure_phase_only(sci, ins, window_um=window, nbin=args.bin)
                    for ins in sci.insnames
                ]
            elif args.match_gravity_viscal:
                viscal_params = {
                    "delta_time": args.delta_time_calib,
                    "smoothing": not args.no_smoothing,
                    "maxdeg": args.maxdeg_tfvis,
                    "nsmooth": args.nsmooth_tfvis,
                    "force_calib": args.force_calib,
                }
                parity_diameters = {
                    c.target: diameter_of(diameters, c.target) for c in cals
                }  # resolved per calibrator, so viscal looks them up by exact name
                results = [
                    calibrate_like_viscal(sci, cals, ins, parity_diameters, **viscal_params)
                    for ins in sci.insnames
                ]
            else:
                results = [
                    calibrate_exposure(
                        sci,
                        cals,
                        ins,
                        strict=not args.allow_mismatch,
                        allow_unrecorded=args.allow_unrecorded,
                        tf_modes=tf_modes,
                        tf_mode_shape=args.tf_mode_shape or "poly1",
                        **options,
                    )
                    for ins in sci.insnames
                ]
            if cp_only and cals:  # calibrated closure phases, V² set aside
                results = [mark_closure_phase_only(r, sci) for r in results]
            modes = {}  # the pipeline does not propagate diameter errors (manual §9.20)
            if any(errors.values()) and not args.match_gravity_viscal and not cp_only:
                modes = {
                    ins: diameter_modes(sci, cals, ins, diameter_errors=errors, **options)
                    for ins in sci.insnames
                }
                relative = {  # the same parameters, relative to the model
                    ins: {
                        f"{RELATIVE_PREFIX}{mode}": by_obs
                        for mode, by_obs in diameter_modes(
                            sci, cals, ins, diameter_errors=errors, relative=True, **options
                        ).items()
                    }
                    for ins in sci.insnames
                }
            os.makedirs(os.path.join(args.out, night), exist_ok=True)
            stamp = _stamp(sci)
            name = f"{args.prefix or _slug(sci.target)}_{stamp}.fits"
            path = write_calibrated(
                sci,
                results,
                os.path.join(args.out, night, name),
                CP_ONLY_UNCALIBRATED.format(window=window or "all", bin=args.bin)
                if not cals
                else history + (CP_ONLY_CALIBRATED if cp_only else ""),
                provenance=_provenance(args, cals, diameters, errors),
            )
            if modes:
                both = {ins: {**modes[ins], **relative[ins]} for ins in modes}
                write_modes(path[: -len(".fits")] + ".modes.npz", both)
            if tf_modes:
                write_modes(
                    path[: -len(".fits")] + ".tf_modes.npz",
                    {r.insname: transfer_function_modes(sci, r) for r in results},
                )
            if template is not None:
                write_modes(
                    path[: -len(".fits")] + ".template_modes.npz",
                    {
                        r.insname: _closure_only(
                            template_modes(
                                template,
                                sci,
                                r.insname,
                                wave=r.wave,
                                strict=not args.any_template_config,
                            ),
                            cp_only,
                        )  # fmt: skip
                        for r in results
                    },  # on the written (windowed, binned) channels
                )
            good = 1 - np.mean([np.mean(r["vis2"].flag) for r in results])
            # The largest diameter mode on V², in units of the V² error.
            worst = max(
                (
                    float(np.max(np.abs(by_obs["vis2"])[ok] / r["vis2"].error[ok]))
                    for r in results
                    for by_obs in modes.get(r.insname, {}).values()
                    if (ok := ~r["vis2"].flag & (r["vis2"].error > 0)).any()
                ),
                default=0.0,
            )
            aniso = None
            if args.dual_field:
                losses = [gravity.anisoplanatic_v2_loss(sci.header, r.wave) for r in results]
                losses = [np.max(x) for x in losses if x is not None]
                aniso = round(float(max(losses)), 4) if losses else None
                if aniso is not None and aniso > args.aniso_warn:
                    warnings.warn(
                        f"{sci.date}: predicted anisoplanatic V² loss up to {aniso:.3f} at "
                        f"{gravity.fibre_separation_mas(sci.header) / 1000:.2f} arcsec "
                        "from the fringe tracker; the on-axis calibrator does not remove it",
                        stacklevel=2,
                    )
            rows.append(
                {
                    "night": night,
                    "date_obs": sci.date,
                    "source": os.path.basename(sci.path),
                    "calibrated": os.path.relpath(path, args.out),
                    "n_cal_frames": len(cals),
                    "cp_only": cp_only,
                    "good_v2_fraction": round(float(good), 3),
                    "max_diameter_mode_sigma": round(worst, 3),
                    "anisoplanatic_v2_loss": "" if aniso is None else aniso,
                }
            )
            print(f"{sci.date}  {len(cals)} CAL  -> {rows[-1]['calibrated']}  (good V2 {good:.2f})")
    unmatched = [n for n in args.cp_only if n not in used]
    if unmatched:
        raise SystemExit(
            f"--cp-only {unmatched} marked no science file (matched against input file names, "
            "output stems <target>_<YYYYMMDDTHHMMSS>, and DATE-OBS); nothing was marked"
        )
    if not rows:
        raise SystemExit("Nothing calibrated: no night had both science and calibrator data")
    if not any(r["cp_only"] for r in rows):  # the column appears only when it matters
        for r in rows:
            del r["cp_only"]
    with open(os.path.join(args.out, "index.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    insnames = exposures[0].insnames
    for ins in insnames:
        try:
            columns = flux_ratio(
                exposures,
                ins,
                window_um=window,
                nbin=args.bin,
                norm_um=tuple(args.flux_norm) if args.flux_norm else None,
            )
        except KeyError:  # no OI_FLUX table
            continue
        if set(columns) == {"wavelength_um", "mean", "std"}:  # no night with science and calibrator
            continue
        suffix = "" if len(insnames) == 1 else f"_{ins}"
        with open(os.path.join(args.out, f"flux_ratio{suffix}.csv"), "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(columns.keys())
            w.writerows(np.round(np.array(list(columns.values())).T, 6).tolist())


if __name__ == "__main__":
    main()
