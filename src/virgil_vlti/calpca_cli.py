"""``virgil-vlti-calpca``: the calibrator PCA on a folder of calibrator products.

    virgil-vlti-calpca reduced/ --out pca/ --k-max 6
    virgil-vlti-calpca reduced/ --out pca_dual/ --field DUAL_OFF --references same-star

Builds leave-one-out calibration residuals (:mod:`virgil_vlti.calpca`), fits
weighted PCAs to V² and closure phases, cross-validates the number of
components by night, measures held-out widths, the telescope fraction and the
dependence on seeing, τ₀ and airmass, and writes

* ``<out>/calpca.npz``: wavelengths, components and widths per observable;
* ``<out>/summary.json``: everything else, what went in, and the calibrators
  excluded by screening (:func:`~virgil_vlti.calpca.screen_calibrators`);
* ``<out>/residuals.npz``: the residuals, errors, masks and per-row metadata;
* ``<out>/components.png`` if matplotlib is available.
"""

import argparse
import json
import os
import time

import numpy as np

from .calpca import (
    REFERENCES,
    Residuals,
    calibrator_residuals,
    clip_rows,
    condition_correlations,
    cross_validate,
    heldout_widths,
    screen_calibrators,
    target_scores,
    telescope_fraction,
    weighted_pca,
)
from .products import expand_paths, instrument_of, read_exposure


def analyse(res, k_max, folds, min_gain=0.02, clip=5.0):
    """The PCA, its cross-validation and diagnostics, for one observable.

    Outlier rows are clipped first (:func:`~virgil_vlti.calpca.clip_rows`,
    with two components, at ``clip`` robust deviations; 0 turns it off). A
    component is kept if it lowers the held-out χ² by more than fitting noise
    would, plus ``min_gain`` of what remains.
    """
    n_rows = res.r.shape[0]
    dropped = np.array([], int)
    if clip > 0:
        res, dropped = clip_rows(res, k=2, threshold=clip)
    cv = cross_validate(res, k_max, folds=folds)
    k = 0
    for (_, c0, e0), (k1, c1, e1) in zip(cv, cv[1:], strict=False):
        if c1 / c0 < e1 / e0 - min_gain:
            k = k1
        else:
            break
    per_row = res.valid.sum(1).mean()
    info = {
        "rows_in": int(n_rows),
        "rows_clipped": int(len(dropped)),
        "rows": int(res.r.shape[0]),
        "nights": int(len(np.unique(res.meta["night"]))),
        "files": int(len(np.unique(res.meta["file"]))),
        "chi2_per_sample_no_components": float(cv[0][1]),
        "cross_validation": [
            {"k": int(a), "heldout_chi2_per_sample": float(b), "noise_only_ratio": float(c)}
            for a, b, c in cv
        ],
        "components_kept": int(k),
    }
    if k == 0:  # nothing beats noise: report no components
        empty = np.zeros((0, res.r.shape[1]))
        info.update(heldout_widths=[], heldout_rms_per_channel=[], telescope_fraction=None,
                    condition_rank_correlations={})  # fmt: skip
        return empty, info
    pca = weighted_pca(res, k)
    widths = heldout_widths(res, k, folds=folds)
    info.update(
        heldout_widths=widths.tolist(),
        heldout_rms_per_channel=(widths / np.sqrt(per_row)).tolist(),
        telescope_fraction=telescope_fraction(res, pca.amplitudes).tolist()
        if res.kind == "vis2"
        else None,
        condition_rank_correlations=condition_correlations(res, pca.amplitudes),
    )
    return pca.components, info


def load_residuals(path):
    """Residuals saved by an earlier run (``residuals.npz``), as ``{kind: Residuals}``."""
    d = np.load(path, allow_pickle=False)
    out = {}
    for kind in ("vis2", "t3phi"):
        prefix = f"{kind}_meta_"
        meta = {f[len(prefix) :]: d[f] for f in d.files if f.startswith(prefix)}
        out[kind] = Residuals(kind, d["wave"], d[f"{kind}_r"], d[f"{kind}_sigma"],
                              d[f"{kind}_valid"], meta)  # fmt: skip
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument(
        "inputs",
        nargs="+",
        help="calibrator products, folders of *.fits, or one residuals.npz to re-analyse",
    )
    p.add_argument("--out", required=True)
    p.add_argument("--insname", default="GRAVITY_SC")
    p.add_argument(
        "--field",
        default=None,
        type=str.upper,
        choices=("SINGLE", "DUAL", "DUAL_ON", "DUAL_OFF", "DUAL_WIDE"),
        help="keep only products in this field configuration (virgil_vlti.gravity.field_mode)",
    )
    p.add_argument(
        "--references",
        default="other-star",
        choices=REFERENCES,
        help="which exposures calibrate each calibrator: the night's other stars (default), "
        "the star's own other exposures, or any",
    )
    p.add_argument(
        "--allow-unrecorded",
        action="store_true",
        help="warn, rather than stop, when a header key or pipeline parameter that must match "
        "between a calibrator and its references is missing",
    )
    p.add_argument("--k-max", type=int, default=6)
    p.add_argument("--folds", type=int, default=5)
    p.add_argument(
        "--window",
        nargs=2,
        type=float,
        default=(2.0, 2.45),
        metavar=("LO", "HI"),
        help="µm; GRAVITY's band edges are noisy (default 2.0 2.45)",
    )
    p.add_argument(
        "--min-gain",
        type=float,
        default=0.02,
        help="keep a component if it lowers the held-out χ² by this fraction beyond noise",
    )
    p.add_argument(
        "--clip",
        type=float,
        default=5.0,
        help="clip outlier rows beyond this many robust deviations (0: off)",
    )
    p.add_argument(
        "--screen-factor",
        type=float,
        default=5.0,
        help="exclude calibrators scoring above this × the median star (0: off)",
    )
    args = p.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    start = time.time()
    window = tuple(args.window)
    if len(args.inputs) == 1 and args.inputs[0].endswith(".npz"):
        # Re-analyse saved residuals (from an earlier run's residuals.npz).
        residuals = load_residuals(args.inputs[0])
        excluded = []
        summary = {"from_residuals": args.inputs[0]}
        print(f"residuals loaded from {args.inputs[0]}", flush=True)
    else:
        paths = expand_paths(args.inputs)
        exposures = []
        for path in paths:
            try:
                e = read_exposure(path, args.insname)
                if instrument_of(e.header) == "MATISSE":
                    # needs calpca(group_key=matisse.group_key) (design/matisse_plan.md)
                    raise ValueError("MATISSE calpca is Stage 5")
                exposures.append(e)
            except (KeyError, ValueError) as err:
                print(f"skipped {os.path.basename(path)}: {err}")
        select = {
            "window_um": window,
            "field": args.field,
            "references": args.references,
            "allow_unrecorded": args.allow_unrecorded,
        }
        if args.screen_factor > 0:
            residuals, excluded = screen_calibrators(
                exposures, args.insname, factor=args.screen_factor, **select
            )
        else:
            residuals = calibrator_residuals(exposures, args.insname, **select)
            excluded = []
        if not residuals:
            raise SystemExit("No calibrator has company to be scored against; nothing to do")
        print(f"{len(exposures)} files read; residuals in {time.time() - start:.0f} s; "
              f"excluded {[e['target'] for e in excluded]}", flush=True)  # fmt: skip
        summary = {
            "inputs": len(paths),
            "read": len(exposures),
            "insname": args.insname,
            "window_um": list(window),
            "field": args.field,
            "references": args.references,
            "screen_factor": args.screen_factor,
            "excluded_calibrators": excluded,
        }
    summary["target_scores"] = {k: target_scores(r) for k, r in residuals.items()}
    # The residuals themselves, for diagnostics.
    np.savez_compressed(
        os.path.join(args.out, "residuals.npz"),
        **{
            f"{k}_{a}": getattr(r, a) for k, r in residuals.items() for a in ("r", "sigma", "valid")
        },
        **{f"{k}_meta_{m}": v for k, r in residuals.items() for m, v in r.meta.items()},
        wave=next(iter(residuals.values())).wave,
    )
    arrays = {}
    for kind, res in residuals.items():
        components, info = analyse(res, args.k_max, args.folds, args.min_gain, args.clip)
        summary[kind] = info
        arrays[f"{kind}_wave"] = res.wave
        arrays[f"{kind}_components"] = components
        arrays[f"{kind}_heldout_widths"] = np.asarray(info["heldout_widths"])
        print(f"{kind}: {json.dumps(info)}", flush=True)
    np.savez(os.path.join(args.out, "calpca.npz"), **arrays)
    summary["elapsed_s"] = round(time.time() - start, 1)
    with open(os.path.join(args.out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=float)
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
        for ax, kind, unit in zip(axes, ("vis2", "t3phi"), ("log V²", "deg"), strict=True):
            w = arrays[f"{kind}_heldout_widths"]
            for j, comp in enumerate(arrays[f"{kind}_components"]):
                ax.plot(arrays[f"{kind}_wave"] * 1e6, w[j] * comp, label=f"PC{j + 1}")
            ax.set(xlabel="wavelength (µm)", ylabel=f"1σ mode ({unit})", title=kind)
            if len(w):
                ax.legend(fontsize=8)
            else:
                ax.text(0.5, 0.5, "no component beats noise", ha="center", transform=ax.transAxes)
        fig.tight_layout()
        fig.savefig(os.path.join(args.out, "components.png"), dpi=120)
    except ImportError:
        print("matplotlib is not installed: components.png not written")
    print(f"done in {summary['elapsed_s']} s; {args.out}/summary.json")


if __name__ == "__main__":
    main()
