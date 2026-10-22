"""Calibration residuals of calibrators, and their principal components.

The plan is design/gravity_calibrator_pca.md. What a science fit suffers from
is not the transfer function itself but what calibration leaves behind. So
each calibrator is calibrated, as if it were a science target, by the other
calibrators of its night (same setup and stations, another star), and
compared with its own uniform-disk model:

* V²: r = log(V²_calibrated / V²_UD) per (frame, baseline, channel), with
  σ_r = σ(V²_calibrated)/V²_calibrated;
* closure phases: r = T3PHI_calibrated (degrees; a uniform disk's is zero).

These residuals hold the noise and every calibration systematic the
transfer function does not remove. A **weighted PCA** (alternating weighted
least squares, weights σ⁻², after Bailey 2012) finds their dominant spectral
shapes in physical units (fractional V², degrees). Cross-validation by night
says how many shapes are real: a component counts only if it lowers the
held-out χ² by more than fitting noise would. The amplitudes' spread, minus
their noise, gives each shape's width: the templates and priors virgil's
correlated nuisances (its Stage 6d) need. Finally, the per-frame amplitudes
of the six baselines are tested for telescope-based structure (a_kl =
t_k + t_l), as Kammerer et al. 2020 found for GRAVITY's noise.

σ holds only statistical noise: the calibrator's own reported errors and the
other calibrators', carried through the interpolation
(``error="propagated"``). Anything beyond that is a calibration systematic,
which is what this measures.
"""

import dataclasses
import warnings

import numpy as np
from scipy.optimize import linear_sum_assignment

from . import gravity
from .calibrate import calibrate_exposure, uniform_disk_v2
from .products import normalise_name
from .viscal import setup_string


class NoResiduals(ValueError):
    """No calibrator has company to calibrate it: there is nothing to score."""


@dataclasses.dataclass
class Residuals:
    """Calibration residuals, one row per (frame, baseline) or (frame, triangle).

    ``r`` and ``sigma`` are (rows, channels); ``valid`` marks usable samples;
    ``meta`` holds one array per key, one value per row.
    """

    kind: str  # "vis2" or "t3phi"
    wave: np.ndarray
    r: np.ndarray
    sigma: np.ndarray
    valid: np.ndarray
    meta: dict

    def select(self, rows):
        return Residuals(self.kind, self.wave, self.r[rows], self.sigma[rows], self.valid[rows],
                         {k: v[rows] for k, v in self.meta.items()})  # fmt: skip


def _group_key(e):
    """Exposures that may calibrate each other: one night, setup, stations and field.

    The field configuration (:func:`~virgil_vlti.gravity.field_mode`) is part
    of the key, so single- and dual-field products, or dual-field pointings
    on and off the fringe-tracking star, never calibrate each other.
    """
    h = e.header
    stations = tuple(h.get(f"ESO ISS CONF STATION{k}") for k in range(1, 5))
    return e.night, setup_string(h), stations, gravity.field_mode(h)


# Which exposures of a group calibrate each one (see calibrator_residuals).
REFERENCES = ("other-star", "same-star", "any")


def _references(e, group, references):
    star = normalise_name(e.target)
    if references == "other-star":
        return [c for c in group if normalise_name(c.target) != star]
    if references == "same-star":
        return [c for c in group if c is not e and normalise_name(c.target) == star]
    return [c for c in group if c is not e]


def _conditions(h):
    def mean(key):
        vals = [h.get(f"{key} {end}") for end in ("START", "END")]
        vals = [float(v) for v in vals if v is not None and float(v) > 0]
        return float(np.mean(vals)) if vals else np.nan

    return {
        "seeing": mean("ESO ISS AMBI FWHM"),
        "tau0": mean("ESO ISS AMBI TAU0"),
        "airmass": mean("ESO ISS AIRM"),
        "telescope": "UT"
        if str(h.get("ESO ISS CONF T1NAME", "AT")).upper().startswith("UT")
        else "AT",
    }


def _on_grid(values, wave, grid):
    """Interpolate rows onto ``grid`` in 1/λ (or return them, if the grids agree)."""
    if wave.shape == grid.shape and np.allclose(wave, grid, rtol=1e-7):
        return values
    x, xg = 1 / wave, 1 / grid
    order = np.argsort(x)
    return np.array(
        [np.interp(xg, x[order], row[order], left=np.nan, right=np.nan) for row in values]
    )


def calibrator_residuals(
    exposures,
    insname="GRAVITY_SC",
    grid=None,
    *,
    exclude=(),
    window_um=None,
    field=None,
    references="other-star",
    **calibrate_kwargs,
):
    """Leave-one-out calibration residuals of every calibrator with company.

    Parameters
    ----------
    exposures : list of Exposure
        Calibrator products (``*_CAL_VIS``) with ``pipeline_diameter``s.
        They are grouped by night, setup and stations; each is calibrated by
        the others of its group that observed another star.
    insname : str
    grid : array, optional
        Wavelengths to put every residual on (default: the first file's).
    exclude : iterable of str
        Calibrator stars to leave out entirely, as targets and as references
        (names matched as in :func:`~virgil_vlti.products.classify`).
    window_um : (float, float), optional
        Keep only residual channels inside this window (µm).
    field : str, optional
        Keep only products in this field configuration
        (:func:`~virgil_vlti.gravity.field_matches`: ``"SINGLE"``,
        ``"DUAL"``, ``"DUAL_ON"``, ``"DUAL_OFF"`` or ``"DUAL_WIDE"``). Whatever
        is kept, configurations never calibrate each other.
    references : {"other-star", "same-star", "any"}
        Which exposures of the group calibrate each one. ``"other-star"``
        (default): the night's other stars, as a science target is
        calibrated, so the residuals hold everything calibration leaves
        behind. ``"same-star"``: the star's own other exposures, so diameter
        errors and anything else constant over the night cancel, and what is
        left is the exposure-to-exposure variation (as in the pair test,
        :mod:`virgil_vlti.pairs`). Dual-field calibrators are rarely more than
        one star a night, so this is often the only choice there.
    **calibrate_kwargs
        For :func:`~virgil_vlti.calibrate.calibrate_exposure`; the default is
        Gaussian time weighting (Δ = 1 h) and statistical errors only
        (``error="propagated"``). Channels are not binned: ``nbin`` other
        than 1 is refused, as the residuals are taken against the
        full-resolution uniform-disk model and put on ``grid``. Header
        keys missing on either side are refused; pass ``allow_unrecorded=True``
        to warn instead.

    Returns
    -------
    dict
        ``{"vis2": Residuals, "t3phi": Residuals}``.
    """
    if references not in REFERENCES:
        raise ValueError(f"references must be one of {REFERENCES}, not {references!r}")
    if calibrate_kwargs.get("nbin", 1) != 1:
        raise ValueError(
            "calibrator_residuals works at full spectral resolution: "
            f"nbin={calibrate_kwargs['nbin']!r} is not supported"
        )
    kw = {"method": "gaussian", "error": "propagated", "strict": False, **calibrate_kwargs}
    dropped = {normalise_name(t) for t in exclude}
    groups = {}
    for e in exposures:
        if e.pipeline_diameter is None or insname not in e.wave:
            continue
        if normalise_name(e.target) in dropped:
            continue
        if field is not None and not gravity.field_matches(field, gravity.field_mode(e.header)):
            continue
        groups.setdefault(_group_key(e), []).append(e)
    if not groups:
        where = "" if field is None else f" in field {field}"
        raise NoResiduals(
            f"No calibrator product with a pipeline diameter and {insname} tables{where} "
            "among the inputs"
        )
    if grid is None:
        grid = next(iter(groups.values()))[0].wave[insname]
    out = {"vis2": [], "t3phi": []}
    for key, group in groups.items():
        for e in group:
            others = _references(e, group, references)
            if not others:
                continue
            diameters = {c.target: c.pipeline_diameter[0] for c in group}
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                cal = calibrate_exposure(e, others, insname, diameters=diameters, **kw)
            wave = e.wave[insname]
            cond = _conditions(e.header)
            common = {"file": str(e.path), "night": key[0], "target": e.target,
                      "n_cal": len(others), "field": str(key[3]),
                      "separation_mas": _separation(e.header), **cond}  # fmt: skip
            # V²: log ratio to the calibrator's own uniform-disk model.
            v2 = e.v2[insname]
            model = uniform_disk_v2(v2, wave, e.pipeline_diameter[0])
            c = cal["vis2"]
            ok = ~c.flag & (c.value > 0) & (c.error > 0)
            with np.errstate(invalid="ignore", divide="ignore"):
                r = np.where(ok, np.log(c.value / model), np.nan)
                s = np.where(ok, c.error / c.value, np.nan)
            rows = _rows_meta(v2, common, frames=True)
            out["vis2"].append((_on_grid(r, wave, grid), _on_grid(s, wave, grid), rows))
            t3 = e.t3[insname]
            c = cal["t3phi"]
            ok = ~c.flag & (c.error > 0)
            r = np.where(ok, c.value, np.nan)
            s = np.where(ok, c.error, np.nan)
            out["t3phi"].append(
                (
                    _on_grid(r, wave, grid),
                    _on_grid(s, wave, grid),
                    _rows_meta(t3, common, frames=True),
                )
            )
    result = {}
    for kind, parts in out.items():
        if not parts:
            raise NoResiduals(
                {"other-star": "No night has two calibrators of different stars in one setup",
                 "same-star": "No star has two exposures on one night in one setup",
                 "any": "No night has two calibrator exposures in one setup"}[references]
            )  # fmt: skip
        r = np.vstack([p[0] for p in parts])
        s = np.vstack([p[1] for p in parts])
        meta = {k: np.concatenate([p[2][k] for p in parts]) for k in parts[0][2]}
        valid = np.isfinite(r) & np.isfinite(s) & (s > 0)
        if window_um is not None:
            grid_um = np.asarray(grid) * 1e6
            valid &= ((grid_um > window_um[0]) & (grid_um < window_um[1]))[None, :]
        result[kind] = Residuals(kind, np.asarray(grid), np.where(valid, r, 0.0),
                                 np.where(valid, s, np.inf), valid, meta)  # fmt: skip
    return result


def _separation(header):
    """The science fibre's offset from the fringe tracker (mas), nan if not recorded."""
    try:
        return gravity.fibre_separation_mas(header)
    except (KeyError, TypeError, ValueError):
        return np.nan


def _rows_meta(table, common, frames):
    d = table.data
    stations = [tuple(int(x) for x in s) for s in d["STA_INDEX"]]
    count, frame = {}, []
    for s in stations:  # the k-th row of a station tuple is frame k
        frame.append(count.get(s, 0))
        count[s] = frame[-1] + 1
    n = len(stations)
    meta = {k: np.array([v] * n) for k, v in common.items()}
    meta["stations"] = np.array(["-".join(map(str, s)) for s in stations])
    meta["frame"] = np.array(frame)
    if "UCOORD" in d.names:
        meta["baseline_m"] = np.hypot(d["UCOORD"], d["VCOORD"])
    return meta


def target_scores(res):
    """Per calibrator star: the median over its rows of the rms whitened residual."""
    z2 = np.where(res.valid, (res.r / res.sigma) ** 2, 0.0)
    n = res.valid.sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        rms = np.sqrt(z2.sum(1) / n)
    names = np.array([normalise_name(t) for t in res.meta["target"]])
    out = {}
    for t in np.unique(names):
        rows = rms[(names == t) & (n > 0)]
        out[t] = float(np.median(rows)) if rows.size else np.nan  # nan: nothing usable
    return out


def screen_calibrators(exposures, insname="GRAVITY_SC", *, factor=5.0, max_rounds=50, **kwargs):
    """Leave out calibrators that do not behave like their uniform disks, one at a time.

    A binary or a badly modelled star spoils its own residuals and those of
    every star it calibrates. Each round scores every star
    (:func:`target_scores`) in V² and in closure phase, relative to the
    median star; if the worst exceeds ``factor`` in either, it is excluded
    (as target and as reference) and the residuals are recomputed.

    Screening uses every night, so it is a quality cut made before any
    cross-validation: the held-out χ² and widths that follow are conditional
    on it (an excluded star's held-out nights also helped exclude it).

    When no star has company to be scored against (no night has two
    calibrator stars in one setup, or no calibrator product at all), nothing
    can be screened: a warning is issued and ``({}, [])`` is returned. It never
    raises for that reason, so a calibrator-less set of nights is "not
    screened", not an error.

    Returns
    -------
    residuals : dict
        As :func:`calibrator_residuals`, without the excluded stars; empty
        when nothing could be screened.
    excluded : list of dict
        Each excluded star, its relative scores and the round.
    """
    excluded = []

    def residuals():
        return calibrator_residuals(
            exposures, insname, exclude=[e["target"] for e in excluded], **kwargs
        )

    for round_ in range(max_rounds):
        try:
            res = residuals()
        except NoResiduals as err:
            warnings.warn(f"calibrators not screened: {err}", stacklevel=2)
            return {}, excluded
        rel = {}
        for kind in ("vis2", "t3phi"):
            scores = {t: v for t, v in target_scores(res[kind]).items() if np.isfinite(v)}
            ref = np.median(list(scores.values())) if scores else np.nan
            if not np.isfinite(ref) or ref <= 0:
                continue  # no usable scores for this observable
            for t, v in scores.items():
                rel.setdefault(t, {})[kind] = v / ref
        if not rel:
            warnings.warn(
                "calibrators not screened: no star has usable V² or closure-phase residuals "
                "to score",
                stacklevel=2,
            )
            return res, excluded
        worst = max(rel, key=lambda t: max(rel[t].values()))
        if max(rel[worst].values()) <= factor:
            return res, excluded
        targets = res["vis2"].meta["target"]
        name = next(t for t in targets if normalise_name(t) == worst)  # as the file names it
        scores = {f"{k}_score": v for k, v in rel[worst].items()}
        excluded.append({"target": str(name), "round": round_, **scores})
    try:
        return residuals(), excluded  # with every exclusion applied
    except NoResiduals as err:  # the exclusions left no star with company
        warnings.warn(f"calibrators not screened: {err}", stacklevel=2)
        return {}, excluded


# ---------------------------------------------------------------------------
# Weighted PCA
# ---------------------------------------------------------------------------


@dataclasses.dataclass
class Components:
    """Weighted principal components of residual spectra.

    ``components`` (k, channels) are orthonormal shapes in the residuals'
    units; ``amplitudes`` (rows, k) and their noise variances
    ``amplitude_var`` (rows, k) come from weighted least squares per row.
    """

    components: np.ndarray
    amplitudes: np.ndarray
    amplitude_var: np.ndarray

    @property
    def widths(self):
        """Each component's intrinsic rms amplitude: total spread minus noise."""
        excess = np.mean(self.amplitudes**2, 0) - np.mean(self.amplitude_var, 0)
        return np.sqrt(np.clip(excess, 0, None))


def _amplitudes(r, w, C):
    """Weighted least-squares amplitudes of every row on components C (k, n)."""
    k = C.shape[0]
    A = np.zeros((r.shape[0], k))
    V = np.full((r.shape[0], k), np.inf)
    for i in range(r.shape[0]):
        M = (C * w[i]) @ C.T
        if np.linalg.matrix_rank(M) < k:
            continue
        Minv = np.linalg.inv(M)
        A[i] = Minv @ (C * w[i]) @ r[i]
        V[i] = np.diag(Minv)
    return A, V


def weighted_pca(res, k, iterations=50, seed=0):
    """Weighted PCA of residual spectra by alternating weighted least squares.

    Weights are σ⁻² (zero where invalid). With uniform weights this reduces
    to an ordinary (uncentred) PCA; the mean residual, if any, is a component
    like any other.
    """
    r = np.where(res.valid, res.r, 0.0)
    w = np.where(res.valid, res.sigma**-2.0, 0.0)
    n = r.shape[1]
    # Start from the SVD of the weighted residuals.
    _, _, vt = np.linalg.svd(r * np.sqrt(w), full_matrices=False)
    C = vt[:k] if k else np.zeros((0, n))
    if k == 0:
        return Components(C, np.zeros((r.shape[0], 0)), np.zeros((r.shape[0], 0)))
    for _ in range(iterations):
        A, _ = _amplitudes(r, w, C)
        # Components: per channel, weighted least squares on the amplitudes.
        for c in range(n):
            M = (A * w[:, c : c + 1]).T @ A
            if np.linalg.matrix_rank(M) == k:
                C[:, c] = np.linalg.solve(M, (A * w[:, c : c + 1]).T @ r[:, c])
        q, _ = np.linalg.qr(C.T)
        C = q.T
    # Order by the variance each component carries.
    A, V = _amplitudes(r, w, C)
    order = np.argsort(-np.mean(np.where(np.isfinite(V), A**2, 0), 0))
    C, A, V = C[order], A[:, order], V[:, order]
    signs = np.sign(C[np.arange(len(C)), np.argmax(np.abs(C), 1)])  # a fixed sign
    return Components(C * signs[:, None], A * signs, V)


def _mad(x, axis=0):
    """The normal-equivalent median absolute deviation (1.4826 × MAD)."""
    return 1.4826 * np.median(np.abs(x - np.median(x, axis=axis)), axis=axis)


def clip_rows(res, k=2, threshold=5.0, iterations=5):
    """Drop outlier rows (frame–baseline or frame–triangle spectra), iteratively.

    A least-squares PCA is pulled hardest by its worst rows, and a handful of
    them (a bad frame, a near-zero V² whose log diverges, a phase jump) can
    dominate both the components and an rms width. Each iteration fits ``k``
    components, then drops rows whose amplitude on any component, or whose
    rms whitened residual after the fit, lies more than ``threshold`` robust
    deviations (:func:`_mad`) from the median; until nothing changes.

    Returns
    -------
    Residuals
        The kept rows.
    numpy.ndarray
        Indices (into ``res``) of the dropped rows.
    """
    keep = np.arange(res.r.shape[0])
    for _ in range(iterations):
        sub = res.select(keep)
        pca = weighted_pca(sub, k)
        A, V = pca.amplitudes, pca.amplitude_var
        ok = np.all(np.isfinite(V), 1)
        bad = ~ok
        for j in range(A.shape[1]):
            a = A[ok, j]
            dev = np.abs(A[:, j] - np.median(a)) / max(_mad(a), 1e-12)
            bad |= ok & (dev > threshold)
        c2, n = chi2_after(sub, pca.components)
        with np.errstate(invalid="ignore", divide="ignore"):
            rms = np.sqrt(c2 / n)
        good = np.isfinite(rms)
        dev = np.abs(rms - np.median(rms[good])) / max(_mad(rms[good]), 1e-12)
        bad |= ~good | (dev > threshold)
        if not bad.any():
            break
        keep = keep[~bad]
    dropped = np.setdiff1d(np.arange(res.r.shape[0]), keep)
    return res.select(keep), dropped


def chi2_after(res, C):
    """χ² per valid sample of each row after removing its best fit on C."""
    r = np.where(res.valid, res.r, 0.0)
    w = np.where(res.valid, res.sigma**-2.0, 0.0)
    A, _ = _amplitudes(r, w, C) if len(C) else (np.zeros((r.shape[0], 0)), None)
    model = A @ C if len(C) else 0.0
    return np.sum(w * (r - model) ** 2, 1), res.valid.sum(1)


def cross_validate(res, k_max, folds=5, seed=0):
    """Held-out χ² per sample against the number of components, by night.

    Nights are split into ``folds``; components fitted on the others are
    fitted to each held-out night's rows. Returns ``(k, chi2_per_sample,
    expected_for_noise)``: fitting k components to pure noise lowers χ² per
    sample by about k / (channels per row), so a component is real when
    the held-out drop exceeds that.
    """
    fold = _folds(res, folds, seed)
    folds = fold.max() + 1
    out = []
    for k in range(k_max + 1):
        chi2, dof = 0.0, 0
        for f in range(folds):
            train, test = res.select(fold != f), res.select(fold == f)
            if not len(test.r):
                continue
            C = weighted_pca(train, k).components
            c2, n = chi2_after(test, C)
            chi2 += c2.sum()
            dof += n.sum()
        per_row = res.valid.sum(1).mean()
        out.append((k, chi2 / dof, (1 - k / per_row)))
    return out


def _folds(res, folds, seed):
    """A fold per row, splitting nights; at most one fold per night."""
    nights = np.unique(res.meta["night"])
    if len(nights) < 2:
        raise ValueError("Cross-validation by night needs residuals from at least two nights")
    folds = min(folds, len(nights))
    rng = np.random.default_rng(seed)
    fold_of = dict(zip(rng.permutation(nights), np.arange(len(nights)) % folds, strict=True))
    return np.array([fold_of[n] for n in res.meta["night"]])


def heldout_widths(res, k, folds=5, seed=0, robust=True):
    """Component widths measured on held-out nights, free of selection bias.

    In-sample, the PCA picks the directions of largest variance, so even pure
    noise gets a positive "width". Here components are fitted without each
    fold's nights and that fold's rows are projected on them (matched to the
    full fit's components by overlap, sign included); the width is the
    held-out amplitudes' spread minus their noise, zero if negative. With
    ``robust`` (default) the spread is the normal-equivalent MAD about the
    median, and the noise the median noise variance; otherwise the rms. Its own noise floor
    is ≈ σ √(2 · δ), where δ is the scatter of the variance ratio for a fixed
    direction (≈ √(2/rows)): widths below that are not detected.
    """
    full = weighted_pca(res, k).components
    fold = _folds(res, folds, seed)
    folds = fold.max() + 1
    amps, noise = [], []
    r = np.where(res.valid, res.r, 0.0)
    w = np.where(res.valid, res.sigma**-2.0, 0.0)
    for f in range(folds):
        test = fold == f
        if not test.any():
            continue
        C = weighted_pca(res.select(~test), k).components
        # One-to-one matching to the full fit's components, by |overlap|.
        _, match = linear_sum_assignment(-np.abs(full @ C.T))
        C = C[match] * np.sign(np.sum(full * C[match], 1))[:, None]
        A, V = _amplitudes(r[test], w[test], C)
        ok = np.all(np.isfinite(V), 1)
        amps.append(A[ok])
        noise.append(V[ok])
    A, V = np.vstack(amps), np.vstack(noise)
    if robust:
        spread2 = _mad(A) ** 2  # insensitive to the few rows that dominate an rms
        return np.sqrt(np.clip(spread2 - np.median(V, 0), 0, None))
    return np.sqrt(np.clip(np.mean(A**2, 0) - np.mean(V, 0), 0, None))


def telescope_fraction(res, amplitudes):
    """How much of each component's per-frame baseline amplitudes is telescope-based.

    For each (file, frame) with all six baselines of four telescopes, fit
    a_kl = t_k + t_l (four parameters for six numbers, not centred, so a
    common gain counts as telescope-based) and return the fraction of the
    amplitudes' sum of squares it explains, per component. With no telescope
    structure the expectation is 4/6 ≈ 0.67; 1 means purely telescope-based.
    """
    keys = list(zip(res.meta["file"], res.meta["frame"], strict=True))
    groups = {}
    for i, key in enumerate(keys):
        groups.setdefault(key, []).append(i)
    explained, total = np.zeros(amplitudes.shape[1]), np.zeros(amplitudes.shape[1])
    for rows in groups.values():
        stations = [tuple(res.meta["stations"][i].split("-")) for i in rows]
        tels = sorted({t for s in stations for t in s})
        if len(rows) != 6 or len(tels) != 4:
            continue
        X = np.array([[t in s for t in tels] for s in stations], float)  # rank 4
        a = amplitudes[rows]  # not centred: a common gain is telescope-based too
        coef, *_ = np.linalg.lstsq(X, a, rcond=None)
        fit = X @ coef
        explained += np.sum(fit**2, 0)
        total += np.sum(a**2, 0)
    with np.errstate(invalid="ignore"):
        return explained / total


def condition_correlations(res, amplitudes, keys=("seeing", "tau0", "airmass")):
    """Spearman-like rank correlation of |amplitude| with each observing condition."""
    out = {}
    for key in keys:
        x = np.asarray(res.meta[key], float)
        ok = np.isfinite(x)
        if ok.sum() < 10:
            continue
        rx = np.argsort(np.argsort(x[ok]))
        out[key] = [
            float(np.corrcoef(rx, np.argsort(np.argsort(np.abs(a[ok]))))[0, 1])
            for a in amplitudes.T
        ]
    return out


__all__ = [
    "Components",
    "NoResiduals",
    "Residuals",
    "calibrator_residuals",
    "chi2_after",
    "condition_correlations",
    "clip_rows",
    "cross_validate",
    "heldout_widths",
    "screen_calibrators",
    "target_scores",
    "telescope_fraction",
    "weighted_pca",
]
