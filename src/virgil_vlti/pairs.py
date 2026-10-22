"""The pair test: the structure of the errors, from consecutive exposures.

Two frames of one target in one configuration, taken minutes apart, see
almost the same uv points, so their difference is almost pure noise. Fitting
the covariance of many such differences measures what the reported errors
leave out (design/pionier_systematics.md):

* **V²**: independent noise, ``white`` times the reported σ in each
  channel; a gain per (frame, baseline), shared by all channels, of
  fractional rms ``baseline``; and a gain per (frame, telescope), which
  enters the three baselines through that telescope, of rms ``telescope``.
* **Closure phases**: independent noise, ``white`` times the reported σ,
  with the ±1/3 correlation between triangles that share a baseline that
  baseline-phase noise implies (as in virgil since its #120; only
  independent triangles are kept); and an offset per (frame, triangle),
  shared by all channels, of rms ``offset`` (degrees).

A frame is one exposure's set of baselines: a PIONIER product usually holds
one, a GRAVITY product several.

Nested models are fitted by maximum likelihood, and their log-likelihood
gains say which terms the data need. Consecutive exposures usually share
calibrators, so calibration errors common to both cancel: the correlated
amplitudes are **lower bounds**.

If the target's visibilities change appreciably between the two exposures
(a resolved target and rotating baselines), subtract a model's prediction
with ``model``; otherwise the change is mistaken for noise.
"""

import dataclasses
import itertools

import numpy as np
from scipy.optimize import minimize

V2_MODELS = {
    "independent only": ("white",),
    "+ baseline": ("white", "baseline"),
    "+ telescope": ("white", "telescope"),
    "+ baseline + telescope": ("white", "baseline", "telescope"),
}
CP_MODELS = {
    "independent only": ("white",),
    "+ triangle offset": ("white", "offset"),
}


@dataclasses.dataclass
class Block:
    """One pair's differences for one observable, and its covariance terms.

    ``terms[name]`` is the covariance that term contributes at unit width,
    so the total is Σ width² · terms[name].
    """

    diff: np.ndarray
    terms: dict


@dataclasses.dataclass
class Frame:
    """One frame of one exposure: its rows in ``OI_VIS2`` and ``OI_T3`` by stations."""

    exposure: object
    insname: str
    mjd: float
    v2_rows: dict  # station pair -> row
    t3_rows: dict  # station triple -> row
    span: tuple = (np.nan, np.nan)  # MJD range that the frame's rows cover, ± the gap


def _ins(exposure, insname):
    """``insname``, or the exposure's only instrument when it is None."""
    if insname is not None:
        return insname
    if len(exposure.insnames) != 1:
        raise ValueError(f"{exposure.path} has several instruments; choose insname")
    return exposure.insnames[0]


FRAME_GAP_S = 60.0  # rows further apart in time than this belong to different frames


def frames_of(exposure, insname, frame_gap_s=FRAME_GAP_S):
    """Split an exposure into frames by time.

    ``OI_VIS2`` rows are sorted by MJD and a new frame starts wherever two
    consecutive rows are more than ``frame_gap_s`` apart (GRAVITY stamps the
    baselines of one frame up to ~12 s apart; frames are minutes apart). A
    frame that holds a baseline twice is malformed and dropped. Each
    ``OI_T3`` row joins the frame whose time span (± ``frame_gap_s``) holds
    its MJD. A PIONIER product usually holds one frame; a GRAVITY product
    holds several.
    """
    insname = _ins(exposure, insname)
    gap = frame_gap_s / 86400
    v2, t3 = exposure.v2[insname].data, exposure.t3[insname].data
    order = np.argsort(v2["MJD"], kind="stable")
    groups, last = [], None
    for row in order:
        t = v2["MJD"][row]
        if last is None or t - last > gap:
            groups.append([])
        groups[-1].append(row)
        last = t
    out = []
    for rows in groups:
        stations = [tuple(v2["STA_INDEX"][r]) for r in rows]
        if len(set(stations)) != len(stations):
            continue  # a baseline twice in one frame: misaligned
        times = v2["MJD"][rows]
        span = (times.min() - gap, times.max() + gap)
        rows_by_sta = dict(zip(stations, rows, strict=True))
        out.append(Frame(exposure, insname, float(np.mean(times)), rows_by_sta, {}, span))
    doubled = set()
    for row, sta in enumerate(map(tuple, t3["STA_INDEX"])):
        frame = next((f for f in out if f.span[0] <= t3["MJD"][row] <= f.span[1]), None)
        if frame is None:
            continue
        if sta in frame.t3_rows:
            doubled.add(id(frame))  # a triangle twice: misaligned
        frame.t3_rows[sta] = row
    for f in out:
        if id(f) in doubled:
            f.t3_rows = {}
    return out


def _compatible(a, b, max_gap_min):
    """Whether two frames can be differenced: one target, close in time, one wavelength grid."""
    if a.exposure.target != b.exposure.target or (b.mjd - a.mjd) * 1440 > max_gap_min:
        return False
    wa, wb = a.exposure.wave[a.insname], b.exposure.wave[b.insname]
    return wa.shape == wb.shape and np.allclose(wa, wb, rtol=1e-4)


def find_pairs(exposures, insname=None, max_gap_min=30.0):
    """Disjoint pairs of consecutive, compatible frames, in time order.

    The pairs do not overlap: (A, B) and (B, C) would share B's noise, and
    their differences would not be independent, so their likelihoods could
    not simply be added.
    """
    frames = sorted((f for e in exposures for f in frames_of(e, insname)), key=lambda f: f.mjd)
    pairs, k = [], 0
    while k < len(frames) - 1:
        if _compatible(frames[k], frames[k + 1], max_gap_min):
            pairs.append((frames[k], frames[k + 1]))
            k += 2
        else:
            k += 1
    return pairs


def _matched(rows_a, rows_b):
    """Stations present in both frames, in a fixed order, and their rows."""
    common = sorted(set(rows_a) & set(rows_b))
    return common, [rows_a[s] for s in common], [rows_b[s] for s in common]


def _valid(table, rows, value, error):
    d = table.data
    err = d[error][rows]
    return ~d["FLAG"][rows] & np.isfinite(d[value][rows]) & np.isfinite(err) & (err > 0)


def v2_block(a, b, model=None):
    """The V² differences of one pair of frames, and their covariance terms."""
    ta, tb = a.exposure.v2[a.insname], b.exposure.v2[b.insname]
    stations, ia, ib = _matched(a.v2_rows, b.v2_rows)
    if not stations:
        return None
    va, vb = ta.data["VIS2DATA"][ia], tb.data["VIS2DATA"][ib]
    ea, eb = ta.data["VIS2ERR"][ia], tb.data["VIS2ERR"][ib]
    ok = _valid(ta, ia, "VIS2DATA", "VIS2ERR") & _valid(tb, ib, "VIS2DATA", "VIS2ERR")
    diff = va - vb
    ref = (va + vb) / 2
    if model is not None:
        pa = model(a.exposure, a.insname)["vis2"][ia]
        pb = model(b.exposure, b.insname)["vis2"][ib]
        diff = diff - (pa - pb)
        ref = (pa + pb) / 2
    base, chan = np.nonzero(ok)
    if base.size == 0:
        return None
    # A gain g in each frame changes the difference by (g_a − g_b)·V². The
    # covariance is built from a reference V² (the model, else the pair's
    # mean), not from each frame's noisy V²: otherwise the gain terms are
    # nearly full rank and absorb the white noise (and data-built covariances
    # bias multiplicative fits; Lachaume 2021).
    m = ref[base, chan]
    vv = 2 * np.outer(m, m)
    same = (base[:, None] == base[None, :]).astype(float)
    telescopes = sorted({t for s in stations for t in s})
    on = np.array([[t in stations[k] for k in base] for t in telescopes], float)
    shared = on.T @ on  # how many telescopes two samples' baselines share
    return Block(
        diff=diff[base, chan],
        terms={
            "white": np.diag(ea[base, chan] ** 2 + eb[base, chan] ** 2),
            "baseline": same * vv,
            "telescope": shared * vv,
        },
    )


def closure_matrix(triangles):
    """T, mapping baseline phases to closure phases, and its baselines.

    Triangles are (a, b, c) station tuples, with closure phase
    φ_ab + φ_bc − φ_ac; baselines are (low, high) station pairs, and a
    baseline taken the other way round enters with the opposite sign.
    """
    baselines = sorted({tuple(sorted(p)) for t in triangles for p in itertools.combinations(t, 2)})
    index = {b: k for k, b in enumerate(baselines)}
    T = np.zeros((len(triangles), len(baselines)))
    for r, (p, q, s) in enumerate(triangles):
        for (x, y), sign in (((p, q), 1), ((q, s), 1), ((p, s), -1)):
            T[r, index[tuple(sorted((x, y)))]] += sign if x < y else -sign
    return T, baselines


def closure_correlation(triangles):
    """R = T Tᵀ/3: with equal, independent baseline-phase noise, triangles that
    share a baseline correlate at ±1/3."""
    T, _ = closure_matrix(triangles)
    return T @ T.T / 3


def independent_triangles(triangles):
    """Indices of a maximal set of linearly independent triangles."""
    R = closure_correlation(triangles)
    keep = []
    for k in range(len(triangles)):
        trial = keep + [k]
        if np.linalg.matrix_rank(R[np.ix_(trial, trial)], tol=1e-9) == len(trial):
            keep.append(k)
    return keep


def cp_block(a, b, model=None):
    """The closure-phase differences of one pair of frames (degrees), and their terms.

    In each channel, a maximal set of linearly independent triangles is kept
    from those valid in both frames, so a flagged triangle is replaced by
    another where possible.
    """
    ta, tb = a.exposure.t3[a.insname], b.exposure.t3[b.insname]
    triangles, ia, ib = _matched(a.t3_rows, b.t3_rows)
    if not triangles:
        return None
    pa, pb = ta.data["T3PHI"][ia], tb.data["T3PHI"][ib]
    ea, eb = ta.data["T3PHIERR"][ia], tb.data["T3PHIERR"][ib]
    ok = _valid(ta, ia, "T3PHI", "T3PHIERR") & _valid(tb, ib, "T3PHI", "T3PHIERR")
    diff = pa - pb
    if model is not None:
        ma, mb = model(a.exposure, a.insname)["t3phi"], model(b.exposure, b.insname)["t3phi"]
        diff = diff - (ma[ia] - mb[ib])
    diff = (diff + 180) % 360 - 180
    tri, chan = [], []
    for c in range(ok.shape[1]):
        valid = np.flatnonzero(ok[:, c])
        for k in independent_triangles([triangles[v] for v in valid]):
            tri.append(valid[k])
            chan.append(c)
    if not tri:
        return None
    tri, chan = np.array(tri), np.array(chan)
    R = closure_correlation(triangles)[np.ix_(tri, tri)]
    same_chan = (chan[:, None] == chan[None, :]).astype(float)
    sa, sb = ea[tri, chan], eb[tri, chan]
    white = R * same_chan * (np.outer(sa, sa) + np.outer(sb, sb))
    same = (tri[:, None] == tri[None, :]).astype(float)
    return Block(diff=diff[tri, chan], terms={"white": white, "offset": 2.0 * same})


def _neg_loglike(log_widths, names, blocks):
    if np.any(np.abs(log_widths) > 30):  # keep the optimiser away from overflow
        return np.inf
    total = 0.0
    for block in blocks:
        cov = sum(np.exp(2 * w) * block.terms[n] for w, n in zip(log_widths, names, strict=True))
        try:
            L = np.linalg.cholesky(cov)
        except np.linalg.LinAlgError:
            return np.inf
        z = np.linalg.solve(L, block.diff)
        total += 0.5 * z @ z + np.log(np.diag(L)).sum()
    return total


@dataclasses.dataclass
class Fit:
    """One nested noise model: its widths and its log-likelihood gain."""

    widths: dict
    delta_lnl: float  # over the "independent only" model


def fit_noise(blocks, models, starts):
    """Fit each nested model by maximum likelihood.

    ``starts`` gives a starting width per term name.
    """
    out, base = {}, None
    for label, names in models.items():
        x0 = np.log([starts[n] for n in names])
        best = minimize(_neg_loglike, x0, args=(names, blocks), method="Nelder-Mead",
                        options={"xatol": 1e-4, "fatol": 1e-4, "maxiter": 4000})  # fmt: skip
        if base is None:
            base = best.fun
        out[label] = Fit(dict(zip(names, np.exp(best.x), strict=True)), float(base - best.fun))
    return out


@dataclasses.dataclass
class PairTest:
    """The result of :func:`pair_test`."""

    n_pairs: int  # disjoint pairs of consecutive frames
    v2: dict  # label -> Fit
    cp: dict  # label -> Fit

    def summary(self):
        lines = [f"{self.n_pairs} pairs"]
        for kind, fits in (("V²", self.v2), ("CP", self.cp)):
            for label, f in fits.items():
                widths = "  ".join(f"{k} {v:.4g}" for k, v in f.widths.items())
                lines.append(f"{kind} {label:24s} ΔlnL {f.delta_lnl:8.1f}   {widths}")
        return "\n".join(lines)


def pair_test(exposures, insname=None, *, max_gap_min=30.0, model=None):
    """Fit the error structure of consecutive exposures.

    Parameters
    ----------
    exposures : list of Exposure
        Calibrated products of one target (or several; pairs never mix
        targets).
    insname : str, optional
        The instrument; by default each exposure's only one (PIONIER, for
        instance, puts the wavelength solution in ``INSNAME``, so it changes
        from file to file).
    max_gap_min : float
        Longest gap between the two frames of a pair.
    model : callable, optional
        ``model(exposure, insname)`` returning ``{"vis2": ..., "t3phi": ...}``
        arrays shaped like that exposure's tables: the target's predicted
        values, whose change between the two exposures is subtracted.

    Returns
    -------
    PairTest
        V² widths: ``white`` (× σ), ``baseline`` and ``telescope``
        (fractional gains). Closure-phase widths: ``white`` (× σ) and
        ``offset`` (degrees).
    """
    pairs = find_pairs(exposures, insname, max_gap_min)
    v2 = [b for a, c in pairs if (b := v2_block(a, c, model)) is not None]
    cp = [b for a, c in pairs if (b := cp_block(a, c, model)) is not None]
    if not v2:
        raise ValueError("No pairs of consecutive exposures to compare")
    starts = {"white": 0.7, "baseline": 0.01, "telescope": 0.01, "offset": 0.3}
    return PairTest(
        n_pairs=len(pairs),
        v2=fit_noise(v2, V2_MODELS, starts),
        cp=fit_noise(cp, CP_MODELS, starts) if cp else {},
    )
