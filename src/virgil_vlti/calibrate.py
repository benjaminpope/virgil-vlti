"""Calibrate interferometric observables with a transfer function from calibrator exposures.

For each night, every calibrator (``CAL``) exposure gives one transfer-function
estimate per baseline (or triangle) and channel, averaged over the frames of
the file. Amplitudes (V², T3AMP and, optionally, VISAMP) are divided by the
calibrator's own model (a uniform disk, or a point source); phases (T3PHI and,
optionally, VISPHI) are averaged as unit vectors and estimate the
instrumental phase. These estimates are interpolated in time to each science
(``SCI``) frame: science amplitudes are divided by the transfer function and
the instrumental phases are subtracted.

The interpolation is either linear in time (clamped at the ends), or
Gaussian weights exp(−2 ΔT²/Δ²)/σ² with Δ = 1 h, the formula of the GRAVITY
pipeline manual (§9.20, eq. 79), with σ the estimate's median error floored at
2% (amplitudes) or 0.2° (phases). The floors are adapted from the 1.9.7 source
but applied to the error of the file mean, not the per-frame error, so they
bite sooner than in the pipeline. The source itself differs from the manual
(exp(−2|ΔT|/Δ)): ``--match-gravity-viscal`` reproduces the pipeline. Its
uncertainty enters the errors in quadrature, as one of:

* ``"loo"`` (default): the leave-one-out prediction error. Each calibrator
  estimate is predicted from the others with the same interpolator, and the
  r.m.s. error is used. This is the error of *predicting* the transfer
  function at another time, which is what matters for the science frames
  (design/gravity_calibration_review.md, §7);
* ``"scatter"``: the standard deviation of the calibrator estimates over the
  night.

The term is computed on the binned channels, from the binned calibrator
estimates, because it is common to neighbouring channels: averaging it in
quadrature would shrink it by up to √nbin (:func:`binned_estimates`). It is
then added independently to each binned channel, although it is common to all
channels of a frame and baseline, and to every science exposure of the night
on that baseline. With ``tf_modes=True`` it is instead returned as nuisance
modes, one per night and baseline (or triangle), for virgil's correlated
calibration nuisances (its Stage 6d): :func:`transfer_function_modes`.
Errors in the calibrators' diameters are not added to the errors: they are
common to everything a calibrator calibrates, and :func:`diameter_modes`
returns them as nuisance modes instead.

Channels can be restricted to a window and averaged in groups (``nbin``).
GRAVITY spectra are about twofold oversampled, so neighbouring channels are
correlated and binning by 2 or more is recommended when the correlations are
not modelled.
"""

import contextlib
import dataclasses
import os
import warnings

import numpy as np
from scipy.special import j1

from .gravity import NOT_RECORDED
from .products import instrument_module, instrument_of, normalise_name, read_flux

MAS = np.pi / 180 / 3600e3
PIPELINE_TIMESCALE_S = 3600.0


# ---------------------------------------------------------------------------
# Calibrator estimates
# ---------------------------------------------------------------------------


def uniform_disk_amplitude(u, v, wave, diameter_mas):
    """|V| of a uniform disk at baselines ``u, v`` (m, per row) and ``wave`` (m)."""
    shape = (np.size(u), wave.size)
    if not diameter_mas:
        return np.ones(shape)
    b = np.hypot(u, v)[:, None]
    x = np.pi * diameter_mas * MAS * b / wave[None, :]
    with np.errstate(invalid="ignore", divide="ignore"):
        amp = np.abs(np.where(x > 0, 2 * j1(x) / x, 1.0))
    if np.any(x > 3.0):
        warnings.warn(
            f"A calibrator of {diameter_mas} mas is near its first null "
            "on these baselines; it is a poor calibrator",
            stacklevel=3,
        )
    return amp


def uniform_disk_v2(table, wave, diameter_mas):
    """V² of a uniform disk at every row and channel of an ``OI_VIS2`` table."""
    return uniform_disk_visamp(table, wave, diameter_mas) ** 2


def uniform_disk_visamp(table, wave, diameter_mas):
    """|V| of a uniform disk at every row and channel of an ``OI_VIS`` table."""
    d = table.data
    return uniform_disk_amplitude(d["UCOORD"], d["VCOORD"], wave, diameter_mas)


def uniform_disk_t3amp(table, wave, diameter_mas):
    """|V₁V₂V₃| of a uniform disk at every row and channel of an ``OI_T3`` table."""
    d = table.data
    u1, v1, u2, v2 = d["U1COORD"], d["V1COORD"], d["U2COORD"], d["V2COORD"]
    return (
        uniform_disk_amplitude(u1, v1, wave, diameter_mas)
        * uniform_disk_amplitude(u2, v2, wave, diameter_mas)
        * uniform_disk_amplitude(u1 + u2, v1 + v2, wave, diameter_mas)
    )


@dataclasses.dataclass(frozen=True)
class Observable:
    """Where an observable lives in an :class:`~virgil_vlti.products.Exposure`.

    ``model(table, wave, diameter_mas)`` gives a calibrator's expected value
    for an amplitude; phases have none (a calibrator's phases are taken as
    zero).
    """

    table: str  # Exposure attribute: "v2", "t3" or "vis"
    value: str
    error: str
    ucol: str
    model: object = None

    @property
    def phase(self):
        return self.model is None


OBSERVABLES = {
    "vis2": Observable("v2", "VIS2DATA", "VIS2ERR", "UCOORD", uniform_disk_v2),
    "t3amp": Observable("t3", "T3AMP", "T3AMPERR", "U1COORD", uniform_disk_t3amp),
    "t3phi": Observable("t3", "T3PHI", "T3PHIERR", "U1COORD"),
    "visamp": Observable("vis", "VISAMP", "VISAMPERR", "UCOORD", uniform_disk_visamp),
    "visphi": Observable("vis", "VISPHI", "VISPHIERR", "UCOORD"),
}


def _bad(table, i, value, error, ucol):
    """Flagged, non-finite or bad-error channels, or a whole row with no data.

    A bad error is non-finite, zero or negative (an infinite weight downstream), or huge.

    An empty row has (u, v) = (0, 0): u = 0 alone is a genuine north–south
    baseline. ``ucol`` is ``UCOORD`` or ``U1COORD``, and v is read from the
    matching ``VCOORD`` or ``V1COORD``.
    """
    d = table.data
    err = d[error][i]
    bad = d["FLAG"][i] | ~np.isfinite(d[value][i]) | ~np.isfinite(err) | (err <= 0) | (err > 1e3)
    vcol = "V" + ucol[1:]
    return bad | ((d[ucol][i] == 0) & (d[vcol][i] == 0))


@dataclasses.dataclass
class Estimates:
    """Calibrator estimates for one baseline or triangle, in time order.

    ``values`` are transfer-function amplitudes or unit phasors, shape
    (n_calibrators, n_channels), each the mean over one file's frames;
    ``errors`` are the errors of those means.
    """

    times: np.ndarray
    values: np.ndarray
    errors: np.ndarray


def fit_channel_polynomial(values, errors, flags, maxdeg):
    """``gravi_vis_fit_amp`` (GRAVITY pipeline 1.9.7, ``gravi_vis.c``): a polynomial in the
    channel index, each residual weighted by σ⁻² (so σ⁻⁴ on its square), flagged
    channels and errors below 1e-5 excluded, evaluated at every channel.
    With no more usable channels than the degree, ``values`` is returned.
    """
    nwave = values.size
    w = np.where(flags, 0.0, errors**-2.0)
    w = np.where(w > 1e10, 0.0, w)  # "faulty uncertainty" in the source
    if np.count_nonzero(w) <= maxdeg:
        return values
    x = np.arange(nwave, dtype=float)
    x = 2 * x / max(nwave - 1, 1) - 1  # rescaled for conditioning; the fit is the same
    A = np.vander(x, maxdeg + 1, increasing=True) * w[:, None]
    coef, *_ = np.linalg.lstsq(A, np.where(flags, 0.0, values) * w, rcond=None)
    return np.vander(x, maxdeg + 1, increasing=True) @ coef


def smooth_spectrum(wave, values, errors, degree):
    """Fit a polynomial across channels, as ``gravity_viscal`` smooths the transfer function.

    :func:`fit_channel_polynomial` (``smoothing``, ``maxdeg-tfvis-sc``;
    manual §6.13), excluding invalid channels and evaluated where ``values``
    is finite. ``wave`` only sets the number of channels.
    """
    if int(degree) != degree or degree < 0:
        raise ValueError(f"smoothing degree must be a non-negative integer, not {degree!r}")
    values, errors = np.asarray(values, float), np.asarray(errors, float)
    bad = ~np.isfinite(values) | ~np.isfinite(errors) | (errors <= 0)
    fitted = fit_channel_polynomial(values, errors, bad, int(degree))
    return np.where(np.isfinite(values), fitted, np.nan)


# Calibrator channels whose model amplitude is at most this (near a null of
# the uniform disk) give no transfer-function estimate.
MIN_CALIBRATOR_MODEL = 1e-6


def calibrator_estimates(cals, insname, observable, diameters=None, smooth_degree=None):
    """Transfer-function estimates per station tuple, one per calibrator file.

    Parameters
    ----------
    cals : list of Exposure
    insname : str
    observable : str
        A key of :data:`OBSERVABLES`.
    diameters : dict, optional
        Uniform-disk diameter (mas) by target name; the key None gives a
        default for every calibrator. Absent means unresolved.
    smooth_degree : int, optional
        Smooth each amplitude estimate across wavelength with a polynomial of
        this degree (:func:`smooth_spectrum`), as ``gravity_viscal`` does with
        degree 5. Phases are not smoothed.
    """
    diameters = diameters or {}
    spec = OBSERVABLES[observable]
    value, error, ucol = spec.value, spec.error, spec.ucol
    rows = {}
    for cal in cals:
        table = getattr(cal, spec.table)[insname]
        d = table.data
        if not spec.phase:
            model = spec.model(table, cal.wave[insname], diameter_of(diameters, cal.target))
        per = {}
        for i, sta in enumerate(map(tuple, d["STA_INDEX"])):
            bad = _bad(table, i, value, error, ucol)
            if spec.phase:
                val, err = np.exp(1j * np.deg2rad(d[value][i])), d[error][i]
            else:  # amplitude and its error, both in transfer-function units
                # A model at (or near) a null would divide by ~0: no estimate.
                bad = bad | ~(np.abs(model[i]) > MIN_CALIBRATOR_MODEL)
                with np.errstate(invalid="ignore", divide="ignore"):
                    val, err = d[value][i] / model[i], d[error][i] / model[i]
            per.setdefault(sta, []).append(
                (d["MJD"][i], np.where(bad, np.nan, val), np.where(bad, np.nan, err))
            )
        # Each baseline gets the mean time of its own frames: pipelines may
        # stamp the baselines of one frame with different times.
        for sta, items in per.items():
            t = float(np.mean([t for t, _, _ in items]))
            mean = np.nanmean([v for _, v, _ in items], 0)
            if spec.phase:
                mean = mean / np.abs(mean)
            # The error of the mean over the file's frames.
            errs = np.array([e for _, _, e in items])
            with np.errstate(invalid="ignore", divide="ignore"):
                err = np.sqrt(np.nansum(errs**2, 0)) / np.sum(np.isfinite(errs), 0)
            if smooth_degree is not None and not spec.phase:
                mean = smooth_spectrum(cal.wave[insname], mean, err, smooth_degree)
            rows.setdefault(sta, []).append((t, mean, err))
    out = {}
    for sta, items in rows.items():
        items.sort(key=lambda r: r[0])
        out[sta] = Estimates(
            times=np.array([r[0] for r in items]),
            values=np.array([r[1] for r in items]),
            errors=np.array([r[2] for r in items]),
        )
    return out


# ---------------------------------------------------------------------------
# Interpolation in time, and its uncertainty
# ---------------------------------------------------------------------------


def _linear_weights(t, times, valid):
    """Linear-interpolation weights (estimates, channels) at ``t``, per channel.

    Each channel interpolates between its nearest *valid* estimates either
    side of ``t`` (the nearest one beyond the ends); a channel with none has
    zero weight everywhere. ``times`` must be sorted.
    """
    n, nchan = valid.shape
    w = np.zeros((n, nchan))
    patterns, which = np.unique(valid.T, axis=0, return_inverse=True)
    for p, pattern in enumerate(patterns):
        idx = np.flatnonzero(pattern)
        if idx.size == 0:
            continue
        cols = np.flatnonzero(which.ravel() == p)
        tv = times[idx]
        if idx.size == 1 or t <= tv[0]:
            w[idx[0], cols] = 1
        elif t >= tv[-1]:
            w[idx[-1], cols] = 1
        else:
            k = np.searchsorted(tv, t) - 1
            f = (t - tv[k]) / (tv[k + 1] - tv[k])
            w[idx[k], cols], w[idx[k + 1], cols] = 1 - f, f
    return w


# Floors on the error that weights a calibrator estimate, adapted from
# gravi_tf.c (1.9.7; see viscal.py): 2% of the transfer function for
# amplitudes, 0.2 degrees for phases. The source floors the per-frame error;
# here the estimate's error is that of the file mean, so the floors bite sooner.
SIGMA_FLOOR_REL = 0.02
SIGMA_FLOOR_PHASE_DEG = 0.2


def _gaussian_weights(t, times, values, errors, timescale_s):
    """Weights (estimates, channels) exp(-2 dT²/Δ²)/σ², zero where ``values`` is not finite.

    σ is each estimate's median error, floored as the pipeline does, so that a
    calibrator with an implausibly small error cannot dominate.
    """
    dt = (t - times) * 86400.0
    with np.errstate(all="ignore"):
        sigma = np.nanmedian(errors, axis=1)
        if np.iscomplexobj(values):
            floor = SIGMA_FLOOR_PHASE_DEG
        else:
            floor = SIGMA_FLOOR_REL * np.abs(np.nanmedian(values, axis=1))
        sigma = np.fmax(sigma, floor)
        w = np.exp(-2 * (dt / timescale_s) ** 2) / sigma**2
    return np.where(np.isfinite(values), w[:, None], 0.0)


def predict(t, est, method="linear", timescale_s=PIPELINE_TIMESCALE_S, keep=None):
    """The transfer function at time ``t`` (MJD), per channel.

    ``keep`` optionally selects which estimates to use (for leave-one-out).
    Phasor estimates are interpolated as complex numbers and renormalised.
    """
    times, values, errors = est.times, est.values, est.errors
    if keep is not None:
        times, values, errors = times[keep], values[keep], errors[keep]
    if method == "linear":
        w = _linear_weights(t, times, np.isfinite(values))
        with np.errstate(invalid="ignore", divide="ignore"):
            out = np.sum(w * np.nan_to_num(values), 0) / np.sum(w, 0)
        return out / np.abs(out) if np.iscomplexobj(values) else out
    if method == "gaussian":
        w = _gaussian_weights(t, times, values, errors, timescale_s)
        out = np.nansum(w * np.nan_to_num(values), 0) / np.sum(w, 0)
        return out / np.abs(out) if np.iscomplexobj(values) else out
    raise ValueError(f"Unknown interpolation {method!r}")


def interpolation_weights(t, est, method="linear", timescale_s=PIPELINE_TIMESCALE_S):
    """The weights (estimates, channels), summing to 1, that :func:`predict` uses at ``t``."""
    times, values, errors = est.times, est.values, est.errors
    valid = np.isfinite(values)
    if method == "linear":
        w = _linear_weights(t, times, valid)
    elif method == "gaussian":
        w = _gaussian_weights(t, times, values, errors, timescale_s)
    else:
        raise ValueError(f"Unknown interpolation {method!r}")
    with np.errstate(invalid="ignore", divide="ignore"):
        return w / np.sum(w, 0)


def propagated_error(t, est, method="linear", timescale_s=PIPELINE_TIMESCALE_S):
    """The interpolated transfer function's statistical error at ``t``: √(Σ w² σ²).

    Only the calibrators' reported errors, carried through the interpolation:
    no calibration systematic. In V² units for V², degrees for phases.
    """
    w = interpolation_weights(t, est, method, timescale_s)
    return np.sqrt(np.nansum(w**2 * est.errors**2, 0))


def _residual(value, prediction):
    if np.iscomplexobj(value):
        return np.angle(value * np.conj(prediction), deg=True)
    return value - prediction


def prediction_error(est, mode="loo", method="linear", timescale_s=PIPELINE_TIMESCALE_S, t=None):
    """The uncertainty of the interpolated transfer function, per channel.

    In the observable's units for amplitudes, in degrees for phases. The
    systematic is the leave-one-out prediction error (``"loo"``) or the
    night's scatter (``"scatter"``); it needs two calibrator estimates and
    is zero with fewer. Given the time ``t`` of the datum, the error is floored
    at the calibrators' own statistical error carried through the
    interpolation (:func:`propagated_error`), which the systematic already
    contains for two or more estimates, so it is not added to it: with one
    calibrator estimate it is all there is, and a
    warning says so, since nothing measures the transfer function's
    variation. Without ``t`` only the systematic is returned. ``"none"``
    is zero.
    """
    n = len(est.times)
    nchan = est.values.shape[1]
    if mode == "none":
        return np.zeros(nchan)
    if mode not in ("loo", "scatter"):
        raise ValueError(f"Unknown error mode {mode!r}")
    if n < 2:
        systematic = np.zeros(nchan)
    elif mode == "scatter":
        if np.iscomplexobj(est.values):
            mean = np.nanmean(est.values, 0)
            systematic = np.nanstd(np.angle(est.values * np.conj(mean), deg=True), 0, ddof=1)
        else:
            systematic = np.nanstd(est.values, 0, ddof=1)
    else:
        res = []
        for c in range(n):
            keep = np.arange(n) != c
            res.append(
                _residual(
                    est.values[c],
                    predict(est.times[c], est, method, timescale_s, keep=keep),
                )
            )
        systematic = np.sqrt(np.nanmean(np.square(res), 0))
    if t is None:
        return systematic
    if n == 1:
        warnings.warn(
            "One calibrator estimate: the transfer function's variation cannot be measured, "
            "so its error is only the calibrator's statistical error and is likely too small",
            stacklevel=3,
        )
    # a floor, not a sum: the leave-one-out residuals and the scatter already contain
    # the calibrators' statistical noise
    return np.maximum(systematic, propagated_error(t, est, method, timescale_s))


def prediction_residuals(est, mode="loo", method="linear", timescale_s=PIPELINE_TIMESCALE_S):
    """The residuals behind :func:`prediction_error`, (estimates, channels).

    ``"loo"``: each estimate minus its prediction from the others;
    ``"scatter"``: each estimate's deviation from the night's mean (phases as
    angles from the mean phasor). In the observable's units, degrees for
    phases. Empty with fewer than two estimates.
    """
    n, nchan = est.values.shape
    if n < 2 or mode not in ("loo", "scatter"):
        return np.zeros((0, nchan))
    if mode == "loo":
        rows = []
        for c in range(n):
            keep = np.arange(n) != c
            rows.append(
                _residual(est.values[c], predict(est.times[c], est, method, timescale_s, keep=keep))
            )
        return np.array(rows)
    mean = np.nanmean(est.values, 0)
    if np.iscomplexobj(est.values):
        mean = mean / np.abs(mean)
    return _residual(est.values, mean)


def channel_correlation(residuals):
    """How much of each channel's transfer-function error is common to all channels.

    Fits the channels' correlation matrix of ``residuals`` (estimates ×
    channels, from :func:`prediction_residuals`) by one shared mode: its
    leading eigenvector, scaled by √λ. Returns r (one value per channel, in
    [−1, 1]): a fraction r² of channel j's error variance τ_j² is common to
    all channels (a mode r τ), and 1 − r² is independent between channels.
    The signs are relative: anticorrelated channels get opposite signs, so
    the mode reproduces the fitted covariance; only the global sign is
    arbitrary, and it is chosen to make Σ r ≥ 0.
    Only estimates valid in every channel are used; with fewer than two, or a
    single channel, everything is common (r = 1). With n estimates the
    correlation matrix has rank ≤ n, so few calibrators overstate the common
    part: the safe side, since a common error cannot be averaged down.
    """
    residuals = np.asarray(residuals, float)
    nchan = residuals.shape[1]
    rows = residuals[np.all(np.isfinite(residuals), 1)]
    if len(rows) < 2 or nchan < 2:
        return np.ones(nchan)
    rms = np.sqrt(np.mean(rows**2, 0))
    live = rms > 0
    r = np.ones(nchan)
    if live.sum() < 2:
        return r
    z = rows[:, live] / rms[live]
    corr = z.T @ z / len(z)
    lam, vec = np.linalg.eigh(corr)
    v = vec[:, -1] if vec[:, -1].sum() >= 0 else -vec[:, -1]
    r[live] = np.clip(np.sqrt(max(lam[-1], 0.0)) * v, -1.0, 1.0)
    return r


TF_MODE_SHAPES = ("flat", "poly1", "poly2", "raw")


def regularize_shape(common, total, wave, shape="poly1"):
    """Smooth a transfer-function mode's channel shape, keeping each channel's variance.

    ``common`` is the mode's 1σ size per channel (from :func:`channel_correlation`,
    r τ) and ``total`` the channel's whole transfer-function error τ, both in
    the units the mode is stored in. Both are estimated channel by channel
    from a handful of calibrators, so the shape r τ scatters from channel to
    channel (a coefficient of variation of about 0.16 for V² and 0.25 for
    closure phases with seven calibrators). A mode spanning hundreds of data
    has its amplitude pinned by them, so marginalizing it removes only the
    transfer-function error along that noisy shape; the rest, the shape's
    scatter times the realized error, is left unmodelled and the fitted
    parameters' errors come out about twofold too small. Replacing the shape
    by a smooth one fixes this:

    * ``"flat"``: the mean of ``|common|`` over its finite channels;
    * ``"poly1"``, ``"poly2"``: a least-squares polynomial of that degree in
      wavelength (for a transfer-function error with a colour);
    * ``"raw"``: ``common`` unchanged (for comparison only).

    The smoothing acts on the size ``|common|`` and each channel keeps its
    sign: :func:`channel_correlation` gives anticorrelated channels opposite
    signs (a closure-phase offset that changes sign across the band, a V²
    tilt about its centre), and a signed mean would shrink such a mode
    towards zero and make it one-signed, sending its variance to the
    independent remainder, which averages down over channels and exposures.

    Returns ``(shape, remainder, added)``: the smoothed shape, the
    independent per-channel remainder √(τ² − shape²), so that shape² +
    remainder² = τ² in every channel, and the variance added where shape² >
    τ² (zero elsewhere). There the remainder is clipped to zero and the
    channel's total error grows to |shape|, the conservative side. Channels
    where ``common`` or ``total`` is not finite keep NaN.
    """
    common = np.asarray(common, float)
    total = np.asarray(total, float)
    if shape not in TF_MODE_SHAPES:
        raise ValueError(f"shape must be one of {TF_MODE_SHAPES}, not {shape!r}")
    ok = np.isfinite(common) & np.isfinite(total)
    fitted = np.full(common.shape, np.nan)
    sign = np.where(np.signbit(common), -1.0, 1.0)  # the mode's sign in each channel
    if shape == "raw":
        fitted[ok] = common[ok]
    elif ok.any() and shape == "flat":
        fitted[ok] = np.mean(np.abs(common[ok])) * sign[ok]
    elif ok.any():
        x = np.asarray(wave, float)
        ok = ok & np.isfinite(x)
        x = (x - np.mean(x[ok])) / max(np.ptp(x[ok]), np.finfo(float).tiny)
        degree = min(int(shape[4:]), ok.sum() - 1)
        coef = np.polynomial.polynomial.polyfit(x[ok], np.abs(common[ok]), degree)
        fitted[ok] = np.polynomial.polynomial.polyval(x[ok], coef) * sign[ok]
    with np.errstate(invalid="ignore"):
        gap = np.where(ok, total**2 - fitted**2, np.nan)
        remainder = np.where(ok, np.sqrt(np.clip(gap, 0, None)), total)
        added = np.where(ok, np.clip(-gap, 0, None), 0.0)
        added = np.where(added > 1e-9 * total**2, added, 0.0)  # not rounding
    return fitted, remainder, added


# |slope / its standard error| above which a flat mode is questioned: with ~10 binned
# channels a pure-noise slope exceeds 5 about once in a thousand times per baseline
SLOPE_WARN = 5.0


def shape_slope_significance(common, wave):
    """The significance of a linear colour in a mode's channel shape.

    Least-squares fits ``|common|`` (the mode's size per channel) as a line in
    wavelength and returns slope / (its standard error from the scatter about
    the line), or 0 with fewer than five finite channels or no scatter. A
    large value says the calibrators' common error has a colour that the
    ``"flat"`` shape cannot follow (use ``"poly1"``).
    """
    y = np.abs(np.asarray(common, float))
    x = np.asarray(wave, float)
    ok = np.isfinite(y) & np.isfinite(x)
    if ok.sum() < 5 or np.ptp(x[ok]) == 0:
        return 0.0
    x = (x[ok] - np.mean(x[ok])) / np.ptp(x[ok])
    y = y[ok]
    sxx = x @ x
    slope = (x @ (y - np.mean(y))) / sxx
    resid = y - np.mean(y) - slope * x
    se = np.sqrt(resid @ resid / (len(y) - 2) / sxx)
    return float(abs(slope) / se) if se > 1e-12 * np.mean(y) else 0.0


# ---------------------------------------------------------------------------
# Calibration of one science exposure
# ---------------------------------------------------------------------------


def bin_channels(x, keep, n):
    """Group the kept channels in groups of ``n`` along a new last axis."""
    x = x[..., keep]
    m = x.shape[-1] // n * n
    return x[..., :m].reshape(*x.shape[:-1], -1, n)


def channel_window(wave, window_um=None):
    """Boolean mask of the channels inside ``window_um`` (all if None)."""
    if window_um is None:
        return np.ones(wave.shape, bool)
    lo, hi = window_um
    return (wave > lo * 1e-6) & (wave < hi * 1e-6)


@dataclasses.dataclass
class Observed:
    """One calibrated, binned observable: values, errors and flags."""

    value: np.ndarray
    error: np.ndarray
    flag: np.ndarray
    tf_mode: np.ndarray = None  # see calibrate_exposure(tf_modes=True)


@dataclasses.dataclass
class Calibrated:
    """Calibrated, binned observables of one instrument in one exposure.

    ``obs`` maps names in :data:`OBSERVABLES` to :class:`Observed`.
    """

    insname: str
    wave: np.ndarray
    band: np.ndarray
    obs: dict
    n_cal: int
    cp_only: bool = False  # closure phases only: see closure_phase_only

    def __getitem__(self, name):
        return self.obs[name]


def binned_estimates(est, keep, nbin):
    """Calibrator estimates on the windowed, binned channel grid.

    Each bin's value is the mean of its valid channels (phasors renormalised)
    and its error that mean's error. The transfer function's uncertainty is
    computed from these, not channel by channel: its errors are common to
    neighbouring channels (on real GRAVITY data the per-channel V² scatter is
    almost fully correlated within a bin of 4), so averaging them in
    quadrature would shrink them by up to √nbin.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # empty bins stay NaN
        values = np.nanmean(bin_channels(est.values, keep, nbin), -1)
        if np.iscomplexobj(values):
            values = values / np.abs(values)
        b = bin_channels(est.errors, keep, nbin)
        n = np.sum(np.isfinite(b), -1)
        errors = np.where(n > 0, np.sqrt(np.nansum(b**2, -1)) / np.maximum(n, 1), np.nan)
    return Estimates(times=est.times, values=values, errors=errors)


def _apply(
    table,
    spec,
    estimates,
    keep,
    nbin,
    method,
    error,
    timescale_s,
    tf_modes=False,
    tf_mode_shape="poly1",
    wave=None,
):
    """Calibrate one table, then window and bin it.

    The errors have two parts, combined in quadrature after binning: the
    science data's own errors, independent between channels (binned in
    quadrature), and the transfer function's uncertainty, computed on the
    binned calibrator estimates (:func:`binned_estimates`), so that what is
    common to a bin's channels is not averaged down.

    With ``tf_modes``, the transfer function's uncertainty is split
    (:func:`channel_correlation`) into a part common to all channels, returned
    as ``Observed.tf_mode`` (relative for amplitudes, degrees for phases), and
    an independent remainder, the only part added to the errors. The common
    part's channel shape is smoothed by :func:`regularize_shape`
    (``tf_mode_shape``, over the binned wavelengths ``wave``) in the units it
    is stored in, and the variance this takes from it is returned to the
    remainder, so each channel's total transfer-function error is unchanged.
    """
    d = table.data
    out, errs, tf_rel, shapes = [], [], [], []
    missing = set()
    split = {}
    excess = []
    colour = {}
    nb = bin_channels(np.zeros(d[spec.value].shape[-1]), keep, nbin).shape[0]
    for i in range(len(d)):
        sta = tuple(d["STA_INDEX"][i])
        nw = d[spec.value][i].size
        if sta not in estimates:
            missing.add(sta)
            out.append(np.full(nw, np.nan))
            errs.append(np.full(nw, np.nan))
            tf_rel.append(np.full(nb, np.nan))
            shapes.append(np.zeros(nb))
            continue
        est = estimates[sta]
        g = predict(d["MJD"][i], est, method, timescale_s)
        coarse = binned_estimates(est, keep, nbin)
        if error == "propagated":
            tau = propagated_error(d["MJD"][i], coarse, method, timescale_s)
        else:
            tau = prediction_error(coarse, error, method, timescale_s, d["MJD"][i])
        common = np.zeros_like(tau)
        if tf_modes:  # the same for every row of this baseline or triangle (loo or scatter)
            systematic = prediction_error(coarse, error, method, timescale_s)
            if sta not in split:
                split[sta] = channel_correlation(
                    prediction_residuals(coarse, error, method, timescale_s)
                )
            r = split[sta]
            # The mode is the systematic's common part. The calibrators' statistical error
            # (the floor in tau) is independent between channels: it stays on the diagonal.
            common = systematic * r
            if tf_mode_shape == "raw":  # the channel-by-channel estimate, unsmoothed
                tau = np.sqrt(np.maximum(tau**2 - common**2, 0.0))
            # else: split below, in the stored units
        bad = _bad(table, i, spec.value, spec.error, spec.ucol)
        if spec.phase:
            bad = bad | ~np.isfinite(g)
            cal = np.exp(1j * np.deg2rad(d[spec.value][i])) * np.conj(g)
            e = d[spec.error][i]
            total = tau  # degrees
        else:
            bad = bad | ~(g > 0)
            cal = d[spec.value][i] / g
            e = d[spec.error][i] / g
            with np.errstate(invalid="ignore", divide="ignore"):
                level = predict(d["MJD"][i], coarse, method, timescale_s)
                total, common = tau / level, common / level  # relative
        if tf_modes and tf_mode_shape != "raw":
            if sta not in colour:
                colour[sta] = shape_slope_significance(common, wave)
            variance = total**2
            common, total, added = regularize_shape(common, total, wave, tf_mode_shape)
            excess.append((added > 0, added, np.where(np.isfinite(variance), variance, 0.0)))
        tf_rel.append(total)
        shapes.append(common)
        out.append(np.where(bad, np.nan, cal))
        errs.append(np.where(bad, np.nan, e))
    if missing:
        warnings.warn(f"No calibrator data on stations {sorted(missing)}; flagged", stacklevel=3)
    steep = {sta: t for sta, t in colour.items() if t > SLOPE_WARN}
    if tf_mode_shape == "flat" and steep:
        warnings.warn(
            f"{spec.value}: the calibrators' common transfer-function error has a significant "
            f"colour on {len(steep)} of {len(colour)} baselines or triangles (linear slope up "
            f"to {max(steep.values()):.1f} standard errors); the flat mode shape cannot follow "
            "it and its errors may be too small. Use tf_mode_shape='poly1'",
            stacklevel=3,
        )
    over = np.array([x[0] for x in excess])
    if over.any():
        added = sum(np.sum(x[1]) for x in excess) / sum(np.sum(x[2]) for x in excess)
        warnings.warn(
            f"{spec.value}: the {tf_mode_shape} transfer-function mode is larger than the "
            f"channel's transfer-function error in {over.sum()} of {over.size} binned channels "
            "(all rows); their independent remainder is clipped to zero, so their total error "
            f"grows to the mode's size, adding {100 * added:.1f}% to the transfer-function "
            "variance overall",
            stacklevel=3,
        )
    b = bin_channels(np.array(out), keep, nbin)
    be = bin_channels(np.array(errs), keep, nbin)
    n = np.sum(np.isfinite(b), -1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # empty bins stay flagged
        mean = np.nanmean(b, -1)
        stat = np.sqrt(np.nansum(be**2, -1)) / np.maximum(n, 1)
        if spec.phase:
            tf = np.array(tf_rel)
            mean = np.angle(mean, deg=True)
        else:
            tf = np.abs(mean) * np.array(tf_rel)
        err = np.hypot(stat, np.nan_to_num(tf))
    flag = n < (nbin + 1) // 2
    mode = np.where(flag, 0.0, np.nan_to_num(np.array(shapes))) if tf_modes else None
    return Observed(mean, err, flag, mode)


def _vis_definition(exposure, insname):
    """The ``OI_VIS`` table of one exposure, checked, and its phase definition."""
    table = exposure.vis.get(insname)
    if table is None:
        raise ValueError(f"{exposure.path} has no OI_VIS table for {insname}")
    amptyp = str(table.header.get("AMPTYP", "absolute")).strip().lower()
    if amptyp != "absolute":
        raise ValueError(
            f"{exposure.path}: OI_VIS has AMPTYP = {amptyp!r}; only visibility amplitudes "
            "can be calibrated by a transfer-function ratio"
        )
    if "VISREFMAP" in table.columns.names:
        raise ValueError(
            f"{exposure.path}: OI_VIS has a VISREFMAP (the reference channels of its "
            "differential phases), which windowing and binning would invalidate; not supported"
        )
    phityp = str(table.header.get("PHITYP", "absolute")).strip().lower()
    return phityp, table.header.get("PHIORDER")


def observables_of(exposure, insname, vis=False, cals=()):
    """The observables to calibrate: V², T3AMP and T3PHI, plus OI_VIS if asked.

    ``OI_VIS`` is calibrated only when the science exposure and every
    calibrator have one with visibility amplitudes (``AMPTYP`` absolute, the
    default; correlated fluxes would need the calibrator's flux), no
    ``VISREFMAP``, and the same phase definition (``PHITYP``, ``PHIORDER``),
    since a calibrator's phase is subtracted from the target's. The output
    keeps that definition.
    """
    names = ["vis2", "t3amp", "t3phi"]
    if vis:
        definition = _vis_definition(exposure, insname)
        for cal in cals:
            other = _vis_definition(cal, insname)
            if other != definition:
                raise ValueError(
                    f"{cal.path}: OI_VIS phases are (PHITYP, PHIORDER) = {other}, but the "
                    f"science exposure's are {definition}; they cannot calibrate each other"
                )
        names += ["visamp", "visphi"]
    return names


def calibrate_exposure(
    sci,
    cals,
    insname,
    *,
    window_um=None,
    nbin=1,
    method="linear",
    error="loo",
    timescale_s=PIPELINE_TIMESCALE_S,
    diameters=None,
    strict=True,
    allow_unrecorded=False,
    vis=False,
    smooth_degree=None,
    tf_modes=False,
    tf_mode_shape="poly1",
):
    """Calibrate the observables of one science exposure.

    Parameters
    ----------
    sci : Exposure
    cals : list of Exposure
        The calibrator exposures to build the transfer function from.
    insname : str
    window_um : (float, float), optional
        Keep only channels inside this window (µm).
    nbin : int
        Average the kept channels in groups of this many.
    method : {"linear", "gaussian"}
        Interpolation in time; see the module docstring.
    error : {"loo", "scatter", "propagated", "none"}
        The transfer-function uncertainty added to the errors: a leave-one-out
        prediction error, the night's scatter, the calibrators' reported
        errors carried through the interpolation (statistical only, no
        systematic), or none.
    timescale_s : float
        Δ of the ``"gaussian"`` interpolator.
    diameters : dict, optional
        Calibrator uniform-disk diameters (mas) by target name.
    strict : bool
        Refuse calibrators reduced with other pipeline settings
        (:func:`check_compatible`); if False, only warn. Calibrators in
        another instrument mode are always refused.
    allow_unrecorded : bool
        A header key or pipeline parameter missing on either side cannot be
        compared; by default that is refused. If True, only warn.
    vis : bool
        Also calibrate ``OI_VIS`` (VISAMP and VISPHI); see :func:`observables_of`.
    smooth_degree : int, optional
        Smooth the amplitude transfer function across wavelength with a
        polynomial of this degree before interpolating it in time; 5 matches
        ``gravity_viscal``'s default (:func:`smooth_spectrum`).
    tf_modes : bool
        For V² and closure phases, return the transfer function's uncertainty
        (``"loo"`` or ``"scatter"``) as a nuisance mode per baseline or
        triangle instead of adding it to the errors; only its part independent
        between channels (:func:`channel_correlation`) is still added. Each
        :class:`Observed` then has ``tf_mode``, the mode's 1σ size per datum:
        relative (δV²/V²) for V², degrees for closure phases;
        :func:`transfer_function_modes` names and lays them out. Other
        observables keep the diagonal term.
    tf_mode_shape : {"flat", "poly1", "poly2", "raw"}
        How each mode's channel shape is smoothed (:func:`regularize_shape`):
        its channel mean, a polynomial in wavelength, or ``"raw"``,
        the channel-by-channel estimate (for comparison only: its scatter
        leaves part of the error unmodelled). Each channel's total
        transfer-function error is the same for all four.
        Default ``poly1``: ``flat`` under-estimates errors at realistic chromatic
        transfer functions (pull rms 3.1 at colour 0.36, i.e. λ⁻² across K,
        against about 1.3 for ``poly1`` at every colour; ``flat`` is calibrated
        only for achromatic errors, pull 0.98). This is from our own Monte
        Carlo and has not yet been checked against ESO/VLTI documentation of
        transfer-function chromaticity.
    """
    if not cals:
        raise ValueError(f"No calibrator exposures for {sci.path}")
    if tf_mode_shape not in TF_MODE_SHAPES:
        raise ValueError(f"tf_mode_shape must be one of {TF_MODE_SHAPES}, not {tf_mode_shape!r}")
    if tf_modes and error not in ("loo", "scatter"):
        raise ValueError(
            f"tf_modes needs error='loo' or 'scatter' (the systematic to model), not {error!r}"
        )
    check_compatible(sci, cals, insname, strict=strict, allow_unrecorded=allow_unrecorded)
    keep = channel_window(sci.wave[insname], window_um)
    kw = {"keep": keep, "nbin": nbin, "method": method, "error": error, "timescale_s": timescale_s}
    wave = np.mean(bin_channels(sci.wave[insname], keep, nbin), -1)
    obs = {}
    with _quiet_empty_channels():
        for name in observables_of(sci, insname, vis, cals):
            spec = OBSERVABLES[name]
            table = getattr(sci, spec.table)[insname]
            estimates = calibrator_estimates(cals, insname, name, diameters, smooth_degree)
            obs[name] = _apply(
                table,
                spec,
                estimates,
                **kw,
                tf_modes=tf_modes and name in TF_MODE_OBSERVABLES,
                tf_mode_shape=tf_mode_shape,
                wave=wave,
            )
    return Calibrated(
        insname=insname,
        wave=wave,
        band=np.sum(bin_channels(sci.band[insname], keep, nbin), -1),
        obs=obs,
        n_cal=len(cals),
    )


def _unusable(n, nb):
    """An observable of ``n`` rows and ``nb`` channels that is NaN and flagged everywhere."""
    nan = np.full((n, nb), np.nan)
    return Observed(nan, nan.copy(), np.ones((n, nb), bool))


def mark_closure_phase_only(calibrated, sci):
    """Flag the V² and T3AMP of a calibrated exposure and mark it closure-phase only.

    For a file whose calibrators are fine for closure phases but whose V² the
    caller does not trust: the T3PHI stays as calibrated (transfer-function
    error included), and nothing downstream can use V², T3AMP or OI_VIS.
    """
    nb = calibrated.wave.size
    calibrated.obs = {
        "vis2": _unusable(len(sci.v2[calibrated.insname].data), nb),
        "t3amp": _unusable(len(sci.t3[calibrated.insname].data), nb),
        "t3phi": calibrated.obs["t3phi"],
    }
    calibrated.cp_only = True
    return calibrated


def closure_phase_only(sci, insname, *, window_um=None, nbin=1):
    """The science closure phases, windowed and binned but **not calibrated**.

    For a night with no calibrator: V² cannot be used without a transfer
    function, so only closure phases are kept. They go through the same
    window, binning and flagging as :func:`calibrate_exposure`
    (:func:`_apply`, with a unit transfer function and no error from it).

    The result is uncalibrated: a calibrator's *measured* closure phase is
    its instrumental one (non-zero and chromatic), which calibration
    subtracts and this does not, so the T3PHI still holds it. The errors are
    the file's statistical ones only, with nothing for the instrumental
    offset or its drift. Fit it with the template's closure-phase offset
    modes (``--template``, written beside the file) or another nuisance
    offset per triangle, and treat the errors as a lower bound. V² and T3AMP
    are NaN and flagged; the returned :class:`Calibrated` has ``cp_only``
    set. Uncalibrated V² is never offered.
    """
    warnings.warn(
        f"{os.path.basename(str(sci.path))}: closure phases without a calibrator keep the "
        "instrumental closure phase and carry statistical errors only; fit them with the "
        "template's closure-phase offset modes (--template)",
        stacklevel=2,
    )
    keep = channel_window(sci.wave[insname], window_um)
    wave = np.mean(bin_channels(sci.wave[insname], keep, nbin), -1)
    table = sci.t3[insname]
    nw = sci.wave[insname].size
    unit = {  # a calibrator of unit transfer function (phasor 1 + 0j, no error)
        sta: Estimates(
            times=np.array([sci.mjd]), values=np.ones((1, nw), complex), errors=np.zeros((1, nw))
        )
        for sta in {tuple(s) for s in table.data["STA_INDEX"]}
    }
    with _quiet_empty_channels():
        phi = _apply(
            table,
            OBSERVABLES["t3phi"],
            unit,
            keep=keep,
            nbin=nbin,
            method="linear",
            error="none",
            timescale_s=PIPELINE_TIMESCALE_S,
            wave=wave,
        )
    return Calibrated(
        insname=insname,
        wave=wave,
        band=np.sum(bin_channels(sci.band[insname], keep, nbin), -1),
        obs={
            "vis2": _unusable(len(sci.v2[insname].data), wave.size),
            "t3amp": _unusable(len(table.data), wave.size),
            "t3phi": phi,
        },
        n_cal=0,
        cp_only=True,
    )


TF_MODE_OBSERVABLES = ("vis2", "t3phi")


def transfer_function_modes(sci, calibrated):
    """The transfer function's error as nuisance modes, one per night and baseline or triangle.

    The transfer function interpolated to a science exposure comes from the
    night's calibrators, so its error is shared by every science exposure of
    that night on that baseline (or triangle), and across its channels. As
    independent errors per datum it would average down by √N over N
    exposures; as one nuisance parameter per (night, baseline) it does not.

    ``calibrated`` is a :class:`Calibrated` from
    ``calibrate_exposure(..., tf_modes=True)``. Returns ``{mode: {observable:
    array}}`` as :func:`~virgil_vlti.write.write_modes` takes, each array
    shaped like the calibrated table (rows × binned channels), zero outside
    its baseline or triangle and on flagged data:

    * ``tf-vis2-base<a>-<b>@<night>``: on V², d log|V| per 1σ, i.e. half the
      relative V² error (``OIData.with_gains(modes=...)``);
    * ``tf-t3phi-tri<a>-<b>-<c>@<night>``: on closure phases, degrees per 1σ
      (convert to radians for ``OIData.with_closure_offsets(modes=...)``).

    Station numbers are the file's ``STA_INDEX``. Files of one night share
    these names, and a mode with the same name in several files is the same
    parameter. Every exposure of a night is treated as sharing one transfer
    function error per baseline: exact when all are interpolated between the
    same calibrators, an overstatement of the correlation when different
    exposures lean on different calibrators.
    """
    modes = {}
    for name in TF_MODE_OBSERVABLES:
        observed = calibrated.obs.get(name)
        if observed is None or observed.tf_mode is None:
            continue
        spec = OBSERVABLES[name]
        stations = getattr(sci, spec.table)[calibrated.insname].data["STA_INDEX"]
        scale = 1.0 if spec.phase else 0.5  # δV²/V² -> d log|V|
        kind = "tri" if spec.phase else "base"
        for sta in sorted({tuple(int(x) for x in s) for s in stations}):
            rows = np.all(stations == np.array(sta), 1)
            arr = np.zeros_like(observed.tf_mode)
            arr[rows] = scale * observed.tf_mode[rows]
            if np.any(arr):
                label = "-".join(str(x) for x in sta)
                modes[f"tf-{name}-{kind}{label}@{sci.night}"] = {name: arr}
    return modes


def diameter_of(diameters, target):
    """A calibrator's diameter (mas): by name, else the default (key None), else 0.

    Names match as in :func:`~virgil_vlti.products.classify` (case, spaces,
    underscores and a leading ``NAME`` ignored).
    """
    if not diameters:
        return 0.0
    key = normalise_name(target)
    for name, value in diameters.items():
        if name is not None and normalise_name(name) == key:
            return value or 0.0
    return diameters.get(None, 0.0) or 0.0


def diameter_modes(sci, cals, insname, diameters, diameter_errors, relative=False, **kwargs):
    """How each calibrator's diameter error moves the calibrated amplitudes.

    A calibrator's diameter error is not noise: it scales the transfer
    function of every frame it calibrates, in every channel, so it must not
    be added to each datum's error. It is one nuisance per calibrator and
    night (design/gravity_calibration_review.md, §7c.3), entering every
    calibrated amplitude through ∂(datum)/∂θ. This returns those columns,
    scaled to 1σ of the diameter, for virgil's low-rank nuisances (its Stage
    6d) to marginalise.

    The derivative is a secant over [θ − σ, θ + σ], clipped at zero, so it
    is not zero for an unresolved calibrator (θ = 0), whose effect is
    quadratic in θ.

    Parameters
    ----------
    sci, cals, insname
        As for :func:`calibrate_exposure`.
    diameters, diameter_errors : dict
        Diameters and their 1σ errors (mas) by calibrator name, the key None
        giving a default (:func:`diameter_of`).
    relative : bool
        Return model-relative modes instead of absolute ones (see below).
    **kwargs
        Passed to :func:`calibrate_exposure` (window, binning, interpolation).

    Returns
    -------
    dict
        ``{"<calibrator>@<night>": {observable: array}}`` for the amplitude
        observables (V², T3AMP and, if calibrated, VISAMP), each shaped like
        the calibrated values.

    Notes
    -----
    By default the modes are **absolute**: the change of each calibrated
    value (V², |V| or T3AMP) per 1σ of the diameter. They are built from
    the data, so they carry the data's noise, and where the calibrated value
    is near zero they are near zero whatever the model says. As a covariance
    this is the data-built form that biases fits low (Peelle's Pertinent
    Puzzle, Lachaume 2021; review §0 item 11).

    With ``relative=True`` they are **model-relative**, as virgil's
    ``OIData.with_gains(modes=...)`` and the transfer-function and template
    modes are. A calibrated amplitude is the raw one divided by the
    interpolated transfer function, so the diameter changes it by the same
    factor whatever its value: d ln(calibrated) = −d ln(TF). That factor is
    computed from unit data (every unflagged value replaced by 1), so it
    does not depend on the science measurement at all, and returned as
    d log|V| per 1σ for V² (half of d ln V²) and VISAMP, and d ln T3AMP for
    T3AMP. The absolute mode for a model is then 2 V²_model × the V² mode,
    |V|_model × the VISAMP mode, and T3AMP_model × the T3AMP mode. With
    binning, the factor is that of the bin's mean of 1/TF.
    """
    modes = {}
    kwargs = {**kwargs, "strict": False}
    for target in sorted({c.target for c in cals}):
        sigma = diameter_of(diameter_errors, target)
        theta = diameter_of(diameters, target)
        if not (np.isfinite(sigma) and np.isfinite(theta) and sigma >= 0 and theta >= 0):
            raise ValueError(f"{target}: diameter {theta} ± {sigma} mas is not physical")
        if not sigma:
            continue
        lo, hi = max(theta - sigma, 0.0), theta + sigma
        data = _unit_amplitudes(sci, insname) if relative else sci
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # already reported by the main calibration
            plus = calibrate_exposure(
                data, cals, insname, diameters={**(diameters or {}), target: hi}, **kwargs
            )
            minus = calibrate_exposure(
                data, cals, insname, diameters={**(diameters or {}), target: lo}, **kwargs
            )
        scale = sigma / (hi - lo)
        amplitudes = [name for name in plus.obs if not OBSERVABLES[name].phase]
        if relative:
            with np.errstate(invalid="ignore", divide="ignore"):
                modes[f"{target}@{sci.night}"] = {
                    name: np.nan_to_num(
                        (0.5 if name == "vis2" else 1.0)
                        * np.log(plus[name].value / minus[name].value)
                        * scale
                    )
                    for name in amplitudes
                }
        else:
            modes[f"{target}@{sci.night}"] = {
                name: np.nan_to_num((plus[name].value - minus[name].value) * scale)
                for name in amplitudes
            }
    return modes


def _unit_amplitudes(exposure, insname):
    """A copy of ``exposure`` whose amplitudes (V², T3AMP, VISAMP) of ``insname`` are 1.

    Flags and non-finite values are kept, so the same data are calibrated.
    """
    tables = {}
    for spec in OBSERVABLES.values():
        if spec.phase:
            continue
        source = getattr(exposure, spec.table)
        if insname not in source:
            continue
        table = tables.setdefault(spec.table, {**source})
        hdu = table[insname] = table[insname].copy()
        values = hdu.data[spec.value]
        hdu.data[spec.value] = np.where(np.isfinite(values), 1.0, values)
    return dataclasses.replace(exposure, **tables)


def check_compatible(sci, cals, insname, strict=True, allow_unrecorded=False):
    """Check that the calibrators were observed and reduced like the target.

    A calibrator in another instrument mode (spectral resolution,
    polarisation), without this instrument's tables, or on another
    wavelength grid cannot calibrate the target: that always raises.
    Different pipeline settings that shape V² raise too, or with
    ``strict=False`` only warn. Differences in DIT or pipeline version are
    warned about. A header key or pipeline parameter missing on either side
    (``NOT_RECORDED``) cannot be compared: that raises too, unless
    ``allow_unrecorded``, when it only warns. The rules are the instrument's
    (:data:`virgil_vlti.products.INSTRUMENT_MODULES`): see
    :func:`virgil_vlti.gravity.compatibility` and
    :func:`virgil_vlti.matisse.compatibility`. A calibrator from another
    instrument always raises.
    """
    name = os.path.basename(sci.path)
    fatal, biased, unrecorded, notes = [], [], [], set()
    rules = instrument_module(sci.header)
    for cal in cals:
        tag = os.path.basename(cal.path)
        if instrument_of(cal.header) != instrument_of(sci.header):
            fatal.append(
                f"{tag}: INSTRUME {instrument_of(cal.header)!r}, "
                f"science {instrument_of(sci.header)!r}"
            )
            continue
        mode, reduction, minor = rules.compatibility(sci.header, cal.header)
        for m in (*mode, *reduction):
            # no check could be made: a named outcome, with its own opt-out
            (unrecorded if NOT_RECORDED in m else fatal if m in mode else biased).append(
                f"{tag}: {m}"
            )
        notes.update(minor)
        if insname not in cal.wave:
            fatal.append(f"{tag}: no {insname} tables")
        elif cal.wave[insname].shape != sci.wave[insname].shape or not np.allclose(
            cal.wave[insname], sci.wave[insname], rtol=1e-6
        ):
            fatal.append(f"{tag}: another wavelength grid")
    if fatal:
        raise ValueError(f"Calibrators incompatible with {name}: " + "; ".join(sorted(fatal)))
    if unrecorded:
        message = f"Compatibility with {name} cannot be checked: " + "; ".join(sorted(unrecorded))
        if not allow_unrecorded:
            raise ValueError(message + " (allow_unrecorded=True or --allow-unrecorded to proceed)")
        warnings.warn(message, stacklevel=3)
    if biased:
        message = f"Calibrators incompatible with {name}: " + "; ".join(sorted(biased))
        if strict:
            raise ValueError(message + " (strict=False or --allow-mismatch to proceed)")
        warnings.warn(message, stacklevel=3)
    for note in sorted(notes):
        warnings.warn(f"{name}: calibrators differ in {note}", stacklevel=3)


@contextlib.contextmanager
def _quiet_empty_channels():
    """Silence numpy's warnings about all-NaN channels, which stay flagged."""
    with warnings.catch_warnings():
        for message in (
            "Mean of empty slice",
            "All-NaN slice encountered",
            "Degrees of freedom <= 0",
            "invalid value encountered in divide",
        ):
            warnings.filterwarnings("ignore", message=message, category=RuntimeWarning)
        yield


def frame_times(v2, t3):
    """One MJD per frame, shared by the V² and closure-phase rows.

    Some pipelines (GRAVITY's among them) stamp each baseline of a frame with
    its own time, up to ~12 s apart, beyond the tolerance virgil uses to match
    triangles to baselines. When the tables are frame-ordered (the station
    pattern repeats every frame, in both tables, with the same number of
    frames), each frame gets the mean of its baseline times. Otherwise the
    times are returned unchanged.
    """
    sv, st = v2.data["STA_INDEX"], t3.data["STA_INDEX"]
    nb = len({tuple(s) for s in sv})
    nt = len({tuple(s) for s in st})
    nf = len(sv) // nb
    ordered = (
        len(sv) == nf * nb
        and len(st) == nf * nt
        and np.all(sv == np.tile(sv[:nb], (nf, 1)))
        and np.all(st == np.tile(st[:nt], (nf, 1)))
    )
    if not ordered:
        warnings.warn("Tables are not frame-ordered; times left as stored", stacklevel=2)
        return v2.data["MJD"], t3.data["MJD"]
    mjd = v2.data["MJD"].reshape(nf, nb).mean(1)
    return np.repeat(mjd, nb), np.repeat(mjd, nt)


# ---------------------------------------------------------------------------
# Whole datasets
# ---------------------------------------------------------------------------


def group_by_night(exposures):
    """``{night: exposures}``, in time order."""
    groups = {}
    for e in sorted(exposures, key=lambda e: e.date):
        groups.setdefault(e.night, []).append(e)
    return groups


def flux_ratio(exposures, insname, *, window_um=None, nbin=1, norm_um=None, calibrator_slopes=None):
    """The science-to-calibrator total-flux ratio (``OI_FLUX``), night by night.

    Fluxes are averaged over telescopes, frames and exposures. Dividing by the
    calibrator removes the tellurics and the instrument's response, but the
    calibrator's own spectral lines remain. Fibre injection changes the level
    from night to night, so each night is normalised to its median over
    ``norm_um`` (default: the whole window) and only the shape is meaningful.

    ``calibrator_slopes`` restores the calibrator's own spectral slope
    (:func:`virgil_vlti.gravity.calibrator_sed_slope`): the ratio is multiplied
    by (λ/λ_ref)^slope, since sci/cal = F_sci / F_cal and F_cal ∝ λ^slope. It is
    a number (all nights) or a dict night -> slope; a dict must have every
    night that has both SCI and CAL exposures (``KeyError`` otherwise, so a
    typo cannot leave one night tilted against the rest). With slopes set, a
    night with more than one calibrator target is refused (``ValueError``):
    one slope per night would fit none of them. The slope's own uncertainty
    (``calibrator_sed_slope``'s ``err``) is a fully correlated tilt across
    all nights and is **not** in ``std``, which is only the night-to-night
    scatter; add it as a nuisance tilt mode if it matters. λ_ref is the middle of the
    normalisation window, which only fixes the level that is divided out anyway. The default, None,
    changes nothing.

    Returns a dict of columns: ``wavelength_um``, one per night, ``mean`` and
    ``std``.
    """
    wave = exposures[0].wave[insname]
    keep = channel_window(wave, window_um)
    norm = channel_window(wave, norm_um if norm_um is not None else window_um)
    columns = {"wavelength_um": np.mean(bin_channels(wave, keep, nbin), -1) * 1e6}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for night, group in group_by_night(exposures).items():
            mean = {}
            for kind in ("SCI", "CAL"):
                fluxes = [np.nanmean(read_flux(e, insname), 0) for e in group if e.kind == kind]
                mean[kind] = np.nanmean(fluxes, 0) if fluxes else None
            if mean["SCI"] is None or mean["CAL"] is None:
                continue
            full = mean["SCI"] / mean["CAL"]
            slope = calibrator_slopes
            if isinstance(calibrator_slopes, dict):
                if night not in calibrator_slopes:
                    raise KeyError(f"calibrator_slopes has no entry for night {night!r}")
                slope = calibrator_slopes[night]
            if slope is not None:
                targets = {e.target for e in group if e.kind == "CAL"}
                if len(targets) > 1:
                    raise ValueError(
                        f"night {night}: calibrator_slopes needs one calibrator per night, "
                        f"got {sorted(targets)}"
                    )
                full = full * (wave / np.median(wave[norm])) ** slope
            full = full / np.nanmedian(full[norm])
            columns[night] = np.nanmean(bin_channels(full, keep, nbin), -1)
        nights = [k for k in columns if k != "wavelength_um"]
        stack = np.array([columns[k] for k in nights])
        columns["mean"] = np.nanmean(stack, 0)
        columns["std"] = (
            np.nanstd(stack, 0, ddof=1) if len(nights) > 1 else np.zeros_like(columns["mean"])
        )
    return columns


def injection_scatter(exposures, insname, *, window_um=None):
    """The spread of per-telescope ``OI_FLUX`` levels: a width for virgil's per-spectrum flux scale.

    The flux coupled into each fibre changes from frame to frame and
    telescope to telescope (fibre injection, tip-tilt, transmission), so
    virgil can fit one grey scale per (frame, telescope) spectrum
    (``OIData.with_flux_scale(per="row")``), which needs a prior
    (design/gravity_calibration_review.md, §1c.1). This measures it: the
    band-averaged flux of every ``OI_FLUX`` row (one telescope in one
    frame) of ``exposures``, which should be one target (or one target's
    calibrator) on one night.

    Which width to use. The overall width mixes fixed throughput
    differences between stations with frame-to-frame injection scatter.
    When the station means differ (by more than their own widths), use each
    station's mean with its own width, one prior per telescope; use the
    overall mean and width only when the station means agree.

    The review recommends a log-normal prior. virgil's is Gaussian, which is
    the small-width approximation of it, adequate for a log-width ≲ 0.2.
    Above that, check first for lost-fibre or outlying rows, then either
    sample the scales under a log-normal prior or treat the Gaussian
    marginal as approximate (a wider, symmetric prior than the data show).

    Rows whose band mean is not positive (a lost fibre, over-subtracted
    background) have no logarithm and are dropped, with a warning. They are
    the worst-injected rows, so the widths are then **censored**: biased
    low. ``n_dropped`` says how many.

    Returns
    -------
    dict
        ``mean`` and ``sd`` of the band-averaged row fluxes, in the file's
        flux units (virgil's ``scale=(mean, sd)`` for a model spectrum
        normalized to a band mean of 1); ``log_sd``, the standard deviation
        of their natural logarithm (the log-normal width); ``n``, the
        number of rows used; ``n_dropped``, the rows dropped; and
        ``per_telescope``, ``{STA_INDEX: {"mean", "sd", "log_sd", "n"}}``
        from each station's rows alone (widths NaN with fewer than two).
    """
    keep = channel_window(exposures[0].wave[insname], window_um)
    rows, stations = [], []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for e in exposures:
            flux, sta = read_flux(e, insname, stations=True)
            rows.append(np.nanmean(flux[:, keep], -1))
            stations.append(sta)
    level, sta = np.concatenate(rows), np.concatenate(stations)
    finite = np.isfinite(level)
    good = finite & (level > 0)
    n_dropped = int(np.sum(finite & ~good))
    if n_dropped:
        warnings.warn(
            f"{n_dropped} OI_FLUX rows of {insname} have a non-positive band mean and are "
            "dropped: the injection widths are censored (biased low)",
            stacklevel=2,
        )
    level, sta = level[good], sta[good]
    if level.size < 2:
        raise ValueError(f"Fewer than two positive OI_FLUX rows for {insname}")
    log = np.log(level)

    def _spread(x):
        many = x.size > 1
        return {
            "mean": float(np.mean(x)),
            "sd": float(np.std(x, ddof=1)) if many else np.nan,
            "log_sd": float(np.std(np.log(x), ddof=1)) if many else np.nan,
            "n": int(x.size),
        }

    per_telescope = {int(k): _spread(level[sta == k]) for k in np.unique(sta)}
    return {
        "mean": float(np.mean(level)),
        "sd": float(np.std(level, ddof=1)),
        "log_sd": float(np.std(log, ddof=1)),
        "n": int(level.size),
        "n_dropped": n_dropped,
        "per_telescope": per_telescope,
    }
