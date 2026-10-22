"""A port of ESO's ``gravity_viscal`` calibration, for parity with ESO's products.

This reproduces what the GRAVITY pipeline (``gravity/1.9.7``) does in
``gravity_viscal``, read from its source (``gravi/gravi_tf.c``:
``gravi_compute_tf``, ``gravi_apply_tf_amp``, ``gravi_apply_tf_phi``,
``gravi_calibrate_vis``; ``gravi/gravi_vis.c``: ``gravi_vis_smooth`` and its
helpers). Where the source and the pipeline manual (1.9.7, §6.13, §9.20)
disagree, it follows the source, and says so. It exists to match ESO's
calibrated products (``*_SCI_VIS_CALIBRATED``) channel by channel; the
default calibration in :mod:`virgil_vlti.calibrate` is a separate, documented
design.

What it does, per instrument (``INSNAME``) and for every calibrator frame
(row) separately:

1. **Transfer function** (``gravi_compute_tf``): VIS2 and its error are
   divided by the uniform-disk V², VISAMP and its error by |V|. T3AMP is not
   corrected for the calibrator's size, and phases are used as they are.
2. **Spectral smoothing** (``gravi_vis_smooth``; not in LOW resolution): the
   frame's flags are extended (NaN in any column, amplitude error ≤ 0, and,
   with ≥ 100 channels, errors above 5 times a 9-channel running median of
   the *first* row's errors — the phase errors for OI_VIS and OI_T3), and
   are shared by a table's amplitude and phase (one FLAG column); with
   ``nsmooth > 0`` a running average over ±``nsmooth`` channels (amplitudes
   linearly, phases as unit vectors, both weighted by median-filtered
   errors); then a polynomial of degree ``maxdeg`` in the
   **channel index**, weighted by σ⁻² on the residual (so σ⁻⁴ on its
   square), fitted to VIS2, VISAMP, T3AMP, VISPHI and T3PHI, the phases in
   degrees as they are, and evaluated at every channel.
3. **Weights** (``gravi_apply_tf_amp``, ``_phi``): each calibrator frame of
   the same baseline or triangle gets W = exp(−2|ΔT|/Δ)/σ², with ΔT from
   the tables' ``TIME`` column (seconds from 0 h UT of each file's own
   date). σ is the frame's median error, floored at 2% of its median
   transfer function for amplitudes and at 0.2° for phases. The manual
   (§9.20, eq. 79) gives exp(−2ΔT²/Δ²); the source uses the exponential in
   |ΔT|. Flags are not used here. Across UT midnight ``TIME`` jumps, so
   calibrators on the other side of midnight are effectively ignored;
   that is reproduced.
4. **Calibration**: amplitude/TF, error/TF (no transfer-function error;
   the source marks it ``FIXME``); phase − arg(Σ W exp(iφ)).
5. **Setups** (``gravi_calib_setupstring``): a calibrator is used only if its
   spectral resolution, both polarisation modes and both DITs match the
   science, unless ``force_calib``.

The calibrator diameter is the one ``gravity_vis`` recorded
(``ESO QC TF CAL_SC DIAM``), unless given. Its error is not propagated
(manual §9.20).
"""

import numpy as np
from scipy.ndimage import median_filter

from .calibrate import Calibrated, Observed, fit_channel_polynomial, uniform_disk_amplitude
from .products import kind_from_category

# What gravity_viscal calibrates, with its table and columns, and the
# calibrator model it divides by (None: none).
_OBSERVABLES = {
    "vis2": ("v2", "VIS2DATA", "VIS2ERR", "v2"),
    "visamp": ("vis", "VISAMP", "VISAMPERR", "v"),
    "t3amp": ("t3", "T3AMP", "T3AMPERR", None),
    "visphi": ("vis", "VISPHI", "VISPHIERR", None),
    "t3phi": ("t3", "T3PHI", "T3PHIERR", None),
}
_PHASES = ("visphi", "t3phi")

# Defaults of gravity_viscal (manual 1.9.7 §6.13).
DEFAULTS = {
    "delta_time": 3600.0,
    "smoothing": True,
    "maxdeg": 5,
    "nsmooth": 0,
    "force_calib": False,
}


def setup_string(header):
    """What ``gravi_calib_setupstring`` compares: resolution, polarisations, DITs.

    DITs are rounded as the pipeline prints them (FT in ms to 0.01, SC in s to
    0.01); a missing card is None on both sides.
    """

    def dit(key, scale):
        value = header.get(key)
        return None if value is None else round(float(value) * scale, 2)

    return (
        header.get("ESO INS SPEC RES"),
        header.get("ESO INS POLA MODE"),
        header.get("ESO FT POLA MODE"),
        dit("ESO DET3 SEQ1 DIT", 1e3),
        dit("ESO DET2 SEQ1 DIT", 1.0),
    )


def _running_median_errors(errors):
    """``cpl_vector_filter_median_create(…, 4)``: a 9-channel running median."""
    return median_filter(errors, size=9, mode="nearest")


def _smooth_running(values, errors, flags, nsamp, phase):
    """``gravi_vis_smooth_amp`` / ``_phi``: running mean over ±nsamp channels."""
    out, err = values.copy(), errors.copy()
    nwave = values.size
    w_all = median_filter(errors, size=2 * nsamp + 1, mode="nearest") ** -2.0
    w_all = np.where(flags, 1e-19, w_all)
    v = np.where(flags, 0.0, np.exp(1j * np.deg2rad(values)) if phase else values)
    for k in range(nwave):
        # The source's loop is ``samp < CPL_MIN(nwave, wave + nsamp)``: channel
        # k + nsamp is left out, so the window is k − nsamp … k + nsamp − 1,
        # off-centre by half a channel. Kept as written, to match ESO's output.
        sl = slice(max(0, k - nsamp), min(nwave, k + nsamp))
        w = w_all[sl]
        total = np.sum(v[sl] * w)
        out[k] = np.angle(total, deg=True) if phase else total / np.sum(w)
        err[k] = np.sum(w) ** -0.5
    return out, err


# Per table, as gravi_vis_smooth treats it: the amplitude and phase columns
# (which share one FLAG column), the model the amplitude is divided by, and
# the error column whose running median extends the flags.
_TABLES = {
    "v2": ("VIS2DATA", "VIS2ERR", None, None, "v2", "VIS2ERR"),
    "vis": ("VISAMP", "VISAMPERR", "VISPHI", "VISPHIERR", "v", "VISPHIERR"),
    "t3": ("T3AMP", "T3AMPERR", "T3PHI", "T3PHIERR", None, "T3PHIERR"),
}
# gravi_vis_flag_nan checks all of these that a table has.
_NAN_COLUMNS = ("VIS2DATA", "VIS2ERR", "VISAMP", "VISAMPERR", "VISPHI", "VISPHIERR",
                "T3PHI", "T3PHIERR", "T3AMP", "T3AMPERR")  # fmt: skip


def _table_tf(cal, insname, table_name, diameter, *, smoothing, maxdeg, nsmooth):
    """One calibrator table's transfer function: ``{column: (values, errors)}``.

    ``gravi_compute_tf`` then ``gravi_vis_smooth``. The amplitude and phase
    of one table share their flags, which come from both: NaN in any column,
    the amplitude error ≤ 0, and (with ≥ 100 channels) the median-flag
    column's errors above 5× a running median of its **first** row. The
    amplitude is fitted first; its "faulty uncertainty" flags (error < 1e-5)
    carry over to the phase fit.
    """
    amp, amp_err, phi, phi_err, model, median_col = _TABLES[table_name]
    table = getattr(cal, table_name)[insname]
    d = table.data
    a = np.array(d[amp], float)
    ae = np.array(d[amp_err], float)
    if model is not None:
        v = uniform_disk_amplitude(d["UCOORD"], d["VCOORD"], cal.wave[insname], diameter)
        m = v**2 if model == "v2" else v
        a, ae = a / m, ae / m
    out = {amp: (a, ae)}
    if phi is not None:
        out[phi] = (np.array(d[phi], float), np.array(d[phi_err], float))
    spec_res = str(cal.header.get("ESO INS SPEC RES", "")).upper()
    if not smoothing or spec_res == "LOW":
        return out
    flags = np.array(d["FLAG"], bool)
    for col in _NAN_COLUMNS:
        if col in d.names:
            flags |= ~np.isfinite(np.asarray(d[col], float))
    flags |= ~(ae > 0)  # gravi_vis_flag_lower on the amplitude error
    if a.shape[1] >= 100:  # gravi_vis_flag_median, against the first row
        errs = ae if median_col == amp_err else out[phi][1]
        flags |= errs > 5.0 * _running_median_errors(errs[0])
    for r in range(a.shape[0]):
        if nsmooth > 0:
            a[r], ae[r] = _smooth_running(a[r], ae[r], flags[r], nsmooth, False)
            if phi is not None:
                p, pe = out[phi]
                p[r], pe[r] = _smooth_running(p[r], pe[r], flags[r], nsmooth, True)
        a[r] = fit_channel_polynomial(a[r], ae[r], flags[r], maxdeg)
        flags[r] |= ~flags[r] & (ae[r] ** -2.0 > 1e10)  # "faulty uncertainty" in the source
        if phi is not None:
            p, pe = out[phi]
            p[r] = fit_channel_polynomial(p[r], pe[r], flags[r], maxdeg)
    return out


def transfer_function(cal, insname, name, diameter, *, smoothing, maxdeg, nsmooth):
    """One calibrator's transfer-function rows for one observable.

    Returns ``(stations, times, values, errors)`` per row, ``times`` in
    seconds from the table's ``TIME`` column.
    """
    table_name, value, _, _ = _OBSERVABLES[name]
    tf = _table_tf(cal, insname, table_name, diameter,
                   smoothing=smoothing, maxdeg=maxdeg, nsmooth=nsmooth)  # fmt: skip
    values, errors = tf[value]
    d = getattr(cal, table_name)[insname].data
    return [tuple(s) for s in d["STA_INDEX"]], np.array(d["TIME"], float), values, errors


def _apply(sci, cals, insname, name, diameters, params):
    table_name, value, error, _ = _OBSERVABLES[name]
    table = getattr(sci, table_name)[insname]
    d = table.data
    phase = name in _PHASES
    rows = []  # (stations, time, values, errors) over every calibrator frame
    for cal in cals:
        diameter = diameters.get(cal.target, 0.0) if diameters else 0.0
        sta, times, vals, errs = transfer_function(
            cal, insname, name, diameter,
            smoothing=params["smoothing"], maxdeg=params["maxdeg"], nsmooth=params["nsmooth"],
        )  # fmt: skip
        rows += list(zip(sta, times, vals, errs, strict=True))
    out = np.array(d[value], float)
    err = np.array(d[error], float)
    for r, (sta, t) in enumerate(zip(map(tuple, d["STA_INDEX"]), d["TIME"], strict=True)):
        mine = [(tc, v, e) for s, tc, v, e in rows if s == sta]
        if not mine:
            out[r] = np.nan
            continue
        total, norm = 0.0, 0.0
        for tc, v, e in mine:
            if phase:
                sigma = max(np.nanmedian(e), 0.2)
                v = np.exp(1j * np.deg2rad(v))
            else:
                sigma = max(np.nanmedian(e), 0.02 * abs(np.nanmedian(v)), 1e-10)
            w = np.exp(-2 * abs(t - tc) / params["delta_time"]) / sigma**2
            total = total + w * v
            norm += w
        if phase:
            out[r] = np.angle(np.exp(1j * np.deg2rad(out[r])) / total, deg=True)
        else:
            tf = total / norm
            out[r], err[r] = out[r] / tf, err[r] / tf
    flag = np.array(d["FLAG"], bool) | ~np.isfinite(out) | ~np.isfinite(err) | (err < 0)
    return Observed(out, err, flag)


def calibrate_like_viscal(sci, cals, insname, diameters=None, **params):
    """Calibrate one science exposure as ``gravity_viscal`` does.

    Parameters
    ----------
    sci : Exposure
    cals : list of Exposure
        Calibrator exposures (``*_CAL_VIS``). Those whose setup differs from
        the science are dropped, as by the pipeline, unless ``force_calib``.
    insname : str
    diameters : dict, optional
        Uniform-disk diameters (mas) by calibrator name; by default each
        calibrator's ``pipeline_diameter``.
    **params
        ``delta_time`` (s), ``smoothing``, ``maxdeg``, ``nsmooth``,
        ``force_calib``; defaults in :data:`DEFAULTS`.

    Returns
    -------
    Calibrated
        At the science exposure's own channels (no window or binning), with
        ``vis2``, ``t3amp``, ``t3phi`` and, when the science has ``OI_VIS``,
        ``visamp`` and ``visphi``.
    """
    params = {**DEFAULTS, **params}
    unknown = set(params) - set(DEFAULTS)
    if unknown:
        raise TypeError(f"Unknown gravity_viscal parameters: {sorted(unknown)}")
    if diameters is None:
        diameters = {c.target: c.pipeline_diameter[0] for c in cals if c.pipeline_diameter}
    setup = setup_string(sci.header)
    used = [c for c in cals if params["force_calib"] or setup_string(c.header) == setup]
    if not used:
        raise ValueError(f"No calibrator with the science setup {setup} (force_calib to override)")
    for c in used:
        if kind_from_category(c.header) == "SCI":
            raise ValueError(f"{c.path} is a science product, not a calibrator")
        if c.wave[insname].shape != sci.wave[insname].shape:
            raise ValueError(f"{c.path}: another wavelength grid")
    names = ["vis2", "t3amp", "t3phi"]
    if insname in sci.vis and all(insname in c.vis for c in used):
        names += ["visamp", "visphi"]
    obs = {n: _apply(sci, used, insname, n, diameters, params) for n in names}
    return Calibrated(
        insname=insname,
        wave=np.array(sci.wave[insname]),
        band=np.array(sci.band[insname]),
        obs=obs,
        n_cal=len(used),
    )
