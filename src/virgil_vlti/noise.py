"""Night-shared calibration terms for virgil ``OIData``: a gain per (night, baseline) and a
closure-phase offset per (night, triangle).

A night's calibrators and transfer function are common to every frame of that night, so the V²
gain on a baseline and the closure phase on a triangle share one error per night.
:func:`with_night_terms` attaches them to an ``OIData`` as virgil's spanning modes (marginalized
analytically by the likelihood), together with the per-frame telescope, baseline and triangle
terms that go in the same ``with_gains`` / ``with_closure_offsets`` call (a second call would
replace the modes). It works from the data alone (MJD and station indices), so it applies to any
``OIData``, including simulated.

The night is the Paranal observing night of :class:`virgil_vlti.products.Exposure`: it starts at
12 UT (the same boundary), so it labels files by the date at the start of the night.
"""

from __future__ import annotations

import numpy as np

from .products import NIGHT_BOUNDARY_UT_HOURS

__all__ = ["night_labels", "with_night_terms"]


def night_labels(data, boundary_ut: float = NIGHT_BOUNDARY_UT_HOURS) -> dict:
    """Night, (night, baseline) cell and frame of every sample.

    ``u_*`` is per sample of ``data.u``, ``v2_*`` per V² (``data.vis``), ``t3_*`` per closure
    phase (``data.phi``; its night and frame are those of its first baseline, its triangle the
    set of its three stations). A night is an integer: ``floor(MJD - boundary_ut / 24)``.
    """
    mjd = float(data.t_ref) + np.asarray(data.dt, float)
    night = np.floor(mjd - boundary_ut / 24.0).astype(int)
    st = np.sort(np.asarray(data.stations, int), axis=1)
    frame = np.asarray(data.frame, int)
    base = np.unique(st, axis=0, return_inverse=True)[1].ravel()
    cell_u = np.unique(np.stack([night, base], 1), axis=0, return_inverse=True)[1].ravel()
    vi = data.vis_index
    vi = np.arange(mjd.size) if vi is None or np.ndim(vi) == 0 else np.asarray(vi)
    i1, i2, i3 = (np.asarray(i, int) for i in (data.i_cps1, data.i_cps2, data.i_cps3))
    tri = np.unique(
        np.sort(np.concatenate([st[i1], st[i2], st[i3]], 1), axis=1), axis=0, return_inverse=True
    )[1].ravel()
    cell_t = np.unique(np.stack([night[i1], tri], 1), axis=0, return_inverse=True)[1].ravel()
    return {
        "u_night": night,
        "u_cell": cell_u,
        "v2_night": night[vi],
        "v2_cell": cell_u[vi],
        "v2_frame": frame[vi],
        "t3_night": night[i1],
        "t3_cell": cell_t,
        "t3_frame": frame[i1],
    }


def night_modes(data, night_gain: float = 0.0, night_t3_deg: float = 0.0, boundary_ut=None):
    """``(g, o)``: one row per (night, baseline) on log|V| of width ``night_gain``, one per
    (night, triangle) on the closure phases of width ``night_t3_deg`` (degrees, returned in
    radians). A zero width gives ``None`` for its group."""
    kw = {} if boundary_ut is None else {"boundary_ut": boundary_ut}
    lab = night_labels(data, **kw)
    g = o = None
    if night_gain > 0:
        g = np.zeros((lab["u_cell"].max() + 1, lab["u_cell"].size))
        g[lab["u_cell"], np.arange(lab["u_cell"].size)] = night_gain
    if night_t3_deg > 0:
        o = np.zeros((lab["t3_cell"].max() + 1, lab["t3_cell"].size))
        o[lab["t3_cell"], np.arange(lab["t3_cell"].size)] = np.radians(night_t3_deg)
    return g, o


def with_night_terms(
    data,
    *,
    night_gain: float = 0.0,
    night_t3_deg: float = 0.0,
    telescope: float = 0.0,
    baseline: float = 0.0,
    triangle_deg: float = 0.0,
    boundary_ut=None,
):
    """``data`` with night-shared and per-frame calibration terms.

    ``night_gain`` (log|V|) per (night, baseline) and ``night_t3_deg`` (degrees) per (night,
    triangle) are spanning modes; ``telescope`` and ``baseline`` (log|V|) and ``triangle_deg``
    (degrees) are virgil's per-frame gains and closure offsets, one per frame. A zero width leaves
    its term out. The night closure-phase modes need a virgil whose closure modes may span frames
    (virgil#708, after 0.4.1); on an older virgil this raises ``ValueError``.

    Raises
    ------
    ValueError
        If a requested group would replace gains or closure offsets already on ``data`` (for
        example the template modes ``to_oidata(..., template=...)`` attaches, or its generic
        gains): virgil's ``with_gains`` / ``with_closure_offsets`` replace, never add, so the
        earlier terms would be dropped silently and the error understated. Apply this function
        to data without those terms, or leave its group out.
    """
    g, o = night_modes(data, night_gain, night_t3_deg, boundary_ut)
    want_gains = g is not None or bool(telescope) or bool(baseline)
    want_offsets = o is not None or bool(triangle_deg)
    if want_gains and getattr(data, "gains", None) is not None:
        raise ValueError(
            "data already has gains (e.g. template or generic gains from to_oidata); "
            "with_night_terms would replace them. Leave the V² groups out or start from data "
            "without them."
        )
    if want_offsets and getattr(data, "phase_offsets", None) is not None:
        raise ValueError(
            "data already has closure offsets; with_night_terms would replace them. Leave the "
            "closure groups out or start from data without them."
        )
    if want_gains:
        data = data.with_gains(telescope=telescope or None, baseline=baseline or None, modes=g)
    if want_offsets:
        data = data.with_closure_offsets(
            triangle=np.radians(triangle_deg) if triangle_deg else None, modes=o
        )
    return data
