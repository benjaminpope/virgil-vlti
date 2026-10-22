"""Simulated GRAVITY products with a known truth, for examples and tests.

:func:`write_product` writes one noise-free product file whose V² and closure
phases are functions you give, so that a calibration can be checked exactly,
or noise and systematics of known size added and then measured. The tutorials
and the test suite use it.

The layout follows an archival ``DUAL_SCI_VIS`` product
(design/gravity_calibration_review.md, §10): baselines in beam pairs 12, 13,
14, 23, 24, 34 with beam 1 on the highest station index; triangles (1,2,3),
(1,2,4), (1,3,4), (2,3,4); a fringe-tracker table beside the science
channel; and each baseline of a frame stamped with its own time.
"""

import numpy as np
from astropy.io import fits

STATIONS = (28, 23, 18, 1)  # beams 1-4
BEAM_PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
TRIANGLES = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
# Telescope positions (m, east and north) giving 30-130 m baselines.
POSITIONS = np.array([[0.0, 0.0], [60.0, 20.0], [-30.0, 90.0], [100.0, -60.0]])


def _table(extname, insname, columns, **header):
    hdu = fits.BinTableHDU.from_columns(columns)
    hdu.header.update(EXTNAME=extname, INSNAME=insname, OI_REVN=2, **header)
    return hdu


def write_product(
    path,
    *,
    mjd,
    vis2,
    t3phi,
    wave=None,
    sobj=(0.0, 0.0),
    target="STAR",
    catg="DUAL_SCI_VIS",
    nframes=3,
    insnames=("GRAVITY_SC",),
    flux=1.0,
    time_jitter_s=8.0,
    vis2err=0.01,
    t3phierr=1.0,
    header=None,
    t3amp=None,
    vis=None,
):
    """Write one product whose V² and closure phases are given functions.

    ``vis2(u, v, wave, mjd)`` and ``t3phi(u1, v1, u2, v2, wave, mjd)`` give the
    observed values (degrees for phases) of each row at its own time; the noise
    is zero, so tests can check the calibration exactly. ``t3amp`` (same
    arguments; default 1) gives T3AMP. ``vis(u, v, wave, mjd)``, if given,
    returns (amplitude, phase in degrees) for an ``OI_VIS`` table in GRAVITY's
    layout, with its extra pipeline columns. ``header`` overrides
    or adds primary-header cards (the defaults are a MEDIUM, combined-
    polarisation, VFACTOR-corrected product).
    """
    wave = np.linspace(2.0e-6, 2.4e-6, 40) if wave is None else wave
    primary = fits.PrimaryHDU()
    primary.header.update(
        {
            "DATE-OBS": _date(mjd),
            "MJD-OBS": mjd,
            "OBJECT": target,
            "HIERARCH ESO PRO CATG": catg,
            "HIERARCH ESO INS SOBJ X": sobj[0],
            "HIERARCH ESO INS SOBJ Y": sobj[1],
            "HIERARCH ESO INS SPEC RES": "MEDIUM",
            "HIERARCH ESO INS POLA MODE": "COMBINED",
            "HIERARCH ESO DET2 SEQ1 DIT": 30.0,
            "HIERARCH ESO PRO REC1 PARAM1 NAME": "vis-correction-sc",
            "HIERARCH ESO PRO REC1 PARAM1 VALUE": "VFACTOR",
            **{f"HIERARCH {k}": v for k, v in (header or {}).items()},
        }
    )
    array = fits.BinTableHDU.from_columns(
        [
            fits.Column("TEL_NAME", "16A", array=[f"AT{k + 1}" for k in range(4)]),
            fits.Column("STA_NAME", "16A", array=[f"S{s}" for s in STATIONS]),
            fits.Column("STA_INDEX", "1I", array=STATIONS),
            fits.Column("DIAMETER", "1E", array=np.full(4, 1.8)),
            fits.Column("STAXYZ", "3D", array=np.c_[POSITIONS, np.zeros(4)]),
        ]
    )
    array.header.update(EXTNAME="OI_ARRAY", ARRNAME="VLTI", OI_REVN=2)
    targets = fits.BinTableHDU.from_columns(
        [fits.Column("TARGET_ID", "1I", array=[1]), fits.Column("TARGET", "16A", array=[target])]
    )
    targets.header.update(EXTNAME="OI_TARGET", OI_REVN=2)
    hdus = [primary, array, targets]
    # A fringe-tracker table that must never be calibrated.
    for ins in ("GRAVITY_FT", *insnames):
        w = wave if ins != "GRAVITY_FT" else np.linspace(2.0e-6, 2.4e-6, 6)
        funcs = {"vis2": vis2, "t3phi": t3phi, "t3amp": t3amp, "vis": vis}
        hdus += _instrument(ins, w, mjd, nframes, funcs, flux, time_jitter_s, vis2err, t3phierr)
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return path


def _instrument(ins, wave, mjd, nframes, funcs, flux, jitter, vis2err, t3phierr):
    nw = wave.size
    vis2, t3phi, vis = funcs["vis2"], funcs["t3phi"], funcs["vis"]
    t3amp = funcs["t3amp"] or (lambda u1, v1, u2, v2, w, t: np.ones_like(w))
    wtab = _table(
        "OI_WAVELENGTH",
        ins,
        [
            fits.Column("EFF_WAVE", "1E", "m", array=wave),
            fits.Column("EFF_BAND", "1E", "m", array=np.full(nw, wave[1] - wave[0])),
        ],
    )
    rows, trows = [], []
    for f in range(nframes):
        t = mjd + f * 300 / 86400  # 5-min frames
        rot = 0.02 * f  # the uv tracks rotate a little between frames
        c, s = np.cos(rot), np.sin(rot)
        pos = POSITIONS @ np.array([[c, -s], [s, c]])
        for k, (a, b) in enumerate(BEAM_PAIRS):
            u, v = pos[b] - pos[a]
            tk = t + k * jitter / 6 / 86400
            vrow = vis(u, v, wave, tk) if vis else (None, None)
            rows.append((tk, (STATIONS[a], STATIONS[b]), u, v, vis2(u, v, wave, tk), *vrow))
        for a, b, c3 in TRIANGLES:
            u1, v1 = pos[b] - pos[a]
            u2, v2 = pos[c3] - pos[b]
            trows.append(
                (
                    t,
                    (STATIONS[a], STATIONS[b], STATIONS[c3]),
                    u1,
                    v1,
                    u2,
                    v2,
                    t3phi(u1, v1, u2, v2, wave, t),
                    t3amp(u1, v1, u2, v2, wave, t),
                )
            )
    n = len(rows)
    v2tab = _table(
        "OI_VIS2",
        ins,
        [
            fits.Column("TARGET_ID", "1I", array=np.ones(n)),
            fits.Column("TIME", "1D", "s", array=np.zeros(n)),
            fits.Column("MJD", "1D", "d", array=[r[0] for r in rows]),
            fits.Column("INT_TIME", "1D", "s", array=np.full(n, 180.0)),
            fits.Column("VIS2DATA", f"{nw}D", array=np.array([r[4] for r in rows])),
            fits.Column("VIS2ERR", f"{nw}D", array=np.full((n, nw), vis2err)),
            fits.Column("UCOORD", "1D", "m", array=[r[2] for r in rows]),
            fits.Column("VCOORD", "1D", "m", array=[r[3] for r in rows]),
            fits.Column("STA_INDEX", "2I", array=np.array([r[1] for r in rows])),
            fits.Column("FLAG", f"{nw}L", array=np.zeros((n, nw), bool)),
        ],
        ARRNAME="VLTI",
    )
    m = len(trows)
    t3tab = _table(
        "OI_T3",
        ins,
        [
            fits.Column("TARGET_ID", "1I", array=np.ones(m)),
            fits.Column("TIME", "1D", "s", array=np.zeros(m)),
            fits.Column("MJD", "1D", "d", array=[r[0] for r in trows]),
            fits.Column("INT_TIME", "1D", "s", array=np.full(m, 180.0)),
            fits.Column("T3AMP", f"{nw}D", array=np.array([r[7] for r in trows])),
            fits.Column("T3AMPERR", f"{nw}D", array=np.full((m, nw), 0.01)),
            fits.Column("T3PHI", f"{nw}D", array=np.array([r[6] for r in trows])),
            fits.Column("T3PHIERR", f"{nw}D", array=np.full((m, nw), t3phierr)),
            fits.Column("U1COORD", "1D", "m", array=[r[2] for r in trows]),
            fits.Column("V1COORD", "1D", "m", array=[r[3] for r in trows]),
            fits.Column("U2COORD", "1D", "m", array=[r[4] for r in trows]),
            fits.Column("V2COORD", "1D", "m", array=[r[5] for r in trows]),
            fits.Column("STA_INDEX", "3I", array=np.array([r[1] for r in trows])),
            fits.Column("FLAG", f"{nw}L", array=np.zeros((m, nw), bool)),
        ],
        ARRNAME="VLTI",
    )
    vtab = []
    if vis:
        common = [
            fits.Column("TARGET_ID", "1I", array=np.ones(n)),
            fits.Column("TIME", "1D", "s", array=np.zeros(n)),
            fits.Column("MJD", "1D", "d", array=[r[0] for r in rows]),
            fits.Column("INT_TIME", "1D", "s", array=np.full(n, 180.0)),
        ]
        amp = np.array([r[5] for r in rows])
        phi = np.array([r[6] for r in rows])
        vtab = [
            _table(
                "OI_VIS",
                ins,
                [
                    *common,
                    # A pipeline extra that the writer must drop.
                    fits.Column("VISDATA", f"{nw}C", array=amp * np.exp(1j * np.deg2rad(phi))),
                    fits.Column("VISAMP", f"{nw}D", array=amp),
                    fits.Column("VISAMPERR", f"{nw}D", array=np.full((n, nw), 0.01)),
                    fits.Column("VISPHI", f"{nw}D", "deg", array=phi),
                    fits.Column("VISPHIERR", f"{nw}D", "deg", array=np.full((n, nw), 1.0)),
                    fits.Column("UCOORD", "1D", "m", array=[r[2] for r in rows]),
                    fits.Column("VCOORD", "1D", "m", array=[r[3] for r in rows]),
                    fits.Column("STA_INDEX", "2I", array=np.array([r[1] for r in rows])),
                    fits.Column("FLAG", f"{nw}L", array=np.zeros((n, nw), bool)),
                ],
                ARRNAME="VLTI",
                AMPTYP="absolute",
                PHITYP="absolute",
            )
        ]
    nf = nframes * 4
    ftab = _table(
        "OI_FLUX",
        ins,
        [
            fits.Column("TARGET_ID", "1I", array=np.ones(nf)),
            fits.Column("MJD", "1D", "d", array=np.full(nf, mjd)),
            fits.Column(
                "FLUXDATA",
                f"{nw}D",
                array=np.full((nf, nw), flux(wave) if callable(flux) else flux) * (wave / 2.2e-6),
            ),
            fits.Column("FLUXERR", f"{nw}D", array=np.full((nf, nw), 0.01)),
            fits.Column("STA_INDEX", "1I", array=np.tile(STATIONS, nframes)),
            fits.Column("FLAG", f"{nw}L", array=np.zeros((nf, nw), bool)),
        ],
    )
    return [wtab, *vtab, v2tab, t3tab, ftab]


def _date(mjd):
    from astropy.time import Time

    return Time(mjd, format="mjd").isot[:19]


def add_systematics(
    path,
    insname="GRAVITY_SC",
    *,
    rng,
    white=0.0,
    telescope=0.0,
    baseline=0.0,
    shape=None,
    cp_white=0.0,
    cp_offset=0.0,
):
    """Add noise and calibration systematics of known size to a product, in place.

    Per frame (rows grouped by time, as :func:`virgil_vlti.pairs.frames_of`
    does), each V² is multiplied by ``exp((t_a + t_b + b_ab) · shape(λ))``.
    Here ``t ~ N(0, telescope)`` is a gain per telescope and
    ``b ~ N(0, baseline)`` a gain per baseline, both on log V² (so 0.01 is a
    1% gain on V²), and ``shape`` (one value per channel; default 1, an
    achromatic gain) their spectral shape. Independent noise of
    ``white × VIS2ERR`` is then added.

    Closure phases get baseline-phase noise, ``cp_white × T3PHIERR`` per
    closure phase with the ±1/3 correlation between triangles that share a
    baseline. They also get an offset ``N(0, cp_offset)`` (degrees) per
    triangle, shared by its channels.

    Returns the path.
    """
    from .pairs import closure_matrix, frames_of
    from .products import read_exposure

    e = read_exposure(str(path))
    v2, t3 = e.v2[insname].data, e.t3[insname].data
    nw = v2["VIS2DATA"].shape[1]
    shape = np.ones(nw) if shape is None else np.asarray(shape, float)
    for f in frames_of(e, insname):
        stations = sorted({s for sta in f.v2_rows for s in sta})
        gain = dict(zip(stations, rng.normal(0, telescope, len(stations)), strict=True))
        for sta, row in f.v2_rows.items():
            log_gain = gain[sta[0]] + gain[sta[1]] + rng.normal(0, baseline)
            noise = rng.normal(0, 1, nw) * white * v2["VIS2ERR"][row]
            v2["VIS2DATA"][row] = v2["VIS2DATA"][row] * np.exp(log_gain * shape) + noise
        triangles = list(f.t3_rows)
        if triangles:
            T, baselines = closure_matrix(triangles)
            err = np.mean([t3["T3PHIERR"][f.t3_rows[t]] for t in triangles], axis=0)
            phi = rng.normal(0, 1, (len(baselines), nw)) * cp_white * err / np.sqrt(3)
            noise = T @ phi
            for r, tri in enumerate(triangles):
                t3["T3PHI"][f.t3_rows[tri]] += noise[r] + rng.normal(0, cp_offset)
    with fits.open(path, mode="update") as h:
        for hdu in h:
            if hdu.header.get("INSNAME") == insname and hdu.name == "OI_VIS2":
                hdu.data["VIS2DATA"][:] = v2["VIS2DATA"]
            if hdu.header.get("INSNAME") == insname and hdu.name == "OI_T3":
                hdu.data["T3PHI"][:] = t3["T3PHI"]
    return path
