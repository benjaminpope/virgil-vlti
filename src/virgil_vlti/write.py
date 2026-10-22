"""Write calibrated exposures as OIFITS files that virgil can read."""

import json
import re
from importlib.metadata import PackageNotFoundError, version

import numpy as np
from astropy.io import fits

from .calibrate import frame_times
from .products import OUTPUT_FLAG, product_category


def _array_column(table, name, values, fmt):
    col = table.columns[name]
    nw = np.asarray(values).shape[1]
    return fits.Column(name=name, format=f"{nw}{fmt}", unit=col.unit, array=values)


def _rebuild(table, replace, keep):
    """A copy of ``table`` with the columns in ``replace`` swapped in.

    Only the OIFITS 2 columns in ``keep`` are written: pipeline extras are
    per-channel arrays of the uncalibrated data, whose width the binning
    would leave wrong.
    """
    cols = [replace.get(col.name, col) for col in table.columns if col.name in keep]
    return fits.BinTableHDU.from_columns(cols, header=table.header)


# The OIFITS 2 columns of OI_VIS; pipeline extras (GRAVITY's VISDATA, RVIS,
# IVIS, GDELAY and so on) describe the uncalibrated data and are dropped.
_OI_VIS_COLUMNS = (
    "TARGET_ID", "TIME", "MJD", "INT_TIME", "VISAMP", "VISAMPERR", "VISPHI", "VISPHIERR",
    "UCOORD", "VCOORD", "STA_INDEX", "FLAG",
)  # fmt: skip


_OI_VIS2_COLUMNS = (
    "TARGET_ID", "TIME", "MJD", "INT_TIME", "VIS2DATA", "VIS2ERR",
    "UCOORD", "VCOORD", "STA_INDEX", "FLAG",
)  # fmt: skip
_OI_T3_COLUMNS = (
    "TARGET_ID", "TIME", "MJD", "INT_TIME", "T3AMP", "T3AMPERR", "T3PHI", "T3PHIERR",
    "U1COORD", "V1COORD", "U2COORD", "V2COORD", "STA_INDEX", "FLAG",
)  # fmt: skip


def _finish(hdu, nwave, extver=1):
    hdu.header["NWAVE"] = nwave
    hdu.header["EXTVER"] = extver
    for key in ("CHECKSUM", "DATASUM"):
        hdu.header.remove(key, ignore_missing=True)
    return hdu


def _nan_where(x, flag):
    return np.where(flag, np.nan, x)


def calibrated_hdus(sci, cal, extver=1):
    """``OI_WAVELENGTH``, ``OI_VIS2``, ``OI_T3`` and, if calibrated, ``OI_VIS``.

    ``OI_T3`` and ``OI_VIS`` each have one ``FLAG`` column for an amplitude
    and a phase. It follows the phase, so that a bad amplitude never costs a
    good phase; amplitudes flagged on their own are written as NaN.
    ``extver`` numbers this instrument's tables (the same ``EXTVER`` on its
    wavelength and data tables, unique among the file's tables of one name).
    """
    v2, t3 = sci.v2[cal.insname], sci.t3[cal.insname]
    mjd_v2, mjd_t3 = frame_times(v2, t3)
    nw = cal.wave.size
    w = fits.BinTableHDU.from_columns(
        [
            fits.Column("EFF_WAVE", "1E", "m", array=cal.wave),
            fits.Column("EFF_BAND", "1E", "m", array=cal.band),
        ]
    )
    w.header.update(EXTNAME="OI_WAVELENGTH", OI_REVN=2, INSNAME=cal.insname, EXTVER=extver)
    vis2, amp, phi = cal["vis2"], cal["t3amp"], cal["t3phi"]
    v2 = _rebuild(
        v2,
        {
            "VIS2DATA": _array_column(v2, "VIS2DATA", vis2.value, "D"),
            "VIS2ERR": _array_column(v2, "VIS2ERR", vis2.error, "D"),
            "FLAG": _array_column(v2, "FLAG", vis2.flag, "L"),
            "MJD": fits.Column("MJD", "1D", "d", array=mjd_v2),
        },
        _OI_VIS2_COLUMNS,
    )
    t3 = _rebuild(
        t3,
        {
            "T3AMP": _array_column(t3, "T3AMP", _nan_where(amp.value, amp.flag), "D"),
            "T3AMPERR": _array_column(t3, "T3AMPERR", _nan_where(amp.error, amp.flag), "D"),
            "T3PHI": _array_column(t3, "T3PHI", phi.value, "D"),
            "T3PHIERR": _array_column(t3, "T3PHIERR", phi.error, "D"),
            "FLAG": _array_column(t3, "FLAG", phi.flag, "L"),
            "MJD": fits.Column("MJD", "1D", "d", array=mjd_t3),
        },
        _OI_T3_COLUMNS,
    )
    hdus = [w, _finish(v2, nw, extver), _finish(t3, nw, extver)]
    if "visamp" in cal.obs:
        vis = sci.vis[cal.insname]
        mjd_vis, _ = frame_times(vis, sci.t3[cal.insname])
        a, p = cal["visamp"], cal["visphi"]
        replace = {
            "VISAMP": _array_column(vis, "VISAMP", _nan_where(a.value, a.flag), "D"),
            "VISAMPERR": _array_column(vis, "VISAMPERR", _nan_where(a.error, a.flag), "D"),
            "VISPHI": _array_column(vis, "VISPHI", p.value, "D"),
            "VISPHIERR": _array_column(vis, "VISPHIERR", p.error, "D"),
            "FLAG": _array_column(vis, "FLAG", p.flag, "L"),
            "MJD": fits.Column("MJD", "1D", "d", array=mjd_vis),
        }
        cols = [replace.get(c.name, c) for c in vis.columns if c.name in _OI_VIS_COLUMNS]
        hdus.append(_finish(fits.BinTableHDU.from_columns(cols, header=vis.header), nw, extver))
    return hdus


def flux_hdus(sci, calibrated):
    """The ``OI_FLUX`` tables of ``sci`` for the calibrated instruments, as stored on disk.

    Flux is not calibrated (it is the observed spectrum, as ESO writes it), so
    the tables are copied unchanged; their channels must be those of the
    calibrated tables (no windowing or binning).
    """
    wanted = {cal.insname: cal.wave.size for cal in calibrated}
    extvers = {cal.insname: k for k, cal in enumerate(calibrated, 1)}
    out = []
    with fits.open(sci.path, memmap=False) as h:
        for hdu in h:
            name = str(hdu.header.get("INSNAME", "")).strip()
            if hdu.name != "OI_FLUX" or name not in wanted:
                continue
            nw = np.shape(hdu.data["FLAG"])[-1]
            if nw != wanted[name]:
                raise ValueError(
                    f"{sci.path}: OI_FLUX has {nw} channels but the calibrated {name} has "
                    f"{wanted[name]}; write flux only for unbinned, unwindowed calibrations"
                )
            copy = fits.BinTableHDU(data=hdu.data.copy(), header=hdu.header.copy())
            copy.header["EXTVER"] = extvers[name]
            for key in ("CHECKSUM", "DATASUM"):
                copy.header.remove(key, ignore_missing=True)
            out.append(copy)
    return out


CP_ONLY_KEY = "CPONLY"  # primary-header keyword of a closure-phase-only product


def _package_version():
    try:
        return version("virgil-vlti")
    except PackageNotFoundError:
        return "unknown"


_PRO_CARD = re.compile(r"ESO PRO (CATG|REC\d+ .*)")


def _provenance_header(header, provenance):
    """``header`` as an output product's: ESO's processing cards removed, ours added.

    The input's ``PRO CATG`` and ``PRO RECn`` cards describe ESO's pipeline
    run, not this one, so they are not inherited; the input's category and
    pipeline versions are recorded under ``HIERARCH VIRGIL INPUT``. The
    product has no ESO category: ``HIERARCH VIRGIL CALIBRATED`` marks it.
    """
    pipe = [
        f"{value}" for key, value in header.items() if re.fullmatch(r"ESO PRO REC\d+ PIPE ID", key)
    ]
    inherited = {
        "HIERARCH VIRGIL INPUT CATG": product_category(header),
        "HIERARCH VIRGIL INPUT PIPE ID": ",".join(dict.fromkeys(pipe)) or "not recorded",
    }
    for key in [k for k in header if _PRO_CARD.fullmatch(k)]:
        del header[key]
    header[OUTPUT_FLAG] = (True, "calibrated by virgil-vlti")
    header["HIERARCH VIRGIL VERSION"] = _package_version()
    ours = {f"HIERARCH VIRGIL {k}": v for k, v in provenance.items()}
    for key, value in {**inherited, **ours}.items():
        header[key] = value if isinstance(value, (bool, int, float)) else str(value)


def write_calibrated(sci, calibrated, path, history=None, flux=False, provenance=None):
    """Write one science exposure, calibrated, to ``path``.

    Parameters
    ----------
    sci : Exposure
    calibrated : list of Calibrated
        One per instrument (e.g. both GRAVITY polarisations).
    history : str, optional
        Added to the primary header's HISTORY.
    flux : bool
        Also copy the exposure's ``OI_FLUX`` tables (uncalibrated, as
        observed; :func:`flux_hdus`), so that readers asking for flux find it.
    provenance : dict, optional
        Extra ``HIERARCH VIRGIL <key>`` cards (calibrators, diameters, error
        mode): the package version and the input's category and pipeline
        versions are always recorded.
    """
    primary = fits.PrimaryHDU(header=sci.header)
    _provenance_header(primary.header, provenance or {})
    if any(cal.cp_only for cal in calibrated):
        if not all(cal.cp_only for cal in calibrated):
            raise ValueError("a file is closure-phase only for every instrument or for none")
        primary.header[CP_ONLY_KEY] = (True, "closure phases only; V2 not usable")
    if history:
        primary.header.add_history(history)
    hdus = [primary, sci.array, sci.targets]
    for k, cal in enumerate(calibrated, 1):
        hdus += calibrated_hdus(sci, cal, k)
    if flux:
        hdus += flux_hdus(sci, calibrated)
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return path


SCHEMA = "virgil_vlti.modes/1"
SCHEMA_KEY = "__schema__"
UNITS_KEY = "__units__"


def _mode_unit(key):
    """Unit of one mode array, from its ``"<insname>/<observable>/<mode>"`` name.

    The unit is inferred from the mode-name prefixes, which must stay in step with
    ``cli.RELATIVE_PREFIX`` (``rel-``), ``calibrate.transfer_function_modes`` (``tf-vis2-``,
    ``tf-t3phi-``) and ``templates.template_modes`` (``vis2-``, ``t3phi-``): a new prefix
    needs a line here, or its modes are labelled absolute.
    """
    _, obs, mode = key.split("/", 2)
    if obs == "t3phi":
        return "deg per 1 sigma"
    if obs in ("vis2", "visamp") and mode.startswith(("rel-", "vis2-", "tf-vis2-")):
        return "d log|V| per 1 sigma"
    if obs == "t3amp" and mode.startswith("rel-"):
        return "d ln T3AMP per 1 sigma"
    return "absolute change in the observable per 1 sigma"


def write_modes(path, modes):
    """Write nuisance modes beside a calibrated file, as ``.npz``.

    Each array is named ``"<insname>/<observable>/<mode>"`` and has the shape
    of that observable's table in the calibrated OIFITS file (rows ×
    channels, in the file's order), so it lines up with the data datum by
    datum. Modes with the same name in different files (one calibrator on one
    night) are the same nuisance parameter.

    Parameters
    ----------
    path : str
        The ``.npz`` file to write. Besides the mode arrays it holds
        ``__schema__`` (``virgil_vlti.modes/1``) and ``__units__`` (a JSON
        object, array name to unit; closure-phase modes are in degrees).
    modes : dict
        ``{insname: {mode: {observable: array}}}``, as from
        :func:`virgil_vlti.calibrate.diameter_modes` per instrument.
    """
    arrays = {
        f"{ins}/{obs}/{mode}": values
        for ins, by_mode in modes.items()
        for mode, by_obs in by_mode.items()
        for obs, values in by_obs.items()
    }
    units = {key: _mode_unit(key) for key in arrays}
    np.savez_compressed(
        path,
        **arrays,
        **{SCHEMA_KEY: np.array(SCHEMA), UNITS_KEY: np.array(json.dumps(units))},
    )
    return path
