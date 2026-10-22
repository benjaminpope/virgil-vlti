"""Read VLTI OIFITS products into exposures, and label them science or calibrator.

An :class:`Exposure` holds one product file: its primary header, its array and
target tables, and, per instrument (``INSNAME``), the wavelength table and the
``OI_VIS2`` and ``OI_T3`` tables. Nothing here changes the data.

Each exposure is labelled ``"SCI"`` or ``"CAL"`` by one of three rules, chosen
by :func:`classify`:

* **dual field** (GRAVITY only): the science fibre's offset from the fringe
  tracker, ``ESO INS SOBJ X/Y``. Frames with the fibre on the fringe-tracking
  star are an in-field calibrator, if that star is unresolved
  (:func:`virgil_vlti.gravity.is_on_axis`);
* **named calibrators**: exposures whose target is in a given list;
* **the product category**: ESO's ``ESO PRO CATG`` (GRAVITY ``*_CAL_VIS``
  against ``*_SCI_VIS``; MATISSE ``CALIB_RAW_INT`` against ``TARGET_RAW_INT``).

Products that are already calibrated (``*_CALIBRATED``, MATISSE's
``*_CAL_INT``, or ``VISCAL = CALIBRATED``) are refused.
"""

import contextlib
import dataclasses
import glob
import os
import re
from datetime import datetime, timedelta

import numpy as np
from astropy.io import fits

from . import gravity, matisse

# The UT hour at which one observing night ends and the next begins. Paranal
# observes from about 23:00 to 10:00 UT, so noon UT never splits a night.
NIGHT_BOUNDARY_UT_HOURS = 12.0


@dataclasses.dataclass
class Exposure:
    """One product file, read into memory."""

    path: str
    date: str  # DATE-OBS
    mjd: float  # MJD-OBS
    target: str
    header: fits.Header
    array: fits.BinTableHDU
    targets: fits.BinTableHDU
    wave: dict  # INSNAME -> EFF_WAVE (m)
    band: dict  # INSNAME -> EFF_BAND (m)
    v2: dict  # INSNAME -> OI_VIS2 table
    t3: dict  # INSNAME -> OI_T3 table
    kind: str = "SCI"
    vis: dict = dataclasses.field(default_factory=dict)  # INSNAME -> OI_VIS, if any

    @property
    def pipeline_diameter(self):
        """(diameter, error) in mas that the pipeline used for a calibrator, or None.

        GRAVITY's ``gravity_vis`` writes the uniform-disk diameter it took from
        its catalogue as ``ESO QC TF CAL_SC DIAM`` (and ``… DIAM ERR``) into
        calibrator products (``*_CAL_VIS``). That is the diameter
        ``gravity_viscal`` calibrates with; the observer's
        ``ESO INS SOBJ DIAMETER`` can differ. Science products may carry the
        same keyword, which is ignored here.
        """
        # Our own output has no ESO category; a null (a calibrator calibrated by
        # others) keeps the diameter card it inherited.
        if not (
            product_category(self.header).endswith("_CAL_VIS")
            or (
                self.header.get(OUTPUT_FLAG) is True
                and str(self.header.get("HIERARCH VIRGIL INPUT CATG", "")).endswith("_CAL_VIS")
            )
        ):
            return None
        diam = self.header.get("ESO QC TF CAL_SC DIAM")
        if diam is None or not np.isfinite(float(diam)) or float(diam) < 0:
            return None
        err = self.header.get("ESO QC TF CAL_SC DIAM ERR")
        ok = err is not None and np.isfinite(float(err)) and float(err) >= 0
        return float(diam), float(err) if ok else None

    @property
    def night(self):
        """The UT date at the start of the night (YYYY-MM-DD)."""
        t = datetime.fromisoformat(self.date[:19])
        return (t - timedelta(hours=NIGHT_BOUNDARY_UT_HOURS)).date().isoformat()

    @property
    def insnames(self):
        return list(self.v2)


def _insname(hdu):
    name = hdu.header.get("INSNAME")
    return None if name is None else str(name).strip()


def instrument_names(hdul):
    """The ``INSNAME`` of every ``OI_VIS2`` table, in file order."""
    return [_insname(h) for h in hdul if h.name == "OI_VIS2"]


def default_insnames(names):
    """The instruments to calibrate when none are named.

    A file with one instrument uses it. A GRAVITY file uses its science
    channel (``GRAVITY_SC``, or both polarisations ``_P1``/``_P2``), never
    the fringe tracker. Anything else must be chosen explicitly.
    """
    names = list(dict.fromkeys(names))
    if len(names) == 1:
        return names
    science = gravity.science_insnames(names)
    if science:
        return science
    raise ValueError(f"This file holds several instruments, {names}; choose with insname=")


def _table(hdul, extname, insname):
    for h in hdul:
        if h.name == extname and _insname(h) == insname:
            return h
    raise KeyError(f"No {extname} table with INSNAME = {insname!r}")


def read_exposure(path, insname=None):
    """Read one product file.

    Parameters
    ----------
    path : str
        An OIFITS file.
    insname : str or sequence of str, optional
        The instruments to keep; by default, :func:`default_insnames`.
    """
    with fits.open(path, memmap=False) as h:
        p = h[0].header
        names = instrument_names(h)
        if insname is None:
            wanted = default_insnames(names)
        else:
            wanted = [insname] if isinstance(insname, str) else list(insname)
        wave, band, v2, t3, vis = {}, {}, {}, {}, {}
        for name in wanted:
            w = _table(h, "OI_WAVELENGTH", name).data
            wave[name] = np.asarray(w["EFF_WAVE"], float)
            band[name] = np.asarray(w["EFF_BAND"], float)
            v2[name] = _table(h, "OI_VIS2", name).copy()
            t3[name] = _table(h, "OI_T3", name).copy()
            with contextlib.suppress(KeyError):  # OI_VIS is optional
                vis[name] = _table(h, "OI_VIS", name).copy()
        target = p.get("ESO OBS TARG NAME") or p.get("OBJECT")
        if target is None:
            target = str(h["OI_TARGET"].data["TARGET"][0])
        return Exposure(
            path=path,
            date=p["DATE-OBS"],
            mjd=float(p.get("MJD-OBS", np.nan)),
            target=str(target).strip(),
            header=p.copy(),
            array=h["OI_ARRAY"].copy(),
            targets=h["OI_TARGET"].copy(),
            wave=wave,
            band=band,
            v2=v2,
            t3=t3,
            vis=vis,
        )


def read_flux(exposure, insname, stations=False):
    """The ``OI_FLUX`` data of one instrument, read from disk, flagged values NaN.

    With ``stations``, also return each row's ``STA_INDEX``.
    """
    with fits.open(exposure.path, memmap=False) as h:
        d = _table(h, "OI_FLUX", insname).data
        column = "FLUXDATA" if "FLUXDATA" in d.names else "FLUX"
        flux = np.where(d["FLAG"], np.nan, d[column])
        return (flux, np.array(d["STA_INDEX"]).ravel()) if stations else flux


def classify(exposures, *, dual_field=False, calibrators=None, onaxis_mas=None):
    """Label each exposure ``"SCI"`` or ``"CAL"`` in place, and return them.

    With neither ``dual_field`` nor ``calibrators``, the product category decides.

    Parameters
    ----------
    exposures : list of Exposure
    dual_field : bool
        Use the GRAVITY science-fibre offset: on the fringe-tracking star
        means calibrator.
    calibrators : iterable of str, optional
        Target names of calibrator stars. Matching ignores case, spaces,
        underscores and a leading ``NAME``.
    onaxis_mas : float, optional
        For ``dual_field``: the largest fibre offset still counted as on
        axis (default :data:`virgil_vlti.gravity.ONAXIS_MAS`).
    """
    names = {normalise_name(c) for c in calibrators or ()}
    for e in exposures:
        instrument_module(e.header)  # refuses an unsupported INSTRUME
        if is_transfer_function_product(e.header):
            raise ValueError(
                f"{e.path} is a transfer-function product ({product_category(e.header)!r}), "
                "not raw visibilities; give the *_SCI_VIS and *_CAL_VIS products"
            )
        if is_calibrated(e.header):
            raise ValueError(
                f"{e.path} is already calibrated ({product_category(e.header)!r}); "
                "give the uncalibrated products"
            )
        if dual_field:
            kwargs = {} if onaxis_mas is None else {"onaxis_mas": onaxis_mas}
            cal = gravity.is_on_axis(e.header, **kwargs)
        elif names:
            cal = normalise_name(e.target) in names
        else:
            kind = kind_from_category(e.header)
            if kind is None:
                raise ValueError(
                    f"{e.path}: cannot classify ESO PRO CATG "
                    f"{product_category(e.header)!r}; name the calibrators or use dual_field=True"
                )
            cal = kind == "CAL"
        e.kind = "CAL" if cal else "SCI"
    return exposures


# Product categories (ESO PRO CATG), matched exactly or by suffix: GRAVITY
# writes SINGLE_/DUAL_ SCI_VIS and CAL_VIS (and *_CALIBRATED after
# gravity_viscal); MATISSE writes TARGET_RAW_INT and CALIB_RAW_INT and, once
# calibrated, *_CAL_INT (design/matisse_calibration_review.md, footgun 13).
_SCI_CATEGORIES = ("TARGET_RAW_INT",)
_CAL_CATEGORIES = ("CALIB_RAW_INT",)


# The module holding each instrument's header conventions (compatibility,
# setup, pipeline parameters). Any other named instrument (PIONIER) is
# refused: GRAVITY's keys and rules do not describe it. A header without
# INSTRUME is read with GRAVITY's.
INSTRUMENT_MODULES = {"GRAVITY": gravity, "MATISSE": matisse}


def instrument_of(header):
    """The product's ``INSTRUME``, stripped and upper-case (``""`` if absent)."""
    return str(header.get("INSTRUME", "")).strip().upper()


def instrument_module(header):
    """:data:`INSTRUMENT_MODULES` entry for a header; :mod:`gravity` if INSTRUME is absent.

    Raises ValueError for any other instrument (PIONIER, say): calibrating it
    with GRAVITY's rules would apply the wrong checks and conventions.
    """
    name = instrument_of(header)
    if name and name not in INSTRUMENT_MODULES:
        raise ValueError(
            f"INSTRUME {name!r} is not supported: virgil-vlti calibrates "
            f"{' and '.join(INSTRUMENT_MODULES)} products only"
        )
    return INSTRUMENT_MODULES.get(name, gravity)


_PARAM_NAME = re.compile(r"ESO PRO REC(\d+) PARAM\d+ NAME")


def pipeline_params(header, recipe=None):
    """Pipeline recipe parameters recorded in a product header, by name.

    ``recipe`` is the index i of ``ESO PRO RECi`` (GRAVITY's ``gravity_vis``
    is REC1), or None for every recipe, a later card overriding an earlier
    one of the same name. Values are returned as stored.
    """
    out = {}
    for key, name in header.items():
        m = _PARAM_NAME.fullmatch(key)
        if m and (recipe is None or int(m.group(1)) == recipe):
            out[str(name)] = header.get(key[: -len("NAME")] + "VALUE")
    return out


def product_category(header):
    return str(header.get("ESO PRO CATG", "")).strip().upper()


def is_calibrated(header):
    """Whether a product is already calibrated, so it cannot be calibrated again.

    Reads ESO's categories, this package's output (``HIERARCH VIRGIL CALIBRATED``)
    and the ``VISCAL`` card it used to write.
    """
    catg = product_category(header)
    viscal = str(header.get("VISCAL", "")).strip().upper()
    return (
        catg.endswith(("_CALIBRATED", "_CAL_INT"))
        or viscal == "CALIBRATED"
        or header.get(OUTPUT_FLAG) is True
    )


# Header cards of the calibrated files this package writes (write.py).
OUTPUT_FLAG = "HIERARCH VIRGIL CALIBRATED"


def is_transfer_function_product(header):
    """ESO's ``*_CAL_TF`` / ``*_SCI_TF`` products (``gravity_viscal``'s transfer functions)."""
    return product_category(header).endswith(("_CAL_TF", "_SCI_TF"))


def kind_from_category(header):
    """``"SCI"``, ``"CAL"`` or None, from the product category alone."""
    catg = product_category(header)
    if catg.endswith("_SCI_VIS") or catg in _SCI_CATEGORIES:
        return "SCI"
    if catg.endswith("_CAL_VIS") or catg in _CAL_CATEGORIES:
        return "CAL"
    return None


def normalise_name(name):
    """A target name for matching: case, spaces, underscores and a leading ``NAME`` ignored."""
    name = re.sub(r"[\s_]+", "", str(name)).upper()
    return name.removeprefix("NAME")


def expand_paths(inputs):
    """Expand directories to their sorted ``*.fits`` files; pass other paths through."""
    out = []
    for item in inputs:
        out += sorted(glob.glob(os.path.join(item, "*.fits"))) if os.path.isdir(item) else [item]
    return out
