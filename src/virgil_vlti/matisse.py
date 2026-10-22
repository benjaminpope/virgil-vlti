"""MATISSE-specific conventions: headers, recipes, refusals and BCD facts.

The analogue of :mod:`virgil_vlti.gravity` for MATISSE ``RAW_INT`` products
(``TARGET_RAW_INT``, ``CALIB_RAW_INT``). Stage 1 of design/matisse_plan.md:
read every fact the calibration needs from a header, and refuse what must not
be calibrated.

MATISSE is outside the maintainer's core expertise, so every convention here
cites its source, and a convention not yet verified refuses rather than
guesses (AGENTS.md, "MATISSE work rules"). Sources are cited as in
design/matisse_plan.md, §9:

* [PM]  ESO MATISSE Pipeline User Manual, VLT-MAN-MAT-15860-9310, issue 2.5.0;
* [UM]  ESO MATISSE User Manual, VLT-MAN-ESO-289268, issue 1.0;
* [Lo22] Lopez et al. 2022, A&A 659, A192;
* [GR22] Gámez Rosas et al. 2022, Nature 602, 403;
* [MT]  the MATISSE consortium tools at b44f965;
* [ESO-news] the ESO MATISSE news page, read 2026-10-07;
* "plan §n": design/matisse_plan.md, and its §7.1 list of unverified items
  for a MATISSE expert.

Unverified conventions, each refused or warned about until checked:

* which header key is which BCD device (:data:`BCD_KEY_ORDER`, plan §3.1,
  §7.1 item 1): unrecorded, so no position is named in device order unless
  both devices are in the same state;
* the N-band phase sign before 1.5.2 (:func:`phase_convention`, plan §7.1
  item 7): ``"unknown"``, refused by default;
* the chopping keyword (``ESO ISS CHOP ST``, from the ν Hor headers only,
  plan §7.1 item 5);
* the spatial-filter keywords (``ESO INS SFL/SFN NAME``, from the ν Hor
  headers only, plan §7.1 item 11);
* the ``OI_VIS2``/``OI_T3`` row order (:data:`VIS2_PEAKS`, :data:`T3_PEAKS`,
  from one ν Hor ``OUT-OUT`` product, plan §3.1).
"""

import dataclasses
import os
import re
import warnings

import numpy as np

# ---------------------------------------------------------------------------
# Bands, positions and the BCD devices
# ---------------------------------------------------------------------------

# ``ESO DET NAME`` -> band. L and M share one detector; N has its own
# (plan §1). The bands are always separate datasets.
BANDS = {"MATISSE-LM": "LM", "MATISSE-N": "N"}

# BCD positions in *device order*: "<IP1-IP3 device>-<IP5-IP7 device>", so
# ``IN-OUT`` means the device that swaps inputs IP1 and IP3 is in (plan §3.1;
# [PM] §4.1 and [UM] §2.2.1 for the swaps). This is not necessarily
# ``<BCD1>-<BCD2>`` of the header: see :data:`BCD_KEY_ORDER`.
POSITIONS = ("OUT-OUT", "IN-OUT", "OUT-IN", "IN-IN")
CODES = {"OUT-OUT": "OO", "IN-OUT": "IO", "OUT-IN": "OI", "IN-IN": "II"}

# Which header key is which device: ``"BCD1=IP13"`` (``ESO INS BCD1`` swaps
# IP1 and IP3) or ``"BCD1=IP57"``. Neither manual says. The [PM] §11.9
# closure-phase formula implies the first word of a position name is the
# IP1-IP3 device; the consortium tools ([MT] ``libPostTools.py``,
# ``mat_removeBCD``) treat ``BCD1 = IN`` as the IP5-IP7 swap (plan §3.1). It
# is a Stage 2 fact, decided on a raw product by the [PM] Table 4.1 shutter
# sequence, and until it is recorded here it is None: a position whose two
# devices differ cannot be named in device order, and nothing that needs the
# name (position codes, mode names, the BCD combination) may proceed
# (plan §7.1 item 1).
BCD_KEY_ORDER = None
BCD_KEY_ORDERS = ("BCD1=IP13", "BCD1=IP57")

# The interferometric inputs (``ESO ISS CONF INPUTk``) and the input swap
# each BCD device makes when in ([PM] §4.1; [UM] §2.2.1 gives the same swap
# as beam order 1-2-3-4 -> 2-1-4-3).
INPUTS = (1, 3, 5, 7)
DEVICE_SWAPS = {"IP13": (1, 3), "IP57": (5, 7)}

# The fringe peaks of the ``OI_VIS2`` and ``OI_T3`` rows, in file order,
# labelled by the inputs they combine. Read from one ν Hor ``OUT-OUT``
# product (``STA_INDEX`` K0=28, G2=19, D0=13, J3=24 at inputs 1, 3, 5, 7):
# plan §3.1. Unverified for other products and pipeline versions (Stage 2).
VIS2_PEAKS = ((5, 7), (1, 3), (3, 5), (3, 7), (1, 5), (1, 7))
T3_PEAKS = ((3, 5, 7), (1, 3, 5), (1, 3, 7), (1, 5, 7))

# The baseline seen on each ``OI_VIS2`` peak at each position, labelled by
# the ``OUT-OUT`` input of each telescope and in the order the telescopes
# reach the peak: plan §3.1, first table. A reversed pair is a reversed
# baseline. Rows in the order of :data:`VIS2_PEAKS`.
PEAK_BASELINES = {
    "OUT-OUT": ((5, 7), (1, 3), (3, 5), (3, 7), (1, 5), (1, 7)),
    "IN-OUT": ((5, 7), (3, 1), (1, 5), (1, 7), (3, 5), (3, 7)),
    "OUT-IN": ((7, 5), (1, 3), (3, 7), (3, 5), (1, 7), (1, 5)),
    "IN-IN": ((7, 5), (3, 1), (1, 7), (1, 5), (3, 7), (3, 5)),
}

# The physical triangle on each ``OI_T3`` peak at each position (sorted
# ``OUT-OUT`` inputs), and the sign s of the permutation that brings it to
# the order its telescopes reach the peak: the closure phase on the peak is
# s times that of the sorted triangle, plus the peak's instrumental term.
# Plan §3.1, second table; [PM] §11.9 writes the resulting combination,
# (φ_OO + φ_OI − φ_IO − φ_II)/4, for the IP1-IP3-IP5 triangle only, and the
# other three are derived (plan §7.1 item 3). Rows in the order of
# :data:`T3_PEAKS`.
PEAK_TRIANGLES = {
    "OUT-OUT": (((3, 5, 7), +1), ((1, 3, 5), +1), ((1, 3, 7), +1), ((1, 5, 7), +1)),
    "IN-OUT": (((1, 5, 7), +1), ((1, 3, 5), -1), ((1, 3, 7), -1), ((3, 5, 7), +1)),
    "OUT-IN": (((3, 5, 7), -1), ((1, 3, 7), +1), ((1, 3, 5), +1), ((1, 5, 7), -1)),
    "IN-IN": (((1, 5, 7), -1), ((1, 3, 7), -1), ((1, 3, 5), -1), ((3, 5, 7), -1)),
}

# Chopped L-band exposures are taken only at these positions ([PM] Table 4.1:
# 4 chopped at IN-IN and 4 at OUT-OUT; [Lo22] §3.1 words it differently,
# plan §7.1 item 19). A chopped LM exposure with the two devices in
# different states is outside that sequence and is warned about.
LM_CHOPPED_POSITIONS = ("OUT-OUT", "IN-IN")

# [ESO-news], October 2024: N-band frames observed in this date range were
# tagged with the wrong chopping state; the affected keyword is not named
# (plan §7.1 item 5). N-band exposures in it are warned about.
N_CHOP_TAG_SUSPECT = ("2024-05-16", "2024-08-12")


class UnverifiedConvention(ValueError):
    """A MATISSE convention this package has not verified, so it refuses."""


def _str(header, key):
    value = header.get(key)
    return None if value is None else str(value).strip().upper()


def require_bcd_key_order(bcd_key_order=None):
    """The recorded ``bcd_key_order`` fact, or a refusal naming it.

    ``bcd_key_order`` overrides the recorded :data:`BCD_KEY_ORDER` (tests,
    and Stage 2 once it has measured it). Everything that needs a position
    in device order calls this: position codes, mode names, the BCD
    combination (Stage 4).
    """
    order = BCD_KEY_ORDER if bcd_key_order is None else bcd_key_order
    if order is None:
        raise UnverifiedConvention(
            "bcd_key_order is not recorded: which of ESO INS BCD1/BCD2 is the IP1-IP3 "
            "device is unverified (design/matisse_plan.md §3.1, §7.1 item 1), so a position "
            "with one device in cannot be named, combined or calibrated across positions"
        )
    if order not in BCD_KEY_ORDERS:
        raise ValueError(f"bcd_key_order must be one of {BCD_KEY_ORDERS}, not {order!r}")
    return order


def device_position(bcd_keys, bcd_key_order=None):
    """The position in device order from the raw ``(BCD1, BCD2)`` key values.

    Both devices in the same state (``OUT-OUT``, ``IN-IN``) need no key
    order. Otherwise the ``bcd_key_order`` fact decides, and its absence
    raises :class:`UnverifiedConvention`.
    """
    b1, b2 = bcd_keys
    for value in bcd_keys:
        if value not in ("IN", "OUT"):
            raise ValueError(f"BCD state {value!r} is neither IN nor OUT")
    if b1 == b2:
        return f"{b1}-{b2}"
    order = require_bcd_key_order(bcd_key_order)
    return f"{b1}-{b2}" if order == "BCD1=IP13" else f"{b2}-{b1}"


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Setup:
    """The instrument setup of one MATISSE product, from its primary header."""

    band: str  # "LM" or "N", from ESO DET NAME
    resolution: str  # ESO INS DIL NAME (LM) or ESO INS DIN NAME (N)
    bcd: str | None  # position in device order; None if bcd_key_order is needed and absent
    bcd_keys: tuple  # raw (ESO INS BCD1 NAME, ESO INS BCD2 NAME); never used for signs
    chopped: bool  # ESO ISS CHOP ST == "T"; never ESO DET CHOP ST
    dit: float | None  # ESO DET SEQ1 DIT [s]
    telescopes: str | None  # "AT" or "UT", from ESO ISS CONF T1NAME
    spatial_filter: str | None  # ESO INS SFL NAME (LM) or SFN NAME (N); unverified keys
    fringe_tracker: str | None  # ESO DEL FT SENSOR
    template: str | None  # ESO TPL ID
    tpl_start: str | None  # ESO TPL START
    inputs: tuple  # ESO ISS CONF INPUT1..4


def setup(header, bcd_key_order=None):
    """Read the :class:`Setup` of a MATISSE product.

    Refuses (``ValueError``) a product whose band, resolution, BCD state or
    chopping state cannot be read, and one whose ``ESO CFG BCD MODE``
    disagrees with ``<BCD1>-<BCD2>`` of the INS keys (which key gives which
    device is undocumented, so a disagreement is an unresolved convention;
    plan §4, Stage 1, step 1).

    ``bcd`` is the position in device order (plan §3.1). With the two
    devices in different states it needs the ``bcd_key_order`` fact
    (:data:`BCD_KEY_ORDER`, or the argument); without it ``bcd`` is None,
    and :func:`position_code` refuses. Comparing ``bcd_keys`` needs no
    order, since equal keys are the same physical state.
    """
    det = _str(header, "ESO DET NAME")
    if det not in BANDS:
        raise ValueError(f"ESO DET NAME {det!r} is not a MATISSE detector {tuple(BANDS)}")
    band = BANDS[det]
    # The resolution of the band's own dispersive element ([PM] §9.1; plan §4,
    # Stage 1): DIL for L/M, DIN for N. The ν Hor headers carry both.
    res_key = "ESO INS DIL NAME" if band == "LM" else "ESO INS DIN NAME"
    resolution = _str(header, res_key)
    if not resolution:
        raise ValueError(f"{res_key} is missing: the {band}-band resolution is unknown")
    bcd_keys = (_str(header, "ESO INS BCD1 NAME"), _str(header, "ESO INS BCD2 NAME"))
    if None in bcd_keys:
        raise ValueError("ESO INS BCD1 NAME or BCD2 NAME is missing: the BCD state is unknown")
    for value in bcd_keys:
        if value not in ("IN", "OUT"):
            raise ValueError(f"BCD state {value!r} is neither IN nor OUT")
    cfg = _str(header, "ESO CFG BCD MODE")
    if cfg is not None and cfg != f"{bcd_keys[0]}-{bcd_keys[1]}":
        raise UnverifiedConvention(
            f"ESO CFG BCD MODE {cfg!r} disagrees with ESO INS BCD1/BCD2 NAME "
            f"{bcd_keys[0]!r}/{bcd_keys[1]!r}; which key gives which device is unverified "
            "(design/matisse_plan.md §3.1, §7.1 item 1)"
        )
    try:
        bcd = device_position(bcd_keys, bcd_key_order)
    except UnverifiedConvention:
        bcd = None
    # Chopping: ESO ISS CHOP ST, never ESO DET CHOP ST, which is False in the
    # ν Hor chopped products (plan §4, Stage 1, step 1; keys from the ν Hor
    # headers only, plan §7.1 item 5).
    chop = header.get("ESO ISS CHOP ST")
    if chop is True or str(chop).strip().upper() == "T":
        chopped = True
    elif chop is False or str(chop).strip().upper() == "F":
        chopped = False
    else:
        raise ValueError(f"ESO ISS CHOP ST is {chop!r}, not T or F: the chopping state is unknown")
    dit = header.get("ESO DET SEQ1 DIT")
    t1 = _str(header, "ESO ISS CONF T1NAME")
    telescopes = None if t1 is None else "UT" if t1.startswith("UT") else "AT"
    # Spatial filter: ESO INS SFL NAME (LM) and SFN NAME (N) in the ν Hor
    # headers ('1.50' and '2.00', λ/D, as [Lo22] Table 2's standard pinholes);
    # neither manual names the keys (plan §7.1 item 11). Compared only.
    spatial = _str(header, "ESO INS SFL NAME" if band == "LM" else "ESO INS SFN NAME")
    inputs = tuple(header.get(f"ESO ISS CONF INPUT{k}") for k in range(1, 5))
    s = Setup(
        band=band,
        resolution=resolution,
        bcd=bcd,
        bcd_keys=bcd_keys,
        chopped=chopped,
        dit=None if dit is None else float(dit),
        telescopes=telescopes,
        spatial_filter=spatial,
        fringe_tracker=_str(header, "ESO DEL FT SENSOR"),
        template=_str(header, "ESO TPL ID"),
        tpl_start=None if header.get("ESO TPL START") is None else str(header["ESO TPL START"]),
        inputs=inputs,
    )
    for message in setup_warnings(header, s):
        warnings.warn(message, stacklevel=2)
    return s


def setup_warnings(header, s):
    """Warnings about a product's setup that do not refuse it."""
    out = []
    if s.band == "LM" and s.chopped and s.bcd_keys[0] != s.bcd_keys[1]:
        out.append(
            f"chopped LM exposure at BCD {s.bcd_keys[0]}/{s.bcd_keys[1]}: [PM] Table 4.1 "
            f"chops only at {' and '.join(LM_CHOPPED_POSITIONS)}"
        )
    date = str(header.get("DATE-OBS", ""))[:10]
    if s.band == "N" and N_CHOP_TAG_SUSPECT[0] <= date <= N_CHOP_TAG_SUSPECT[1]:
        out.append(
            f"N-band exposure of {date}: ESO reports wrongly tagged chopping states between "
            f"{N_CHOP_TAG_SUSPECT[0]} and {N_CHOP_TAG_SUSPECT[1]} ([ESO-news], October 2024)"
        )
    return out


def position_code(s, bcd_key_order=None):
    """``"OO"``, ``"IO"``, ``"OI"`` or ``"II"`` for a :class:`Setup`.

    Refuses (:class:`UnverifiedConvention`) a position with one device in
    when ``bcd_key_order`` is not recorded.
    """
    return CODES[s.bcd if s.bcd is not None else device_position(s.bcd_keys, bcd_key_order)]


def peak_baseline(position, peak):
    """The baseline seen on ``OI_VIS2`` row ``peak`` at ``position`` (device order).

    In ``OUT-OUT`` inputs, in the order the telescopes reach the peak
    (plan §3.1). ``peak`` is a row index or an input pair of
    :data:`VIS2_PEAKS`.
    """
    row = VIS2_PEAKS.index(tuple(peak)) if not isinstance(peak, int) else peak
    return PEAK_BASELINES[position][row]


def peak_triangle(position, peak):
    """(sorted triangle, sign) seen on ``OI_T3`` row ``peak`` at ``position``.

    The closure phase on the peak is ``sign`` times that of the sorted
    triangle (in ``OUT-OUT`` inputs) plus the peak's instrumental term
    (plan §3.1). ``peak`` is a row index or an input triple of
    :data:`T3_PEAKS`.
    """
    row = T3_PEAKS.index(tuple(peak)) if not isinstance(peak, int) else peak
    return PEAK_TRIANGLES[position][row]


# ---------------------------------------------------------------------------
# Recipes and versions
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Recipe:
    """One ``ESO PRO RECi`` entry of a product header."""

    index: int  # i in ESO PRO RECi
    name: str  # ESO PRO RECi ID
    pipe_id: str | None  # ESO PRO RECi PIPE ID, e.g. "matisse/2.0.0"
    version: tuple | None  # (major, minor, patch) from pipe_id
    params: dict  # PARAMn NAME -> PARAMn VALUE, as stripped strings


_REC = re.compile(r"ESO PRO REC(\d+) ID")
_PARAM = re.compile(r"ESO PRO REC(\d+) PARAM(\d+) NAME")


def parse_version(pipe_id):
    """``"matisse/1.7.9"`` -> (1, 7, 9); anything without x.y.z -> None."""
    m = re.search(r"(\d+)\.(\d+)\.(\d+)", str(pipe_id or ""))
    return None if m is None else tuple(int(g) for g in m.groups())


def recipes(header):
    """Every ``ESO PRO RECi`` recipe, by name; the lowest i wins a repeated name.

    MATISSE products carry several recipes (``mat_raw_estimates``, then
    ``mat_cal_oifits`` and others), each with its own parameters (plan §2).
    """
    found = {}
    for key in header:
        m = _REC.fullmatch(key)
        if m:
            found[int(m.group(1))] = str(header[key]).strip()
    params = {i: {} for i in found}
    for key in header:
        m = _PARAM.fullmatch(key)
        if m and int(m.group(1)) in params:
            value = header.get(key[: -len("NAME")] + "VALUE")
            params[int(m.group(1))][str(header[key]).strip()] = (
                None if value is None else str(value).strip()
            )
    out = {}
    for i in sorted(found):
        if found[i] in out:
            continue
        pipe = header.get(f"ESO PRO REC{i} PIPE ID")
        pipe = None if pipe is None else str(pipe).strip()
        out[found[i]] = Recipe(i, found[i], pipe, parse_version(pipe), params[i])
    return out


def raw_estimates(header):
    """The ``mat_raw_estimates`` recipe of a product, or None."""
    return recipes(header).get("mat_raw_estimates")


# ---------------------------------------------------------------------------
# The N-band phase sign
# ---------------------------------------------------------------------------

# N-band phases and closure phases of mat_raw_estimates in [1.5.2, 2.0.0) are
# rotated by 180° relative to L and M; 2.0.0 corrected it ([PM] §2.5, the
# 2.0.0 release notes; [GR22] Methods §4.3, §4.8; [ESO-news], December 2023).
# The lower bound is not confirmed by the 1.5.2 notes ([PM] §2.11), so
# earlier versions are "unknown", not "standard" (plan §7, decision 2,
# revised; §7.1 item 7).
N_PHASE_FLIPPED = ((1, 5, 2), (2, 0, 0))

# Above this median wavelength (m) a product without ESO DET NAME is N band.
N_BAND_MIN_WAVE = 7e-6


def phase_convention(header, wave=None):
    """``"standard"``, ``"flipped"`` or ``"unknown"`` for a product's phases.

    L and M are always ``"standard"``. In N, a ``mat_raw_estimates`` version
    in [1.5.2, 2.0.0) is ``"flipped"``, one at or after 2.0.0
    ``"standard"``, and anything else (no version, or before 1.5.2)
    ``"unknown"``. Without ``ESO DET NAME``, N band is a median wavelength
    above 7 µm (plan §4, Stage 1, step 3).
    """
    det = _str(header, "ESO DET NAME")
    if det in BANDS:
        n_band = BANDS[det] == "N"
    elif wave is not None and np.size(wave):
        n_band = float(np.median(wave)) > N_BAND_MIN_WAVE
    else:
        return "unknown"
    if not n_band:
        return "standard"
    rec = raw_estimates(header)
    version = None if rec is None else rec.version
    if version is None or version < N_PHASE_FLIPPED[0]:
        return "unknown"
    return "flipped" if version < N_PHASE_FLIPPED[1] else "standard"


N_PHASE_CHOICES = ("refuse", "negate", "keep")


def n_phase_action(convention, n_phase="refuse", n_phase_unknown="refuse"):
    """What to do with a product's phases: ``"keep"``, ``"negate"`` or ``"refuse"``.

    ``n_phase`` applies only to ``"flipped"`` products and ``n_phase_unknown``
    only to ``"unknown"`` ones, so negating flipped data never rotates older
    data unasked (plan §4, Stage 1, CLI). ``"standard"`` is always kept.
    """
    for name, value in (("n_phase", n_phase), ("n_phase_unknown", n_phase_unknown)):
        if value not in N_PHASE_CHOICES:
            raise ValueError(f"{name} must be one of {N_PHASE_CHOICES}, not {value!r}")
    if convention == "flipped":
        return n_phase
    if convention == "unknown":
        return n_phase_unknown
    return "keep"


def negate_phases(exposure, version=None):
    """Negate ``T3PHI`` and ``VISPHI`` of an exposure in place, with a HISTORY card.

    The same as (u, v) -> (−u, −v): the image rotated by 180° (plan §4,
    Stage 1, step 3). That VISPHI flips with T3PHI follows from the rotation
    ([PM] §2.5 names phases and closure phases; plan §7.1 item 7).

    Refuses a product that already carries this HISTORY card, so phases are
    never negated twice (e.g. a virgil-vlti output read back).
    """
    if was_negated(exposure.header):
        raise ValueError(f"{exposure.path}: N-band phases were already negated by virgil-vlti")
    for table in exposure.t3.values():
        table.data["T3PHI"] = -table.data["T3PHI"]
    for table in exposure.vis.values():
        table.data["VISPHI"] = -table.data["VISPHI"]
    exposure.header.add_history(f"{NEGATED} (mat_raw_estimates {version or 'unknown version'})")
    return exposure


NEGATED = "virgil-vlti: N-band phases negated"


def was_negated(header):
    """Whether virgil-vlti has negated this product's N-band phases (its HISTORY says so)."""
    history = header.get("HISTORY")
    return history is not None and any(str(h).startswith(NEGATED) for h in history)


# ---------------------------------------------------------------------------
# Merged products
# ---------------------------------------------------------------------------

# A table INT_TIME longer than this many exposure times is a merge (plan §4,
# Stage 1, step 4; ν Hor: 120.46 s against ESO DET SEQ1 EXPTIME = 61.8 s).
MERGE_INT_TIME_FACTOR = 1.5

# Whether one baseline at several MJDs refuses a product. Off until Stage 2
# has counted rows and MJDs per baseline in a single RAW_INT product (plan
# §4, Stage 1, step 4): until then it is evidence, reported, not a refusal.
MULTI_MJD_IS_FATAL = False

# File names of consortium merges and calibrated products ([MT]
# mat_mergeAllOiFits.py output names; ESO *_CAL_INT): plan §4, Stage 1, step 4.
_MERGE_NAMES = ("_CAL_OIFITS_", "_CHOP", "_NOCHOP", "_CAL_INT")
_MULTI_MJD = "one baseline at several MJDs"


def merge_evidence(header, tables=None, path=None):
    """Every sign that a product merges several exposures; empty if none.

    ``tables`` are ``OI_VIS2`` tables (FITS records with ``MJD``,
    ``STA_INDEX`` and ``INT_TIME``); ``path`` is the file name. A count of
    ``ESO PRO REC1 RAWn`` cannot detect a merge: a single ν Hor product
    lists all 12 template exposures and 2 skies (plan §4, Stage 1, step 4).
    """
    out = []
    exptime = header.get("ESO DET SEQ1 EXPTIME")
    for t in tables or ():
        names = t.columns.names if hasattr(t, "columns") else list(t)
        if exptime is not None and "INT_TIME" in names:
            longest = float(np.max(t["INT_TIME"]))
            if longest > MERGE_INT_TIME_FACTOR * float(exptime):
                out.append(
                    f"INT_TIME {longest:.2f} s > {MERGE_INT_TIME_FACTOR} x "
                    f"ESO DET SEQ1 EXPTIME {float(exptime):.2f} s"
                )
        mjds = {}
        for sta, mjd in zip(np.asarray(t["STA_INDEX"]).tolist(), t["MJD"], strict=True):
            mjds.setdefault(tuple(sorted(sta)), set()).add(float(mjd))
        repeated = sorted(k for k, v in mjds.items() if len(v) > 1)
        if repeated:
            out.append(f"{_MULTI_MJD}: {repeated}")
    for key, value in (
        ("PIPEFILE", header.get("PIPEFILE")),
        ("ARCFILE", header.get("ARCFILE")),
        ("file name", None if path is None else os.path.basename(str(path))),
    ):
        text = str(value or "").upper()
        hits = [n for n in _MERGE_NAMES if n in text]
        if hits:
            out.append(f"{key} {value!s} contains {hits[0].lower()}")
    history = header.get("HISTORY")
    if history is not None and any("mat_mergeOifits" in str(h) for h in history):
        out.append("HISTORY names mat_mergeOifits")
    return list(dict.fromkeys(out))


def _fatal(evidence):
    return [e for e in evidence if MULTI_MJD_IS_FATAL or not e.startswith(_MULTI_MJD)]


def is_merged(header, tables=None, path=None):
    """Whether a product merges several exposures (see :func:`merge_evidence`)."""
    return bool(_fatal(merge_evidence(header, tables, path)))


# ---------------------------------------------------------------------------
# Compatibility and groups
# ---------------------------------------------------------------------------

# Setup fields that must agree between science and calibrator, for every
# observable (plan §4, Stage 1, step 5). ``chopped`` may be pooled for
# closure phases only with ``pool_chopping_t3`` (plan §3.3, revised).
# ``dit`` is also a hard split in :func:`group_key`; it is compared here too
# because the CLI does not group by :func:`group_key` until Stage 4, and DIT
# changes V² through jitter and coherence loss.
MODE_KEYS = ("band", "resolution", "bcd_keys", "telescopes", "spatial_filter", "chopped", "dit")

# mat_raw_estimates parameters that must agree. 2.5.0 renamed spectralBinning
# -> spectralAverage and coherentAlgo -> cpAlgo ([PM] §2.1), and cpAlgo=2
# means coherent closure phases, not the bispectrum. A product with the old
# names and one with the new are a reduction mismatch, never silently
# equated. replaceTel replaces one telescope's photometry ([PM] §9.1.4).
REDUCTION_PARAMS = (
    "corrFlux",
    "coherentAlgo",
    "cpAlgo",
    "spectralBinning",
    "spectralAverage",
    "replaceTel",
)

# The version of mat_raw_estimates, read as gravity.WARN_KEYS reads GRAVITY's.
WARN_KEYS = ("ESO PRO REC1 PIPE ID",)

OBSERVABLES = ("vis2", "t3amp", "t3phi", "visamp", "visphi")


def _check_observable(observable):
    if observable not in OBSERVABLES:
        raise ValueError(f"observable must be one of {OBSERVABLES}, not {observable!r}")


def _pipe_id(header):
    """The mat_raw_estimates version string: REC1's, or its own entry's."""
    rec = raw_estimates(header)
    if rec is not None and rec.index != 1:
        return rec.pipe_id
    return header.get(WARN_KEYS[0])


def compatibility(sci_header, cal_header, observable="vis2", pool_chopping_t3=False):
    """How a calibrator's header differs from the science exposure's.

    Returns three lists of messages, as :func:`virgil_vlti.gravity.compatibility`:
    mode differences (band, resolution, BCD state, telescopes, spatial
    filter, chopping, DIT; and, in N, the phase convention, since calibration
    subtracts phases), mat_raw_estimates parameter differences, and the
    pipeline version.

    DIT is also a hard split in :func:`group_key`, but ``virgil-vlti-calibrate``
    groups by night only until Stage 4 wires :func:`group_key` into the CLI,
    so it is compared here as a mode difference (fatal). Until then a night
    with several BCD positions fails compatibility on ``bcd_keys`` for every
    science exposure: safe, but such nights cannot be calibrated yet.
    """
    _check_observable(observable)
    a, b = setup(sci_header), setup(cal_header)
    mode, reduction, minor = [], [], []
    for key in MODE_KEYS:
        if key == "chopped" and observable == "t3phi" and pool_chopping_t3:
            continue
        x, y = getattr(a, key), getattr(b, key)
        if x != y:
            mode.append(f"{key}: science {x!r}, calibrator {y!r}")
    if "N" in (a.band, b.band):
        x, y = phase_convention(sci_header), phase_convention(cal_header)
        if x != y:
            mode.append(f"N-band phase convention: science {x!r}, calibrator {y!r}")
    p, q = raw_estimates(sci_header), raw_estimates(cal_header)
    p, q = (p.params if p else {}), (q.params if q else {})
    for name in REDUCTION_PARAMS:
        if p.get(name) != q.get(name):
            reduction.append(
                f"mat_raw_estimates {name}: science {p.get(name)!r}, calibrator {q.get(name)!r}"
            )
    x, y = _pipe_id(sci_header), _pipe_id(cal_header)
    if x != y:
        minor.append(f"mat_raw_estimates PIPE ID: science {x!r}, calibrator {y!r}")
    return mode, reduction, minor


def group_key(exposure, observable="vis2", pool_chopping_t3=False):
    """What may calibrate what: (night, band, resolution, DIT, BCD, chopped).

    Stricter than the consortium's ``mat_autoCalib`` ([MT] ``findClosestCal``,
    which matches BCD, band, DIT and chopping, within 3 h, but not
    resolution): plan §3.3. ``chopped`` is in the key for every observable,
    and is dropped for ``"t3phi"`` only with ``pool_chopping_t3`` (plan §3.3,
    revised). The BCD state is the raw ``(BCD1, BCD2)`` pair, which needs no
    ``bcd_key_order``: equal pairs are the same physical state.

    Not yet used by ``virgil-vlti-calibrate``, ``calpca`` or ``pairs``: Stage 4
    wires it into the CLI (and Stage 5 into ``calpca``). Until then
    :func:`compatibility` refuses calibrators that differ in any of these.
    """
    _check_observable(observable)
    s = setup(exposure.header)
    pooled = observable == "t3phi" and pool_chopping_t3
    return (
        exposure.night,
        s.band,
        s.resolution,
        s.dit,
        s.bcd_keys,
        None if pooled else s.chopped,
    )


# ---------------------------------------------------------------------------
# Refusals before calibration
# ---------------------------------------------------------------------------


def _truthy(value):
    text = str(value).strip().upper()
    if text in ("TRUE", "T", "1"):
        return True
    if text in ("FALSE", "F", "0"):
        return False
    return None


def screen(exposure, *, n_phase="refuse", n_phase_unknown="refuse"):
    """Refuse a MATISSE exposure that must not be calibrated; fix its phase sign.

    Raises ``ValueError`` for a merged product, and for one whose
    ``corrFlux`` is true or not recorded (its ``VIS2DATA`` may hold squared
    correlated fluxes, [PM] §7; review footgun 1; only Stage 7 reads them).
    N-band phases follow :func:`n_phase_action`: refused, negated in place
    (:func:`negate_phases`), or kept. Returns the exposure.
    """
    h = exposure.header
    name = str(exposure.path)
    evidence = merge_evidence(h, [t.data for t in exposure.v2.values()], name)
    fatal = _fatal(evidence)
    if fatal:
        raise ValueError(f"{name} is a merged product ({'; '.join(fatal)}); give RAW_INT products")
    for note in set(evidence) - set(fatal):
        warnings.warn(f"{name}: possible merge, {note}", stacklevel=2)
    s = setup(h)
    rec = raw_estimates(h)
    corr = None if rec is None else _truthy(rec.params.get("corrFlux"))
    if corr is not False:
        recorded = None if rec is None else rec.params.get("corrFlux")
        raise ValueError(
            f"{name}: mat_raw_estimates corrFlux is {recorded!r}; "
            "corrFlux=true products hold correlated fluxes, not visibilities, and are refused "
            "for V², T3AMP and VISAMP"
        )
    wave = next(iter(exposure.wave.values()), None)
    convention = phase_convention(h, wave)
    action = n_phase_action(convention, n_phase, n_phase_unknown)
    if action == "refuse":
        option = "--n-phase" if convention == "flipped" else "--n-phase-unknown"
        raise UnverifiedConvention(
            f"{name}: {s.band}-band phases of mat_raw_estimates "
            f"{rec.pipe_id if rec else 'of unknown version'} are {convention!r} "
            f"([PM] §2.5; design/matisse_plan.md §7, decision 2); {option} negate or keep "
            "to proceed"
        )
    if action == "negate":
        negate_phases(exposure, rec.pipe_id if rec else None)
    return exposure


# ---------------------------------------------------------------------------
# Inspection
# ---------------------------------------------------------------------------


def describe(header, tables=None, path=None, wave=None):
    """A summary of one product, and whether it may be calibrated here.

    Returns a dict with the setup, pipeline version, phase convention,
    merge evidence, ``replaceTel`` and a ``verdict``: ``"ok"`` or the
    reasons it is refused (as the defaults of :func:`screen` would).
    """
    from .products import is_calibrated, product_category

    out = {
        "file": path,
        "instrument": _str(header, "INSTRUME"),
        "category": product_category(header),
        "calibrated": is_calibrated(header),
    }
    refusals, notes = [], []
    if out["calibrated"]:
        refusals.append("calibrated")
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            s = setup(header)
        notes += [str(w.message) for w in caught]
        out.update(
            band=s.band,
            resolution=s.resolution,
            bcd=s.bcd or f"BCD1={s.bcd_keys[0]} BCD2={s.bcd_keys[1]} (device order unverified)",
            chopped=s.chopped,
        )
    except ValueError as err:
        refusals.append(str(err))
    rec = raw_estimates(header)
    out["pipeline"] = None if rec is None else rec.pipe_id
    out["phase_convention"] = phase_convention(header, wave)
    if n_phase_action(out["phase_convention"]) == "refuse":
        refusals.append(f"N-band phases {out['phase_convention']}")
    evidence = merge_evidence(header, tables, path)
    out["merged"] = bool(_fatal(evidence))
    out["merge_evidence"] = evidence
    if out["merged"]:
        refusals.append("merged")
    corr = None if rec is None else rec.params.get("corrFlux")
    if _truthy(corr) is not False:
        refusals.append(f"corrFlux={corr}")
    replace = None if rec is None else rec.params.get("replaceTel")
    out["replaceTel"] = replace
    if replace not in (None, "0"):
        # [PM] §9.1.4: a two-digit code (UT digit, AT digit); a digit k != 0
        # replaces telescope k's photometry by the mean of the other three.
        notes.append(f"replaceTel={replace}: one telescope's photometry replaced ([PM] §9.1.4)")
    out["notes"] = notes
    out["verdict"] = "ok" if not refusals else "refused: " + "; ".join(refusals)
    return out
