"""MATISSE headers, refusals and BCD facts (design/matisse_plan.md, §4, Stage 1).

Headers only: the ν Hor header is a dict literal of cards read from the
consortium-merged ν Hor products (2023-12-12, LM LOW, ATs K0 G2 D0 J3); no
data file is committed. The BCD tests check the §3.1 tables against the
input swaps of [PM] §4.1 and against a binary model, computed without the
tables.
"""

import itertools
import types
import warnings

import numpy as np
import pytest
from astropy.io import fits

from virgil_vlti import matisse
from virgil_vlti.calibrate import check_compatible
from virgil_vlti.inspect_cli import main as inspect_main
from virgil_vlti.products import (
    Exposure,
    instrument_module,
    instrument_of,
    pipeline_params,
)

# Cards of the ν Hor product 2023-12-12T004705 (merged by the consortium
# tools, so INT_TIME is 120.46 s and the category TARGET_CAL_INT).
NU_HOR = {
    "INSTRUME": "MATISSE",
    "OBJECT": "nu. Hor",
    "MJD-OBS": 60290.03948737,
    "DATE-OBS": "2023-12-12T00:56:51.7088",
    "PIPEFILE": "TARGET_CAL_INT_0006.fits",
    "ARCFILE": "MATIS.2023-12-12T00:56:51.709.fits",
    "ESO TPL ID": "MATISSE_hyb_obs_ft",
    "ESO TPL START": "2023-12-12T00:47:05",
    "ESO INS BCD1 NAME": "OUT",
    "ESO INS BCD2 NAME": "OUT",
    "ESO INS DIL NAME": "LOW",
    "ESO INS DIN NAME": "LOW",
    "ESO INS SFL NAME": "1.50",
    "ESO INS SFN NAME": "2.00",
    "ESO DET CHOP ST": False,
    "ESO DET NAME": "MATISSE-LM",
    "ESO DET SEQ1 DIT": 0.0751997,
    "ESO DET SEQ1 EXPTIME": 61.828,
    "ESO PRO CATG": "TARGET_CAL_INT",
    "ESO PRO REC1 ID": "mat_raw_estimates",
    "ESO PRO REC1 PIPE ID": "matisse/2.0.0",
    "ESO PRO REC1 PARAM9 NAME": "hampelFilterKernel",
    "ESO PRO REC1 PARAM9 VALUE": "10",
    "ESO PRO REC1 PARAM10 NAME": "replaceTel",
    "ESO PRO REC1 PARAM10 VALUE": "3",
    "ESO PRO REC1 PARAM13 NAME": "corrFlux",
    "ESO PRO REC1 PARAM13 VALUE": "false",
    "ESO PRO REC1 PARAM14 NAME": "cumulBlock",
    "ESO PRO REC1 PARAM14 VALUE": "false",
    "ESO PRO REC1 PARAM15 NAME": "coherentAlgo",
    "ESO PRO REC1 PARAM15 VALUE": "2",
    "ESO PRO REC1 PARAM16 NAME": "spectralBinning",
    "ESO PRO REC1 PARAM16 VALUE": "5",
    "ESO PRO REC2 ID": "mat_cal_oifits",
    "ESO PRO REC2 PIPE ID": "matisse/1.7.9",
    "ESO PRO REC2 PARAM1 NAME": "tfKeep",
    "ESO PRO REC2 PARAM1 VALUE": "0",
    "ESO PRO REC2 PARAM2 NAME": "tfInterp",
    "ESO PRO REC2 PARAM2 VALUE": "0",
    "ESO PRO REC2 PARAM3 NAME": "cumulBlock",
    "ESO PRO REC2 PARAM3 VALUE": "true",
    "ESO PRO REC2 PARAM4 NAME": "bcdMode",
    "ESO PRO REC2 PARAM4 VALUE": "5",
    "ESO PRO REC3 ID": "mat_cal_oifits",
    "ESO PRO REC3 PIPE ID": "matisse/1.7.9",
    "ESO PRO REC3 PARAM1 NAME": "tfKeep",
    "ESO PRO REC3 PARAM1 VALUE": "1",
    "ESO DEL FT SENSOR": "GRAVITY",
    "ESO ISS CHOP ST": "T",
    "ESO ISS CONF INPUT1": 1,
    "ESO ISS CONF INPUT2": 3,
    "ESO ISS CONF INPUT3": 5,
    "ESO ISS CONF INPUT4": 7,
    "ESO ISS CONF STATION1": "K0",
    "ESO ISS CONF STATION2": "G2",
    "ESO ISS CONF STATION3": "D0",
    "ESO ISS CONF STATION4": "J3",
    "ESO ISS CONF T1NAME": "AT1",
    "ESO CFG BCD MODE": "OUT-OUT",
}
NU_HOR_INT_TIME = 120.45799994
NU_HOR_STA = ((13, 24), (28, 19), (19, 13), (19, 24), (28, 13), (28, 24))

# A single-exposure TARGET_RAW_INT product, as Stage 2 will make them.
RAW = {
    **NU_HOR,
    "ESO PRO CATG": "TARGET_RAW_INT",
    "PIPEFILE": "TARGET_RAW_INT_0001.fits",
    "ESO PRO REC1 PARAM10 VALUE": "0",
}
for _k in [k for k in RAW if k.startswith(("ESO PRO REC2", "ESO PRO REC3"))]:
    del RAW[_k]


def header(base=RAW, drop=(), **cards):
    """A FITS header from ``base``, with ``cards`` (spaces as ``_``) set and ``drop`` removed."""
    h = fits.Header()
    for k, v in {**base, **{k.replace("_", " "): v for k, v in cards.items()}}.items():
        if k not in drop:
            h[f"HIERARCH {k}" if k.startswith("ESO ") else k] = v
    return h


def n_band(version, base=RAW, **cards):
    """An N-band product header of a given mat_raw_estimates version (None: no PIPE ID)."""
    extra = {"ESO DET NAME": "MATISSE-N", **cards}
    if version is None:
        return header(base, drop=("ESO PRO REC1 PIPE ID",), **extra)
    return header(base, **{"ESO PRO REC1 PIPE ID": f"matisse/{version}"}, **extra)


def keyed(**cards):
    """A RAW header with cards given by their real names."""
    h = header(RAW)
    for k, v in cards.items():
        h[f"HIERARCH {k}"] = v
    return h


def exposure(h, *, mjds=None, int_time=61.828, wave=None, path="TARGET_RAW_INT_0001.fits"):
    """A minimal in-memory exposure: OI_VIS2 (6 rows), OI_T3 (4 rows), OI_VIS."""
    wave = np.linspace(4.2e-6, 2.8e-6, 8) if wave is None else wave
    nw = wave.size
    mjds = np.full(6, 60290.04) if mjds is None else np.asarray(mjds, float)
    v2 = fits.BinTableHDU.from_columns(
        [
            fits.Column("MJD", "1D", array=mjds),
            fits.Column("INT_TIME", "1D", array=np.full(mjds.size, int_time)),
            fits.Column("STA_INDEX", "2I", array=np.resize(NU_HOR_STA, (mjds.size, 2))),
            fits.Column("VIS2DATA", f"{nw}D", array=np.full((mjds.size, nw), 0.5)),
        ]
    )
    t3 = fits.BinTableHDU.from_columns(
        [fits.Column("T3PHI", f"{nw}D", array=np.tile(np.arange(nw, dtype=float), (4, 1)))]
    )
    vis = fits.BinTableHDU.from_columns(
        [fits.Column("VISPHI", f"{nw}D", array=np.full((6, nw), 7.0))]
    )
    return Exposure(
        path=path,
        date=str(h["DATE-OBS"]),
        mjd=float(h["MJD-OBS"]),
        target="nu. Hor",
        header=h,
        array=None,
        targets=None,
        wave={"MATISSE": wave},
        band={"MATISSE": np.full(nw, 2e-8)},
        v2={"MATISSE": v2},
        t3={"MATISSE": t3},
        vis={"MATISSE": vis},
    )


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------


def test_setup_reads_lm_low_from_dil():
    s = matisse.setup(keyed(**{"ESO INS DIL NAME": "LOW", "ESO INS DIN NAME": "HIGH"}))
    assert (s.band, s.resolution) == ("LM", "LOW")
    s = matisse.setup(header(NU_HOR))
    assert s.telescopes == "AT"
    assert s.spatial_filter == "1.50"
    assert s.fringe_tracker == "GRAVITY"
    assert s.inputs == (1, 3, 5, 7)
    assert s.dit == pytest.approx(0.0751997)
    assert s.tpl_start == "2023-12-12T00:47:05"


def test_setup_reads_n_resolution_from_din():
    h = keyed(
        **{"ESO DET NAME": "MATISSE-N", "ESO INS DIL NAME": "LOW", "ESO INS DIN NAME": "HIGH"}
    )
    s = matisse.setup(h)
    assert (s.band, s.resolution, s.spatial_filter) == ("N", "HIGH", "2.00")


def test_setup_refuses_unknown_detector_and_missing_cards():
    with pytest.raises(ValueError, match="DET NAME"):
        matisse.setup(keyed(**{"ESO DET NAME": "HAWAII-2RG"}))
    with pytest.raises(ValueError, match="DIL NAME"):
        matisse.setup(header(RAW, drop=("ESO INS DIL NAME",)))
    with pytest.raises(ValueError, match="CHOP ST"):
        matisse.setup(header(RAW, drop=("ESO ISS CHOP ST",)))
    with pytest.raises(ValueError, match="BCD"):
        matisse.setup(header(RAW, drop=("ESO INS BCD2 NAME",)))


def test_chopping_read_from_iss_not_det():
    s = matisse.setup(keyed(**{"ESO ISS CHOP ST": "T", "ESO DET CHOP ST": False}))
    assert s.chopped is True
    s = matisse.setup(keyed(**{"ESO ISS CHOP ST": "F", "ESO DET CHOP ST": True}))
    assert s.chopped is False


def _bcd(b1, b2, **cards):
    """A RAW header at BCD keys (b1, b2); not chopped, as [PM] Table 4.1 has single-device
    positions."""
    return keyed(**{"ESO INS BCD1 NAME": b1, "ESO INS BCD2 NAME": b2,
                    "ESO CFG BCD MODE": f"{b1}-{b2}", "ESO ISS CHOP ST": "F", **cards})  # fmt: skip


def test_bcd_position_from_ins_keys():
    s = matisse.setup(_bcd("IN", "OUT"), bcd_key_order="BCD1=IP13")
    assert s.bcd == "IN-OUT"
    assert matisse.position_code(s) == "IO"
    for b in ("IN", "OUT"):  # equal states need no key order
        s = matisse.setup(_bcd(b, b))
        assert s.bcd == f"{b}-{b}"
        assert matisse.position_code(s) == matisse.CODES[f"{b}-{b}"]


def test_bcd_position_in_device_order():
    s = matisse.setup(_bcd("IN", "OUT"), bcd_key_order="BCD1=IP57")
    assert (s.bcd, s.bcd_keys) == ("OUT-IN", ("IN", "OUT"))
    s = matisse.setup(_bcd("IN", "OUT"), bcd_key_order="BCD1=IP13")
    assert (s.bcd, s.bcd_keys) == ("IN-OUT", ("IN", "OUT"))
    with pytest.raises(ValueError, match="bcd_key_order must be"):
        matisse.device_position(("IN", "OUT"), "BCD1=IN")


def test_cfg_bcd_mode_disagreement_refused():
    h = keyed(
        **{"ESO INS BCD1 NAME": "IN", "ESO INS BCD2 NAME": "OUT", "ESO CFG BCD MODE": "OUT-IN"}
    )
    with pytest.raises(matisse.UnverifiedConvention, match="'OUT-IN'.*'IN'/'OUT'"):
        matisse.setup(h)


def test_bcd_key_order_required():
    assert matisse.BCD_KEY_ORDER is None  # a Stage 2 fact, not yet recorded
    s = matisse.setup(_bcd("OUT", "IN"))
    assert s.bcd is None and s.bcd_keys == ("OUT", "IN")
    with pytest.raises(matisse.UnverifiedConvention, match="bcd_key_order"):
        matisse.position_code(s)
    with pytest.raises(matisse.UnverifiedConvention, match="bcd_key_order"):
        matisse.require_bcd_key_order()
    # Calibration across positions is a mode mismatch, with or without the fact.
    mode, _, _ = matisse.compatibility(_bcd("IN", "OUT"), _bcd("OUT", "IN"))
    assert any(m.startswith("bcd_keys") for m in mode)
    mode, _, _ = matisse.compatibility(_bcd("IN", "IN"), _bcd("OUT", "OUT"))
    assert any(m.startswith("bcd_keys") for m in mode)
    assert matisse.compatibility(_bcd("IN", "OUT"), _bcd("IN", "OUT"))[0] == []


def test_lm_chopping_outside_the_sequence_warns():
    chopped = {"ESO ISS CHOP ST": "T"}
    with pytest.warns(UserWarning, match="Table 4.1"):
        matisse.setup(_bcd("IN", "OUT", **chopped))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        matisse.setup(_bcd("IN", "IN", **chopped))
        matisse.setup(_bcd("IN", "OUT"))


def test_n_band_chopping_tag_period_warns():
    h = n_band("2.0.0", **{"DATE-OBS": "2024-06-01T03:00:00"})
    with pytest.warns(UserWarning, match="ESO-news"):
        matisse.setup(h)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        matisse.setup(n_band("2.0.0", **{"DATE-OBS": "2024-09-01T03:00:00"}))


# ---------------------------------------------------------------------------
# Recipes, versions and the N-band phase sign
# ---------------------------------------------------------------------------


def test_recipes_parsed_per_recipe():
    r = matisse.recipes(header(NU_HOR))
    assert set(r) == {"mat_raw_estimates", "mat_cal_oifits"}
    raw, cal = r["mat_raw_estimates"], r["mat_cal_oifits"]
    assert (raw.index, raw.pipe_id, raw.version) == (1, "matisse/2.0.0", (2, 0, 0))
    assert raw.params["spectralBinning"] == "5"
    assert raw.params["corrFlux"] == "false"
    assert "bcdMode" not in raw.params
    assert (cal.index, cal.version) == (2, (1, 7, 9))  # REC2, not REC3: the lowest i wins
    assert cal.params == {"tfKeep": "0", "tfInterp": "0", "cumulBlock": "true", "bcdMode": "5"}
    assert matisse.raw_estimates(header(NU_HOR)) == raw
    # the generic reader, as GRAVITY's gravity.pipeline_params uses it
    assert pipeline_params(header(NU_HOR), recipe=1)["replaceTel"] == "3"
    assert pipeline_params(header(NU_HOR), recipe=2)["cumulBlock"] == "true"


def test_parse_version():
    assert matisse.parse_version("matisse/1.7.9") == (1, 7, 9)
    assert matisse.parse_version("matisse/2.0.0") == (2, 0, 0)
    for text in ("", None, "matisse", "matisse/2.0"):
        assert matisse.parse_version(text) is None


@pytest.mark.parametrize(
    ("version", "n_expected"),
    [("1.5.1", "unknown"), ("1.5.2", "flipped"), ("1.7.9", "flipped"),
     ("1.9.99", "flipped"), ("2.0.0", "standard"), ("2.5.0", "standard")],
)  # fmt: skip
def test_phase_convention_by_version(version, n_expected):
    assert matisse.phase_convention(n_band(version)) == n_expected
    lm = header(RAW, **{"ESO PRO REC1 PIPE ID": f"matisse/{version}"})
    assert matisse.phase_convention(lm) == "standard"


def test_unknown_version_in_n_band():
    assert matisse.phase_convention(n_band(None)) == "unknown"
    assert matisse.phase_convention(header(RAW, drop=("ESO PRO REC1 PIPE ID",))) == "standard"
    # without ESO DET NAME the wavelengths decide
    no_det = header(RAW, drop=("ESO DET NAME", "ESO PRO REC1 PIPE ID"))
    assert matisse.phase_convention(no_det, wave=np.linspace(8e-6, 13e-6, 5)) == "unknown"
    assert matisse.phase_convention(no_det, wave=np.linspace(3e-6, 4e-6, 5)) == "standard"
    assert matisse.phase_convention(no_det) == "unknown"


def test_negate_does_not_touch_unknown():
    old = n_band("1.5.1")
    with pytest.raises(matisse.UnverifiedConvention, match="n-phase-unknown"):
        matisse.screen(exposure(old.copy()))
    with pytest.raises(matisse.UnverifiedConvention, match="unknown"):
        matisse.screen(exposure(old.copy()), n_phase="negate")
    e = matisse.screen(exposure(old.copy()), n_phase_unknown="keep")
    np.testing.assert_array_equal(e.t3["MATISSE"].data["T3PHI"][0], np.arange(8.0))
    e = matisse.screen(exposure(old.copy()), n_phase_unknown="negate")
    np.testing.assert_array_equal(e.t3["MATISSE"].data["T3PHI"][0], -np.arange(8.0))
    np.testing.assert_array_equal(e.vis["MATISSE"].data["VISPHI"], -7.0)
    assert "N-band phases negated (mat_raw_estimates matisse/1.5.1)" in str(e.header["HISTORY"])
    with pytest.raises(ValueError, match="already negated"):  # never twice
        matisse.screen(exposure(e.header), n_phase_unknown="negate")


def test_flipped_n_phases_refused_then_negated():
    flipped = n_band("1.7.9")
    with pytest.raises(matisse.UnverifiedConvention, match="'flipped'.*--n-phase negate"):
        matisse.screen(exposure(flipped.copy()))
    with pytest.raises(matisse.UnverifiedConvention):  # the unknown option does not apply
        matisse.screen(exposure(flipped.copy()), n_phase_unknown="negate")
    e = matisse.screen(exposure(flipped.copy()), n_phase="negate")
    np.testing.assert_array_equal(e.t3["MATISSE"].data["T3PHI"][0], -np.arange(8.0))
    e = matisse.screen(exposure(flipped.copy()), n_phase="keep")
    np.testing.assert_array_equal(e.t3["MATISSE"].data["T3PHI"][0], np.arange(8.0))
    # standard N and all LM data are never touched
    for h in (n_band("2.0.0"), header(RAW)):
        e = matisse.screen(exposure(h), n_phase="negate", n_phase_unknown="negate")
        np.testing.assert_array_equal(e.t3["MATISSE"].data["T3PHI"][0], np.arange(8.0))


def test_corrflux_products_refused_for_vis2():
    with pytest.raises(ValueError, match="corrFlux is 'true'"):
        matisse.screen(exposure(header(RAW, **{"ESO PRO REC1 PARAM13 VALUE": "true"})))
    with pytest.raises(ValueError, match="corrFlux is None"):  # not recorded
        matisse.screen(exposure(header(RAW, drop=("ESO PRO REC1 PARAM13 NAME",))))
    matisse.screen(exposure(header(RAW)))  # corrFlux=false is read


# ---------------------------------------------------------------------------
# Merges
# ---------------------------------------------------------------------------


def _tables(e):
    return [t.data for t in e.v2.values()]


def test_merge_from_int_time():
    e = exposure(header(NU_HOR), int_time=NU_HOR_INT_TIME)
    evidence = matisse.merge_evidence(e.header, _tables(e))
    assert any(x.startswith("INT_TIME 120.46 s > 1.5 x") for x in evidence)
    assert matisse.is_merged(e.header, _tables(e))
    # the INT_TIME rule alone, on an otherwise clean header
    e = exposure(header(RAW), int_time=NU_HOR_INT_TIME)
    assert matisse.is_merged(e.header, _tables(e))
    with pytest.raises(ValueError, match="merged product.*INT_TIME"):
        matisse.screen(e)


def test_merge_from_repeated_baselines():
    mjds = [60290.04] * 6 + [60290.05] * 6
    e = exposure(header(RAW), mjds=mjds)
    evidence = matisse.merge_evidence(e.header, _tables(e))
    assert evidence and evidence[0].startswith("one baseline at several MJDs")
    # Evidence only, until Stage 2 counts MJDs per baseline in a RAW_INT product.
    assert not matisse.MULTI_MJD_IS_FATAL
    assert not matisse.is_merged(e.header, _tables(e))
    with pytest.warns(UserWarning, match="possible merge"):
        matisse.screen(e)


def test_merge_from_consortium_names():
    h = header(RAW, PIPEFILE="nu_Hor_IR-LM_LOW_noChop_cal_oifits_0.fits")
    assert matisse.is_merged(h)
    assert matisse.is_merged(header(RAW), path="/data/x_IR-LM_LOW_Chop_cal_oifits_0.fits")
    assert not matisse.is_merged(header(RAW), path="/data/chopped_runs/TARGET_RAW_INT_0001.fits")
    h = header(RAW)
    h.add_history("produced by mat_mergeOifits")
    assert matisse.is_merged(h)


def test_single_exposure_not_merged():
    e = exposure(header(RAW))
    assert matisse.merge_evidence(e.header, _tables(e), e.path) == []
    assert not matisse.is_merged(e.header, _tables(e), e.path)


# ---------------------------------------------------------------------------
# Compatibility, groups and dispatch
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("cards", "where", "what"),
    [
        ({"ESO INS BCD1 NAME": "IN", "ESO INS BCD2 NAME": "IN", "ESO CFG BCD MODE": "IN-IN"},
         0, "bcd_keys"),
        ({"ESO INS DIL NAME": "MED"}, 0, "resolution"),
        ({"ESO DET NAME": "MATISSE-N"}, 0, "band"),
        ({"ESO ISS CONF T1NAME": "UT1"}, 0, "telescopes"),
        ({"ESO INS SFL NAME": "1.00"}, 0, "spatial_filter"),
        ({"ESO DET SEQ1 DIT": 0.111}, 0, "dit"),
        ({"ESO PRO REC1 PARAM16 VALUE": "7"}, 1, "spectralBinning"),
        ({"ESO PRO REC1 PARAM10 VALUE": "3"}, 1, "replaceTel"),
        ({"ESO PRO REC1 PARAM16 NAME": "spectralAverage"}, 1, "spectralBinning"),
        ({"ESO PRO REC1 PARAM15 NAME": "cpAlgo"}, 1, "coherentAlgo"),
    ],
)  # fmt: skip
def test_compatibility_refuses_other_setup(cards, where, what):
    """The plan's _other_bcd, _other_resolution, _other_band, _other_spectral_binning, and more."""
    lists = matisse.compatibility(header(RAW), keyed(**cards))
    assert any(m.split(":")[0].endswith(what) for m in lists[where]), lists
    assert lists[1 - where] == []


def test_compatibility_refuses_other_chopping():
    sci, cal = header(RAW), keyed(**{"ESO ISS CHOP ST": "F"})
    for observable in ("vis2", "t3phi"):
        mode, _, _ = matisse.compatibility(sci, cal, observable)
        assert any(m.startswith("chopped") for m in mode)
    assert matisse.compatibility(sci, cal, "vis2", pool_chopping_t3=True)[0]
    assert matisse.compatibility(sci, cal, "t3phi", pool_chopping_t3=True)[0] == []
    with pytest.raises(ValueError, match="observable"):
        matisse.compatibility(sci, cal, "t3")


def test_group_key_dit_and_chopping_split_groups():
    base = exposure(header(RAW))
    other_dit = exposure(keyed(**{"ESO DET SEQ1 DIT": 0.111}))
    no_chop = exposure(keyed(**{"ESO ISS CHOP ST": "F"}))
    for obs in ("vis2", "t3phi"):
        k = matisse.group_key(base, obs)
        assert k == ("2023-12-11", "LM", "LOW", 0.0751997, ("OUT", "OUT"), True)
        assert matisse.group_key(other_dit, obs) != k
        assert matisse.group_key(no_chop, obs) != k
    pooled = {"pool_chopping_t3": True}
    assert matisse.group_key(no_chop, "t3phi", **pooled) == matisse.group_key(
        base, "t3phi", **pooled
    )
    assert matisse.group_key(no_chop, "vis2", **pooled) != matisse.group_key(base, "vis2", **pooled)
    assert matisse.group_key(other_dit, "t3phi", **pooled) != matisse.group_key(
        base, "t3phi", **pooled
    )


def test_compatibility_warns_on_pipeline_version():
    mode, reduction, minor = matisse.compatibility(
        header(RAW), keyed(**{"ESO PRO REC1 PIPE ID": "matisse/2.5.0"})
    )
    assert (mode, reduction) == ([], [])
    assert minor == [
        "mat_raw_estimates PIPE ID: science 'matisse/2.0.0', calibrator 'matisse/2.5.0'"
    ]
    # mat_raw_estimates in another slot than REC1: its own PIPE ID is used
    moved = fits.Header()
    for k, v in header(RAW).items():
        moved[f"HIERARCH {k.replace('REC1', 'REC2')}" if k.startswith("ESO ") else k] = v
    moved["HIERARCH ESO PRO REC1 ID"] = "mat_other"
    moved["HIERARCH ESO PRO REC1 PIPE ID"] = "matisse/9.9.9"
    assert matisse.compatibility(header(RAW), moved) == ([], [], [])


def test_n_band_mixed_phase_convention_refused():
    mode, _, minor = matisse.compatibility(n_band("1.7.9"), n_band("2.0.0"))
    assert any("phase convention" in m and "'flipped'" in m for m in mode)
    assert minor  # and the version differs
    assert matisse.compatibility(n_band("2.0.0"), n_band("2.5.0"))[0] == []


def _fake(h, path):
    return types.SimpleNamespace(header=h, path=path, wave={"X": np.ones(3)})


def test_check_compatible_dispatches_to_matisse():
    sci = _fake(header(RAW), "sci.fits")
    check_compatible(sci, [_fake(header(RAW), "cal.fits")], "X")
    with pytest.raises(ValueError, match="cal.fits: bcd_keys"):
        check_compatible(sci, [_fake(_bcd("IN", "IN"), "cal.fits")], "X")
    with pytest.raises(ValueError, match="spectralBinning.*--allow-mismatch"):
        check_compatible(sci, [_fake(keyed(**{"ESO PRO REC1 PARAM16 VALUE": "1"}), "c")], "X")
    gravity = header({"INSTRUME": "GRAVITY", "DATE-OBS": "2023-12-12T00:00:00"})
    with pytest.raises(ValueError, match="INSTRUME 'GRAVITY', science 'MATISSE'"):
        check_compatible(sci, [_fake(gravity, "grav.fits")], "X")


def test_check_compatible_refuses_other_dit():
    """The CLI groups by night only until Stage 4, so a DIT mismatch must be fatal here."""
    sci = _fake(header(RAW), "sci.fits")
    other = _fake(keyed(**{"ESO DET SEQ1 DIT": 0.111}), "cal.fits")
    with pytest.raises(ValueError, match="cal.fits: dit"):
        check_compatible(sci, [other], "X")
    with pytest.raises(ValueError, match="cal.fits: dit"):  # strict=False only relaxes reduction
        check_compatible(sci, [other], "X", strict=False)


class _Stop(Exception):
    pass


def _refused_by_cli(monkeypatch, tmp_path, capsys, module, downstream, argv):
    """Run a CLI on one MATISSE file; return the exposures it passed on, and its output."""
    path = tmp_path / "TARGET_RAW_INT_0001.fits"
    path.write_bytes(b"")
    monkeypatch.setattr(module, "read_exposure", lambda p, insname=None: exposure(header(RAW)))
    passed = []

    def capture(exposures, *args, **kwargs):
        passed.extend(exposures)
        raise _Stop

    monkeypatch.setattr(module, downstream, capture)
    with pytest.raises(_Stop):
        module.main([str(path), *argv])
    return passed, capsys.readouterr().out


def test_calpca_refuses_matisse(monkeypatch, tmp_path, capsys):
    from virgil_vlti import calpca_cli

    passed, out = _refused_by_cli(
        monkeypatch, tmp_path, capsys, calpca_cli, "calibrator_residuals",
        ["--out", str(tmp_path / "out"), "--screen-factor", "0"],
    )  # fmt: skip
    assert passed == []
    assert "skipped TARGET_RAW_INT_0001.fits: MATISSE calpca is Stage 5" in out


def test_pairs_cli_skips_matisse(monkeypatch, tmp_path, capsys):
    from virgil_vlti import pairs_cli

    passed, out = _refused_by_cli(monkeypatch, tmp_path, capsys, pairs_cli, "pair_test", [])
    assert passed == []
    assert "skipped TARGET_RAW_INT_0001.fits: MATISSE pair test needs find_pairs" in out


def _gravity_header(instrument, **cards):
    base = {
        "ESO INS SPEC RES": "MEDIUM",
        "ESO INS POLA MODE": "COMBINED",
        "ESO DET2 SEQ1 DIT": 30.0,
        "ESO PRO REC1 PARAM1 NAME": "vis-correction-sc",
        "ESO PRO REC1 PARAM1 VALUE": "VFACTOR",
    }
    if instrument:
        base["INSTRUME"] = instrument
    return header(base, **cards)


def test_gravity_dispatch_unchanged():
    for instrument in ("GRAVITY", None):
        sci = _fake(_gravity_header(instrument), "sci.fits")
        cal = _fake(_gravity_header(instrument, **{"ESO INS POLA MODE": "SPLIT"}), "cal.fits")
        assert instrument_module(sci.header).__name__ == "virgil_vlti.gravity"
        with pytest.raises(ValueError, match="POLA MODE"):
            check_compatible(sci, [cal], "X")
        check_compatible(sci, [_fake(_gravity_header(instrument), "c2.fits")], "X")


def test_pionier_is_refused():
    sci = header({"INSTRUME": "PIONIER", "ESO PRO REC1 PIPE ID": "pionier/3.1"})
    assert instrument_of(sci) == "PIONIER"
    with pytest.raises(ValueError, match="'PIONIER' is not supported"):
        instrument_module(sci)


# ---------------------------------------------------------------------------
# virgil-vlti-inspect
# ---------------------------------------------------------------------------


def _write(path, h, int_time):
    nw = 4
    w = fits.BinTableHDU.from_columns(
        [fits.Column("EFF_WAVE", "1E", array=np.linspace(4.2e-6, 2.8e-6, nw)),
         fits.Column("EFF_BAND", "1E", array=np.full(nw, 2e-8))]
    )  # fmt: skip
    w.header.update(EXTNAME="OI_WAVELENGTH", INSNAME="MATISSE")
    v2 = fits.BinTableHDU.from_columns(
        [
            fits.Column("MJD", "1D", array=np.full(6, 60290.04)),
            fits.Column("INT_TIME", "1D", array=np.full(6, int_time)),
            fits.Column("STA_INDEX", "2I", array=NU_HOR_STA),
        ]
    )
    v2.header.update(EXTNAME="OI_VIS2", INSNAME="MATISSE")
    fits.HDUList([fits.PrimaryHDU(header=h), w, v2]).writeto(path)
    return str(path)


def test_inspect_cli_verdicts(tmp_path, capsys):
    nu_hor = _write(tmp_path / "nuhor.fits", header(NU_HOR), NU_HOR_INT_TIME)
    assert inspect_main([nu_hor]) == 1
    out = capsys.readouterr().out
    for text in ("band              LM", "resolution        LOW", "bcd               OUT-OUT",
                 "chopped           True", "pipeline          matisse/2.0.0",
                 "merged            True", "calibrated        True",
                 "verdict           refused: calibrated; merged", "replaceTel=3"):  # fmt: skip
        assert text in out
    raw = _write(tmp_path / "raw.fits", header(RAW), 61.828)
    assert inspect_main([raw]) == 0
    assert "verdict           ok" in capsys.readouterr().out
    assert inspect_main([raw, nu_hor, "--json"]) == 1
    assert '"verdict": "ok"' in capsys.readouterr().out
    flipped = _write(tmp_path / "n.fits", n_band("1.7.9"), 61.828)
    assert inspect_main([flipped]) == 1
    assert "N-band phases flipped" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# BCD facts (plan §3.1), checked without the tables
# ---------------------------------------------------------------------------


# Which telescope (labelled by its OUT-OUT input) reaches each input, per
# device state, from the swaps alone ([PM] §4.1).
def _arrivals(position):
    ip13, ip57 = position.split("-")
    arrive = {i: i for i in matisse.INPUTS}
    if ip13 == "IN":
        arrive[1], arrive[3] = 3, 1
    if ip57 == "IN":
        arrive[5], arrive[7] = 7, 5
    return arrive


def _parity(seq):
    """+1 for an even permutation of sorted(seq), −1 for an odd one."""
    inversions = sum(a > b for a, b in itertools.combinations(seq, 2))
    return -1 if inversions % 2 else 1


def test_bcd_tables_follow_from_input_swaps():
    for position in matisse.POSITIONS:
        arrive = _arrivals(position)
        for row, (i, j) in enumerate(matisse.VIS2_PEAKS):
            assert matisse.peak_baseline(position, row) == (arrive[i], arrive[j])
            assert matisse.peak_baseline(position, (i, j)) == (arrive[i], arrive[j])
        for row, peak in enumerate(matisse.T3_PEAKS):
            seen = tuple(arrive[i] for i in peak)
            assert matisse.peak_triangle(position, row) == (tuple(sorted(seen)), _parity(seen))
    # [PM] §11.9: (φ_OO + φ_OI − φ_IO − φ_II)/4 for triangle IP1-IP3-IP5, on its own peak.
    signs = {p: matisse.peak_triangle(p, (1, 3, 5)) for p in matisse.POSITIONS}
    assert {p: s for p, (_, s) in signs.items()} == {
        "OUT-OUT": 1, "OUT-IN": 1, "IN-OUT": -1, "IN-IN": -1
    }  # fmt: skip
    assert signs["OUT-OUT"][0] == signs["IN-OUT"][0] == (1, 3, 5)


# A binary of flux ratio 0.3 seen on stations K0 G2 D0 J3 (STA_INDEX 28, 19,
# 13, 24, at inputs 1, 3, 5, 7 as in ν Hor). Station positions (m, east and
# north) are illustrative; only their geometry matters.
STATION_OF_INPUT = {1: 28, 3: 19, 5: 13, 7: 24}
STATION_XY = {28: (106.0, -14.0), 19: (-20.0, 48.0), 13: (-60.0, -32.0), 24: (82.0, 70.0)}
BINARY_MAS = (0.3, 4.3)  # every |closure phase| between 19° and 161°
FLUX_RATIO = 0.3
WAVE_LM = np.linspace(3.0e-6, 4.0e-6, 5)
MAS_RAD = np.pi / 180 / 3600e3


def _vis(sta_a, sta_b):
    """Complex visibility of the binary on baseline a -> b, u = (x_b − x_a)/λ."""
    d = np.subtract(STATION_XY[sta_b], STATION_XY[sta_a])
    phase = 2 * np.pi * (d[0] * BINARY_MAS[0] + d[1] * BINARY_MAS[1]) * MAS_RAD / WAVE_LM
    return (1 + FLUX_RATIO * np.exp(-1j * phase)) / (1 + FLUX_RATIO)


def _closure(a, b, c):
    """Closure phase (deg) of stations in the order given: arg V_ab V_bc V_ac*."""
    return np.degrees(np.angle(_vis(a, b) * _vis(b, c) * np.conj(_vis(a, c))))


# A different instrumental closure phase on each peak (deg), ±40° at most.
PHI_INST = {
    (3, 5, 7): 23.0 + 0 * WAVE_LM,
    (1, 3, 5): -37.0 + 4 * np.arange(5),
    (1, 3, 7): 11.0 - 2 * np.arange(5),
    (1, 5, 7): -8.0 + 0 * WAVE_LM,
}


def _measure(position):
    """T3PHI per peak at a position, from the swaps and the model alone."""
    arrive = _arrivals(position)
    return [
        _closure(*(STATION_OF_INPUT[arrive[i]] for i in peak)) + PHI_INST[peak]
        for peak in matisse.T3_PEAKS
    ]


def _keys(position, truth):
    """The header keys (BCD1, BCD2) of a device-order position, if ``truth`` holds."""
    ip13, ip57 = position.split("-")
    return (ip13, ip57) if truth == "BCD1=IP13" else (ip57, ip13)


def _combined(truth, assumed):
    """Equal-weight BCD combination, through the headers, setup and §3.1 tables."""
    sums = {}
    for position in matisse.POSITIONS:
        b1, b2 = _keys(position, truth)
        s = matisse.setup(_bcd(b1, b2), bcd_key_order=assumed)
        for row, phi in enumerate(_measure(position)):
            triangle, sign = matisse.peak_triangle(s.bcd, row)
            stations = [STATION_OF_INPUT[i] for i in triangle]
            key = tuple(sorted(stations))
            value = _parity(stations) * sign * phi  # to sorted-station orientation
            sums[key] = sums.get(key, 0) + np.exp(1j * np.radians(value)) / 4
    return {k: np.degrees(np.angle(v)) for k, v in sums.items()}


def _wrap(x):
    return (np.asarray(x) + 180) % 360 - 180


@pytest.mark.parametrize("truth", matisse.BCD_KEY_ORDERS)
def test_closure_phase_sign_against_model(truth):
    """The §3.1 tables combine simulated BCD measurements into the model's closure phases.

    This catches **transcription errors only** (see
    ``test_closure_phase_sign_test_catches_table_errors``). The expected
    closure phases come from station positions, but the simulated
    measurements (``_measure``, via ``_arrivals``) use the same §3.1
    input-swap model the tables were derived from: if the real BCD did
    something else (reversed the peak orientation as well, say, or another
    ``T3_PEAKS`` row order), tables and simulator would be wrong together and
    this would still pass. It is not the plan's sign gate. The only
    independent anchor here is the [PM] §11.9 sign row
    (``test_bcd_tables_follow_from_input_swaps``); real validation needs the
    pipeline's combined closure phase on D1 (plan Stage 4).
    """
    expected = {
        tri: _closure(*tri) for tri in itertools.combinations(sorted(STATION_OF_INPUT.values()), 3)
    }
    for v in expected.values():  # a sign error would show
        assert np.all((np.abs(_wrap(v)) > 10) & (np.abs(_wrap(v)) < 170))
    got = _combined(truth, assumed=truth)
    assert sorted(got) == sorted(expected)
    for tri, want in expected.items():
        np.testing.assert_allclose(_wrap(got[tri] - want), 0, atol=1e-6, err_msg=str(tri))
    # The reversed key order mislabels IN-OUT and OUT-IN and fails.
    other = next(o for o in matisse.BCD_KEY_ORDERS if o != truth)
    wrong = _combined(truth, assumed=other)
    assert max(np.max(np.abs(_wrap(wrong[t] - expected[t]))) for t in expected) > 1


@pytest.mark.parametrize(
    "mutate",
    ["swap rows", "flip sign", "swap positions"],
)
def test_closure_phase_sign_test_catches_table_errors(monkeypatch, mutate):
    """A swapped row, a flipped sign or swapped positions in the §3.1 table fail the model test."""
    table = dict(matisse.PEAK_TRIANGLES)
    if mutate == "swap rows":
        rows = list(table["IN-OUT"])
        rows[1], rows[2] = rows[2], rows[1]
        table["IN-OUT"] = tuple(rows)
    elif mutate == "flip sign":
        rows = list(table["IN-IN"])
        rows[0] = (rows[0][0], -rows[0][1])
        table["IN-IN"] = tuple(rows)
    else:
        table["IN-OUT"], table["OUT-IN"] = table["OUT-IN"], table["IN-OUT"]
    monkeypatch.setattr(matisse, "PEAK_TRIANGLES", table)
    got = _combined("BCD1=IP13", assumed="BCD1=IP13")
    errors = [np.max(np.abs(_wrap(got[t] - _closure(*t)))) for t in got]
    assert max(errors) > 1


def test_unequal_positions_leave_instrumental_phase():
    """OO plus II alone (no single-device positions) leaves (φ_i(135) − φ_i(137))/2."""
    sums = {}
    for position in ("OUT-OUT", "IN-IN"):
        for row, phi in enumerate(_measure(position)):
            triangle, sign = matisse.peak_triangle(position, row)
            stations = [STATION_OF_INPUT[i] for i in triangle]
            key = tuple(sorted(stations))
            value = _parity(stations) * sign * phi
            sums[key] = sums.get(key, 0) + np.exp(1j * np.radians(value))
    errors = [
        np.max(np.abs(_wrap(np.degrees(np.angle(v)) - _closure(*k)))) for k, v in sums.items()
    ]
    assert max(errors) > 1
