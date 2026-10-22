import argparse
import dataclasses
import json
import warnings

import numpy as np
import pytest
from astropy.io import fits
from scipy.special import j1

from virgil_vlti import calibrate_exposure, classify, frame_times, read_exposure, write_calibrated
from virgil_vlti.calibrate import (
    MAS,
    bin_channels,
    calibrator_estimates,
    channel_correlation,
    diameter_modes,
    flux_ratio,
    predict,
    prediction_error,
    propagated_error,
    regularize_shape,
    shape_slope_significance,
    transfer_function_modes,
    uniform_disk_v2,
)
from virgil_vlti.cli import _calibrator, _diameter, main
from virgil_vlti.gravity import (
    anisoplanatic_v2_factor,
    anisoplanatic_v2_loss,
    telescope_diameter,
)
from virgil_vlti.products import kind_from_category
from virgil_vlti.synthetic import add_systematics, write_product

MJD0 = 60400.1
WAVE = np.linspace(2.0e-6, 2.4e-6, 40)


def gaussian_v2(fwhm_mas):
    def v2(u, v, wave):
        rho = np.hypot(u, v) / wave
        return np.exp(-((np.pi * fwhm_mas * MAS * rho) ** 2) / (4 * np.log(2))) ** 2

    return v2


def transfer(t, wave):
    """A chromatic transfer function drifting linearly in time."""
    return (0.8 - 0.3 * (t - MJD0)) * (wave / 2.2e-6) ** 0.5


def instrumental_cp(wave):
    """A constant, chromatic instrumental closure phase (degrees)."""
    return 4.0 + 3.0 * (wave - 2.2e-6) / 0.2e-6


def science_cp(u1, v1, u2, v2, wave):
    return 20.0 * np.sin(u1 / 40.0) + 0.0 * wave


def amplitude_transfer(t, wave):
    """A transfer function for |V| and T3AMP, also linear in time."""
    return (0.9 - 0.2 * (t - MJD0)) * (wave / 2.2e-6) ** 0.25


def instrumental_visphi(wave):
    return 2.0 - 5.0 * (wave - 2.2e-6) / 0.2e-6


def science_visphi(u, v, wave):
    return 10.0 * np.sin(u / 30.0) + 0.0 * wave


def uniform_disk(diameter):
    def amp(u, v, w):
        if not diameter:
            return np.ones_like(w)
        x = np.pi * diameter * MAS * np.hypot(u, v) / w
        return np.abs(2 * j1(x) / x)

    return amp


def product(tmp_path, name, t, *, cal, diameter=0.0, sobj=None, **kw):
    """A calibrator or science exposure seen through the instrument at time t.

    Writes V², T3AMP, T3PHI and OI_VIS, each through its own transfer
    function (linear in time) or instrumental phase (constant).
    """
    if cal:
        amp = uniform_disk(diameter)

        def cp(u1, v1, u2, v2, w):
            return 0.0 * w

        def visphi(u, v, w):
            return 0.0 * w
    else:

        def amp(u, v, w):
            return np.sqrt(gaussian_v2(2.0)(u, v, w))

        cp, visphi = science_cp, science_visphi

    def t3amp(u1, v1, u2, v2, w, tk):
        triple = amp(u1, v1, w) * amp(u2, v2, w) * amp(u1 + u2, v1 + v2, w)
        return triple * amplitude_transfer(tk, w)

    def vis(u, v, w, tk):
        return amp(u, v, w) * amplitude_transfer(tk, w), visphi(u, v, w) + instrumental_visphi(w)

    return write_product(
        tmp_path / f"{name}.fits",
        mjd=t,
        wave=WAVE,
        vis2=lambda u, v, w, tk: amp(u, v, w) ** 2 * transfer(tk, w),
        t3phi=lambda u1, v1, u2, v2, w, tk: cp(u1, v1, u2, v2, w) + instrumental_cp(w),
        t3amp=t3amp,
        vis=vis,
        sobj=sobj if sobj is not None else ((0.0, 0.0) if cal else (210.0, -710.0)),
        catg="DUAL_CAL_VIS" if cal else "DUAL_SCI_VIS",
        target="CALSTAR" if cal else "SCISTAR",
        **kw,
    )


@pytest.fixture
def night(tmp_path):
    """Calibrators before and after a science exposure, in dual-field layout."""
    paths = [
        product(tmp_path, "cal1", MJD0, cal=True),
        product(tmp_path, "sci", MJD0 + 0.02, cal=False),
        product(tmp_path, "cal2", MJD0 + 0.05, cal=True),
    ]
    exposures = classify([read_exposure(p) for p in paths], dual_field=True)
    return exposures, tmp_path


def _truth_v2(sci):
    t = sci.v2["GRAVITY_SC"].data
    return np.array(
        [gaussian_v2(2.0)(u, v, WAVE) for u, v in zip(t["UCOORD"], t["VCOORD"], strict=True)]
    )


def _truth_cp(sci):
    t = sci.t3["GRAVITY_SC"].data
    return np.array([science_cp(t["U1COORD"][i], 0, 0, 0, WAVE) for i in range(len(t))])


def test_fringe_tracker_is_never_read(night):
    exposures, _ = night
    assert all(e.insnames == ["GRAVITY_SC"] for e in exposures)


def test_dual_field_classification(night):
    exposures, _ = night
    assert [e.kind for e in exposures] == ["CAL", "SCI", "CAL"]


def test_category_and_name_classification(night):
    exposures, _ = night
    assert [e.kind for e in classify(exposures)] == ["CAL", "SCI", "CAL"]
    assert [e.kind for e in classify(exposures, calibrators=["cal star"])] == ["CAL", "SCI", "CAL"]


def test_point_source_calibration_is_exact(night):
    exposures, _ = night
    sci = exposures[1]
    cals = [e for e in exposures if e.kind == "CAL"]
    out = calibrate_exposure(sci, cals, "GRAVITY_SC")
    # A transfer function linear in time is interpolated exactly.
    np.testing.assert_allclose(out["vis2"].value, _truth_v2(sci), rtol=1e-6)
    # A constant instrumental closure phase is removed exactly.
    np.testing.assert_allclose(out["t3phi"].value, _truth_cp(sci), atol=1e-6)
    assert not out["vis2"].flag.any() and not out["t3phi"].flag.any()


def test_resolved_calibrator_needs_its_diameter(tmp_path):
    paths = [
        product(tmp_path, "cal1", MJD0, cal=True, diameter=1.2),
        product(tmp_path, "sci", MJD0 + 0.02, cal=False),
        product(tmp_path, "cal2", MJD0 + 0.05, cal=True, diameter=1.2),
    ]
    exposures = classify([read_exposure(p) for p in paths], dual_field=True)
    sci, cals = exposures[1], [exposures[0], exposures[2]]
    truth = _truth_v2(sci)
    right = calibrate_exposure(sci, cals, "GRAVITY_SC", diameters={"CALSTAR": 1.2})
    wrong = calibrate_exposure(sci, cals, "GRAVITY_SC")
    np.testing.assert_allclose(right["vis2"].value, truth, rtol=1e-6)
    assert np.all(wrong["vis2"].value > truth * 1.001)  # a resolved calibrator inflates V²


def test_window_and_binning(night):
    exposures, _ = night
    out = calibrate_exposure(
        exposures[1], [exposures[0], exposures[2]], "GRAVITY_SC", window_um=(2.1, 2.4), nbin=4
    )
    keep = (WAVE > 2.1e-6) & (WAVE < 2.4e-6)
    n = keep.sum() // 4
    assert out["vis2"].value.shape == (18, n)
    np.testing.assert_allclose(out.wave, WAVE[keep][: 4 * n].reshape(n, 4).mean(1), rtol=1e-6)
    expect = _truth_v2(exposures[1])[:, keep][:, : 4 * n].reshape(18, n, 4).mean(-1)
    np.testing.assert_allclose(out["vis2"].value, expect, rtol=1e-6)


def test_split_polarisation_calibrated_separately(tmp_path):
    ins = ("GRAVITY_SC_P1", "GRAVITY_SC_P2")
    paths = [
        product(tmp_path, "cal1", MJD0, cal=True, insnames=ins),
        product(tmp_path, "sci", MJD0 + 0.02, cal=False, insnames=ins),
        product(tmp_path, "cal2", MJD0 + 0.05, cal=True, insnames=ins),
    ]
    exposures = classify([read_exposure(p) for p in paths], dual_field=True)
    assert exposures[1].insnames == list(ins)
    for name in ins:
        out = calibrate_exposure(exposures[1], [exposures[0], exposures[2]], name)
        np.testing.assert_allclose(
            out["vis2"].value, _truth_v2_named(exposures[1], name), rtol=1e-6
        )


def test_split_polarisation_round_trip(tmp_path):
    """P1 and P2 are written with unique EXTVERs matching their data tables, and the
    extra per-channel columns of a pipeline table are dropped, not mis-sized."""
    ins = ("GRAVITY_SC_P1", "GRAVITY_SC_P2")
    paths = [
        product(tmp_path, "cal1", MJD0, cal=True, insnames=ins),
        product(tmp_path, "sci", MJD0 + 0.02, cal=False, insnames=ins),
        product(tmp_path, "cal2", MJD0 + 0.05, cal=True, insnames=ins),
    ]
    exposures = classify([read_exposure(p) for p in paths], dual_field=True)
    sci = exposures[1]
    for name in ins:  # a pipeline-style extra column with one value per channel
        for tables in (sci.v2, sci.t3):
            t = tables[name]
            extra = fits.Column("EXTRA", f"{WAVE.size}D", array=np.zeros((len(t.data), WAVE.size)))
            tables[name] = fits.BinTableHDU(
                data=fits.BinTableHDU.from_columns(t.columns + extra).data, header=t.header
            )
    results = [calibrate_exposure(sci, [exposures[0], exposures[2]], name, nbin=4) for name in ins]
    path = write_calibrated(sci, results, str(tmp_path / "out.fits"))
    with fits.open(path) as h:
        for extname in ("OI_WAVELENGTH", "OI_VIS2", "OI_T3"):
            tables = [x for x in h if x.name == extname and x.header["INSNAME"] in ins]
            assert [x.header["INSNAME"] for x in tables] == list(ins)
            assert [x.ver for x in tables] == [1, 2]
        for x in h:
            if x.name in ("OI_VIS2", "OI_T3"):
                assert "EXTRA" not in x.columns.names
                assert x.header["NWAVE"] == WAVE.size // 4
                assert np.shape(x.data["FLAG"])[1] == x.header["NWAVE"]
        assert h["OI_WAVELENGTH", 2].header["INSNAME"] == ins[1]


def test_gaussian_phase_floor(night):
    exposures, _ = night
    est = calibrator_estimates([exposures[0], exposures[2]], "GRAVITY_SC", "t3phi")
    e = next(iter(est.values()))
    t = MJD0 + 0.02
    dt = (t - e.times) * 86400
    tiny = dataclasses.replace(e, errors=np.full_like(e.errors, 1e-6))
    floored = np.exp(-2 * (dt / 3600) ** 2) / 0.2**2  # both errors floored at 0.2 degrees
    w = floored[:, None] * np.ones(e.values.shape)
    expect = (w * e.values).sum(0)
    np.testing.assert_allclose(predict(t, tiny, "gaussian"), expect / np.abs(expect), rtol=1e-6)


def test_many_long_calibrator_names_fit_the_header(night, tmp_path):
    exposures, _ = night
    sci, cals = exposures[1], [exposures[0], exposures[2]]
    out = calibrate_exposure(sci, cals, "GRAVITY_SC")
    prov = {"NCAL": 4}
    for k in range(1, 5):
        prov |= {f"CAL{k} NAME": f"HD {123450 + k} (NAME LONG CALIBRATOR)", f"CAL{k} DIAM": 0.0123,
                 f"CAL{k} DIAMERR": 0.0012}  # fmt: skip
    path = write_calibrated(sci, [out], str(tmp_path / "out.fits"), provenance=prov)
    with fits.open(path) as h:
        assert h[0].header["HIERARCH VIRGIL CAL4 NAME"] == "HD 123454 (NAME LONG CALIBRATOR)"
        assert h[0].header["HIERARCH VIRGIL CAL4 DIAMERR"] == 0.0012


@pytest.mark.parametrize("catg", ["DUAL_SCI_TF", "DUAL_CAL_TF"])
def test_transfer_function_products_are_refused(night, catg):
    exposures, _ = night
    exposures[1].header["HIERARCH ESO PRO CATG"] = catg
    with pytest.raises(ValueError, match="transfer-function product"):
        classify(exposures, dual_field=True)
    with pytest.raises(ValueError, match="transfer-function product"):
        classify(exposures, calibrators=["CALSTAR"])


def _truth_v2_named(sci, name):
    t = sci.v2[name].data
    return np.array(
        [gaussian_v2(2.0)(u, v, WAVE) for u, v in zip(t["UCOORD"], t["VCOORD"], strict=True)]
    )


def test_frame_times_are_unified(night):
    exposures, _ = night
    sci = exposures[1]
    v2, t3 = sci.v2["GRAVITY_SC"], sci.t3["GRAVITY_SC"]
    assert np.ptp(v2.data["MJD"][:6]) > 1e-5  # the baselines of a frame differ
    mjd_v2, mjd_t3 = frame_times(v2, t3)
    assert np.ptp(mjd_v2[:6]) == 0 and np.all(mjd_t3[:4] == mjd_v2[0])
    np.testing.assert_allclose(mjd_v2[0], v2.data["MJD"][:6].mean())


def test_prediction_errors(night):
    exposures, _ = night
    cals = [exposures[0], exposures[2]]
    est = calibrator_estimates(cals, "GRAVITY_SC", "vis2")[(28, 23)]
    # Two estimates: leave-one-out predicts each from the other (a constant).
    loo = prediction_error(est, "loo")
    np.testing.assert_allclose(loo, np.abs(est.values[1] - est.values[0]), rtol=1e-6)
    scatter = prediction_error(est, "scatter")
    np.testing.assert_allclose(scatter, np.std(est.values, 0, ddof=1), rtol=1e-6)
    one = calibrator_estimates(cals[:1], "GRAVITY_SC", "vis2")[(28, 23)]
    assert np.all(prediction_error(one, "loo") == 0)


def test_gaussian_interpolation_weights(night):
    exposures, _ = night
    est = calibrator_estimates([exposures[0], exposures[2]], "GRAVITY_SC", "vis2")[(28, 23)]
    t = MJD0 + 0.02
    sigma = np.fmax(np.nanmedian(est.errors, 1), 0.02 * np.abs(np.nanmedian(est.values, 1)))
    w = np.exp(-2 * ((t - est.times) * 86400 / 3600) ** 2) / sigma**2
    expect = (w[:, None] * est.values).sum(0) / w.sum()
    np.testing.assert_allclose(predict(t, est, "gaussian"), expect, rtol=1e-6)
    # An implausibly small error is floored at 2% of the transfer function, so it
    # cannot dominate (as in the pipeline: gravi_tf.c).
    tiny = dataclasses.replace(est, errors=est.errors * np.array([[1e-6], [1.0]]))
    sigma = np.fmax(np.nanmedian(tiny.errors, 1), 0.02 * np.abs(np.nanmedian(est.values, 1)))
    assert sigma[0] > np.nanmedian(tiny.errors, 1)[0]
    w = np.exp(-2 * ((t - est.times) * 86400 / 3600) ** 2) / sigma**2
    np.testing.assert_allclose(
        predict(t, tiny, "gaussian"), (w[:, None] * est.values).sum(0) / w.sum(), rtol=1e-6
    )
    # Closure-phase estimates stay unit phasors.
    cp = calibrator_estimates([exposures[0], exposures[2]], "GRAVITY_SC", "t3phi")
    for e in cp.values():
        np.testing.assert_allclose(np.abs(predict(t, e, "gaussian")), 1.0)


def test_missing_stations_are_flagged(tmp_path, night):
    exposures, _ = night
    sci, cals = exposures[1], [exposures[0], exposures[2]]
    for cal in cals:  # drop one baseline from the calibrators
        d = cal.v2["GRAVITY_SC"].data
        d["STA_INDEX"][d["STA_INDEX"][:, 0] == 18] = (99, 98)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        out = calibrate_exposure(sci, cals, "GRAVITY_SC")
    assert any("No calibrator data" in str(x.message) for x in w)
    assert out["vis2"].flag.any() and not out["vis2"].flag.all()


def test_cli_writes_files_virgil_can_read(night):
    _, folder = night
    out = folder / "calibrated"
    main([str(folder), "--out", str(out), "--dual-field", "--bin", "2"])
    files = sorted(out.glob("*/*.fits"))
    assert len(files) == 1
    assert (out / "index.csv").exists() and (out / "flux_ratio.csv").exists()
    with fits.open(files[0]) as h:
        head = h[0].header
        assert head["HIERARCH VIRGIL CALIBRATED"] is True
        assert "VISCAL" not in head and "ESO PRO CATG" not in head
        assert not [k for k in head if k.startswith("ESO PRO REC")]
        assert head["HIERARCH VIRGIL INPUT CATG"] == "DUAL_SCI_VIS"
        assert head["HIERARCH VIRGIL INPUT PIPE ID"] == "not recorded"
        assert head["HIERARCH VIRGIL CAL1 NAME"] == "CALSTAR"
        assert head["HIERARCH VIRGIL NCAL"] == 1
        assert head["HIERARCH VIRGIL TF ERROR"] == "loo"
        assert head["HIERARCH VIRGIL VERSION"]
        assert {x.header.get("INSNAME") for x in h[3:]} == {"GRAVITY_SC"}
    oifits = pytest.importorskip("virgil.oifits")
    data = oifits.read_oifits(str(files[0]))
    assert np.size(data["vis"]) == 18 * 20


def _with_cal_header(tmp_path, header):
    paths = [
        product(tmp_path, "cal1", MJD0, cal=True, header=header),
        product(tmp_path, "sci", MJD0 + 0.02, cal=False),
        product(tmp_path, "cal2", MJD0 + 0.05, cal=True),
    ]
    exposures = classify([read_exposure(p) for p in paths], dual_field=True)
    return exposures[1], [exposures[0], exposures[2]]


@pytest.mark.parametrize("header", [{"ESO INS SPEC RES": "HIGH"}, {"ESO INS POLA MODE": "SPLIT"}])
def test_calibrators_in_another_mode_are_always_refused(tmp_path, header):
    sci, cals = _with_cal_header(tmp_path, header)
    for strict in (True, False):
        with pytest.raises(ValueError, match="incompatible"):
            calibrate_exposure(sci, cals, "GRAVITY_SC", strict=strict)


def test_calibrators_on_another_wavelength_grid_are_refused(tmp_path):
    paths = [
        product(tmp_path, "cal1", MJD0, cal=True),
        product(tmp_path, "sci", MJD0 + 0.02, cal=False),
    ]
    exposures = classify([read_exposure(p) for p in paths], dual_field=True)
    exposures[0].wave["GRAVITY_SC"] = exposures[0].wave["GRAVITY_SC"][:-1]
    with pytest.raises(ValueError, match="wavelength grid"):
        calibrate_exposure(exposures[1], exposures[:1], "GRAVITY_SC", strict=False)


def test_other_pipeline_settings_refused_unless_not_strict(tmp_path):
    sci, cals = _with_cal_header(tmp_path, {"ESO PRO REC1 PARAM1 VALUE": "NONE"})
    with pytest.raises(ValueError, match="vis-correction-sc"):
        calibrate_exposure(sci, cals, "GRAVITY_SC")
    with pytest.warns(UserWarning, match="vis-correction-sc"):
        calibrate_exposure(sci, cals, "GRAVITY_SC", strict=False)


def test_different_dit_is_only_warned_about(tmp_path):
    sci, cals = _with_cal_header(tmp_path, {"ESO DET2 SEQ1 DIT": 10.0})
    with pytest.warns(UserWarning, match="DIT"):
        calibrate_exposure(sci, cals, "GRAVITY_SC")


def _truth_amp(u, v):
    return np.sqrt(gaussian_v2(2.0)(u, v, WAVE))


@pytest.mark.parametrize("diameter", [0.0, 1.2])
def test_t3amp_and_oi_vis_are_calibrated(tmp_path, diameter):
    paths = [
        product(tmp_path, "cal1", MJD0, cal=True, diameter=diameter),
        product(tmp_path, "sci", MJD0 + 0.02, cal=False),
        product(tmp_path, "cal2", MJD0 + 0.05, cal=True, diameter=diameter),
    ]
    exposures = classify([read_exposure(p) for p in paths], dual_field=True)
    sci, cals = exposures[1], [exposures[0], exposures[2]]
    out = calibrate_exposure(sci, cals, "GRAVITY_SC", vis=True, diameters={"CALSTAR": diameter})
    t = sci.t3["GRAVITY_SC"].data
    u1, v1, u2, v2 = t["U1COORD"], t["V1COORD"], t["U2COORD"], t["V2COORD"]
    triple = np.array(
        [
            _truth_amp(a, b) * _truth_amp(c, d) * _truth_amp(a + c, b + d)
            for a, b, c, d in zip(u1, v1, u2, v2, strict=True)
        ]
    )
    np.testing.assert_allclose(out["t3amp"].value, triple, rtol=1e-6)
    d = sci.vis["GRAVITY_SC"].data
    pairs = zip(d["UCOORD"], d["VCOORD"], strict=True)
    np.testing.assert_allclose(
        out["visamp"].value, np.array([_truth_amp(a, b) for a, b in pairs]), rtol=1e-6
    )
    expect = np.array([science_visphi(a, 0, WAVE) for a in d["UCOORD"]])
    np.testing.assert_allclose(out["visphi"].value, expect, atol=1e-6)


def test_oi_vis_is_optional_and_needs_visibility_amplitudes(night):
    exposures, _ = night
    sci, cals = exposures[1], [exposures[0], exposures[2]]
    assert set(calibrate_exposure(sci, cals, "GRAVITY_SC").obs) == {"vis2", "t3amp", "t3phi"}
    sci.vis["GRAVITY_SC"].header["AMPTYP"] = "correlated flux"
    with pytest.raises(ValueError, match="AMPTYP"):
        calibrate_exposure(sci, cals, "GRAVITY_SC", vis=True)


def test_written_oi_vis_has_standard_columns_and_phase_flags(night):
    _, folder = night
    out = folder / "calibrated"
    main([str(folder), "--out", str(out), "--dual-field", "--vis"])
    (path,) = out.glob("*/*.fits")
    with fits.open(path) as h:
        vis = h["OI_VIS"]
        assert "VISDATA" not in vis.columns.names
        assert {"VISAMP", "VISPHI", "FLAG"} <= set(vis.columns.names)
        assert vis.header["PHITYP"] == "absolute"
        assert np.isfinite(h["OI_T3"].data["T3AMP"]).all()


def test_a_bad_amplitude_never_costs_a_good_phase(night, tmp_path):
    exposures, _ = night
    sci, cals = exposures[1], [exposures[0], exposures[2]]
    for cal in cals:  # T3AMP unusable in channel 5 of every calibrator
        cal.t3["GRAVITY_SC"].data["T3AMP"][:, 5] = np.nan
    out = calibrate_exposure(sci, cals, "GRAVITY_SC")
    assert out["t3amp"].flag[:, 5].all() and not out["t3phi"].flag[:, 5].any()
    path = write_calibrated(sci, [out], str(tmp_path / "out.fits"))
    with fits.open(path) as h:
        t3 = h["OI_T3"].data
        assert not t3["FLAG"][:, 5].any()
        assert np.isnan(t3["T3AMP"][:, 5]).all() and np.isfinite(t3["T3PHI"][:, 5]).all()


@pytest.mark.parametrize(
    ("catg", "kind"),
    [
        ("DUAL_SCI_VIS", "SCI"),
        ("SINGLE_CAL_VIS", "CAL"),
        ("TARGET_RAW_INT", "SCI"),
        ("CALIB_RAW_INT", "CAL"),
        ("SOMETHING_ELSE", None),
        ("TARGET_OTHER", None),
        ("CALIB_OTHER", None),
    ],
)
def test_kind_from_category(catg, kind):
    assert kind_from_category(fits.Header({"HIERARCH ESO PRO CATG": catg})) == kind


@pytest.mark.parametrize(
    "card",
    [
        {"HIERARCH ESO PRO CATG": "TARGET_CAL_INT"},
        {"HIERARCH ESO PRO CATG": "CALIB_CAL_INT"},
        {"HIERARCH ESO PRO CATG": "DUAL_SCI_VIS_CALIBRATED"},
        {"VISCAL": "CALIBRATED"},
        {"HIERARCH VIRGIL CALIBRATED": True},
    ],
)
def test_calibrated_products_are_refused(night, card):
    exposures, _ = night
    exposures[1].header.update(card)
    with pytest.raises(ValueError, match="already calibrated"):
        classify(exposures, dual_field=True)


def test_diameter_modes_match_a_finite_difference(tmp_path):
    paths = [
        product(tmp_path, "cal1", MJD0, cal=True, diameter=1.0),
        product(tmp_path, "sci", MJD0 + 0.02, cal=False),
        product(tmp_path, "cal2", MJD0 + 0.05, cal=True, diameter=1.0),
    ]
    exposures = classify([read_exposure(p) for p in paths], dual_field=True)
    sci, cals = exposures[1], [exposures[0], exposures[2]]
    diam, sigma = {"CALSTAR": 1.0}, {"CALSTAR": 0.1}
    modes = diameter_modes(sci, cals, "GRAVITY_SC", diam, sigma)
    (name,) = modes
    assert name == f"CALSTAR@{sci.night}"
    assert set(modes[name]) == {"vis2", "t3amp"}
    # A larger calibrator would have had a smaller model V², so a larger
    # transfer function and a smaller calibrated V²: the mode is negative.
    base = calibrate_exposure(sci, cals, "GRAVITY_SC", diameters=diam)["vis2"].value
    bigger = calibrate_exposure(sci, cals, "GRAVITY_SC", diameters={"CALSTAR": 1.001})
    slope = (bigger["vis2"].value - base) / 0.001
    np.testing.assert_allclose(modes[name]["vis2"], 0.1 * slope, rtol=2e-2)
    assert np.all(modes[name]["vis2"] < 0)


def test_unresolved_calibrator_still_has_a_mode(night):
    exposures, _ = night
    sci, cals = exposures[1], [exposures[0], exposures[2]]
    modes = diameter_modes(sci, cals, "GRAVITY_SC", {}, {None: 0.5})
    (mode,) = modes.values()
    assert np.all(mode["vis2"] <= 0) and np.any(mode["vis2"] < 0)


def test_cli_writes_diameter_modes(night):
    _, folder = night
    out = folder / "calibrated"
    main([str(folder), "--out", str(out), "--dual-field", "--diameter", "0.3±0.05"])
    (path,) = out.glob("*/*.fits")
    modes = np.load(str(path)[: -len(".fits")] + ".modes.npz")
    files = [k for k in modes.files if not k.startswith("__")]
    key = next(k for k in files if k.startswith("GRAVITY_SC/vis2/"))
    with fits.open(path) as h:
        assert modes[key].shape == h["OI_VIS2"].data["VIS2DATA"].shape
    # The model-relative form of the same mode, under its own key.
    rel = key.replace("/vis2/", "/vis2/rel-")
    absolute = [k for k in files if "/rel-" not in k]
    assert rel in files
    assert sorted(files) == sorted(
        absolute + [k.rsplit("/", 1)[0] + "/rel-" + k.rsplit("/", 1)[1] for k in absolute]
    )
    assert np.all(modes[rel] <= 0) and np.any(modes[rel] < 0)
    index = (out / "index.csv").read_text()
    assert "max_diameter_mode_sigma" in index


@pytest.mark.parametrize(
    ("card", "match"),
    [
        ({"AMPTYP": "correlated flux"}, "AMPTYP"),
        ({"PHITYP": "differential"}, "PHITYP"),
    ],
)
def test_every_calibrator_oi_vis_is_checked(night, card, match):
    exposures, _ = night
    sci, cals = exposures[1], [exposures[0], exposures[2]]
    cals[1].vis["GRAVITY_SC"].header.update(card)
    with pytest.raises(ValueError, match=match):
        calibrate_exposure(sci, cals, "GRAVITY_SC", vis=True)
    del cals[1].vis["GRAVITY_SC"]
    with pytest.raises(ValueError, match="no OI_VIS"):
        calibrate_exposure(sci, cals, "GRAVITY_SC", vis=True)


def test_resolved_calibrator_errors_are_in_transfer_function_units(tmp_path):
    paths = [product(tmp_path, "cal1", MJD0, cal=True, diameter=1.2)]
    (cal,) = classify([read_exposure(p) for p in paths], dual_field=True)
    est = calibrator_estimates([cal], "GRAVITY_SC", "vis2", {"CALSTAR": 1.2})
    d = cal.v2["GRAVITY_SC"].data
    model = uniform_disk_v2(cal.v2["GRAVITY_SC"], WAVE, 1.2)
    rows = [i for i, s in enumerate(map(tuple, d["STA_INDEX"])) if s == (28, 23)]
    frames = d["VIS2ERR"][rows] / model[rows]  # the error of the mean over frames
    expect = np.sqrt(np.sum(frames**2, 0)) / len(rows)
    np.testing.assert_allclose(est[(28, 23)].errors[0], expect, rtol=1e-6)


def test_anisoplanatism_reproduces_gravity_wide():
    # GRAVITY+ Collaboration 2022, §4.2: HD 48543B, 7.90 arcsec from the
    # fringe tracker, seeing 0.41 and theta0 1.7 arcsec at 500 nm, ATs: a
    # visibility loss of about 0.2-0.4 across K, theta0 = 10.0 arcsec in K.
    factor = anisoplanatic_v2_factor(np.array([2.0e-6, 2.2e-6, 2.4e-6]), 7.90, 0.41, 1.7, 1.8)
    loss = 1 - np.sqrt(factor)
    assert 0.2 < loss[-1] < loss[1] < loss[0] < 0.4
    # The exponent goes as 1/λ².
    ratio = np.log(factor[0]) / np.log(factor[2])
    np.testing.assert_allclose(ratio, (2.4 / 2.0) ** 2, rtol=1e-10)
    assert anisoplanatic_v2_factor(2.2e-6, 0.0, 0.41, 1.7, 1.8) == 1.0


def test_anisoplanatic_loss_from_headers(night):
    exposures, _ = night
    sci = exposures[1]
    assert anisoplanatic_v2_loss(sci.header, WAVE) is None  # no seeing recorded
    sci.header.update(
        {
            "HIERARCH ESO ISS AMBI FWHM START": 0.7,
            "HIERARCH ESO ISS AMBI FWHM END": 0.7,
            "HIERARCH ESO ISS AMBI THETA0 START": 1.8,
            "HIERARCH ESO ISS AMBI THETA0 END": -1.0,  # missing value, ignored
        }
    )
    loss = anisoplanatic_v2_loss(sci.header, WAVE)
    assert np.all((loss > 0) & (loss < 0.05)) and loss[0] > loss[-1]


@pytest.mark.parametrize("text", ["1±-0.1", "-1", "nan", "1±inf", "big"])
def test_unphysical_diameters_are_rejected(text):
    with pytest.raises(argparse.ArgumentTypeError):
        _diameter(text)
    with pytest.raises(argparse.ArgumentTypeError):
        _calibrator(f"HD 1={text}")


def test_diameter_parsing():
    assert _diameter("0.45±0.03") == (0.45, 0.03)
    assert _diameter("0.45+-0.03") == (0.45, 0.03)
    assert _diameter("0.45") == (0.45, None)
    assert _calibrator("HD 678") == ("HD 678", None, None)


@pytest.mark.parametrize(("name", "diameter"), [("UT1", 8.2), ("AT1", 1.8), (None, 1.8)])
def test_telescope_diameter_from_header(name, diameter):
    header = fits.Header()
    if name:
        header["HIERARCH ESO ISS CONF T1NAME"] = name
    assert telescope_diameter(header) == diameter


def test_uts_lose_less_than_ats():
    # GRAVITY+ Collaboration 2022, §4.3: the UTs keep more visibility at a given
    # separation (D/r0 enters as its −1/3 power in the exponent).
    wave = np.array([2.2e-6])
    at = anisoplanatic_v2_factor(wave, 10.0, 0.7, 2.0, 1.8)
    ut = anisoplanatic_v2_factor(wave, 10.0, 0.7, 2.0, 8.2)
    np.testing.assert_allclose(np.log(ut) / np.log(at), (8.2 / 1.8) ** (-1 / 3), rtol=1e-12)


def test_smooth_spectrum_fits_a_polynomial():
    from virgil_vlti.calibrate import smooth_spectrum

    rng = np.random.default_rng(3)
    truth = 0.8 + 0.1 * ((WAVE - 2.2e-6) / 0.2e-6) ** 2
    noisy = truth + rng.normal(0, 0.01, WAVE.size)
    noisy[4] = np.nan
    smooth = smooth_spectrum(WAVE, noisy, np.full(WAVE.size, 0.01), 2)
    assert np.isnan(smooth[4])
    ok = np.isfinite(smooth)
    assert np.std(smooth[ok] - truth[ok]) < 0.4 * np.std(noisy[ok] - truth[ok])
    # Too few channels for the degree: unchanged.
    few = np.full(WAVE.size, np.nan)
    few[:3] = 1.0
    assert np.array_equal(smooth_spectrum(WAVE, few, np.ones(WAVE.size), 5), few, equal_nan=True)


def test_smoothing_removes_calibrator_channel_noise(tmp_path):
    rng = np.random.default_rng(4)
    paths = [
        product(tmp_path, "cal1", MJD0, cal=True),
        product(tmp_path, "sci", MJD0 + 0.02, cal=False),
    ]
    sci_cal = classify([read_exposure(p) for p in paths], dual_field=True)
    cal, sci = sci_cal
    data = cal.v2["GRAVITY_SC"].data
    data["VIS2DATA"] *= 1 + rng.normal(0, 0.05, data["VIS2DATA"].shape)  # channel noise
    truth = _truth_v2(sci)
    raw = calibrate_exposure(sci, [cal], "GRAVITY_SC")["vis2"].value
    smooth = calibrate_exposure(sci, [cal], "GRAVITY_SC", smooth_degree=5)["vis2"].value
    assert np.std(smooth / truth - 1) < 0.5 * np.std(raw / truth - 1)


def test_pipeline_diameter_only_from_calibrator_products(night):
    exposures, _ = night
    cal, sci = exposures[0], exposures[1]
    for e in (cal, sci):
        e.header["HIERARCH ESO QC TF CAL_SC DIAM"] = 0.365
        e.header["HIERARCH ESO QC TF CAL_SC DIAM ERR"] = 0.005
    assert cal.pipeline_diameter == (0.365, 0.005)
    assert sci.pipeline_diameter is None  # a science product's own diameter is not used


def test_cli_uses_pipeline_diameters_unless_told_not_to(night, capsys):
    _, folder = night
    with fits.open(folder / "cal1.fits", mode="update") as h:
        h[0].header["HIERARCH ESO QC TF CAL_SC DIAM"] = 1.2
        h[0].header["HIERARCH ESO QC TF CAL_SC DIAM ERR"] = 0.05
    main([str(folder), "--out", str(folder / "a"), "--dual-field", "--smooth-tf", "3"])
    main([str(folder), "--out", str(folder / "b"), "--dual-field", "--no-pipeline-diameters"])
    (a,) = (folder / "a").glob("*/*.fits")
    (b,) = (folder / "b").glob("*/*.fits")
    assert a.with_suffix(".modes.npz").exists()  # the pipeline's diameter error gives a mode
    assert not b.with_suffix(".modes.npz").exists()
    with fits.open(a) as x, fits.open(b) as y:
        assert np.nanmedian(x["OI_VIS2"].data["VIS2DATA"] / y["OI_VIS2"].data["VIS2DATA"]) != 1
        assert "smoothed in wavelength" in str(x[0].header["HISTORY"])


def test_cal_names_match_loosely_for_diameters(night):
    exposures, _ = night
    sci, cals = exposures[1], [exposures[0], exposures[2]]
    loose = calibrate_exposure(sci, cals, "GRAVITY_SC", diameters={"cal star": 1.2})
    exact = calibrate_exposure(sci, cals, "GRAVITY_SC", diameters={"CALSTAR": 1.2})
    np.testing.assert_array_equal(loose["vis2"].value, exact["vis2"].value)


def test_transfer_function_error_is_not_averaged_down_by_binning(night):
    """The transfer function's error is common to a bin's channels.

    The synthetic transfer function drifts by the same fraction in every
    channel, so its relative scatter is the same after binning; averaging it in
    quadrature over four channels would halve it.
    """
    exposures, _ = night
    sci, cals = exposures[1], [exposures[0], exposures[2]]
    kw = {"window_um": (2.1, 2.4)}
    rel, stat = {}, {}
    for nbin in (1, 4):
        with_tf = calibrate_exposure(sci, cals, "GRAVITY_SC", nbin=nbin, error="scatter", **kw)
        without = calibrate_exposure(sci, cals, "GRAVITY_SC", nbin=nbin, error="none", **kw)
        tf = np.sqrt(with_tf["vis2"].error ** 2 - without["vis2"].error ** 2)
        rel[nbin] = tf / np.abs(with_tf["vis2"].value)
        stat[nbin] = np.median(without["vis2"].error)
    # The data's own errors are independent: binning still averages them down.
    assert 0.4 < stat[4] / stat[1] < 0.6
    est = calibrator_estimates(cals, "GRAVITY_SC", "vis2")[(28, 23)]
    # the scatter, and the calibrators' statistical error through the interpolation
    tau = np.maximum(np.std(est.values, 0, ddof=1), propagated_error(MJD0 + 0.02, est))
    expect = tau / predict(MJD0 + 0.02, est)
    np.testing.assert_allclose(rel[1][0], expect[(WAVE > 2.1e-6) & (WAVE < 2.4e-6)], rtol=1e-3)
    np.testing.assert_allclose(rel[4], rel[1][:, : rel[4].shape[1]], rtol=0.02)


def test_unbinned_errors_are_unchanged(night):
    """With nbin=1 the error is hypot(σ/g, V·τ/g), channel by channel, as before."""
    exposures, _ = night
    sci, cals = exposures[1], [exposures[0], exposures[2]]
    out = calibrate_exposure(sci, cals, "GRAVITY_SC", error="scatter")
    est = calibrator_estimates(cals, "GRAVITY_SC", "vis2")
    d = sci.v2["GRAVITY_SC"].data
    for i in range(len(d)):
        e = est[tuple(d["STA_INDEX"][i])]
        g = predict(d["MJD"][i], e)
        tau = prediction_error(e, "scatter", t=d["MJD"][i])
        cal = d["VIS2DATA"][i] / g
        expect = np.hypot(d["VIS2ERR"][i] / g, cal * tau / g)
        np.testing.assert_allclose(out["vis2"].error[i], expect, rtol=1e-10)


def test_cli_tf_error_propagated_is_statistical_only(tmp_path):
    for k, t in enumerate((MJD0, MJD0 + 0.05)):
        product(tmp_path, f"cal{k}", t, cal=True)
    product(tmp_path, "sci", MJD0 + 0.02, cal=False)
    errors = {}
    for mode in ("loo", "propagated"):
        out = tmp_path / mode
        main([str(tmp_path), "--out", str(out), "--tf-error", mode])
        (path,) = out.glob("*/*.fits")
        with fits.open(path) as h:
            errors[mode] = h["OI_VIS2"].data["VIS2ERR"].copy()
    # The statistical-only error is never larger than the leave-one-out one.
    assert (errors["propagated"] <= errors["loo"] + 1e-12).all()


# ---------------------------------------------------------------------------
# The transfer function's error as nuisance modes (--tf-modes)
# ---------------------------------------------------------------------------


def _tf_night(tmp_path, rng):
    """Three calibrators and two science exposures; calibrator closure phases offset."""
    paths = [product(tmp_path, f"cal{k}", MJD0 + 0.05 * k, cal=True) for k in range(3)] + [
        product(tmp_path, f"sci{k}", MJD0 + 0.02 + 0.05 * k, cal=False) for k in range(2)
    ]
    for p in paths[:3]:
        add_systematics(p, rng=rng, cp_offset=0.5)
    return classify([read_exposure(p) for p in paths], dual_field=True)


def test_channel_correlation():
    rng = np.random.default_rng(1)
    common = rng.normal(size=(2000, 1)) * np.linspace(1, 3, 8)  # one shape, varying size
    np.testing.assert_allclose(channel_correlation(common), 1.0, atol=1e-12)
    # Independent channels, or half common and half independent: a rank-one fit
    # puts a fraction ρ + (1 − ρ)/m of the variance in the mode (on average).
    assert np.mean(channel_correlation(rng.normal(size=(2000, 8))) ** 2) < 1 / 8 + 0.05
    mixed = rng.normal(size=(2000, 1)) + rng.normal(size=(2000, 8))
    r2 = channel_correlation(mixed) ** 2
    assert np.all((r2 > 0.5 - 0.05) & (r2 < 0.5 + 0.5 / 8 + 0.05))
    # A single channel, or a single calibrator, is all common.
    np.testing.assert_array_equal(channel_correlation(rng.normal(size=(5, 1))), [1.0])
    np.testing.assert_array_equal(channel_correlation(rng.normal(size=(1, 4))), np.ones(4))
    # Anticorrelated channels keep opposite signs, so the mode reproduces −1.
    x = rng.normal(size=500)
    r = channel_correlation(np.stack([x, -x], 1))
    np.testing.assert_allclose(r[0] * r[1], -1.0, atol=1e-9)


@pytest.mark.parametrize("error", ["scatter", "loo"])
def test_tf_modes_replace_the_diagonal_term(tmp_path, error):
    """With modes, the transfer function's error leaves the diagonal and becomes the mode.

    The synthetic transfer function drifts by one fraction in every channel and
    the calibrators' closure phases by one offset per triangle, so the error is
    fully common to the channels: the errors are those without any transfer
    function term, and the modes carry all of it.
    """
    exposures = _tf_night(tmp_path, np.random.default_rng(2))
    sci, cals = exposures[3], exposures[:3]
    kw = {"window_um": (2.1, 2.4), "nbin": 4, "error": error}
    diag = calibrate_exposure(sci, cals, "GRAVITY_SC", **kw)
    modal = calibrate_exposure(sci, cals, "GRAVITY_SC", tf_modes=True, **kw)
    none = calibrate_exposure(sci, cals, "GRAVITY_SC", **{**kw, "error": "none"})
    assert diag["vis2"].tf_mode is None
    for name in ("vis2", "t3phi"):
        np.testing.assert_allclose(modal[name].value, diag[name].value)
        # the diagonal keeps only what is not in the mode: the calibrators' statistical
        # floor where it exceeds the systematic's common part
        assert np.all(modal[name].error >= none[name].error - 1e-12)
        assert np.all(modal[name].error <= diag[name].error + 1e-12)
        assert np.all(diag[name].error >= none[name].error)
        assert np.median(diag[name].error / none[name].error) > 1.01
        tf = np.sqrt(diag[name].error ** 2 - none[name].error ** 2)
        if name == "vis2":
            tf = tf / np.abs(diag[name].value)  # relative
        # the mode is the systematic's common part: never more than the whole term
        assert np.all(modal[name].tf_mode <= tf * (1 + 1e-4))
    # Other observables keep the diagonal term.
    np.testing.assert_allclose(modal["t3amp"].error, diag["t3amp"].error)
    with pytest.raises(ValueError, match="tf_modes"):
        calibrate_exposure(sci, cals, "GRAVITY_SC", tf_modes=True, error="propagated")


def test_tf_modes_are_named_by_night_and_baseline(tmp_path):
    exposures = _tf_night(tmp_path, np.random.default_rng(3))
    cals = exposures[:3]
    by_file = {}
    for sci in exposures[3:]:
        out = calibrate_exposure(sci, cals, "GRAVITY_SC", nbin=4, error="scatter", tf_modes=True)
        by_file[sci.path] = (sci, out, transfer_function_modes(sci, out))
    names = [set(m) for _, _, m in by_file.values()]
    assert names[0] == names[1]  # one parameter per (night, baseline), shared by the files
    night = exposures[3].night
    stations = {tuple(s) for s in exposures[3].v2["GRAVITY_SC"].data["STA_INDEX"]}
    triangles = {tuple(s) for s in exposures[3].t3["GRAVITY_SC"].data["STA_INDEX"]}
    assert {n for n in names[0] if n.startswith("tf-vis2")} == {
        f"tf-vis2-base{a}-{b}@{night}" for a, b in stations
    }
    assert {n for n in names[0] if n.startswith("tf-t3phi")} == {
        f"tf-t3phi-tri{a}-{b}-{c}@{night}" for a, b, c in triangles
    }
    for sci, out, modes in by_file.values():
        sta = sci.v2["GRAVITY_SC"].data["STA_INDEX"]
        a, b = sta[0]
        arr = modes[f"tf-vis2-base{a}-{b}@{night}"]["vis2"]
        rows = np.all(sta == sta[0], 1)
        np.testing.assert_allclose(arr[rows], 0.5 * out["vis2"].tf_mode[rows])  # d log|V|
        assert not np.any(arr[~rows])
        tri = sci.t3["GRAVITY_SC"].data["STA_INDEX"]
        name = "tf-t3phi-tri" + "-".join(str(x) for x in tri[0]) + f"@{night}"
        cp = modes[name]["t3phi"]
        rows = np.all(tri == tri[0], 1)
        np.testing.assert_allclose(cp[rows], out["t3phi"].tf_mode[rows])  # degrees
        assert 0.1 < np.median(cp[rows]) < 2.0  # the injected 0.5° offsets


def test_cli_writes_tf_modes_and_drops_the_diagonal_term(tmp_path):
    _tf_night(tmp_path, np.random.default_rng(4))
    common = [str(tmp_path), "--dual-field", "--bin", "4", "--tf-error", "scatter"]
    main([*common, "--out", str(tmp_path / "modal"), "--tf-modes"])
    main([*common, "--out", str(tmp_path / "diag")])
    files = sorted((tmp_path / "modal").glob("*/*.fits"))
    assert len(files) == 2
    keys = []
    for path in files:
        with np.load(path.with_suffix(".tf_modes.npz")) as d:
            names = {k for k in d.files if not k.startswith("__")}
            units = json.loads(str(d["__units__"]))
            assert all(units[k].startswith("deg" if "/t3phi/" in k else "d log|V|") for k in names)
            keys.append(names)
            assert all(k.startswith(("GRAVITY_SC/vis2/tf-vis2-", "GRAVITY_SC/t3phi/tf-t3phi-"))
                       for k in names)  # fmt: skip
        with (
            fits.open(path) as h,
            fits.open(tmp_path / "diag" / path.relative_to(tmp_path / "modal")) as g,
        ):
            assert np.all(h["OI_VIS2"].data["VIS2ERR"] < g["OI_VIS2"].data["VIS2ERR"])
            assert "modes" in str(h[0].header["HISTORY"])
            assert "shape poly1" in str(h[0].header["HISTORY"])
    assert keys[0] == keys[1]
    assert not list((tmp_path / "diag").glob("*/*.tf_modes.npz"))
    with pytest.raises(SystemExit, match="loo or scatter"):
        main([*common[:-2], "--out", str(tmp_path / "x"), "--tf-modes", "--tf-error", "none"])
    with pytest.raises(SystemExit, match="not both"):
        main([*common, "--out", str(tmp_path / "y"), "--tf-modes", "--template", "any.npz"])
    with pytest.raises(SystemExit, match="only with --tf-modes"):
        main([*common, "--out", str(tmp_path / "z"), "--tf-mode-shape", "poly1"])


def test_regularize_shape_is_smooth_and_conserves_variance():
    wave = np.linspace(2.0e-6, 2.4e-6, 12)
    rng = np.random.default_rng(5)
    total = 1 + 0.2 * rng.normal(size=12)
    common = 0.8 * total
    common[3] = np.nan  # a flagged channel stays NaN and out of the fit
    total[3] = np.nan
    flat, rem, added = regularize_shape(common, total, wave, "flat")
    ok = np.isfinite(total)
    np.testing.assert_allclose(flat[ok], np.nanmean(common))
    assert np.isnan(flat[3]) and np.isnan(rem[3]) and added[3] == 0
    # shape² + remainder² is τ², or shape² where that is larger (remainder clipped)
    np.testing.assert_allclose(flat[ok] ** 2 + rem[ok] ** 2, np.maximum(total, flat)[ok] ** 2)
    np.testing.assert_allclose(added[ok], np.clip(flat**2 - total**2, 0, None)[ok])
    assert np.any(added > 0) and np.any(rem[ok] > 0)  # both branches exercised
    # A polynomial shape reproduces a polynomial exactly, and fits by least squares.
    line = 0.5 + 0.3 * (wave - 2.2e-6) / 0.4e-6
    for degree in (1, 2):
        fitted, rem, added = regularize_shape(line, np.full(12, 2.0), wave, f"poly{degree}")
        np.testing.assert_allclose(fitted, line, rtol=1e-10)
        np.testing.assert_allclose(fitted**2 + rem**2, 4.0)
        assert not np.any(added)
    quad = line + 0.1 * rng.normal(size=12)
    fitted, _, _ = regularize_shape(quad, np.full(12, 2.0), wave, "poly1")
    np.testing.assert_allclose(fitted, np.polyval(np.polyfit(wave, quad, 1), wave), rtol=1e-8)
    raw, rem, _ = regularize_shape(common, total, wave, "raw")
    np.testing.assert_array_equal(raw, common)
    np.testing.assert_allclose(rem[ok], 0.6 * total[ok])
    with pytest.raises(ValueError, match="shape"):
        regularize_shape(common, total, wave, "cubic")


def test_regularize_shape_keeps_the_sign_of_anticorrelated_channels():
    """A mode that changes sign across the band keeps its size and its signs.

    A signed mean of r·τ shrinks such a mode towards zero and makes it
    one-signed, sending its variance to the independent remainder, which
    averages down over channels and exposures: errors too small.
    """
    wave = np.linspace(2.0e-6, 2.4e-6, 12)
    total = np.full(12, 1.0)
    common = np.where(np.arange(12) < 6, 0.9, -0.9) * (1 + 0.1 * np.sin(np.arange(12)))
    for shape in ("flat", "poly1"):
        fitted, rem, added = regularize_shape(common, total, wave, shape)
        np.testing.assert_array_equal(np.sign(fitted), np.sign(common))
        assert np.mean(np.abs(fitted)) > 0.8  # not shrunk towards zero
        # the mode reproduces the anticorrelation, and the variance is conserved
        assert fitted[0] * fitted[-1] < -0.5
        np.testing.assert_allclose(fitted**2 + rem**2, np.maximum(total, np.abs(fitted)) ** 2)
    # the same through channel_correlation: a tilt about the band centre
    rng = np.random.default_rng(3)
    res = rng.normal(size=(300, 1)) * np.linspace(-1, 1, 12)[None] + 0.2 * rng.normal(
        size=(300, 12)
    )
    r = channel_correlation(res)
    tau = np.sqrt(np.mean(res**2, 0))
    fitted, _, _ = regularize_shape(r * tau, tau, wave, "flat")
    assert np.all(np.sign(fitted[np.abs(r) > 0.3]) == np.sign(r[np.abs(r) > 0.3]))
    assert np.mean(np.abs(fitted)) > 0.5 * np.mean(np.abs(r * tau))


def test_shape_slope_significance():
    wave = np.linspace(2.0e-6, 2.4e-6, 20)
    rng = np.random.default_rng(4)
    x = (wave - wave.mean()) / np.ptp(wave)
    assert shape_slope_significance(1 + 0.03 * rng.normal(size=20), wave) < 3
    assert shape_slope_significance(1 + 0.6 * x + 0.03 * rng.normal(size=20), wave) > 10
    # the sign of a channel's mode does not matter, only its size
    signed = np.where(np.arange(20) < 10, 1.0, -1.0) * (1 + 0.6 * x + 0.03 * rng.normal(size=20))
    assert shape_slope_significance(signed, wave) > 10
    assert shape_slope_significance(np.ones(3), wave[:3]) == 0.0  # too few channels
    assert shape_slope_significance(np.ones(20), wave) == 0.0  # no scatter


@pytest.mark.parametrize("shape", ["flat", "poly1", "raw"])
def test_tf_mode_shapes_through_calibration(tmp_path, shape):
    """The stored mode is flat (or linear) in each row, and no channel loses variance.

    Per-channel noise on the calibrators makes τ and r vary from channel to
    channel, so the raw shape is not flat and the smoothing has work to do.
    """
    paths = [product(tmp_path, f"cal{k}", MJD0 + 0.05 * k, cal=True) for k in range(5)] + [
        product(tmp_path, "sci", MJD0 + 0.07, cal=False)
    ]
    rng = np.random.default_rng(6)
    for p in paths[:5]:
        add_systematics(p, rng=rng, baseline=0.03, white=2.0, cp_offset=1.0, cp_white=2.0)
    exposures = classify([read_exposure(p) for p in paths], dual_field=True)
    sci, cals = exposures[5], exposures[:5]
    kw = {"nbin": 4, "error": "loo"}
    diag = calibrate_exposure(sci, cals, "GRAVITY_SC", **kw)
    none = calibrate_exposure(sci, cals, "GRAVITY_SC", **{**kw, "error": "none"})
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        modal = calibrate_exposure(
            sci, cals, "GRAVITY_SC", tf_modes=True, tf_mode_shape=shape, **kw
        )
    clipped = [w for w in caught if "clipped to zero" in str(w.message)]
    assert bool(clipped) == (shape != "raw")  # said so when a channel's variance grows
    wave = modal.wave
    for name in ("vis2", "t3phi"):
        m, ok = modal[name], ~modal[name].flag
        mode = m.tf_mode if name == "t3phi" else m.tf_mode * np.abs(m.value)  # absolute
        tf2 = diag[name].error ** 2 - none[name].error ** 2
        rem2 = m.error**2 - none[name].error ** 2
        # Each channel's total transfer-function variance is kept, or grows where clipped.
        np.testing.assert_allclose(
            (rem2 + mode**2)[ok], np.maximum(tf2, mode**2)[ok], rtol=1e-6, atol=1e-12
        )
        for row, keep in zip(m.tf_mode, ok, strict=True):
            row = row[keep]
            spread = np.std(row) / np.mean(np.abs(row))
            if shape != "raw":  # smooth in size; channels keep their sign (anticorrelation)
                row = np.abs(row)
            if shape == "flat":
                np.testing.assert_allclose(row, row[0], rtol=1e-12)
            elif shape == "poly1":
                fit = np.polyval(np.polyfit(wave[keep], row, 1), wave[keep])
                np.testing.assert_allclose(row, fit, rtol=1e-8)
                assert spread > 1e-6
            else:
                assert spread > 0.02  # the raw shape scatters from channel to channel
        if shape == "flat":  # both kept (remainder > 0) and clipped channels occur
            assert np.any((mode**2 > 1.01 * tf2)[ok]) and np.any((rem2 > 0.01 * tf2)[ok])
        if shape == "raw":  # the old split: nothing clipped, variance exactly conserved
            np.testing.assert_allclose((rem2 + mode**2)[ok], tf2[ok], rtol=1e-6, atol=1e-12)


def _gls_pulls(shape, colour, rng, n_real=1000, nchan=24, ncal=7, nbase=6, nrow=40):
    """Pulls of a GLS fit through transfer-function modes, numpy only.

    Per baseline, ``ncal`` calibrators each draw a transfer-function error
    (1σ = 1, with spectral shape 1 + colour·x) plus channel noise; τ, r and
    the mode come from them as in :func:`calibrate_exposure`. The science data
    (``nrow`` rows sharing one realized error per baseline, white noise 1)
    are fitted for two linear parameters with C = diag(σ² + remainder²) +
    mode modeᵀ. Returns the pulls' rms per parameter (1 if the errors are right).
    """
    wave = np.linspace(2.0e-6, 2.4e-6, nchan)
    x = (wave - wave.mean()) / np.ptp(wave)
    colour_shape = 1 + colour * x
    out = []
    for _ in range(n_real):
        F, g = np.zeros((2, 2)), np.zeros(2)
        for k in range(nbase):
            res = rng.normal(size=(ncal, 1)) * colour_shape + rng.normal(0, 0.3, (ncal, nchan))
            tau = np.sqrt(np.mean(res**2, 0))
            mode, rem, _ = regularize_shape(tau * channel_correlation(res), tau, wave, shape)
            J = np.tile(np.stack([np.cos(0.3 * k + 2 * x), np.sin(0.7 * k + 3 * x)], -1), (nrow, 1))
            y = rng.normal() * np.tile(colour_shape, nrow) + rng.normal(size=nrow * nchan)
            d, u = 1 + np.tile(rem**2, nrow), np.tile(mode, nrow)
            Ju, uy = J.T @ (u / d), u @ (y / d)
            s = 1 + u @ (u / d)  # Woodbury, one mode per baseline
            F += J.T @ (J / d[:, None]) - np.outer(Ju, Ju) / s
            g += J.T @ (y / d) - Ju * uy / s
        out.append(np.linalg.solve(F, g) / np.sqrt(np.diag(np.linalg.inv(F))))
    return np.sqrt(np.mean(np.square(out), 0))


def test_smoothed_tf_modes_give_calibrated_pulls():
    """Achromatic transfer-function errors: flat modes give pulls of 1, raw ones do not.

    Seven calibrators estimate each channel's τ and r with about 25% scatter;
    a raw mode leaves (shape scatter × realized error) unmodelled, which the
    hundreds of data per mode then fit as signal.
    """
    flat = _gls_pulls("flat", 0.0, np.random.default_rng(7))
    raw = _gls_pulls("raw", 0.0, np.random.default_rng(7))
    assert np.all((flat > 0.88) & (flat < 1.1)), flat  # about 0.95: clipping is conservative
    assert np.all(raw > 1.2) and np.max(raw) > 1.5, raw  # about 1.6 and 1.35


def test_poly1_tf_modes_follow_a_coloured_transfer_function():
    """A transfer-function error with a linear colour needs poly1: flat misses its slope."""
    poly1 = _gls_pulls("poly1", 1.5, np.random.default_rng(8))
    flat = _gls_pulls("flat", 1.5, np.random.default_rng(8))
    assert np.all(poly1 < 1.4), poly1  # about 1.23: the fitted slope's own noise
    assert np.all(flat > 3.0), flat  # about 11


def test_flat_tf_modes_at_a_realistic_colour_are_not_calibrated():
    """At colour 0.36 (about ±18 per cent across K, as λ⁻² losses give) flat is too small.

    Pull rms (sep-like, PA-like parameters), 400 realizations: flat ~3.1,
    poly1 ~1.3. The default ``flat`` is validated only on achromatic errors
    (PROVENANCE.md); this documents where it stops being calibrated.
    """
    flat = {
        c: _gls_pulls("flat", c, np.random.default_rng(11), n_real=200) for c in (0.0, 0.2, 0.36)
    }
    poly1 = {c: _gls_pulls("poly1", c, np.random.default_rng(11), n_real=200) for c in (0.0, 0.36)}
    assert np.all(flat[0.0] < 1.15)
    assert np.all(flat[0.2] > 1.4), flat[0.2]  # about 1.9
    assert np.all(flat[0.36] > 2.2), flat[0.36]  # about 3.1
    assert np.all(poly1[0.36] < 1.5), poly1[0.36]  # about 1.3, at every colour
    assert np.all(np.diff([flat[0.0], flat[0.2], flat[0.36]], axis=0) > 0)  # grows with colour


def _night_pulls(tmp_path, shape, colour, rng, n_real=120, ncal=12, nbin=4):
    """Pulls through :func:`calibrate_exposure`, on synthetic nights with an injected error.

    ``ncal`` calibrators and one science exposure, written once. In each
    realization every calibrator and the science exposure get their own
    gain per baseline (1σ = 0.1 on log V², spectral shape 1 + colour·x) plus
    channel noise, so the transfer function interpolated to the science
    exposure is wrong by a known, realized amount. The calibrated V² are
    fitted for two parameters (a scale and a scale tilt of the true V²),
    generalized least squares with C = diag(error²) + Σ mode modeᵀ, the modes
    along the true V² (as :class:`virgil.OIData` takes model-relative modes).
    Returns the pulls' rms per parameter.
    """
    ins = "GRAVITY_SC"
    paths = [
        product(tmp_path, f"cal{k}", MJD0 + 0.01 * k, cal=True, nframes=1) for k in range(ncal)
    ] + [product(tmp_path, "sci", MJD0 + 0.01 * (ncal // 2) + 0.005, cal=False, nframes=1)]
    *cals, sci = classify([read_exposure(p) for p in paths], dual_field=True)
    x = (WAVE - WAVE.mean()) / np.ptp(WAVE)
    gain_shape = 1 + colour * x
    base = {id(e): e.v2[ins].data["VIS2DATA"].copy() for e in [*cals, sci]}
    stations = {id(e): [tuple(s) for s in e.v2[ins].data["STA_INDEX"]] for e in [*cals, sci]}
    names = sorted(set(stations[id(sci)]))
    truth = bin_channels(_truth_v2(sci), np.ones(len(WAVE), bool), nbin).mean(-1)
    pulls = []
    for _ in range(n_real):
        for e, white in [*[(e, 0.03) for e in cals], (sci, 0.01)]:
            gain = dict(zip(names, rng.normal(0, 0.1, len(names)), strict=True))
            log_gain = np.array([gain[s] for s in stations[id(e)]])[:, None]
            e.v2[ins].data["VIS2DATA"][:] = base[id(e)] * np.exp(log_gain * gain_shape) + (
                rng.normal(0, white, base[id(e)].shape)
            )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cal = calibrate_exposure(
                sci, cals, ins, nbin=nbin, error="loo", tf_modes=True, tf_mode_shape=shape
            )
        o, ok = cal["vis2"], ~cal["vis2"].flag
        modes = np.stack(
            [
                np.where(np.array([s == n for s in stations[id(sci)]])[:, None], o.tf_mode, 0)
                * truth
                for n in names
            ],
            -1,
        )[ok]
        C = np.diag(o.error[ok] ** 2) + modes @ modes.T
        xs = np.tile((cal.wave - cal.wave.mean()) / np.ptp(cal.wave), (len(o.value), 1))[ok]
        J = np.stack([truth[ok], truth[ok] * xs], -1)
        F = J.T @ np.linalg.solve(C, J)
        theta = np.linalg.solve(F, J.T @ np.linalg.solve(C, (o.value - truth)[ok]))
        pulls.append(theta / np.sqrt(np.diag(np.linalg.inv(F))))
    return np.sqrt(np.mean(np.square(pulls), 0))


def test_tf_mode_pulls_through_calibrate_exposure(tmp_path):
    """The shipped path (`_apply`: relative units, loo residuals, binned estimates) is calibrated.

    Pull rms with 12 calibrators, 120 realizations (the band is loose: the
    pulls inherit the small-sample error of the scatter, 1.1 at 30
    calibrators): flat ~1.35 / 1.2 and raw ~1.9 / 0.95 at colour 0; at
    colour 0.36 flat's slope pull is ~1.9 and poly1's ~1.4. A mis-wired split
    (swapped total and common, a wrong level, the wrong wavelength grid)
    moves them.
    """
    rng = np.random.default_rng(5)
    for name in ("flat", "raw", "coloured"):
        (tmp_path / name).mkdir()
    flat = _night_pulls(tmp_path / "flat", "flat", 0.0, rng)
    raw = _night_pulls(tmp_path / "raw", "raw", 0.0, rng)
    assert np.all((flat > 0.8) & (flat < 1.9)), flat
    assert raw[0] > 1.25 * flat[0], (raw, flat)  # the noisy raw shape is worse
    coloured = _night_pulls(tmp_path / "coloured", "flat", 0.36, rng)
    assert coloured[1] > flat[1] + 0.3, (coloured, flat)  # the slope pull grows with colour


def test_flat_mode_warns_about_a_coloured_transfer_function(tmp_path):
    """Calibrators whose common error has a colour warn that ``flat`` cannot follow it."""
    ins = "GRAVITY_SC"

    def run(colour, shape="flat"):
        sub = tmp_path / f"{colour}{shape}"
        sub.mkdir()
        rng = np.random.default_rng(8)
        paths = [product(sub, f"cal{k}", MJD0 + 0.01 * k, cal=True, nframes=1) for k in range(8)]
        paths.append(product(sub, "sci", MJD0 + 0.035, cal=False, nframes=1))
        x = (WAVE - WAVE.mean()) / np.ptp(WAVE)
        for p in paths[:-1]:
            add_systematics(p, rng=rng, baseline=0.1, white=2.0, shape=1 + colour * x)
        *cals, sci = classify([read_exposure(p) for p in paths], dual_field=True)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            calibrate_exposure(
                sci, cals, ins, nbin=2, error="loo", tf_modes=True, tf_mode_shape=shape
            )
        return [str(w.message) for w in caught if "significant colour" in str(w.message)]

    assert not run(0.0)
    hit = run(1.0)
    assert len(hit) == 1 and "VIS2DATA" in hit[0] and "poly1" in hit[0]
    assert not run(1.0, "poly1")  # a shape that can follow it does not warn


def test_relative_diameter_modes_are_model_relative(tmp_path):
    from virgil_vlti.calibrate import _unit_amplitudes

    paths = [
        product(tmp_path, "cal1", MJD0, cal=True, diameter=1.0),
        product(tmp_path, "sci", MJD0 + 0.02, cal=False),
        product(tmp_path, "cal2", MJD0 + 0.05, cal=True, diameter=1.0),
    ]
    exposures = classify([read_exposure(p) for p in paths], dual_field=True)
    sci, cals = exposures[1], [exposures[0], exposures[2]]
    diam, sigma = {"CALSTAR": 1.0}, {"CALSTAR": 0.1}
    raw = sci.v2["GRAVITY_SC"].data["VIS2DATA"].copy()
    (absolute,) = diameter_modes(sci, cals, "GRAVITY_SC", diam, sigma).values()
    (relative,) = diameter_modes(sci, cals, "GRAVITY_SC", diam, sigma, relative=True).values()
    np.testing.assert_array_equal(sci.v2["GRAVITY_SC"].data["VIS2DATA"], raw)  # not modified
    base = calibrate_exposure(sci, cals, "GRAVITY_SC", diameters=diam)
    # d V² = 2 V² d log|V|; d T3AMP = T3AMP d ln T3AMP (second order in the secant).
    np.testing.assert_allclose(
        2 * base["vis2"].value * relative["vis2"], absolute["vis2"], rtol=1e-2
    )
    np.testing.assert_allclose(
        base["t3amp"].value * relative["t3amp"], absolute["t3amp"], rtol=1e-2
    )
    # A transfer function only: the same for any science data, unlike the
    # absolute mode, even when binning averages channels of different values
    # and signs (where a ratio of calibrated data would not be).
    other = _unit_amplitudes(sci, "GRAVITY_SC")
    rng = np.random.default_rng(0)
    other.v2["GRAVITY_SC"].data["VIS2DATA"] = rng.uniform(-0.2, 1.0, raw.shape)
    t3amp = other.t3["GRAVITY_SC"].data["T3AMP"]
    other.t3["GRAVITY_SC"].data["T3AMP"] = rng.uniform(-0.2, 1.0, t3amp.shape)
    binned = {}
    for name, data in (("sci", sci), ("other", other)):
        for rel in (False, True):
            (binned[name, rel],) = diameter_modes(
                data, cals, "GRAVITY_SC", diam, sigma, relative=rel, nbin=2
            ).values()
    for obs in ("vis2", "t3amp"):
        np.testing.assert_allclose(binned["other", True][obs], binned["sci", True][obs], rtol=1e-12)
    assert not np.allclose(binned["other", False]["vis2"], binned["sci", False]["vis2"], rtol=1e-2)
    assert np.all(relative["vis2"] < 0)


def test_injection_scatter(tmp_path):
    from virgil_vlti.calibrate import injection_scatter
    from virgil_vlti.products import read_flux

    rng = np.random.default_rng(2)
    levels = []
    paths = []
    for k in range(3):
        level = np.exp(rng.normal(0.0, 0.1, (12, 1))) * np.array([1.0, 2.0, 1.0, 1.0] * 3)[:, None]
        levels.append(level[:, 0])
        paths.append(
            write_product(
                tmp_path / f"f{k}.fits",
                mjd=MJD0 + 0.01 * k,
                vis2=lambda u, v, w, t: np.ones_like(w),
                t3phi=lambda u1, v1, u2, v2, w, t: 0 * w,
                wave=WAVE,
                flux=level,
            )
        )
    exposures = [read_exposure(p) for p in paths]
    out = injection_scatter(exposures, "GRAVITY_SC")
    # The synthetic flux is level × (λ/2.2 µm): its band mean is level × mean(λ)/2.2 µm.
    band = np.concatenate(levels) * WAVE.mean() / 2.2e-6
    assert out["n"] == 36
    np.testing.assert_allclose(out["mean"], band.mean())
    np.testing.assert_allclose(out["sd"], band.std(ddof=1))
    np.testing.assert_allclose(out["log_sd"], np.log(band).std(ddof=1))
    assert out["n_dropped"] == 0
    _, sta = read_flux(exposures[0], "GRAVITY_SC", stations=True)
    assert set(out["per_telescope"]) == set(sta.tolist())
    # Per telescope, the doubled station's offset drops out of the width and
    # shows in its mean.
    second = band[np.tile(np.arange(12) % 4 == 1, 3)]
    station = out["per_telescope"][int(sta[1])]
    np.testing.assert_allclose(station["mean"], second.mean())
    np.testing.assert_allclose(station["sd"], second.std(ddof=1))
    np.testing.assert_allclose(station["log_sd"], np.log(second).std(ddof=1))
    assert station["n"] == 9
    assert station["log_sd"] < out["log_sd"]
    assert station["mean"] > 1.5 * out["per_telescope"][int(sta[0])]["mean"]
    # A lost fibre (band mean ≤ 0) is dropped, with a warning and a count.
    lost = levels[0][:, None].copy()
    lost[2] = -0.01
    path = write_product(
        tmp_path / "lost.fits",
        mjd=MJD0 + 0.05,
        vis2=lambda u, v, w, t: np.ones_like(w),
        t3phi=lambda u1, v1, u2, v2, w, t: 0 * w,
        wave=WAVE,
        flux=lost,
    )
    with pytest.warns(UserWarning, match="1 OI_FLUX rows .* censored"):
        out = injection_scatter([*exposures, read_exposure(path)], "GRAVITY_SC")
    assert out["n_dropped"] == 1 and out["n"] == 47


def test_bad_keeps_north_south_baselines():
    """A row is empty only when u and v are both zero; u = 0 alone is a N–S baseline."""
    from types import SimpleNamespace

    from virgil_vlti.calibrate import _bad

    n = 4
    for ucol, vcol in (("UCOORD", "VCOORD"), ("U1COORD", "V1COORD")):
        data = {
            "FLAG": np.zeros((2, n), bool),
            "VIS2DATA": np.ones((2, n)),
            "VIS2ERR": np.full((2, n), 0.01),
            ucol: np.array([0.0, 0.0]),
            vcol: np.array([50.0, 0.0]),
        }
        table = SimpleNamespace(data=data)
        assert not _bad(table, 0, "VIS2DATA", "VIS2ERR", ucol).any()
        assert _bad(table, 1, "VIS2DATA", "VIS2ERR", ucol).all()


def test_linear_interpolation_skips_invalid_calibrators():
    """A NaN calibrator channel is bridged by the valid ones, as in the Gaussian branch."""
    from virgil_vlti.calibrate import Estimates, interpolation_weights, predict, propagated_error

    times = np.array([0.0, 1.0, 2.0])
    values = np.array([[1.0, 1.0], [np.nan, 2.0], [3.0, 3.0]])
    errors = np.array([[0.1, 0.1], [np.nan, 0.2], [0.3, 0.3]])
    est = Estimates(times, values, errors)
    t = 1.5
    np.testing.assert_allclose(predict(t, est), [2.5, 2.5])
    w = interpolation_weights(t, est)
    np.testing.assert_allclose(w[:, 0], [0.25, 0.0, 0.75])
    np.testing.assert_allclose(w[:, 1], [0.0, 0.5, 0.5])
    np.testing.assert_allclose(np.nansum(w * values, 0), predict(t, est))
    np.testing.assert_allclose(
        propagated_error(t, est), [np.hypot(0.25 * 0.1, 0.75 * 0.3), np.hypot(0.1, 0.15)]
    )
    # Phasors too, and a channel with no valid calibrator stays NaN.
    phasors = np.exp(1j * np.deg2rad(np.array([[0.0, np.nan], [np.nan, np.nan], [20.0, np.nan]])))
    out = predict(t, Estimates(times, phasors, errors))
    np.testing.assert_allclose(
        np.angle(out[0], deg=True),
        np.degrees(np.angle(0.25 + 0.75 * np.exp(1j * np.deg2rad(20.0)))),
    )
    assert np.isnan(out[1])


def test_calibrator_model_null_is_masked(night, monkeypatch):
    """A calibrator model of zero (a uniform-disk null) masks the channel, not inf."""
    import dataclasses

    from virgil_vlti import calibrate as cal_mod

    def model(table, wave, diameter):
        m = np.ones((len(table.data), wave.size))
        m[:, 5] = 0.0
        m[:, 6] = np.nan
        return m

    spec = dataclasses.replace(cal_mod.OBSERVABLES["vis2"], model=model)
    monkeypatch.setitem(cal_mod.OBSERVABLES, "vis2", spec)
    exposures, _ = night
    est = calibrator_estimates([exposures[0], exposures[2]], "GRAVITY_SC", "vis2")[(28, 23)]
    assert np.isnan(est.values[:, 5:7]).all() and np.isnan(est.errors[:, 5:7]).all()
    assert np.isfinite(np.delete(est.values, [5, 6], axis=1)).all()


def test_readme_states_tf_mode_shape_default():
    """The README's table and caveat both give calibrate_exposure's --tf-mode-shape default."""
    import inspect
    from pathlib import Path

    from virgil_vlti.calibrate import TF_MODE_SHAPES

    default = inspect.signature(calibrate_exposure).parameters["tf_mode_shape"].default
    readme = (Path(__file__).parents[1] / "README.md").read_text()
    assert f"| `--tf-mode-shape` | `{default}` |" in readme
    assert f"smooths the shape (default `{default}`)" in readme
    for other in set(TF_MODE_SHAPES) - {default}:
        assert f"(default `{other}`)" not in readme


def test_flux_ratio_calibrator_slopes(night):
    exposures, _ = night
    base = flux_ratio(exposures, "GRAVITY_SC")
    np.testing.assert_array_equal(  # default, and a zero slope, change nothing
        flux_ratio(exposures, "GRAVITY_SC", calibrator_slopes=None)["mean"], base["mean"]
    )
    np.testing.assert_allclose(
        flux_ratio(exposures, "GRAVITY_SC", calibrator_slopes=0.0)["mean"], base["mean"]
    )
    wave = base["wavelength_um"]
    out = flux_ratio(exposures, "GRAVITY_SC", calibrator_slopes=-3.5)["mean"]
    expect = base["mean"] * (wave / np.median(wave)) ** -3.5
    np.testing.assert_allclose(out / np.nanmedian(out), expect / np.nanmedian(expect), rtol=1e-6)
    assert not np.allclose(out, base["mean"])


def test_flux_ratio_calibrator_slopes_physics(tmp_path):
    """CAL = F_true λ^s T, SCI = flat T: slope s flattens the ratio, -s does not; also binned."""
    s = 3.0
    shape = lambda w: (w / 2.2e-6) ** s  # noqa: E731
    paths = []
    for name, t, cal in (("c", MJD0, True), ("s", MJD0 + 0.02, False)):
        paths.append(
            write_product(
                tmp_path / f"{name}.fits",
                mjd=t,
                vis2=lambda u, v, w, tk: np.ones_like(w),
                t3phi=lambda u1, v1, u2, v2, w, tk: 0 * w,
                wave=WAVE,
                flux=shape if cal else 1.0,
                sobj=(0.0, 0.0) if cal else (210.0, -710.0),
                catg="DUAL_CAL_VIS" if cal else "DUAL_SCI_VIS",
                target="CALSTAR" if cal else "SCISTAR",
            )
        )
    exposures = classify([read_exposure(p) for p in paths], dual_field=True)
    for nbin in (1, 2):
        good = flux_ratio(exposures, "GRAVITY_SC", nbin=nbin, calibrator_slopes=s)["mean"]
        bad = flux_ratio(exposures, "GRAVITY_SC", nbin=nbin, calibrator_slopes=-s)["mean"]
        np.testing.assert_allclose(good, 1.0, rtol=0, atol=1e-5 if nbin == 1 else 0.02)
        assert np.ptp(bad) > 0.3


def test_flux_ratio_calibrator_slopes_guards(night):
    import dataclasses

    exposures, _ = night
    nights = {e.night for e in exposures}
    with pytest.raises(KeyError, match="no entry"):
        flux_ratio(exposures, "GRAVITY_SC", calibrator_slopes={"1999-01-01": -3.0})
    flux_ratio(exposures, "GRAVITY_SC", calibrator_slopes=dict.fromkeys(nights, -3.0))
    other = [
        dataclasses.replace(e, target="OTHER") if i == 2 else e for i, e in enumerate(exposures)
    ]
    with pytest.raises(ValueError, match="one calibrator per night"):
        flux_ratio(other, "GRAVITY_SC", calibrator_slopes=-3.0)
    flux_ratio(other, "GRAVITY_SC")  # fine without slopes


# ---------------------------------------------------------------------------
# Review fixes: unrecorded keys, one calibrator, bad errors
# ---------------------------------------------------------------------------


def test_unrecorded_key_is_refused_unless_allowed(tmp_path):
    sci, cals = _with_cal_header(tmp_path, {})
    del cals[0].header["ESO INS POLA MODE"]
    with pytest.raises(ValueError, match=r"cannot be checked.*\[not recorded\]"):
        calibrate_exposure(sci, cals, "GRAVITY_SC")
    with pytest.raises(ValueError, match="cannot be checked"):  # strict=False is not the opt-out
        calibrate_exposure(sci, cals, "GRAVITY_SC", strict=False)
    with pytest.warns(UserWarning, match=r"\[not recorded\]"):
        calibrate_exposure(sci, cals, "GRAVITY_SC", allow_unrecorded=True)


def test_unrecorded_pipeline_parameter_is_refused_on_either_side(tmp_path):
    sci, cals = _with_cal_header(tmp_path, {"ESO PRO REC1 PARAM1 NAME": "other"})
    with pytest.raises(ValueError, match=r"vis-correction-sc.*\[not recorded\]"):
        calibrate_exposure(sci, cals, "GRAVITY_SC")
    with pytest.raises(ValueError, match="cannot be checked"):  # science side
        calibrate_exposure(cals[0], [sci], "GRAVITY_SC")


def test_gravity_vis_recipe_found_among_several_recs():
    from virgil_vlti.gravity import compatibility, pipeline_params

    def header(value, first="gravity_dualscivis"):
        return fits.Header(
            {
                "ESO PRO REC1 ID": first,
                "ESO PRO REC1 PARAM1 NAME": "vis-correction-sc",
                "ESO PRO REC1 PARAM1 VALUE": "decoy",
                "ESO PRO REC2 ID": "gravity_vis",
                "ESO PRO REC2 PARAM1 NAME": "vis-correction-sc",
                "ESO PRO REC2 PARAM1 VALUE": value,
                "ESO INS SPEC RES": "MEDIUM",
                "ESO INS POLA MODE": "COMBINED",
            }
        )

    assert pipeline_params(header("VFACTOR")) == {"vis-correction-sc": "VFACTOR"}
    _, reduction, _ = compatibility(header("VFACTOR"), header("NONE"))
    assert len(reduction) == 1 and "'VFACTOR'" in reduction[0] and "'NONE'" in reduction[0]
    assert compatibility(header("VFACTOR"), header("VFACTOR")) == ([], [], [])
    # no ID card: REC1 as before
    assert pipeline_params(
        fits.Header({"ESO PRO REC1 PARAM1 NAME": "a", "ESO PRO REC1 PARAM1 VALUE": 1})
    ) == {"a": 1}


def _estimates(n, sigma=0.01):
    from virgil_vlti.calibrate import Estimates

    rng = np.random.default_rng(1)
    return Estimates(
        times=MJD0 + 0.01 * np.arange(n),
        values=0.8 + 0.005 * rng.normal(size=(n, 4)),
        errors=np.full((n, 4), sigma),
    )


@pytest.mark.parametrize("mode", ["loo", "scatter"])
def test_one_calibrator_estimate_keeps_its_statistical_error(mode):
    from virgil_vlti.calibrate import propagated_error

    est = _estimates(1)
    t = MJD0 + 0.003
    with pytest.warns(UserWarning, match="One calibrator estimate"):
        tau = prediction_error(est, mode, t=t)
    np.testing.assert_allclose(tau, propagated_error(t, est))
    assert np.all(tau > 0)


@pytest.mark.parametrize("mode", ["loo", "scatter"])
def test_prediction_error_is_never_below_the_statistical_error(mode):
    from virgil_vlti.calibrate import propagated_error

    est = _estimates(4, sigma=0.05)
    t = MJD0 + 0.015
    tau = prediction_error(est, mode, t=t)
    assert np.all(tau >= propagated_error(t, est))
    assert np.all(tau >= prediction_error(est, mode))  # without t: systematic only
    assert np.all(prediction_error(est, "none", t=t) == 0)


def test_one_calibrator_calibration_has_a_transfer_function_error(tmp_path):
    paths = [
        product(tmp_path, "cal1", MJD0, cal=True),
        product(tmp_path, "sci", MJD0 + 0.02, cal=False),
    ]
    sci, cal = classify([read_exposure(p) for p in paths], dual_field=True)[::-1]
    with pytest.warns(UserWarning, match="One calibrator estimate"):
        loo = calibrate_exposure(sci, [cal], "GRAVITY_SC", error="loo")
    none = calibrate_exposure(sci, [cal], "GRAVITY_SC", error="none")
    good = ~loo["vis2"].flag
    assert np.all(loo["vis2"].error[good] > none["vis2"].error[good])
    with pytest.warns(UserWarning, match="One calibrator estimate"):
        modes = calibrate_exposure(sci, [cal], "GRAVITY_SC", error="loo", tf_modes=True)
    assert np.all(modes["vis2"].tf_mode == 0)  # no systematic to split
    assert np.all(modes["vis2"].error[~modes["vis2"].flag] > none["vis2"].error[~none["vis2"].flag])


@pytest.mark.parametrize("bad", [np.nan, 0.0, -1.0])
def test_non_finite_zero_or_negative_errors_are_flagged(tmp_path, bad):
    sci, cals = _with_cal_header(tmp_path, {})
    clean = calibrate_exposure(sci, cals, "GRAVITY_SC")
    sci.v2["GRAVITY_SC"].data["VIS2ERR"][0, 3] = bad
    out = calibrate_exposure(sci, cals, "GRAVITY_SC")
    n = out["vis2"].value.shape[-1]
    assert np.isnan(out["vis2"].value[0, 3]) and np.isfinite(clean["vis2"].value[0, 3])
    assert np.all(np.isfinite(out["vis2"].error[~out["vis2"].flag]))
    assert n == clean["vis2"].value.shape[-1]
