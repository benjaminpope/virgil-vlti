"""The gravity_viscal port: the pipeline's own rules, on synthetic products."""

import numpy as np
import pytest

from virgil_vlti import classify, read_exposure
from virgil_vlti.calibrate import fit_channel_polynomial, uniform_disk_v2
from virgil_vlti.cli import main
from virgil_vlti.viscal import calibrate_like_viscal, setup_string

from .test_calibrate import MJD0, _truth_v2, product

INS = "GRAVITY_SC"


def _night(tmp_path, cal_header=None, n_cal=1):
    # One frame per file, so each calibrator row is that baseline's only TF.
    paths = [product(tmp_path, "sci", MJD0 + 0.02, cal=False, nframes=1)]
    paths += [
        product(tmp_path, f"cal{k}", MJD0 + 0.04 * k, cal=True, header=cal_header, nframes=1)
        for k in range(n_cal)
    ]
    exposures = classify([read_exposure(p) for p in paths])
    return exposures[0], exposures[1:]


def test_fit_matches_the_pipeline_normal_equations():
    # gravi_vis_fit_amp solves the normal equations of rows scaled by w = σ⁻²
    # (not σ⁻¹), in the channel index.
    rng = np.random.default_rng(5)
    n, deg = 60, 3
    values = 1 + 0.1 * np.sin(np.arange(n) / 9) + rng.normal(0, 0.02, n)
    errors = rng.uniform(0.01, 0.05, n)
    flags = np.zeros(n, bool)
    flags[[3, 40]] = True
    x = np.arange(n, dtype=float)
    w = np.where(flags, 0, errors**-2.0)
    A = np.vander(x, deg + 1, increasing=True) * w[:, None]
    coef = np.linalg.solve(A.T @ A, A.T @ (np.where(flags, 0, values) * w))
    expect = np.vander(x, deg + 1, increasing=True) @ coef
    np.testing.assert_allclose(
        fit_channel_polynomial(values, errors, flags, deg), expect, rtol=1e-9
    )


def test_point_source_without_smoothing_is_exact(tmp_path):
    sci, cals = _night(tmp_path)
    out = calibrate_like_viscal(sci, cals, INS, smoothing=False)
    # One calibrator: the TF is that calibrator's, with no interpolation.
    cal = cals[0].v2[INS].data["VIS2DATA"]
    np.testing.assert_allclose(out["vis2"].value, sci.v2[INS].data["VIS2DATA"] / cal, rtol=1e-12)
    np.testing.assert_allclose(out["vis2"].error, sci.v2[INS].data["VIS2ERR"] / cal, rtol=1e-12)


def test_t3amp_is_not_corrected_for_the_calibrator_size(tmp_path):
    paths = [
        product(tmp_path, "sci", MJD0 + 0.02, cal=False, nframes=1),
        product(tmp_path, "cal", MJD0, cal=True, diameter=1.2, nframes=1),
    ]
    sci, cal = classify([read_exposure(p) for p in paths])
    out = calibrate_like_viscal(sci, [cal], INS, {cal.target: 1.2}, smoothing=False)
    d = cal.v2[INS]
    tf = d.data["VIS2DATA"] / uniform_disk_v2(d, cal.wave[INS], 1.2)  # the file's float32 grid
    np.testing.assert_allclose(out["vis2"].value, sci.v2[INS].data["VIS2DATA"] / tf, rtol=1e-12)
    raw = sci.t3[INS].data["T3AMP"] / cal.t3[INS].data["T3AMP"]
    np.testing.assert_allclose(out["t3amp"].value, raw, rtol=1e-12)


def test_weights_are_exponential_in_time(tmp_path):
    sci, cals = _night(tmp_path, n_cal=2)
    times = (100.0, 2000.0)
    for cal, t in zip(cals, times, strict=True):
        cal.v2[INS].data["TIME"][:] = t
        cal.v2[INS].data["VIS2DATA"][:] = 0.5 if t == 100.0 else 0.8
        cal.v2[INS].data["VIS2ERR"][:] = 0.01
    sci.v2[INS].data["TIME"][:] = 600.0
    out = calibrate_like_viscal(sci, cals, INS, smoothing=False, delta_time=3600.0)
    w = np.exp(-2 * np.abs(600.0 - np.array(times)) / 3600.0) / 0.016**2  # 2% of the TF floors σ
    w[1] = np.exp(-2 * 1400 / 3600.0) / max(0.01, 0.02 * 0.8) ** 2
    w[0] = np.exp(-2 * 500 / 3600.0) / max(0.01, 0.02 * 0.5) ** 2
    tf = (w[0] * 0.5 + w[1] * 0.8) / w.sum()
    np.testing.assert_allclose(out["vis2"].value, sci.v2[INS].data["VIS2DATA"] / tf, rtol=1e-12)


def test_setup_mismatch_dropped_unless_forced(tmp_path):
    sci, cals = _night(tmp_path, cal_header={"ESO DET2 SEQ1 DIT": 10.0})
    assert setup_string(sci.header) != setup_string(cals[0].header)
    with pytest.raises(ValueError, match="setup"):
        calibrate_like_viscal(sci, cals, INS)
    assert calibrate_like_viscal(sci, cals, INS, force_calib=True).n_cal == 1


def test_smoothing_follows_a_smooth_truth(tmp_path):
    rng = np.random.default_rng(6)
    sci, cals = _night(tmp_path)
    d = cals[0].v2[INS].data
    d["VIS2DATA"] *= 1 + rng.normal(0, 0.03, d["VIS2DATA"].shape)
    raw = calibrate_like_viscal(sci, cals, INS, smoothing=False)["vis2"].value
    smooth = calibrate_like_viscal(sci, cals, INS)["vis2"].value
    truth = _truth_v2(sci)
    assert np.std(smooth / truth - 1) < 0.5 * np.std(raw / truth - 1)


def test_unknown_parameters_are_refused(tmp_path):
    sci, cals = _night(tmp_path)
    with pytest.raises(TypeError, match="Unknown"):
        calibrate_like_viscal(sci, cals, INS, maxdegree=3)


def test_cli_parity_mode(tmp_path):
    _night(tmp_path)
    out = tmp_path / "out"
    with pytest.raises(SystemExit):
        main([str(tmp_path), "--out", str(out), "--match-gravity-viscal", "--bin", "2"])
    main([str(tmp_path), "--out", str(out), "--match-gravity-viscal"])
    assert len(list(out.glob("*/*.fits"))) == 1
    for bad in (["--smooth-tf", "-1"], ["--maxdeg-tfvis", "x"]):
        with pytest.raises(SystemExit):
            main([str(tmp_path), "--out", str(tmp_path / "o2"), *bad])


def test_smooth_running_window_matches_pipeline():
    """``gravi_vis_smooth_amp`` averages channels k − nsamp … k + nsamp − 1 (exclusive bound)."""
    from virgil_vlti.viscal import _smooth_running

    nwave, nsamp = 20, 2
    values = np.zeros(nwave)
    values[10] = 1.0  # an impulse at channel 10
    errors = np.ones(nwave)
    flags = np.zeros(nwave, bool)
    out, _ = _smooth_running(values, errors, flags, nsamp, phase=False)
    # Channels whose window k−2 … k+1 contains channel 10: k = 9 … 12, not 8.
    assert np.flatnonzero(out).tolist() == [9, 10, 11, 12]
    np.testing.assert_allclose(out[9:13], 1 / (2 * nsamp))
