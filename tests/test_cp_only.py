"""Closure-phase-only products, and screening without calibrator company."""

import numpy as np
import pytest
from astropy.io import fits

from virgil_vlti import calibrate_exposure, classify, read_exposure
from virgil_vlti.calibrate import closure_phase_only
from virgil_vlti.cli import main
from virgil_vlti.templates import gain_template_from_pairs

from .test_calibrate import MJD0, product
from .test_calpca import INS
from .test_load import assert_same
from .test_templates import _pairs_json


def zero_closure_calibrators(tmp_path):
    """Calibrators whose closure phases are zero, so that calibrating leaves T3PHI as observed."""
    cals = []
    for name, t in (("cal1", MJD0), ("cal2", MJD0 + 0.05)):
        path = product(tmp_path, name, t, cal=True)
        with fits.open(path, mode="update") as h:
            for hdu in h:
                if hdu.name == "OI_T3":
                    hdu.data["T3PHI"][:] = 0.0
        cals.append(path)
    return cals


def science_night(tmp_path, name="sci", t=MJD0 + 0.02):
    return product(tmp_path, name, t, cal=False)


def read(path, ext="OI_VIS2"):
    with fits.open(path) as h:
        hdu = next(x for x in h if x.name == ext and x.header["INSNAME"] == INS)
        return h[0].header.copy(), hdu.data.copy()


@pytest.fixture
def uncalibrated(tmp_path):
    """A folder with one science file and no calibrator, and the output folder's CLI arguments."""
    folder = tmp_path / "data"
    folder.mkdir()
    science_night(folder)
    template = gain_template_from_pairs(
        _pairs_json(0.06, 0.0, 1.4), "gains", config={"field": "DUAL"}
    )
    path = tmp_path / "template.npz"
    template.save(path)
    return folder, tmp_path / "out", str(path)


def run(folder, out, *extra):
    main([str(folder), "--out", str(out), "--dual-field", "--bin", "2",
          "--window", "2.1", "2.3", *extra])  # fmt: skip
    (path,) = out.rglob("*.fits")
    return path


def test_no_calibrator_is_skipped_unless_asked(uncalibrated):
    folder, out, _ = uncalibrated
    with pytest.raises(SystemExit, match="Nothing calibrated"):
        main([str(folder), "--out", str(out), "--dual-field"])


def test_closure_phase_only_output_has_no_usable_v2(uncalibrated):
    folder, out, template = uncalibrated
    path = run(folder, out, "--cp-only-when-uncalibrated", "--template", template)
    header, v2 = read(path)
    assert header["CPONLY"] is True
    assert v2["FLAG"].all() and np.isnan(v2["VIS2DATA"]).all() and np.isnan(v2["VIS2ERR"]).all()
    _, t3 = read(path, "OI_T3")
    assert np.isnan(t3["T3AMP"]).all()
    assert (~t3["FLAG"]).any() and np.isfinite(t3["T3PHI"][~t3["FLAG"]]).all()
    modes = np.load(str(path)[: -len(".fits")] + ".template_modes.npz")
    names = [k for k in modes.files if not k.startswith("__")]
    kinds = {key.split("/")[1] for key in names}
    assert kinds == {"t3phi"}
    assert all(modes[k].shape == t3["T3PHI"].shape for k in names)


def test_closure_phase_binning_matches_the_calibrated_path(tmp_path):
    """Window, binning and flags only: the calibrators here have zero closure phase, so
    this does not test calibration (the science file's instrumental offset stays in both)."""
    cals = zero_closure_calibrators(tmp_path)
    sci_path = science_night(tmp_path)
    sci, *cal = classify([read_exposure(p) for p in (sci_path, *cals)], dual_field=True)
    kw = {"window_um": (2.1, 2.3), "nbin": 3}
    calibrated = calibrate_exposure(sci, cal, INS, error="none", **kw)
    with pytest.warns(UserWarning, match="instrumental closure phase"):
        only = closure_phase_only(sci, INS, **kw)
    assert only.cp_only and only.n_cal == 0
    np.testing.assert_array_equal(only.wave, calibrated.wave)
    np.testing.assert_array_equal(only["t3phi"].flag, calibrated["t3phi"].flag)
    ok = ~only["t3phi"].flag
    np.testing.assert_allclose(only["t3phi"].value[ok], calibrated["t3phi"].value[ok], atol=1e-9)
    np.testing.assert_allclose(only["t3phi"].error[ok], calibrated["t3phi"].error[ok], rtol=1e-9)
    assert only["vis2"].flag.all() and only["t3amp"].flag.all()


def test_cp_only_marks_a_file_that_has_calibrators(tmp_path):
    folder = tmp_path / "data"
    folder.mkdir()
    for p in zero_closure_calibrators(folder):
        assert p
    science_night(folder, "suspect_sci")
    science_night(folder, "good_sci", MJD0 + 0.03)
    out = tmp_path / "out"
    main([str(folder), "--out", str(out), "--dual-field", "--cp-only", "suspect"])
    flagged = {p.name: read(p)[0].get("CPONLY", False) for p in out.rglob("*.fits")}
    assert sorted(flagged.values()) == [False, True]


def _needs_new_virgil():
    oidata = pytest.importorskip("virgil.oidata")
    if not hasattr(oidata, "_valid_samples"):
        pytest.skip("needs a virgil newer than 0.3.0")


def test_loader_reads_closure_phase_only_files(uncalibrated):
    _needs_new_virgil()
    from virgil.oidata import OIData
    from virgil.oifits import read_oifits

    from virgil_vlti import to_oidata

    folder, out, template = uncalibrated
    path = run(folder, out, "--cp-only-when-uncalibrated", "--template", template)
    _, t3 = read(path, "OI_T3")
    data = to_oidata(path, INS)
    assert data.phase_offsets is not None
    assert data.gains is None
    plain = read_oifits(str(path), insname=INS)
    assert int(np.sum(np.isfinite(np.asarray(plain["vis"])))) == 0
    n_cp = int(np.sum(~t3["FLAG"] & np.isfinite(t3["T3PHI"])))
    assert np.size(data.phi) == n_cp
    # the template's closure modes, aligned cell for cell with the kept closure phases
    modes = np.load(str(path)[: -len(".fits")] + ".template_modes.npz")
    cells = np.flatnonzero(~t3["FLAG"].reshape(-1) & np.isfinite(t3["T3PHI"].reshape(-1)))
    expected = np.deg2rad(
        [
            modes[k].reshape(-1)[cells]
            for k in sorted(k for k in modes.files if not k.startswith("__"))
        ]
    )
    expected = expected[np.any(expected != 0, axis=1)]
    bare = OIData(plain).with_closure_offsets(modes=expected)
    assert_same(data.phase_offsets, bare.phase_offsets)


def test_loader_drops_v2_modes_beside_a_closure_phase_only_file(uncalibrated):
    _needs_new_virgil()
    from virgil_vlti import to_oidata
    from virgil_vlti.write import write_modes

    folder, out, template = uncalibrated
    path = run(folder, out, "--cp-only-when-uncalibrated", "--template", template)
    v2 = read(path)[1]
    write_modes(
        str(path)[: -len(".fits")] + ".modes.npz",
        {INS: {"rel-STAR@night": {"vis2": np.ones(v2["VIS2DATA"].shape)}}},
    )
    data = to_oidata(path, INS)
    assert data.gains is None and data.phase_offsets is not None


def test_screening_without_company_is_not_screened_and_never_raises(tmp_path):
    from virgil_vlti.calpca import screen_calibrators

    from .test_calpca import calibrator_nights

    exposures = calibrator_nights(tmp_path, n_nights=1)[:1]  # one calibrator star: no company
    with pytest.warns(UserWarning, match="not screened"):
        res, excluded = screen_calibrators(exposures, INS)
    assert res == {} and excluded == []
    with pytest.warns(UserWarning, match="not screened"):  # no calibrator product at all
        assert screen_calibrators([], INS) == ({}, [])


def test_closure_phase_only_without_calibrator_keeps_the_instrumental_offset(tmp_path):
    """Against truth: the uncalibrated T3PHI is the science phase plus the instrumental one."""
    from virgil_vlti.calibrate import bin_channels, channel_window

    from .test_calibrate import WAVE, _truth_cp, instrumental_cp

    (sci,) = classify([read_exposure(science_night(tmp_path))], dual_field=True)
    with pytest.warns(UserWarning, match="instrumental closure phase"):
        only = closure_phase_only(sci, INS, window_um=(2.1, 2.3), nbin=2)
    keep = channel_window(sci.wave[INS], (2.1, 2.3))
    truth = np.mean(bin_channels(_truth_cp(sci) + instrumental_cp(WAVE), keep, 2), -1)
    bare = np.mean(bin_channels(_truth_cp(sci), keep, 2), -1)
    ok = ~only["t3phi"].flag
    np.testing.assert_allclose(only["t3phi"].value[ok], truth[ok], atol=0.5)
    assert np.abs(only["t3phi"].value[ok] - bare[ok]).mean() > 2.0  # the offset is still in


def test_cp_only_on_a_file_with_calibrators_calibrates_its_closure_phases(tmp_path):
    folder = tmp_path / "data"
    folder.mkdir()
    cals = zero_closure_calibrators(folder)
    sci_path = science_night(folder, "suspect_sci")
    out = tmp_path / "out"
    main([str(folder), "--out", str(out), "--dual-field", "--cp-only", "suspect"])
    (path,) = out.rglob("*.fits")
    header, v2 = read(path)
    assert header["CPONLY"] is True and v2["FLAG"].all()
    sci, *cal = classify([read_exposure(p) for p in (sci_path, *cals)], dual_field=True)
    expected = calibrate_exposure(sci, cal, INS)["t3phi"]
    _, t3 = read(path, "OI_T3")
    np.testing.assert_allclose(t3["T3PHIERR"], expected.error)  # transfer-function error kept


def test_cp_only_name_that_marks_nothing_is_refused(uncalibrated):
    folder, out, _ = uncalibrated
    with pytest.raises(SystemExit, match="marked no science file"):
        main([str(folder), "--out", str(out), "--dual-field", "--cp-only", "no_such_star"])


def test_cp_only_matches_the_output_stem(uncalibrated):
    folder, out, _ = uncalibrated
    main([str(folder), "--out", str(out), "--dual-field", "--cp-only", "scistar_20240331T025248"])
    (path,) = out.rglob("*.fits")
    assert read(path)[0]["CPONLY"] is True
