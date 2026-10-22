"""Dual-field calibration templates: field configurations, selection and refusal."""

import dataclasses
import json
import warnings

import numpy as np
import pytest
from astropy.io import fits

from virgil_vlti import classify, read_exposure
from virgil_vlti.calpca import calibrator_residuals
from virgil_vlti.calpca_cli import main as calpca_main
from virgil_vlti.cli import main as calibrate_main
from virgil_vlti.gravity import field_matches, field_mode
from virgil_vlti.pairs_cli import main as pairs_main
from virgil_vlti.templates import build_template, template_modes, template_name
from virgil_vlti.templates_cli import main as template_main

from .test_calibrate import MJD0, product
from .test_calpca import INS, calibrator_nights

OFF = (700.0, -500.0)  # mas: the science fibre on another star in the field


def header(**cards):
    return fits.Header({f"HIERARCH {k}" if " " in k else k: v for k, v in cards.items()})


@pytest.mark.parametrize(
    ("cards", "mode"),
    [
        ({"ESO PRO CATG": "SINGLE_CAL_VIS", "ESO INS SOBJ X": 700.0, "ESO INS SOBJ Y": 0.0},
         "SINGLE"),
        ({"ESO PRO CATG": "DUAL_CAL_VIS", "ESO INS SOBJ X": 0.01, "ESO INS SOBJ Y": 0.01},
         "DUAL_ON"),
        ({"ESO PRO CATG": "DUAL_CAL_VIS", "ESO INS SOBJ X": 700.0, "ESO INS SOBJ Y": -500.0},
         "DUAL_OFF"),
        ({"ESO PRO CATG": "DUAL_SCI_VIS", "ESO INS SOBJ X": 9000.0, "ESO INS SOBJ Y": 0.0},
         "DUAL_WIDE"),
        ({"ESO PRO CATG": "DUAL_CAL_VIS"}, "DUAL"),  # no fibre offset recorded
        ({"INSMODE": "SINGLE,MEDIUM,COMBINED,COMBINED"}, "SINGLE"),
        ({"ESO QC FIELD MODE": "DUAL", "ESO PRO CATG": "X", "ESO INS SOBJ X": 0.0,
          "ESO INS SOBJ Y": 0.0}, "DUAL_ON"),
        ({}, None),
    ],
)  # fmt: skip
def test_field_mode(cards, mode):
    assert field_mode(header(**cards)) == mode


def test_field_matches():
    assert field_matches("DUAL", "DUAL_OFF") and field_matches("dual_off", "DUAL_OFF")
    assert not field_matches("DUAL_OFF", "DUAL_ON")
    assert not field_matches("DUAL_OFF", "DUAL")  # geometry unknown: cannot match
    assert not field_matches("SINGLE", None)
    with pytest.raises(ValueError, match="Unknown field"):
        field_matches("TRIPLE", "SINGLE")


def test_template_name():
    config = {"field": "DUAL_OFF", "resolution": "MEDIUM", "polarisation": "COMBINED"}
    assert template_name(config, "2023-2025") == "gravity_dual-off_medium_combined_2023-2025"
    single = {**config, "field": "SINGLE"}
    assert template_name(single, "2023-2025") == "gravity_single_medium_combined_2023-2025"
    with pytest.raises(ValueError, match="polarisation"):
        template_name({"field": "DUAL", "resolution": "MEDIUM"}, "2024")


@pytest.fixture(scope="module")
def mixed(tmp_path_factory):
    """Single-field and off-axis dual-field calibrators on the same nights, setup and stations."""
    path = tmp_path_factory.mktemp("mixed")
    single = calibrator_nights(path, n_nights=4, per_night=3, prefix="s")
    dual = calibrator_nights(path, n_nights=4, per_night=3, seed=1, catg="DUAL_CAL_VIS",
                             sobj=OFF, prefix="d")  # fmt: skip
    return path, single, dual


def test_field_selection_and_no_mixing(mixed):
    _, single, dual = mixed
    res = calibrator_residuals(single + dual, INS)
    v2 = res["vis2"]
    assert set(v2.meta["field"]) == {"SINGLE", "DUAL_OFF"}
    assert np.all(v2.meta["n_cal"] == 2)  # each calibrated by its own field's two others
    off = calibrator_residuals(single + dual, INS, field="DUAL_OFF")["vis2"]
    assert set(off.meta["field"]) == {"DUAL_OFF"} and off.r.shape[0] == v2.r.shape[0] // 2
    assert np.allclose(off.meta["separation_mas"], np.hypot(*OFF))
    assert set(calibrator_residuals(single + dual, INS, field="DUAL")["vis2"].meta["field"]) == {
        "DUAL_OFF"
    }
    with pytest.raises(ValueError, match="in field DUAL_ON"):
        calibrator_residuals(single + dual, INS, field="DUAL_ON")


def test_same_star_references(tmp_path):
    # Dual-field programmes usually observe one calibrator star a night.
    cals = calibrator_nights(tmp_path, n_nights=3, per_night=3, catg="DUAL_CAL_VIS", sobj=OFF,
                             stars_per_night=1)  # fmt: skip
    with pytest.raises(ValueError, match="two calibrators of different stars"):
        calibrator_residuals(cals, INS)
    res = calibrator_residuals(cals, INS, references="same-star")["vis2"]
    assert res.r.shape[0] == 3 * 3 * 3 * 6 and np.all(res.meta["n_cal"] == 2)
    assert calibrator_residuals(cals, INS, references="any")["vis2"].r.shape == res.r.shape
    with pytest.raises(ValueError, match="references"):
        calibrator_residuals(cals, INS, references="nearest")


def test_calpca_cli_records_field_and_references(mixed, tmp_path):
    path, _, _ = mixed
    out = tmp_path / "pca"
    calpca_main([str(path), "--out", str(out), "--field", "dual_off", "--references", "any",
                 "--k-max", "1", "--folds", "2", "--screen-factor", "0",
                 "--clip", "0"])  # fmt: skip
    summary = json.loads((out / "summary.json").read_text())
    assert summary["field"] == "DUAL_OFF" and summary["references"] == "any"
    assert summary["vis2"]["files"] == 12
    res = np.load(out / "residuals.npz")
    assert set(res["vis2_meta_field"]) == {"DUAL_OFF"}


@pytest.fixture(scope="module")
def dual_template(mixed):
    _, _, dual = mixed
    res = calibrator_residuals(dual, INS)
    config = {"field": "DUAL_OFF", "resolution": "MEDIUM", "polarisation": "COMBINED"}
    return build_template(res, "dual", config, 1, 0, folds=2)


def test_template_refuses_other_fields(dual_template, tmp_path):
    off = classify([read_exposure(product(tmp_path, "off", MJD0, cal=False, sobj=OFF))])[0]
    assert dual_template.mismatches(off.header) == []
    assert template_modes(dual_template, off, INS)  # matching: modes, no warning
    on = classify([read_exposure(product(tmp_path, "on", MJD0, cal=True))])[0]
    with pytest.raises(ValueError, match="field: data 'DUAL_ON', template 'DUAL_OFF'"):
        template_modes(dual_template, on, INS)
    with pytest.warns(UserWarning, match="another configuration"):
        template_modes(dual_template, on, INS, strict=False)
    single = product(tmp_path, "single", MJD0, cal=False, header={"ESO PRO CATG": "SINGLE_SCI_VIS"})
    (single,) = classify([read_exposure(single)])
    with pytest.raises(ValueError, match="field: data 'SINGLE'"):
        template_modes(dual_template, single, INS)
    split = product(tmp_path, "split", MJD0, cal=False, sobj=OFF,
                    header={"ESO INS POLA MODE": "SPLIT"})  # fmt: skip
    with pytest.raises(ValueError, match="polarisation"):
        template_modes(dual_template, classify([read_exposure(split)])[0], INS)


def test_template_cli_names_from_config(mixed, tmp_path):
    path, _, _ = mixed
    pca = tmp_path / "pca"
    calpca_main([str(path), "--out", str(pca), "--field", "DUAL_OFF", "--k-max", "1",
                 "--folds", "2", "--screen-factor", "0"])  # fmt: skip
    config = ["field=dual_off", "resolution=MEDIUM", "polarisation=COMBINED"]
    template_main([str(pca / "residuals.npz"), "--k-vis2", "1", "--folds", "2", "--epoch",
                   "2023-2025", "--config", *config, "--out", str(tmp_path)])  # fmt: skip
    assert (tmp_path / "gravity_dual-off_medium_combined_2023-2025.npz").exists()
    with pytest.raises(SystemExit):  # an unknown field
        template_main([str(pca / "residuals.npz"), "--k-vis2", "1", "--epoch", "2024",
                       "--config", "field=WIDE", "resolution=MEDIUM",
                       "polarisation=COMBINED"])  # fmt: skip


def test_calibrate_cli_refuses_a_template_of_another_field(dual_template, tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    for k, t in enumerate((MJD0, MJD0 + 0.05)):
        product(data, f"cal{k}", t, cal=True)
    product(data, "sci", MJD0 + 0.02, cal=False)  # off axis: DUAL_OFF
    config = {**dual_template.config, "field": "SINGLE"}
    single = dataclasses.replace(dual_template, name="single", config=config)
    path = str(single.save(tmp_path / "single.npz"))
    with pytest.raises(SystemExit, match="another configuration"):
        calibrate_main([str(data), "--out", str(tmp_path / "out"), "--dual-field",
                        "--template", path])  # fmt: skip
    assert not list((tmp_path / "out").glob("*/*.fits"))  # refused before writing anything
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        calibrate_main([str(data), "--out", str(tmp_path / "out"), "--dual-field",
                        "--template", path, "--any-template-config"])  # fmt: skip
    assert list((tmp_path / "out").glob("*/*.template_modes.npz"))


def test_pairs_cli_by_star_night(tmp_path, capsys):
    calibrator_nights(tmp_path, n_nights=2, per_night=3, catg="DUAL_CAL_VIS", sobj=OFF,
                      stars_per_night=1)  # fmt: skip
    calibrator_nights(tmp_path, n_nights=2, per_night=2, prefix="s")  # single field
    out = tmp_path / "pairs.json"
    pairs_main([str(tmp_path), "--field", "DUAL_OFF", "--by-star-night", "--json", str(out)])
    result = json.loads(out.read_text())
    assert len(result) == 2 and all(r["field"] == ["DUAL_OFF"] for r in result)
    assert all(len(r["files"]) == 3 and r["n_pairs"] >= 1 for r in result)
    assert "telescope" in result[0]["v2"]["+ baseline + telescope"]["widths"]
