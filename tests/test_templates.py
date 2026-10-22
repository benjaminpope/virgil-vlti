"""Calibration templates: building, saving, shipping and laying out modes."""

import dataclasses
import re
import warnings
from importlib import resources

import numpy as np
import pytest

from virgil_vlti import classify, read_exposure
from virgil_vlti.calpca import calibrator_residuals
from virgil_vlti.cli import main
from virgil_vlti.templates import Template, build_template, shipped, template_modes

from .test_calibrate import MJD0, product
from .test_calpca import INS, SHAPE, WAVE, WIDTH, calibrator_nights

PILOT = "gravity_single_medium_combined_2023-2025"


@pytest.fixture(scope="module")
def synthetic_template(tmp_path_factory):
    res = calibrator_residuals(calibrator_nights(tmp_path_factory.mktemp("t")), INS)
    return build_template(res, "synthetic", {"resolution": "MEDIUM"}, 1, 0, folds=4)


def test_telescope_mode_gets_the_telescope_share(synthetic_template):
    t = synthetic_template
    share = t.provenance["vis2_telescope_share"][0]
    assert share > 0.85  # the injected mode is purely telescope-based
    shape = SHAPE / np.linalg.norm(SHAPE)
    tel = t.vis2_tel[0]
    assert abs(tel @ shape) / np.linalg.norm(tel) > 0.98
    # Per telescope on log|V|: ½ × WIDTH (on log V²) × the leave-one-out factor.
    expect = 0.5 * WIDTH * np.linalg.norm(SHAPE) * np.sqrt(1 + 1 / 3)
    assert np.linalg.norm(tel) == pytest.approx(expect, rel=0.3)
    assert np.linalg.norm(t.vis2_base[0]) < 0.5 * np.linalg.norm(tel)


def test_save_and_load(synthetic_template, tmp_path):
    path = synthetic_template.save(tmp_path / "t.npz")
    back = Template.load(str(path))
    np.testing.assert_array_equal(back.vis2_tel, synthetic_template.vis2_tel)
    assert back.config == synthetic_template.config and back.name == "synthetic"


def test_shipped_pilot_template():
    assert PILOT in shipped()
    t = Template.load(PILOT)
    assert t.config["resolution"] == "MEDIUM" and t.config["field"] == "SINGLE"
    assert t.vis2_tel.shape[0] == 4 and t.t3phi.shape[0] == 2
    lo, hi = t.provenance["valid_um"]
    assert 1.99 < lo < 2.01 and 2.44 < hi < 2.46
    inside = (t.wave * 1e6 >= lo) & (t.wave * 1e6 <= hi)
    assert np.all(t.vis2_tel[:, ~inside] == 0) and np.any(t.vis2_tel[:, inside] != 0)


def test_mode_layout(synthetic_template, tmp_path):
    (sci,) = classify([read_exposure(product(tmp_path, "sci", MJD0, cal=False))])
    modes = template_modes(synthetic_template, sci, INS)
    frames = 3
    assert len(modes) == frames * (4 + 6)  # one component: 4 telescopes, 6 baselines per frame
    v2 = sci.v2[INS].data
    tel = next(m["vis2"] for name, m in modes.items() if "-tel28@" in name and "#frame0" in name)
    rows = np.flatnonzero(np.any(tel != 0, 1))
    assert len(rows) == 3 and all(28 in v2["STA_INDEX"][r] for r in rows)


def test_mismatch_refused_and_range_warnings(synthetic_template, tmp_path):
    path = product(tmp_path, "sci", MJD0, cal=False, header={"ESO INS SPEC RES": "HIGH"})
    (sci,) = classify([read_exposure(path)])
    with pytest.raises(ValueError, match="another configuration"):
        template_modes(synthetic_template, sci, INS)
    with pytest.warns(UserWarning, match="another configuration"):
        template_modes(synthetic_template, sci, INS, strict=False)
    t = Template.load(PILOT)  # valid 2.0-2.45 µm; the synthetic grid starts at 2.0
    (sci,) = classify([read_exposure(product(tmp_path, "sci2", MJD0, cal=False))])
    sci.wave[INS] = np.linspace(1.95e-6, 2.4e-6, WAVE.size)  # starts below 2.0 µm
    with warnings.catch_warnings(record=True) as caught:  # it also differs in field
        warnings.simplefilter("always")
        template_modes(t, sci, INS, strict=False)
    assert any("outside the template" in str(w.message) for w in caught)


def test_cli_writes_template_modes(tmp_path):
    for k, t in enumerate((MJD0, MJD0 + 0.05)):
        product(tmp_path, f"cal{k}", t, cal=True)
    product(tmp_path, "sci", MJD0 + 0.02, cal=False)
    out = tmp_path / "out"
    # The synthetic products are dual field; the pilot template is single field.
    with pytest.raises(SystemExit, match="another configuration"):
        main([str(tmp_path), "--out", str(out), "--dual-field", "--template", PILOT])
    with warnings.catch_warnings(record=True) as caught:  # and its grid starts at 2.0 µm
        warnings.simplefilter("always")
        main([str(tmp_path), "--out", str(out), "--dual-field", "--template", PILOT,
              "--any-template-config"])  # fmt: skip
    assert any("another configuration" in str(w.message) for w in caught)
    (fits_path,) = out.glob("*/*.fits")
    modes = np.load(str(fits_path)[: -len(".fits")] + ".template_modes.npz")
    assert any("vis2-pc1-tel" in k for k in modes.files)
    assert any("t3phi-pc2-tri" in k for k in modes.files)
    assert WAVE.size in modes[modes.files[0]].shape


def test_no_modes_outside_the_valid_range(synthetic_template):
    t = dataclasses.replace(  # a copy: the fixture is shared
        synthetic_template, provenance={**synthetic_template.provenance, "valid_um": [2.1, 2.3]}
    )
    um = np.linspace(2.0, 2.4, 81)
    tel, base, phi = t.on_grid(um * 1e-6)
    outside = (um < 2.1) | (um > 2.3)
    assert np.all(tel[:, outside] == 0) and np.any(tel[:, ~outside] != 0)


def test_closure_phase_only_template(tmp_path):
    res = calibrator_residuals(calibrator_nights(tmp_path, n_nights=4), INS)
    t = build_template(res, "cp-only", {}, 0, 1, folds=4)
    assert t.vis2_tel.shape[0] == 0 and t.t3phi.shape[0] == 1


def test_undefined_telescope_fraction_is_refused(tmp_path):
    res = calibrator_residuals(calibrator_nights(tmp_path, n_nights=4), INS)
    v2 = res["vis2"]
    keep = np.flatnonzero(v2.meta["stations"] != v2.meta["stations"][0])  # one baseline gone
    res["vis2"] = v2.select(keep)
    with pytest.raises(ValueError, match="six-baseline"):
        build_template(res, "broken", {}, 1, 0, folds=4, clip=0)


def test_frames_by_time_survive_a_missing_row(synthetic_template, tmp_path):
    (sci,) = classify([read_exposure(product(tmp_path, "sci", MJD0, cal=False))])
    v2 = sci.v2[INS]
    v2.data = v2.data[np.arange(len(v2.data)) != 7]  # frame 1 loses one baseline
    modes = template_modes(synthetic_template, sci, INS)
    for name, m in modes.items():
        if "-tel" in name:
            rows = np.flatnonzero(np.any(m["vis2"] != 0, 1))
            times = v2.data["MJD"][rows]
            assert np.ptp(times) * 86400 < 60  # all from one frame


def test_cli_modes_match_windowed_binned_tables(tmp_path):
    from astropy.io import fits

    for k, t in enumerate((MJD0, MJD0 + 0.05)):
        product(tmp_path, f"cal{k}", t, cal=True)
    product(tmp_path, "sci", MJD0 + 0.02, cal=False)
    out = tmp_path / "out"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        main([str(tmp_path), "--out", str(out), "--dual-field", "--template", PILOT,
              "--any-template-config", "--window", "2.1", "2.35", "--bin", "3"])  # fmt: skip
    (fits_path,) = out.glob("*/*.fits")
    with fits.open(fits_path) as h:
        shape_v2 = h["OI_VIS2"].data["VIS2DATA"].shape
        shape_t3 = h["OI_T3"].data["T3PHI"].shape
    modes = np.load(str(fits_path)[: -len(".fits")] + ".template_modes.npz")
    for key in (k for k in modes.files if not k.startswith("__")):
        assert modes[key].shape == (shape_v2 if "/vis2/" in key else shape_t3)


def _pairs_json(telescope=0.06, baseline=0.0, offset=1.4):
    return {
        "n_pairs": 24,
        "v2": {
            "independent only": {"widths": {"white": 1.37}, "delta_lnl": 0.0},
            "+ baseline + telescope": {
                "widths": {"white": 0.92, "baseline": baseline, "telescope": telescope},
                "delta_lnl": 898.0,
            },
        },
        "cp": {"+ triangle offset": {"widths": {"white": 1.29, "offset": offset}, "delta_lnl": 40}},
    }


def test_gain_template_reproduces_the_pair_test_covariance(tmp_path):
    """Flat telescope gains g give each baseline V² a variance 2 g² V², as in the pair test."""
    from virgil_vlti.templates import gain_template_from_pairs

    sci = classify([read_exposure(product(tmp_path, "sci", MJD0, cal=False))], dual_field=True)[0]
    g, o = 0.06, 1.4
    t = gain_template_from_pairs(_pairs_json(g, 0.0, o), "gains", config={"field": "DUAL"})
    assert t.mismatches(sci.header) == []
    modes = template_modes(t, sci, "GRAVITY_SC")
    v2 = [m["vis2"] for k, m in modes.items() if "#frame0" in k and "vis2" in m]
    assert len(v2) == 4  # one per telescope; no baseline modes when that width is zero
    assert not any("-base" in k for k in modes)
    # dV²/V² = 2 d log|V|: the variance per V² row, over V², is Σ (2 mode)².
    var = np.sum([(2 * m) ** 2 for m in v2], 0)
    rows = [r for r in range(var.shape[0]) if var[r].any()]
    np.testing.assert_allclose(var[rows], 2 * g**2)
    cp = [m["t3phi"] for k, m in modes.items() if "#frame0" in k and "t3phi" in m]
    assert len(cp) == 4
    np.testing.assert_allclose(np.sum(cp, 0)[np.sum(cp, 0) != 0], o)
    assert t.provenance["pair_test"]["v2_white"] == 0.92


def test_gain_template_cli_from_pair_test_json(tmp_path):
    import json

    from virgil_vlti.templates_cli import main as template_main

    path = tmp_path / "pairs.json"
    path.write_text(json.dumps(_pairs_json(0.06, 0.04, 1.4)))
    template_main(["--from-pairs", str(path), "--name", "gains", "--cp-model", "none",
                   "--config", "field=DUAL", "--out", str(tmp_path)])  # fmt: skip
    t = Template.load(str(tmp_path / "gains.npz"))
    np.testing.assert_allclose(t.vis2_tel, 0.03)
    np.testing.assert_allclose(t.vis2_base, 0.02)
    assert t.t3phi.shape[0] == 0
    with pytest.raises(SystemExit):
        template_main(["--name", "x"])  # neither residuals nor a pair test


def test_gain_widths_must_be_physical():
    from virgil_vlti.templates import gain_template

    with pytest.raises(ValueError):
        gain_template("bad", telescope=-0.1)


NOWR = "gravity_single_medium_combined_2023-2025_nowr"
HIGH = "gravity_single_high_combined_2023-2025"
SINGLE_HEADER = {
    "ESO PRO CATG": "SINGLE_SCI_VIS",
    "ESO INS POLA MODE": "COMBINED",
    "ESO ISS CONF T1NAME": "AT1",
}


def test_shipped_wr_templates_load():
    assert {HIGH, NOWR} <= set(shipped())
    high, nowr = Template.load(HIGH), Template.load(NOWR)
    assert high.config["resolution"] == "HIGH" and nowr.config["resolution"] == "MEDIUM"
    assert high.wave.size == 1628 and nowr.wave.size == 233
    assert high.vis2_tel.shape == high.vis2_base.shape == (3, 1628) and high.t3phi.shape == (
        2,
        1628,
    )
    assert nowr.vis2_tel.shape == (4, 233) and nowr.t3phi.shape == (2, 233)
    assert high.config["telescopes"] == "AT" and nowr.config["telescopes"] == "mostly-AT"


def test_shipped_widths_match_the_stored_modes(tmp_path):
    for name in (HIGH, NOWR, PILOT):
        t = Template.load(name)
        # build_template: closure shapes are W c; V² tel = ½ W √(share/2) c and
        # base = ½ W √(1 − share) c, c orthonormal, so the pieces add back to W².
        w = np.asarray(t.provenance["t3phi_heldout_widths"])
        np.testing.assert_allclose(np.linalg.norm(t.t3phi, axis=1), w, rtol=1e-3)
        wv = np.asarray(t.provenance["vis2_heldout_widths"])
        total = 4 * (2 * np.sum(t.vis2_tel**2, 1) + np.sum(t.vis2_base**2, 1))
        np.testing.assert_allclose(total, wv**2, rtol=1e-3)
        back = Template.load(str(t.save(tmp_path / f"{name}.npz")))
        assert back.provenance["t3phi_heldout_widths"] == w.tolist()
    assert np.allclose(
        Template.load(HIGH).provenance["t3phi_heldout_widths"], (15.0, 7.0), atol=0.1
    )


def _principal_cosines(a, b):
    qa, qb = np.linalg.qr(a.T)[0], np.linalg.qr(b.T)[0]
    return np.linalg.svd(qa.T @ qb, compute_uv=False)


def test_nowr_rebuild_has_the_pilot_shapes():
    pilot, nowr = Template.load(PILOT), Template.load(NOWR)
    np.testing.assert_array_equal(pilot.wave, nowr.wave)
    assert _principal_cosines(pilot.t3phi, nowr.t3phi).min() > 0.99
    assert _principal_cosines(pilot.vis2_tel, nowr.vis2_tel).min() > 0.97


def test_mostly_at_and_missing_telescopes():
    pilot, high = Template.load(PILOT), Template.load(HIGH)
    assert pilot.config["telescopes"] == "mostly-AT"
    spec = {"ESO INS SPEC RES": "MEDIUM", "ESO INS POLA MODE": "COMBINED",
            "ESO PRO CATG": "SINGLE_SCI_VIS"}  # fmt: skip
    assert pilot.mismatches({**spec, "ESO ISS CONF T1NAME": "AT1"}) == []
    (problem,) = pilot.mismatches({**spec, "ESO ISS CONF T1NAME": "UT1"})
    assert problem.startswith("telescopes")
    assert pilot.mismatches(spec) != []  # unrecorded: cannot be shown to match
    # The HIGH template is AT only: a UT exposure is refused.
    hspec = {**spec, "ESO INS SPEC RES": "HIGH"}
    assert high.mismatches({**hspec, "ESO ISS CONF T1NAME": "AT1"}) == []
    (problem,) = high.mismatches({**hspec, "ESO ISS CONF T1NAME": "UT1"})
    assert problem.startswith("telescopes")
    # A template without the key does not compare telescopes.
    free = dataclasses.replace(high, config={"resolution": "HIGH"})
    assert free.mismatches({**hspec, "ESO ISS CONF T1NAME": "UT1"}) == []


@pytest.mark.parametrize("name", [HIGH, NOWR])
def test_template_modes_with_shipped_wr_templates(name, tmp_path):
    t = Template.load(name)
    res = t.config["resolution"]
    (sci,) = classify([read_exposure(product(tmp_path, "sci", MJD0, cal=False))])
    sci.header.update({**SINGLE_HEADER, "ESO INS SPEC RES": res})
    assert t.mismatches(sci.header) == []
    modes = template_modes(t, sci, INS)  # strict: raises on a mismatch
    frames, k_v2, k_cp = 3, t.vis2_tel.shape[0], t.t3phi.shape[0]
    assert len(modes) == frames * (k_v2 * 10 + k_cp * 4)
    tel = next(m["vis2"] for n, m in modes.items() if "-tel28@" in n and "pc1" in n)
    assert tel.shape[1] == WAVE.size and np.any(tel != 0)
    sci.header["ESO INS SPEC RES"] = "LOW"
    with pytest.raises(ValueError, match="resolution"):
        template_modes(t, sci, INS)


def test_mostly_at_template_matches_at_data(synthetic_template, tmp_path):
    import dataclasses

    (sci,) = classify([read_exposure(product(tmp_path, "sci", MJD0, cal=False))])
    sci.header["ESO ISS CONF T1NAME"] = "AT1"
    config = {**synthetic_template.config, "telescopes": "mostly-AT"}
    t = dataclasses.replace(synthetic_template, config=config)
    assert not any("telescopes" in p for p in t.mismatches(sci.header))
    sci.header["ESO ISS CONF T1NAME"] = "U1"
    assert any("telescopes" in p for p in t.mismatches(sci.header))


@pytest.mark.parametrize("name", shipped())
def test_shipped_template_metadata_has_no_private_paths(name):
    path = resources.files("virgil_vlti") / "templates" / f"{name}.npz"
    with np.load(path, allow_pickle=False) as data:
        meta = str(data["meta"])
    found = re.search(
        r"/Users/|/home/|/fred/|/scratch/|~/|[A-Za-z]:\\\\|ozstar|oz\d{3}|bpope",
        meta,
        re.IGNORECASE,
    )
    assert not found, found.group()
