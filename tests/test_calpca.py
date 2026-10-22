"""The calibrator PCA on synthetic calibrators with a known systematic."""

import json

import numpy as np
import pytest
from astropy.io import fits

from virgil_vlti import read_exposure
from virgil_vlti.calibrate import uniform_disk_amplitude
from virgil_vlti.calpca import (
    calibrator_residuals,
    cross_validate,
    heldout_widths,
    telescope_fraction,
    weighted_pca,
)
from virgil_vlti.calpca_cli import main
from virgil_vlti.synthetic import STATIONS, write_product

INS = "GRAVITY_SC"
WAVE = np.linspace(2.0e-6, 2.4e-6, 40)
SHAPE = (2.2e-6 / WAVE) ** 2 - 1  # a coherence-loss-like chromatic shape
WIDTH, NOISE = 0.03, 0.004  # per-telescope rms amplitude; fractional V² noise


def calibrator_nights(
    tmp_path,
    n_nights=8,
    per_night=4,
    seed=0,
    width=WIDTH,
    *,
    catg="SINGLE_CAL_VIS",
    sobj=(0.0, 0.0),
    stars_per_night=None,
    diameter=0.0,
    prefix="",
    disk=False,
):
    """Point-source calibrators whose V² carry a telescope-based chromatic mode.

    ``stars_per_night`` (default ``per_night``: every exposure another star)
    sets how many distinct stars share a night's exposures; ``diameter`` is
    the catalogue diameter the products record. By default the stars are
    point sources; with ``disk=True`` their true V² is the uniform disk of
    that diameter, and the noise is additive, matching the recorded VIS2ERR.
    """

    def truth(u, v, w, t):
        if not disk:
            return np.full_like(w, 0.8)
        return uniform_disk_amplitude(np.atleast_1d(u), np.atleast_1d(v), w, diameter)[0] ** 2

    rng = np.random.default_rng(seed)
    stars = per_night if stars_per_night is None else stars_per_night
    out = []
    for night in range(n_nights):
        for k in range(per_night):
            mjd = 60400.1 + night + k * 0.03
            path = write_product(
                tmp_path / f"{prefix}n{night}c{k}.fits",
                mjd=mjd,
                wave=WAVE,
                vis2=truth,
                t3phi=lambda *a: np.zeros_like(a[4]),
                catg=catg,
                sobj=sobj,
                target=f"{prefix}STAR{night}_{k % stars}",
                vis2err=0.8 * NOISE,
                t3phierr=0.5,
                header={"ESO QC TF CAL_SC DIAM": diameter, "ESO QC TF CAL_SC DIAM ERR": 0.01},
            )
            e = read_exposure(str(path))
            d = e.v2[INS].data
            seen = {}
            gains = {}
            for row, sta in enumerate(map(tuple, d["STA_INDEX"])):
                frame = seen.get(sta, 0)
                seen[sta] = frame + 1
                if frame not in gains:
                    gains[frame] = dict(zip(STATIONS, rng.normal(0, width, 4), strict=True))
                a = gains[frame][sta[0]] + gains[frame][sta[1]]
                if disk:
                    d["VIS2DATA"][row] = d["VIS2DATA"][row] * np.exp(a * SHAPE) + rng.normal(
                        0, 0.8 * NOISE, WAVE.size
                    )
                else:
                    d["VIS2DATA"][row] *= np.exp(a * SHAPE) * (1 + rng.normal(0, NOISE, WAVE.size))
            t3 = e.t3[INS].data
            t3["T3PHI"] += rng.normal(0, 0.5, t3["T3PHI"].shape)
            with fits.open(path, mode="update") as h:  # for readers of the file too
                for hdu in h:
                    if hdu.header.get("INSNAME") == INS and hdu.name == "OI_VIS2":
                        hdu.data["VIS2DATA"][:] = d["VIS2DATA"]
                    if hdu.header.get("INSNAME") == INS and hdu.name == "OI_T3":
                        hdu.data["T3PHI"][:] = t3["T3PHI"]
            out.append(e)
    return out


@pytest.fixture(scope="module")
def residuals(tmp_path_factory):
    return calibrator_residuals(calibrator_nights(tmp_path_factory.mktemp("cal")), INS)


def test_residual_rows(residuals):
    v2, cp = residuals["vis2"], residuals["t3phi"]
    assert v2.r.shape == (8 * 4 * 3 * 6, WAVE.size)  # nights × files × frames × baselines
    assert cp.r.shape[0] == 8 * 4 * 3 * 4
    assert set(np.unique(v2.meta["night"])) and v2.valid.all()


def test_pca_finds_the_injected_shape_and_width(residuals):
    v2 = residuals["vis2"]
    pca = weighted_pca(v2, 2)
    shape = SHAPE / np.linalg.norm(SHAPE)
    assert abs(pca.components[0] @ shape) > 0.98
    # Per-baseline amplitude: two telescopes, each N(0, WIDTH), minus what the
    # other calibrators' average shares; on its unit-norm component the rms is
    # ~ sqrt(2)·WIDTH·|SHAPE| times (1 + 1/n_others)^½.
    expect = np.sqrt(2) * WIDTH * np.linalg.norm(SHAPE) * np.sqrt(1 + 1 / 3)
    assert pca.widths[0] == pytest.approx(expect, rel=0.25)
    assert pca.widths[1] < 0.3 * pca.widths[0]


def test_cross_validation_keeps_one_component(residuals):
    cv = cross_validate(residuals["vis2"], 3, folds=4)
    chi2 = [c for _, c, _ in cv]
    noise = [e for _, _, e in cv]
    assert chi2[1] < 0.5 * chi2[0]  # the real component
    assert chi2[2] / chi2[1] > noise[2] / noise[1] - 0.02  # the next is only noise


def test_the_mode_is_telescope_based(residuals):
    v2 = residuals["vis2"]
    pca = weighted_pca(v2, 1)
    assert telescope_fraction(v2, pca.amplitudes)[0] > 0.9


def test_no_systematic_no_component(tmp_path):
    res = calibrator_residuals(calibrator_nights(tmp_path, n_nights=4, width=0.0), INS)
    cv = cross_validate(res["vis2"], 1, folds=4)
    assert cv[1][1] > cv[0][1] * (cv[1][2] - 0.02)  # no better than fitting noise


def test_heldout_widths_are_unbiased(residuals, tmp_path):
    expect = np.sqrt(2) * WIDTH * np.linalg.norm(SHAPE) * np.sqrt(1 + 1 / 3)
    widths = heldout_widths(residuals["vis2"], 2, folds=4)
    assert widths[0] == pytest.approx(expect, rel=0.25) and widths[1] < 0.1 * widths[0]
    # Closure phases carry no systematic here. With ~380 rows the held-out
    # variance ratio scatters by ~0.09, so a noise "width" up to
    # ~√(2 × 0.09) × 0.6° ≈ 0.25° is not significant.
    assert heldout_widths(residuals["t3phi"], 1, folds=4)[0] < 0.35


def test_one_night_cannot_be_cross_validated(tmp_path):
    res = calibrator_residuals(calibrator_nights(tmp_path, n_nights=1), INS)
    with pytest.raises(ValueError, match="two nights"):
        cross_validate(res["vis2"], 1)


def test_cli_end_to_end(tmp_path):
    (tmp_path / "data").mkdir()
    calibrator_nights(tmp_path / "data", n_nights=6)
    out = tmp_path / "out"
    main([str(tmp_path / "data"), "--out", str(out), "--k-max", "2", "--folds", "3",
          "--window", "1.9", "2.5"])  # fmt: skip
    summary = json.loads((out / "summary.json").read_text())
    assert summary["vis2"]["components_kept"] == 1
    assert summary["vis2"]["telescope_fraction"][0] > 0.9
    arrays = np.load(out / "calpca.npz")
    assert abs(arrays["vis2_components"][0] @ (SHAPE / np.linalg.norm(SHAPE))) > 0.98


def test_no_eligible_calibrators_is_a_clear_error():
    with pytest.raises(ValueError, match="pipeline diameter"):
        calibrator_residuals([], INS)


def test_telescope_fraction_null_and_common_mode():
    from virgil_vlti.calpca import Residuals

    rng = np.random.default_rng(9)
    pairs = [(a, b) for a in range(4) for b in range(a + 1, 4)]
    n_frames = 400
    meta = {
        "file": np.repeat([f"f{i}" for i in range(n_frames)], 6),
        "frame": np.zeros(6 * n_frames, int),
        "stations": np.tile([f"{a}-{b}" for a, b in pairs], n_frames),
    }
    res = Residuals("vis2", np.zeros(1), np.zeros((6 * n_frames, 1)), np.ones((6 * n_frames, 1)),
                    np.ones((6 * n_frames, 1), bool), meta)  # fmt: skip
    noise = rng.normal(size=(6 * n_frames, 1))
    assert telescope_fraction(res, noise)[0] == pytest.approx(4 / 6, abs=0.03)
    common = np.repeat(rng.normal(size=n_frames), 6)[:, None]  # the same gain on every baseline
    assert telescope_fraction(res, common)[0] == pytest.approx(1.0)


def test_heldout_widths_with_several_components_stay_finite(residuals):
    widths = heldout_widths(residuals["vis2"], 3, folds=4)
    assert np.all(np.isfinite(widths)) and widths[0] > 5 * widths[1:].max()


def test_cli_reports_no_components_when_none_beat_noise(tmp_path):
    (tmp_path / "data").mkdir()
    calibrator_nights(tmp_path / "data", n_nights=4, width=0.0)
    out = tmp_path / "out"
    main([str(tmp_path / "data"), "--out", str(out), "--k-max", "1", "--folds", "4"])
    summary = json.loads((out / "summary.json").read_text())
    arrays = np.load(out / "calpca.npz")
    for kind in ("vis2", "t3phi"):
        k = summary[kind]["components_kept"]
        assert arrays[f"{kind}_components"].shape[0] == k
        assert len(summary[kind]["heldout_widths"]) == k


def test_screening_removes_a_binary_calibrator(tmp_path):
    from virgil_vlti.calpca import screen_calibrators

    exposures = calibrator_nights(tmp_path, n_nights=6)
    binary = exposures[5]  # one star of the second night
    t3 = binary.t3[INS].data
    t3["T3PHI"] += 25.0 * np.sin(2 * np.pi * (2.2e-6 / WAVE))[None, :]  # a companion's phase
    v2 = binary.v2[INS].data
    v2["VIS2DATA"] *= 1 + 0.3 * np.cos(2 * np.pi * (2.2e-6 / WAVE))[None, :]
    res, excluded = screen_calibrators(exposures, INS)
    assert [e["target"] for e in excluded] == ["STAR1_1"]
    assert excluded[0]["t3phi_score"] > 5
    pca = weighted_pca(res["vis2"], 1)
    assert abs(pca.components[0] @ (SHAPE / np.linalg.norm(SHAPE))) > 0.98


def test_window_masks_channels(tmp_path):
    res = calibrator_residuals(calibrator_nights(tmp_path, n_nights=2), INS, window_um=(2.1, 2.3))
    inside = (WAVE > 2.1e-6) & (WAVE < 2.3e-6)
    assert not res["vis2"].valid[:, ~inside].any() and res["vis2"].valid[:, inside].all()


def test_screening_rounds_and_missing_scores(tmp_path):
    from virgil_vlti.calpca import screen_calibrators

    exposures = calibrator_nights(tmp_path, n_nights=6)
    binary = exposures[5]
    binary.t3[INS].data["T3PHI"] += 25.0 * np.sin(2 * np.pi * (2.2e-6 / WAVE))[None, :]
    # The last allowed round's exclusion is applied to what is returned.
    res, excluded = screen_calibrators(exposures, INS, max_rounds=1)
    assert [e["target"] for e in excluded] == ["STAR1_1"]
    assert "STAR1_1" not in set(res["t3phi"].meta["target"])
    # A star without usable V² is scored on closure phase alone (no NaN ranking).
    binary.v2[INS].data["FLAG"][:] = True
    res, excluded = screen_calibrators(exposures, INS)
    assert [e["target"] for e in excluded] == ["STAR1_1"]
    res, excluded = screen_calibrators(exposures, INS, max_rounds=0)
    assert excluded == [] and "vis2" in res


def test_clipping_drops_outlier_rows_and_widths_stay_robust(residuals):
    from virgil_vlti.calpca import Residuals, clip_rows

    v2 = residuals["vis2"]
    r = v2.r.copy()
    bad = np.arange(0, r.shape[0], 40)  # a few corrupted frames
    r[bad] += 0.8 * np.sin(np.arange(r.shape[1]) / 3.0)
    corrupt = Residuals(v2.kind, v2.wave, r, v2.sigma, v2.valid, v2.meta)
    clean, dropped = clip_rows(corrupt, k=1)
    assert set(bad) <= set(dropped) and len(dropped) < 2 * len(bad)
    expect = np.sqrt(2) * WIDTH * np.linalg.norm(SHAPE) * np.sqrt(1 + 1 / 3)
    assert heldout_widths(clean, 1, folds=4)[0] == pytest.approx(expect, rel=0.3)
    assert abs(weighted_pca(clean, 1).components[0] @ (SHAPE / np.linalg.norm(SHAPE))) > 0.98


def test_calibrator_residuals_refuses_binning(tmp_path):
    """The residuals are against the full-resolution model, so nbin > 1 is refused."""
    with pytest.raises(ValueError, match="nbin=3"):
        calibrator_residuals(calibrator_nights(tmp_path), INS, nbin=3)
