import json

import numpy as np
import pytest
from astropy.io import fits

from virgil_vlti import read_exposure
from virgil_vlti.pairs import (
    closure_correlation,
    closure_matrix,
    cp_block,
    find_pairs,
    frames_of,
    independent_triangles,
    pair_test,
)
from virgil_vlti.pairs_cli import main
from virgil_vlti.synthetic import write_product

INS = "GRAVITY_SC"
WAVE = np.linspace(2.0e-6, 2.4e-6, 30)
TRUTH = {"white": 0.6, "baseline": 0.01, "telescope": 0.015, "cp_white": 0.5, "offset": 0.4}
V2, V2ERR, CPERR = 0.5, 0.01, 1.0


def noisy_exposures(tmp_path, n_files=12, seed=1):
    """Flat V² and zero closure phase, with the structured noise of TRUTH injected."""
    rng = np.random.default_rng(seed)
    out = []
    for k in range(n_files):
        path = write_product(
            tmp_path / f"e{k}.fits",
            mjd=60400.1 + k * 15 / 1440,  # 15 min apart, 3 frames 5 min apart
            wave=WAVE,
            vis2=lambda u, v, w, t: np.full_like(w, V2),
            t3phi=lambda *a: np.zeros_like(a[4]),
            vis2err=V2ERR,
            t3phierr=CPERR,
        )
        e = read_exposure(path)
        v2, t3 = e.v2[INS].data, e.t3[INS].data
        for f in frames_of(e, INS):
            stations = sorted({s for sta in f.v2_rows for s in sta})
            gain_t = dict(zip(stations, rng.normal(0, TRUTH["telescope"], 4), strict=True))
            for sta, row in f.v2_rows.items():
                gain = gain_t[sta[0]] + gain_t[sta[1]] + rng.normal(0, TRUTH["baseline"])
                white = rng.normal(0, TRUTH["white"] * V2ERR, WAVE.size)
                v2["VIS2DATA"][row] = V2 * (1 + gain) + white
            triangles = list(f.t3_rows)
            T, baselines = closure_matrix(triangles)
            # Baseline-phase noise gives closure-phase noise of std cp_white·σ
            # with the ±1/3 structure between triangles.
            phi = rng.normal(0, TRUTH["cp_white"] * CPERR / np.sqrt(3), (len(baselines), WAVE.size))
            noise = T @ phi
            for r, tri in enumerate(triangles):
                t3["T3PHI"][f.t3_rows[tri]] = noise[r] + rng.normal(0, TRUTH["offset"])
        with fits.open(path, mode="update") as h:  # so that readers of the files see it too
            for hdu in h:
                if hdu.header.get("INSNAME") == INS and hdu.name in ("OI_VIS2", "OI_T3"):
                    col = "VIS2DATA" if hdu.name == "OI_VIS2" else "T3PHI"
                    hdu.data[col][:] = (v2 if hdu.name == "OI_VIS2" else t3)[col]
        out.append(e)
    return out


def test_closure_structure():
    triangles = [(1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4)]
    R = closure_correlation(triangles)
    np.testing.assert_allclose(np.diag(R), 1.0)
    assert np.allclose(np.abs(R[~np.eye(4, dtype=bool)]), 1 / 3)
    assert np.linalg.matrix_rank(R) == 3
    assert len(independent_triangles(triangles)) == 3


def test_frames_split_by_occurrence(tmp_path):
    (e,) = noisy_exposures(tmp_path, n_files=1)
    frames = frames_of(e, INS)
    assert len(frames) == 3
    assert all(len(f.v2_rows) == 6 and len(f.t3_rows) == 4 for f in frames)
    assert frames[0].mjd < frames[1].mjd < frames[2].mjd


def test_pair_test_recovers_injected_noise(tmp_path):
    result = pair_test(noisy_exposures(tmp_path), INS)
    assert result.n_pairs == 12 * 3 // 2  # disjoint pairs of the 36 frames
    full = result.v2["+ baseline + telescope"].widths
    assert full["white"] == pytest.approx(TRUTH["white"], rel=0.05)
    assert full["baseline"] == pytest.approx(TRUTH["baseline"], rel=0.35)
    assert full["telescope"] == pytest.approx(TRUTH["telescope"], rel=0.3)
    assert result.v2["+ telescope"].delta_lnl > 20
    cp = result.cp["+ triangle offset"].widths
    assert cp["white"] == pytest.approx(TRUTH["cp_white"], rel=0.05)
    assert cp["offset"] == pytest.approx(TRUTH["offset"], rel=0.3)
    print(result.summary())


def test_white_noise_alone_needs_no_correlated_terms(tmp_path, monkeypatch):
    monkeypatch.setitem(TRUTH, "baseline", 0.0)
    monkeypatch.setitem(TRUTH, "telescope", 0.0)
    monkeypatch.setitem(TRUTH, "offset", 0.0)
    result = pair_test(noisy_exposures(tmp_path, n_files=6, seed=2), INS)
    assert result.v2["+ baseline + telescope"].delta_lnl < 5
    assert result.cp["+ triangle offset"].delta_lnl < 3


def test_cli_writes_json(tmp_path):
    noisy_exposures(tmp_path, n_files=3)
    out = tmp_path / "pairs.json"
    main([str(tmp_path), "--json", str(out)])
    result = json.loads(out.read_text())
    assert result["n_pairs"] == 4
    widths = result["v2"]["+ baseline + telescope"]["widths"]
    assert widths["white"] == pytest.approx(TRUTH["white"], rel=0.2)


def test_frames_by_time_survive_a_missing_row(tmp_path):
    (e,) = noisy_exposures(tmp_path, n_files=1)
    v2 = e.v2[INS]
    v2.data = v2.data[np.arange(len(v2.data)) != 7]  # frame 1 loses one baseline
    frames = frames_of(e, INS)
    assert [len(f.v2_rows) for f in frames] == [6, 5, 6]
    for f in frames:  # every row of a frame is from that frame's time
        times = v2.data["MJD"][list(f.v2_rows.values())]
        assert np.ptp(times) * 86400 < 60


def test_pairs_are_disjoint_and_need_one_wavelength_grid(tmp_path):
    exposures = noisy_exposures(tmp_path, n_files=2)
    pairs = find_pairs(exposures, INS)
    frames = [f for pair in pairs for f in pair]
    assert len({id(f) for f in frames}) == len(frames)
    exposures[1].wave[INS] = exposures[1].wave[INS] * 1.01
    pairs = find_pairs(exposures, INS)
    assert all(a.exposure is b.exposure for a, b in pairs)


def test_unusable_samples_are_left_out(tmp_path):
    exposures = noisy_exposures(tmp_path, n_files=2)
    (a, b), *_ = find_pairs(exposures, INS)
    t3 = a.exposure.t3[INS].data
    first = list(a.t3_rows.values())
    t3["T3PHIERR"][first[0], :] = np.inf  # one triangle unusable in frame a
    block = cp_block(a, b)
    nchan = WAVE.size
    assert block.diff.size == 3 * nchan  # another triangle takes its place
    assert np.all(np.isfinite(block.terms["white"]))
    for row in first:
        t3["FLAG"][row, :] = True
    assert cp_block(a, b) is None


def test_add_systematics_is_recovered_by_the_pair_test(tmp_path):
    from virgil_vlti.synthetic import add_systematics

    rng = np.random.default_rng(3)
    paths = []
    for k in range(12):
        path = write_product(
            tmp_path / f"s{k}.fits",
            mjd=60400.1 + k * 15 / 1440,
            wave=WAVE,
            vis2=lambda u, v, w, t: np.full_like(w, V2),
            t3phi=lambda *a: np.zeros_like(a[4]),
            vis2err=V2ERR,
            t3phierr=CPERR,
        )
        add_systematics(path, INS, rng=rng, white=0.6, telescope=0.015, baseline=0.01,
                        cp_white=0.5, cp_offset=0.4)  # fmt: skip
        paths.append(path)
    result = pair_test([read_exposure(str(p)) for p in paths], INS)
    full = result.v2["+ baseline + telescope"].widths
    assert full["white"] == pytest.approx(0.6, rel=0.05)
    assert full["telescope"] == pytest.approx(0.015, rel=0.3)
    assert full["baseline"] == pytest.approx(0.01, rel=0.35)
    cp = result.cp["+ triangle offset"].widths
    assert cp["white"] == pytest.approx(0.5, rel=0.05)
    assert cp["offset"] == pytest.approx(0.4, rel=0.3)
