"""virgil_vlti.noise: night-shared gains and closure offsets on a virgil OIData."""

import numpy as np
import pytest

from virgil_vlti.noise import night_labels, night_modes, with_night_terms

virgil = pytest.importorskip("virgil")

PAIRS = [(1, 2), (2, 3), (1, 3)]
# two nights of two frames; MJD fraction .1/.2 is before 12 UT: the night that began the day before
MJDS = [100.1, 100.2, 101.1, 101.2]


def toy_data():
    """Three telescopes, two nights x two frames, three baselines, two channels."""
    rows = [(f, m, a, c, w) for f, m in enumerate(MJDS) for a, c in PAIRS for w in (1.6e-6, 1.7e-6)]
    r = np.array(rows)
    n = len(rows)
    d = {
        "u": np.linspace(10, 60, n),
        "v": np.linspace(-20, 30, n),
        "vis": np.full(n, 0.8),
        "d_vis": np.full(n, 0.02),
        "vis_flag": np.zeros(n, bool),
        "wavel": r[:, 4],
        "mjd": r[:, 1],
        "stations": r[:, 2:4].astype(int),
        "frame": r[:, 0].astype(int),
    }
    i1 = [6 * f + k for f in range(4) for k in range(2)]
    d.update(
        phi=np.zeros(len(i1)),
        d_phi=np.full(len(i1), 0.02),
        phi_flag=np.zeros(len(i1), bool),
        i_cps1=np.array(i1),
        i_cps2=np.array(i1) + 2,
        i_cps3=np.array(i1) + 4,
    )
    return virgil.OIData(d)


def test_labels_group_samples_by_night_and_baseline_or_triangle():
    lab = night_labels(toy_data())
    assert len(np.unique(lab["u_night"])) == 2
    assert len(np.unique(lab["v2_cell"])) == 6  # 2 nights x 3 baselines
    assert len(np.unique(lab["t3_cell"])) == 2  # 2 nights x 1 triangle
    for c in np.unique(lab["v2_cell"]):  # a cell spans both frames of its night, on one baseline
        assert len(np.unique(lab["v2_frame"][lab["v2_cell"] == c])) == 2
        assert len(np.unique(lab["v2_night"][lab["v2_cell"] == c])) == 1


def test_the_night_boundary_is_noon_ut_by_default():
    d = toy_data()
    assert len(np.unique(night_labels(d)["u_night"])) == 2
    assert len(np.unique(night_labels(d, boundary_ut=12.0)["u_night"])) == 2
    # a boundary at 3.6 UT (MJD .15) falls between the first two frames of a night's morning
    assert len(np.unique(night_labels(d, boundary_ut=3.6)["u_night"])) == 3


def test_modes_have_one_row_per_cell_and_zero_width_gives_none():
    d = toy_data()
    g, o = night_modes(d, 0.02, 1.0)
    assert g.shape == (6, d.u.size) and o.shape == (2, np.size(d.phi))
    assert set(np.unique(g[g > 0])) == {0.02}
    assert np.allclose(np.unique(o[o > 0]), np.radians(1.0))
    assert np.all((g > 0).sum(axis=0) == 1) and np.all((o > 0).sum(axis=0) == 1)
    assert night_modes(d, 0.0, 0.0) == (None, None)


def test_with_night_terms_keeps_frame_and_night_terms_in_one_call():
    d = toy_data()
    try:
        d.with_closure_offsets(modes=np.ones((1, d.phi.size)))
    except ValueError as e:  # virgil before closure modes could span frames
        pytest.skip(f"this virgil has no frame-spanning closure modes: {e}")
    nd = with_night_terms(
        d, night_gain=0.02, night_t3_deg=1.0, telescope=0.005, baseline=0.004, triangle_deg=0.2
    )
    assert nd.gains.spanning.shape == (d.u.size, 6)
    assert np.shape(nd.phase_offsets.spanning) == (d.phi.size, 2)
    assert int((np.asarray(nd.gains.group) >= 0).sum()) > 0  # the per-frame gains are still there
    n0 = with_night_terms(d, telescope=0.005, baseline=0.004, triangle_deg=0.2)
    assert n0.gains.spanning is None and n0.phase_offsets.spanning is None
    assert with_night_terms(d).gains is None
