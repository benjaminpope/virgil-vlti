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


P4 = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
TRI4 = [
    (0, 3, 1),
    (0, 4, 2),
    (1, 5, 2),
    (3, 5, 4),
]  # baseline indices of triangles 123, 124, 134, 234


def toy4_data():
    """Four telescopes, two nights x two frames, six baselines, four triangles, one channel."""
    rows = [(f, m, a, c) for f, m in enumerate(MJDS) for a, c in P4]
    r = np.array(rows)
    n = len(rows)
    d = {
        "u": np.linspace(10, 60, n),
        "v": np.linspace(-20, 30, n),
        "vis": np.full(n, 0.8),
        "d_vis": np.full(n, 0.02),
        "vis_flag": np.zeros(n, bool),
        "wavel": np.full(n, 1.65e-6),
        "mjd": r[:, 1],
        "stations": r[:, 2:4].astype(int),
        "frame": r[:, 0].astype(int),
    }
    ids = np.array([[6 * f + k for k in t] for f in range(4) for t in TRI4])
    d.update(
        phi=np.zeros(len(ids)),
        d_phi=np.full(len(ids), 0.02),
        phi_flag=np.zeros(len(ids), bool),
        i_cps1=ids[:, 0],
        i_cps2=ids[:, 1],
        i_cps3=ids[:, 2],
    )
    return virgil.OIData(d)


def test_cells_hold_one_baseline_or_triangle_in_one_night():
    d = toy4_data()
    lab = night_labels(d)
    st = np.sort(np.asarray(d.stations), axis=1)
    pair = [tuple(x) for x in st]
    for c in np.unique(lab["v2_cell"]):
        m = lab["v2_cell"] == c
        assert len({pair[i] for i in np.flatnonzero(m)}) == 1  # one baseline
        assert len(np.unique(lab["v2_night"][m])) == 1  # one night
        assert len(np.unique(lab["v2_frame"][m])) == 2  # both frames of that night
    assert len(np.unique(lab["v2_cell"])) == 12  # 2 nights x 6 baselines
    tri = [
        frozenset(np.concatenate([st[a], st[b], st[c]]).tolist())
        for a, b, c in zip(d.i_cps1, d.i_cps2, d.i_cps3, strict=True)
    ]
    assert len(np.unique(lab["t3_cell"])) == 8  # 2 nights x 4 triangles
    for c in np.unique(lab["t3_cell"]):
        m = np.flatnonzero(lab["t3_cell"] == c)
        assert len({tri[i] for i in m}) == 1  # one triangle
        assert len(np.unique(lab["t3_night"][m])) == 1
        assert len(np.unique(lab["t3_frame"][m])) == 2


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
    # widths: per-frame telescope and baseline (log|V|), triangle in radians; modes carry theirs
    assert np.allclose(nd.gains.widths[:2], [0.005, 0.004])
    assert np.allclose(nd.phase_offsets.widths[0], np.radians(0.2))
    assert np.allclose(np.unique(np.asarray(nd.gains.spanning))[-1], 0.02)
    assert np.allclose(np.unique(np.asarray(nd.phase_offsets.spanning))[-1], np.radians(1.0))
    n0 = with_night_terms(d, telescope=0.005, baseline=0.004, triangle_deg=0.2)
    assert n0.gains.spanning is None and n0.phase_offsets.spanning is None
    assert with_night_terms(d).gains is None
    # only night_gain: no zero-width per-frame groups
    assert list(with_night_terms(d, night_gain=0.02).gains.widths) == [1.0]


def test_with_night_terms_refuses_to_replace_existing_terms():
    d = toy_data()
    with_gains = d.with_gains(modes=np.full((1, d.u.size), 0.01))
    with pytest.raises(ValueError, match="already has gains"):
        with_night_terms(with_gains, night_gain=0.02)
    with pytest.raises(ValueError, match="already has gains"):
        with_night_terms(with_gains, telescope=0.005)
    kept = with_night_terms(with_gains, triangle_deg=0.2)  # a group not requested is untouched
    assert np.array_equal(np.asarray(kept.gains.spanning), np.asarray(with_gains.gains.spanning))
    with_offsets = d.with_closure_offsets(triangle=0.01)
    with pytest.raises(ValueError, match="already has closure offsets"):
        with_night_terms(with_offsets, triangle_deg=0.2)
