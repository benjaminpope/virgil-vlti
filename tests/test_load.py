"""The loader from calibrated files and their mode files to a virgil OIData."""

import warnings

import jax
import numpy as np
import pytest
from astropy.io import fits

pytest.importorskip("virgil")

import virgil.oidata as _oidata  # noqa: E402
from virgil.oidata import OIData  # noqa: E402

# to_oidata uses virgil API (OIData.select(ranges=), _valid_samples) that the 0.3.0 release
# lacks; the "against virgil main" CI job runs these.
pytestmark = pytest.mark.skipif(
    not hasattr(_oidata, "_valid_samples"), reason="needs a virgil newer than 0.3.0"
)

from virgil_vlti import to_oidata  # noqa: E402
from virgil_vlti.templates import gain_template_from_pairs  # noqa: E402
from virgil_vlti.write import write_modes  # noqa: E402

from .test_calibrate import MJD0, product  # noqa: E402
from .test_calpca import INS  # noqa: E402
from .test_templates import _pairs_json  # noqa: E402

WINDOW = [(2.1e-6, 2.3e-6)]


def flagged_file(tmp_path, name="sci", t=MJD0):
    """A calibrated-looking file with some V² and closure phases flagged."""
    path = product(tmp_path, name, t, cal=False)
    with fits.open(path, mode="update") as h:
        for ext, col, cells in (("OI_VIS2", "FLAG", [(0, 3), (2, 10), (5, 39)]),
                                ("OI_T3", "FLAG", [(1, 4), (3, 20)])):  # fmt: skip
            hdu = next(x for x in h if x.name == ext and x.header["INSNAME"] == INS)
            for row, chan in cells:
                hdu.data[col][row, chan] = True
    return path


def tables(path):
    with fits.open(path) as h:
        v2 = next(x for x in h if x.name == "OI_VIS2" and x.header["INSNAME"] == INS).data
        t3 = next(x for x in h if x.name == "OI_T3" and x.header["INSNAME"] == INS).data
    return v2, t3


def beside(path, suffix, modes):
    """Write ``modes`` as ``<file>.<suffix>.npz`` beside a ``.fits`` file."""
    write_modes(str(path)[: -len(".fits")] + f".{suffix}.npz", {INS: modes})


def wave_of(path):
    with fits.open(path) as h:
        table = next(x for x in h if x.name == "OI_WAVELENGTH" and x.header["INSNAME"] == INS)
        return table.data["EFF_WAVE"]


def cell_ids(shape, base):
    """Each cell of a table gets its own value, so any misalignment shows."""
    return base + np.arange(np.prod(shape), dtype=float).reshape(shape) / 1e4


def write_pattern_modes(path, pattern=1.0):
    v2, t3 = tables(path)
    modes = {
        "vis2-pc1-tel28@x": {"vis2": cell_ids(v2["VIS2DATA"].shape, pattern)},
        "vis2-pc1-base28-23@x": {"vis2": cell_ids(v2["VIS2DATA"].shape, 3 * pattern)},
    }
    ids = cell_ids(t3["T3PHI"].shape, 7 * pattern)
    for row in range(len(ids)):  # one mode per triangle row, so each lies in one frame
        only = np.zeros_like(ids)
        only[row] = ids[row]
        modes[f"t3phi-pc1-row{row:02d}@x"] = {"t3phi": only}
    beside(path, "template_modes", modes)
    return modes


def expected_gain_modes(plain, v2, wave, patterns):
    """The modes by hand: the cells of the table (row by row) inside the window, in order.

    Checked against the samples' own baseline and wavelength, so that a
    misaligned layout cannot pass.
    """
    nrow, nchan = v2["VIS2DATA"].shape
    cells = np.flatnonzero(((wave >= WINDOW[0][0]) & (wave <= WINDOW[0][1]))[None].repeat(nrow, 0))
    assert cells.size == plain.u.size
    row, chan = np.divmod(cells, nchan)
    np.testing.assert_allclose(plain.u, v2["UCOORD"][row])
    np.testing.assert_allclose(plain.v, v2["VCOORD"][row])
    np.testing.assert_allclose(np.broadcast_to(plain.wavel, plain.u.shape), wave[chan])
    return np.array([p.reshape(-1)[cells] for p in patterns])


def leaves(x):
    return [np.asarray(leaf) for leaf in jax.tree_util.tree_leaves(x)]


def assert_same(a, b):
    la, lb = leaves(a), leaves(b)
    assert len(la) == len(lb)
    for x, y in zip(la, lb, strict=True):
        np.testing.assert_array_equal(x, y)


def test_gain_modes_align_after_window_and_flags(tmp_path):
    path = flagged_file(tmp_path)
    modes = write_pattern_modes(path)
    v2, _ = tables(path)
    data = to_oidata(path, INS, ranges=WINDOW)
    assert data.gains is not None
    # the same selection without modes, with the modes laid out by hand
    from virgil.oifits import read_oifits

    bare = OIData(read_oifits(str(path), insname=INS)).select(ranges=WINDOW)
    assert data.u.size == bare.u.size
    patterns = [modes[k]["vis2"] for k in sorted(modes) if "vis2" in modes[k]]
    ref = bare.with_gains(modes=expected_gain_modes(bare, v2, wave_of(path), patterns))
    assert_same(data.gains, ref.gains)


def test_closure_modes_align_after_window_and_flags(tmp_path):
    path = flagged_file(tmp_path)
    modes = write_pattern_modes(path)
    _, t3 = tables(path)
    data = to_oidata(path, INS, ranges=WINDOW)
    assert data.phase_offsets is not None
    from virgil.oifits import read_oifits

    bare = OIData(read_oifits(str(path), insname=INS)).select(ranges=WINDOW)
    nrow, nchan = t3["T3PHI"].shape
    wave = wave_of(path)
    inside = ((wave >= WINDOW[0][0]) & (wave <= WINDOW[0][1]))[None].repeat(nrow, 0)
    cells = np.flatnonzero(inside & ~t3["FLAG"].astype(bool))  # row by row: unflagged, in window
    np.testing.assert_allclose(np.asarray(bare.phi), np.deg2rad(t3["T3PHI"].reshape(-1)[cells]))
    rows = sorted(k for k in modes if "t3phi" in modes[k])
    expected = np.array([np.deg2rad(modes[k]["t3phi"].reshape(-1)[cells]) for k in rows])
    assert expected.shape == (nrow, cells.size)
    ref = bare.with_closure_offsets(modes=expected)
    assert_same(data.phase_offsets, ref.phase_offsets)


def test_zero_modes_reproduce_plain_oidata(tmp_path):
    from virgil.oifits import read_oifits

    path = flagged_file(tmp_path)
    plain = OIData(read_oifits(str(path), insname=INS))
    assert_same(to_oidata(path, INS), plain)  # no mode files at all
    v2, t3 = tables(path)
    zero = {"vis2-pc1-tel28@x": {"vis2": np.zeros(v2["VIS2DATA"].shape)},
            "t3phi-pc1-row00@x": {"t3phi": np.zeros(t3["T3PHI"].shape)}}  # fmt: skip
    beside(path, "template_modes", zero)
    assert_same(to_oidata(path, INS), plain)
    # a mode that is nonzero only on flagged samples is zero on every sample virgil keeps
    only_flagged = {"vis2": np.where(v2["FLAG"], 1.0, 0.0)}
    beside(path, "template_modes", {"vis2-flagged@x": only_flagged})
    assert_same(to_oidata(path, INS), plain)


def test_two_files_do_not_share_template_modes_but_share_diameter_modes(tmp_path):
    from virgil.oifits import read_oifits

    dirs = [tmp_path / "n1", tmp_path / "n2"]
    paths = []
    for d, t in zip(dirs, (MJD0, MJD0 + 0.1), strict=True):
        d.mkdir()
        paths.append(flagged_file(d, "same", t))  # the same basename on two nights
    shape = tables(paths[0])[0]["VIS2DATA"].shape
    for path, scale in zip(paths, (1.0, 2.0), strict=True):
        beside(
            path, "template_modes", {"vis2-pc1-tel28@same.fits#frame0": {"vis2": np.ones(shape)}}
        )
        diameter = {"rel-cal@night": {"vis2": scale * np.ones(shape)},
                    "cal@night": {"vis2": np.ones(shape)}}  # fmt: skip
        beside(path, "modes", diameter)
    with pytest.warns(UserWarning, match="ignored modes"):
        data = to_oidata(paths, INS)
    bare = OIData(read_oifits([str(p) for p in paths], insname=INS))
    n = shape[0] * shape[1]
    assert bare.u.size == 2 * n
    good = [~tables(p)[0]["FLAG"].astype(bool).reshape(-1) for p in paths]
    # own modes of the two files (full path first), then the shared rel- diameter mode;
    # the absolute one is ignored
    expected = np.zeros((3, 2 * n))
    expected[0, :n], expected[1, n:] = good[0], good[1]
    expected[2, :n], expected[2, n:] = good[0] * 1.0, good[1] * 2.0
    assert_same(data.gains, bare.with_gains(modes=expected).gains)


def test_t3amp_and_absolute_modes_are_ignored_with_a_warning(tmp_path):
    path = flagged_file(tmp_path)
    v2, t3 = tables(path)
    beside(path, "modes", {"cal@night": {"vis2": np.ones(v2["VIS2DATA"].shape),
                                         "t3amp": np.ones(t3["T3PHI"].shape)}})  # fmt: skip
    with pytest.warns(UserWarning, match="ignored modes.*cal@night"):
        data = to_oidata(path, INS)
    assert data.gains is None and data.phase_offsets is None


def test_use_template_false_and_generic_fallback(tmp_path):
    path = flagged_file(tmp_path)
    write_pattern_modes(path)
    assert to_oidata(path, INS, use_template=False).gains is None
    widths = {"telescope": 0.02, "chromatic": 0.01}
    bare = to_oidata(path, INS, use_template=False)
    assert_same(to_oidata(path, INS, use_template=False, generic_gains=widths).gains,
                bare.with_gains(**widths).gains)  # fmt: skip
    with pytest.warns(UserWarning, match="counted twice"):
        with_modes = to_oidata(path, INS, generic_gains=widths)
    assert_same(with_modes.gains, to_oidata(path, INS).gains)
    offsets = to_oidata(path, INS, use_template=False, cp_offsets={"triangle": 0.01})
    assert offsets.phase_offsets is not None


def test_layout_is_checked_against_the_tables(tmp_path):
    path = flagged_file(tmp_path)
    write_pattern_modes(path)
    from virgil_vlti import load

    original = load._cells
    load._cells = lambda shape: tuple(np.asarray(i)[::-1] for i in original(shape))
    try:
        with pytest.raises(RuntimeError, match="not laid out"):
            to_oidata(path, INS)
    finally:
        load._cells = original


def test_round_trip_with_real_virgil(tmp_path):
    """Modes in, finite whitened residuals and likelihood out, with and without the modes."""
    from virgil.likelihood import loglike, whitened_residuals
    from virgil.models import GaussianDisk

    path = flagged_file(tmp_path)
    write_pattern_modes(path, pattern=0.02)
    model = GaussianDisk(sigma=2.0)
    for data in (to_oidata(path, INS, ranges=WINDOW), to_oidata(path, INS, use_template=False)):
        assert np.isfinite(
            np.concatenate(
                [np.ravel(r) for r in jax.tree_util.tree_leaves(whitened_residuals(model, data))]
            )
        ).all()
        assert np.isfinite(float(loglike([2.0], ["sigma"], model, data)))
    with_modes = to_oidata(path, INS, ranges=WINDOW)
    assert with_modes.gains is not None and with_modes.phase_offsets is not None


def test_flat_template_matches_its_written_modes(tmp_path):
    """A pair-test template through the loader gives the modes ``virgil-vlti-template`` writes."""
    from virgil_vlti import classify, read_exposure
    from virgil_vlti.templates import template_modes

    path = flagged_file(tmp_path)
    t = gain_template_from_pairs(_pairs_json(0.06, 0.0, 1.4), "gains", config={"field": "DUAL"})
    direct = to_oidata(path, INS, template=t)
    (sci,) = classify([read_exposure(path)], dual_field=True)
    beside(path, "template_modes", template_modes(t, sci, INS))
    assert direct.gains is not None and direct.phase_offsets is not None
    # the widths by hand, as the pair test gives them: g = 0.06 per telescope is 2 g² of V² variance
    # (relative) on every baseline, and the closure offset 1.4° per triangle
    frame0 = [
        m["vis2"] for k, m in template_modes(t, sci, INS).items() if "#frame0" in k and "vis2" in m
    ]
    variance = np.sum([(2 * m) ** 2 for m in frame0], axis=0)
    np.testing.assert_allclose(variance[variance != 0], 2 * 0.06**2)
    assert_same(direct, to_oidata(path, INS))
    flat = template_modes(t, sci, INS)  # 4 telescope gains in each of 3 frames; no baseline gains
    assert sum("vis2" in by_obs for by_obs in flat.values()) == 12


def test_observables_leave_out_what_they_drop(tmp_path):
    path = flagged_file(tmp_path)
    write_pattern_modes(path)
    vis = to_oidata(path, INS, observables="vis")
    phi = to_oidata(path, INS, observables="phi")
    assert vis.gains is not None and vis.phase_offsets is None
    assert phi.gains is None and phi.phase_offsets is not None


def test_missing_modes_warn_and_wrong_size_refused(tmp_path):
    a, b = flagged_file(tmp_path, "a"), flagged_file(tmp_path, "b", MJD0 + 0.1)
    write_pattern_modes(a)
    with pytest.warns(UserWarning, match="no nuisance modes"):
        to_oidata([a, b], INS)
    beside(b, "template_modes", {"vis2-m": {"vis2": np.ones((2, 2))}})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with pytest.raises(ValueError, match="were the modes written for this file"):
            to_oidata([a, b], INS)


def test_shipped_high_template_through_the_loader(tmp_path):
    header = {
        "ESO INS SPEC RES": "HIGH",
        "ESO QC FIELD MODE": "SINGLE",
        "ESO ISS CONF T1NAME": "AT",
    }
    path = product(tmp_path, "sci", MJD0, cal=False, header=header)
    data = to_oidata(path, INS, template="gravity_single_high_combined_2023-2025")
    assert data.gains is not None and data.phase_offsets is not None


def test_transfer_function_modes_are_read_and_shared_by_night(tmp_path):
    from virgil.oifits import read_oifits

    paths = [flagged_file(tmp_path, "a"), flagged_file(tmp_path, "b", MJD0 + 0.1)]
    v2, t3 = tables(paths[0])
    for path in paths:  # the same night: the same parameters
        tf = {"tf-vis2-base28-23@night": {"vis2": np.ones(v2["VIS2DATA"].shape)},
              "tf-t3phi-tri28-23-18@night": {"t3phi": np.ones(t3["T3PHI"].shape)}}  # fmt: skip
        beside(path, "tf_modes", tf)
    with pytest.warns(UserWarning, match="span several frames"):
        data = to_oidata(paths, INS)
    bare = OIData(read_oifits([str(p) for p in paths], insname=INS))
    n = v2["VIS2DATA"].size
    good = np.concatenate([~tables(p)[0]["FLAG"].astype(bool).reshape(-1) for p in paths])
    assert bare.u.size == 2 * n
    assert_same(data.gains, bare.with_gains(modes=good[None] * 1.0).gains)
    assert data.phase_offsets is None  # the night-long closure mode cannot be used


def test_no_modes_at_all_warns(tmp_path):
    path = flagged_file(tmp_path)
    with pytest.warns(UserWarning, match="no nuisance modes beside"):
        to_oidata(path, INS)
