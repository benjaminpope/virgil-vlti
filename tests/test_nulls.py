"""Empirical null datasets and block resampling, on synthetic calibrator nights."""

import json
import os

import numpy as np
import pytest
from astropy.io import fits

from virgil_vlti import read_exposure
from virgil_vlti.calibrate import uniform_disk_v2
from virgil_vlti.calpca import calibrator_residuals
from virgil_vlti.nulls import block_resample, null_datasets, resampled_null_datasets
from virgil_vlti.products import read_flux

from .test_calpca import INS, WAVE, calibrator_nights


@pytest.fixture(scope="module")
def exposures(tmp_path_factory):
    return calibrator_nights(tmp_path_factory.mktemp("cal"), n_nights=3, per_night=4)


@pytest.fixture(scope="module")
def nulls(exposures, tmp_path_factory):
    return null_datasets(exposures, tmp_path_factory.mktemp("null"), INS, vis=True)


def test_one_file_per_calibrator_exposure(nulls, exposures):
    assert len(nulls.files) == len(exposures) == 12
    assert not nulls.dropped and not nulls.skipped
    with open(nulls.files[0].rsplit("/", 2)[0] + "/null_manifest.json") as fh:
        manifest = json.load(fh)
    assert len(manifest["files"]) == 12


def test_tables_and_shapes(nulls, exposures):
    src = {str(e.path).split("/")[-1].replace(".fits", ""): e for e in exposures}
    for path in nulls.files:
        out = read_exposure(path)
        e = src[path.split("/")[-1].replace("_null.fits", "")]
        assert out.v2[INS].data["VIS2DATA"].shape == e.v2[INS].data["VIS2DATA"].shape
        assert out.t3[INS].data["T3PHI"].shape == e.t3[INS].data["T3PHI"].shape
        np.testing.assert_allclose(out.wave[INS], WAVE)
        with fits.open(path) as h:
            names = [x.name for x in h]
        assert {"OI_VIS2", "OI_T3", "OI_FLUX"} <= set(names)
        assert read_flux(out, INS).shape[1] == WAVE.size


def test_null_is_the_catalogue_uniform_disk(tmp_path):
    """Truth-based: stars that are 0.5 mas disks, with no systematic, calibrate to that disk."""
    (tmp_path / "in").mkdir()
    exps = calibrator_nights(
        tmp_path / "in", n_nights=2, per_night=4, width=0.0, diameter=0.5, disk=True
    )
    got = null_datasets(exps, tmp_path / "out", INS, screen=False, pair_factor=None, vis=False)
    src = {str(e.path).split("/")[-1].replace(".fits", ""): e for e in exps}
    chi2 = n = 0
    for path in got.files:
        out = read_exposure(path)
        e = src[path.split("/")[-1].replace("_null.fits", "")]
        v2 = out.v2[INS]
        model = uniform_disk_v2(v2, WAVE, 0.5)
        err = v2.data["VIS2ERR"]
        ok = ~v2.data["FLAG"] & (err > 0)
        chi2 += np.sum(((v2.data["VIS2DATA"] - model) / err)[ok] ** 2)
        n += ok.sum()
        # the propagated error includes the source's own
        assert np.all(err[ok] >= e.v2[INS].data["VIS2ERR"][ok] * (1 - 1e-9))
        assert abs(np.mean(out.t3[INS].data["T3PHI"])) < 0.5
    assert 0.7 < chi2 / n < 1.4
    assert all(m["diameter_mas"] == 0.5 for m in got.manifest)
    assert all(m["diameter_source"] == "ESO QC TF CAL_SC DIAM" for m in got.manifest)


def test_exclude_is_target_and_reference(exposures, tmp_path):
    got = null_datasets(exposures, tmp_path, INS, exclude=["STAR0_1"], screen=False)
    assert len(got.files) == 11
    assert not any("n0c1_" in f for f in got.files)
    assert all("STAR0_1" not in " ".join(m["references"]) for m in got.manifest)


def test_min_calibrators_per_night(exposures, tmp_path):
    with pytest.raises(ValueError, match="enough other calibrators"):
        null_datasets(exposures, tmp_path, INS, min_calibrators_per_night=4, screen=False)
    got = null_datasets(exposures, tmp_path, INS, exclude=["STAR0_1"], screen=False)
    assert len(got.manifest[0]["references"]) >= 2


def test_mode_filters(exposures, tmp_path):
    with pytest.raises(ValueError, match="passes the filters"):
        null_datasets(exposures, tmp_path, INS, resolution="HIGH")
    with pytest.raises(ValueError, match="passes the filters"):
        null_datasets(exposures, tmp_path, INS, telescope="UT")
    with pytest.raises(ValueError, match="polarisation must be"):
        null_datasets(exposures, tmp_path, INS, polarisation="BOTH")
    got = null_datasets(
        exposures, tmp_path, INS, resolution="MEDIUM", polarisation="COMBINED", telescope="AT"
    )
    assert len(got.files) == 12


def test_binary_calibrator_is_dropped_and_reported(tmp_path):
    (tmp_path / "in").mkdir()
    exposures = calibrator_nights(tmp_path / "in", n_nights=4, per_night=6, width=0.0)
    star = "STAR1_2"
    for e in exposures:
        if e.target == star:
            e.v2[INS].data["VIS2DATA"][:] *= 0.5  # a bad star
            e.t3[INS].data["T3PHI"][:] += 40.0
    got = null_datasets(exposures, tmp_path / "out", INS)
    # A bad reference spoils its night's other stars too, so others may go with it. With
    # the error floors the weights are near-equal, so six stars a night are used: the
    # bad one is then still found, after good ones.
    assert star in [d["target"] for d in got.dropped]
    assert {d["reason"] for d in got.dropped} == {"screen"}
    # flagged stars are kept out as references only: their own nulls are written
    flagged = {m["target"] for m in got.manifest if m["flagged"]}
    assert star in flagged
    assert any("n1c2_" in f for f in got.files)
    assert all(star not in " ".join(m["references"]).upper() for m in got.manifest)
    assert all(m["flag_reasons"] == ["screen"] for m in got.manifest if m["flagged"])


def test_pair_test_drops_a_star_that_changes(tmp_path):
    (tmp_path / "in").mkdir()
    exposures = calibrator_nights(
        tmp_path / "in", n_nights=3, per_night=6, width=0.0, stars_per_night=3
    )
    # Three stars a night, two exposures each. One star's V² scatters between its exposures
    # by far more than its errors, channel by channel (a gain per baseline would be absorbed
    # by the fitted correlated terms).
    rng = np.random.default_rng(0)
    for e in (e for e in exposures if e.target == "STAR1_0"):
        d = e.v2[INS].data["VIS2DATA"]
        d *= 1.0 + 0.1 * rng.normal(size=d.shape)
    got = null_datasets(
        exposures, tmp_path / "out", INS, screen=False, pair_factor=3.0, pair_max_gap_min=300.0
    )
    assert [d["target"] for d in got.dropped] == ["STAR1_0"]
    assert got.dropped[0]["reason"] == "pair_test"
    flagged = [m for m in got.manifest if m["flagged"]]
    assert {m["target"] for m in flagged} == {"STAR1_0"}
    assert flagged[0]["flag_reasons"] == ["pair_test"]
    assert any("n1c0_" in f for f in got.files)  # its nulls are kept, flagged
    assert not any("STAR1_0" in " ".join(m["references"]).upper() for m in got.manifest)


def test_window_is_refused(exposures, tmp_path):
    with pytest.raises(ValueError, match="full-resolution"):
        null_datasets(exposures, tmp_path, INS, window_um=(2.0, 2.3))


# --- block resampling -------------------------------------------------------


@pytest.fixture(scope="module")
def residuals(exposures):
    return calibrator_residuals(exposures, INS)


def test_block_resample_shapes_and_seed(residuals):
    a = block_resample(residuals, 3, "exposure", seed=1)
    b = block_resample(residuals, 3, "exposure", seed=1)
    c = block_resample(residuals, 3, "exposure", seed=2)
    assert len(a) == 3
    for kind in ("vis2", "t3phi"):
        assert a[0][kind].r.shape == residuals[kind].r.shape
        np.testing.assert_array_equal(a[0][kind].r, b[0][kind].r)
        assert a[0][kind].r is not residuals[kind].r
    assert not np.array_equal(a[0]["vis2"].r, c[0]["vis2"].r)
    assert not np.array_equal(a[0]["vis2"].r, a[1]["vis2"].r)


def test_exposure_blocks_stay_intact(residuals):
    v2 = residuals["vis2"]
    real = block_resample(v2, 2, "exposure", seed=0)[0]
    files = np.asarray(real.meta["file"])
    for f in np.unique(files):
        rows = np.flatnonzero(files == f)
        source = np.unique(real.meta["source_file"][rows])
        assert len(source) == 1
        src = np.flatnonzero(np.asarray(v2.meta["file"]) == source[0])
        np.testing.assert_array_equal(real.r[rows], v2.r[src])
        np.testing.assert_array_equal(real.sigma[rows], v2.sigma[src])


def test_dict_keeps_vis2_and_t3_together(residuals):
    real = block_resample(residuals, 1, "exposure", seed=3)[0]
    assert (real["vis2"].meta["source_file"][0] ==
            real["t3phi"].meta["source_file"][0])  # fmt: skip
    f = real["vis2"].meta["file"][0]
    assert real["t3phi"].meta["file"][np.flatnonzero(real["t3phi"].meta["file"] == f)[0]] == f


def test_chunk_blocks_are_whole_chunks_of_one_exposure(residuals):
    v2 = residuals["vis2"]
    real = block_resample(v2, 1, "chunk", chunk_channels=10, seed=4)[0]
    assert real.r.shape == v2.r.shape
    files = np.asarray(v2.meta["file"])
    by_file = {f: np.flatnonzero(files == f) for f in np.unique(files)}
    rows = by_file[files[0]]
    for lo in range(0, WAVE.size, 10):
        chunk = real.r[rows, lo : lo + 10]
        matches = [
            g for g, idx in by_file.items() if np.array_equal(chunk, v2.r[idx, lo : lo + 10])
        ]
        assert matches  # a chunk of some exposure, all its rows together
    assert not np.array_equal(real.r, v2.r)


def test_block_resample_rejects_bad_block(residuals):
    with pytest.raises(ValueError, match="block must be"):
        block_resample(residuals, 1, "row")


def test_block_resample_pools_by_stations_and_warns_when_alone(residuals):
    v2 = residuals["vis2"]
    stations = np.asarray(v2.meta["stations"]).copy()
    first = np.flatnonzero(np.asarray(v2.meta["file"]) == v2.meta["file"][0])
    stations[first] = stations[first][::-1]  # a different station order: its own layout
    meta = {**v2.meta, "stations": stations}
    odd = type(v2)(v2.kind, v2.wave, v2.r, v2.sigma, v2.valid, meta)
    with pytest.warns(UserWarning, match="no other exposure"):
        real = block_resample(odd, 1, "exposure", seed=0)[0]
    f = v2.meta["file"][0]
    rows = np.flatnonzero(np.asarray(real.meta["file"]) == f)
    assert set(real.meta["source_file"][rows]) == {f}


def test_reference_subsets_make_more_nulls(exposures, tmp_path):
    """Four stars a night: three others per exposure, three distinct pairs of them."""
    kw = {"screen": False, "pair_factor": None, "vis": False}
    one = null_datasets(exposures, tmp_path / "one", INS, **kw)
    assert len(one.files) == 12 and all(m["subset"] == 0 for m in one.manifest)
    two = null_datasets(exposures, tmp_path / "two", INS, reference_subsets=2, **kw)
    assert len(two.files) == 24
    again = null_datasets(exposures, tmp_path / "again", INS, reference_subsets=2, **kw)
    assert [m["references"] for m in two.manifest] == [m["references"] for m in again.manifest]
    allsets = null_datasets(exposures, tmp_path / "all", INS, reference_subsets=9, **kw)
    assert len(allsets.files) == 36
    for m in allsets.manifest:
        assert len(m["references"]) == 2 and m["target"] not in m["references"]
    by_source = {}
    for m in allsets.manifest:
        by_source.setdefault(m["source"], set()).add(tuple(m["references"]))
    assert all(len(v) == 3 for v in by_source.values())
    for f, m in zip(allsets.files, allsets.manifest, strict=True):
        assert f.endswith(f"_null_s{m['subset']}.fits")


def test_reference_subsets_need_enough_stars(exposures, tmp_path):
    with pytest.raises(ValueError, match="enough other calibrators"):
        null_datasets(exposures, tmp_path, INS, reference_subsets=2, subset_size=4, screen=False)
    with pytest.raises(ValueError, match="at least 1"):
        null_datasets(exposures, tmp_path, INS, reference_subsets=0)


# --- resampled nulls as OIFITS ----------------------------------------------


def test_resampled_nulls_have_the_layout_of_null_datasets(exposures, residuals, tmp_path):
    got = resampled_null_datasets(residuals, 2, tmp_path / "rs", INS, seed=3)
    assert len(got.files) == 2 * len(exposures)
    real = null_datasets(exposures, tmp_path / "real", INS, screen=False, pair_factor=None,
                         vis=False)  # fmt: skip
    with fits.open(real.files[0]) as h:
        want = [x.name for x in h if x.name != "OI_VIS"]
    src = {str(e.path): e for e in exposures}
    for m, path in zip(got.manifest, got.files, strict=True):
        with fits.open(path) as h:
            assert [x.name for x in h] == want
            assert "seed=3" in "".join(str(x) for x in h[0].header["HISTORY"])
        out, e = read_exposure(path), src[m["source"]]
        for a, b in ((out.v2[INS], e.v2[INS]), (out.t3[INS], e.t3[INS])):
            assert a.data.shape == b.data.shape
            np.testing.assert_array_equal(a.data["STA_INDEX"], b.data["STA_INDEX"])
        np.testing.assert_allclose(out.wave[INS], WAVE)
        assert read_flux(out, INS).shape[1] == WAVE.size
        assert m["seed"] == 3 and m["block"] == "exposure" and m["sources"]
    listed = json.loads((tmp_path / "rs" / "null_manifest.json").read_text())
    assert listed["files"][0]["realisation"] == 0


def test_resampled_nulls_carry_the_residuals(tmp_path):
    """Data are the slot's disk times exp(r), and the donor's closure phases, flags included."""
    (tmp_path / "in").mkdir()
    exps = calibrator_nights(tmp_path / "in", n_nights=3, per_night=4, diameter=0.5, disk=True)
    res = calibrator_residuals(exps, INS)
    got = resampled_null_datasets(res, 1, tmp_path / "out", INS, seed=5)
    (real,) = block_resample(res, 1, "exposure", seed=5)
    for m in got.manifest:
        out = read_exposure(m["file"])
        for kind in ("vis2", "t3phi"):
            idx = np.flatnonzero(np.asarray(real[kind].meta["file"]) == m["source"])
            ok = real[kind].valid[idx]
            if kind == "t3phi":
                d = out.t3[INS].data
                np.testing.assert_allclose(d["T3PHI"][ok], real[kind].r[idx][ok])
                np.testing.assert_allclose(d["T3PHIERR"][ok], real[kind].sigma[idx][ok])
            else:
                d = out.v2[INS].data
                model = uniform_disk_v2(out.v2[INS], WAVE, m["diameter_mas"])
                assert np.all(model < 1)  # a real disk, so a wrong diameter would show
                np.testing.assert_allclose(
                    d["VIS2DATA"][ok], (model * np.exp(real[kind].r[idx]))[ok]
                )
                np.testing.assert_allclose(
                    d["VIS2ERR"][ok], d["VIS2DATA"][ok] * real[kind].sigma[idx][ok]
                )
            np.testing.assert_array_equal(d["FLAG"], ~ok)


def test_resampled_nulls_refuse_excluded_stars_and_record_pools(residuals, tmp_path):
    with pytest.raises(ValueError, match="excluded star"):
        resampled_null_datasets(residuals, 1, tmp_path / "a", INS, exclude=["STAR0_1"])
    got = resampled_null_datasets(residuals, 1, tmp_path / "b", INS, exclude=["NOT_THERE"])
    assert all(m["donor_pool"] >= 1 for m in got.manifest)
    assert got.manifest[0]["excluded"] == ["NOT_THERE"]


def test_subset_draw_does_not_depend_on_the_directory(exposures, tmp_path):
    import shutil

    moved = tmp_path / "elsewhere"
    moved.mkdir()
    copies = []
    for e in exposures:
        shutil.copy(e.path, moved)
        copies.append(read_exposure(str(moved / os.path.basename(str(e.path)))))
    kw = {"screen": False, "pair_factor": None, "vis": False, "reference_subsets": 2}
    a = null_datasets(exposures, tmp_path / "a", INS, **kw)
    b = null_datasets(copies, tmp_path / "b", INS, **kw)
    key = lambda m: (os.path.basename(m["source"]), m["subset"], tuple(m["references"]))  # noqa: E731
    assert sorted(map(key, a.manifest)) == sorted(map(key, b.manifest))


def test_resampled_nulls_chunk_and_errors(residuals, tmp_path):
    got = resampled_null_datasets(residuals, 1, tmp_path, INS, block="chunk", chunk_channels=4)
    assert got.manifest[0]["chunk_channels"] == 4 and got.manifest[0]["sources"] == []
    with pytest.raises(ValueError, match="both"):
        resampled_null_datasets({"vis2": residuals["vis2"]}, 1, tmp_path, INS)
    with pytest.raises(ValueError, match="block must be"):
        resampled_null_datasets(residuals, 1, tmp_path, INS, block="night")


def test_null_errors_default_to_the_science_settings(exposures, nulls, tmp_path):
    """Nulls are calibrated with calibrate_exposure's defaults (loo), as science is."""
    propagated = null_datasets(exposures, tmp_path, INS, method="gaussian", error="propagated")
    explicit = null_datasets(exposures, tmp_path / "e", INS, method="linear", error="loo")
    for path, wide, narrow, same in zip(
        nulls.files, nulls.files, propagated.files, explicit.files, strict=True
    ):
        err = lambda f: read_exposure(f).v2[INS].data["VIS2ERR"]  # noqa: E731
        np.testing.assert_allclose(err(path), err(same))
        assert np.nanmean(err(wide)) > np.nanmean(err(narrow))


def test_manifest_lists_reference_diameter_errors(nulls):
    with open(nulls.files[0].rsplit("/", 2)[0] + "/null_manifest.json") as fh:
        entry = json.load(fh)["files"][0]
    assert entry["reference_diameter_errs_mas"]
