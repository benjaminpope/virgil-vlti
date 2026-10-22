"""Guard G2 (design/matisse_plan.md, §4): GRAVITY outputs do not change.

A synthetic GRAVITY night is calibrated by ``virgil-vlti-calibrate``, and
every output array is compared with ``tests/golden/gravity_golden.npz``,
written by ``tests/golden/make_gravity_golden.py`` before any MATISSE code.
``allclose`` rather than exact equality, since interpolation and ``j1`` may
differ in the last bits across numpy, scipy and platforms; the exact
before/after comparison is made by hand in each refactoring PR.
"""

import importlib.util
import pathlib

import numpy as np
import pytest

GOLDEN_DIR = pathlib.Path(__file__).with_name("golden")

# Arrays added since the golden file was written, by deliberate PRs. Each is
# checked only for presence; every array of the golden file must still match.
# The file was regenerated when the transfer-function error gained a floor at the
# calibrators' statistical error (the *ERR arrays and .tf_modes.npz changed), and
# now includes #30's model-relative diameter modes (``rel-<cal>@<night>``).
ADDED = set()


def _maker():
    spec = importlib.util.spec_from_file_location(
        "make_gravity_golden", GOLDEN_DIR / "make_gravity_golden.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_gravity_outputs_unchanged(tmp_path):
    new = _maker().outputs(tmp_path)
    # The sidecars' schema and units entries (added in 0.1.0) are not part of the golden set.
    new = {k: v for k, v in new.items() if not k.endswith(("|__schema__", "|__units__"))}
    with np.load(GOLDEN_DIR / "gravity_golden.npz") as d:
        old = {k: d[k] for k in d.files}
    assert sorted(set(new) - ADDED) == sorted(old)
    assert set(new) >= ADDED
    for key, want in old.items():
        got = new[key]
        assert got.shape == want.shape, key
        if want.dtype.kind in "fc":
            np.testing.assert_allclose(got, want, rtol=1e-9, atol=1e-12, err_msg=key)
        elif ".csv|row" in key:
            # CSV cells are text; compare numbers as numbers, the rest exactly.
            for a, b in zip(got.tolist(), want.tolist(), strict=True):
                try:
                    np.testing.assert_allclose(float(a), float(b), rtol=1e-9, atol=1e-12)
                except ValueError:
                    assert a == b, key
        else:
            np.testing.assert_array_equal(got, want, err_msg=key)


@pytest.mark.parametrize("shape", ["poly1", "flat"])
def test_smoothed_tf_mode_shapes_change_only_the_tf_terms(tmp_path, shape):
    """The default (poly1) and flat shapes change the modes and their errors, nothing else."""
    raw = _maker().outputs(tmp_path / "raw")
    flat = _maker().outputs(tmp_path / shape, shape=shape)
    assert sorted(raw) == sorted(flat)
    changed = set()
    for key, want in raw.items():
        got = flat[key]
        if want.dtype.kind in "fc" and not np.allclose(got, want, rtol=1e-9, atol=1e-12):
            changed.add(key)
        elif key == "index.csv|row0":
            # Expected knock-on of the changed V2 errors: the "good V2" fraction
            # (column 6) moves 0.273 -> 0.272; every other cell is unchanged.
            assert got.tolist()[:6] == want.tolist()[:6]
            assert got.tolist()[7:] == want.tolist()[7:]
            assert float(want[6]) == 0.273 and float(got[6]) == 0.272
        elif want.dtype.kind not in "fc":
            np.testing.assert_array_equal(got, want, err_msg=key)
    assert any(".tf_modes.npz|" in k for k in changed)
    assert {k for k in changed if ".tf_modes.npz|" not in k} == {
        "2024-03-30/scistar_20240331T025248.fits|4:OI_VIS2|VIS2ERR",
        "2024-03-30/scistar_20240331T025248.fits|5:OI_T3|T3PHIERR",
    }
