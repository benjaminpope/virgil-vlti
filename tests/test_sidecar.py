"""The modes sidecar's schema and units (no virgil needed)."""

import json

import numpy as np
import pytest

from virgil_vlti.load import _npz
from virgil_vlti.write import SCHEMA, SCHEMA_KEY, UNITS_KEY, write_modes

INS = "GRAVITY_SC"


def test_schema_and_units_are_written(tmp_path):
    path = str(tmp_path / "m.npz")
    write_modes(path, {INS: {"a": {"t3phi": np.ones((2, 3))}, "rel-b": {"vis2": np.ones((2, 3))}}})
    with np.load(path) as d:
        assert str(d[SCHEMA_KEY]) == SCHEMA
        units = json.loads(str(d[UNITS_KEY]))
    assert units[f"{INS}/t3phi/a"].startswith("deg")
    assert units[f"{INS}/vis2/rel-b"].startswith("d log|V|")
    assert set(_npz(path, INS)) == {"t3phi", "vis2"}


def test_legacy_file_without_schema_is_read(tmp_path):
    path = tmp_path / "old.npz"
    np.savez(path, **{f"{INS}/vis2/rel-a": np.ones((2, 3))})
    assert list(_npz(str(path), INS)["vis2"]) == ["rel-a"]


@pytest.mark.parametrize(
    ("schema", "message"),
    [
        ("virgil_vlti.modes/99", "newer"),
        ("other.modes/1", "unknown"),
        ("virgil_vlti.modes/1.1", "unknown"),
    ],
)
def test_unreadable_schema_is_refused(tmp_path, schema, message):
    path = tmp_path / "x.npz"
    np.savez(path, **{SCHEMA_KEY: np.array(schema)})
    with pytest.raises(ValueError, match=message):
        _npz(str(path), INS)
