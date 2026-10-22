"""Build the golden synthetic GRAVITY night, and (when run) write its outputs.

Guard G2 of design/matisse_plan.md, §4: ``virgil-vlti-calibrate`` on a
synthetic GRAVITY night (two calibrators, one science file, ``--tf-modes``,
``--tf-error loo``). Every array of the ``.fits`` and ``.npz`` outputs and
every cell of ``index.csv`` is collected into one flat dictionary, which
``tests/test_gravity_golden.py`` compares with ``gravity_golden.npz``.

Regenerate the golden file only deliberately, and give the reason in the
PR::

    .venv/bin/python tests/golden/make_gravity_golden.py

The first version was written from origin/main at 4e8e808 (no package code
changed since 11b4ada, the commit the plan names), before any MATISSE code.
The transfer-function modes are written with ``--tf-mode-shape raw``, the
only shape there was then; the default shape since (``poly1``, and ``flat``) changes only
VIS2ERR, T3PHIERR and the ``.tf_modes.npz`` arrays, which
``tests/test_gravity_golden.py`` also checks.
"""

import csv
import pathlib
import sys
import tempfile

import numpy as np
from astropy.io import fits

from virgil_vlti.cli import main
from virgil_vlti.synthetic import add_systematics, write_product

GOLDEN = pathlib.Path(__file__).with_name("gravity_golden.npz")
MJD0 = 60400.1
WAVE = np.linspace(2.0e-6, 2.4e-6, 40)
MAS = np.pi / 180 / 3600 / 1000


def _night(folder):
    """Two calibrators bracketing one science exposure, in dual-field layout."""

    def transfer(t, w):
        return (0.8 - 0.3 * (t - MJD0)) * (w / 2.2e-6) ** 0.5

    def gaussian_v2(u, v, w):
        rho = np.hypot(u, v) / w
        return np.exp(-((np.pi * 2.0 * MAS * rho) ** 2) / (2 * np.log(2)))

    def write(name, t, cal):
        def v2(u, v, w, tk):
            return (1.0 if cal else gaussian_v2(u, v, w)) * transfer(tk, w)

        def t3phi(u1, v1, u2, v2_, w, tk):
            science = 0.0 if cal else 20.0 * np.sin(u1 / 40.0)
            return science + 4.0 + 3.0 * (w - 2.2e-6) / 0.2e-6

        return write_product(
            folder / f"{name}.fits",
            mjd=t,
            wave=WAVE,
            vis2=v2,
            t3phi=t3phi,
            sobj=(0.0, 0.0) if cal else (210.0, -710.0),
            catg="DUAL_CAL_VIS" if cal else "DUAL_SCI_VIS",
            target="CALSTAR" if cal else "SCISTAR",
        )

    rng = np.random.default_rng(20261007)
    for name, t, cal in (
        ("cal1", MJD0, True),
        ("sci", MJD0 + 0.02, False),
        ("cal2", MJD0 + 0.05, True),
    ):
        path = write(name, t, cal)
        add_systematics(path, rng=rng, white=1.0, telescope=0.02, cp_white=1.0, cp_offset=0.5)


def outputs(folder, shape="raw"):
    """Calibrate the golden night in ``folder``; every output array, by a flat key."""
    folder = pathlib.Path(folder)
    inputs, out = folder / "in", folder / "out"
    inputs.mkdir(parents=True)
    _night(inputs)
    main(
        [
            str(inputs),
            "--out",
            str(out),
            "--dual-field",
            "--diameter",
            "0.5±0.05",
            "--tf-modes",
            "--tf-mode-shape",  # raw: the mode shapes as written when the golden file was made
            shape,
            "--tf-error",
            "loo",
            "--bin",
            "4",
        ]
    )
    arrays = {}
    for path in sorted(out.rglob("*")):
        rel = path.relative_to(out).as_posix()
        if path.suffix == ".fits":
            with fits.open(path) as h:
                for k, hdu in enumerate(h):
                    if isinstance(hdu, fits.BinTableHDU):
                        for col in hdu.columns.names:
                            arrays[f"{rel}|{k}:{hdu.name}|{col}"] = np.asarray(hdu.data[col])
        elif path.suffix == ".npz":
            with np.load(path) as d:
                for key in d.files:
                    arrays[f"{rel}|{key}"] = d[key]
        elif path.suffix == ".csv":
            with open(path, newline="") as f:
                rows = list(csv.reader(f))
            arrays[f"{rel}|header"] = np.array(rows[0])
            for r, row in enumerate(rows[1:]):
                arrays[f"{rel}|row{r}"] = np.array(row)
    return arrays


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as tmp:
        result = outputs(tmp)
    target = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else GOLDEN
    np.savez_compressed(target, **result)
    print(f"{len(result)} arrays -> {target}")
