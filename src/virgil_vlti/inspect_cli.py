"""``virgil-vlti-inspect``: what is a VLTI product, and may virgil-vlti calibrate it?

Examples::

    virgil-vlti-inspect reduced/*.fits
    virgil-vlti-inspect reduced/*.fits --json

Prints one row per file: instrument, product category, and for MATISSE the
band, resolution, BCD position, chopping, ``mat_raw_estimates`` version,
N-band phase convention, merge evidence, ``replaceTel`` and a verdict: ``ok``,
or the reasons the product is refused as a calibration input (by the
defaults of ``virgil-vlti-calibrate``). Exit status 0 if every file is
``ok``, else 1. See design/matisse_plan.md, §4, Stage 1.
"""

import argparse
import json
import os
import sys

from astropy.io import fits

from . import matisse
from .products import instrument_of, is_calibrated, product_category

COLUMNS = (
    "instrument",
    "category",
    "band",
    "resolution",
    "bcd",
    "chopped",
    "pipeline",
    "phase_convention",
    "merged",
    "calibrated",
    "replaceTel",
    "verdict",
)


def inspect(path):
    """The :func:`virgil_vlti.matisse.describe` summary of one file (any instrument)."""
    try:
        with fits.open(path, memmap=False) as h:
            header = h[0].header.copy()
            tables = [x.data.copy() for x in h if x.name == "OI_VIS2"]
            waves = [x.data["EFF_WAVE"] for x in h if x.name == "OI_WAVELENGTH"]
    except (OSError, KeyError, ValueError) as err:
        return {"file": path, "verdict": f"refused: unreadable ({err})"}
    if instrument_of(header) == "MATISSE":
        return matisse.describe(header, tables, path, waves[0] if waves else None)
    calibrated = is_calibrated(header)
    return {
        "file": path,
        "instrument": instrument_of(header),
        "category": product_category(header),
        "calibrated": calibrated,
        "verdict": "refused: calibrated" if calibrated else "ok",
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("paths", nargs="+", help="product files")
    p.add_argument("--json", action="store_true", help="print JSON instead of a table")
    args = p.parse_args(argv)
    rows = [inspect(path) for path in args.paths]
    if args.json:
        print(json.dumps(rows, indent=1, default=str))
    else:
        for row in rows:
            print(os.path.basename(row["file"]))
            for key in COLUMNS:
                if key in row:
                    print(f"  {key:17s} {row[key]}")
            for e in row.get("merge_evidence", ()):
                print(f"  {'merge evidence':17s} {e}")
            for note in row.get("notes", ()):
                print(f"  {'warning':17s} {note}")
    ok = all(row["verdict"] == "ok" for row in rows)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
