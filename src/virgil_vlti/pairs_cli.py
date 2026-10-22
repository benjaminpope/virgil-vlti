"""``virgil-vlti-pairs``: the pair test on a set of products.

Example::

    virgil-vlti-pairs data/*.fits --max-gap 30 --json pairs.json

With ``--by-star-night``, each star's frames of each night are tested on
their own (as on archival calibrators, to see how the widths vary from star
to star), and ``--json`` gets a list of results, one per star and night::

    virgil-vlti-pairs cals/ --field DUAL_OFF --by-star-night --json pairs.json
"""

import argparse
import dataclasses
import json
import os

from .gravity import field_matches, field_mode
from .pairs import pair_test
from .products import expand_paths, instrument_of, normalise_name, read_exposure


def by_star_night(exposures, insname=None, max_gap_min=30.0):
    """The pair test per (star, night): a list of dicts, groups without pairs left out.

    Each holds the star, the night, the files, the field configuration(s)
    and the :class:`~virgil_vlti.pairs.PairTest` as a dict.
    """
    groups = {}
    for e in exposures:
        groups.setdefault((normalise_name(e.target), e.night), []).append(e)
    out = []
    for (star, night), group in sorted(groups.items()):
        try:
            result = pair_test(group, insname, max_gap_min=max_gap_min)
        except ValueError:  # no pair of consecutive frames
            continue
        out.append(
            {
                "star": star,
                "night": night,
                "files": [os.path.basename(str(e.path)) for e in group],
                "field": sorted({str(field_mode(e.header)) for e in group}),
                **dataclasses.asdict(result),
            }
        )
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("inputs", nargs="+", help="product files, or folders of *.fits")
    p.add_argument("--insname", default=None, help="instrument (default: each file's only one)")
    p.add_argument("--max-gap", type=float, default=30.0, help="minutes between paired frames")
    p.add_argument(
        "--field",
        default=None,
        type=str.upper,
        choices=("SINGLE", "DUAL", "DUAL_ON", "DUAL_OFF", "DUAL_WIDE"),
        help="keep only products in this field configuration (virgil_vlti.gravity.field_mode)",
    )
    p.add_argument(
        "--by-star-night", action="store_true", help="one test per star and night, not pooled"
    )
    p.add_argument("--json", default=None, help="also write the fitted widths here")
    args = p.parse_args(argv)
    paths = expand_paths(args.inputs)
    exposures = []
    for path in paths:
        try:
            e = read_exposure(path, args.insname)
            if instrument_of(e.header) == "MATISSE":
                # neighbours differ in BCD and chopping (design/matisse_plan.md)
                raise ValueError("MATISSE pair test needs find_pairs(same=group_key), Stage 4/5")
            exposures.append(e)
        except (KeyError, ValueError) as err:  # e.g. another instrument's tables
            print(f"skipped {os.path.basename(path)}: {err}")
    if args.field:
        exposures = [e for e in exposures if field_matches(args.field, field_mode(e.header))]
    if args.by_star_night:
        result = by_star_night(exposures, args.insname, max_gap_min=args.max_gap)
        print(f"{len(result)} star-nights with pairs, of {len(exposures)} files")
        for r in result:
            fit = r["v2"]["+ baseline + telescope"]["widths"]
            print(f"{r['night']} {r['star']:24s} {r['n_pairs']:3d} pairs  "
                  + "  ".join(f"{k} {v:.4g}" for k, v in fit.items()))  # fmt: skip
        if args.json:
            with open(args.json, "w") as f:
                json.dump(result, f, indent=1, default=float)
        return
    result = pair_test(exposures, args.insname, max_gap_min=args.max_gap)
    print(result.summary())
    if args.json:
        with open(args.json, "w") as f:
            json.dump(dataclasses.asdict(result), f, indent=1, default=float)


if __name__ == "__main__":
    main()
