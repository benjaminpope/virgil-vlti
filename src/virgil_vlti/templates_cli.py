"""``virgil-vlti-template``: build a calibration template from saved residuals.

    virgil-vlti-template residuals.npz --epoch 2023-2025 \\
        --k-vis2 4 --k-t3phi 2 --config field=SINGLE resolution=MEDIUM polarisation=COMBINED \\
        --out src/virgil_vlti/templates/

writes ``gravity_single_medium_combined_2023-2025.npz`` (the name follows
:func:`~virgil_vlti.templates.template_name`; ``--name`` overrides it). The
config says what data the template describes: ``field`` is ``SINGLE``,
``DUAL`` (any dual field), ``DUAL_ON``, ``DUAL_OFF`` or ``DUAL_WIDE``
(:func:`~virgil_vlti.gravity.field_mode`), and data in another configuration
are refused when the modes are laid out.

The residuals come from ``virgil-vlti-calpca`` (its ``residuals.npz``); choose
the numbers of components from its cross-validation (``summary.json``).

Or make a template of flat gains from a pair test of the target's own data
(``virgil-vlti-pairs --json``), e.g. a gain per frame and telescope::

    virgil-vlti-template --from-pairs pairs.json --name target_gains --config field=DUAL
"""

import argparse
import json
import os

from .calpca_cli import load_residuals
from .gravity import field_matches
from .templates import build_template, gain_template_from_pairs, template_name


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("residuals", nargs="?", help="residuals.npz from virgil-vlti-calpca")
    p.add_argument(
        "--from-pairs",
        metavar="JSON",
        help="instead, flat gains from a pair test (virgil-vlti-pairs --json)",
    )
    p.add_argument(
        "--v2-model",
        default="+ baseline + telescope",
        help="with --from-pairs: the V² model whose widths to use",
    )
    p.add_argument(
        "--cp-model",
        default="+ triangle offset",
        help="with --from-pairs: the closure-phase model ('none' for no closure offsets)",
    )
    p.add_argument("--name", default=None, help="default: from --config and --epoch")
    p.add_argument("--epoch", default=None, help="years the data span, e.g. 2023-2025")
    p.add_argument("--k-vis2", type=int, default=None)
    p.add_argument("--k-t3phi", type=int, default=0)
    p.add_argument(
        "--config",
        nargs="*",
        default=[],
        metavar="KEY=VALUE",
        help="what data the template describes, e.g. field=SINGLE resolution=MEDIUM",
    )
    p.add_argument(
        "--provenance",
        nargs="*",
        default=[],
        metavar="KEY=VALUE",
        help="where the residuals came from (job ids, commits, data)",
    )
    p.add_argument("--folds", type=int, default=5)
    p.add_argument("--clip", type=float, default=5.0)
    p.add_argument("--out", default=".")
    args = p.parse_args(argv)
    pairs = lambda items: dict(item.split("=", 1) for item in items)  # noqa: E731
    if (args.residuals is None) == (args.from_pairs is None):
        p.error("give either residuals or --from-pairs")
    config = pairs(args.config)
    if "field" in config:
        config["field"] = config["field"].upper()
        try:
            field_matches(config["field"], None)  # refuses an unknown field
        except ValueError as err:
            p.error(str(err))
    if args.name is None:
        if args.epoch is None:
            p.error("give --name, or --epoch to name the template from its config")
        try:
            args.name = template_name(config, args.epoch)
        except ValueError as err:
            p.error(str(err))
    if args.from_pairs:
        with open(args.from_pairs) as f:
            result = json.load(f)
        template = gain_template_from_pairs(
            result, args.name, v2_model=args.v2_model,
            cp_model=None if args.cp_model.lower() == "none" else args.cp_model,
            config=config,
            provenance={"pairs": os.path.abspath(args.from_pairs), **pairs(args.provenance)},
        )  # fmt: skip
        os.makedirs(args.out, exist_ok=True)
        path = template.save(os.path.join(args.out, f"{args.name}.npz"))
        print(f"wrote {path}")
        print(f"  widths: {template.provenance['widths']}")
        print(f"  pair test: {template.provenance['pair_test']}")
        return
    if args.k_vis2 is None:
        p.error("--k-vis2 is required with residuals")
    template = build_template(
        load_residuals(args.residuals), args.name, config, args.k_vis2, args.k_t3phi,
        folds=args.folds, clip=args.clip,
        provenance={"residuals": os.path.abspath(args.residuals), **pairs(args.provenance)},
    )  # fmt: skip
    os.makedirs(args.out, exist_ok=True)
    path = template.save(os.path.join(args.out, f"{args.name}.npz"))
    print(f"wrote {path}")
    for key in ("vis2_heldout_widths", "vis2_telescope_share", "t3phi_heldout_widths", "valid_um"):
        print(f"  {key}: {template.provenance.get(key)}")


if __name__ == "__main__":
    main()
