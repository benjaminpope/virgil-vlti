# AGENTS.md

Guidance for AI coding agents working in this repository.

## What this is
Instrument-specific calibration and error models for VLTI data (GRAVITY calibration; a PIONIER pair test; MATISSE header handling), for use with [virgil](https://github.com/benjaminpope/virgil). Read [PLAN.md](PLAN.md) first: it sets the boundary with virgil and the order of work. The design notes are in `design/`; the GRAVITY review (`design/gravity_calibration_review.md`) is the main reference.

## Rules
- **virgil never imports this package.** Generic machinery (readers, likelihoods, nuisance modes) belongs in virgil, on a branch and PR there. Code here may import virgil only where it hands data or nuisance modes to it, and the calibration itself needs only numpy, scipy and astropy.
- **Design notes are general.** A science target appears only as a worked example. Target-specific scripts and science stay in that target's data folder.
- **No data in the repository.** `*.fits` and `data/` are git-ignored. Tests use synthetic products (`tests/synthetic.py`); checks against real data are examples that take a path.
- **Record provenance.** Every new method gets a row in [PROVENANCE.md](PROVENANCE.md): where it comes from (a paper with section and equation, a pipeline's manual or code, our own earlier code, or new here), how it is checked, and whether an independent check is possible.
- Merge, never rebase or force-push, pushed branches.
- Maintainers keep private workflow notes elsewhere.

## MATISSE work rules
MATISSE is outside the maintainer's core expertise, so its code and notes are held to a higher standard.
- **Conservative defaults.** Where a convention is ambiguous (BCD order, phase sign, station orientation, units), refuse or warn by default, and require an explicit option or a verified instrument fact to proceed. Prefer a larger, stated error to a silently small one.
- **Cite every convention and number.** Each MATISSE convention, keyword or number in code, facts tables or design notes cites the literature or ESO documentation (paper and section, equation or table; manual issue and section). An uncited one is marked unverified and listed for the MATISSE expert.
- **Expert review.** Send MATISSE conventions (BCD order, phase sign) to the MATISSE expert named in the maintainer's private notes before merging.
- **Adversarial review first.** A MATISSE PR is not marked ready until it has had an adversarial review of its MATISSE content against those sources.

## Commands
| Task | Command |
|---|---|
| Set up | `uv venv --python 3.11 .venv && uv pip install --python .venv/bin/python -e ".[dev]"` |
| Tests | `.venv/bin/python -m pytest` |
| Lint | `.venv/bin/ruff check . && .venv/bin/ruff format --check .` |

`.[dev]` installs virgil-astro>=0.4 from PyPI, which has the API `to_oidata` needs (`OIData.select(ranges=)`, `_valid_samples`), so the loader tests run. To test against virgil's main branch instead: `uv pip install --python .venv/bin/python "virgil-astro @ git+https://github.com/benjaminpope/virgil"` (CI's `virgil-main` job does this).
