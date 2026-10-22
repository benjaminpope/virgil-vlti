"""Calibrated OIFITS to a virgil ``OIData`` with its nuisance modes attached.

.. note::
   This is a **bridge**. The modes travel as ``.npz`` sidecars whose units and
   layout are conventions (below); the architecture review
   (``design/architecture_review_2026-10.md``) proposes a versioned
   ``virgil.modes/1`` OIFITS extension that would carry them, with units and
   layout declared in the file, and this loader would then read that.

:func:`to_oidata` is the one step from the files ``virgil-vlti-calibrate``
writes (and their sidecars) to data that virgil fits with gains and closure
offsets marginalised. virgil never imports this package, so the hand-off is
arrays: one value per sample of the ``OIData`` (gains) or per closure phase
(offsets).

Units, by sidecar (converted here, once, explicitly):

* V² modes (template ``vis2-*`` modes, and the ``rel-`` diameter modes of
  ``.modes.npz``): d log|V| per 1σ, which is what ``with_gains(modes=)``
  takes. Absolute diameter modes (dV², without the ``rel-`` prefix) are not
  model-relative and are ignored, with a warning.
* transfer-function modes (``tf-vis2-*`` in ``.tf_modes.npz``, written by
  ``virgil-vlti-calibrate --tf-modes``, which takes their share out of the
  diagonal errors): d log|V| per 1σ and, for ``tf-t3phi-*``, degrees per 1σ.
  They are shared by name (one parameter per night). virgil's closure offsets
  are per frame, so ``tf-t3phi-*`` modes, which span a night, are ignored
  with a warning.
* ``t3phi`` template modes: degrees per 1σ, converted to radians for
  ``with_closure_offsets(modes=)``.
* files with ``CPONLY`` in the header (closure phases only, no usable calibrator:
  ``virgil-vlti-calibrate --cp-only-when-uncalibrated``) have V² and T3AMP flagged, so
  they give no visibility samples; any V² modes beside them are dropped.
* ``t3amp`` and ``visamp`` modes: virgil's likelihood does not use these
  observables here; they are ignored, with a warning.

Layout: every mode is a table (rows × channels, the file's order), placed on
the samples of virgil's data by (row, channel), and the placement is checked
against the data (baseline coordinates and wavelength for V², closure-phase
values for T3PHI), so that a change of virgil's ordering fails loudly.
"""

import os
import warnings

import numpy as np

from .products import read_exposure
from .templates import Template, template_modes
from .write import CP_ONLY_KEY, SCHEMA, SCHEMA_KEY

# Mode-name prefixes: the model-relative diameter modes (cli.RELATIVE_PREFIX), the
# transfer-function modes (calibrate.transfer_function_modes) and the template modes
# (templates.template_modes).
RELATIVE = "rel-"
TF = "tf-"
TEMPLATE = {"vis2": "vis2-", "t3phi": "t3phi-"}
IGNORED = ("t3amp", "visamp")


def _stem(path):
    return os.path.splitext(path)[0]


def _npz(path, insname):
    """``{observable: {mode: array}}`` of one ``.npz`` written by ``write_modes``."""
    out = {}
    with np.load(path) as d:
        # Files from before 0.1.0 have no schema: read as version 1.
        if SCHEMA_KEY in d.files:
            schema = str(d[SCHEMA_KEY])
            family, _, number = schema.rpartition("/")
            if family != SCHEMA.rpartition("/")[0] or not number.isdigit():
                raise ValueError(f"{path}: unknown modes schema {schema!r} (expected {SCHEMA!r})")
            if int(number) > int(SCHEMA.rpartition("/")[2]):
                raise ValueError(
                    f"{path}: modes schema {schema!r} is newer than this virgil-vlti "
                    f"understands ({SCHEMA!r}); upgrade virgil-vlti"
                )
        prefix = f"{insname}/"
        for key in d.files:
            if key.startswith(prefix):
                obs, _, mode = key[len(prefix) :].partition("/")
                out.setdefault(obs, {})[mode] = d[key]
    return out


def _file_modes(path, exposure, insname, template, use_template, strict):
    """The V² and closure-phase modes of one calibrated file, and what was ignored.

    From ``template`` if given, else from ``<file>.template_modes.npz`` (both
    skipped when ``use_template`` is False); and from ``<file>.modes.npz`` the
    ``rel-`` diameter modes and any template modes (``vis2-``, ``t3phi-``).
    Template modes are keyed by the file's full path (a parameter of one file
    and frame, never shared); diameter modes by their name (one calibrator on
    one night is one parameter across files). Returns ``(found, ignored)``,
    ``found`` ``None`` if the file has none.
    """
    full = os.path.realpath(path)
    found = {"vis2": {}, "t3phi": {}}
    ignored = []
    sources = []
    if template is not None:
        sources.append(
            {
                obs: {
                    m: by_obs[obs]
                    for m, by_obs in template_modes(
                        template, exposure, insname, strict=strict
                    ).items()
                    if obs in by_obs
                }
                for obs in ("vis2", "t3phi")
            }
            if use_template
            else {}
        )
    elif os.path.exists(_stem(path) + ".template_modes.npz"):
        sources.append(_npz(_stem(path) + ".template_modes.npz", insname) if use_template else {})
    for suffix in (".modes.npz", ".tf_modes.npz"):
        if os.path.exists(_stem(path) + suffix):
            sources.append(_npz(_stem(path) + suffix, insname))
    for source in sources:
        for obs, by_mode in source.items():
            for name, values in by_mode.items():
                if obs in IGNORED:
                    ignored.append(f"{obs}/{name}")
                elif obs in found and (
                    name.startswith(TF) or (obs == "vis2" and name.startswith(RELATIVE))
                ):
                    found[obs][name] = values  # shared by name: one parameter per night
                elif obs in found and name.startswith(TEMPLATE[obs]):
                    if use_template:
                        found[obs][f"{full}::{name}"] = values
                else:
                    ignored.append(f"{obs}/{name}")
    found = {obs: by_mode for obs, by_mode in found.items() if by_mode}
    return found or None, ignored


def _cells(shape):
    """(row, channel) of every cell of a table, in the order virgil lays its samples out."""
    return np.unravel_index(np.arange(int(np.prod(shape))), shape)


def _check_layout(path, record, exposure, insname):
    """Assert that virgil's samples are the file's cells, row by row, as the modes assume."""
    v2 = exposure.v2[insname].data
    wave = np.asarray(exposure.wave[insname], float)
    rows, chans = _cells(v2["VIS2DATA"].shape)
    n = rows.size
    wavel = np.broadcast_to(np.asarray(record["wavel"], float), np.shape(record["u"]))
    ok = (
        np.size(record["u"]) >= n
        and np.allclose(np.asarray(record["u"])[:n], v2["UCOORD"][rows], rtol=1e-5, atol=1e-6)
        and np.allclose(np.asarray(record["v"])[:n], v2["VCOORD"][rows], rtol=1e-5, atol=1e-6)
        and np.allclose(wavel[:n], wave[chans], rtol=1e-6)
    )
    t3 = exposure.t3[insname].data
    if ok and np.size(record["phi"]):
        rows, chans = _cells(t3["T3PHI"].shape)
        phi = np.deg2rad(np.asarray(t3["T3PHI"], float)[rows, chans])
        ok = np.size(record["phi"]) == phi.size and np.allclose(
            record["phi"], phi, rtol=1e-4, atol=1e-6, equal_nan=True
        )
    if not ok:
        raise RuntimeError(
            f"{os.path.basename(path)}: virgil's samples are not laid out as the file's tables "
            "(rows × channels); the modes cannot be aligned"
        )


def _stack(per_file, observable, usable, blocks=None):
    """Modes of several files in one array: ``(n_mode, sum(blocks))``.

    File i's modes are tables (rows × channels); ``usable[i]`` marks, per cell
    in row-major order, those virgil keeps (the rest get zero). They fill the
    start of the file's room ``blocks[i]`` (default: its cells). Same-named
    modes are one parameter.
    """
    blocks = [m.size for m in usable] if blocks is None else blocks
    names = sorted({m for found in per_file if found for m in found.get(observable, {})})
    offsets = np.concatenate([[0], np.cumsum(blocks)])
    modes = np.zeros((len(names), offsets[-1]))
    for i, found in enumerate(per_file):
        for name, values in (found or {}).get(observable, {}).items():
            values = np.asarray(values, float)
            if values.ndim != 2 or values.size != usable[i].size:
                raise ValueError(
                    f"mode {name!r} has shape {values.shape}, not the {observable} table's "
                    f"({usable[i].size} cells): were the modes written for this file?"
                )
            rows, chans = _cells(values.shape)
            modes[names.index(name), offsets[i] : offsets[i] + values.size] = np.where(
                usable[i], values[rows, chans], 0.0
            )
    return names, modes


def to_oidata(
    paths,
    insname,
    *,
    ranges=None,
    exclude=None,
    observables=None,
    template=None,
    use_template=True,
    generic_gains=None,
    cp_offsets=None,
    strict=True,
    target=None,
    extras=(),
):
    """Calibrated OIFITS files as a virgil ``OIData`` with gains and closure offsets.

    The files are read with virgil, restricted with ``OIData.select``, and
    given the nuisance modes written beside them (see the module docstring
    for units and layout): ``with_gains(modes=...)`` from the template's V²
    modes and the calibrators' model-relative diameter modes, and
    ``with_closure_offsets(modes=...)`` from the template's closure-phase
    modes. They are aligned sample for sample with the selected data and
    are zero on flagged samples (flagged V² samples stay in ``u``, flagged
    closure phases are dropped). External modes are virgil's
    ``modes`` group, with the widths in the data (1σ each); they are never
    mixed into its telescope and baseline groups. Files without modes
    contribute zeros (with a warning); if none has modes, the result is the plain
    ``OIData`` (or has ``generic_gains``), again with a warning.

    Parameters
    ----------
    paths : str, path-like, or list of them
        Calibrated ``.fits`` files (the sidecars are found by replacing the
        extension).
    insname : str
        The instrument table to read, e.g. ``"GRAVITY_SC"``.
    ranges, exclude : sequence of (lo, hi), optional
        Wavelength windows to keep, and to drop (metres), as for
        ``OIData.select``.
    observables : str or sequence of str, optional
        ``"vis"``, ``"phi"`` (and virgil's extras) to keep, as for
        ``OIData.select``. Gains go with the visibilities, offsets with the
        closure phases.
    template : Template or str, optional
        A calibration template (or its name or path), laid out on each file
        (:func:`~virgil_vlti.templates.template_modes`) instead of reading
        ``.template_modes.npz``. For a flat one from a pair test
        (:func:`~virgil_vlti.templates.gain_template_from_pairs`, e.g.
        PIONIER) this puts its widths in the data.
    use_template : bool
        False leaves the template's modes out (``template=`` and the
        ``.template_modes.npz`` and ``vis2-``/``t3phi-`` modes of
        ``.modes.npz``); the diameter modes stay.
    generic_gains : dict, optional
        Widths for virgil's own gains (``telescope``, ``baseline``,
        ``chromatic``), a fallback for data **without** V² modes. When any
        file has V² modes they are not added (a template already holds the
        telescope and baseline gains, which would be counted twice), and
        files without modes then get none.
    cp_offsets : dict, optional
        Widths (radians) of virgil's own closure-phase offsets
        (``baseline``, ``triangle``), added to the modes. Give them only if
        the template does not already hold closure offsets.
    strict : bool
        With ``template``, refuse files in another configuration (see
        ``template_modes``).
    target, extras
        As for ``virgil.oifits.read_oifits``.

    Returns
    -------
    virgil.oidata.OIData
    """
    import inspect

    import virgil.oidata as virgil_oidata
    from virgil.observables import in_ranges
    from virgil.oidata import OIData
    from virgil.oifits import read_oifits

    if not hasattr(virgil_oidata, "_valid_samples") or (
        "ranges" not in inspect.signature(OIData.select).parameters
    ):
        raise ImportError(
            "to_oidata needs OIData.select(ranges=...) and virgil.oidata._valid_samples, which "
            "virgil-astro 0.3.0 lacks: install virgil-astro>=0.4.0 (until it is released, "
            "virgil's main branch: pip install git+https://github.com/benjaminpope/virgil)"
        )

    if isinstance(paths, (str, os.PathLike)):
        paths = [paths]
    paths = [os.fspath(p) for p in paths]
    if isinstance(template, str):
        template = Template.load(template)
    # One record per file gives each file's block of samples; the joint read
    # numbers the frames across files.
    records = [read_oifits(p, target=target, insname=insname) for p in paths]
    record = read_oifits(paths, target=target, insname=insname, extras=extras)
    data = OIData(record)
    selected = data
    window = ranges is not None or exclude is not None
    if window or observables is not None:
        selected = data.select(
            ranges=ranges if ranges is not None else [(-np.inf, np.inf)],
            exclude=exclude,
            observables=observables,
        )
    exposures = [read_exposure(p, insname=insname) for p in paths]
    found, ignored = [], []
    for p, record_i, exposure in zip(paths, records, exposures, strict=True):
        _check_layout(p, record_i, exposure, insname)
        modes, skipped = _file_modes(p, exposure, insname, template, use_template, strict)
        if exposure.header.get(CP_ONLY_KEY, False):  # closure phases only: V² is not to be used
            modes = {obs: m for obs, m in (modes or {}).items() if obs != "vis2"} or None
        found.append(modes)
        ignored += skipped
    if ignored:
        warnings.warn(
            "ignored modes that are not model-relative V² or closure-phase modes: "
            f"{sorted(set(ignored))}",
            stacklevel=2,
        )

    names = (
        {"vis", "phi", *extras}
        if observables is None
        else {observables}
        if isinstance(observables, str)
        else set(observables)
    )
    out, has_gains = selected, False
    if not any(found):
        warnings.warn(
            f"no nuisance modes beside {paths[:3]}{'…' if len(paths) > 3 else ''} "
            "(.template_modes.npz, .modes.npz, .tf_modes.npz): the systematics are not "
            "marginalised, only the calibrated errors are used",
            stacklevel=2,
        )
    if any(found):
        missing = [p for p, f in zip(paths, found, strict=True) if not f]
        if missing:
            warnings.warn(f"no nuisance modes for {missing}: they get zero modes", stacklevel=2)
        # The samples select kept: it keeps by wavelength alone (not by observable).
        wavel = np.broadcast_to(np.asarray(data.wavel), np.shape(data.u))
        keep = np.ones(wavel.shape, bool)
        if window:
            keep = in_ranges(wavel, ranges if ranges is not None else [(-np.inf, np.inf)])
            if exclude is not None:
                keep &= ~in_ranges(wavel, exclude)
        if keep.sum() != np.size(selected.u):
            raise RuntimeError("the selection differs from virgil's: cannot align the modes")
        if "vis" in names:
            # Flagged V² samples stay in ``u`` (their baselines may serve closure phases),
            # so they have a place in the modes: zero.
            usable = [_usable(r["vis"], r["d_vis"], r.get("vis_flag")) for r in records]
            # Each file's room is all its samples: virgil may append closure legs to a table's.
            _, vis_modes = _stack(found, "vis2", usable, [np.size(r["u"]) for r in records])
            vis_modes = _nonzero(vis_modes[:, keep])
            if len(vis_modes):
                out = out.with_gains(modes=vis_modes)
                has_gains = True
        if "phi" in names and out.cp_flag:
            usable = [_usable(r["phi"], r["d_phi"], r.get("phi_flag")) for r in records]
            _, cp_modes = _stack(found, "t3phi", usable)
            # The closure phases virgil keeps, and of those the ones in the selection
            # (a triangle goes with its first leg). Degrees in the file, radians in virgil.
            kept = np.flatnonzero(np.concatenate(usable))
            inside = keep[np.asarray(data.i_cps1)]
            cp_modes = np.deg2rad(cp_modes[:, kept])[:, inside]
            # virgil's closure offsets are per frame: a mode that spans frames (a night-long
            # transfer-function mode) cannot be marginalised there.
            frame = np.asarray(data.frame)[np.asarray(data.i_cps1)][inside]
            spans = [len(np.unique(frame[m != 0])) > 1 for m in cp_modes]
            if any(spans):
                warnings.warn(
                    f"{sum(spans)} closure-phase modes span several frames (virgil's closure "
                    "offsets are per frame) and are ignored",
                    stacklevel=2,
                )
                cp_modes = cp_modes[~np.array(spans)]
            cp_modes = _nonzero(cp_modes)
            if len(cp_modes):
                out = out.with_closure_offsets(modes=cp_modes, **(cp_offsets or {}))
                cp_offsets = None
    if generic_gains and has_gains:
        warnings.warn(
            "V² modes present: generic_gains not added (they would be counted twice)", stacklevel=2
        )
    elif generic_gains and "vis" in names:
        out = out.with_gains(**generic_gains)
    if cp_offsets and "phi" in names and out.cp_flag:
        out = out.with_closure_offsets(**cp_offsets)
    return out


def _nonzero(modes):
    """The modes that are not zero everywhere (virgil refuses those)."""
    return modes[np.any(modes != 0, axis=1)]


def _usable(values, errors, flag):
    """Which samples virgil keeps: finite value and positive error, and unflagged if ``flag``."""
    from virgil.oidata import _valid_samples

    usable = _valid_samples(values, errors)
    return usable if flag is None else usable & ~np.asarray(flag, bool).reshape(-1)
