"""Empirical null datasets: calibrators calibrated by the night's other calibrators.

A science fit's false-alarm rate depends on what calibration leaves behind:
correlated chromatic and telescope-based systematics, the real uv and
wavelength coverage, the real error bars. None of it is white, so a χ²
threshold from the reported errors is wrong (χ²_red of 2–10 is common). A
calibrator is a star with no companion to find, so calibrating it as if it
were a science target, with the other calibrators of its night, gives data of
exactly the science kind whose truth is "nothing there". Searching these for
companions, over the same grid, gives the distribution of the maximum Δχ² under
the null (the false-alarm probability, look-elsewhere effect included), and
injecting companions into them gives recovery and limits. See
design/null_false_alarm.md.

* :func:`null_datasets` writes one calibrated OIFITS file per calibrator
  exposure (``OI_VIS``, ``OI_VIS2``, ``OI_T3`` and ``OI_FLUX``), through
  :func:`~virgil_vlti.write.write_calibrated`, using the same calibration as
  :func:`~virgil_vlti.calpca.calibrator_residuals`. Calibrators that look
  like binaries are kept out as *references* only: their own nulls are still
  written, flagged in the manifest, so that the false-alarm rate can be
  computed with and without them.
* :func:`block_resample` makes further realisations from the residuals,
  by redrawing whole exposures or channel chunks of exposures, so that the
  correlations inside a block survive; :func:`resampled_null_datasets` writes
  them as OIFITS files of the same layout, so there can be as many as wanted
  where nights have too few calibrator stars for :func:`null_datasets`.
"""

import dataclasses
import itertools
import json
import os
import tempfile
import warnings
import zlib

import numpy as np

from . import gravity
from .calibrate import Calibrated, Observed, calibrate_exposure, uniform_disk_v2
from .calpca import Residuals, _conditions, _group_key, screen_calibrators
from .pairs import pair_test
from .products import default_insnames, normalise_name, read_exposure
from .write import write_calibrated

RESOLUTIONS = ("LOW", "MEDIUM", "HIGH")
POLARISATIONS = ("COMBINED", "SPLIT")
TELESCOPES = ("AT", "UT")
SCREEN_WINDOW_UM = (2.0, 2.45)  # as the calibrator PCA's default


@dataclasses.dataclass
class NullDatasets:
    """What :func:`null_datasets` wrote and what it left out.

    ``files`` are the calibrated OIFITS paths; ``dropped`` the calibrators
    flagged as possible binaries and kept out of the reference sets (``target``,
    ``reason`` ``"screen"`` or ``"pair_test"``, and the scores); their own
    null files are written too, with ``flagged: true`` in the manifest.
    ``skipped`` the exposures that could not be calibrated (``file``,
    ``reason``). ``manifest`` is one dict per file.
    """

    files: list
    dropped: list
    skipped: list
    manifest: list

    def summary(self):
        lines = [f"{len(self.files)} null files; {len(self.skipped)} exposures skipped"]
        lines += [f"flagged {d['target']}: {d['reason']}" for d in self.dropped]
        return "\n".join(lines)


def _passes(e, resolution, polarisation, telescope):
    h = e.header
    if resolution and str(h.get("ESO INS SPEC RES", "")).strip().upper() != resolution:
        return False
    if polarisation and str(h.get("ESO INS POLA MODE", "")).strip().upper() != polarisation:
        return False
    return not (telescope and _conditions(h)["telescope"] != telescope)


def _check_choice(value, allowed, name):
    if value is None:
        return None
    value = str(value).upper()
    if value not in allowed:
        raise ValueError(f"{name} must be one of {allowed}, not {value!r}")
    return value


def _reference_sets(others, subsets, size, rng):
    """The reference sets of one exposure: all ``others`` (``subsets=1``), or distinct star subsets.

    With ``subsets > 1``, the sets are subsets of exactly ``size`` of the
    other *stars* (every exposure of a chosen star goes in); all of them if
    there are no more than ``subsets``, otherwise ``subsets`` drawn without
    replacement. Sets share stars and so their nulls are not independent.
    """
    if subsets == 1:
        return [others]
    by_star = {}
    for c in others:
        by_star.setdefault(normalise_name(c.target), []).append(c)
    stars = sorted(by_star)
    combos = list(itertools.combinations(stars, size))
    if len(combos) > subsets:
        combos = [combos[i] for i in sorted(rng.choice(len(combos), subsets, replace=False))]
    return [[c for s in combo for c in by_star[s]] for combo in combos]


def _calibrate_all(
    exposures, insname, no_reference, min_cal, calibrate_kwargs, vis, subsets=1, size=1, seed=0
):
    """Calibrate every calibrator by its group's other stars: ``[(exposure, results, others, k)]``.

    Stars in ``no_reference`` (flagged as possible binaries) are calibrated
    like the rest but never serve as references. ``k`` numbers the reference
    sets of an exposure (:func:`_reference_sets`); the draw of the sets is
    seeded by ``seed`` and the file's base name, so repeated rounds, and the same
    products in another directory, redraw the same sets.
    """
    gone = {normalise_name(t) for t in no_reference}
    groups = {}
    for e in exposures:
        groups.setdefault(_group_key(e), []).append(e)
    done, skipped = [], []
    for group in groups.values():
        diameters = {c.target: c.pipeline_diameter[0] for c in group}
        for e in group:
            star = normalise_name(e.target)
            others = [
                c
                for c in group
                if normalise_name(c.target) != star and normalise_name(c.target) not in gone
            ]
            n_stars = len({normalise_name(c.target) for c in others})
            need = max(min_cal, size) if subsets > 1 else min_cal
            if n_stars < need:
                skipped.append({"file": str(e.path), "reason": f"{n_stars} other calibrator "
                                f"star(s) in the night's setup; need {need}"})  # fmt: skip
                continue
            rng = np.random.default_rng([seed, zlib.crc32(os.path.basename(str(e.path)).encode())])
            for k, refs in enumerate(_reference_sets(others, subsets, size, rng)):
                insnames = default_insnames(list(e.v2))
                results = []
                try:
                    for ins in insnames:
                        both = vis and all(ins in x.vis for x in (e, *refs))
                        kw = {"diameters": diameters, **calibrate_kwargs}
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            try:
                                results.append(calibrate_exposure(e, refs, ins, vis=both, **kw))
                            except ValueError:
                                if not both:
                                    raise
                                results.append(calibrate_exposure(e, refs, ins, vis=False, **kw))
                except ValueError as err:
                    skipped.append({"file": str(e.path), "reason": str(err)})
                    continue
                done.append((e, results, refs, k))
    return done, skipped


def _write_all(done, out_dir, suffixed=False):
    files, manifest = [], []
    for e, results, others, k in done:
        stem = os.path.splitext(os.path.basename(str(e.path)))[0]
        tag = f"_s{k}" if suffixed else ""
        path = os.path.join(out_dir, e.night, f"{stem}_null{tag}.fits")
        os.makedirs(os.path.dirname(path), exist_ok=True)
        stars = sorted({normalise_name(c.target) for c in others})
        diam, diam_err = e.pipeline_diameter
        history = (f"virgil_vlti.nulls: empirical null, {e.target} calibrated as a science "
                   f"target by {', '.join(stars)} (same night, setup and stations); its truth is a "
                   f"uniform disk of {diam} mas")  # fmt: skip
        write_calibrated(e, results, path, history, flux=True)
        files.append(path)
        manifest.append({"file": path, "source": str(e.path), "target": str(e.target),
                         "night": e.night, "n_cal": len(others), "references": stars,
                         "subset": k,
                         "reference_diameter_errs_mas": {
                             str(c.target): c.pipeline_diameter[1] for c in others},
                         "diameter_mas": diam, "diameter_err_mas": diam_err,
                         "diameter_source": "ESO QC TF CAL_SC DIAM",
                         "flagged": False, "flag_reasons": [],
                         "insnames": [r.insname for r in results],
                         "vis": [("visamp" in r.obs) for r in results],
                         **_conditions(e.header)})  # fmt: skip
    return files, manifest


def _model(exposure, insname):
    """A calibrator's own expectation: its uniform-disk V², closure phases zero."""
    d = exposure.pipeline_diameter[0]
    return {"vis2": uniform_disk_v2(exposure.v2[insname], exposure.wave[insname], d),
            "t3phi": np.zeros(np.shape(exposure.t3[insname].data["T3PHI"]))}  # fmt: skip


def _pair_flags(files, factor, max_gap_min):
    """Stars whose consecutive null exposures differ by more than ``factor`` × their errors.

    The pair test (:func:`~virgil_vlti.pairs.pair_test`) is run on each star's
    calibrated null products, with the star's uniform disk subtracted; its
    fitted ``white`` widths (× the reported σ) should be about 1. ``files``
    must carry *statistical* errors only (calibrated with ``error="none"``):
    consecutive exposures share their references, so the transfer-function
    error cancels in their difference, and counting it in σ would bias the
    widths low.
    """
    stars = {}
    for f in files:
        e = read_exposure(f)
        stars.setdefault((e.night, normalise_name(e.target)), []).append(e)
    flagged = {}
    for (_, star), exps in stars.items():
        if len(exps) < 2:
            continue
        scores = {}
        for ins in default_insnames(list(exps[0].v2)):
            try:
                pt = pair_test(exps, ins, max_gap_min=max_gap_min, model=_model)
            except ValueError:  # no pairs of consecutive exposures
                continue
            v2 = pt.v2["+ baseline + telescope"].widths["white"]
            cp = pt.cp["+ triangle offset"].widths["white"] if pt.cp else np.nan
            scores["v2_white"] = max(scores.get("v2_white", 0.0), float(v2))
            scores["cp_white"] = max(scores.get("cp_white", 0.0), float(np.nan_to_num(cp)))
        if scores and max(scores.values()) > factor:
            old = flagged.get(star, {})
            flagged[star] = {k: max(scores[k], old.get(k, 0.0)) for k in scores}
    return flagged


def null_datasets(
    exposures,
    out_dir,
    insname="GRAVITY_SC",
    *,
    exclude=(),
    resolution=None,
    polarisation=None,
    telescope=None,
    field=None,
    min_calibrators_per_night=2,
    reference_subsets=1,
    subset_size=None,
    subset_seed=0,
    screen=True,
    screen_factor=5.0,
    screen_window_um=SCREEN_WINDOW_UM,
    pair_factor=3.0,
    pair_max_gap_min=30.0,
    vis=True,
    **calibrate_kwargs,
):
    """Calibrated calibrator-on-calibrator OIFITS files: the empirical null.

    Each calibrator exposure is calibrated as a science target by the other
    stars of its group (night, setup, stations and field, as in
    :func:`~virgil_vlti.calpca.calibrator_residuals`, by
    :func:`~virgil_vlti.calibrate.calibrate_exposure`), and written with
    :func:`~virgil_vlti.write.write_calibrated`: ``OI_VIS2``, ``OI_T3``,
    ``OI_VIS`` (amplitude and phase, where every exposure involved has it)
    and the unchanged ``OI_FLUX``. Channels are not windowed or binned.

    Parameters
    ----------
    exposures : list of Exposure
        Calibrator products (``*_CAL_VIS``) with pipeline diameters.
    out_dir : str
        Files go to ``<out_dir>/<night>/<name>_null.fits``; a
        ``null_manifest.json`` lists them.
    insname : str
        The instrument used for screening; each exposure is calibrated and
        written in all its science instruments (:func:`default_insnames`).
    exclude : iterable of str
        Calibrator stars to leave out entirely, as targets and as references.
    resolution, polarisation, telescope : str, optional
        Keep only ``"LOW"``/``"MEDIUM"``/``"HIGH"``, ``"COMBINED"``/
        ``"SPLIT"`` and ``"AT"``/``"UT"`` products.
    field : str, optional
        As in :func:`~virgil_vlti.calpca.calibrator_residuals`.
    min_calibrators_per_night : int
        Fewest *other* calibrator stars in an exposure's night and setup for
        it to give a null (default 2, as a science target is calibrated by
        several); exposures with fewer are reported in ``skipped``.
    reference_subsets : int
        Nulls per calibrator exposure. 1 (default): calibrate by all the
        night's other stars, one null. More: each null uses a different subset of
        exactly ``subset_size`` of the other stars (all the subsets if there
        are no more than ``reference_subsets``, else that many drawn without
        replacement, seeded by ``subset_seed`` and the file name), written as
        ``<name>_null_s<k>.fits`` with ``subset`` and ``references`` in the
        manifest. Subsets of one exposure share its data and possibly stars,
        so their nulls are correlated, not independent. Exposures with fewer
        than ``subset_size`` other stars are skipped.
    subset_size : int, optional
        Stars per reference set when ``reference_subsets > 1`` (default
        ``min_calibrators_per_night``).
    subset_seed : int
        Seed of the subset draw.
    screen, screen_factor, screen_window_um
        Flag stars found by :func:`~virgil_vlti.calpca.screen_calibrators`
        as possible binaries. A flagged star is kept out of the reference
        sets but is still written as a null target, with ``flagged: true``
        and ``flag_reasons`` in the manifest: a science target is never
        screened on its own scatter, so removing such stars as targets would
        bias the tail of the null low. Compute the false-alarm rate with
        and without the flagged nulls.
    pair_factor : float or None
        Then flag stars whose consecutive exposures differ by more than this
        multiple of their *statistical* errors (the test runs on a second
        calibration with ``error="none"``, so that the transfer-function error
        shared by consecutive exposures is not counted twice;
        :func:`~virgil_vlti.pairs.pair_test`'s fitted ``white`` width, in V²
        or closure phase; about 1 for a star with no company). The worst star
        is flagged each round and the references recomputed. None skips the
        test.
    vis : bool
        Calibrate ``OI_VIS`` too when possible.
    **calibrate_kwargs
        For :func:`~virgil_vlti.calibrate.calibrate_exposure`. The defaults
        are the science run's (linear weighting, leave-one-out
        transfer-function error), so the null's error bars are as wide as a
        science target's; pass ``error=`` and ``method=`` to match a science
        run that used others. ``allow_unrecorded=True`` lets calibrators whose
        header keys are missing on either side through, with a warning.
        The references' diameter errors are not in the
        errors (they are common modes, :func:`~virgil_vlti.calibrate.diameter_modes`):
        the manifest lists them as ``reference_diameter_errs_mas``, and a null
        compared with a science fit that marginalises them is too narrow by
        that part.

    Returns
    -------
    NullDatasets
    """
    for key in ("window_um", "nbin"):
        if calibrate_kwargs.get(key) not in (None, 1):
            raise ValueError(f"null_datasets writes full-resolution files: {key} is not supported")
    if reference_subsets < 1:
        raise ValueError("reference_subsets must be at least 1")
    size = int(subset_size or max(min_calibrators_per_night, 1))
    if size < 1:
        raise ValueError("subset_size must be at least 1")
    resolution = _check_choice(resolution, RESOLUTIONS, "resolution")
    polarisation = _check_choice(polarisation, POLARISATIONS, "polarisation")
    telescope = _check_choice(telescope, TELESCOPES, "telescope")
    kw = {"strict": False, **calibrate_kwargs}  # method, error: calibrate_exposure's, as science
    gone = {normalise_name(t) for t in exclude}
    pool = [
        e
        for e in exposures
        if e.pipeline_diameter is not None
        and normalise_name(e.target) not in gone
        and _passes(e, resolution, polarisation, telescope)
        and (field is None or gravity.field_matches(field, gravity.field_mode(e.header)))
    ]
    if not pool:
        raise ValueError("No calibrator product with a pipeline diameter passes the filters")
    dropped = []
    if screen:
        _, bad = screen_calibrators(
            pool, insname, factor=screen_factor, window_um=screen_window_um, **kw
        )
        dropped += [{**b, "reason": "screen"} for b in bad]
    flagged = [d["target"] for d in dropped]
    seen = set()
    while True:
        for f in seen:  # files of an earlier round
            if os.path.exists(f):
                os.remove(f)
        done, skipped = _calibrate_all(
            pool, insname, flagged, min_calibrators_per_night, kw, vis,
            reference_subsets, size, subset_seed,
        )  # fmt: skip
        files, manifest = _write_all(done, out_dir, suffixed=reference_subsets > 1)
        seen = set(files)
        if pair_factor is None or not files:
            break
        stat = {**kw, "error": "none"}  # statistical σ only, for the pair test
        with tempfile.TemporaryDirectory() as tmp:
            sdone, _ = _calibrate_all(
                pool, insname, flagged, min_calibrators_per_night, stat, False
            )
            flags = _pair_flags(_write_all(sdone, tmp)[0], pair_factor, pair_max_gap_min)
        already = {normalise_name(t) for t in flagged}
        flags = {k: v for k, v in flags.items() if k not in already}
        if not flags:
            break
        names = {normalise_name(e.target): e.target for e in pool}
        # One star a round, the worst: a bad reference spoils the others of its night.
        star = max(flags, key=lambda k: max(flags[k].values()))
        dropped.append({"target": str(names[star]), "reason": "pair_test", **flags[star]})
        flagged.append(names[star])
    if not files:
        raise ValueError("No exposure has enough other calibrators to give a null")
    reasons = {}
    for d in dropped:
        reasons.setdefault(normalise_name(d["target"]), []).append(d["reason"])
    for m in manifest:
        m["flag_reasons"] = reasons.get(normalise_name(m["target"]), [])
        m["flagged"] = bool(m["flag_reasons"])
    with open(os.path.join(out_dir, "null_manifest.json"), "w") as fh:
        json.dump({"files": manifest, "dropped": dropped, "skipped": skipped}, fh, indent=1)
    return NullDatasets(files, dropped, skipped, manifest)


# ---------------------------------------------------------------------------
# Block resampling
# ---------------------------------------------------------------------------

BLOCKS = ("exposure", "chunk")


def _rows_of_files(res):
    files = np.asarray(res.meta["file"])
    return {f: np.flatnonzero(files == f) for f in dict.fromkeys(files)}


def _pools(kinds):
    """Rows of each exposure, the exposures, their row layouts, and the pools of equal layout."""
    rows = {k: _rows_of_files(r) for k, r in kinds.items()}
    files = list(rows[next(iter(rows))])
    for k in rows:
        if list(rows[k]) != files:
            raise ValueError("Every kind must hold the same exposures, in the same order")

    def _layout(f):
        key = []
        for k, r in kinds.items():
            idx = rows[k][f]
            stations = tuple(np.asarray(r.meta["stations"])[idx]) if "stations" in r.meta else ()
            field = str(r.meta["field"][idx[0]]) if "field" in r.meta and len(idx) else ""
            key.append((len(idx), stations, field))
        return tuple(key)

    layout = {f: _layout(f) for f in files}
    pool = {}
    for f, lay in layout.items():
        pool.setdefault(lay, []).append(f)
    return rows, files, layout, pool


def block_resample(data, n, block="exposure", chunk_channels=16, seed=0):
    """Resampled realisations of calibration residuals that keep their blocks.

    Every exposure of ``data`` is a slot to fill, and is refilled from an
    exposure drawn with replacement among those with the same row layout
    (the same field and the same station pair or triple of every row, in
    order, so baseline-length structure stays on its baselines), so each
    realisation has the shape of the original and keeps its uv coverage,
    flags layout and row metadata (``meta``, with ``source_file`` added).

    * ``block="exposure"``: the whole exposure's residuals move together,
      keeping all correlations within it (between baselines and triangles,
      frames, channels).
    * ``block="chunk"``: each contiguous run of ``chunk_channels`` channels of
      each exposure is drawn independently, so correlations survive within
      a chunk (and between baselines in it) but not between chunks: a check that
      results do not rest on long-range chromatic structure.

    Parameters
    ----------
    data : Residuals or dict of Residuals
        As from :func:`~virgil_vlti.calpca.calibrator_residuals`. For a dict
        (``vis2`` and ``t3phi``) one source is drawn for both, so an
        exposure's V² and closure phases stay together.
    n : int
        Number of realisations.
    block : {"exposure", "chunk"}
    chunk_channels : int
    seed : int

    Returns
    -------
    list
        ``n`` objects like ``data``.
    """
    if block not in BLOCKS:
        raise ValueError(f"block must be one of {BLOCKS}, not {block!r}")
    single = isinstance(data, Residuals)
    kinds = {"data": data} if single else dict(data)
    if block == "chunk" and chunk_channels < 1:
        raise ValueError("chunk_channels must be at least 1")
    first = next(iter(kinds.values()))
    nw = first.r.shape[1]
    if any(r.r.shape[1] != nw for r in kinds.values()):
        raise ValueError("Every kind must have the same channel grid")
    rows, files, layout, pool = _pools(kinds)
    lonely = [f for f in files if len(pool[layout[f]]) == 1]
    if lonely:
        warnings.warn(
            f"{len(lonely)} exposure(s) have no other exposure with the same stations and "
            "field to be resampled from, so they are not redrawn and add no spread",
            stacklevel=2,
        )
    edges = list(range(0, nw, chunk_channels)) + [nw] if block == "chunk" else [0, nw]
    spans = list(zip(edges[:-1], edges[1:], strict=True))
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        new = {
            k: {
                "r": r.r.copy(),
                "sigma": r.sigma.copy(),
                "valid": r.valid.copy(),
                "source": np.empty(len(r.r), dtype=object),
            }
            for k, r in kinds.items()
        }
        for f in files:
            candidates = pool[layout[f]]
            for lo, hi in spans:
                g = candidates[rng.integers(len(candidates))]
                for k, r in kinds.items():
                    dst, src = rows[k][f], rows[k][g]
                    for name in ("r", "sigma", "valid"):
                        new[k][name][dst, lo:hi] = getattr(r, name)[src, lo:hi]
                    new[k]["source"][dst] = g
        real = {}
        for k, r in kinds.items():
            meta = dict(r.meta)
            if block == "exposure":
                meta["source_file"] = new[k]["source"].astype(str)
            real[k] = Residuals(r.kind, r.wave, new[k]["r"], new[k]["sigma"], new[k]["valid"], meta)
        out.append(real["data"] if single else real)
    return out


# ---------------------------------------------------------------------------
# Resampled residuals as OIFITS
# ---------------------------------------------------------------------------


def _resampled_exposure(real, path, insname):
    """One null exposure: the slot exposure ``path`` with the rows of ``real`` as its data."""
    e = read_exposure(path)
    if insname not in e.wave:
        raise ValueError(f"{path} has no {insname} tables")
    diam, _ = e.pipeline_diameter
    out = {}
    for kind, res in real.items():
        rows = np.flatnonzero(np.asarray(res.meta["file"]) == path)
        if res.r.shape[1] != e.wave[insname].size or not np.allclose(
            res.wave, e.wave[insname], rtol=1e-7
        ):
            raise ValueError(
                f"{path}: the residuals' channel grid is not that of its {insname} tables; "
                "pass grid=None to calibrator_residuals for files with one wavelength grid"
            )
        out[kind] = (rows, res.r[rows], res.sigma[rows], res.valid[rows])
    rows, r, sigma, ok = out["vis2"]
    v2 = e.v2[insname]
    model = uniform_disk_v2(v2, e.wave[insname], diam)
    if len(rows) != len(v2.data):
        raise ValueError(f"{path}: {len(rows)} residual V² rows for {len(v2.data)} table rows")
    with np.errstate(over="ignore", invalid="ignore"):
        value = np.where(ok, model * np.exp(r), np.nan)
        error = np.where(ok, value * sigma, np.nan)
    vis2 = Observed(value, error, ~ok)
    rows, r, sigma, ok = out["t3phi"]
    if len(rows) != len(e.t3[insname].data):
        raise ValueError(f"{path}: {len(rows)} residual T3 rows for {len(e.t3[insname].data)}")
    phi = Observed(np.where(ok, r, np.nan), np.where(ok, sigma, np.nan), ~ok)
    nan = np.full(phi.value.shape, np.nan)
    amp = Observed(nan, nan.copy(), np.ones(nan.shape, bool))  # not resampled; written as NaN
    cal = Calibrated(insname, e.wave[insname], e.band[insname],
                     {"vis2": vis2, "t3amp": amp, "t3phi": phi}, n_cal=0)  # fmt: skip
    return e, cal


def resampled_null_datasets(
    residuals,
    n,
    out_dir,
    insname="GRAVITY_SC",
    *,
    block="exposure",
    chunk_channels=16,
    seed=0,
    exclude=(),
):
    """Block-resampled null realisations written as OIFITS files.

    Draws ``n`` realisations with :func:`block_resample` and writes each as
    one file per slot exposure, with the layout :func:`null_datasets` writes
    (``OI_ARRAY``, ``OI_TARGET``, ``OI_WAVELENGTH``, ``OI_VIS2``, ``OI_T3`` and the
    slot's unchanged ``OI_FLUX``, through :func:`~virgil_vlti.write.write_calibrated`),
    so that readers of the nulls (virgil's ``OIData``, the WR scripts) need
    no change. In each file the slot exposure supplies the uv coverage,
    times, stations and flux; the resampled residuals supply the data:
    ``VIS2DATA = V²_UD · exp(r)`` with the slot's own uniform disk (the
    truth, as in :func:`null_datasets`), ``VIS2ERR = VIS2DATA · σ``, and
    ``T3PHI = r``, ``T3PHIERR = σ`` (degrees, as in the files the residuals came
    from). Samples invalid in the donor are flagged and NaN. ``T3AMP`` is
    not resampled and is written as flagged NaN, and no ``OI_VIS`` table is
    written (as for ``vis=False`` in :func:`null_datasets`); only ``insname``
    is written.

    Parameters
    ----------
    residuals : dict
        ``{"vis2": Residuals, "t3phi": Residuals}`` from
        :func:`~virgil_vlti.calpca.calibrator_residuals`, called with the
        default ``grid`` (the first file's wavelengths) on files that share one
        wavelength grid, and without ``window_um`` (every channel is written).
        Its ``meta["file"]`` must still point at the (readable) calibrator
        products, which are the slots.
    n : int
        Number of realisations.
    out_dir : str
        Files go to ``<out_dir>/r<i>/<night>/<name>_null.fits`` (``i`` from 0)
        and ``null_manifest.json`` lists them; its entries have the keys of
        :func:`null_datasets`' manifest plus ``realisation``, ``block``,
        ``seed`` and ``sources`` (the donor exposures; for ``block="chunk"``
        the donor of each chunk varies, and ``sources`` is empty). The same
        provenance, and the source files, are in each file's primary header
        HISTORY.
    insname, block, chunk_channels, seed
        As in :func:`block_resample` (``insname`` names the tables written).
    exclude : iterable of str
        Stars that must not be in the residuals: pass the stars
        :func:`null_datasets` flagged (``[d["target"] for d in nulls.dropped]``)
        here **and** to :func:`~virgil_vlti.calpca.calibrator_residuals` as
        ``exclude=``, which keeps them out of the references as well as the
        donors. :func:`calibrator_residuals` does no screening itself, and a
        binary's companion signal in a donor, or in the references of the
        night, is real signal in the "null". A residual set that still has
        rows of an excluded star is refused. Recorded in the manifest
        (``excluded``).

    Notes
    -----
    The realisations are *not* independent: each slot is redrawn from the
    exposures of its layout pool, often few, so a run holds about as many
    independent nulls as there are distinct donors, not ``n``. Each manifest
    entry records ``donor_pool`` (the exposures its slot can be drawn from).
    Use the files for bootstrap errors and checks, and count independent
    nulls by donors and nights (design/null_false_alarm.md).

    Returns
    -------
    NullDatasets
        ``files`` and ``manifest`` are filled; ``dropped`` and ``skipped`` are empty.
    """
    for kind in ("vis2", "t3phi"):
        if kind not in residuals:
            raise ValueError(f"residuals must hold both 'vis2' and 't3phi'; {kind!r} is missing")
    gone = {normalise_name(t) for t in exclude}
    present = sorted(
        {str(t) for t in residuals["vis2"].meta["target"] if normalise_name(t) in gone}
    )
    if present:
        raise ValueError(
            f"residuals contain excluded star(s) {present}: recompute them with "
            "calibrator_residuals(..., exclude=...) so that they are neither donors nor references"
        )
    draws = block_resample(residuals, n, block=block, chunk_channels=chunk_channels, seed=seed)
    _, slots, layout, pool = _pools(residuals)
    files, manifest = [], []
    for i, real in enumerate(draws):
        for path in slots:
            rows = np.flatnonzero(np.asarray(real["vis2"].meta["file"]) == path)
            if block == "exposure":
                sources = sorted(set(np.asarray(real["vis2"].meta["source_file"])[rows]))
            else:
                sources = []
            e, cal = _resampled_exposure(real, path, insname)
            stem = os.path.splitext(os.path.basename(path))[0]
            out = os.path.join(out_dir, f"r{i}", e.night, f"{stem}_null.fits")
            os.makedirs(os.path.dirname(out), exist_ok=True)
            diam, diam_err = e.pipeline_diameter
            donors = ", ".join(os.path.basename(x) for x in sources) or "chunks of many exposures"
            history = (f"virgil_vlti.nulls.resampled_null_datasets: realisation {i} of "
                       f"{n}, block={block}"
                       + (f" ({chunk_channels} channels)" if block == "chunk" else "")
                       + f", seed={seed}; slot {os.path.basename(path)}, residuals from {donors}; "
                       f"its truth is a uniform disk of {diam} mas")  # fmt: skip
            write_calibrated(e, [cal], out, history, flux=True)
            files.append(out)
            manifest.append({"file": out, "source": path, "target": str(e.target),
                             "night": e.night, "n_cal": 0, "references": [],
                             "diameter_mas": diam, "diameter_err_mas": diam_err,
                             "diameter_source": "ESO QC TF CAL_SC DIAM",
                             "flagged": False, "flag_reasons": [], "insnames": [insname],
                             "vis": [False], "realisation": i, "block": block,
                             "chunk_channels": chunk_channels if block == "chunk" else None,
                             "seed": seed, "sources": sources,
                             "donor_pool": len(pool[layout[path]]),
                             "excluded": sorted(str(t) for t in exclude),
                             **_conditions(e.header)})  # fmt: skip
    with open(os.path.join(out_dir, "null_manifest.json"), "w") as fh:
        json.dump({"files": manifest, "dropped": [], "skipped": []}, fh, indent=1)
    return NullDatasets(files, [], [], manifest)


__all__ = ["NullDatasets", "block_resample", "null_datasets", "resampled_null_datasets"]
