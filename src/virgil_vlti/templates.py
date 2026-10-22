"""Calibration templates: the calibrator PCA's systematics, as nuisance modes.

A template holds the residual shapes the calibrator PCA (:mod:`virgil_vlti.calpca`)
found for one instrument configuration, with their widths, in the form virgil's
correlated calibration nuisances (its Stage 6d) are designed to take: gains on
log |V| in low-rank blocks per frame, and optionally closure-phase offsets.

**V².** Each PCA component c(λ) of log V² becomes two kinds of mode on log |V|
(½ log V²):

* a **telescope gain** per (frame, telescope): log |V_kl| += (t_k + t_l) s_T(λ);
* a **baseline gain** per (frame, baseline): log |V_kl| += b_kl s_B(λ),

with t, b ~ N(0, 1). The split follows the component's telescope fraction f
(:func:`~virgil_vlti.calpca.telescope_fraction`): with no telescope structure f
is 4/6, so the telescope share of the per-baseline variance is 3f − 2 (clipped
to [0, 1]). With W the component's held-out width (per-baseline rms
amplitude on the unit-norm shape, in log V²),
s_T = ½ W √((3f − 2)/2) c and s_B = ½ W √(3 − 3f) c.

**Closure phases.** Each component becomes an offset per (frame, triangle):
φ_t += o_t s(λ), o ~ N(0, 1), s = W c in degrees. virgil's plan prefers no
closure offsets by default (they are optional), so these are marked optional.

**Flat gains from the pair test.** :func:`gain_template` makes a template
of achromatic (flat) shapes from widths the pair test measured on the target's
own data (:func:`~virgil_vlti.pairs.pair_test`: a fractional V² gain per
frame and telescope g_T, per frame and baseline g_B, and a closure offset per
frame and triangle o). A V² gain g is d log|V| = g/2, so s_T = g_T/2,
s_B = g_B/2 and s = o.

**Configurations.** A template's ``config`` says which data it describes:
``field`` (``SINGLE``; ``DUAL`` for any dual field, or one of ``DUAL_ON``,
``DUAL_OFF``, ``DUAL_WIDE``, as :func:`~virgil_vlti.gravity.field_mode`
reads the science fibre's offset), ``resolution``, ``polarisation`` and
optionally ``telescopes`` (``AT`` or ``UT``; ``mostly-AT`` matches ATs).
:func:`template_modes` refuses data in another configuration. Shipped templates are named by
:func:`template_name`, e.g. ``gravity_single_medium_combined_2023-2025``.

Modes are **relative**: for V², what is stored is d log|V| per 1σ; a consumer
turns it into dV² = 2 V²_model · d log|V| (covariances built from the model,
not the data; Lachaume 2021). Closure-phase modes are in degrees.
"""

import dataclasses
import json
import os
import warnings
from importlib import resources

import numpy as np

from . import gravity
from .calpca import clip_rows, heldout_widths, telescope_fraction, weighted_pca
from .pairs import frames_of

FORMAT_VERSION = 1

# How each field configuration appears in a template's name.
FIELD_TOKENS = {
    "SINGLE": "single",
    "DUAL": "dual",
    "DUAL_ON": "dual-on",
    "DUAL_OFF": "dual-off",
    "DUAL_WIDE": "dual-wide",
}


def template_name(config, epoch, instrument="gravity"):
    """A template's conventional name, ``<instrument>_<field>_<resolution>_<polarisation>_<epoch>``.

    e.g. ``template_name({"field": "DUAL_OFF", "resolution": "MEDIUM",
    "polarisation": "COMBINED"}, "2023-2025")`` is
    ``gravity_dual-off_medium_combined_2023-2025``.
    """
    missing = [k for k in ("field", "resolution", "polarisation") if not config.get(k)]
    if missing:
        raise ValueError(f"A template name needs {missing} in its config")
    field = str(config["field"]).upper()
    if field not in FIELD_TOKENS:
        raise ValueError(f"Unknown field {field!r}; one of {list(FIELD_TOKENS)}")
    parts = (instrument, FIELD_TOKENS[field], config["resolution"], config["polarisation"], epoch)
    return "_".join(str(x).lower() for x in parts)


@dataclasses.dataclass
class Template:
    """Systematic shapes for one instrument configuration (see the module docstring).

    ``vis2_tel`` and ``vis2_base`` are (k, channels): the 1σ telescope- and
    baseline-gain shapes on log |V|. ``t3phi`` is (k, channels): the 1σ
    closure-phase offset shapes in degrees. ``wave`` (m) is their grid;
    ``config`` says which data they describe; ``provenance`` how they were
    made.
    """

    name: str
    wave: np.ndarray
    vis2_tel: np.ndarray
    vis2_base: np.ndarray
    t3phi: np.ndarray
    config: dict
    provenance: dict

    def save(self, path):
        np.savez_compressed(
            path,
            wave=self.wave,
            vis2_tel=self.vis2_tel,
            vis2_base=self.vis2_base,
            t3phi=self.t3phi,
            meta=json.dumps(
                {
                    "format": FORMAT_VERSION,
                    "name": self.name,
                    "config": self.config,
                    "provenance": self.provenance,
                }  # fmt: skip
            ),
        )
        return path

    @classmethod
    def load(cls, name_or_path):
        """Load a template from a path, or by name from those shipped with the package."""
        path = name_or_path
        if not os.path.exists(path):
            path = resources.files("virgil_vlti") / "templates" / f"{name_or_path}.npz"
            if not path.is_file():
                raise FileNotFoundError(f"No template {name_or_path!r}; shipped: {shipped()}")
        with np.load(path, allow_pickle=False) as d:
            meta = json.loads(str(d["meta"]))
            if meta.get("format") != FORMAT_VERSION:
                raise ValueError(f"Template format {meta.get('format')}, expected {FORMAT_VERSION}")
            return cls(meta["name"], d["wave"], d["vis2_tel"], d["vis2_base"], d["t3phi"],
                       meta["config"], meta["provenance"])  # fmt: skip

    def mismatches(self, header):
        """How a product's configuration differs from the template's (empty if none).

        Keys the template does not set are not compared. A product that does
        not record a key the template sets is a mismatch: it cannot be
        shown to match. A ``telescopes`` of ``"mostly-AT"`` (a template built
        mostly from AT data) matches AT data.
        """
        out = []
        tel = header.get("ESO ISS CONF T1NAME")
        checks = {
            "resolution": header.get("ESO INS SPEC RES"),
            "polarisation": header.get("ESO INS POLA MODE"),
            "telescopes": None if tel is None else str(tel).strip().upper()[:2],
        }
        for key, value in checks.items():
            want = self.config.get(key)
            if want is None:
                continue
            cmp = str(want).strip().upper()
            if key == "telescopes":  # "mostly-AT" matches AT
                cmp = cmp.removeprefix("MOSTLY-")
            if value is None or str(value).strip().upper() != cmp:
                out.append(f"{key}: data {value!r}, template {want!r}")
        want = self.config.get("field")
        if want is not None:
            mode = gravity.field_mode(header)
            if not gravity.field_matches(want, mode):
                out.append(f"field: data {mode!r}, template {want!r}")
        return out

    def on_grid(self, wave):
        """The template's shapes interpolated in 1/λ onto ``wave``; zero outside its range."""
        x, xg = 1 / self.wave, 1 / np.asarray(wave, float)
        order = np.argsort(x)

        def interp(shapes):
            return np.array(
                [np.interp(xg, x[order], s[order], left=0.0, right=0.0) for s in shapes]
            ).reshape(len(shapes), xg.size)

        lo, hi = self.provenance.get("valid_um", [-np.inf, np.inf])
        um = np.asarray(wave, float) * 1e6
        inside = (um >= lo) & (um <= hi)  # no modes outside the fitted window
        return tuple(interp(a) * inside for a in (self.vis2_tel, self.vis2_base, self.t3phi))


def shipped():
    """Names of the templates shipped with the package."""
    folder = resources.files("virgil_vlti") / "templates"
    return sorted(p.name[: -len(".npz")] for p in folder.iterdir() if p.name.endswith(".npz"))


def build_template(residuals, name, config, k_vis2, k_t3phi, *, folds=5, clip=5.0, provenance=None):
    """Make a template from calibration residuals.

    The residuals are as :func:`~virgil_vlti.calpca.calibrator_residuals` gives.

    Rows are clipped (:func:`~virgil_vlti.calpca.clip_rows`), ``k_vis2`` and
    ``k_t3phi`` components fitted, and their robust held-out widths and (for
    V²) telescope fractions measured, as ``virgil-vlti-calpca`` does.
    """
    provenance = dict(provenance or {})
    v2 = residuals["vis2"]
    if clip > 0:
        v2, dropped = clip_rows(v2, k=2, threshold=clip)
        provenance["vis2_rows_clipped"] = int(len(dropped))
    pca = weighted_pca(v2, k_vis2)
    if k_vis2:
        widths = heldout_widths(v2, k_vis2, folds=folds)
        f = telescope_fraction(v2, pca.amplitudes)
        if not np.all(np.isfinite(f)):
            raise ValueError(
                "No complete six-baseline frame to measure the telescope fraction: "
                "cannot split the V² components into telescope and baseline gains"
            )
    else:
        widths, f = np.zeros(0), np.zeros(0)
    share = np.clip(3 * f - 2, 0, 1)  # telescope share of per-baseline variance
    vis2_tel = 0.5 * (widths * np.sqrt(share / 2))[:, None] * pca.components
    vis2_base = 0.5 * (widths * np.sqrt(1 - share))[:, None] * pca.components
    provenance.update(
        vis2_rows=int(v2.r.shape[0]),
        vis2_heldout_widths=widths.tolist(),
        vis2_telescope_fraction=f.tolist(),
        vis2_telescope_share=share.tolist(),
    )
    cp = residuals["t3phi"]
    if k_t3phi:
        if clip > 0:
            cp, dropped = clip_rows(cp, k=min(2, k_t3phi), threshold=clip)
            provenance["t3phi_rows_clipped"] = int(len(dropped))
        cpca = weighted_pca(cp, k_t3phi)
        cw = heldout_widths(cp, k_t3phi, folds=folds)
        t3phi = cw[:, None] * cpca.components
        provenance.update(t3phi_rows=int(cp.r.shape[0]), t3phi_heldout_widths=cw.tolist())
    else:
        t3phi = np.zeros((0, v2.r.shape[1]))
    wave = np.asarray(v2.wave)
    valid = v2.valid.any(0)
    provenance["valid_um"] = [float(wave[valid].min() * 1e6), float(wave[valid].max() * 1e6)]
    # Outside the fitted window the shapes are not known: zero them.
    vis2_tel[:, ~valid], vis2_base[:, ~valid], t3phi[:, ~valid] = 0.0, 0.0, 0.0
    return Template(name, wave, vis2_tel, vis2_base, t3phi, dict(config), provenance)


def template_modes(template, exposure, insname, wave=None, *, strict=True):
    """The template's nuisance modes for one exposure, laid out like its tables.

    Returns ``{mode: {observable: array}}`` as :func:`~virgil_vlti.write.write_modes`
    takes. V² modes are d log|V| per 1σ (relative: multiply by 2 V²_model for
    dV²), closure-phase modes are degrees per 1σ. Mode names say what they
    are: ``vis2-pc1-tel28@<file>#frame0`` (a telescope gain),
    ``vis2-pc1-base28-23@…`` (a baseline gain), ``t3phi-pc1-tri28-23-18@…``
    (an optional closure offset). Modes with the same name are one parameter.

    ``wave`` is the channel grid of the calibrated tables the modes go with
    (after any window and binning); by default the exposure's own. Frames
    are grouped by time (:func:`~virgil_vlti.pairs.frames_of`), so a missing
    row cannot mix frames.

    A template measured in one configuration (field, resolution,
    polarisation, telescopes: :meth:`Template.mismatches`) says nothing about
    another, so an exposure in another one is refused (``ValueError``);
    ``strict=False`` only warns.
    """
    problems = template.mismatches(exposure.header)
    if problems:
        message = (
            f"{os.path.basename(str(exposure.path))}: template {template.name!r} is for another "
            f"configuration ({'; '.join(problems)})"
        )
        if strict:
            raise ValueError(message + "; pass strict=False to use it anyway")
        warnings.warn(message, stacklevel=2)
    wave = exposure.wave[insname] if wave is None else np.asarray(wave, float)
    lo, hi = template.provenance.get("valid_um", [-np.inf, np.inf])
    outside = (wave * 1e6 < lo) | (wave * 1e6 > hi)
    if outside.any():
        warnings.warn(
            f"{outside.sum()} channels outside the template's range {lo:.2f}–{hi:.2f} µm "
            "get no template modes",
            stacklevel=2,
        )
    tel, base, phi = template.on_grid(wave)
    tag = os.path.basename(str(exposure.path))
    modes = {}
    nv2, nt3 = len(exposure.v2[insname].data), len(exposure.t3[insname].data)
    for frame, f in enumerate(frames_of(exposure, insname)):
        telescopes = sorted({t for sta in f.v2_rows for t in sta})
        for j in range(len(tel)):
            for t in telescopes if np.any(tel[j]) else ():
                arr = np.zeros((nv2, wave.size))
                for sta, r in f.v2_rows.items():
                    if t in sta:
                        arr[r] = tel[j]
                modes[f"vis2-pc{j + 1}-tel{t}@{tag}#frame{frame}"] = {"vis2": arr}
            for (a, b), r in f.v2_rows.items() if np.any(base[j]) else ():
                arr = np.zeros((nv2, wave.size))
                arr[r] = base[j]
                modes[f"vis2-pc{j + 1}-base{a}-{b}@{tag}#frame{frame}"] = {"vis2": arr}
        for j in range(len(phi)):
            for sta, r in f.t3_rows.items():
                arr = np.zeros((nt3, wave.size))
                arr[r] = phi[j]
                name = "-".join(str(int(x)) for x in sta)
                modes[f"t3phi-pc{j + 1}-tri{name}@{tag}#frame{frame}"] = {"t3phi": arr}
    return modes


FLAT_UM = (0.5, 15.0)  # a flat template's range: every VLTI band


def gain_template(name, *, telescope=0.0, baseline=0.0, offset=0.0, config=None, provenance=None):
    """A template of flat (achromatic) gains and closure offsets.

    Parameters
    ----------
    name : str
    telescope, baseline : float
        Fractional rms of the V² gain per (frame, telescope) and per (frame,
        baseline), as the pair test reports them; a telescope gain g enters
        the baselines through that telescope.
    offset : float
        Rms closure-phase offset per (frame, triangle), degrees.
    config, provenance : dict, optional
        What data the widths describe, and where they come from.
    """
    widths = {"telescope": telescope, "baseline": baseline, "offset": offset}
    for key, value in widths.items():
        if not (np.isfinite(value) and value >= 0):
            raise ValueError(f"{key} width must be finite and non-negative, not {value!r}")
    wave = np.array(FLAT_UM) * 1e-6
    flat = np.ones((1, wave.size))
    k_v2 = 1 if telescope or baseline else 0
    return Template(
        name,
        wave,
        (flat * telescope / 2)[:k_v2],
        (flat * baseline / 2)[:k_v2],
        (flat * offset)[: 1 if offset else 0],
        dict(config or {}),
        {"kind": "flat gains", "widths": widths, "valid_um": list(FLAT_UM), **(provenance or {})},
    )


def gain_template_from_pairs(
    result, name, *, v2_model="+ baseline + telescope", cp_model="+ triangle offset", **kwargs
):
    """:func:`gain_template` with the widths of one pair-test model.

    ``result`` is a :class:`~virgil_vlti.pairs.PairTest` or its JSON form (as
    ``virgil-vlti-pairs --json`` writes it). ``v2_model`` and ``cp_model``
    choose the fitted models (None leaves that observable out). The white-noise
    factors are not modes; they are recorded in the provenance.
    """
    if not isinstance(result, dict):
        result = dataclasses.asdict(result)
    v2 = result["v2"][v2_model]["widths"] if v2_model else {}
    cp = result["cp"][cp_model]["widths"] if cp_model and result.get("cp") else {}
    provenance = {
        "pair_test": {
            "n_pairs": result.get("n_pairs"),
            "v2_model": v2_model,
            "cp_model": cp_model if cp else None,
            "v2_white": v2.get("white"),
            "cp_white": cp.get("white"),
        },
        **kwargs.pop("provenance", {}),
    }
    return gain_template(
        name,
        telescope=v2.get("telescope", 0.0),
        baseline=v2.get("baseline", 0.0),
        offset=cp.get("offset", 0.0),
        provenance=provenance,
        **kwargs,
    )
