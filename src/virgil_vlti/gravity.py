"""GRAVITY-specific conventions used by the calibration.

A GRAVITY product holds tables for the fringe tracker (``GRAVITY_FT``, a few
low-resolution channels) beside the science channel (``GRAVITY_SC``, or
``GRAVITY_SC_P1`` and ``_P2`` in split polarisation). The fringe-tracker data
are never calibrated here (design/gravity_calibration_review.md, §12.6).

In dual-field mode the science fibre sits at an offset ``ESO INS SOBJ X/Y``
(mas) from the fringe-tracking star. An offset near zero puts the science
fibre on that star, which is then a calibrator observed in the same mode,
minutes from the target, if it is unresolved.
"""

import re

import numpy as np

# Science-fibre offsets closer than this (mas) count as on the fringe-tracking
# star. GRAVITY's fibre field of view is about 60 mas on the UTs, so 50 mas is
# well inside an AT fibre and keeps every off-axis pointing out. (The fibre's
# response at larger offsets is `fibre_response`.)
ONAXIS_MAS = 50.0


def science_insnames(names):
    """The science-channel ``INSNAME``s among ``names`` (FT tables dropped)."""
    return [n for n in dict.fromkeys(names) if n.upper().startswith("GRAVITY_SC")]


def fibre_offset_mas(header):
    """The science fibre's offset from the fringe tracker, (x, y) in mas."""
    return float(header["ESO INS SOBJ X"]), float(header["ESO INS SOBJ Y"])


def fibre_separation_mas(header):
    """The science fibre's angular separation from the fringe tracker, in mas."""
    return float(np.hypot(*fibre_offset_mas(header)))


def is_on_axis(header, onaxis_mas=ONAXIS_MAS):
    """Whether the science fibre is on the fringe-tracking star."""
    return bool(fibre_separation_mas(header) < onaxis_mas)


# Beyond this fibre offset (mas) a dual-field pointing is outside the field of
# non-Wide dual field (2" on the UTs, 4" on the ATs): GRAVITY Wide.
WIDE_MAS = 4000.0

# Field configurations, as :func:`field_mode` names them. A template or a
# selection may also say just "DUAL", which accepts every dual-field one.
FIELD_MODES = ("SINGLE", "DUAL_ON", "DUAL_OFF", "DUAL_WIDE")


def field_mode(header, onaxis_mas=ONAXIS_MAS):
    """The product's field configuration: where the science fibre is.

    * ``"SINGLE"``: single field, one fibre for both;
    * ``"DUAL_ON"``: dual field, science fibre on the fringe-tracking star
      (offset below ``onaxis_mas``): an in-field calibrator, no anisoplanatism;
    * ``"DUAL_OFF"``: dual field, science fibre on another object within the
      non-Wide field (up to :data:`WIDE_MAS`);
    * ``"DUAL_WIDE"``: dual field beyond it (GRAVITY Wide);
    * ``"DUAL"`` if the fibre offset is not recorded, None if not even the
      field is.

    The field comes from ``ESO QC FIELD MODE``, else the product category
    (``SINGLE_*``/``DUAL_*``), else ``INSMODE``. Every dual-field calibrator
    product checked in the archive (2026-10-07) records the same optics
    (``ESO INS OPTI1 NAME = STAR``), so only the offset tells these apart.
    """
    field = str(header.get("ESO QC FIELD MODE", "")).strip().upper()
    if field not in ("SINGLE", "DUAL"):
        catg = str(header.get("ESO PRO CATG") or header.get("HIERARCH VIRGIL INPUT CATG") or "")
        catg = catg.strip().upper()
        insmode = str(header.get("INSMODE", "")).strip().upper().split(",")[0]
        field = next((f for f in ("SINGLE", "DUAL") if catg.startswith(f) or insmode == f), None)
    if field != "DUAL":
        return field
    try:
        sep = fibre_separation_mas(header)
    except (KeyError, TypeError, ValueError):
        return "DUAL"
    if sep < onaxis_mas:
        return "DUAL_ON"
    return "DUAL_OFF" if sep <= WIDE_MAS else "DUAL_WIDE"


def field_matches(wanted, mode):
    """Whether data in field configuration ``mode`` belong to ``wanted``.

    ``wanted`` is one of :data:`FIELD_MODES` or ``"DUAL"`` (any dual field).
    Data whose configuration is unknown (None, or ``"DUAL"`` without a fibre
    offset when ``wanted`` is a specific dual-field one) do not match.
    """
    wanted = str(wanted).strip().upper()
    if wanted not in (*FIELD_MODES, "DUAL"):
        raise ValueError(f"Unknown field {wanted!r}; one of {(*FIELD_MODES, 'DUAL')}")
    if mode is None:
        return False
    return mode == wanted or (wanted == "DUAL" and mode.startswith("DUAL"))


# Header keys that must agree between science and calibrator exposures: a
# transfer function from another spectral resolution or polarisation mode
# does not apply.
MODE_KEYS = ("ESO INS SPEC RES", "ESO INS POLA MODE")

# Pipeline (``gravity_vis``) parameters that must agree. ``vis-correction-sc``
# (VFACTOR, FORCE or NONE) changes the chromatic shape of V² (pipeline
# manual §6, §9.10), so mixing it between target and calibrator leaves a
# coloured error (design/gravity_calibration_review.md, §7c.5).
REDUCTION_PARAMS = ("vis-correction-sc",)

# Keys whose differences are only warned about. The science-channel DIT
# changes the coherence loss that the vFactor correction removes.
WARN_KEYS = ("ESO DET2 SEQ1 DIT", "ESO PRO REC1 PIPE ID")


NOT_RECORDED = "[not recorded]"


def pipeline_params(header):
    """The pipeline recipe parameters of ``gravity_vis``, by name.

    The recipe is found among the ``ESO PRO RECn ID`` cards (the one whose ID
    is ``gravity_vis``), as the card order depends on the pipeline chain.
    Where no ID card names it, ``ESO PRO REC1`` is used.
    """
    from .products import pipeline_params as params

    for key, value in header.items():
        m = re.fullmatch(r"ESO PRO REC(\d+) ID", key)
        if m and str(value).strip().lower() == "gravity_vis":
            return params(header, recipe=int(m.group(1)))
    return params(header, recipe=1)


def compatibility(sci_header, cal_header):
    """How a calibrator's header differs from the science exposure's.

    Returns three lists of messages: mode differences (which make the
    transfer function inapplicable), pipeline-setting differences (which
    bias it), and minor differences (DIT, pipeline version). A key or
    parameter absent on either side is a named outcome, ``"[not recorded]"``
    (in the message), reported with the same weight as a difference: no check
    could be made. Only when absent from both headers is a minor key skipped.
    """
    mode, reduction, minor = [], [], []

    def compare(out, label, a, b):
        if a is None or b is None:
            out.append(f"{label}: {NOT_RECORDED} (science {a!r}, calibrator {b!r})")
        elif a != b:
            out.append(f"{label}: science {a!r}, calibrator {b!r}")

    for key in MODE_KEYS:
        compare(mode, key, sci_header.get(key), cal_header.get(key))
    p, q = pipeline_params(sci_header), pipeline_params(cal_header)
    for name in REDUCTION_PARAMS:
        compare(reduction, f"pipeline {name}", p.get(name), q.get(name))
    for key in WARN_KEYS:
        a, b = sci_header.get(key), cal_header.get(key)
        if a != b:
            minor.append(f"{key}: science {a!r}, calibrator {b!r}")
    return mode, reduction, minor


# ---------------------------------------------------------------------------
# Anisoplanatism in dual-field mode
# ---------------------------------------------------------------------------

ARCSEC = np.pi / 180 / 3600
MAS = ARCSEC / 1000


def _mean_card(header, key):
    """Mean of a ``<key> START`` / ``<key> END`` pair, ignoring missing or non-positive."""
    values = [header.get(f"{key} {end}") for end in ("START", "END")]
    values = [float(v) for v in values if v is not None and float(v) > 0]
    return float(np.mean(values)) if values else None


def telescope_diameter(header):
    """Aperture diameter (m): 8.2 for the UTs, 1.8 for the ATs."""
    name = str(header.get("ESO ISS CONF T1NAME", "AT")).strip().upper()
    return 8.2 if name.startswith("UT") else 1.8


def anisoplanatic_v2_factor(wave, separation_arcsec, seeing_arcsec, theta0_arcsec, diameter_m):
    """The fraction of V² kept by off-axis fringe tracking, per wavelength.

    GRAVITY+ Collaboration 2022 (GRAVITY Wide), eqs. 1, 2, 4 and 5, after
    Elhalkouj et al. 2008 (Kolmogorov turbulence, infinite outer scale,
    small apertures): V(θ) = V exp(−2π² σ_p²/λ²), with
    σ_p/λ = 0.12 π^{1/3} (D/r₀)^{−1/6} θ/θ₀. The seeing ε and isoplanatic
    angle θ₀ are the MASS-DIMM values at 500 nm, scaled as λ^{6/5}, and
    r₀ = 0.98 λ/ε. The loss goes as exp(−a/λ²), the chromatic coherence loss
    of design/gravity_calibration_review.md, §4.

    The paper finds the model follows the trend with separation, but single
    observations scatter about it (its §4.3): use it to judge the size of
    the effect, not to correct for it.
    """
    wave = np.asarray(wave, float)
    scale = (wave / 500e-9) ** 1.2
    r0 = 0.98 * 500e-9 / (seeing_arcsec * ARCSEC) * scale
    theta0 = theta0_arcsec * scale
    sigma_over_lambda = (
        0.12 * np.pi ** (1 / 3) * (diameter_m / r0) ** (-1 / 6) * (separation_arcsec / theta0)
    )
    return np.exp(-4 * np.pi**2 * sigma_over_lambda**2)  # V² = V·V


def anisoplanatic_v2_loss(header, wave):
    """1 − the predicted V² factor for this exposure's fibre offset, or None.

    None when the header lacks the fibre offset, seeing or isoplanatic angle.
    """
    try:
        separation = fibre_separation_mas(header) / 1000.0
    except KeyError:
        return None
    seeing = _mean_card(header, "ESO ISS AMBI FWHM")
    theta0 = _mean_card(header, "ESO ISS AMBI THETA0")
    if seeing is None or theta0 is None:
        return None
    factor = anisoplanatic_v2_factor(wave, separation, seeing, theta0, telescope_diameter(header))
    return 1 - factor


# ---------------------------------------------------------------------------
# Fibre injection (the primary beam)
# ---------------------------------------------------------------------------

# Coupling into a single-mode fibre falls off with a component's offset from
# the fibre centre. Its field of view is about 60 mas at K on the UTs
# (Lacour et al. 2019, §2, via design/gravity_calibration_review.md, §0
# item 8), close to λ/D = 55 mas at 2.2 µm; the review (§1c.2) models the
# coupled power as a Gaussian of FWHM ≈ λ/D. The factor below scales that.
FIBRE_FWHM_LAMBDA_OVER_D = 1.0

# A component farther than this fraction of the fibre FWHM from its centre
# loses more than a quarter of its coupled power (exp(−4 ln 2/9) = 0.73), so
# the flat primary beam that virgil assumes no longer holds (review §1c.2;
# the threshold is the review's estimate).
FIBRE_OFFSET_LIMIT = 1 / 3

FWHM_PER_SIGMA = 2 * np.sqrt(2 * np.log(2))


def fibre_fwhm_mas(wave, diameter_m):
    """FWHM (mas) of the fibre's coupled-power beam on the sky: λ/D.

    ``wave`` in metres; ``diameter_m`` 8.2 for the UTs and 1.8 for the ATs
    (:func:`telescope_diameter`), which makes the AT beam 4.6 times wider.
    """
    return FIBRE_FWHM_LAMBDA_OVER_D * np.asarray(wave, float) / diameter_m / MAS


def fibre_coupling(offset_mas, wave, diameter_m, jitter_mas=0.0, amplitude=False):
    """Fibre coupling at ``offset_mas`` from the fibre centre, relative to on-axis.

    The coupled power η is a Gaussian of FWHM :func:`fibre_fwhm_mas`. With
    tip-tilt jitter of ``jitter_mas`` r.m.s. per axis (the Galactic Centre
    binary fits use about 15 mas on the UTs, GRAVITY Collaboration 2020,
    App. A.3, via the review §1c.2), the Gaussian is averaged over the
    jitter: a Gaussian of width s becomes s²/(s² + j²) exp(−r²/2(s² + j²)).
    The result is relative to a perfectly pointed fibre, so jitter lowers
    even the on-axis value.

    Parameters
    ----------
    offset_mas : float or array
        Distance of the component from the fibre centre (mas).
    wave : float or array
        Wavelengths (m); broadcast against ``offset_mas``.
    diameter_m : float
        Telescope diameter.
    jitter_mas : float
        Tip-tilt jitter, r.m.s. per axis (mas).
    amplitude : bool
        Return the coupled field amplitude √η, averaged over the jitter,
        instead of the power η. A component's flux at telescope k is
        f ⟨η_k⟩; its coherent flux on baseline kl is f ⟨√η_k⟩⟨√η_l⟩ when the
        telescopes jitter independently, so the two differ once there is
        jitter.

    Returns
    -------
    array
        The coupling, 1 on axis without jitter.
    """
    if jitter_mas < 0:
        raise ValueError(f"jitter_mas must be non-negative, not {jitter_mas}")
    s = fibre_fwhm_mas(wave, diameter_m) / FWHM_PER_SIGMA  # σ of the power beam
    if amplitude:
        s = np.sqrt(2.0) * s  # √η is a Gaussian √2 wider
    var = s**2 + jitter_mas**2
    r = np.asarray(offset_mas, float)
    return s**2 / var * np.exp(-(r**2) / (2 * var))


AT_DIAMETER_M = 1.8
UT_DIAMETER_M = 8.2


def fibre_response(sep_mas, wave_m, telescope="AT", jitter_mas=0.0, fwhm_scale=1.0):
    """Fractional flux injected for a source ``sep_mas`` from the fibre centre.

    A Gaussian of FWHM ``fwhm_scale`` · λ/D (D = 1.8 m for ``"AT"``, 8.2 m for
    ``"UT"``) with isotropic tip-tilt jitter of ``jitter_mas`` r.m.s. per axis
    added in quadrature to the Gaussian's width (review §1c.2). The value is
    relative to a perfectly pointed fibre on axis, so it is 1 at zero
    separation without jitter; with jitter it carries the s²/(s² + j²)
    normalisation of :func:`fibre_coupling`, which this wraps. At 2.2 µm on the
    ATs: 0.64 at 100 mas, 0.37 at 150 mas, 0.17 at 200 mas. Q12 to the GRAVITY
    team (the true coupled mode) is open, so ``fwhm_scale`` is a free knob.

    ``sep_mas`` and ``wave_m`` (metres) broadcast against each other; the
    function is plain NumPy, not JAX.
    """
    try:
        diameter = {"AT": AT_DIAMETER_M, "UT": UT_DIAMETER_M}[telescope]
    except KeyError:
        raise ValueError(f"telescope must be 'AT' or 'UT', not {telescope!r}") from None
    if fwhm_scale <= 0:
        raise ValueError(f"fwhm_scale must be positive, not {fwhm_scale}")
    # fibre_coupling's FWHM is FIBRE_FWHM_LAMBDA_OVER_D * lambda/D; rescale its wavelength.
    wave = np.asarray(wave_m, float) * fwhm_scale / FIBRE_FWHM_LAMBDA_OVER_D
    return fibre_coupling(sep_mas, wave, diameter, jitter_mas=jitter_mas)


# ---------------------------------------------------------------------------
# Calibrator spectral slope and K-band masks
# ---------------------------------------------------------------------------

# Anchor effective temperatures (K) by spectral type, at subtype 0 and 5 of each
# letter, for each luminosity class. Dwarfs follow Pecaut & Mamajek 2013 (ApJS
# 208, 9, Table 5); giants and supergiants are rounded values in the spirit of
# the standard calibrations (Cox 2000, Allen's Astrophysical Quantities; de
# Jager & Nieuwenhuijzen 1987). They were entered from memory, not re-read
# against those tables, so they are good to a few per cent (unverified); the
# slope depends on T_eff only weakly (see ``calibrator_sed_slope``).
_LETTERS = "OBAFGKM"
_TEFF_ANCHORS = {
    "V": [41400, 31400, 15700, 9700, 8080, 7220, 6510, 5920, 5660, 5280, 4410, 3850, 3050],
    "III": [39000, 29000, 15000, 10000, 8000, 7000, 6500, 5850, 5050, 4800, 4000, 3850, 3470],
    "I": [40000, 26000, 14000, 9700, 8500, 7500, 6600, 5500, 4800, 4500, 3900, 3700, 3450],
}
# Index of each anchor on the continuous type scale 10*letter + subtype; O5 is
# the first anchor, M5 the last.
_ANCHOR_INDEX = [10 * i + sub for i in range(len(_LETTERS)) for sub in (0, 5)][1:]
_K_BAND_UM = (2.0, 2.4)
_HC_OVER_K = 1.438777e-2  # m K


def _parse_spectral_type(spectral_type):
    # the class may be followed by a suffix ('IIIb', 'III/IV'); only the first is used
    m = re.match(r"\s*([OBAFGKM])\s*(\d(?:\.\d)?)?\s*(IV|I{1,3}|V)?", str(spectral_type))
    if m is None:
        raise ValueError(f"cannot parse spectral type {spectral_type!r}")
    letter, sub, lum = m.groups()
    return letter, float(sub) if sub else 5.0, lum


def calibrator_sed_slope(spectral_type):
    """K-band slope d ln F_λ / d ln λ of a calibrator, and its error.

    Returns ``(slope, err)``. The slope is that of a blackbody at the star's
    effective temperature, averaged over 2.0-2.4 µm:
    d ln B_λ / d ln λ = −5 + x / (1 − e^−x), x = hc/λkT, so it tends to −4
    (Rayleigh-Jeans) for hot stars and rises to about −2.6 for a 3000 K dwarf (a 4800 K
    K0 giant gives −3.2; the −3.5 once hard-coded for K giants is steeper than
    this blackbody value by 0.3, which is within the error returned). T_eff comes
    from a small O-M × V/III/I table (anchors at subtypes 0 and 5, linearly
    interpolated; see the comments at ``_TEFF_ANCHORS`` for the sources). Luminosity
    classes II and IV map to III and V; a missing class is taken as III with a
    larger error.

    ``err`` is a judgement, not a fit: 0.1 for hot stars, and 0.3 below 4600 K
    (where CO and H2O bands and H⁻ opacity make the continuum depart from a
    blackbody) or when the luminosity class is missing. It should be read as a
    prior width on a nuisance slope.
    """
    letter, sub, lum = _parse_spectral_type(spectral_type)
    guessed = lum is None
    lum = {"II": "III", "IV": "V", None: "III"}.get(lum, lum)
    idx = 10 * _LETTERS.index(letter) + sub
    idx = float(np.clip(idx, _ANCHOR_INDEX[0], _ANCHOR_INDEX[-1]))
    teff = float(np.interp(idx, _ANCHOR_INDEX, _TEFF_ANCHORS[lum]))
    wave = np.linspace(*_K_BAND_UM, 41) * 1e-6
    x = _HC_OVER_K / (wave * teff)
    slope = float(np.mean(-5 + x / (-np.expm1(-x))))
    err = 0.3 if (guessed or teff < 4600) else 0.1
    return slope, err


BRGAMMA_UM = 2.1661  # H I 7-4, vacuum wavelength
TELLURIC_UM = (2.06, 2.40)


def brgamma_mask(wave_um, resolving_power=500.0, nfwhm=1.5):
    """True where a calibrator's Brγ absorption affects a channel (to be masked).

    The half-width is ``nfwhm`` resolution FWHMs, nfwhm · λ/R, about 6.5 nm for
    MEDIUM (R ≈ 500) and 0.6 nm for HIGH (R ≈ 4000) with the defaults. Those
    are our choices; widen ``nfwhm`` for broad lines in hot stars.
    """
    half = nfwhm * BRGAMMA_UM / resolving_power
    return np.abs(np.asarray(wave_um, float) - BRGAMMA_UM) <= half


def telluric_mask(wave_um, lo_um=TELLURIC_UM[0], hi_um=TELLURIC_UM[1]):
    """True where a channel is in the CO2/H2O-contaminated edges, below 2.06 or above 2.40 µm.

    The bands are masked rather than modelled (review §1c.3 stays open).
    """
    wave = np.asarray(wave_um, float)
    return (wave < lo_um) | (wave > hi_um)


# ---------------------------------------------------------------------------
# Smearing: how far a scene is from needing it
# ---------------------------------------------------------------------------

# GRAVITY's spectral response per channel, as a Gaussian FWHM in channel
# spacings: PMOIRED's GRAVITY kernel (Mérand 2022, code, via the review §8a).
# MEDIUM samples about 2 channels per resolution element (GRAVITY
# Collaboration 2018b, Methods), so 1.4 is the less conservative choice.
SPECTRAL_KERNEL_PX = 1.4

OMEGA_EARTH = 7.2921150e-5  # sidereal rotation, rad/s


def spectral_fwhm(spacing_m, nbin=1, kernel_px=SPECTRAL_KERNEL_PX):
    """FWHM (m) of the spectral response of one channel after binning by ``nbin``.

    The second moments add: a Gaussian kernel of ``kernel_px`` channel
    spacings, plus the mean of ``nbin`` channels (a discrete uniform
    distribution, variance (nbin² − 1)/12 spacings²). The result is the
    FWHM of a Gaussian of the same variance.
    """
    var = (kernel_px / FWHM_PER_SIGMA) ** 2 + (nbin**2 - 1) / 12
    return FWHM_PER_SIGMA * np.sqrt(var) * np.asarray(spacing_m, float)


def bandwidth_smearing_v2_factor(baseline_m, separation_mas, wave, fwhm_m):
    """Bandwidth smearing: γ², the V² factor of a point ``separation_mas`` from the phase centre.

    The fringe phase 2π B θ/λ of a point θ from the phase centre (along a
    projected baseline B) changes across a channel. Averaging the coherent
    flux over a Gaussian response of ``fwhm_m`` (:func:`spectral_fwhm`)
    multiplies that point's complex visibility by γ = exp(−φ²/2), with
    φ = 2π B θ σ_λ/λ², and so its V² by γ² = exp(−φ²), which this returns.
    In a binary with the phase centre on the primary, the companion's term
    of V² gets γ² and the cross term γ = √(returned value). It is a
    diagnostic: whether 1 − γ² exceeds the V² errors says whether the model
    needs smearing (design/gravity_calibration_review.md, §8; the model
    itself belongs in virgil). Arguments broadcast.
    """
    wave = np.asarray(wave, float)
    sigma = np.asarray(fwhm_m, float) / FWHM_PER_SIGMA
    phi = 2 * np.pi * np.asarray(baseline_m) * np.asarray(separation_mas) * MAS * sigma / wave**2
    return np.exp(-(phi**2))


def time_smearing_v2_factor(baseline_m, separation_mas, wave, int_time_s):
    """Time smearing over a coherent integration of ``int_time_s``: a lower bound on γ².

    During an integration the projected baseline turns with the Earth, by at
    most Δu = B ω⊕ T, so the fringe phase of a point θ from the phase centre
    sweeps at most Φ = 2π B ω⊕ T θ/λ. A uniform sweep multiplies its complex
    visibility by γ = sinc(Φ/2) and its V² by γ² = sinc²(Φ/2), which this
    returns (the cross term of a binary gets γ, as in
    :func:`bandwidth_smearing_v2_factor`): a V² loss of about 4% for 30 mas
    on 130 m over 180 s at 2.2 µm (review §8b). The rate is that of the worst
    geometry, so the true sweep is smaller. Below the first null (Φ < 2π)
    sinc² falls monotonically, so a smaller sweep keeps more coherence and
    sinc²(Φ/2) is a lower bound. Beyond it the side lobes rise again, and
    some smaller true sweep (Φ = 2π) gives γ² = 0, so the bound over every
    sweep up to Φ is 0, which is returned there: the scene then certainly
    needs time smearing modelled. Whether the pipeline averages
    the coherent flux over the whole exposure, or only per DIT, is
    unverified (review §8b): with ``int_time_s`` the DIT the bound is the
    optimistic case. Arguments broadcast.
    """
    phi = (
        2
        * np.pi
        * np.asarray(baseline_m)
        * OMEGA_EARTH
        * np.asarray(int_time_s)
        * np.asarray(separation_mas)
        * MAS
        / np.asarray(wave, float)
    )
    # numpy's sinc(x) = sin(πx)/(πx); past the first null the bound is 0.
    return np.where(phi < 2 * np.pi, np.sinc(phi / (2 * np.pi)) ** 2, 0.0)


# ---------------------------------------------------------------------------
# Wavelength calibration
# ---------------------------------------------------------------------------

# Provisional priors on the wavelength scale, per spectral resolution:
# (σ of the scale s, σ of the offset δ in m), for λ' = λ(1 + s) + δ. The
# pipeline manual quotes 0.5 nm absolute accuracy (≈2.3×10⁻⁴ at 2.2 µm) and
# 0.1 nm between baselines (GRAVITY pipeline manual 1.11.0, §5.1); the
# review (§6c.1, §12.5) adopts N(0, 2×10⁻⁴) in scale and N(0, 0.5 nm) in
# offset until the GRAVITY team answers (design/questions_for_gravity_team.md).
# LOW is left out: its effective wavelengths depend on the source's colour
# (review §6c.2).
WAVELENGTH_PRIORS = {"MEDIUM": (2e-4, 0.5e-9), "HIGH": (2e-4, 0.5e-9)}


def wavelength_prior(header):
    """The provisional wavelength-calibration prior for a product, as 1σ widths.

    Returns ``{"scale": σ_s, "offset_m": σ_δ}`` for λ' = λ(1 + s) + δ, from
    :data:`WAVELENGTH_PRIORS` by ``ESO INS SPEC RES``. In virgil they are the
    noise terms ``wavel_scale`` (``dist.Normal(1, σ_s)``, a factor about 1)
    and ``wavel_offset`` (``dist.Normal(0, σ_δ)``, in metres). Fix the offset
    at zero unless spectral lines constrain it, and note that the scale is
    degenerate with a fitted line velocity (review §6c). The 0.1 nm
    mismatch between baselines is a closure-phase error, not a wavelength
    one (§6c.3). These widths are provisional (review §12.5).

    Raises ``ValueError`` for a resolution without a prior (LOW, or none
    recorded).
    """
    res = str(header.get("ESO INS SPEC RES", "")).strip().upper()
    if res not in WAVELENGTH_PRIORS:
        raise ValueError(
            f"No wavelength prior for spectral resolution {res or None!r}: "
            f"only {sorted(WAVELENGTH_PRIORS)} (LOW depends on the source's colour)"
        )
    scale, offset = WAVELENGTH_PRIORS[res]
    return {"scale": scale, "offset_m": offset}
