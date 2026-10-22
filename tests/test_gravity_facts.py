"""GRAVITY instrument facts: fibre injection, smearing, wavelength priors."""

import numpy as np
import pytest
from astropy.io import fits

from virgil_vlti.gravity import (
    FIBRE_OFFSET_LIMIT,
    MAS,
    OMEGA_EARTH,
    bandwidth_smearing_v2_factor,
    brgamma_mask,
    calibrator_sed_slope,
    fibre_coupling,
    fibre_fwhm_mas,
    fibre_response,
    spectral_fwhm,
    telluric_mask,
    time_smearing_v2_factor,
    wavelength_prior,
)


def test_fibre_beam_matches_the_quoted_field():
    # Lacour et al. 2019 §2: about 60 mas at K on the UTs; the ATs are wider
    # by the aperture ratio.
    ut = fibre_fwhm_mas(2.2e-6, 8.2)
    assert 50 < ut < 65
    np.testing.assert_allclose(fibre_fwhm_mas(2.2e-6, 1.8) / ut, 8.2 / 1.8)
    np.testing.assert_allclose(fibre_fwhm_mas(2.4e-6, 8.2) / fibre_fwhm_mas(2.0e-6, 8.2), 1.2)


def test_fibre_coupling_without_jitter():
    fwhm = fibre_fwhm_mas(2.2e-6, 8.2)
    np.testing.assert_allclose(fibre_coupling(0.0, 2.2e-6, 8.2), 1.0)
    np.testing.assert_allclose(fibre_coupling(fwhm / 2, 2.2e-6, 8.2), 0.5)
    np.testing.assert_allclose(fibre_coupling(fwhm / 2, 2.2e-6, 8.2, amplitude=True), np.sqrt(0.5))
    # At the offset limit a component keeps about three quarters of its power.
    kept = fibre_coupling(FIBRE_OFFSET_LIMIT * fwhm, 2.2e-6, 8.2)
    np.testing.assert_allclose(kept, np.exp(-4 * np.log(2) / 9))
    with pytest.raises(ValueError, match="jitter"):
        fibre_coupling(0.0, 2.2e-6, 8.2, jitter_mas=-1.0)


@pytest.mark.parametrize("amplitude", [False, True])
def test_fibre_coupling_with_jitter_is_the_jitter_average(amplitude):
    # Monte Carlo over independent Gaussian pointing errors per axis, against
    # the closed form.
    rng = np.random.default_rng(1)
    jitter, wave, diameter = 15.0, 2.2e-6, 8.2
    pointing = rng.normal(0.0, jitter, (400_000, 2))
    for offset in (0.0, 20.0, 45.0):
        r = np.hypot(offset + pointing[:, 0], pointing[:, 1])
        eta = fibre_coupling(r, wave, diameter)
        sample = np.sqrt(eta) if amplitude else eta
        closed = fibre_coupling(offset, wave, diameter, jitter, amplitude=amplitude)
        np.testing.assert_allclose(closed, sample.mean(), rtol=3e-3)
    assert fibre_coupling(0.0, wave, diameter, jitter) < 1.0


def test_spectral_fwhm_adds_second_moments():
    spacing = 2e-9
    np.testing.assert_allclose(spectral_fwhm(spacing), 1.4 * spacing)
    # Numerically: a Gaussian kernel at each of nbin channel centres, averaged.
    x = np.linspace(-20, 20, 40001)
    sigma = 1.4 / (2 * np.sqrt(2 * np.log(2)))
    for nbin in (2, 4):
        centres = np.arange(nbin) - (nbin - 1) / 2
        profile = sum(np.exp(-((x - c) ** 2) / (2 * sigma**2)) for c in centres)
        var = np.sum(profile * x**2) / np.sum(profile)
        expected = 2 * np.sqrt(2 * np.log(2)) * np.sqrt(var) * spacing
        np.testing.assert_allclose(spectral_fwhm(spacing, nbin), expected, rtol=1e-6)


def test_bandwidth_smearing_against_integration():
    wave, fwhm = 2.2e-6, 4e-9
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    lam = wave + np.linspace(-6, 6, 20001) * sigma
    weight = np.exp(-((lam - wave) ** 2) / (2 * sigma**2))
    for baseline, sep in ((130.0, 30.0), (100.0, 200.0)):
        phase = 2 * np.pi * baseline * sep * MAS / lam
        coherence = np.abs(np.sum(weight * np.exp(1j * phase)) / np.sum(weight))
        factor = bandwidth_smearing_v2_factor(baseline, sep, wave, fwhm)
        np.testing.assert_allclose(factor, coherence**2, rtol=2e-3)
    assert bandwidth_smearing_v2_factor(130.0, 0.0, wave, fwhm) == 1.0


def test_time_smearing_against_integration_and_the_review():
    wave, baseline, sep, t = 2.2e-6, 130.0, 30.0, 180.0
    # The phase sweeps uniformly at the largest rate, B ω T θ/λ in fringes.
    phase = 2 * np.pi * baseline * OMEGA_EARTH * np.linspace(0, t, 100001) * sep * MAS / wave
    coherence = np.abs(np.mean(np.exp(1j * phase)))
    factor = time_smearing_v2_factor(baseline, sep, wave, t)
    np.testing.assert_allclose(factor, coherence**2, rtol=1e-6)
    # design/gravity_calibration_review.md §8b: about 4%.
    assert 0.035 < 1 - factor < 0.045
    assert time_smearing_v2_factor(baseline, sep, wave, 0.0) == 1.0


def test_time_smearing_bound_is_zero_past_the_first_null():
    wave, baseline, sep = 2.2e-6, 130.0, 30.0
    per_second = 2 * np.pi * baseline * OMEGA_EARTH * sep * MAS / wave  # Φ per second
    # Φ_max = 3π: a slower true sweep with Φ = 2π gives γ² = 0, so the bound
    # is 0, not the side lobe sinc²(1.5π) ≈ 0.045.
    assert time_smearing_v2_factor(baseline, sep, wave, 3 * np.pi / per_second) == 0.0
    # A bound over all sweeps up to Φ never rises as Φ grows.
    factors = time_smearing_v2_factor(
        baseline, sep, wave, np.linspace(0, 6 * np.pi, 400) / per_second
    )
    assert np.all(np.diff(factors) <= 0)


def test_wavelength_prior_by_resolution():
    header = fits.Header({"HIERARCH ESO INS SPEC RES": "MEDIUM"})
    assert wavelength_prior(header) == {"scale": 2e-4, "offset_m": 0.5e-9}
    header["HIERARCH ESO INS SPEC RES"] = "LOW"
    with pytest.raises(ValueError, match="LOW"):
        wavelength_prior(header)
    with pytest.raises(ValueError, match="None"):
        wavelength_prior(fits.Header())


def test_fibre_response_lambda_over_d_limit_and_at_values():
    wave = 2.2e-6
    for tel, d in (("AT", 1.8), ("UT", 8.2)):
        fwhm = wave / d / MAS
        np.testing.assert_allclose(fibre_response(0.0, wave, tel), 1.0)
        np.testing.assert_allclose(fibre_response(fwhm / 2, wave, tel), 0.5)
        np.testing.assert_allclose(fibre_response(fwhm, wave, tel, fwhm_scale=2.0), 0.5)
    np.testing.assert_allclose(
        fibre_response(np.array([100.0, 150.0, 200.0]), wave, "AT"),
        [0.64, 0.37, 0.17],
        atol=0.01,
    )
    with pytest.raises(ValueError):
        fibre_response(1.0, wave, "VLT")


def test_fibre_response_jitter_limit_and_monotonic():
    sep = np.linspace(0, 400, 81)
    np.testing.assert_allclose(
        fibre_response(sep, 2.2e-6, "AT", jitter_mas=1e-9), fibre_response(sep, 2.2e-6, "AT")
    )
    for jitter in (0.0, 15.0):
        assert np.all(np.diff(fibre_response(sep, 2.2e-6, "AT", jitter_mas=jitter)) < 0)
    # jitter widens the beam in quadrature: a UT fibre at 60 mas loses more of the on-axis peak,
    # but gains in the wings
    assert fibre_response(200.0, 2.2e-6, "UT", jitter_mas=15.0) > fibre_response(
        200.0, 2.2e-6, "UT"
    )
    # chromatic: the beam is wider at longer wavelengths
    assert fibre_response(100.0, 2.4e-6, "AT") > fibre_response(100.0, 2.0e-6, "AT")


def test_calibrator_sed_slope_table():
    for lum in ("V", "III", "I"):
        types = [f"{c}{n}{lum}" for c in "OBAFGKM" for n in (0, 5)][1:]
        slopes = [calibrator_sed_slope(t)[0] for t in types]
        assert np.all(np.diff(slopes) >= -1e-12), lum  # hotter is steeper, toward -4
        assert -4.0 < slopes[0] < slopes[-1] < -2.0
    slope, err = calibrator_sed_slope("K0III")
    assert -3.3 < slope < -3.0 and err > 0
    assert calibrator_sed_slope("F5III")[0] < calibrator_sed_slope("K0III")[0]
    assert calibrator_sed_slope("B2V")[0] < -3.8
    assert calibrator_sed_slope("K1IIIb") == calibrator_sed_slope("K1III")
    assert calibrator_sed_slope("K1")[1] >= 0.3  # no class: wider error
    with pytest.raises(ValueError):
        calibrator_sed_slope("DA2")


def test_brgamma_and_telluric_masks():
    wave = np.array([2.00, 2.10, 2.1661, 2.1700, 2.30, 2.45])
    assert brgamma_mask(wave, 500).tolist() == [False, False, True, True, False, False]
    assert brgamma_mask(wave, 4000).tolist() == [False, False, True, False, False, False]
    assert telluric_mask(wave).tolist() == [True, False, False, False, False, True]
