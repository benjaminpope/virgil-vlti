# MATISSE calibration: a literature review of products, error models and pitfalls

Status: **review**, 2026-10-04. Nothing here is implemented. This note answers the questions in [`matisse.md`](matisse.md), in the style of [`gravity_calibration_review.md`](gravity_calibration_review.md) (cited below as "the GRAVITY review"). It compares what the MATISSE pipeline and the MATISSE consortium do with what virgil reads and fits, and with what `virgil_vlti.calibrate` does. It respects the boundary in [`PLAN.md`](../PLAN.md): generic machinery is recommended for virgil core, instrument facts and calibration for virgil-vlti. It came with pointers to it in `PLAN.md` and `matisse.md`, and otherwise changes no other note.

For each design element: **(a)** what the literature and the pipeline do, **(b)** whether our design agrees, is incomplete, or is wrong or risky, and **(c)** the change I recommend. Citation keys are in §11, with how far each source was read. **[unverified]** marks a claim I could not check against a primary source. "Our files" are the three calibrated ν Hor products (local only, not in this repository) (L band, LOW, chopped); I read their headers and table layouts, plus a few single-number summaries computed from their arrays (median VISAMP/√V², median T3AMP, median relative errors; §1, §11).

## 0. Footguns to avoid, in priority order

Each item is something the MATISSE literature, the pipeline manual, the consortium tools or our own files show to be wrong or misleading. Items 1–5 affect virgil **today**, if a user hands it a MATISSE file.

1. **Correlated flux stored where a visibility is expected.** With `corrFlux=TRUE` (the usual choice in N), the pipeline writes the *squared correlated flux* into `OI_VIS2.VIS2DATA`, and the correlated flux into `OI_VIS.VISAMP` ([PM] §7, product table; §9.1, §9.11). OIFITS gives `OI_VIS2` no keyword to say so. `read_oifits` reads `VIS2DATA` as V², and falls back to `VISAMP` without looking at `AMPTYP` (`src/virgil/oifits.py`, `_read_visibilities`). The units are then ADU² or Jy², not V². Even `AMPTYP` is unreliable: the manual says the default `VISAMP` is a *differential* visibility with `AMPTYPE = 'differential'` ([PM] §7), but our files say `AMPTYP = 'absolute'`, and their `VISAMP` is close to √V² (median ratio 0.98–1.00). **Fix:** virgil-vlti detects `corrFlux` from the provenance keywords (`ESO PRO REC1 PARAMn NAME = corrFlux`) and refuses, or relabels, before anything reaches virgil. virgil core checks `AMPTYP` on `OI_VIS` and refuses anything but `absolute` as a visibility amplitude (§1, §9).
2. **Trusting pipeline error bars.** They are statistical: the dispersion of sub-sets of frames within one exposure ("chunks" [Va21 App. B]; a bootstrap [PM §2.6]). They contain no calibration term. The dominant errors are broadband: the transfer function moves with the coherence time (by a factor of 2 in L between τ₀ > 5 ms and τ₀ < 3 ms [Lo22 §4.1.2]), giving visibility errors of 0.015–0.08 common to all channels ([Lo22] Table 5). In L, chop residuals add a flux-dependent photometric error ([Lo22] §4.1.4, Table 6). The consortium replaces the errors: Varga et al. add the scatter between exposures and set floors of 0.03 in L visibility, 1° in L closure phase and 8% in N correlated flux ([Va21] §3.1, App. B); Isbell et al. use the scatter of exposures within a cycle and between BCD positions ([Is22] §2.2). Worse, the consortium merge tool *divides* errors by √N when it averages exposures and BCD positions (inverse-variance weights for V², VISAMP and flux; σ/√N for T3PHI and VISPHI; the scatter terms are commented out) ([MT] `libPostTools.mat_mergeOifits`). Our files appear to be such merges (§1).
3. **N-band phases rotated by 180°.** From pipeline 1.5.2 until 2.0.0, N-band phases and closure phases had the opposite sign convention to L and M, which rotates N-band images and models by 180° ([PM] §2.5; [GR22] §4.3, §4.8, who found and corrected it on a known binary). Archival N-band products may carry either convention. **Fix:** virgil-vlti reads `ESO PRO RECi PIPE ID` and, for N-band products from 1.5.2 ≤ version < 2.0.0, negates `T3PHI` and `VISPHI` or refuses them **[the exact first affected version is as the two sources state it; I have not tested a product]**.
4. **Channels that are not independent.** In LOW the spectra are oversampled: our files hold 64 channels from 2.82 to 4.23 µm, 0.022 µm apart, while R ≈ 30 gives a resolution element of about 0.12 µm, i.e. about 5 channels. On top of that the pipeline applies a 5-pixel (L) or 7-pixel (N) spectral average by default in the consortium scripts and recommends it ([PM] §2.14; [MT] `mat_autoPipeline.py`); since 2.5.0 it is a *sliding* average ([PM] §2.1). Our files used `spectralBinning = 5`. Treating 64 channels as independent counts each piece of information about 5 times. `EFF_BAND` in our files is the channel spacing (0.0219 µm), not the resolution width, so a smearing kernel built from `EFF_BAND` is about 5 times too narrow (§5, §7).
5. **Calibration errors common to all channels, and to baselines sharing a telescope.** The broadband transfer-function error is shared by all channels of an exposure and baseline ([Lo22] §4.1.2 calls these "broadband calibration errors"). The broadband photometric error ε_p of telescope i enters the V² of its three baselines through ε_p/F_i ([Lo22] App. A.2, eq. A.9). Diagonal errors, even inflated, miss both. This is the same structure as the GRAVITY review's §4, and the same low-rank modes fit it (§5).
6. **BCD positions mixed in one transfer function.** The transfer function depends on the BCD configuration, so calibration must be done per BCD position, before the positions are combined ([Lo22] §3.3; [GR22] §4.5). The pipeline's BCD combination (`cumulBlock=TRUE` in `mat_cal_oifits`) and the consortium merge both average positions after calibration ([PM] §11.9; [MT]). `virgil_vlti.calibrate` would mix positions, since it matches on `STA_INDEX` only (§8).
7. **Misleading headers in merged products.** The consortium merge copies the primary header of its *first* input file ([MT] `mat_mergeOifits`: the output is a copy of `data[0]`). Our files therefore say `ESO CFG BCD MODE = 'OUT-OUT'` and `ESO DRS BCD COMBINED = F`, although each holds six V² rows built from several exposures (`INT_TIME` = 120.5 s; the header lists 12 raw target files, 8 of them chopped by their times and the standard sequence). Never infer BCD state, exposure time or calibrators from the primary header of a merged file without checking the history **[my reading of the merge code against our files; the exact set of merged exposures is not recorded in the file]**.
8. **Differential phases of unknown order.** `OI_VIS` is `PHITYP = 'differential'` in our files, with no `PHIORDER` keyword. The DRS design removes "any wavelength-linear phase term" ([Mi16] §2.5.2). The chromatic-OPD model changed in 2.5.0, from Mathar's refractive index to Voronin & Zheltikov's ([PM] §2.1). Only two of the six differential phases can be corrected by the BCDs; the other four are only averaged ([PM] §11.9). virgil already refuses differential `VISPHI` as an absolute phase. Fitting it needs the normalisation operator of the GRAVITY review §2 (§2 here).
9. **T3AMP is not data.** Before 2.5.0 the pipeline wrote zeros (our files: median T3AMP = 0) and the manual calls the column "not used" ([PM] §7). From 2.5.0, T3AMP is *derived from VIS2* ([PM] §2.1). Fitting it alongside V² counts V² twice.
10. **Smearing, in wavelength and in time.** At LOW resolution a companion 30 mas away loses about 5% of its fringe contrast on a 130 m baseline in L (§7) **[my estimate]**. Merged products average u, v over about 10 minutes: on a 130 m baseline that sweeps a component 30 mas from the phase centre through up to 1.5 rad of fringe phase in L (§7) **[my estimate]**. virgil models neither.
11. **Total flux through a pinhole.** `OI_FLUX` is uncalibrated (`CALSTAT = 'U'`, ADU in our files) and is the flux through a spatial filter of 1.5λ/D (L) or 2λ/D (N) ([Lo22] Table 2 notes). For an extended source it depends on that filter and on the extraction software: Isbell et al. found an N-band flux three times larger from the DRS than from EWS for the Circinus nucleus, and attribute it to the spatial filter each uses ([Is22] §2.2). The default `replaceTel = 3` replaces the photometry of AT3 by the mean of the other three ([PM] §9.1.4), and our files used it. In N the photometry is taken *after* the fringes, not simultaneously ([Lo22] §3.1).
12. **Coherent-flux bias near the sensitivity limit.** Errors in the per-frame OPD estimate bias the coherent flux low for faint N-band sources, below about 5 Jy at 8.5 µm and 8 Jy at 11 µm on the ATs ([Lo22] §4.2, Fig. 18). Varga et al. dropped such data because the bias could not be given an error ([Va21] App. B). The error bars do not show it.
13. **`virgil_vlti.classify` mislabels MATISSE products.** MATISSE product categories are `TARGET_RAW_INT`, `CALIB_RAW_INT` and `TARGET_CAL_INT` ([PM] §7.5, §9.1.2). `classify` tests for the substrings `_CAL_` and `_SCI_`: `TARGET_CAL_INT` (calibrated *science*) contains `_CAL_` and is labelled a calibrator; the two raw categories match neither and raise (§8).

## 1. Products and the OIFITS layout

**(a) Literature and inspection.**
- **Observables.** MATISSE measures coherent (correlated) flux, visibility, closure phase and differential phase, and derives differential visibility ([PM] §4.1). The ESO products store: `OI_VIS2` (squared visibility, or squared correlated flux with `corrFlux`), `OI_T3` (closure phase; T3AMP "not used"), `OI_VIS` (`VISAMP` = differential visibility by default, or coherent correlated flux with `corrFlux` and `coherentAlgo`; `VISPHI` = differential phase), `OI_FLUX` (total spectrum, ADU) and, in calibrator products, a non-standard `TF2` table holding the transfer function ([PM] §7, product table).
- **Categories.** `mat_raw_estimates` writes `TARGET_RAW_INT` or `CALIB_RAW_INT`, one per exposure. `mat_cal_oifits` writes `TARGET_CAL_INT` ([PM] §7.5, §9.1.2, §9.5).
- **Bands.** L&M and N are reduced independently, into separate files ([Lo22] §3.3). The consortium tools tell them apart by `ESO DET NAME` (`MATISSE-LM` or `MATISSE-N`) ([MT] `libPostTools.mat_mergeByTplStart`). Our L files carry `INSNAME = 'MATISSE'`. I did not see an N-band product, so I do not know whether its `INSNAME` differs **[unverified]**. If both bands share `INSNAME = 'MATISSE'`, virgil's `insname=` selector cannot separate them when files are merged.
- **Chopped and non-chopped.** In the standard hybrid sequence the L band gets 4 non-chopped exposures (one per BCD position) and 8 chopped ones (4 at IN-IN, 4 at OUT-OUT); the N band gets 4 non-chopped fringe exposures and 8 chopped photometric exposures, one shutter open at a time ([PM] §4.2, Table 4.1). Lopez et al. describe the chopped L block as "2 × 4" exposures "for each of the BCD configurations" ([Lo22] §3.1); I follow the manual's explicit table **[the two descriptions differ]**. The pipeline uses only chopped frames for L&M visibilities and fluxes when it combines BCD positions ([PM] §11.9). The consortium merge keeps chopped and non-chopped L data in separate files on request (`--separateChopping`; the `_Chop` and `_noChop` file names) ([MT] `mat_mergeAllOiFits.py`). They are not separate tables or `INSNAME`s within one file.
- **BCD.** Two beam-commuting devices swap input beams IP1/IP3 and IP5/IP7, giving four configurations ([PM] §4.1; [Lo22] §2.5). Each exposure is reduced separately ([PM] §11.9). The consortium tool reorders rows, negates u, v and phases, and flips `STA_INDEX` to bring IN configurations to the OUT-OUT layout before merging ([MT] `mat_removeBCD`). So per-exposure products of different BCD positions list some baselines and triangles in different orientations **[my inference from that code; not checked on raw products]**.
- **Our files.** Each holds `OI_TARGET`, `OI_ARRAY` (4 ATs, K0 G2 D0 J3), `OI_WAVELENGTH` (64 channels, 4.23 → 2.82 µm, decreasing), `OI_VIS2` (6 rows), `OI_T3` (4 rows), `OI_VIS` (6 rows, `AMPTYP = 'absolute'`, `PHITYP = 'differential'`, no `PHIORDER`) and `OI_FLUX` (4 rows, `CALSTAT = 'U'`, `FLUXDATA` in ADU, `FOV = 1.0` arcsec). All are `INSNAME = 'MATISSE'`, OIFITS2 (`OI_REVN = 2`, but 1 on `OI_FLUX`). `OI_VIS2` rows are (D0,J3), (K0,G2), (G2,D0), (G2,J3), (K0,D0), (K0,J3); every triangle's (a,b), (b,c), (a,c) appears in `OI_VIS2` in that orientation, as `read_oifits` assumes. All rows share one MJD and `INT_TIME` = 120.46 s. Provenance: `mat_raw_estimates` from pipeline 2.0.0 (`corrFlux=false`, `coherentAlgo=2`, `spectralBinning=5`, `replaceTel=3`, `hampelFilterKernel=10`), then `mat_cal_oifits` from 1.7.9 (`tfInterp=0`, `cumulBlock=true`), with one or two calibrators per file listed in `ESO PRO JSDC CALn`. The GRAVITY fringe tracker was on (`ESO DEL FT SENSOR = 'GRAVITY'`: GRA4MAT), with an L-band DIT of 75 ms. The file names follow the consortium scripts (`mat_autoCalib` appends `_cal_oifits_<n>`; the merge keeps `_Chop`).

**(b) Verdict.** virgil's reader **agrees** on the parts it reads: `OI_VIS2` plus `OI_T3`, the triangle orientation, and the refusal of differential `VISPHI`. It is **wrong** for correlated-flux products (footgun 1), and **incomplete**: it ignores `AMPTYP`, and nothing stops a user mixing L and N files, or chopped and non-chopped L files, in one `OIData`.

**(c) Changes.**
1. virgil core: refuse `OI_VIS.VISAMP` as a visibility unless `AMPTYP = 'absolute'`; refuse `'correlated flux'` and `'differential'` with a message naming them. This is generic OIFITS2 ([D17] §5), so it belongs in core.
2. virgil-vlti: a MATISSE product reader that records band (`ESO DET NAME`), resolution (`ESO INS DIL/DIN NAME`), BCD (`ESO INS BCD1/2 NAME`), chopping (`ESO ISS CHOP ST`), `corrFlux`, `spectralBinning`, `replaceTel`, the pipeline versions and whether the file is a merge. It refuses `corrFlux` products as V² (footgun 1).
3. Give each band its own dataset in virgil, even when the `INSNAME`s agree: different noise, smearing and calibration.
4. Use L-band V² from chopped exposures only, as the pipeline does ([PM] §11.9); use closure and differential phases from all exposures, since chopping is "unnecessary" for them ([PM] §4.1).

## 2. Incoherent and coherent processing, and their error models

**(a) Literature.**
- **Incoherent** (default, `corrFlux=FALSE`): V² is the frame-averaged power spectrum at each fringe peak, minus a bias measured between the peaks, divided by the product of the photometries ([Lo22] eqs. 3–4; [PM] §11.5; [Mi16] §2.5.1). The closure phase is the argument of the frame-averaged bispectrum, corrected for additive and multiplicative biases ([PM] §11.5). Nothing needs the OPD.
- **Coherent:** the residual atmospheric phase, modelled as φ₀(t) + 2πδ(t)/λ, is estimated per frame and removed before averaging ([Lo22] eqs. 6–9). `coherentAlgo=1` is the AMBER-like differential estimator; `coherentAlgo=2` is the MIDI-like coherent integration, giving phase and coherent flux, or visibility with `corrFlux=FALSE` ([PM] §11.6). Since 2.5.0 the parameter is `cpAlgo`, and its 2 means closure phases computed coherently rather than from the bispectrum ([PM] §2.1) **[the renaming changes what "2" means; I did not check which estimator our files' closure phases used]**.
- **Noise.** The coherent-flux variance per frame is N_T(n* + n_B) + N_pix σ²_RON ([Lo22] eq. A.1). Coherent integration gains √(N_λ N_F) in S/N; incoherent integration loses fast once S/N per frame < 1 ([Lo22] eqs. A.3–A.4). The phase error is 1/(√2 S/N), and the closure-phase error ≈ √3 times the baseline phase error ([Lo22] eqs. A.5–A.6). The V² error adds the coherent-flux and photometric terms in quadrature ([Lo22] eq. A.8), and the broadband photometric error ε_p/F_i on top ([Lo22] eq. A.9).
- **Practice.** In L, both estimators perform alike, and the incoherent V² is preferred because photometry is simultaneous ([Va21] §3.1; [GR22] §4.3). In N the consortium uses coherent correlated fluxes: photometry is not simultaneous, the coherent estimator has a "significant gain" in S/N, and many targets are below the N-band photometric limit ([Va21] §3.1; [Is22] §2.2). The coherent flux is biased low near the sensitivity limit (footgun 12).
- **Differential phase.** The coherent estimator yields "Φdiff, with any wavelength-linear phase term set to zero" ([Mi16] §2.5.2). The OPD term is 2πδ/λ, so "linear" means linear in 1/λ, as in OIFITS2 ([D17] §3.2) **[my reading]**. The chromatic water-vapour OPD is modelled from ambient data and delay-line positions ([Mi16] §2.5.2; [PM] §11.6). Lopez et al. removed residual chromatic OPD by a low-order polynomial in their performance analysis ([Lo22] §4.2).

**(b) Verdict.** For the products a virgil user will meet first (L-band V² and closure phases), virgil's model of V² and closure phases **agrees** with the incoherent estimator. Correlated flux and differential phase are **not supported**, and must not be read as V² or absolute phase (footguns 1, 8).

**(c) Changes.**
1. Correlated flux as an observable (virgil core, later): the model is F_tot(λ) |V(u,v,λ)|, which needs the total spectrum of the GRAVITY review §1. It is the natural N-band observable.
2. Differential phase: reuse the GRAVITY review's normalisation operator N (a fit in a basis polynomial in 1/λ, degree ≥ 1, applied to data and model alike). Its default degree for MATISSE LOW should be 1 (constant plus 1/λ), raised when chromatic OPD residuals show **[my recommendation; no MATISSE paper states a degree]**.
3. Document that coherent N-band fluxes below the bias limit are not trustworthy, whatever their error bars.

## 3. Calibration: transfer function, photometry, chopping, flux calibration and calibrators

**(a) Literature.**
- **Transfer function, V².** The transfer function is the calibrator's raw V² divided by its model V², with diameters from JSDC ([PM] §9.21; [Va21] eq. 1). `mat_cal_vis` interpolates it to the science time and divides ([PM] §9.22). The recipe option `tfInterp` offers only an *average* (0, the default) or a linear function (2) ([PM] §9.5, §9.22). The design paper planned a Gaussian-weighted average with FWHM about 1 h, as for AMBER ([Mi16] §2.6). The consortium script groups calibrators by BCD, detector, DIT and chopping status, and by default takes every matching calibrator within ±3 h (mode `MEAN`), or the nearest, or the bracketing two; it disables the linear interpolation because it "give[s] weird results" ([MT] `mat_autoCalib.py`). Our files used the average of two calibrators.
- **BCD.** The transfer function depends on the BCD position and must be calibrated first ([Lo22] §3.3). After calibration, V² and fluxes are averaged over positions; closure phases and two differential phases are combined with signs, e.g. (φ_OUT-OUT + φ_OUT-IN − φ_IN-OUT − φ_IN-IN)/4 for the first triangle ([PM] §11.9). For closure phases the BCD combination alone is enough; adding a calibrator gave "no additional improvement" ([Lo22] §4.1.3, Fig. 14). Isbell et al. do both: subtract the calibrator's closure phase, then average the BCD positions ([Is22] §2.2).
- **Photometry and chopping.** In L, photometry is simultaneous (SiPhot) because Strehl and tip-tilt make the flux fluctuate; in N it is sequential (HighSens) because background noise dominates and the photometric channels cannot be calibrated on the N detector ([PM] §4.2; [Lo22] §3.1). Chopping at 0.5–1 Hz removes the background from the photometry only; spatial and temporal modulation remove it from the fringes ([Lo22] §2.4; [Mi16] §2.4). Chop residuals are "the most important contribution to the photometric error" ([PM] §4.1). The broadband photometric error is 0.11 / 0.19 / 2.3 Jy (L/M/N) on the ATs and 0.008 / 0.016 / 0.08 Jy on the UTs ([Lo22] Table 6).
- **Correlated-flux calibration** (N band; not in the DRS): F_corr = F_corr,raw / T, with T = F_corr,raw^cal / (F_tot^cal V^cal), the calibrator's model spectrum times its uniform-disk visibility ([Va21] eqs. 2–3; [Is22] eq. 1, with V^cal = 1). The consortium's `mat_fluxCal` does this, optionally with an airmass correction from an ESO SkyCalc transmission model; its error propagation states that the calibrator diameter "is not taken into account", nor is the model spectrum's error ([MT] `libFluxCal_STARSFLUX_new.py`).
- **Total flux.** Cross-calibrate the photometric spectra with spectrophotometric standards (Cohen et al.), and take errors from the scatter between spectra ([GR22] §4.6).
- **Calibrators.** Similar airmass and direction (Varga et al.: azimuth difference < 30°, airmass difference < 0.2; Gámez Rosas et al.: within 10°), and, for faint targets, within a factor 2 in flux ([Va21] §2; [GR22] §3). Diameters below about 3 mas for L and 9 mas for N on medium and long AT baselines; most such stars are too faint in N, so L and N often need different calibrators ([Va21] §2). Catalogues: MDFC (Cruzalèbes et al. 2019) and the MIDI calibrator list (van Boekel 2004) ([Va21] §2; [Lo22] §4.1.1 footnote). An N-band calibrator also serves as the spectral calibrator, so it needs a reliable model spectrum ([Is22] §2.1).
- **A calibration term in the errors.** Neither the manual nor the papers say that `mat_cal_vis` adds the spread of the transfer function to the errors, and Varga et al. treat the pipeline errors as statistical ([Va21] App. B). The consortium therefore adds its own: σ_sys,s, the scatter of the 4 or 8 exposures (minutes), and σ_sys,l, the scatter between calibrations with different calibrators (hours). Only σ_sys,s went into their error bars; the floors cover the rest ([Va21] App. B). **[Whether `mat_cal_vis` propagates any transfer-function error is unverified; I did not read the pipeline source.]**

**(b) Verdict.** The transfer-function idea of `virgil_vlti.calibrate` **agrees**, and its leave-one-out prediction error is closer to what the consortium does by hand (σ_sys,l) than anything the pipeline offers. It is **incomplete** for MATISSE: no BCD grouping, no chopping selection, no correlated-flux calibration, no calibrator spectra (§8).

**(c) Changes.** See §8 for `calibrate.py`, and §10 for the list. The key one: calibrate per BCD position and chopping status, then combine.

## 4. The N band: thermal background and transfer-function variability

**(a) Literature.**
- The background exceeds the source "by far in most cases" ([Lo22] §2.4, Fig. 3). It forces DITs of 20 ms (LOW) or 75 ms (HIGH) to avoid saturation ([Lo22] §3.1); 25 ms in the manual's example ([PM] §11.8).
- N-band visibilities need photometry taken minutes later, so they carry the change in transfer function between the two ([Va21] §3.1). That is one reason the consortium prefers correlated flux in N.
- Transfer-function errors in N (DIT × modulation = 240 ms): visibility 0.045 / 0.02 / 0.015 and closure phase 1.75° / 0.49° / 0.29° for τ₀ = 3.2 / 6.8 / 7.5 ms ([Lo22] Table 5). Bad coherence time hurts N closure phases much more than L closure phases (0.30° / 0.26° / 0.16°).
- Varga et al. find N-band correlated-flux errors < 8% in good conditions and about 20% in poor ones, dominated by systematics ([Va21] App. B).
- Water vapour makes the N-band OPD fluctuate by microns even with fringe tracking in K; GRA4MAT data allow "a posteriori co-phasing" from the K-band delays ([PM] §11.8). GRA4MAT also stabilises the transfer function by reducing piston jitter ([Lo22] §6).
- The FLAG column marks 9.18–9.56 µm (the ozone band) as bad in N ([PM] §2.12).

**(b) Verdict.** No virgil design element addresses the N band yet. The nuisance structure is the GRAVITY review's §4 one, with larger widths: an achromatic per-(exposure, baseline) gain on correlated flux of several per cent, plus a per-calibrator spectral-template error common to everything it calibrates.

**(c) Changes.**
1. In virgil-vlti, default nuisance widths per band and condition, from [Lo22] Table 5 (as relative errors where possible) and the [Va21] floors, in the versioned table of PLAN item 8. Make τ₀ an input: the start-of-exposure value is in the header (`ESO ISS AMBI TAU0 START`; our three files: 7.3, 8.6 and 3.6 ms).
2. A calibrator-spectrum nuisance for correlated flux: one multiplicative, smooth-in-λ factor per calibrator and night, shared by every target frame it calibrates. It is the N-band analogue of PLAN item 1 (diameter errors).

## 5. Correlations between channels, BCD positions and baselines

**(a) Literature.**
- **Channels.** The detector samples the spectral direction with 3 pixels per λ/D ([Lo22] §2.3), and the pipeline then averages 5 (L) or 7 (N) pixels ([PM] §2.14; [GR22] §4.3.1), as a sliding average since 2.5.0 ([PM] §2.1). Gámez Rosas et al. call both detectors "oversampled" in LOW ([GR22] §4.3.1). Our files show about 5 channels per resolution element (footgun 4). Lopez et al. estimate the fundamental noise as the scatter about a low-order polynomial in λ ([Lo22] §4.1.1), which assumes channel-independent noise and cannot see a broadband component.
- **Broadband (all channels).** Transfer-function errors and photometric errors are broadband by construction ([Lo22] §4.1.2, §4.1.4). Varga et al. call the transfer-function error "systematic in nature", shifting "the overall visibility or flux level of the spectrum" ([Va21] App. B).
- **Baselines.** The photometric error of telescope i enters the V² of its three baselines ([Lo22] eq. A.9), and so do injection and Strehl changes. With `replaceTel` one telescope's photometry is the mean of the others ([PM] §9.1.4), which couples all four.
- **BCD positions.** Varga et al. find the exposure-to-exposure scatter of closure phases is the largest closure-phase error, and that it "reflects instrumental effects (revealed by the beam commutation)"; averaging exposures reduces it ([Va21] App. B). The V² transfer function differs between positions ([Lo22] §3.3).
- No MATISSE paper I read measured inter-channel or inter-baseline covariances from frames, as Kammerer et al. did for GRAVITY [K20]. OIFITS2 `OI_CORR` is not written by the pipeline **[I found no mention of it in [PM]]**.

**(b) Verdict.** Diagonal errors are **wrong** for MATISSE in the same way as for GRAVITY, only more so: the oversampling is about 5 rather than 2, and the broadband terms are larger relative to the statistical ones.

**(c) Changes.**
1. Until channel covariances are modelled, bin to one resolution element before fitting (about 5 channels in L LOW for our files), and say why. `calibrate.py`'s `nbin` already does this; its docstring speaks of GRAVITY's twofold oversampling and needs a MATISSE line.
2. Use virgil's planned low-rank nuisances (Stage 6d; GRAVITY review §4c) with these modes per exposure: an achromatic gain per baseline (transfer function); a per-telescope mode on its three baselines, with m ∝ V²_model ε_p/F_i (photometry, [Lo22] eq. A.9); a chromatic per-telescope mode for coherence loss. All are generic and fit the existing design. virgil-vlti supplies the widths.
3. Treat the BCD positions of one cycle as separate frames that share atmosphere (same nuisances on their gains) but not instrument (separate instrumental closure phases). Then the BCD combination is a linear operator, and it can be propagated and whitened like any other (GRAVITY review footgun 3).

## 6. Closure phases

**(a) Literature.**
- Four triangles are written, three independent (our files: 4 rows). The pipeline's closure phase is the argument of the bias-corrected bispectrum ([PM] §11.5), or, from 2.5.0 with `cpAlgo=2`, built coherently ([PM] §2.1).
- BCD combination removes the instrumental closure phase "from several degrees to less than one degree"; residual calibration errors are 0.16–0.30° in L and 0.29–1.75° in N ([Lo22] §2.4, §4.1.3, Table 5). The floor at high flux is 0.1° (L, N) and 0.2° (M) ([Lo22] Table 4).
- A strong deviation of calibrator closure phases beyond 3.9 µm comes from the chromatic OPD in the VLTI tunnels, and the BCD calibration removes it ([Lo22] §4.1.3).
- Errors: the closure-phase error ≈ √3 σ_φ for S/N > 1 ([Lo22] eq. A.6), i.e. built from three baseline phase errors. That is the structure behind the ±1/3 correlation between triangles that share a baseline ([K20], for GRAVITY) **[the analogy is my inference; no MATISSE paper states the correlation]**.
- Isbell et al. used the standard deviation over the four BCD positions as the error (about 1° on the calibrator, 15° on Circinus) ([Is22] §2.2).

**(b) Verdict.** virgil's whitened, independent closure phases (virgil#120) **agree** for the noise part. Two caveats. The BCD residual is triangle-based by nature (it is what does not close), so it need not follow the ±1/3 structure. And if a chopped-only merge contains only IN-IN and OUT-OUT (footgun 7; [PM] Table 4.1), the four-position combination of [PM] §11.9 cannot be formed, and the instrumental closure phase is only partly removed **[my inference from the manual's sequence and formula; not checked on data]**.

**(c) Changes.**
1. For closure phases, prefer the non-chopped exposures with all four BCD positions, or all exposures together, over a chopped-only merge.
2. Report the closure-phase BCD scatter as a diagnostic, per triangle and channel, and use it to seed a closure-phase floor (virgil-vlti), as Isbell et al. did.
3. Leave closure offsets out of the default model, as the GRAVITY review §13 recommends; the BCD is MATISSE's own non-closing-error correction.

## 7. Spectral resolution modes and smearing

**(a) Literature.** Resolutions: LOW-LM and LOW-N R ≈ 31.5; MED-LM 499; HIGH-L 979; HIGH-N 218; VHIGH-L/M 3370 ([Lo22] Table 2). Varga et al. quote λ/δλ ≈ 34 in L and 30 in N ([Va21] §2). Without a fringe tracker, MED and HIGH read only 0.2 and 0.1 µm of the L band ([Lo22] §2.1). Spatial filters: a pinhole of 1.5λ/D in L&M and 2λ/D in N ([Lo22] §2.4, Table 2 note). No paper I read models bandwidth or time smearing for MATISSE; Gámez Rosas et al. note that spectral binning "can reduce the visibility contrast when the interferometric phase varies rapidly with wavelength" ([GR22] §4.3.1).

**My numbers** (top-hat bandpass, point companion; contrast loss sinc(πBθ/(λR)); 5% loss at θ ≈ 0.18 Rλ/B) **[my arithmetic, not from the literature]**:

| Mode | λ | Rλ/B at B = 30 m | Rλ/B at B = 130 m | θ for 5% loss at 130 m | Pinhole (AT / UT) |
|---|---|---|---|---|---|
| LOW L | 3.5 µm | 720 mas | 170 mas | ≈ 30 mas | 600 / 130 mas |
| MED L | 3.5 µm | 12″ | 2.8″ | ≈ 0.5″ | 600 / 130 mas |
| LOW N | 10 µm | 2.1″ | 480 mas | ≈ 85 mas | 2.3″ / 0.5″ |
| HIGH N | 10 µm | 15″ | 3.5″ | ≈ 0.6″ | 2.3″ / 0.5″ |

The pinhole column is the full diameter 1.5λ/D (L) or 2λ/D (N) with D the telescope diameter **[my reading of [Lo22] Table 2's note]**. The spectral average of 5 pixels makes the effective bandwidth wider than R ≈ 30 implies **[by how much, I have not computed]**.

Time smearing: a coherent average over T moves u by up to BωT. For a merged product spanning 10 min on a 130 m baseline, Δu ≈ 5.7 m, and a component 30 mas away sweeps up to 1.5 rad of phase at 3.5 µm. For one 1-min exposure it is 0.15 rad **[my estimate]**.

**(b) Verdict.** virgil's B·θ/λ-versus-R rule (imaging plan Stage 4) **agrees** in spirit. In LOW L on ATs it matters for anything wider than about 30 mas on the long baselines, inside the pinhole field. The two MATISSE-specific traps: `EFF_BAND` is the channel spacing, not the resolution (footgun 4), and merged products average over many minutes and over BCD positions (footgun 10).

**(c) Changes.**
1. virgil core's smearing (Stage 6a) takes an explicit spectral response per dataset, not `EFF_BAND`. virgil-vlti supplies it for MATISSE: R per mode, widened by the spectral average.
2. virgil core adds time smearing for long or merged integrations (GRAVITY review §8c.2); for MATISSE, fit per-exposure products when wide components matter.
3. virgil-vlti records the pinhole size per band and telescope type, for the generic primary beam (GRAVITY review §1c.2); MATISSE's beam is a pinhole, not a single-mode fibre.

## 8. `virgil_vlti.calibrate`: is it appropriate for MATISSE?

**(a) What it does.** Per night, each calibrator file gives a transfer function per baseline and channel: V²/V²_UD, averaged over the file's rows; closure phases as unit phasors. It interpolates in time (linear, or GRAVITY's Gaussian weights), adds a leave-one-out or scatter error in quadrature per channel, and classifies files by `ESO PRO CATG` or by name (`src/virgil_vlti/calibrate.py`, `products.py`).

**(b) Verdict.** The core idea **is appropriate** for MATISSE L&M V² and closure phases, *if* it is fed per-exposure `TARGET_RAW_INT` and `CALIB_RAW_INT` products. It must not be fed our files, which are already calibrated (`TARGET_CAL_INT`) and merged. As written it is **wrong** for MATISSE in five places:
1. **Classification** (footgun 13): `TARGET_CAL_INT` is labelled a calibrator, and both raw categories raise.
2. **BCD and chopping** (footgun 6): estimates are keyed on the `STA_INDEX` tuple only, so a calibrator at one BCD position calibrates science at another. Where a position reverses a baseline's orientation, the key differs, and the baseline would be flagged "no calibrator data" by accident rather than by design.
3. **Combination**: no BCD combination of closure phases or differential phases ([PM] §11.9).
4. **N band**: no correlated-flux calibration, and so no use of calibrator spectra or airmass.
5. **Defaults**: `nbin` advice assumes GRAVITY's twofold oversampling; there is no default wavelength window (the pipeline flags outside 2.8–4.2 and 4.48–5.2 µm in L&M, and outside 8–9.18 and 9.56–13 µm in N [PM §2.12]).

The interpolation choice is less critical. The pipeline averages by default and the consortium distrusts its linear mode ([PM] §9.22; [MT]); our linear and Gaussian interpolators are reasonable. The leave-one-out error is better than either. It should be computed per BCD group, and labelled as broadband (it is common to all channels), as the module docstring already says.

**(c) Changes.**
1. `classify`: recognise MATISSE categories explicitly (`CALIB_RAW_INT` → CAL, `TARGET_RAW_INT` → SCI) and refuse `*_CAL_INT` products as input, since they are already calibrated. Test this on headers, not on substrings.
2. Group calibrators as `mat_autoCalib` does: same band, resolution, DIT, BCD configuration and chopping status ([MT]). Calibrate within each group.
3. Normalise orientation before matching: sort each baseline's `STA_INDEX` and negate phases and u, v as needed, and do the same for triangles. This is the GRAVITY-independent version of `mat_removeBCD`.
4. Combine BCD positions after calibration, with the manual's signs for closure phases and the two correctable differential phases ([PM] §11.9). Propagate errors through the combination (it is linear), and report the scatter between positions as well.
5. Add correlated-flux calibration for N, following [Va21] eqs. 2–3: calibrator spectra (MDFC, the MIDI list, or Cohen standards), uniform-disk visibility, and an optional airmass correction. Propagate the calibrator's diameter and spectrum errors, or at least report them (PLAN item 1).
6. Read the `TF2` table of calibrator products, to compare our transfer function with the pipeline's.
7. A MATISSE validation, as for GRAVITY: our calibration with `--tf-error scatter` against `mat_cal_oifits` with `tfInterp=0`, per BCD position, on the same raw products (in a private validation suite).

## 9. What virgil core needs

All of these are generic, and most are already planned for GRAVITY.

1. **`AMPTYP` handling** in `read_oifits` (§1c.1). New, and small.
2. **Smearing with an explicit spectral response** and **time smearing** (§7c). Planned as Stage 6a; the MATISSE case makes the response argument essential, since `EFF_BAND` cannot be trusted.
3. **Low-rank calibration nuisances** (Stage 6d), with per-telescope modes. MATISSE's photometric error (eq. A.9 of [Lo22]) is exactly a per-telescope mode on three baselines.
4. **Channel covariance** for oversampled spectra: either a banded (sliding-average) covariance per spectrum, or a documented binning rule. A banded covariance is simple: a moving average of width w over independent noise gives correlation (w − k)/w at lag k **[standard result; it ignores the instrument's own line-spread function]**.
5. **Differential phases** through the normalisation operator of the GRAVITY review §2, with a default basis of constant plus 1/λ for MATISSE (§2c.2).
6. **Correlated flux** as an observable, once the total-spectrum machinery exists (§2c.1).
7. **A generic pinhole primary beam**, beside the fibre beam (§7c.3).

## 10. Recommended work items

**virgil-vlti**
1. Fix `classify` for MATISSE product categories (footgun 13; §8c.1). Small; do it first.
2. A MATISSE header reader: band, resolution, BCD, chopping, `corrFlux`, `spectralBinning`, `replaceTel`, pipeline versions, merge detection (§1c.2). Refuse `corrFlux` products as V², and flag N-band products from pipelines 1.5.2–1.x for the phase sign (footguns 1, 3).
3. Per-BCD, per-chopping-status calibration with orientation-normalised matching, then BCD combination with error propagation (§8c.2–4).
4. N-band correlated-flux calibration with calibrator spectra and airmass, and the calibrator-spectrum nuisance (§3, §4c.2, §8c.5).
5. MATISSE entries in the nuisance-width table (PLAN item 8): transfer-function, photometric and closure-phase widths per band and τ₀, from [Lo22] Tables 5–6 and the [Va21] floors (§4c.1).
6. Instrument facts: R per mode and the effective response after spectral averaging; pinhole sizes; default wavelength windows (§7c.1, §7c.3, §8b.5).
7. Validation against `mat_cal_oifits` in a private validation suite (§8c.7).
8. Update `matisse.md` to point here, and the `nbin` docstring (§5c.1).

**virgil core**
1. `AMPTYP` checks in `read_oifits` (§9.1).
2. Smearing with an explicit spectral response, and time smearing (§9.2).
3. Low-rank nuisances with per-telescope modes (§9.3), then channel covariances or a binning rule (§9.4).
4. Differential phases and correlated flux (§9.5–9.6).
5. A pinhole primary beam (§9.7).

## 11. Sources and how far each was checked

Read in full text (PDF converted to text; the sections listed were read):
- **[Lo22]** Lopez et al. 2022, A&A 659, A192, *MATISSE, the VLTI mid-infrared imaging spectro-interferometer*, arXiv:2110.15556: §2.1–2.5, §3.1–3.3, §4.1–4.2, §6, App. A. Tables 2 and 4–7.
- **[PM]** ESO, *MATISSE Pipeline User Manual*, VLT-MAN-MAT-15860-9310, issue 2.5.0 (2026-08-06), from ftp.eso.org: §2 (release notes 1.4.1–2.5.0), §4.1–4.2, §7 (product table), §7.5, §8.2, §9.1, §9.5, §9.10–9.15, §9.21–9.22, §10, §11.3–11.9, App. D. Our files were reduced with 2.0.0 and 1.7.9; I used the release notes to bridge the difference.
- **[Mi16]** Millour et al. 2016, SPIE 9907, *Data reduction for the MATISSE instrument*, arXiv:1608.01913: §2.4–2.7. This predates commissioning and describes the design ("we envision"); where it and [PM] differ, I follow [PM].
- **[Va21]** Varga et al. 2021, A&A (HD 163296), arXiv:2012.05697: §2, §3.1–3.2, App. B. **[Journal volume not checked; I believe it is A&A 647, A56.]**
- **[GR22]** Gámez Rosas et al. 2022, Nature (NGC 1068), arXiv:2112.13694: Supplementary §2–5 (instrument, observations, flagging, reduction, binning, calibration, photometry, faint-source tests, image orientation). **[Journal volume not checked.]**
- **[Is22]** Isbell et al. 2022, A&A 663, A35 (Circinus, N band), arXiv:2205.01575: §2.1–2.2, §2.5. The journal reference is from the VizieR catalogue record J/A+A/663/A35.
- **[MT]** MATISSE consortium tools, github.com/Matisse-Consortium/tools, read 2026-10-04: `mat_autoCalib.py` (calibrator selection and TF modes), `libPostTools.py` (`mat_mergeOifits`, `mat_removeBCD`, `mat_mergeByTplStart`), `mat_mergeAllOiFits.py`, `mat_autoPipeline.py` (default spectral binning), `mat_fluxCal.py` and `libFluxCal_STARSFLUX_new.py` (flux and correlated-flux calibration). Code, not a publication; its behaviour may change.
- **Our files.** Headers and table layouts (column names, shapes, keywords, row order) of the three ν Hor products; a few single-number summaries (median VISAMP/√V², median T3AMP, median relative errors). No data arrays were inspected beyond that.

From the GRAVITY review, not re-read here: **[D17]** Duvert, Young & Hummel 2017 (OIFITS 2); **[K20]** Kammerer et al. 2020 (GRAVITY correlations); **[La21]** Lachaume 2021.

Cited only through the sources above, not read: Cruzalèbes et al. 2019 (MDFC); van Boekel 2004 (MIDI calibrators); Millour et al. 2008 (AMBER BCD); Petrov et al. 2007 and 2016; Schutz et al. 2016 (OPD estimation); Koresko et al. 2006 (K-to-N co-phasing); Voronin & Zheltikov 2017; Jaffe 2004 (EWS); Hofmann et al. 2014, 2016 (IRBis).

Not covered: Robbe-Dubois et al. and Matter et al. (asked for in the brief; I did not find a paper of theirs that describes calibration in more depth than [Lo22], and did not read one); Berio et al. on N-band coherent processing ("in preparation" in [Lo22]); the pipeline source code, so the exact error propagation in `mat_cal_vis` and `mat_cal_cphase` is unverified; MATISSE wavelength-calibration accuracy (I found no quoted number, so no wavelength prior is recommended); an N-band product, so N-band `INSNAME` and `AMPTYP` values are unverified; GRA4MAT-specific products beyond the manual's §11.7–11.8.
