# MIMM Project Notes

## Overview
Implementation of Şişman et al. 2025 "Microstructure-Informed Myelin Mapping (MIMM)"
(Magn Reson Med 93:1499–1515, doi:10.1002/mrm.30369)

Estimates myelin volume fraction (MVF) from multi-gradient echo (mGRE) MRI data using dictionary matching.

## Environment
- MATLAB R2025b
- Requires Image Processing Toolbox (for `medfilt3` in orientation-informed mode)
- Run scripts via bash: `/path/to/matlab -batch "run('scripts/run_stochastic_all.m')"`

## Known Issues
- `dictionary` is a reserved class name in MATLAB R2022b+. Load `.mat` files into a different variable:
  ```matlab
  stoch = load('MIMM_dictionary_stochastic.mat');
  dict = stoch.dictionary;
  ```

## Code Changes vs Original (Zenodo)
Two files changed since the reference results were computed:
- `compute_field.m`: `phi = imrotate(phi, 90)` was removed — changes field orientation calculation
- `MIMM.m`: `lambda_chi` moved from hardcoded inside the function to an input parameter

Reference results from Zenodo cannot be reproduced exactly with the current code. Differences are expected.

## MUMC Pipeline Findings (2026-05-25)

### TE1 = 6.001 ms (not 6.000 ms)
The Philips MR 7700 reports TE1 = 6.001 ms in the JSON sidecar due to gradient raster clock quantization (1 µs steps). The difference is physically negligible (~0.003% signal error at typical WM R2* = 30 s⁻¹). Fixed in `prepare_mgre.m`: TEs are now read directly from the JSON sidecar per echo instead of hardcoded.

### Chi-separation accuracy limitation
The current chi-sep pipeline (ME-GRE only, no spin echo) uses R2* as a proxy for R2', which introduces ~57% error in separated susceptibility components (PMC12620178). This means chi_neg is spatially informative but quantitatively unreliable.

Root cause: within a 1mm³ voxel, diamagnetic myelin and paramagnetic iron partially cancel in QSM. Chi-separation requires R2' = R2* − R2 to disentangle them, but without a spin echo scan, R2 is unavailable.

**Potential fix:** T2-GRASE data (already collected for MWF) may provide a geometric mean T2 map (T2_gm), giving R2 = 1/T2_gm as input to chi-sep. This would be substantially more accurate than the R2* approximation, though not as clean as a dedicated dual-echo TSE with B1+ correction. Awaiting confirmation from Gerald whether T2_gm is output by the T2-GRASE pipeline.

### Streaking artifacts in QSM / chi-separation (confirmed by Gerald Drenthen, 2026-05-27)
The 'X'-shaped linear artifacts in chi-separation maps are **streaking artifacts** from the dipole inversion.

**Physical cause:** The dipole kernel in Fourier space is `d(k) = 1/3 − k_z²/|k|²`, which equals zero at the magic-angle cones (θ = 54.7° from B0). Inverting the forward model requires dividing by d(k); near-zero values amplify noise → artifacts back-project into real space as diagonal streaks forming an 'X' through high-susceptibility structures.

**What is affected:**
- QSM, chi_neg, chi_pos: directly affected — artifact is in the source map
- MIMM MVF: partially buffered — QSM enters the dictionary cost at λ_χ = 0.015 (small weight), so the magnitude fitting term usually dominates; MVF is less affected than chi_neg
- R2*: not affected — fitted from ME-GRE magnitude, no dipole inversion

**Options to reduce streaking:**

| Option | Effect | Feasibility |
|---|---|---|
| Multi-orientation QSM (COSMOS, ≥3 tilts) | Eliminates zeros in dipole kernel — exact inversion | Not feasible clinically |
| Better regularization (MEDI, iLSQR, morphology-enabled) | Suppresses streaks; doesn't eliminate for single orientation | Already in standard chi-sep pipelines |
| T2-GRASE R2 = 1/T2_gm for chi-sep | Fixes ~57% R2'/R2* error; **does not fix streaking** | Feasible — awaiting T2_gm confirmation from Gerald |

**Implication for Bland-Altman:** Streaking adds spatial bias and scatter to chi_neg, compounding the R2' error. Worst in iron-rich regions (basal ganglia, near veins). Reinforces that MVF vs chi_neg is a consistency check only.

### MVF vs chi_neg Bland-Altman — interpretation
Comparing MIMM MVF against chi_neg (chi-sep) is a **consistency check, not a validation**, because both methods derive from the same QSM input. True independent validation requires FAST-T2 MWF (T2-based, no QSM involvement). The chi_neg comparison is further limited by the ~57% quantitative error above. Expect larger Bland-Altman bias in iron-rich tracts (basal ganglia, deep WM) where iron-myelin cancellation is strongest.

### Subject-level Bland-Altman
The paper (Şişman et al. 2025) pools voxels across subjects, inflating N and producing artificially narrow LoA. The correct approach is one dot per subject per ROI. Implemented in `MUMC_pipeline/analysis/generate_bland_altman.py`.

### MVF → MWF conversion (TODO: implement before Bland-Altman vs T2-GRASE)
When T2-GRASE MWF data arrives, MIMM MVF must be converted to MWF units before
Bland-Altman comparison. The paper (Şişman et al. 2025, Eq. 10) defines the conversion as:

    MWF = (ρ_MW × MVF) / (ρ_IEW × (1 − MVF) + ρ_MW × MVF)

Constants from Table 1 (3T, from Hédouin 2021 and Xu 2018):
  ρ_MW  = 0.5   (myelin water proton density relative to free water)
  ρ_IEW = 1.0   (intra/extracellular water proton density)

Example: MVF = 0.20 → MWF ≈ 0.111 (roughly MVF ≈ 2× MWF, nonlinear at higher values).

This conversion uses fixed literature constants — the paper itself notes it is "imperfect"
and may contribute to the residual bias seen in atlas/DTI MIMM vs FAST-T2 MWF.
Implement in `generate_bland_altman.py` before the OLS scaling step, or as an alternative
to OLS scaling when a principled same-units comparison is preferred.

### DTI data
DTI files (.nii.gz, .bval, .bvec) not yet received. Confirmed b-values: 0, 1000, 2000. Preprocessing script ready at `MUMC_pipeline/preprocess/preprocess_dti.sh`. Awaiting account setup at MUMC.

---

## Experiments & Results

### 2026-05-26 — Per-ROI MVF overexpression: basic vs atlas (Zenodo example subject)

**Data:** Zenodo single example subject (Şişman et al. 2025, DOI 10.5281/zenodo.10019720).
All numbers below are from this one subject — not MUMC patient data.

**Method:**
- FA and brain mask extracted from `Example_Data/FA.mat` and `QSM.mat` (via Python/scipy)
- JHU ICBM-DTI-81 atlas registered to subject space: flirt (affine, 12-DOF, normcorr) → fnirt (nonlinear, bending-energy λ=300) → invwarp → applywarp (--interp=nn)
- MVF compared within each JHU ROI: mean(basic) − mean(atlas), relative to atlas mean

**Whole-brain WM summary (MVF > 0.1 mask):**
| Variant | Mean MVF | vs Basic |
|---|---|---|
| Basic | 0.236 | — |
| Atlas | 0.223 | −6.1% |
| DTI | 0.207 | −13.9% |

**Per-ROI: most overexpressed in basic vs atlas:**
| Tract | Basic | Atlas | Overexpression |
|---|---|---|---|
| Post. thalamic radiation L | 0.424 | 0.276 | +54% |
| Post. thalamic radiation R | 0.397 | 0.248 | +60% |
| Sagittal stratum L | 0.332 | 0.210 | +58% |
| Sagittal stratum R | 0.352 | 0.239 | +47% |
| Splenium CC | 0.324 | 0.216 | +51% |
| Retrolenticular IC R | 0.403 | 0.312 | +29% |
| Pontine crossing tract | 0.345 | 0.248 | +39% |

**Nearly unaffected (B−A ≈ 0):** CST, posterior IC, superior CR, medial lemniscus, tapetum.

**Slightly underestimated by basic:** Cerebral peduncle R/L (−3–6%), medial lemniscus R/L (−3–4%).

**Interpretation:** Basic MIMM overestimates MVF most in tracts running perpendicular to B0.
Without an orientation prior it attributes angle-dependent ME-GRE signal variation to myelin
thickness → inflated MVF. Atlas mode corrects this via the HCP1065 fibre orientation prior
(active only where FA > 0.25 and QSM < 0.1 ppm). DTI mode corrects more strongly because
it uses the actual measured per-voxel fibre orientation rather than a population atlas.

**What "percentage" means:** (basic_mean − atlas_mean) / atlas_mean × 100%.
Not a percentage of true myelin — atlas is assumed better, not ground truth.
True validation requires histology or independent FAST-T2 MWF.

**Note on FSL on this Mac:** FSL installed via conda (FSL conda channel + conda-forge).
`register_atlas.sh` updated to use `$FSLDIR` env variable, so it works on both this Mac
and the MUMC Linux server without code changes.

### 2026-04-20 — All three MIMM variants, stochastic dictionary
- **Dictionary:** stochastic (20,000 entries)
- **Modes:** basic, DTI orientation informed, atlas orientation informed
- **Output:** `Example_Results/stochastic_MIMM_results_new.mat`
- **Comparison vs Zenodo reference:** max MVF difference = 0.26
- **Reason:** Stochastic dictionary uses random sampling — different seed from reference

### 2026-04-20 — All three MIMM variants, deterministic dictionary
- **Dictionary:** deterministic (12,540 entries)
- **Modes:** basic, DTI orientation informed, atlas orientation informed
- **Output:** `Example_Results/deterministic_MIMM_results_new.mat`
- **Comparison vs Zenodo reference:** max MVF difference = 0.4075
- **Reason:** `compute_field.m` was updated after reference was computed

### Output map ranges (stochastic, basic MIMM)
| Map | Min | Max | Mean |
|---|---|---|---|
| MVF | 0.0 | 0.556 | 0.140 |
| g_ratio | 0.500 | 1.000 | 0.853 |
| FVF | 0.037 | 0.743 | 0.525 |
| R2s | 14.4 | 83.3 | 21.3 |
| chi_iron_est | 0.0 | 0.289 | 0.019 |
| chi_myelin | -0.056 | 0.0 | -0.014 |
| theta_est | 0.001 | 89.98 | 30.0 |
| error | 0.0 | 0.240 | 0.001 |
