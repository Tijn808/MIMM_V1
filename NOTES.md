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

### JHU atlas registration — critical assessment (2026-05-28)

Both JHU atlases (DTI-81, 48 ROIs; tractography thr25, 20 tracts) have been warped to MUMC subject space using the existing `mni2subj_warp.nii.gz`. Registration geometry is correct — R/L hemisphere assignments are correct and anatomical centroids are in plausible positions (genu anterior to splenium, IC correctly lateralised, paired tracts on correct sides).

**Key limitation — FA_atlas is a population average:**
The HCP1065 FA atlas warped to subject space has mean FA = 0.191 across all brain voxels (FA > 0.1). This is much lower than individual subject FA (~0.4–0.7 in WM) due to inter-subject averaging. Consequence: the MIMM atlas orientation prior (FA > 0.25 threshold) activates in far fewer voxels than intended — the atlas mode runs largely as basic mode in this subject. This will be resolved when subject-specific DTI FA arrives.

**WM masking for figures:** FA_atlas > 0.20 used as WM display mask (180k voxels, primarily WM). Not a perfect segmentation — proper WM mask requires FSL FAST on T1w or subject-specific DTI FA. The chi_myelin map in particular still shows signal in some GM-adjacent voxels because MIMM has no GM class in the dictionary (forced match for all brain voxels).

**Atlas comparison (DTI-81 vs tractography thr25):**
- Both atlases show similar within-ROI CV (0.57) — neither is clearly more homogeneous
- Tractography atlas shows marginally higher Pearson r for MVF vs chi_neg (0.77 vs 0.75) — slightly better spatial agreement between methods
- DTI-81 is preferred for the primary analysis (matches the paper, more granular, bilateral pairs); tractography thr25 used as sensitivity check

**Pending:** Re-run registration with subject-specific DTI FA when available; use DTI FA > 0.25 as definitive WM mask.

### Figure quality fixes (2026-05-28)

Applied three fixes to `generate_figures.py`:
1. **WM mask (FA_atlas > 0.20)** applied to MVF, FVF, g-ratio, chi_myelin, chi_iron, theta — removes GM false-positive signal from microstructure maps
2. **R2* colorscale** reduced from 0–50 to 0–35 s⁻¹ — better WM contrast, basal ganglia iron visible without saturation
3. **FA atlas colorscale** reduced from 0–1 to 0–0.5 — reveals actual WM FA structure in the population-average map

Residual issues (not fixable without better data):
- QSM streaking artifact still visible in MVF difference map (coronal view) — label as artifact in any caption
- chi_myelin still shows some GM coverage — fundamental MIMM limitation (no GM dictionary class)
- FA atlas remains faint — will improve with subject DTI

---

## Experiments & Results

### 2026-05-28 — Per-ROI MVF overestimation: basic vs atlas (MUMC subject)

**Data:** Single MUMC subject (Philips 3T, TEs 6/12/18/24/30 ms).

**Key finding: physics hypothesis confirmed — r=0.66, p<0.0001**

Overestimation (basic − atlas) correlates strongly with mean fibre angle θ per ROI.
Basic MIMM overestimates MVF most in tracts running perpendicular to B0 because without
an orientation prior it attributes angle-dependent ME-GRE signal variation to myelin content.

The regression line crosses zero at θ ≈ 30–35°:
- θ < 30° (CST, medial lemniscus — parallel to B0): near-zero or slight underestimation
- θ > 60° (CC, posterior thalamic radiation — perpendicular): systematic overestimation

**Top overestimated ROIs:**
| Tract | Basic | Atlas | Overest (abs) | Overest (%) | θ (°) |
|---|---|---|---|---|---|
| Splenium CC | 0.359 | 0.262 | +0.097 | +37% | 79.6 |
| Post. thal. radiation L | 0.330 | 0.258 | +0.072 | +28% | 72.3 |
| Post. thal. radiation R | 0.336 | 0.278 | +0.058 | +21% | 76.0 |
| Retrolenticular IC R | 0.329 | 0.277 | +0.052 | +19% | 62.6 |
| Body CC | 0.234 | 0.184 | +0.050 | +27% | 68.6 |
| Genu CC | 0.247 | 0.198 | +0.049 | +25% | 76.8 |
| Ant. limb IC L | 0.305 | 0.260 | +0.045 | +17% | 75.2 |

**Only underestimated:** Medial lemniscus L (−7%) — runs parallel to B0 (θ=18°).

**Notable outlier:** Splenium sits above the regression line — the combination of very high
MVF (~0.36) AND near-perpendicular orientation (~80°) amplifies the absolute error beyond
what theta alone predicts.

**Scatter at high theta:** Several tracts at θ>70° show lower overestimation than expected.
These are regions where FA_atlas < 0.25 — the atlas orientation prior never activates,
so basic and atlas give identical results. This is a limitation of the atlas mode, not a feature.

**Additional correlation:** overestimation vs mean FA: r=0.58, p<0.0001.
High-FA tracts tend to be coherently oriented and thus most affected by orientation error.

Figures: `23_overestimation_ranked.png`, `24_overestimation_vs_theta.png`, `25_overestimation_vs_FA.png`

---

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
