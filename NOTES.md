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

## Experiments & Results

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
