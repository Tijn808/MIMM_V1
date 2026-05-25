# MUMC MIMM Pipeline Documentation

**Project:** Myelin Imaging using Multi-echo GRE and MIMM  
**Dataset:** MUMC IMPROMYMS (MS patients, Philips MR 7700, 3T)  
**Last updated:** 2026-05-13

---

## Acquisition Parameters

| Parameter | Value |
|---|---|
| Sequence | 3D ME-GRE, monopolar readout |
| Scanner | Philips MR 7700, 3T |
| Echo times | 6.001, 12, 18, 24, 30 ms |
| TR | 35 ms |
| Voxel size | 1 × 1 × 1 mm isotropic |
| Matrix | 224 × 224 × 154 |
| B0 direction | [0, 0, 1] |

---

## Full Pipeline Overview

```
Raw Philips DICOMs
        │
        ▼
preprocess/prepare_mgre.m        — Philips RWV phase scaling → NIfTI
        │
        ▼
qsm/run_QSM_chisep.m             — ROMEO + V-SHARP + iLSQR → QSM, R2*, brain mask
        │
        ▼
atlas/register_atlas.sh          — HCP1065 atlas → subject space (fnirt)
        │
        ▼
mimm/run_MIMM_MUMC.m             — Dictionary matching → MVF, FVF, g-ratio, ...
```

---

## Step 1: Preprocessing (`preprocess/prepare_mgre.m`)

Converts raw Philips DICOM data to NIfTI format.

**Key operation — Philips RWV phase scaling:**  
Philips stores phase as integers. The real-world value (RWV) transformation applies:
```
phase_rad = raw_value × RWV_slope + RWV_intercept
```
This converts the stored integers to radians in [−π, π], read from the JSON sidecar.

**Outputs (written to `subj_dir/`):**
- `magnitude.nii.gz` — 4D ME-GRE magnitude [224 × 224 × 154 × 5]
- `phase.nii.gz` — 4D ME-GRE phase in radians [224 × 224 × 154 × 5]

---

## Step 2: QSM Reconstruction

### Pipeline v1 — SEPIA-based (`qsm/run_QSM.m`) *(deprecated)*

| Step | Method |
|---|---|
| Brain mask | FSL BET (f=0.4) |
| Field map | Phase-difference between adjacent echoes, weighted average |
| Background removal | SEPIA BKGRemovalVSHARP (radii 1–10 mm) |
| Dipole inversion | SEPIA qsmIterativeLSQR |

### Pipeline v2 — Chi-sep toolbox (`qsm/run_QSM_chisep.m`) *(current)*

Supervisor's recommended pipeline. Uses the chi-separation toolbox (Shin et al. 2021).

| Step | Method | Toolbox |
|---|---|---|
| Brain mask | MEDI BET | MEDI toolbox |
| R2* fitting | ARLO | MEDI toolbox |
| Phase unwrapping | ROMEO + weighted echo averaging | mritools (Ubuntu 24.04) |
| Background removal | V-SHARP (smvsize=12) | STI Suite V3.0 |
| Dipole inversion | iLSQR (padsize=[12,12,12]) | STI Suite V3.0 |

**Toolbox setup:** Run `Chisep_Toolbox_v1.2.1/setup_toolbox_paths.m` once before first use.

**Key fixes applied during setup:**
- Tukey windowing disabled (`Tukey=0`) — Philips data, Signal Processing Toolbox not available
- NIfTI header: clean dummy header used (not inherited from magnitude) to avoid Philips `scl_slope=739` corrupting FSL display

**Outputs (written to `subj_dir/qsm/`):**
- `QSM.nii.gz` — susceptibility map (ppm)
- `R2star.nii.gz` — R2* map (s⁻¹)
- `brain_mask.nii.gz` — binary brain mask
- `mag_e1.nii.gz` — echo 1 magnitude (used by atlas registration)

### QSM Pipeline Comparison (sub-01, brain mask, P5/P50/P95)

| Metric | Pipeline v1 (SEPIA) | Pipeline v2 (chi-sep) |
|---|---|---|
| QSM P5 | −0.091 ppm | −0.054 ppm |
| QSM P50 | −0.001 ppm | −0.002 ppm |
| QSM P95 | +0.096 ppm | +0.060 ppm |

Both pipelines produce values within expected range for 3T brain (WM: −0.03 to −0.07 ppm, iron-rich structures: +0.05 to +0.15 ppm). Pipeline v2 gives a slightly narrower distribution, consistent with ROMEO producing cleaner phase unwrapping than the phase-difference approach.

---

## Step 3: Atlas Registration (`atlas/register_atlas.sh`)

Registers the HCP1065 population-average DTI atlas into subject ME-GRE space. Used when no subject-specific DTI data is available.

**Atlas source:** FSL standard — `FSL_HCP1065_FA_1mm.nii.gz`, `FSL_HCP1065_V1_1mm.nii.gz`

**Registration pipeline:**
1. Brain-extract echo 1 magnitude using brain mask
2. Affine registration: subject → MNI152 1mm (FSL `flirt`, 12-DOF, normcorr)
3. Nonlinear registration: `fnirt` (bending_energy, λ=300, warpres=10mm, 4 levels)
4. Invert warp: MNI → subject (`invwarp`)
5. Warp FA atlas to subject space (`applywarp`, trilinear)
6. Warp V1 eigenvector atlas with reorientation (`vecreg` — corrects vector directions per local Jacobian)
7. Compute fibre angle: `theta = arccos(|V1_z|) × 180/π`

**Outputs (written to `subj_dir/atlas/`):**
- `FA_atlas.nii.gz` — FA map in subject space
- `V1_atlas.nii.gz` — principal eigenvector [x, y, z] in subject space
- `theta_atlas.nii.gz` — fibre angle relative to B0 (degrees)

---

## Step 4: MIMM Reconstruction (`mimm/run_MIMM_MUMC.m`)

Runs MIMM (Magnetic Imaging of Myelin and Iron Microstructure, Şişman et al. 2025) in two variants.

**Parameters:**
- Dictionary: `MIMM_dictionary_stochastic.mat`
- λ_χ = 0.015 (QSM/magnitude weighting, L-curve optimised)
- Echo times: [6, 12, 18, 24, 30] ms

**Variants:**
- **Basic** — no orientation prior; MIMM estimates fibre angle from data
- **Atlas** — HCP1065 FA + theta used as orientation prior

**Outputs (written to `subj_dir/mimm/`):**

| Map | Unit | Description |
|---|---|---|
| `MVF_{tag}.nii.gz` | fraction | Myelin volume fraction |
| `FVF_{tag}.nii.gz` | fraction | Fibre volume fraction |
| `g_ratio_{tag}.nii.gz` | — | g-ratio (axon/fibre diameter) |
| `R2s_{tag}.nii.gz` | s⁻¹ | Fitted transverse relaxation |
| `chi_myelin_{tag}.nii.gz` | ppm | Myelin susceptibility contribution |
| `chi_iron_est_{tag}.nii.gz` | ppm | Iron susceptibility contribution |
| `theta_est_{tag}.nii.gz` | degrees | Estimated fibre angle |
| `error_{tag}.nii.gz` | — | Dictionary matching residual |

Where `{tag}` is `basic` or `Atlas`.

### MIMM Output Values (sub-01, whole-brain, P5/P50/P95/mean)

| Map | Basic | Atlas |
|---|---|---|
| MVF | 0.002 / 0.144 / 0.443 / 0.159 | 0.002 / 0.143 / 0.428 / 0.157 |
| FVF | 0.037 / 0.640 / 0.742 / 0.498 | 0.037 / 0.640 / 0.742 / 0.496 |
| g-ratio | 0.534 / 0.852 / 0.997 / 0.829 | 0.538 / 0.852 / 0.997 / 0.831 |
| R2s (s⁻¹) | 14.4 / 20.3 / 39.3 / 22.4 | 14.4 / 20.3 / 39.2 / 22.4 |
| chi_myelin (ppm) | −0.044 / −0.014 / 0.000 / −0.016 | −0.043 / −0.014 / 0.000 / −0.016 |
| chi_iron_est (ppm) | 0.000 / 0.011 / 0.078 / 0.019 | 0.000 / 0.011 / 0.078 / 0.019 |

### Sanity Check Against Literature (3T)

| Map | Observed (P50) | Expected | Status |
|---|---|---|---|
| MVF whole-brain | 0.14 | 0.10–0.18 (WM diluted by GM/CSF) | ✓ |
| FVF | 0.64 | 0.60–0.75 (WM-dominant) | ✓ |
| g-ratio | 0.85 | 0.70–0.90 | ✓ |
| R2s | 20 s⁻¹ | 14–40 s⁻¹ | ✓ |
| chi_myelin | −0.014 ppm | negative (diamagnetic) | ✓ |
| chi_iron_est | +0.011 ppm | positive (paramagnetic) | ✓ |

**Basic vs Atlas difference:** The atlas prior (HCP1065 orientation) produces a slightly tighter MVF distribution (P95: 0.428 vs 0.443) but nearly identical medians. FVF, g-ratio, R2s, and susceptibility maps are essentially unchanged — confirming that orientation information primarily affects MVF estimation, as expected from the MIMM biophysical model.

---

## File Structure

```
subj_dir/                       (e.g. ME_GRE/)
├── magnitude.nii.gz            4D ME-GRE magnitude
├── phase.nii.gz                4D ME-GRE phase (radians)
├── qsm/
│   ├── QSM.nii.gz              Susceptibility map (ppm)
│   ├── R2star.nii.gz           R2* map (s⁻¹)
│   ├── brain_mask.nii.gz       Binary brain mask
│   └── mag_e1.nii.gz           Echo 1 magnitude
├── atlas/
│   ├── FA_atlas.nii.gz         HCP1065 FA in subject space
│   ├── V1_atlas.nii.gz         HCP1065 eigenvector in subject space
│   └── theta_atlas.nii.gz      Fibre angle map (degrees)
└── mimm/
    ├── MIMM_results.mat        Full MIMM output struct
    ├── MVF_basic.nii.gz
    ├── MVF_Atlas.nii.gz
    ├── FVF_basic.nii.gz / FVF_Atlas.nii.gz
    ├── g_ratio_basic.nii.gz / g_ratio_Atlas.nii.gz
    ├── R2s_basic.nii.gz / R2s_Atlas.nii.gz
    ├── chi_myelin_basic.nii.gz / chi_myelin_Atlas.nii.gz
    ├── chi_iron_est_basic.nii.gz / chi_iron_est_Atlas.nii.gz
    ├── theta_est_basic.nii.gz / theta_est_Atlas.nii.gz
    └── error_basic.nii.gz / error_Atlas.nii.gz
```

---

## Toolbox Dependencies

| Toolbox | Version | Purpose |
|---|---|---|
| MIMM | — | Dictionary matching reconstruction |
| SEPIA | — | (v1 pipeline only) |
| Chi-separation toolbox | v1.2.1 | QSM preprocessing pipeline |
| MEDI toolbox | 2024.11.26 | BET brain extraction, ARLO R2* |
| STI Suite | V3.0 | V-SHARP background removal, iLSQR |
| SEGUE | 2025-07-03 | Phase unwrapping fallback (Linux) |
| mritools (Ubuntu 24.04) | 4.7.1 | ROMEO phase unwrapping |
| FSL | — | BET, atlas registration (flirt/fnirt/vecreg) |

---

## Known Issues / Notes

- **Philips scl_slope:** The Philips magnitude NIfTI has `scl_slope ≈ 740`. All output NIfTIs use a clean dummy header to prevent FSL from misinterpreting values.
- **Signal Processing Toolbox:** Not available — Tukey windowing is disabled (correct for Philips data, no effect on output).
- **mritools version:** ROMEO uses the Ubuntu 22.04 binary (bundled in Chisep toolbox) even on Ubuntu 24.04 — works correctly.
- **DTI:** No subject DTI available for this dataset. Atlas-based orientation (HCP1065) is used instead.
- **T2-GRASE:** Not yet processed. Required for MWF (myelin water fraction) and chi-separation. Data is available from supervisor.
