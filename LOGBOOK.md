# MIMM Internship Logbook

---

## Input Data Description

All input data comes from the Zenodo example release (record 10019720), accompanying Şişman et al. 2025. It is a single healthy subject scanned on a 3T Siemens Prisma scanner. All volumes share the same space: **208 × 256 × 160 voxels, 1 mm isotropic**.

---

### iField.mat

The core MRI input. Contains the following variables:

| Variable | Content |
|---|---|
| `iField` | 4D array [208 × 256 × 160 × 7] — mGRE magnitude signal, one 3D volume per echo time. Values range 0 to ~1279 (arbitrary scanner units, normalised per voxel before matching). |
| `TE` | 7 echo times in seconds: 3.32, 8.81, 14.30, 19.79, 25.28, 30.77, 36.26 ms |
| `delta_TE` | Echo spacing = 5.49 ms |
| `CF` | Centre frequency = 123,260,884 Hz = 123.26 MHz → confirms 3T scanner (proton Larmor frequency at 3T ≈ 123.2 MHz) |
| `voxel_size` | [1, 1, 1] mm — isotropic 1 mm resolution |
| `matrix_size` | [208, 256, 160] |
| `B0_dir` | [0, 0, 1] — B0 field points along the z-axis (standard axial acquisition) |
| `files` | File path metadata from the original acquisition |

The 7 echoes capture how fast the signal decays in each voxel. A voxel with myelin or iron decays faster (higher R2*) than a voxel with only water. MIMM uses these 7 time points to match against the dictionary signal curves.

**Note:** The MUMC acquisition uses 5 echoes at 6, 12, 18, 24, 30 ms (TE1 = 6 ms, ΔTE = 6 ms). The echo spacing is roughly the same order; the interpolation step in MIMM.m handles the difference automatically.

---

### QSM.mat

Contains two variables:

| Variable | Content |
|---|---|
| `QSM` | 3D array [208 × 256 × 160] — quantitative susceptibility map in ppm. Range: −4.14 to +2.85 ppm. Negative = diamagnetic (myelin), positive = paramagnetic (iron, veins). |
| `Brain_Mask` | 3D binary array [208 × 256 × 160] — 1 inside the brain, 0 outside. Contains 1,150,162 brain voxels (~13% of the full volume). Used by MIMM to restrict processing to brain tissue only. |

QSM is reconstructed from the mGRE phase data (not provided separately). It encodes the magnetic susceptibility of each voxel, which reflects the combined contribution of iron (positive) and myelin (negative). MIMM uses this alongside the magnitude signal to separate the two sources.

The range −4.14 to +2.85 ppm is wider than typical brain tissue values (usually −0.1 to +0.3 ppm) because it includes blood vessels, air cavities, and background field effects at brain edges.

---

### FA.mat

Fractional anisotropy maps derived from DTI. Contains two variables:

| Variable | Content |
|---|---|
| `FA_DTI` | 3D array [208 × 256 × 160] — subject-specific FA from a DTI scan. Range: 0 to 1.2 (values slightly above 1.0 can occur in noisy voxels due to the FA calculation). High FA = coherent fibre structure; FA > 0.25 is the threshold MIMM uses to identify reliable white matter. |
| `FA_atlas` | 3D array [208 × 256 × 160] — FA derived from the ICBM DTI-81 normative atlas, warped into subject space. Range: 0 to 0.769. Smoother than subject DTI due to population averaging. |

FA is used in orientation-informed MIMM to define where the orientation prior is reliable. MIMM only applies the DTI/atlas orientation in voxels where `FA > 0.25` **and** `QSM < 0.1 ppm` (to exclude iron-rich basal ganglia where DTI orientation estimates are unreliable).

---

### theta.mat

Fibre orientation angle maps, one per orientation source:

| Variable | Content |
|---|---|
| `theta_DTI` | 3D array [208 × 256 × 160] — angle in degrees between the principal fibre axis and B0. Derived from the subject's own DTI scan. Range: 0° to 90°. |
| `theta_atlas` | 3D array [208 × 256 × 160] — same angle derived from the ICBM atlas diffusion tensor warped into subject space. Range: 0° to 90°. |

These angles directly parameterise the anisotropic susceptibility and relaxation effects in the HCFM signal model. A fibre running perpendicular to B0 (θ = 90°) has a different magnitude signal evolution than one running parallel (θ = 0°). Without this prior, basic MIMM estimates θ freely, leading to noisy orientation maps and some MVF overestimation in major tracts.

---

## Summary table

| File | Variables | Size | What it represents |
|---|---|---|---|
| `iField.mat` | `iField` (4D), `TE`, `delta_TE`, `CF`, `voxel_size`, `matrix_size`, `B0_dir` | 208×256×160×7 | mGRE magnitude signal at 7 echo times — the primary MRI measurement |
| `QSM.mat` | `QSM`, `Brain_Mask` | 208×256×160 | Susceptibility map (iron/myelin content) + brain mask |
| `FA.mat` | `FA_DTI`, `FA_atlas` | 208×256×160 | Fibre coherence — used to define where orientation prior is reliable |
| `theta.mat` | `theta_DTI`, `theta_atlas` | 208×256×160 | Fibre–B0 angle — the orientation prior for orientation-informed MIMM |

---

## Experiment Log

### 2026-04-20 — Week 16 complete

**Goal:** Run all MIMM variants on the Zenodo example data, verify outputs, set up GitHub repo.

**What was done:**
- Ran stochastic dictionary (20,000 entries): Basic, DTI, Atlas MIMM → saved to `stochastic_MIMM_results_new.mat`
- Ran deterministic dictionary (12,540 entries): Basic, DTI, Atlas MIMM → saved to `deterministic_MIMM_results_new.mat`
- Generated and saved PNG figures for all 8 output maps × 2 dictionaries
- Created private GitHub repo at github.com/Tijn808/MIMM_V1
- Pushed all code, scripts, figures, README, and NOTES

**Output ranges (stochastic, basic MIMM):**

| Map | Min | Max | Mean |
|---|---|---|---|
| MVF | 0.0 | 0.556 | 0.140 |
| g_ratio | 0.500 | 1.000 | 0.853 |
| FVF | 0.037 | 0.743 | 0.525 |
| R2s | 14.4 | 83.3 | 21.3 |
| chi_iron_est | 0.0 | 0.289 | 0.019 |
| chi_myelin | −0.056 | 0.0 | −0.014 |
| theta_est | 0.001 | 89.98 | 30.0 |
| error | 0.0 | 0.240 | 0.001 |

**Validation against paper Table 1:**

| Parameter | Paper range | Our output | Match |
|---|---|---|---|
| MVF | [0, 0.55] | [0.0, 0.556] | ✓ |
| g-ratio | [0.5, 1] | [0.500, 1.000] | ✓ |
| FVF | [0, 0.75] | [0.037, 0.743] | ✓ |
| χ^iron | [0, 300 ppb] | [0, 289 ppb] | ✓ |
| θ | [0°, 90°] | [0.001°, 89.98°] | ✓ |

**Differences vs Zenodo reference results:**
- Stochastic: max MVF difference = 0.26 — expected, different random seed each run
- Deterministic: max MVF difference = 0.4075 — caused by removal of `phi = imrotate(phi, 90)` in `compute_field.m` after reference was computed

**Key observations from figures:**
- MVF: clear WM anatomy, all three modes nearly identical
- theta_est: most striking difference between modes — Basic is random noise, DTI and Atlas show coherent tract structure
- chi_myelin: anatomically sharpest map, clear corpus callosum and internal capsule
- error: near zero everywhere, confirming good dictionary fit
- Quantitative Bland-Altman comparison (paper Figure 4) not possible — requires FAST-T2 MWF reference data not included in Zenodo release

**Known issues / code notes:**
- `dictionary` is a reserved class name in MATLAB R2025b — always load as `stoch = load(...); dict = stoch.dictionary;`
- Image Processing Toolbox required for `medfilt3` in orientation-informed mode
