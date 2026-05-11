# MIMM, Microstructure-Informed Myelin Mapping

Implementation of **Şişman et al. 2025**, *Microstructure-Informed Myelin Mapping (MIMM) from routine multi-echo gradient echo data using multiscale physics modeling of iron and myelin effects and QSM*
(Magn Reson Med 93:1499–1515, doi:10.1002/mrm.30369)

MIMM estimates the **myelin volume fraction (MVF)** from routine multi-gradient echo (mGRE) MRI data by matching measured signals to a biophysical dictionary.

---

## Requirements

- MATLAB R2020a or later (tested on R2025b)
- Image Processing Toolbox (required for orientation-informed modes)
- Input data: mGRE magnitude (`iField`), QSM map, brain mask
- Optional: FA map and fiber orientation map (theta) for orientation-informed modes

Example data available at: https://zenodo.org/records/10019720

---

## Project Structure

```
MIMM/
├── Dictionary/                    # Precomputed dictionaries (.mat)
├── Dictionary_Generation/         # Scripts to generate dictionaries
│   ├── generate_dictionary_stochastic.m
│   ├── generate_dictionary_deterministic.m
│   ├── compute_field.m            # 2D myelin field (Wharton & Bowtell 2012)
│   ├── generate_iron_field.m      # Iron field via dipole convolution
│   ├── generate_iron_volume.m     # Random iron placement
│   ├── dipole_kernel.m            # Magnetic dipole kernel (k-space)
│   ├── polar_mesh.m               # Cartesian to polar grid
│   └── bubblebath.m               # 2D random circle packing
├── Matching_Pursuit/              # Core MIMM algorithm
│   ├── MIMM.m                     # Main matching pursuit function
│   ├── interpolate_dictionary.m   # Polynomial interpolation of dictionary TEs
│   └── bounding_box.m             # Brain mask bounding box
├── scripts/                       # Run and visualization scripts
│   ├── run_stochastic_all.m       # Run all modes with stochastic dictionary
│   ├── run_deterministic_MIMM.m   # Run all modes with deterministic dictionary
│   ├── run_basic_MIMM.m           # Run basic MIMM only
│   └── visualize_results.m        # Save output maps as PNG figures
├── Example_Data/                  # Input MRI data (not tracked in git)
├── Example_Results/               # Output results and figures (not tracked in git)
├── MIMM_set_path.m                # Add project folders to MATLAB path
├── RUNME.m                        # Top-level run script
└── NOTES.md                       # Project notes and experiment log
```

---

## How to Run

### 1. Setup
In MATLAB or via bash:
```matlab
cd('/path/to/MIMM')
MIMM_set_path
```

### 2. Dictionary Generation (optional , ~38h)
Precomputed dictionaries are provided. Skip this step unless you need to regenerate them.
```matlab
% See RUNME.m , Dictionary Generation section
```

### 3. Run MIMM

**Via bash (recommended):**
```bash
# All three modes with stochastic dictionary
/path/to/matlab -batch "run('scripts/run_stochastic_all.m')"

# All three modes with deterministic dictionary
/path/to/matlab -batch "run('scripts/run_deterministic_MIMM.m')"
```

**Via MATLAB:**
```matlab
% Load data and dictionary, then call:
MIMM_basic = MIMM(dict, lambda_chi, QSM, Brain_Mask, iField, TE, 'basic');
MIMM_DTI   = MIMM(dict, lambda_chi, QSM, Brain_Mask, iField, TE, 'orientation_informed', FA_DTI, theta_DTI);
MIMM_Atlas = MIMM(dict, lambda_chi, QSM, Brain_Mask, iField, TE, 'orientation_informed', FA_atlas, theta_atlas);
```

### 4. Visualize Results
```bash
/path/to/matlab -batch "run('scripts/visualize_results.m')"
```
Saves PNG figures to `Example_Results/figures/stochastic/` and `Example_Results/figures/deterministic/`.

---

## MIMM Modes

| Mode | Description |
|---|---|
| `basic` | No fiber orientation prior. Faster, less accurate in white matter tracts. |
| `orientation_informed` (DTI) | Uses subject-specific DTI fiber orientation. Most accurate. |
| `orientation_informed` (Atlas) | Uses atlas-based fiber orientation. No extra scan needed. |

---

## Output Maps

Each MIMM run produces a struct with the following fields:

| Field | Description | Unit |
|---|---|---|
| `MVF` | Myelin volume fraction | fraction [0–0.55] |
| `g_ratio` | Inner/outer axon radius ratio | [0.5–1] |
| `FVF` | Fiber volume fraction | fraction [0–0.75] |
| `R2s` | Transverse relaxation rate | s⁻¹ |
| `chi_iron_est` | Iron susceptibility | ppm |
| `chi_myelin` | Myelin susceptibility contribution | ppm |
| `theta_est` | Estimated fiber orientation | degrees |
| `error` | Dictionary matching error | [0–1] |

---

## Notes
- `lambda_chi = 0.015` is the recommended weighting factor between the magnitude and QSM terms (determined by L-curve analysis in the paper)
- The stochastic dictionary produces smoother MVF maps than the deterministic dictionary
- Processing time: ~4 minutes per subject per mode on a modern desktop CPU

---

## Citation
Şişman M, Nguyen TD, Roberts AG, et al. Microstructure-Informed Myelin Mapping (MIMM) from routine multi-echo gradient echo data using multiscale physics modeling of iron and myelin effects and QSM. *Magn Reson Med*. 2025;93:1499–1515. doi:10.1002/mrm.30369
