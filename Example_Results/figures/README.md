# Figure Interpretation — MIMM Sample Data

All figures show axial slice 80 of the example subject. Each figure contains three panels side by side: **Basic**, **DTI orientation-informed**, and **Atlas orientation-informed**. Figures are saved for both the **stochastic** (20,000 entry) and **deterministic** (12,540 entry) dictionaries.

---

## MVF — Myelin Volume Fraction

**Range displayed:** 0 – 0.4 (fraction)

The MVF map shows the estimated fraction of each voxel occupied by myelin. Bright regions contain more myelin.

**What we see:**
- Clear white matter structures are visible: corpus callosum (the bright arch-shaped structure running left–right), internal capsule, corona radiata, and posterior white matter tracts.
- Gray matter cortex appears uniformly darker, consistent with its lower myelin content.
- Ventricles (CSF-filled cavities in the centre) are black, as expected — no myelin.
- All three modes (Basic, DTI, Atlas) produce very similar MVF maps. DTI is marginally sharper at white matter boundaries.
- Stochastic and deterministic dictionaries give nearly identical MVF results for this subject.

**Expected biology:** White matter MVF is typically 0.10–0.25 in major tracts; gray matter ~0.03–0.10. The observed values are consistent with this.

---

## g_ratio — Inner/Outer Axon Radius Ratio

**Range displayed:** 0.5 – 1.0

The g_ratio describes how thick the myelin sheath is relative to the total fibre radius. A lower g_ratio means a proportionally thicker sheath; g_ratio = 1 means no myelin.

**What we see:**
- The map is generally bright (g_ratio 0.8–1.0), meaning most voxels have a relatively thin myelin sheath, which is physiologically normal.
- Darker patches in deep white matter tracts correspond to regions with thicker myelin sheaths and higher MVF.
- The map is noisy throughout. This is expected: g_ratio is a derived quantity that depends on the ratio of MVF to FVF and is sensitive to dictionary sampling. The stochastic dictionary is slightly smoother; the deterministic dictionary is patchier.
- Gray matter regions correctly approach g_ratio ≈ 1 (effectively unmyelinated).

---

## FVF — Fiber Volume Fraction

**Range displayed:** 0 – 0.8 (fraction)

FVF is the total fraction of the voxel occupied by axon cylinders plus their myelin sheaths combined (i.e. the myelinated fibre cross-section).

**What we see:**
- Most brain parenchyma is very bright (FVF ≈ 0.5–0.7), indicating the model fills a large proportion of each voxel with fibre structure.
- Ventricles are dark (no fibres in CSF).
- There is limited spatial contrast between white and gray matter at this slice. FVF distinguishes tissue from CSF rather than tract-specific anatomy at this scale.
- High FVF values throughout may partly reflect that the dictionary matches by also fitting gray matter voxels with high FVF and low MVF, since it has no explicit GM/WM class constraint.

---

## R2s — Transverse Relaxation Rate

**Range displayed:** 0 – 100 s⁻¹

R2* (R2s) measures the rate of signal decay in a gradient echo sequence. Higher R2* means faster signal decay, associated with iron content, myelin, and susceptibility sources.

**What we see:**
- The map appears dark overall because most brain tissue has R2* in the range 15–25 s⁻¹, while the display scale extends to 100 s⁻¹.
- Very bright spots are visible at venous structures (large draining veins), which have extremely high R2* due to deoxyhemoglobin.
- White matter regions are slightly brighter than gray matter cortex, consistent with the known higher R2* of myelinated tissue.
- All three modes produce almost identical R2* maps, as this parameter is largely determined by the mGRE signal decay curve independent of orientation.

---

## chi_iron_est — Iron Susceptibility Contribution

**Range displayed:** 0 – 0.3 ppm

This map isolates the paramagnetic susceptibility contribution attributed to iron (ferritin/hemosiderin), separated from the myelin diamagnetic contribution.

**What we see:**
- Most of the brain is near zero (dark), consistent with low iron content in white matter and cortex.
- Moderately bright regions are visible around the central brain, corresponding to the basal ganglia area where iron naturally accumulates (globus pallidus, putamen).
- Bright spots also appear at vascular structures.
- DTI mode shows slightly more focal iron contrast compared to Basic, suggesting the orientation correction affects how iron and myelin susceptibility are partitioned.

---

## chi_myelin — Myelin Susceptibility Contribution

**Range displayed:** −0.06 – 0 ppm

Myelin has a negative (diamagnetic) magnetic susceptibility. This map shows the susceptibility contribution from myelin alone, after separating the iron component.

**What we see:**
- The background (outside the brain and in CSF/ventricles) is white (near zero), and myelinated regions are grey-to-dark.
- This is anatomically the clearest map of white matter microstructure: the corpus callosum, internal capsule, and posterior white matter tracts are all sharply delineated as dark structures against the lighter background.
- All three modes show very similar patterns, suggesting chi_myelin is robustly estimated regardless of orientation strategy.
- The deterministic dictionary produces slightly sharper spatial contrast, while the stochastic dictionary appears smoother.

---

## theta_est — Estimated Fiber Orientation Angle

**Range displayed:** 0 – 90°

theta_est is the angle between the dominant fibre direction and the main magnetic field (B0). It affects the anisotropic susceptibility and relaxation of white matter.

**What we see:**
- **This is the map where the three modes differ most dramatically.**
- **Basic mode:** Very noisy — the estimated orientation is essentially random voxel-by-voxel, because no orientation prior is provided. The algorithm selects whatever angle minimises the dictionary error locally, without spatial coherence.
- **DTI mode:** Smooth, spatially coherent orientation map. Clear fibre tract anatomy is visible. The corpus callosum shows low angles (fibres near-perpendicular to B0), while other tracts show higher angles. The structured appearance confirms that the DTI orientation prior is guiding the matching effectively.
- **Atlas mode:** Very similar to DTI, with slightly less fine-grained detail. The atlas-derived orientation is a good approximation of the subject-specific DTI orientation and produces near-identical results in major tracts.
- This map most clearly demonstrates the benefit of orientation-informed MIMM over basic MIMM.

---

## error — Dictionary Matching Error

**Range displayed:** 0 – 0.25

The error is the normalised residual between the measured mGRE signal and the best-matching dictionary entry. Lower is better.

**What we see:**
- The map is almost completely black (error ≈ 0–0.01 everywhere within the brain mask), indicating an excellent dictionary fit throughout.
- Very faint bright spots are visible at the edges of ventricles and at tissue–CSF interfaces, where the biophysical model is less well suited (partial volume effects, CSF pulsation).
- Large vessel locations may also show slightly elevated error due to flow effects.
- No systematic spatial pattern of high error, which confirms that the dictionary adequately samples the range of biophysical parameters present in this subject.

---

## Stochastic vs. Deterministic Dictionary

| Map | Stochastic | Deterministic |
|---|---|---|
| MVF | Smooth, similar values | Slightly patchier, similar values |
| g_ratio | Smoother | Noisier/patchier |
| FVF | Similar | Similar |
| R2s | Nearly identical | Nearly identical |
| chi_iron_est | Similar | Similar |
| chi_myelin | Smooth boundaries | Slightly sharper detail |
| theta_est | Basic is noisy; DTI/Atlas coherent | Basic noisier; DTI/Atlas coherent |
| error | Near zero throughout | Near zero throughout |

The stochastic dictionary (random sampling over parameter space) produces smoother maps because nearby dictionary entries are statistically spread. The deterministic dictionary (systematic grid) has fewer entries in some regions of parameter space, producing slightly more discrete/patchy results. Both are valid — the stochastic dictionary is generally preferred for final results.

---

## Data source

Figures generated from the Şişman et al. 2025 example dataset (Zenodo record 10019720), using the MIMM implementation with `lambda_chi = 0.015`.
