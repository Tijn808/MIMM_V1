# Figure Reference — MUMC MIMM Pipeline

Figures 01–15 are spatial maps (brain slices). This document covers the analytical figures 16–29.

All figures saved to `OUTPUT_DIR/figures/` as defined in `paths.py`.  
Figures 26–29 require re-running `plot_roi_stats.py` against current data.

---

## Figure 16 — Bland-Altman: Paramagnetic Susceptibility (iron)
**File:** `16_BA_chi_iron.png` | **Script:** `generate_figures.py` | **Data:** voxel-level, WM mask (FA > 0.2)

Bland-Altman comparing MIMM `chi_iron_est` vs chi-sep `chi_pos`, both paramagnetic susceptibility in ppm. Density histogram (hot colourmap) used because of the large number of WM voxels.

**Observed values:**
- Bias = **−0.056 ppm** (MIMM assigns less paramagnetic susceptibility than chi-sep)
- LoA: +1.96 SD = −0.022 ppm, −1.96 SD = −0.090 ppm (range ≈ 0.068 ppm wide)

**Observations:**
The bias is entirely negative — MIMM chi_iron is systematically lower than chi-sep chi_pos across all WM voxels. The density cloud shows a clear **funnel/wedge shape**: the spread of differences grows proportionally with the mean value, indicating proportional bias (the discrepancy is larger where iron is higher). The bulk of the WM density sits at mean 0–0.04 ppm, consistent with WM being iron-poor; the long tail extends to ~0.11 ppm where brainstem and cerebellar tracts contribute.

**Interpretation:**
The −0.056 ppm offset is large relative to typical WM chi_pos values (~0.01–0.08 ppm). MIMM may suppress its iron estimate because the biophysical model distributes susceptibility contributions across iron, myelin, and the water pool simultaneously, whereas chi-sep uses a simpler two-component decomposition. Alternatively, this could reflect a pipeline calibration difference. The proportional bias means the agreement is worst in the most iron-rich voxels, which are also the most clinically relevant (deep grey matter edges, rim lesions).

---

## Figure 17 — Bland-Altman: Diamagnetic Susceptibility (myelin)
**File:** `17_BA_chi_myelin.png` | **Script:** `generate_figures.py` | **Data:** voxel-level, WM mask (FA > 0.2)

Bland-Altman comparing MIMM `|chi_myelin|` vs chi-sep `|chi_neg|`, both in ppm. Absolute values used throughout since both quantities are diamagnetic (stored as negative in physics convention; taking abs makes them directly comparable).

**Observed values:**
- Bias = **−0.058 ppm** (MIMM assigns less diamagnetic susceptibility than chi-sep)
- LoA: +1.96 SD = −0.019 ppm, −1.96 SD = −0.097 ppm (range ≈ 0.078 ppm wide)

**Observations:**
Strikingly similar bias magnitude to iron (−0.058 vs −0.056 ppm). The density peak is around mean 0.04–0.06 ppm, consistent with WM myelin susceptibility. The LoA are slightly wider than for iron, reflecting the greater spatial variability of myelin susceptibility across WM. The funnel shape is less pronounced than in figure 16, with the spread remaining relatively constant across the mean range.

**Interpretation:**
The near-identical bias in both channels (−0.056 and −0.058 ppm) is striking and likely not coincidental. It suggests a **pipeline-level systematic offset** rather than random disagreement: MIMM's biophysical matching may systematically compress susceptibility magnitudes in both channels, or there is a global calibration difference between the chi-sep toolbox and MIMM. The practical implication is that absolute susceptibility values from MIMM should not be directly compared to chi-sep values without a correction factor; relative comparisons across ROIs or individuals remain valid.

---

## Figure 30 — Internal Consistency: MIMM χ_total vs Measured QSM
**File:** `30_chi_total_vs_QSM.png` | **Script:** `generate_figures.py` | **Data:** voxel-level, whole-brain mask

Diagnostic for the susceptibility bias seen in figures 16–17. MIMM does not match its two susceptibility components to QSM independently — it softly matches their **sum** (`χ_total = χ_myelin + χ_iron`) to the measured QSM, via the `chi_error` term in `MIMM.m` weighted by `lambda_chi = 0.015`. This figure tests how well that matching actually holds.

Two panels:
- **Panel A (Bland-Altman):** χ_total − QSM difference vs their mean. Bias quantifies systematic under/over-estimation of total susceptibility.
- **Panel B (density scatter):** χ_total vs QSM with identity and regression lines. The **slope** distinguishes additive bias (slope ≈ 1, line offset) from proportional/compressive bias (slope < 1, MIMM compresses the QSM dynamic range).

**Why this matters:**
The diamagnetic and paramagnetic offsets (figures 16–17, both ≈ −0.057 ppm) arise because MIMM's `χ_myelin = χ_iso × MVF` and `χ_iron = IVF × χ_iron_const` model only the myelin- and iron-specific contributions, using fixed physical constants (χ_iso = −0.1 ppm, χ_iron = 0.3 ppm). Chi-sep, by contrast, captures the *total* diamagnetic and paramagnetic susceptibility of each voxel (myelin + water + proteins; ferritin + all paramagnetic sources). So MIMM components are expected to be smaller in magnitude. This figure checks whether their **sum** still reproduces the measured QSM — the quantity MIMM is actually constrained by.

**What to look for:**
- **Slope near 1, bias near 0:** matching is working; the figures 16–17 offsets are purely a component-attribution difference, not a model failure. Absolute MIMM susceptibility components should then be interpreted only relatively (rankings, not magnitudes).
- **Slope < 1 (compression):** χ_total systematically underestimates high-|QSM| voxels. This would mean lambda_chi = 0.015 is too weak to enforce the QSM constraint — the magnitude-decay matching term `(1 − magnitude_correlation)` dominates and the model drifts from the measured susceptibility. **This is a tunable parameter** — raising lambda_chi would tighten QSM agreement (at the cost of magnitude-curve fit), and would be a legitimate, reportable fix.
- **Low r:** the dictionary's chi_total range does not span the measured QSM range — a dictionary coverage problem rather than a weighting problem.

**Note:** Generated by re-running `generate_figures.py`. The terminal prints `bias`, `slope`, and `r` for direct quoting.

---

## λ_chi Sweep — Tuning the QSM-matching weight (diagnostic, not numbered)
**Files:** `lambda_chi_sweep.png`, `lambda_chi_sweep.csv` (in `mimm_dir`) | **Script:** `MUMC_pipeline/mimm/sweep_lambda_chi.m` (MATLAB)

If figure 30 shows χ_total compressing the QSM range (slope < 1), this sweep finds the `lambda_chi` that fixes it. `lambda_chi` weights the QSM term against the magnitude-decay term in MIMM's cost (`MIMM.m` line 123). The paper's default 0.015 was L-curve-optimised on *their* data; this re-derives it on MUMC data.

**Efficiency:** the expensive dictionary×data correlation and `chi_error` are computed once per slice; only the cheap argmin repeats per `lambda_chi`. A full 9-value sweep ≈ one MIMM run. Basic strategy only (`lambda_theta = 0`) to isolate the effect.

Four panels:
- **Slope vs λ_chi** — should rise toward 1 as λ increases (QSM range reproduced).
- **Bias vs λ_chi** — should move toward 0.
- **Magnitude fit vs λ_chi** — the achieved decay-curve correlation at the chosen dictionary entry. *Drops* as λ rises: this is the cost, because MVF/FVF/g-ratio accuracy depends on the magnitude match, not the QSM match.
- **L-curve** — magnitude fit (y) vs χ MAE (x). The **knee** is the data-driven optimal λ_chi: maximal QSM agreement before the magnitude fit collapses. This is the same L-curve method the paper used to pick 0.015.

**Automated knee selection:** the script detects the L-curve knee (utopia-point method: the sampled λ closest to the ideal corner of low χ-MAE and high magnitude fit, both axis-normalised), marks it with a red star on panel 4, and writes it to `mimm_dir/lambda_chi_recommended.txt` (a one-line data file).

**λ_chi precedence in `run_MIMM_MUMC.m`** (each is a data file, not hard-coded; delete to fall back):
1. **Cohort-locked** — `MUMC_pipeline/lambda_chi_cohort.txt` (shared by all subjects)
2. **Per-subject** — `mimm_dir/lambda_chi_recommended.txt` (this subject's own knee)
3. **Paper default** — 0.015, if neither file exists

The source default is never overwritten; `run_MIMM_MUMC.m` prints which value (and source) it used.

**Recommended cohort workflow (avoids circularity):** tuning λ_chi on the same subject's QSM that MIMM is then constrained by is mildly circular. The defensible approach is *derive once, apply fixed to everyone*:
1. Pick one representative healthy reference subject. Point `paths.m` at it.
2. Run with the lock set: `lock_cohort = true; run('sweep_lambda_chi.m')` — this writes `MUMC_pipeline/lambda_chi_cohort.txt`.
3. For every other subject, just run `run_MIMM_MUMC.m` — it reads the cohort-locked value automatically. Do **not** re-run the sweep per subject.
4. Regenerate figures (`generate_figures.py`, `extract_roi_stats.py`, `plot_roi_stats.py`). Check figure 30 slope moved toward 1 and MVF (figure 18) did not degrade unacceptably.

To unlock and revert the whole cohort to the default, delete `lambda_chi_cohort.txt`.

**If the knee lands near 0.015:** the cleanest conclusion is that the figures 16–17/30 bias is an inherent component-attribution difference (MIMM components are myelin-only / iron-only; chi-sep captures total susceptibility), not a tuning error — report it as such and interpret MIMM susceptibilities relatively.

**CSV columns:** `lambda_chi, bias, slope, r, mean_abs_err, magnitude_fit` — for plotting or quoting exact values.

---

## Figure 18 — ROI Bar Chart: MVF Basic vs Atlas (dual atlas)
**File:** `18_ROI_MVF_dual_atlas.png` | **Script:** `plot_roi_stats.py` | **Data:** 48 JHU DTI-81 ROIs (left) + 20 tractography ROIs (right, if CSV exists)

Two-panel bar chart of mean MVF ± SEM per ROI, basic (red) vs atlas (blue). Left panel: JHU DTI-81 (48 ROIs). Right panel: JHU tractography thr25 (20 ROIs). Merges the previously separate figures 18 and 21 into a single figure that directly invites atlas comparison. Error bars are the standard error of the ROI mean (SD/√n_voxels).

**Observations:**
- **Highest MVF:** Posterior limb of internal capsule (L/R, ~0.38–0.40), splenium of corpus callosum (~0.35), posterior thalamic radiation (~0.30). These are densely myelinated projection fibres.
- **Lowest MVF:** Fornix (cres/column, ~0.05–0.08), tapetum, cingulum hippocampus (~0.10). These are limbic structures with sparser myelination.
- **Basic consistently overestimates atlas** in the high-MVF tracts (red bars extend further right). The gap is especially large for the splenium and posterior thalamic radiation.
- **For low-MVF tracts** (bottom of the chart), basic and atlas are nearly identical, confirming that the orientation prior only matters for tracts with a strong orientation signal relative to B0.
- Error bars are large throughout, indicating substantial within-ROI MVF variability — partly genuine (tract margins, partial volume), partly noise.

**Interpretation:**
The anatomical ranking is biologically plausible and matches the established literature on WM myelination. The selective overestimation in projection and callosal fibres, but not limbic fibres, is exactly what the physics predicts: those tracts run perpendicular to B0, creating the largest orientation-dependent signal confound.

---

## Figure 19 — ROI Parameter Heatmap
**File:** `19_ROI_heatmap.png` | **Script:** `plot_roi_stats.py` | **Data:** 48 JHU DTI-81 ROIs × 7 parameters (MVF basic, MVF atlas, FVF, g-ratio, R2*, |χ⁻|, χ⁺), each column normalised 0→1

**Observations:**
- **MVF basic and MVF atlas track extremely closely** — nearly identical colour columns, confirming the atlas prior scales rather than reorders the myelin ranking.
- **FVF column closely parallels MVF** — tracts with high myelin have proportionally high fibre volume fraction, consistent with the biophysical model (more axons = more myelin sheath).
- **g-ratio column shows surprisingly little variation** — most tracts land in a similar mid-range (pale, ~0.6–0.8), suggesting MIMM produces a relatively uniform g-ratio across WM. A few limbic tracts have notably high g-ratio (thin myelin relative to axon calibre).
- **R2* column shows a distinct pattern from MVF** — some brainstem and cerebellar tracts (middle cerebellar peduncle, cerebellar peduncles) are bright in R2* but not in MVF, consistent with iron deposition in these regions without corresponding high myelination.
- **|χ⁻| column broadly tracks MVF** but with exceptions: some tracts with moderate MVF have relatively high |χ⁻| from chi-sep.
- **χ⁺ column** highlights different tracts than MVF — some corticospinal and brainstem tracts show relatively higher paramagnetic signal, consistent with iron content independent of myelin.
- **Fornix, tapetum, and limbic tracts** are consistently dark across all myelin-related columns, confirming they are the least myelinated structures.

**Interpretation:**
The heatmap makes the multi-parametric story visible at once: myelin (MVF, FVF, |χ⁻|) and iron (R2*, χ⁺) show partially independent spatial distributions, validating that MIMM and chi-sep are capturing different tissue properties rather than the same signal in different units.

---

## Figure 20 — ROI Scatter: MVF Basic vs Atlas
**File:** `20_ROI_MVF_scatter.png` | **Script:** `plot_roi_stats.py` | **Data:** 48 JHU DTI-81 ROIs

Scatter of MVF basic (x) vs MVF atlas (y) per ROI, with identity line and ±1 SEM error bars (SD/√n_voxels). SEM replaced the original SD bars, which spanned nearly the whole plot and obscured the relationship. Highlighted-tract labels use leader lines to avoid collisions in the dense 0.18–0.25 cluster.

**Observations:**
- **Most points cluster around the identity line**, confirming overall good agreement between basic and atlas MIMM.
- **High-MVF tracts fall below the identity** (basic > atlas): splenium (~0.35 basic vs ~0.26 atlas), posterior thalamic radiation, body of corpus callosum. Basic overestimates here.
- **Posterior limb of internal capsule sits very close to the identity** (~0.40 on both axes) — despite being a high-MVF tract, the atlas prior does not strongly reduce its estimate. This is because PLIC runs nearly parallel to B0, so the orientation-dependent confound is small.
- **Low-MVF tracts (bottom-left) sit precisely on the identity** — basic and atlas agree exactly for limbic and low-FA tracts.
- **Error bars are large**, particularly in the x-direction (basic), reflecting higher within-ROI variance in the no-prior estimate.

**Interpretation:**
The identity-line deviation pattern directly validates the physics hypothesis: only tracts that are both highly myelinated AND oriented perpendicular to B0 are significantly overestimated by basic MIMM. PLIC being on the identity despite high MVF is a key finding — it shows the atlas prior is not globally suppressing MVF, but selectively correcting orientation-dependent bias.

---

## Figure 21 — Merged into Figure 18

The tractography bar chart (previously figure 21) is now the right panel of figure 18 (`18_ROI_MVF_dual_atlas.png`). Key observations from when it existed standalone: forceps major highest MVF (~0.38 basic, ~0.28 atlas, +37% overestimation); basic overestimates atlas consistently across all 20 tracts; results mirror the DTI-81 splenium finding.

---

## Figure 22 — Atlas Comparison: DTI-81 vs Tractography (CV and r)
**File:** `22_atlas_comparison.png` | **Script:** `plot_roi_stats.py` | **Data:** 48 DTI-81 ROIs vs 20 tractography tracts

Two panels comparing the two atlases on ROI homogeneity (CV) and MIMM–chi-sep spatial agreement (Pearson r).

**Observed values:**
- **CV:** DTI-81 μ = 0.57, Tractography μ = 0.57 — identical
- **Pearson r (MVF vs |χ⁻|):** DTI-81 μ = 0.75, Tractography μ = 0.77

**Observations:**
- **Left panel (CV):** The two distributions are nearly identical in mean (both 0.57). DTI-81 has a right tail with some very heterogeneous ROIs (CV up to 1.3), while tractography terminates around 0.9. Tractography ROIs are not more homogeneous despite covering fewer, larger tracts.
- **Right panel (r):** Both atlases show a broad distribution of within-ROI MIMM–chi-sep correlation, peaking around r = 0.70–0.80. Tractography is marginally better (0.77 vs 0.75) but the difference is negligible.

**Interpretation:**
There is no meaningful advantage to either atlas in terms of ROI homogeneity or MIMM–chi-sep spatial agreement. The DTI-81 atlas is preferred for its finer anatomical parcellation (48 vs 20 ROIs) which allows tract-specific conclusions. The within-ROI Pearson r of 0.75–0.77 is strong and consistent, meaning MVF and |χ⁻| co-vary spatially within tracts — both methods detect the same intra-ROI spatial patterns of myelination.

---

## Figure 23 — MVF Overestimation Ranked by ROI
**File:** `23_overestimation_ranked.png` | **Script:** `analyse_overestimation.py` | **Data:** 48 JHU DTI-81 ROIs

Ranked bar chart of (MVF basic − MVF atlas) per ROI. Red = overestimation, blue = underestimation.

**Observed values (top overestimated):**
- Splenium of corpus callosum: **+37%** (absolute: ~+0.097)
- Posterior thalamic radiation L/R: **+21–22%**
- Retrolenticular internal capsule R: **+22%**
- Body of corpus callosum: ~+18%
- Genu of corpus callosum: ~+17%

**Observed underestimated (blue bars):**
- Medial lemniscus, posterior limb of internal capsule, superior corona radiata, corticospinal tract — these run predominantly parallel to B0.

**Observations:**
The splenium is a clear outlier — its 37% overestimation is nearly double the next highest. This makes physical sense: the splenium runs almost perfectly perpendicular to B0 (θ ≈ 85°), maximising the orientation-dependent signal. The entire corpus callosum body-to-splenium axis shows large overestimation. Projection fibres (thalamic radiations, retrolenticular IC) are also strongly overestimated. The few blue bars are all fibres running roughly parallel to B0 where the atlas slightly raises MVF through regularisation.

**Interpretation:**
This figure is the clearest single demonstration that orientation prior matters. A clinician using basic MIMM on the splenium would overestimate MVF by 37% — a clinically meaningful error for MS lesion characterisation or demyelination tracking. The atlas correction is therefore not optional for callosal and thalamic tracts.

---

## Figure 24 — Overestimation vs Fibre Angle θ
**File:** `24_overestimation_vs_theta.png` | **Script:** `analyse_overestimation.py` | **Data:** 48 JHU DTI-81 ROIs

Scatter of mean ROI fibre angle θ vs overestimation (MVF basic − MVF atlas), colour = MVF atlas (reference myelin content).

**Observed values:**
- **r = 0.66, p = 0.000**

**Observations:**
- Strong, significant positive correlation: higher fibre angle → more overestimation by basic MIMM.
- **Splenium** is the clear top-right outlier (θ ≈ 85°, overestimation ≈ +0.097) and is bright yellow (high MVF atlas ~0.38) — the combination of high myelin AND near-perpendicular orientation maximises the error.
- **Posterior thalamic radiation** (θ ≈ 70°, overestimation +0.07) is the second-highest, also yellow-orange (high MVF).
- **Body of corpus callosum** (θ ≈ 65°, overestimation ~0.05) follows the trend.
- **CST and PLIC** (θ ≈ 25–30°) cluster near zero overestimation, consistent with their more parallel orientation.
- Some scatter around the regression line in the 50–80° range — not all high-angle tracts are equally overestimated, suggesting myelin content also modulates the error.

**Interpretation:**
r = 0.66 confirms the core physics prediction that orientation-dependent signal is the primary driver of MVF overestimation. The residual scatter is partly explained by MVF magnitude (high-MVF tracts show more absolute overestimation for the same θ, visible in the colour coding). This figure is the key physics validation result of the project.

---

## Figure 25 — Overestimation vs FA
**File:** `25_overestimation_vs_FA.png` | **Script:** `analyse_overestimation.py` | **Data:** 48 JHU DTI-81 ROIs

Scatter of mean FA vs overestimation, colour = mean fibre angle θ.

**Observed values:**
- **r = 0.58, p = 0.000**

**Observations:**
- Positive correlation between FA and overestimation (r = 0.58), weaker than the θ correlation (0.66).
- The colour coding reveals the mechanism: **all highly overestimated points (top half) are orange-red = high θ**, and **all underestimated or near-zero points are blue = low θ**. FA is a confound of θ, not an independent predictor.
- **Splenium** (top right, FA ≈ 0.34, overestimation ≈ 0.097) is both the highest-FA and highest-θ callosal segment.
- **Two deep-blue points (θ ≈ 5–10°) near zero overestimation** despite moderate FA (~0.28–0.30) — these are PLIC and CST, parallel to B0, confirming that high FA alone does not cause overestimation without the orientation component.
- Some moderate-FA tracts (0.22–0.26) span a wide range of overestimation — entirely explained by their θ colour.

**Interpretation:**
FA correlates with overestimation only because high-FA tracts tend to have well-defined orientations, making them more susceptible to orientation-dependent signal when also perpendicular to B0. FA is not the mechanistic driver — θ is. The r = 0.58 for FA vs r = 0.66 for θ confirms θ is the stronger predictor, as expected from the physics.

---

## Figure 45 — Spatial Overestimation with JHU Atlas Overlay
**File:** `45_overestimation_spatial_JHU.png` | **Script:** `analyse_overestimation.py` | **Data:** voxel-level + ROI-level, whole brain

Two-row × three-column figure bridging the voxel-level and ROI-level overestimation views.

- **Row 1 (voxel-wise):** MVF basic − MVF atlas shown continuously per voxel on axial/coronal/sagittal slices. Red = basic overestimates, blue = basic underestimates. JHU ROI boundaries overlaid as thin white contour lines. Shows the real spatial extent and magnitude of the overestimation.
- **Row 2 (ROI-mean):** Each voxel coloured by its ROI's mean overestimation value (from `roi_stats.csv`). Removes within-ROI noise and shows the structured, atlas-level result in anatomical context. Top-5 most overestimated ROIs labelled by name and percentage on the axial slice.

Both rows share the same diverging colormap (RdBu_r) and ±99th-percentile scale. Anatomy (echo 1 magnitude) in greyscale underneath.

**What to look for:** The corpus callosum should be prominently red (especially splenium and body), thalamic radiations orange-red, while tracts running parallel to B0 (CST, PLIC) should be near zero or faintly blue. The JHU boundaries show which ROI each coloured region belongs to.

---

## Figures 26–29 — Pending (re-run `plot_roi_stats.py`)

These figures are generated by the updated `plot_roi_stats.py` and require the existing `roi_stats.csv`. They have not yet been generated and will appear in `OUTPUT_DIR/figures/` after the next run.

---

## Figure 26 — ROI Scatter: MIMM MVF vs chi-sep |χ⁻|
**File:** `26_MVF_vs_chineg.png` | **Script:** `plot_roi_stats.py` | **Data:** 48 JHU DTI-81 ROIs, colour by tract family

Scatter of MIMM MVF basic (fraction, y) vs chi-sep |χ⁻| (ppm, x) per ROI, ±1 SEM error bars. Both are myelin markers but different physical quantities — no identity line. Motivated by Şişman et al. 2024 (ISMRM), which showed |χ⁻| correlates with MBP optical density at r = 0.90 (histopathological validation), making χ⁻ a well-validated myelin reference.

**Expected observations:**
Strong positive correlation — both methods should rank WM tracts consistently by myelination. Corpus callosum ROIs (gold dots) expected in the upper-right (high on both axes). Limbic tracts (pink) expected lower-left. Deviations from the regression by tract family reveal where the two methods disagree on relative myelination.

---

## Figure 27 — ROI Scatter: MIMM chi_iron vs chi-sep χ⁺
**File:** `27_iron_chi_comparison.png` | **Script:** `plot_roi_stats.py` | **Data:** 48 JHU DTI-81 ROIs, colour by tract family

Scatter of MIMM `chi_iron_est` (ppm, y) vs chi-sep `chi_pos` (ppm, x) per ROI, ±1 SEM error bars. Both are paramagnetic susceptibility in the same units. ROI-level companion to the voxel Bland-Altman in figure 16. Motivated by Şişman 2024: χ⁺ outperforms QSM for iron (r = 0.67 vs 0.55 with microglia counts), making it a stronger iron reference than plain QSM.

**Expected observations:**
Given the systematic −0.056 ppm bias seen in figure 16, the regression line should sit below an identity line. If the bias is purely additive, the slope should be close to 1 with a negative intercept; proportional bias would give slope > 1.

---

## Figure 28 — ROI Scatter: |χ myelin| MIMM vs |χ⁻| chi-sep (same units, ppm)
**File:** `28_chi_myelin_vs_chineg_ROI.png` | **Script:** `plot_roi_stats.py` | **Data:** 48 JHU DTI-81 ROIs, colour by tract family

Scatter of MIMM `|chi_myelin|` (ppm, x) vs chi-sep `|χ⁻|` (ppm, y) per ROI, ±1 SEM error bars, with identity line. Unlike figure 26 (MVF vs |χ⁻|), both axes are the same physical quantity. The identity line shows absolute numerical agreement — points below it mean chi-sep reports larger diamagnetic susceptibility than MIMM.

**Expected observations:**
Given the voxel-level bias of −0.058 ppm (figure 17), all ROI points should sit above the identity (chi-sep > MIMM). The regression slope relative to 1.0 will indicate whether the discrepancy is additive (parallel to identity) or proportional (slope ≠ 1). Tract family colouring may reveal whether the offset is uniform or structure-dependent.

---

## Figures 49–51 — Three-Way Myelin Comparison: chi-sep | MIMM atlas | T2-GRASE MWF
**Files:** `49_three_way_spatial.png`, `50_three_way_scatter.png`, `51_three_way_ranking.png`
**Script:** `plot_three_way.py` | **Data:** voxel + 48 JHU DTI-81 ROIs | **Dormant until `grase/MWF.nii.gz` exists**

The headline cross-method figure set: the three independent in-vivo myelin
markers side by side — MIMM atlas MVF, chi-separation |χ⁻| (diamagnetic
susceptibility), and T2-GRASE MWF (the closest thing to a myelin gold standard).
Uses MIMM **atlas** (orientation-corrected), not basic.

- **Fig 49 (spatial):** the three maps on axial/coronal/sagittal slices, each on
  its own physical scale with its own colourbar (MVF & MWF as fractions, |χ⁻| in
  ppm). Visual co-localisation check.
- **Fig 50 (pairwise scatter):** three ROI scatters with MWF on the x-axis as the
  reference — MIMM vs MWF (identity line, both fractions), chi-sep vs MWF, and
  MIMM vs chi-sep. Each prints Pearson r and p. This is where the three-way
  agreement is quantified.
- **Fig 51 (ranking):** per-ROI grouped bars, each marker min–max normalised to
  0–1 (so the different units share one axis), ROIs sorted by MWF. Aligned bars =
  the three markers rank that tract consistently; divergent bars flag tracts where
  the methods disagree on relative myelination.

**Expected (once T2 data lands):** strong positive correlations in all three
panels of fig 50. MIMM-vs-MWF and chi-vs-MWF are the genuine external validations
(MWF independent of mGRE); MIMM-vs-chi reproduces the figure-26 result with the
atlas estimate. Tracts where chi-sep diverges from MWF but MIMM tracks it (or vice
versa) are the scientifically interesting cases.

---

## Figure 29 — ROI Scatter: MIMM MVF vs FA atlas
**File:** `29_MVF_vs_FA.png` | **Script:** `plot_roi_stats.py` | **Data:** 48 JHU DTI-81 ROIs, colour by tract family

Scatter of MIMM MVF basic (fraction, y) vs mean FA from HCP1065 atlas (x) per ROI, ±1 SEM error bars on MVF. FA comes from DTI — a completely independent acquisition with no shared signal with mGRE. A positive correlation is expected since both FA and MVF track WM structural integrity. MVF basic is used (not atlas) to keep the comparison independent: MVF atlas incorporates the same FA atlas for orientation weighting, which would artificially inflate the correlation.

**Expected observations:**
Positive correlation with r > 0.5. Corpus callosum and PLIC expected upper-right (high FA, high MVF). Limbic tracts (cingulum hippocampus, fornix) expected lower-left. Outliers from the regression are scientifically interesting — brainstem tracts with moderate FA but variable MVF, or tracts where MIMM and DTI disagree on the myelination ranking.
