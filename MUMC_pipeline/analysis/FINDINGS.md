# Findings — MIMM vs Chi-Separation on the Test Subject

A consolidated write-up of what the figures actually show. All numbers below are
measured on the single test subject (`ME_GRE`), 48 JHU DTI-81 white-matter ROIs
(50 rows incl. two background labels) plus voxel-level analyses on the WM mask.
Companion to `FIGURES.md` (which documents each figure individually) and the
runbook. Source data: `ME_GRE/analysis/roi_stats.csv`.

---

## 1. Headline results

| # | Comparison | Statistic | Figure |
|---|------------|-----------|--------|
| 1 | MIMM MVF vs chi-sep \|χ⁻\| (myelin) | **r = 0.93**, p ≈ 5×10⁻²² | 26 |
| 2 | MIMM χ_iron vs chi-sep χ⁺ (iron) | **r = 0.86**, p ≈ 3×10⁻¹⁵ | 27 |
| 3 | Overestimation vs fibre angle θ | **r = 0.66**, p ≈ 2×10⁻⁷ | 24 |
| 4 | Overestimation vs FA | r = 0.58, p ≈ 1×10⁻⁵ | 25 |
| 5 | MIMM MVF vs DTI FA (independent) | r = 0.51, p ≈ 1.5×10⁻⁴ | 29 |
| 6 | χ_total vs measured QSM (slope) | slope = 0.22 at λ=0.015 | 30 |

Four scientific stories come out of this: **(A)** MIMM's myelin map validates
strongly against an independent, histologically-grounded myelin marker; **(B)**
MIMM correctly separates iron from myelin; **(C)** the orientation prior corrects
a real, physics-predicted overestimation; **(D)** an honest limitation remains in
the absolute susceptibility channel.

---

## A. Myelin validation — MIMM agrees with chi-separation (Fig 26)

**MIMM MVF vs chi-sep |χ⁻| across 48 WM tracts: r = 0.93 (p < 0.001).**

Two independent myelin estimates — one biophysical-dictionary (MIMM), one
susceptibility-source-separation (chi-sep) — rank white-matter tracts almost
identically by myelination. χ⁻ is a meaningful reference because Şişman et al.
(2024, ISMRM) showed diamagnetic susceptibility correlates with myelin basic
protein optical density at r ≈ 0.90 histopathologically. So a cross-method
r = 0.93 at the tract level is a genuine validation of the MIMM myelin map.

The anatomical ranking is biologically correct (Figs 18–19): highest MVF in the
posterior limb of the internal capsule (~0.40), splenium (~0.35), and posterior
thalamic radiation; lowest in fornix, tapetum, and limbic tracts (~0.05–0.10).


Note: |χ_myelin| from MIMM vs |χ⁻| (Fig 28) gives the *same* r = 0.93, as
expected — MIMM defines χ_myelin = χ_iso·MVF, so it is collinear with MVF. Fig 28
adds the absolute-scale comparison (same units), not new ranking information.

## B. Iron/myelin separation works — and it explains the Fig 26 outliers (Fig 27)

**MIMM χ_iron vs chi-sep χ⁺: r = 0.86 (p < 0.001).**

The paramagnetic (iron) channels also agree across tracts. More importantly, this
figure *resolves the outliers in Fig 26*: the two brainstem/corticospinal tracts
that sit **below** the myelin regression in Fig 26 (high |χ⁻| but only moderate
MVF) are the **same** tracts that sit high in **both** iron estimates here. The
"extra" diamagnetic-looking susceptibility in those tracts is not myelin — it is
iron, and MIMM routes it to the iron channel rather than inflating MVF. That
cross-figure consistency is positive evidence that MIMM separates the two
sources the way it is designed to, rather than an artefact.

## C. The orientation prior corrects a real, physics-predicted bias (Figs 23–25, 45)

Basic MIMM (no orientation prior) systematically **overestimates** MVF in tracts
oriented perpendicular to B0:

- **Overestimation scales with fibre angle θ: r = 0.66 (p < 0.001)** (Fig 24) —
  the core physics result. Higher θ (more perpendicular to B0) → larger MVF
  overestimation by basic MIMM.
- **Splenium of the corpus callosum is the worst case: +37%** (θ ≈ 80°), nearly
  double the next tract (Fig 23). Ranked overestimation:

  | Tract | Overestimation | θ |
  |-------|----------------|---|
  | Splenium of corpus callosum | **+37%** | 80° |
  | Posterior thalamic radiation L | +28% | 72° |
  | Body of corpus callosum | +27% | 69° |
  | Posterior thalamic radiation R | +21% | 76° |
  | Retrolenticular internal capsule R | +19% | 63° |

- **FA vs overestimation: r = 0.58** (Fig 25), but the colour coding shows FA is a
  *confound* of θ, not an independent driver — every strongly-overestimated tract
  is also high-θ. θ is the mechanism.
- Tracts running parallel to B0 (CST, PLIC) show near-zero overestimation despite
  high MVF, confirming the prior is *selective*, not a global rescaling (Fig 20:
  PLIC sits on the identity line; Fig 45: corpus callosum lit red, CST/PLIC neutral).

Clinical implication: using basic MIMM on the splenium overestimates myelin by
~37% — a clinically meaningful error for MS demyelination tracking. The atlas
orientation correction is not optional for callosal and thalamic tracts.

## D. Honest limitation — absolute susceptibility (Figs 16–17, 30)

MIMM's absolute susceptibility magnitudes do **not** match chi-sep one-for-one:

- Voxel-level Bland-Altman bias: **−0.056 ppm (iron)**, **−0.058 ppm (myelin)**
  (Figs 16–17). The near-identical offset in both channels points to a systematic
  attribution difference, not random noise.
- This is expected by construction: MIMM models only the myelin-specific
  (χ_iso·MVF) and iron-specific (IVF·χ_iron) contributions with fixed constants,
  whereas chi-sep captures the *total* diamagnetic/paramagnetic susceptibility of
  each voxel (incl. water, proteins, all paramagnetic sources). MIMM components
  are therefore expected to be smaller in magnitude.
- The quantity MIMM is actually constrained by is the **sum** χ_total matched to
  measured QSM (the `chi_error` term, weight λ_chi = 0.015). Fig 30 tests this:
  **χ_total vs QSM slope = 0.22** — MIMM strongly compresses the QSM dynamic range.

### λ_chi sweep — the compression is a dictionary-coverage ceiling, not just tuning

Raising λ_chi (which weights QSM-matching against the magnitude-decay fit) does
improve the slope, but it **saturates** well below 1:

| λ_chi | slope | r (χ_total vs QSM) | magnitude fit |
|-------|-------|--------------------|---------------|
| 0.015 (paper default) | 0.22 | 0.53 | 0.99851 |
| 0.025 (auto-knee) | 0.25 | 0.58 | 0.99848 |
| 0.05 | 0.29 | 0.63 | 0.99840 |
| 0.10 | 0.33 | 0.67 | 0.99826 |
| 0.40 | 0.38 | 0.72 | 0.99779 |

Even at λ = 0.40 (27× the default) the slope only reaches 0.38 while the
magnitude fit barely moves (0.9985 → 0.9978). The QSM agreement is bounded by how
far the dictionary's χ_total range extends, not by the weighting — i.e. a
**dictionary-coverage limitation**. The L-curve auto-knee lands at λ = 0.025,
which is defensibly conservative.

**Conclusion for the report:** MIMM susceptibility components should be
interpreted *relatively* (rankings, cross-ROI and cross-subject comparisons),
not as absolute ppm values directly comparable to chi-sep. Tuning λ_chi cannot
close the gap; the relative agreement (Figs 26–27, r = 0.86–0.93) is what holds.

---

## 2. Secondary observations

- **Atlas choice doesn't matter much (Fig 22):** DTI-81 (48 ROIs) and JHU
  tractography (20 tracts) give identical mean ROI homogeneity (CV = 0.57 both)
  and near-identical MVF–χ⁻ within-ROI correlation (r = 0.75 vs 0.77). DTI-81 is
  preferred only for its finer parcellation.
- **Independent DTI cross-check (Fig 29):** MIMM MVF vs DTI FA gives r = 0.51 —
  a positive correlation across a *completely independent acquisition* (no shared
  signal with mGRE). Basic MVF is used deliberately so the FA atlas doesn't inflate
  the correlation. Confirms MIMM tracks genuine WM structure, not mGRE artefacts.
- **Multiparametric independence (Fig 19):** myelin markers (MVF, FVF, |χ⁻|) and
  iron markers (R2*, χ⁺) show partially independent spatial distributions —
  MIMM and chi-sep capture different tissue properties, not one signal in two units.
- **g-ratio is nearly uniform across WM (Fig 19):** most tracts ~0.6–0.8, a few
  limbic tracts higher (thin myelin relative to axon calibre).

---

## 3. What's still pending (needs patient/cohort data)

These figures are built and validated but dormant until real MUMC data is
available (blocked by data governance, not by the pipeline):

- **MWF validation (Figs 38–41):** MIMM MVF vs T2-GRASE MWF — the strongest
  possible external myelin reference. Needs the GRASE/MWF output.
- **Three-way comparison (Figs 49–51):** chi-sep |χ⁻|, MIMM atlas MVF and T2-GRASE
  MWF side by side — spatial maps, pairwise scatters (MWF as reference), and a
  normalised tract-ranking chart. The figure that ties the whole validation story
  together once T2 is available. Built and tested on synthetic MWF; needs GRASE.
- **Lesion analysis (Figs 46–48):** lesion-vs-NAWM contrast and per-lesion
  MVF-vs-MWF. Needs MS patients with lesion masks.
- **Cohort figures (Figs 31–37, 42–44):** between-subject means with proper SEM.
  Needs N > 1 subjects.

---

*Generated for the internship report. All test-subject statistics reproducible
from `ME_GRE/analysis/roi_stats.csv` and `ME_GRE/mimm/lambda_chi_sweep.csv`.*
