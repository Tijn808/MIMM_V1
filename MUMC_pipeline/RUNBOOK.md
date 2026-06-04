# MUMC MIMM Pipeline — Runbook

Step-by-step for processing MUMC subjects from raw DICOM to cohort figures.
Everything here was validated on the test subject (`ME_GRE`) except the
data-dependent items flagged ⚠️.

Tool locations (this machine; not on PATH — see memory `reference_toolpaths`):
- MATLAB R2025b: `/home/tijn-saes/Desktop/bin/matlab`
- FSL: `FSLDIR=/home/tijn-saes/fsl` (flirt, dcm2niix, …)
- SPM12 + LST: `~/spm12`, `~/spm12/toolbox/LST`
- Chi-sep toolbox: `~/Documents/Internship/Chisep_Toolbox_v1.2.1`

```bash
export FSLDIR=/home/tijn-saes/fsl
export PATH="$FSLDIR/bin:$PATH"
MATLAB=/home/tijn-saes/Desktop/bin/matlab
REPO=~/Documents/Internship/MIMM
```

---

## 0. Cohort layout

```
COHORT_DIR/
  sub-01/  sub-02/  ...        ← one folder per subject
  cohort/                      ← created automatically by aggregation
```

Each `sub-XX/` is built by `sort_dicom.sh` from that subject's DICOM session.

---

## 1. DICOM → NIfTI (per subject)

```bash
bash $REPO/MUMC_pipeline/preprocess/sort_dicom.sh \
  "<DICOM session dir>" \
  COHORT_DIR/sub-01
```

Produces in `sub-01/`:
- `501-ME_GRE_e{1..5}.nii.gz` + `_ph` (+ JSON) — MIMM input
- `grase/grase.nii.gz` — raw GRASE
- `lesion/FLAIR_native.nii.gz` — FLAIR
- `dti/dti.nii.gz` + `.bval`/`.bvec`, `dti/reverse_b0.nii.gz`
- `t1w/T1w_native.nii.gz`

**Verify:** echo times in `501-ME_GRE_e1.json` are 6/12/18/24/30 ms; 5 mag + 5 phase files.

---

## 2. ⚠️ One-time config from the data

**2a. DTI topup params** — read `sub-01/dti/dti.json`:
```bash
grep -E "TotalReadoutTime|PhaseEncodingDirection" COHORT_DIR/sub-01/dti/dti.json
```
Set `PE_DIR` and `READOUT_TIME` in `preprocess/preprocess_dti.sh` accordingly.

**2b. MWF source** — does MUMC's `derivatives/` already have a computed MWF?
- **Yes** → copy it to `sub-01/grase/MWF_native.nii.gz` (or `MWF.nii`).
- **No** → run MUMC's `MWFfit` on `sub-01/grase/grase.nii.gz` (needs ExploreDTI),
  then put the result in `sub-01/grase/`.

**2c. Lesion mask** — does MUMC provide one?
- **Yes** → copy to `sub-01/lesion/lesion_mask_native.nii.gz`.
- **No** (and it's an MS patient) → segment from FLAIR:
  ```bash
  $MATLAB -batch "addpath('$REPO/MUMC_pipeline'); \
    run_subject('COHORT_DIR/sub-01','$REPO','<chisep>',{})"   # sets paths only
  # then:
  $MATLAB -batch "mimm_root='$REPO'; output_dir='COHORT_DIR/sub-01'; \
    run('$REPO/MUMC_pipeline/lesion/run_lst.m')"
  ```
- **Control / no lesions** → skip; lesion figures auto-skip.

After getting a native mask + FLAIR:
```bash
bash $REPO/MUMC_pipeline/lesion/register_flair.sh COHORT_DIR/sub-01
```

---

## 3. Run the pipeline (all subjects + cohort)

```bash
bash $REPO/MUMC_pipeline/run_cohort.sh COHORT_DIR $REPO <chisep_dir>
```

Per subject, in order (each step skips cleanly if its input is absent):
`matlab_pre` (prepare+QSM) → `register` (atlas) → `dti` → `matlab_post`
(chisep+MIMM+grase) → `python` (figs) → then cohort aggregation + figures.

Useful flags:
- `--steps python,cohort` — re-run only analysis after a fix
- `--subjects 'sub-0{1,2,3}'` — subset
- `--skip-done` — skip steps whose output exists
- `--dry-run` — print commands only

---

## 4. ⚠️ λ_chi (optional, do once on a reference subject)

Default `lambda_chi = 0.015`. To tune from data (figure 30 showed it under-weights QSM):
```bash
$MATLAB -batch "lock_cohort=true; run('$REPO/MUMC_pipeline/mimm/sweep_lambda_chi.m')"
```
Writes `MUMC_pipeline/lambda_chi_cohort.txt`; every subject then uses that value.
Derive on ONE subject, apply to all (avoids per-subject circularity).
Note: the sweep on the test subject showed the slope ceilings at ~0.4 (a
dictionary-coverage limit), so don't expect lambda alone to fully fix χ–QSM.

---

## 5. Outputs

- Per subject: `sub-XX/figures/` (01–48), `sub-XX/analysis/*.csv`
- Cohort: `COHORT_DIR/cohort/figures/` (31–44), `cohort/analysis/cohort_roi_stats.csv`
- See `analysis/FIGURES.md` for what every figure shows.

---

## Quick reference — what each figure group is

| Figs | Content |
|------|---------|
| 01–15 | Spatial maps (anatomy, QSM, R2*, MVF, FVF, g-ratio, χ, FA, θ, error) |
| 16–17, 30 | Bland-Altman χ vs χ-sep; χ_total vs QSM consistency |
| 18–22 | ROI tables/heatmap/scatter, dual-atlas, atlas comparison |
| 23–25, 45 | Overestimation (ranked, vs θ, vs FA, spatial+JHU overlay) |
| 26–29 | ROI scatters: MVF/χ⁻, χ_iron/χ⁺, χ_myelin/χ⁻, MVF/FA |
| 38–41 | MWF: spatial, MVF-vs-MWF, Bland-Altman, bar (needs GRASE) |
| 46–48 | Lesion overlay, lesion-vs-NAWM, per-lesion MVF-vs-MWF (needs mask) |
| 31–37, 42–44 | Cohort versions (between-subject SEM) |
