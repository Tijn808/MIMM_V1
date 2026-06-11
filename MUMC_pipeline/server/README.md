# Deploying the MIMM pipeline into the MUMC numbered pipeline

The MIMM processing is split into numbered steps that follow the MUMC template
(`000_scriptTemplate.sh`): each takes `<subjectName>`, sources
`functions.source`/`project.config`, logs with `HeaderLog`, checks its inputs
with `CheckFilesExist`, then `CreateWorkingDir` → process locally → copy results
back to `results/<subject>/` → `GarbageCleanUp`. Every step's output is the next
step's input.

(One adaptation to the template: our outputs are whole folders — `qsm/`,
`atlas/`, `mimm/`, … — so the folder copy in/out uses `cp -r` instead of the
file-only `CopyFilesToDir`, and `FSLOUTPUTTYPE` is `NIFTI_GZ`, not `NIFTI`.)

## The steps

| Step | Script | Reads | Writes |
|------|--------|-------|--------|
| 030 | `030_PrepareQSM.sh` | `nifti/gremag.nii`,`grepha.nii` (010) | `magnitude.nii.gz`,`phase.nii.gz`,`qsm/` |
| 040 | `040_Register.sh` | `qsm/mag_e1`,`brain_mask` | `atlas/` (+`dti/` if a DTI scan exists) |
| 050 | `050_ChiSep.sh` | `magnitude.nii.gz`,`phase.nii.gz` | `chisep/` |
| 060 | `060_MIMM.sh` | `magnitude.nii.gz`,`qsm/`,`atlas/` | `mimm/` |
| 070 | `070_ROIstats.sh` | `mimm/`,`atlas/`,`chisep/`,`qsm/` | `analysis/roi_stats.csv`,`figures/` (6 QC) |
| 080 | `080_cohort.sh` | all `analysis/roi_stats.csv` | `cohort_analysis/` |

`010` (DICOM→NIfTI) and `020` (MWF) are MUMC's own steps and run before these.

## One-time setup on the server

1. **Clone this repo** (carries the toolbox + dictionary):
   ```sh
   git clone https://github.com/Tijn808/MIMM_V1.git "$SCRIPTDIR"/MIMM_V1
   ```
2. **Copy the step scripts into `scripts/`** (file manager is most reliable):
   `030_PrepareQSM.sh`, `040_Register.sh`, `050_ChiSep.sh`, `060_MIMM.sh`,
   `070_ROIstats.sh`, `080_cohort.sh`, `run_cohort_mumc.sh` — all from
   `MIMM_V1/MUMC_pipeline/server/`.
3. **Add to `project.config`** (or leave the defaults):
   - `MIMM_REPO`   (default `$SCRIPTDIR/MIMM_V1`)
   - `CHISEP_DIR`  (default `$SCRIPTDIR/matlab` — searched for a working ROMEO)
   - `MATLAB_BIN`, `FSLDIR` (FSLDIR is usually already set)

## Running

**Whole cohort, one command** (converts with 010 if needed, then 030→070 per
subject, skips finished ones, continues past failures):
```sh
sh run_cohort_mumc.sh --dry-run        # preview the subject list
nohup sh run_cohort_mumc.sh &          # run unattended; progress in cohort_run.log
```
Then once subjects are processed:
```sh
sh 080_cohort.sh                       # cohort myelin-per-ROI flags / patient-vs-HC
```

**A single subject, step by step:**
```sh
sh 030_PrepareQSM.sh IMPROMYMS_002
sh 040_Register.sh   IMPROMYMS_002
sh 050_ChiSep.sh     IMPROMYMS_002
sh 060_MIMM.sh       IMPROMYMS_002
sh 070_ROIstats.sh   IMPROMYMS_002
```
Each step checks its input exists and fails clearly if a previous step is missing.

## Notes
- Per subject, `070` produces only the **6 QC figures** by default. Set
  `MIMM_FULL_FIGURES=1` for the full per-subject figure set.
- The chi-sep toolbox needs a ROMEO binary that runs on this machine; `030`/`050`
  auto-pick a working build from under `CHISEP_DIR`.
- Validated: the chain reproduces the test-subject results (MVF vs chi-sep
  r=0.93, splenium overestimation +37%) and runs end-to-end on the radstation.
