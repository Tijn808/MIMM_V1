# Deploying MIMM into the MUMC pipeline

`030_MIMM.sh` is a single pipeline step that runs the whole MIMM chain for one
subject. It is built from `000_scriptTemplate.sh`, so it behaves like a native
MUMC step: it sources `functions.source` and `project.config`, logs with
`HeaderLog`, copies the ME-GRE input to a local working directory, processes it
there, and copies the results back to the subject's results directory.

## One-time setup on the server

1. **Clone this repo** somewhere on the server (it carries the MIMM toolbox and
   dictionary, so it can't be just the one script). Putting it next to the
   scripts is the default the step expects:
   ```sh
   git clone git@github.com:Tijn808/MIMM_V1.git "$SCRIPTDIR"/MIMM_V1
   ```

2. **Copy `030_MIMM.sh`** into the MUMC `scripts/` folder, next to `010`/`020`,
   so it can source `functions.source` and `project.config` from there.

3. **Add the MIMM settings to `project.config`** (or set them at the top of the
   script):
   - `MIMM_REPO`   - where you cloned the repo (default: `$SCRIPTDIR/MIMM_V1`)
   - `CHISEP_DIR`  - the Chisep_Toolbox_v1.2.1 path
   - `MATLAB_BIN`  - the matlab launcher (default: `matlab` on PATH)
   - `FSLDIR`      - the FSL install (usually already set by the FSL environment)

4. **Add one line to `pipeline.sh`**, after the `020` call:
   ```sh
   sh 030_MIMM.sh "$subjectName"
   ```

## How it runs (per subject)

`$RESULTDIR/<subjectName>/` (where `010` wrote the NIfTIs) -> copied to a local
working dir -> `prepare + QSM (MATLAB)` -> `atlas registration (FSL)` ->
`DTI (FSL, only if a dti/ scan is present)` -> `chi-separation + MIMM (MATLAB)`
-> `ROI stats + figures (Python)` -> results copied back to the subject dir:
- `qsm/`      QSM, R2*, brain mask
- `atlas/`    registered FA / theta / JHU labels
- `chisep/`   chi-separation maps
- `mimm/`     MVF, FVF, g-ratio, chi_myelin, chi_iron (basic + atlas)
- `analysis/` roi_stats.csv
- `figures/`  figures 01-45

If anything fails, the step cleans up the working dir, logs `END`, and exits 1.

## One thing to confirm against a real 010 output

`prepare_mgre.m` expects the ME-GRE NIfTIs named `<prefix>-ME_GRE_e1.nii.gz`
... `_e5.nii.gz` plus `_ph` phase files (dcm2niix `%s-%p`; prefix auto-detected),
as `.nii.gz`. The step's input check globs `*ME_GRE_e1.nii.gz`. If
`010_DicomToNifti.sh` names them differently, or writes uncompressed `.nii`,
the `inputGlob` lines near the top of `030_MIMM.sh` and the patterns in
`prepare_mgre.m` need adjusting to match. Share one converted subject's file
listing and this can be set exactly.

## Validation

The full chain (every script this step calls) was run end-to-end on the test
subject and reproduced the reference results exactly: MVF vs chi-sep |chi-|
r = 0.93, iron r = 0.86, splenium overestimation +37%. The `030_MIMM.sh`
orchestration itself was also dry-run-verified to build the correct calls.
