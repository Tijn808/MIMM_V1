#!/bin/bash
################################################################################
# run_cohort_mumc.sh
#
# Processes the whole cohort: for every subject with raw DICOM
# (sourcedata/dicom/IMPROMYMS_*) it runs 010 (DICOM->NIfTI) first if needed, then
# 030_MIMM.sh, one subject after another. Sequential on purpose (each is a heavy
# MATLAB/FSL job). Subjects that already have MIMM output are skipped, and a
# failure on one subject is logged and does not stop the rest.
#
# Lives in scripts/ next to 030_MIMM.sh (sources functions.source/project.config).
#
# Usage (from scripts/):
#   sh run_cohort_mumc.sh --dry-run        # list subjects it would process
#   nohup sh run_cohort_mumc.sh &          # run unattended; progress in cohort_run.log
#   sh run_cohort_mumc.sh IMPROMYMS_002 IMPROMYMS_004   # only these subjects
################################################################################

source "$( dirname "$0" )/functions.source" 2>/dev/null
source "$( dirname "$0" )/project.config"   2>/dev/null
SELF_DIR=$( dirname "$( readlink -f "$0" )" )
LOG="$SELF_DIR/cohort_run.log"

# Parse args: --dry-run flag, plus any explicit subject names.
DRY=0
explicit=""
for a in "$@"; do
    if [ "$a" = "--dry-run" ]; then DRY=1; else explicit="$explicit $a"; fi
done

# Subject list: explicit args, otherwise every subject that has raw DICOM
# (sourcedata/dicom/IMPROMYMS_*). For each we run 010 (DICOM->NIfTI) first if it
# has not been converted yet, then MIMM. So this processes the whole cohort,
# converting any subject that still needs it.
if [ -n "$explicit" ]; then
    subjects="$explicit"
else
    subjects=""
    for d in "$SOURCEDATADIR"/dicom/IMPROMYMS_*; do
        [ -d "$d" ] || continue
        subjects="$subjects $( basename "$d" )"
    done
    subjects=$( echo "$subjects" | tr ' ' '\n' | sort -u | tr '\n' ' ' )
fi

if [ -z "$( echo "$subjects" | tr -d ' ' )" ]; then
    echo "No subjects found in $SOURCEDATADIR/dicom/IMPROMYMS_*."
    exit 1
fi

echo "Cohort subjects:$subjects"
if [ "$DRY" = "1" ]; then echo "(dry-run) nothing executed."; exit 0; fi

have_gre() { [ -f "$RESULTDIR/$1/nifti/gremag.nii" ] || [ -f "$RESULTDIR/$1/nifti/gremag.nii.gz" ]; }

echo "==== cohort run started $( date ) ====" | tee -a "$LOG"
for s in $subjects; do
    if [ -f "$RESULTDIR/$s/mimm/MVF_basic.nii.gz" ]; then
        echo "[skip] $s (already has MIMM output)" | tee -a "$LOG"
        continue
    fi
    # Convert DICOM -> NIfTI first if this subject has not been converted yet.
    if ! have_gre "$s"; then
        echo "---- $s : 010 DICOM->NIfTI $( date ) ----" | tee -a "$LOG"
        sh "$SELF_DIR/010_DicomToNifti.sh" "$s"
    fi
    if ! have_gre "$s"; then
        echo "[skip] $s (no ME-GRE found after 010)" | tee -a "$LOG"
        continue
    fi
    echo "---- $s : MIMM start $( date ) ----" | tee -a "$LOG"
    sh "$SELF_DIR/030_MIMM.sh" "$s"
    rc=$?
    if [ "$rc" -eq 0 ]; then
        echo "---- $s : DONE $( date ) ----" | tee -a "$LOG"
    else
        echo "---- $s : FAILED (exit $rc) $( date ) ----" | tee -a "$LOG"
    fi
done
echo "==== cohort run finished $( date ) ====" | tee -a "$LOG"
