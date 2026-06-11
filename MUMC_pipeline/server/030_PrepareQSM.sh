#!/bin/bash
################################################################################
# 030_PrepareQSM.sh - prepare ME-GRE and reconstruct QSM + R2* for one subject.
#
# MIMM pipeline step 1 (after 010/020). Output feeds 040.
#   IN : results/<subj>/nifti/gremag.nii + grepha.nii        (from 010)
#   OUT: results/<subj>/qsm/  (QSM, R2star, brain_mask, mag_e1) + magnitude/phase
#
# Usage (from scripts/):  sh 030_PrepareQSM.sh <subjectName>
################################################################################
source "$( dirname "$0" )/functions.source" 2>/dev/null
source "$( dirname "$0" )/project.config"   2>/dev/null
readonly SCRIPTNAME=$( basename "$( readlink -f "$0" )" )
readonly SCRIPTVERSION=1.0
export FSLOUTPUTTYPE=NIFTI_GZ

: "${MIMM_REPO:=$SCRIPTDIR/MIMM_V1}"
: "${CHISEP_DIR:=$SCRIPTDIR/matlab}"
: "${MATLAB_BIN:=matlab}"
: "${FSLDIR:?ERROR - FSLDIR is not set}"
readonly MIMM_PIPE="$MIMM_REPO/MUMC_pipeline"
export FSLDIR
export PATH="$FSLDIR/bin:$PATH"

subjectName="$1"
subjectDir="$RESULTDIR/$subjectName"

HeaderLog "START" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"
fail() { echo "ERROR - $1"; HeaderLog "END" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"; exit 1; }
run_matlab() {
    "$MATLAB_BIN" -nodesktop -nodisplay -nosplash -r \
        "try, $1, catch e, fprintf(2,'MATLAB ERROR: %s\n',e.message); exit(1); end; exit(0)"
}

[ -n "$subjectName" ] || fail "no subjectName (usage: $SCRIPTNAME <subjectName>)"
[ -d "$subjectDir" ]  || fail "subject dir not found: $subjectDir"
ls "$subjectDir"/nifti/gremag.nii* >/dev/null 2>&1 || fail "input nifti/gremag.nii not found in $subjectDir (run 010 first)"

echo "[030] prepare + QSM: $subjectName"
run_matlab "addpath('$MIMM_PIPE'); run_subject('$subjectDir','$MIMM_REPO','$CHISEP_DIR',{'prepare_mgre','qsm'})" \
    || fail "prepare + QSM failed"

HeaderLog "END" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"
exit 0
