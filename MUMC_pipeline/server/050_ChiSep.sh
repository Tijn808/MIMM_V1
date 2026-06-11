#!/bin/bash
################################################################################
# 050_ChiSep.sh - susceptibility source separation (chi-separation).
#
# MIMM pipeline step 3.
#   IN : results/<subj>/qsm/  (QSM, R2star)   (from 030)
#   OUT: results/<subj>/chisep/  (chi_neg, chi_pos, chi_tot)
#
# Usage (from scripts/):  sh 050_ChiSep.sh <subjectName>
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
[ -f "$subjectDir/qsm/QSM.nii.gz" ] || fail "qsm/QSM.nii.gz not found (run 030 first)"

echo "[050] chi-separation: $subjectName"
run_matlab "addpath('$MIMM_PIPE'); run_subject('$subjectDir','$MIMM_REPO','$CHISEP_DIR',{'chisep'})" \
    || fail "chi-separation failed"

HeaderLog "END" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"
exit 0
