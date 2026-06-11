#!/bin/bash
################################################################################
# 060_MIMM.sh - Microstructure-Informed Myelin Mapping (Basic + Atlas).
#
# MIMM pipeline step 4.
#   IN : results/<subj>/qsm/ (QSM, R2star) + results/<subj>/atlas/ (theta)   (030, 040)
#   OUT: results/<subj>/mimm/  (MVF, FVF, g_ratio, chi_myelin, chi_iron, ...)
#
# Usage (from scripts/):  sh 060_MIMM.sh <subjectName>
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
[ -f "$subjectDir/qsm/QSM.nii.gz" ]            || fail "qsm/QSM.nii.gz not found (run 030 first)"
[ -f "$subjectDir/atlas/theta_atlas.nii.gz" ]  || fail "atlas/theta_atlas.nii.gz not found (run 040 first)"

echo "[060] MIMM: $subjectName"
run_matlab "addpath('$MIMM_PIPE'); run_subject('$subjectDir','$MIMM_REPO','$CHISEP_DIR',{'mimm'})" \
    || fail "MIMM failed"

HeaderLog "END" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"
exit 0
