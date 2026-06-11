#!/bin/bash
################################################################################
# 040_Register.sh - register the HCP1065/JHU atlases into subject space.
#
# MIMM pipeline step 2. Output feeds 060 (MIMM).
#   IN : results/<subj>/qsm/mag_e1.nii.gz, brain_mask.nii.gz   (from 030)
#   OUT: results/<subj>/atlas/  (FA_atlas, theta_atlas, JHU labels/tracts)
#        (+ results/<subj>/dti/  if a DTI scan is present)
#
# Usage (from scripts/):  sh 040_Register.sh <subjectName>
################################################################################
source "$( dirname "$0" )/functions.source" 2>/dev/null
source "$( dirname "$0" )/project.config"   2>/dev/null
readonly SCRIPTNAME=$( basename "$( readlink -f "$0" )" )
readonly SCRIPTVERSION=1.0
export FSLOUTPUTTYPE=NIFTI_GZ

: "${MIMM_REPO:=$SCRIPTDIR/MIMM_V1}"
: "${FSLDIR:?ERROR - FSLDIR is not set}"
readonly MIMM_PIPE="$MIMM_REPO/MUMC_pipeline"
export FSLDIR
export PATH="$FSLDIR/bin:$PATH"

subjectName="$1"
subjectDir="$RESULTDIR/$subjectName"

HeaderLog "START" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"
fail() { echo "ERROR - $1"; HeaderLog "END" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"; exit 1; }

[ -n "$subjectName" ] || fail "no subjectName (usage: $SCRIPTNAME <subjectName>)"
[ -f "$subjectDir/qsm/mag_e1.nii.gz" ] || fail "qsm/mag_e1.nii.gz not found (run 030 first)"

echo "[040] atlas registration: $subjectName"
bash "$MIMM_PIPE/atlas/register_atlas.sh" "$subjectDir" || fail "atlas registration failed"

if [ -f "$subjectDir/dti/dti.nii.gz" ]; then
    echo "[040] DTI preprocessing: $subjectName"
    bash "$MIMM_PIPE/preprocess/preprocess_dti.sh" "$subjectDir" || fail "DTI preprocessing failed"
else
    echo "[040] no dti/dti.nii.gz - using atlas orientation"
fi

HeaderLog "END" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"
exit 0
