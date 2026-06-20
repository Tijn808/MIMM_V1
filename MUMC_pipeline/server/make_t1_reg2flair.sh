#!/bin/bash
################################################################################
# make_t1_reg2flair.sh - create t1w_reg2flair.nii for one subject.
#
# 020_segmentation expects t1w_reg2flair.nii (T1 resampled into FLAIR space) as
# input, but that file is produced upstream in the other student's pipeline and
# is missing for subjects outside it (e.g. 026-028). This rigidly co-registers
# the subject's T1 to its FLAIR (within-subject, cross-modal) to make it.
#
# Usage (from scripts/):  sh make_t1_reg2flair.sh <subjectName>
#
# NB: a standard rigid (6-dof, normmi) registration. Confirm with Lis that this
# matches how t1w_reg2flair.nii was made for the other subjects, so 026-028 are
# consistent with the rest.
################################################################################
source "$( dirname "$0" )/project.config" 2>/dev/null
export FSLOUTPUTTYPE=NIFTI          # 020_segmentation expects uncompressed .nii

s="$1"
[ -n "$s" ] || { echo "usage: sh make_t1_reg2flair.sh <subjectName>"; exit 1; }
d="$RESULTDIR/$s"
[ -f "$d/nifti/t1w.nii" ]   || { echo "ERROR: $d/nifti/t1w.nii not found";   exit 1; }
[ -f "$d/nifti/flair.nii" ] || { echo "ERROR: $d/nifti/flair.nii not found"; exit 1; }

echo "Registering T1 -> FLAIR for $s ..."
flirt -in "$d/nifti/t1w.nii" -ref "$d/nifti/flair.nii" -dof 6 -cost normmi \
      -out "$d/t1w_reg2flair" -omat "$d/t1w_reg2flair.mat"
if [ -f "$d/t1w_reg2flair.nii" ]; then
    echo "OK: wrote $d/t1w_reg2flair.nii  (now run 020_segmentation $s)"
else
    echo "ERROR: registration did not produce t1w_reg2flair.nii"; exit 1
fi
