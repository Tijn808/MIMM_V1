#!/bin/bash
################################################################################
# run_mwf_compare.sh - register each subject's MWF map (T2-GRASE space) into
# ME-GRE space so it can be compared per-ROI with MIMM MVF.
#
# Rigid (6-dof) FLIRT of the GRASE first volume to qsm/mag_e1, transform applied
# to MWF.nii -> grase/MWF.nii.gz (1 mm, same grid as the MVF maps). FSL only.
# Runs on whoever has MWF.nii -- no need for the whole cohort.
#
# Usage (from scripts/):
#   sh run_mwf_compare.sh IMPROMYMS_003   # one subject (test)
#   sh run_mwf_compare.sh                 # every subject with MWF.nii
################################################################################
source "$( dirname "$0" )/project.config" 2>/dev/null
export FSLOUTPUTTYPE=NIFTI_GZ
: "${FSLDIR:?ERROR - FSLDIR not set}"
export PATH="$FSLDIR/bin:$PATH"

do_one() {
    local s="$1" d="$RESULTDIR/$1"
    [ -f "$d/MWF.nii" ]            || { echo "[skip] $s (no MWF.nii)";          return; }
    [ -f "$d/nifti/grase.nii" ]    || { echo "[skip] $s (no nifti/grase.nii)";  return; }
    [ -f "$d/qsm/mag_e1.nii.gz" ]  || { echo "[skip] $s (no qsm/mag_e1.nii.gz)"; return; }
    echo "=== MWF -> ME-GRE: $s ==="
    mkdir -p "$d/grase"
    "$FSLDIR/bin/fslroi" "$d/nifti/grase.nii" "$d/grase/grase_vol0" 0 1
    "$FSLDIR/bin/flirt" -in "$d/grase/grase_vol0" -ref "$d/qsm/mag_e1.nii.gz" \
        -omat "$d/grase/grase2mgre.mat" -out "$d/grase/grase_mgre.nii.gz" \
        -dof 6 -cost normmi
    "$FSLDIR/bin/flirt" -in "$d/MWF.nii" -ref "$d/qsm/mag_e1.nii.gz" \
        -applyxfm -init "$d/grase/grase2mgre.mat" -out "$d/grase/MWF.nii.gz" \
        -interp trilinear
    [ -f "$d/grase/MWF.nii.gz" ] && echo "OK $s -> grase/MWF.nii.gz" || echo "FAILED $s"
}

if [ -n "$1" ]; then do_one "$1"; else
    for dd in "$RESULTDIR"/*/; do do_one "$( basename "$dd" )"; done
fi
echo "done"
