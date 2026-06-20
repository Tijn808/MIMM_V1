#!/bin/bash
################################################################################
# run_lesion.sh - lesion analysis on Lis's WMH masks.
#
# Per subject with WMH_mask_copy.nii.gz: stages the FLAIR + mask into lesion/,
# registers the mask into ME-GRE space (lesion/register_flair.sh), and runs the
# lesion-vs-NAWM analysis (analysis/plot_lesion.py).
#
# Usage (from scripts/):
#   sh run_lesion.sh IMPROMYMS_003     # one subject (test this first)
#   sh run_lesion.sh                   # every patient that has a WMH mask
################################################################################
source "$( dirname "$0" )/project.config" 2>/dev/null
: "${MIMM_REPO:=$SCRIPTDIR/MIMM_V1}"
PIPE="$MIMM_REPO/MUMC_pipeline"
export FSLOUTPUTTYPE=NIFTI_GZ
: "${FSLDIR:?ERROR - FSLDIR not set}"
export PATH="$FSLDIR/bin:$PATH"

do_one() {
    local s="$1" d="$RESULTDIR/$1"
    local mask="$d/WMH_mask_copy.nii.gz"
    [ -f "$mask" ]                  || { echo "[skip] $s (no WMH_mask_copy.nii.gz)"; return; }
    [ -f "$d/nifti/flair.nii" ]     || { echo "[skip] $s (no nifti/flair.nii)";     return; }
    [ -f "$d/qsm/mag_e1.nii.gz" ]   || { echo "[skip] $s (no qsm/mag_e1.nii.gz)";    return; }
    echo "=== lesion analysis: $s ==="
    mkdir -p "$d/lesion"
    "$FSLDIR/bin/imcp" "$d/nifti/flair.nii" "$d/lesion/FLAIR_native"
    "$FSLDIR/bin/imcp" "$mask"              "$d/lesion/lesion_mask_native"
    bash "$PIPE/lesion/register_flair.sh" "$d" || { echo "$s: registration failed"; return; }
    MIMM_OUTPUT_DIR="$d" python3 "$PIPE/analysis/plot_lesion.py" || echo "$s: plot_lesion failed"
}

if [ -n "$1" ]; then
    do_one "$1"
else
    for dd in "$RESULTDIR"/*/; do do_one "$( basename "$dd" )"; done
fi
echo "done"
