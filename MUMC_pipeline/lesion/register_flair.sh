#!/bin/bash
# Register T2-FLAIR to ME-GRE space and apply transform to lesion mask.
#
# Inputs expected in <subj_dir>/lesion/:
#   FLAIR_native.nii.gz          T2-FLAIR in native scanner space
#   lesion_mask_native.nii.gz    Binary lesion mask in FLAIR/native space
#                                (output of run_lst.m or pre-computed)
#
# Outputs (written to <subj_dir>/lesion/):
#   FLAIR_mgre.nii.gz            FLAIR registered to ME-GRE space
#   FLAIR2mgre.mat               Registration matrix
#   lesion_mask.nii.gz           Lesion mask in ME-GRE space  ← used by pipeline
#
# Usage:
#   bash register_flair.sh <subj_dir>

set -e
SUBJ_DIR="${1:?Usage: $0 <subj_dir>}"

if [ -z "$FSLDIR" ]; then echo "ERROR: FSLDIR is not set." >&2; exit 1; fi
FSL="${FSLDIR}/bin"

LESION_DIR="$SUBJ_DIR/lesion"
REF="$SUBJ_DIR/qsm/mag_e1.nii.gz"
FLAIR="$LESION_DIR/FLAIR_native.nii.gz"
MASK_IN="$LESION_DIR/lesion_mask_native.nii.gz"

if [ ! -f "$FLAIR" ]; then
    echo "FLAIR_native.nii.gz not found in $LESION_DIR"
    echo "Copy the T2-FLAIR NIfTI there and re-run."
    exit 1
fi
if [ ! -f "$REF" ]; then
    echo "ME-GRE reference not found: $REF"
    echo "Run run_QSM_chisep.m first."
    exit 1
fi

# --- Step 1: register FLAIR to ME-GRE (rigid, 6 DOF) ---
echo "Registering FLAIR to ME-GRE space..."
"$FSL/flirt" \
    -in  "$FLAIR" \
    -ref "$REF" \
    -omat "$LESION_DIR/FLAIR2mgre.mat" \
    -out  "$LESION_DIR/FLAIR_mgre.nii.gz" \
    -dof 6 -cost normcorr
echo "  Written: FLAIR_mgre.nii.gz"

# --- Step 2: apply transform to lesion mask (nearest-neighbour) ---
if [ -f "$MASK_IN" ]; then
    echo "Applying transform to lesion mask..."
    "$FSL/flirt" \
        -in   "$MASK_IN" \
        -ref  "$REF" \
        -applyxfm -init "$LESION_DIR/FLAIR2mgre.mat" \
        -out  "$LESION_DIR/lesion_mask_raw.nii.gz" \
        -interp nearestneighbour
    # Threshold to binary (interpolation artefacts can give non-integer values)
    "$FSL/fslmaths" "$LESION_DIR/lesion_mask_raw.nii.gz" \
        -thr 0.5 -bin "$LESION_DIR/lesion_mask.nii.gz"
    rm -f "$LESION_DIR/lesion_mask_raw.nii.gz"

    # Sanity check
    N_LESION=$("$FSL/fslstats" "$LESION_DIR/lesion_mask.nii.gz" -V | awk '{print $1}')
    echo "  Written: lesion_mask.nii.gz  (${N_LESION} lesion voxels)"
    if [ "$N_LESION" -lt 10 ]; then
        echo "  WARNING: very few lesion voxels — check registration and mask."
    fi
else
    echo "  lesion_mask_native.nii.gz not found — skipping mask registration."
    echo "  Run run_lst.m first, then re-run this script."
fi

echo "Done. Lesion outputs in: $LESION_DIR"
