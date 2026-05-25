#!/bin/bash
# Register HCP1065 FA and V1 atlas into subject ME-GRE space.
# Uses affine (flirt) + nonlinear (fnirt) registration for accurate WM tract alignment.
#
# Outputs (written to <subj_dir>/atlas/):
#   FA_atlas.nii.gz      FA map in subject space
#   V1_atlas.nii.gz      Principal eigenvector [x y z] in subject space
#   theta_atlas.nii.gz   Fibre angle relative to B0 (degrees) in subject space
#
# Usage: bash register_atlas.sh <subj_dir>
#   subj_dir must contain:
#     qsm/brain_mask.nii.gz   — brain mask from BET/VSHARP
#     qsm/mag_e1.nii.gz       — echo 1 magnitude (written by run_QSM.m)

set -e

SUBJ_DIR="${1:?Usage: $0 <subj_dir>}"
FSL=/home/tijn-saes/fsl/bin
STD=/home/tijn-saes/fsl/data/standard
CFG=/home/tijn-saes/fsl/etc/flirtsch
OUT="${SUBJ_DIR}/atlas"
mkdir -p "$OUT"

MAG_E1="${SUBJ_DIR}/qsm/mag_e1.nii.gz"
MASK="${SUBJ_DIR}/qsm/brain_mask.nii.gz"

# --- Step 1: brain-extract echo 1 for registration ---
echo "Brain-extracting echo 1 for registration..."
"$FSL/fslmaths" "$MAG_E1" -mas "$MASK" "${OUT}/mag_e1_brain.nii.gz"

# --- Step 2: affine registration subject → MNI152 1mm ---
echo "Affine registration to MNI152..."
"$FSL/flirt" \
    -in  "${OUT}/mag_e1_brain.nii.gz" \
    -ref "${STD}/MNI152_T1_1mm_brain.nii.gz" \
    -omat "${OUT}/subj2mni_affine.mat" \
    -out  "${OUT}/mag_e1_affine_mni.nii.gz" \
    -dof 12 -cost normcorr

# --- Step 3: nonlinear registration subject → MNI152 (fnirt) ---
echo "Nonlinear registration to MNI152 (fnirt — this takes ~10 min)..."
"$FSL/fnirt" \
    --in="${OUT}/mag_e1_brain.nii.gz" \
    --ref="${STD}/MNI152_T1_1mm_brain.nii.gz" \
    --refmask="${STD}/MNI152_T1_1mm_brain_mask_dil.nii.gz" \
    --aff="${OUT}/subj2mni_affine.mat" \
    --cout="${OUT}/subj2mni_warp.nii.gz" \
    --subsamp=4,2,1,1 \
    --miter=5,5,5,5 \
    --warpres=10,10,10 \
    --infwhm=8,4,2,2 \
    --reffwhm=4,2,0,0 \
    --regmod=bending_energy \
    --lambda=300 \
    --iout="${OUT}/mag_e1_fnirt_mni.nii.gz"

# --- Step 4: invert warp (MNI → subject) ---
echo "Inverting warp field..."
"$FSL/invwarp" \
    -w "${OUT}/subj2mni_warp.nii.gz" \
    -r "$MAG_E1" \
    -o "${OUT}/mni2subj_warp.nii.gz"

# --- Step 5: apply warp to FA atlas ---
echo "Warping FA atlas to subject space..."
"$FSL/applywarp" \
    -i "${STD}/FSL_HCP1065_FA_1mm.nii.gz" \
    -r "$MAG_E1" \
    -w "${OUT}/mni2subj_warp.nii.gz" \
    -o "${OUT}/FA_atlas.nii.gz" \
    --interp=trilinear

# --- Step 6: apply warp to V1 atlas (with vector reorientation) ---
echo "Warping V1 eigenvector atlas to subject space..."
"$FSL/vecreg" \
    -i "${STD}/FSL_HCP1065_V1_1mm.nii.gz" \
    -o "${OUT}/V1_atlas.nii.gz" \
    -r "$MAG_E1" \
    -w "${OUT}/mni2subj_warp.nii.gz"

# --- Step 7: compute theta = acos(|V1_z|) * 180/pi ---
echo "Computing theta map..."
"$FSL/fslsplit" "${OUT}/V1_atlas.nii.gz" "${OUT}/V1_comp_" -t
"$FSL/fslmaths" "${OUT}/V1_comp_0002.nii.gz" -abs "${OUT}/V1_z_abs.nii.gz"
"$FSL/fslmaths" "${OUT}/V1_z_abs.nii.gz" -acos -mul 57.2958 "${OUT}/theta_atlas.nii.gz"

rm -f "${OUT}"/V1_comp_000*.nii.gz "${OUT}/V1_z_abs.nii.gz"

echo "Done. Atlas outputs written to: ${OUT}"
echo "  FA_atlas.nii.gz"
echo "  V1_atlas.nii.gz"
echo "  theta_atlas.nii.gz"
