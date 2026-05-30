#!/bin/bash
# DTI preprocessing for MUMC data: topup (optional) → eddy → dtifit → register to ME-GRE space
#
# b-values: 0, 1000, 2000
#
# Usage:
#   bash preprocess_dti.sh <subj_dir>                        # no distortion correction
#   bash preprocess_dti.sh <subj_dir> <reverse_b0.nii.gz>   # with topup
#
# subj_dir must contain:
#   dti/dti.nii.gz        4D DTI volume (all b-shells concatenated)
#   dti/dti.bval          b-values (space-separated, one per volume)
#   dti/dti.bvec          gradient directions (3 x N matrix)
#   qsm/mag_e1.nii.gz     ME-GRE echo 1 magnitude (registration reference)
#   qsm/brain_mask.nii.gz brain mask
#
# Outputs (written to <subj_dir>/dti/):
#   FA.nii.gz             FA map in ME-GRE space
#   V1.nii.gz             principal eigenvector in ME-GRE space
#   theta.nii.gz          fibre angle relative to B0 (degrees) in ME-GRE space

set -e

SUBJ_DIR="${1:?Usage: $0 <subj_dir> [reverse_b0.nii.gz]}"
REVERSE_B0="${2:-}"          # optional: reverse-PE b=0 for topup

if [ -n "$FSLDIR" ]; then
    FSL="${FSLDIR}/bin"
else
    FSL=/home/tijn-saes/fsl/bin
fi
DTI_DIR="${SUBJ_DIR}/dti"
MAG_E1="${SUBJ_DIR}/qsm/mag_e1.nii.gz"
MASK="${SUBJ_DIR}/qsm/brain_mask.nii.gz"

# ---------------------------------------------------------------------------
# TODO: set these to match the MUMC DTI sequence parameters
#   PE_DIR:        phase-encode direction (y = AP, -y = PA, x = RL, etc.)
#   READOUT_TIME:  total EPI readout time in seconds
#                  = EchoSpacing_s × (NumberOfPhaseEncodeSteps - 1)
#                  Check the JSON sidecar: TotalReadoutTime or compute from
#                  EffectiveEchoSpacing × ReconMatrixPE
PE_DIR="y"
READOUT_TIME=0.0630
# ---------------------------------------------------------------------------

cd "$DTI_DIR"

# --- Step 1: brain extraction of b=0 for mask ---
echo "Extracting b=0 volume..."
"$FSL/fslroi" dti.nii.gz b0.nii.gz 0 1
"$FSL/bet" b0.nii.gz b0_brain.nii.gz -m -f 0.3

# --- Step 2: topup (only if reverse-PE b=0 provided) ---
if [ -n "$REVERSE_B0" ]; then
    echo "Running topup..."

    # Merge forward and reverse b=0
    "$FSL/fslmerge" -t b0s_for_topup.nii.gz b0.nii.gz "$REVERSE_B0"

    # Acquisition parameters: [PE_vec  readout_time]
    # Row 1 = forward PE, row 2 = reverse PE
    printf "0 1 0 ${READOUT_TIME}\n0 -1 0 ${READOUT_TIME}\n" > acqp.txt

    "$FSL/topup" \
        --imain=b0s_for_topup.nii.gz \
        --datain=acqp.txt \
        --config=b02b0.cnf \
        --out=topup_result \
        --iout=b0s_corrected.nii.gz \
        --fout=topup_fieldmap.nii.gz

    TOPUP_ARGS="--topup=topup_result"
    printf "0 1 0 ${READOUT_TIME}\n" > acqp.txt   # eddy needs 1-row acqp when using topup
else
    echo "No reverse b=0 provided — skipping topup."
    printf "0 1 0 ${READOUT_TIME}\n" > acqp.txt
    TOPUP_ARGS=""
fi

# --- Step 3: build eddy index file (all volumes → acquisition 1) ---
N_VOLS=$("$FSL/fslnvols" dti.nii.gz)
python3 -c "print(' '.join(['1']*${N_VOLS}))" > index.txt

# --- Step 4: eddy (motion + eddy-current correction) ---
echo "Running eddy (${N_VOLS} volumes)..."
"$FSL/eddy_openmp" \
    --imain=dti.nii.gz \
    --mask=b0_brain_mask.nii.gz \
    --acqp=acqp.txt \
    --index=index.txt \
    --bvecs=dti.bvec \
    --bvals=dti.bval \
    --out=dti_eddy \
    $TOPUP_ARGS \
    --repol \
    --verbose

# --- Step 5: dtifit ---
echo "Running dtifit..."
"$FSL/dtifit" \
    --data=dti_eddy \
    --out=dtifit \
    --mask=b0_brain_mask.nii.gz \
    --bvecs=dti_eddy.eddy_rotated_bvecs \
    --bvals=dti.bval

# dtifit outputs: dtifit_FA.nii.gz, dtifit_V1.nii.gz, dtifit_MD.nii.gz, ...

# --- Step 6: register b=0 to ME-GRE space ---
echo "Registering DTI to ME-GRE space..."
"$FSL/flirt" \
    -in  b0_brain.nii.gz \
    -ref "$MAG_E1" \
    -omat dti2mgre.mat \
    -out  b0_in_mgre.nii.gz \
    -dof 6 \
    -cost normcorr

# --- Step 7: apply transform to FA and V1 ---
echo "Applying transform to FA and V1..."
"$FSL/flirt" \
    -in  dtifit_FA.nii.gz \
    -ref "$MAG_E1" \
    -applyxfm -init dti2mgre.mat \
    -out FA.nii.gz \
    -interp trilinear

"$FSL/vecreg" \
    -i dtifit_V1.nii.gz \
    -o V1.nii.gz \
    -r "$MAG_E1" \
    -t dti2mgre.mat

# --- Step 8: compute theta = acos(|V1_z|) * 180/pi ---
# Matches the convention used in register_atlas.sh
echo "Computing theta map..."
"$FSL/fslsplit" V1.nii.gz V1_comp_ -t
"$FSL/fslmaths" V1_comp_0002.nii.gz -abs V1_z_abs.nii.gz
"$FSL/fslmaths" V1_z_abs.nii.gz -acos -mul 57.2958 theta.nii.gz

rm -f V1_comp_000*.nii.gz V1_z_abs.nii.gz

echo "Done. DTI outputs written to: ${DTI_DIR}"
echo "  FA.nii.gz      — FA map in ME-GRE space"
echo "  V1.nii.gz      — principal eigenvector in ME-GRE space"
echo "  theta.nii.gz   — fibre angle relative to B0 (degrees)"
