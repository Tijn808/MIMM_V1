#!/bin/bash
# Convert a MUMC DICOM session to NIfTI and sort it into the pipeline layout.
#
# Uses dcm2niix. Series are matched by name (not number) so it works across
# subjects whose series numbers differ. Each modality is optional — missing
# series are skipped with a message.
#
# Usage:
#   bash sort_dicom.sh <dicom_session_dir> <output_subj_dir>
#
#   <dicom_session_dir>  folder containing the per-series subfolders, e.g.
#                        .../IMPROMYMS_001/20251218-.../  with 501-ME_GRE/ etc.
#   <output_subj_dir>    the subject's pipeline dir (e.g. cohort/sub-01)
#
# Series → destination mapping:
#   *ME_GRE*       → <subj>/            501-ME_GRE_e1.nii.gz ...  (prepare_mgre input)
#   *GRASE*        → <subj>/grase/      grase.nii.gz              (→ MWF via MUMC MWFfit)
#   *FLAIR*        → <subj>/lesion/     FLAIR_native.nii.gz       (lesion segmentation)
#   *DTI_tra* (not Reg) → <subj>/dti/   dti.nii.gz + .bval/.bvec  (orientation)
#   *REVERSE*      → <subj>/dti/        reverse_b0.nii.gz         (topup)
#   *3D_ISO_SAG* (not FLAIR/V3D) → <subj>/t1w/  T1w_native.nii.gz (registration/LST)
#
# After this: fill paths.m, then run the pipeline on <output_subj_dir>.

set -e
SESSION="${1:?Usage: $0 <dicom_session_dir> <output_subj_dir>}"
SUBJ="${2:?Usage: $0 <dicom_session_dir> <output_subj_dir>}"

if [ -n "$FSLDIR" ]; then DCM="$FSLDIR/bin/dcm2niix"; else DCM=dcm2niix; fi
if ! command -v "$DCM" >/dev/null 2>&1 && [ ! -x "$DCM" ]; then
    DCM=$(command -v dcm2niix) || { echo "dcm2niix not found. Set FSLDIR or add to PATH."; exit 1; }
fi
echo "Using dcm2niix: $DCM"

mkdir -p "$SUBJ" "$SUBJ/grase" "$SUBJ/lesion" "$SUBJ/dti" "$SUBJ/t1w"

# Find a series subfolder by name pattern (first match), optionally excluding a pattern.
find_series() {  # $1 = include glob, $2 = optional exclude regex
    local inc="$1" exc="$2"
    while IFS= read -r d; do
        [ -n "$exc" ] && echo "$d" | grep -qiE "$exc" && continue
        echo "$d"; return 0
    done < <(find "$SESSION" -maxdepth 1 -type d -iname "*$inc*" | sort)
    return 1
}

# Convert one series dir → dest dir with a given filename stem.
convert() {  # $1 = series dir, $2 = dest dir, $3 = -f format/stem, $4 = label
    local src="$1" dest="$2" fmt="$3" label="$4"
    if [ -z "$src" ]; then echo "  [skip] $label — series not found"; return; fi
    echo "  [convert] $label  ($(basename "$src"))"
    "$DCM" -z y -b y -f "$fmt" -o "$dest" "$src" >/dev/null 2>&1 || \
        echo "    WARNING: dcm2niix returned nonzero for $label"
}

echo ""
echo "Sorting DICOM → $SUBJ"

# 1. ME-GRE → subject root, named <seriesnum>-<protocol> so prepare_mgre finds
#    501-ME_GRE_e1.nii.gz etc. (dcm2niix appends _e<n> per echo and _ph for phase)
convert "$(find_series ME_GRE)"           "$SUBJ"        '%s-%p' 'ME-GRE (MIMM input)'

# 2. GRASE → grase/grase.nii.gz  (raw; MWF computed separately by MUMC MWFfit)
convert "$(find_series GRASE)"            "$SUBJ/grase"  'grase' 'T2-GRASE (raw, for MWF)'

# 3. FLAIR → lesion/FLAIR_native.nii.gz
convert "$(find_series FLAIR)"            "$SUBJ/lesion" 'FLAIR_native' 'T2-FLAIR (lesions)'

# 4. DTI → dti/dti.nii.gz (+ .bval/.bvec). Exclude registered/derived (Reg) series.
convert "$(find_series DTI_tra 'Reg')"    "$SUBJ/dti"    'dti'   'DTI (orientation)'

# 5. Reverse-PE b0 → dti/reverse_b0.nii.gz (for topup distortion correction)
convert "$(find_series REVERSE)"          "$SUBJ/dti"    'reverse_b0' 'reverse-PE b0 (topup)'

# 6. T1w → t1w/T1w_native.nii.gz. 3D_ISO_SAG, excluding FLAIR and V3D variants.
convert "$(find_series 3D_ISO_SAG 'FLAIR|V3D')" "$SUBJ/t1w" 'T1w_native' 'T1w (registration)'

echo ""
echo "Done. Converted NIfTI in: $SUBJ"
echo ""
echo "VERIFY before running the pipeline:"
echo "  - ls $SUBJ/*ME_GRE_e*.nii.gz   (expect 5 magnitude + 5 phase _ph)"
echo "  - check echo times in the JSON sidecars match 6/12/18/24/30 ms"
echo "  - DTI: confirm dti/dti.bval and dti/dti.bvec were written"
echo "  - set PE_DIR and READOUT_TIME in preprocess_dti.sh from dti/dti.json"
echo "  - MWF: run MUMC MWFfit on grase/grase.nii.gz (or drop in their MWF.nii)"
