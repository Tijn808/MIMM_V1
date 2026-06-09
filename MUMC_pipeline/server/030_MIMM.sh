#!/bin/bash
################################################################################
# Script properties
#
# Name        : 030_MIMM.sh
# Description : Microstructure-Informed Myelin Mapping (MIMM) for one subject.
#               Runs the full chain on the ME-GRE data produced by 010:
#                 1. prepare ME-GRE + QSM        (MATLAB)
#                 2. atlas registration          (FSL)
#                 3. DTI preprocessing           (FSL, only if a DTI scan is present)
#                 4. chi-separation + MIMM       (MATLAB)
#                 5. ROI statistics + figures    (Python)
#               Follows the MUMC pipeline template (000_scriptTemplate.sh):
#               inputs are copied to a local working dir, processed there, and
#               the results copied back to the subject's results directory.
#
# Arguments   : $1 = subjectName
# Exit code   : 0 if successfully executed, !0 if there was an error
################################################################################

################################################################################
# Version History
#
# 20260608 1.0 Tijn Saes - Initial MIMM pipeline step
################################################################################

################################################################################
# Includes
source $( dirname "$0" )/functions.source
source $( dirname "$0" )/project.config

################################################################################
# Constants
readonly SCRIPTNAME=$(basename $( readlink -f "$0" ))
readonly SCRIPTVERSION=1.0
# MIMM works in .nii.gz throughout (overrides the template's NIFTI default).
export FSLOUTPUTTYPE=NIFTI_GZ

################################################################################
# MIMM-specific configuration
#  Add these to project.config (preferred) or set them here.
#  - MIMM_REPO  : where the MIMM_V1 repo is cloned (has MIMM_set_path.m +
#                 MUMC_pipeline/). Defaults to a clone next to the scripts.
#  - CHISEP_DIR : the Chisep_Toolbox_v1.2.1 directory
#  - MATLAB_BIN : matlab launcher (default: matlab on PATH)
#  - FSLDIR     : FSL install (normally already set by the FSL environment)
: "${MIMM_REPO:=$SCRIPTDIR/MIMM_V1}"
: "${CHISEP_DIR:?ERROR - set CHISEP_DIR (Chisep_Toolbox path) in project.config}"
: "${MATLAB_BIN:=matlab}"
: "${FSLDIR:?ERROR - FSLDIR is not set (needed for FSL and the atlases)}"
readonly MIMM_PIPE="$MIMM_REPO/MUMC_pipeline"
export FSLDIR
export PATH="$FSLDIR/bin:$PATH"

################################################################################
# Variables
subjectName="$1"
logFile="$LOGDIR"/"$subjectName".log
subjectDir="$RESULTDIR"/"$subjectName"
workingDir=~/"$PROJECTDIR"/"$subjectName"

# ME-GRE inputs written by MUMC's 010_DicomToNifti: two 4D files in nifti/
#   nifti/gremag.nii  (magnitude, all echoes)
#   nifti/grepha.nii  (phase,     all echoes)
# prepare_mgre.m reads this layout directly.
niftiDir="$subjectDir/nifti"
inputMag="$niftiDir/gremag.nii"
inputPha="$niftiDir/grepha.nii"

################################################################################
# Start
HeaderLog "START" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"

echo "Project       : $PROJECTNAME"
echo " - results     : $RESULTDIR"
echo " - subject     : $subjectName"
echo " - subject dir : $subjectDir"
echo " - working dir : $workingDir"
echo " - MIMM repo   : $MIMM_REPO"

# fail() - report an error, clean up, log END, and exit 1
fail() {
    echo "ERROR - $1"
    GarbageCleanUp "$workingDir" 2>/dev/null
    HeaderLog "END" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"
    exit 1
}

[ -n "$subjectName" ] || fail "no subjectName given (usage: $SCRIPTNAME <subjectName>)"
[ -d "$subjectDir" ]  || fail "subject directory not found: $subjectDir"

# Check the ME-GRE input exists (accept .nii or .nii.gz)
{ [ -f "$inputMag" ] || [ -f "$inputMag.gz" ]; } || fail "ME-GRE magnitude (nifti/gremag.nii) not found in $subjectDir"
{ [ -f "$inputPha" ] || [ -f "$inputPha.gz" ]; } || fail "ME-GRE phase (nifti/grepha.nii) not found in $subjectDir"

################################################################################
# Copy inputs to a local working directory
CreateWorkingDir "$workingDir"

echo "Copying ME-GRE input to local working directory..."
cp "$niftiDir"/gremag.nii* "$niftiDir"/grepha.nii* "$workingDir"/ 2>/dev/null || fail "could not copy gremag/grepha NIfTIs"
cp "$niftiDir"/gremag.json "$niftiDir"/grepha.json "$workingDir"/ 2>/dev/null
# Copy a DTI scan too if one exists (enables the DTI-informed variant)
if [ -d "$subjectDir/dti" ]; then
    cp -r "$subjectDir/dti" "$workingDir"/dti
fi

################################################################################
# Actual image processing — the MIMM chain (runs in $workingDir)

# run_matlab <statement> - run MATLAB headless with error handling
run_matlab() {
    "$MATLAB_BIN" -nodesktop -nodisplay -nosplash -r \
        "try, $1, catch e, fprintf(2,'MATLAB ERROR: %s\n',e.message); exit(1); end; exit(0)"
}

echo "1/5 prepare + QSM (MATLAB)..."
run_matlab "addpath('$MIMM_PIPE'); run_subject('$workingDir','$MIMM_REPO','$CHISEP_DIR',{'prepare_mgre','qsm'})" \
    || fail "prepare + QSM step failed"

echo "2/5 atlas registration (FSL)..."
bash "$MIMM_PIPE/atlas/register_atlas.sh" "$workingDir" || fail "atlas registration failed"

if [ -f "$workingDir/dti/dti.nii.gz" ]; then
    echo "3/5 DTI preprocessing (FSL)..."
    bash "$MIMM_PIPE/preprocess/preprocess_dti.sh" "$workingDir" || fail "DTI preprocessing failed"
else
    echo "3/5 DTI - skipped (no dti/dti.nii.gz); using atlas orientation"
fi

echo "4/5 chi-separation + MIMM (MATLAB)..."
run_matlab "addpath('$MIMM_PIPE'); run_subject('$workingDir','$MIMM_REPO','$CHISEP_DIR',{'chisep','mimm'})" \
    || fail "chi-separation + MIMM step failed"

echo "5/5 ROI stats + figures (Python)..."
export MIMM_OUTPUT_DIR="$workingDir"
python3 "$MIMM_PIPE/analysis/extract_roi_stats.py"     || fail "extract_roi_stats.py failed"
python3 "$MIMM_PIPE/analysis/generate_figures.py"      || fail "generate_figures.py failed"
python3 "$MIMM_PIPE/analysis/plot_roi_stats.py"        || fail "plot_roi_stats.py failed"
python3 "$MIMM_PIPE/analysis/analyse_overestimation.py" || fail "analyse_overestimation.py failed"
unset MIMM_OUTPUT_DIR

################################################################################
# End — copy results back to the subject directory, then clean up
echo "Copying results to subject directory..."
for d in qsm atlas chisep mimm analysis figures; do
    if [ -d "$workingDir/$d" ]; then
        rm -rf "$subjectDir/$d"
        cp -r "$workingDir/$d" "$subjectDir/$d" || fail "could not copy results ($d) back to $subjectDir"
    fi
done

GarbageCleanUp "$workingDir"
HeaderLog "END" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"
exit 0
