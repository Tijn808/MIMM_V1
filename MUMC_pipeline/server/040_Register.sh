#!/bin/bash
################################################################################
# Script properties
#
# Name        : 040_Register.sh
# Description : Register the HCP1065/JHU atlases into subject space (MIMM step 2).
# Arguments   : subjectName
# Exit code   : 0 if successfully executed, !0 if there was an error
################################################################################

################################################################################
# Version History
#
# 20260611 1.0 Tijn Saes   MIMM pipeline step 2 (atlas registration)
################################################################################

################################################################################
# Includes
source $( dirname "$0" )/functions.source
source $( dirname "$0" )/project.config

################################################################################
# Constants
readonly SCRIPTNAME=$( basename $( readlink -f "$0" ) )
readonly SCRIPTVERSION=1.0
export FSLOUTPUTTYPE=NIFTI_GZ

: "${MIMM_REPO:=$SCRIPTDIR/MIMM_V1}"
readonly MIMM_PIPE="$MIMM_REPO/MUMC_pipeline"

################################################################################
# Variables
subjectName="$1"
logFile="$LOGDIR"/"$subjectName".log
subjectDir="$RESULTDIR"/"$subjectName"
workingDir=~/"$PROJECTNAME"/"$subjectName"

inputFileList=( qsm/mag_e1.nii.gz qsm/brain_mask.nii.gz )   # checked with CheckFilesExist
outputList=( atlas )                                        # (+ dti/ if a DTI scan exists)

################################################################################
# Start
HeaderLog "START" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"

echo "Project       : ""$PROJECTNAME"
echo " - results     : ""$RESULTDIR"
echo " - subject     : ""$subjectName"
echo " - working dir : ""$workingDir"

if [ ! -d "$subjectDir" ]; then
    echo "ERROR - subject directory not found: $subjectDir"
    HeaderLog "END" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"
    exit 1
fi

CheckFilesExist ${inputFileList[@]/#/"$subjectDir"/}
if [ $? -ne 0 ]; then
    echo "ERROR - missing input files (run 030 first)"
    HeaderLog "END" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"
    exit 1
fi

# Create local working directory and copy the qsm/ folder (+ dti/ if present) in
CreateWorkingDir "$workingDir"
echo "Copying input from subject directory to local working directory..."
cp -r "$subjectDir"/qsm "$workingDir"/qsm
[ -d "$subjectDir"/dti ] && cp -r "$subjectDir"/dti "$workingDir"/dti

################################################################################
# Actual image processing part
echo "Running atlas registration..."
bash "$MIMM_PIPE/atlas/register_atlas.sh" "$workingDir"
if [ $? -ne 0 ]; then
    echo "ERROR - atlas registration failed"
    GarbageCleanUp "$workingDir"
    HeaderLog "END" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"
    exit 1
fi

if [ -f "$workingDir"/dti/dti.nii.gz ]; then
    echo "Running DTI preprocessing..."
    bash "$MIMM_PIPE/preprocess/preprocess_dti.sh" "$workingDir"
    [ $? -eq 0 ] && outputList=( atlas dti )
fi

################################################################################
# End - copy output back to the subject directory and clean up
echo "Copying output from local working directory to subject directory..."
for o in "${outputList[@]}"; do
    rm -rf "$subjectDir"/"$o"
    cp -r "$workingDir"/"$o" "$subjectDir"/"$o"
done
GarbageCleanUp "$workingDir"
HeaderLog "END" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"
exit 0
