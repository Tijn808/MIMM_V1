#!/bin/bash
################################################################################
# Script properties
#
# Name        : 030_PrepareQSM.sh
# Description : Prepare the ME-GRE data and reconstruct QSM + R2* (MIMM step 1).
# Arguments   : subjectName
# Exit code   : 0 if successfully executed, !0 if there was an error
################################################################################

################################################################################
# Version History
#
# 20260611 1.0 Tijn Saes   MIMM pipeline step 1 (prepare ME-GRE + QSM)
################################################################################

################################################################################
# Includes
source $( dirname "$0" )/functions.source
source $( dirname "$0" )/project.config

################################################################################
# Constants
readonly SCRIPTNAME=$( basename $( readlink -f "$0" ) )
readonly SCRIPTVERSION=1.0
export FSLOUTPUTTYPE=NIFTI_GZ          # MIMM works in .nii.gz (template default = NIFTI)

# MIMM configuration (override in project.config if needed)
: "${MIMM_REPO:=$SCRIPTDIR/MIMM_V1}"
: "${CHISEP_DIR:=$SCRIPTDIR/matlab}"
: "${MATLAB_BIN:=matlab}"
readonly MIMM_PIPE="$MIMM_REPO/MUMC_pipeline"

################################################################################
# Variables
subjectName="$1"
logFile="$LOGDIR"/"$subjectName".log
subjectDir="$RESULTDIR"/"$subjectName"
workingDir=~/"$PROJECTNAME"/"$subjectName"

inputFileList=( nifti/gremag.nii nifti/grepha.nii )   # checked with CheckFilesExist
outputList=( magnitude.nii.gz phase.nii.gz qsm )      # copied back (files + folders)

################################################################################
# Start
HeaderLog "START" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"

echo "Project       : ""$PROJECTNAME"
echo " - results     : ""$RESULTDIR"
echo " - subject     : ""$subjectName"
echo " - working dir : ""$workingDir"

# Check the subject directory exists
if [ ! -d "$subjectDir" ]; then
    echo "ERROR - subject directory not found: $subjectDir"
    HeaderLog "END" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"
    exit 1
fi

# Check input files exist
CheckFilesExist ${inputFileList[@]/#/"$subjectDir"/}
if [ $? -ne 0 ]; then
    echo "ERROR - missing input files (run 010 first)"
    HeaderLog "END" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"
    exit 1
fi

# Create local working directory and copy the ME-GRE input into it
CreateWorkingDir "$workingDir"
echo "Copying input from subject directory to local working directory..."
mkdir -p "$workingDir"/nifti
cp "$subjectDir"/nifti/gremag.nii* "$subjectDir"/nifti/grepha.nii* "$workingDir"/nifti/ 2>/dev/null
cp "$subjectDir"/nifti/gremag.json "$subjectDir"/nifti/grepha.json "$workingDir"/nifti/ 2>/dev/null

################################################################################
# Actual image processing part
echo "Running prepare_mgre + QSM in Matlab..."
"$MATLAB_BIN" -nodesktop -nodisplay -nosplash -r \
    "try, addpath('$MIMM_PIPE'); run_subject('$workingDir','$MIMM_REPO','$CHISEP_DIR',{'prepare_mgre','qsm'}); catch e, fprintf(2,'%s\n',e.message); exit(1); end; exit(0)"
if [ $? -ne 0 ]; then
    echo "ERROR - prepare + QSM failed"
    GarbageCleanUp "$workingDir"
    HeaderLog "END" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"
    exit 1
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
