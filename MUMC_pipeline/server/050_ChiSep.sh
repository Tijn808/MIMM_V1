#!/bin/bash
################################################################################
# Script properties
#
# Name        : 050_ChiSep.sh
# Description : Susceptibility source separation / chi-separation (MIMM step 3).
# Arguments   : subjectName
# Exit code   : 0 if successfully executed, !0 if there was an error
################################################################################

################################################################################
# Version History
#
# 20260611 1.0 Tijn Saes   MIMM pipeline step 3 (chi-separation)
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
: "${CHISEP_DIR:=$SCRIPTDIR/matlab}"
: "${MATLAB_BIN:=matlab}"
readonly MIMM_PIPE="$MIMM_REPO/MUMC_pipeline"

################################################################################
# Variables
subjectName="$1"
logFile="$LOGDIR"/"$subjectName".log
subjectDir="$RESULTDIR"/"$subjectName"
workingDir=~/"$PROJECTNAME"/"$subjectName"

inputFileList=( magnitude.nii.gz phase.nii.gz )   # checked with CheckFilesExist
outputList=( chisep )

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

CreateWorkingDir "$workingDir"
echo "Copying input from subject directory to local working directory..."
cp "$subjectDir"/magnitude.nii.gz "$subjectDir"/phase.nii.gz "$workingDir"/

################################################################################
# Actual image processing part
echo "Running chi-separation in Matlab..."
"$MATLAB_BIN" -nodesktop -nodisplay -nosplash -r \
    "try, addpath('$MIMM_PIPE'); run_subject('$workingDir','$MIMM_REPO','$CHISEP_DIR',{'chisep'}); catch e, fprintf(2,'%s\n',e.message); exit(1); end; exit(0)"
if [ $? -ne 0 ]; then
    echo "ERROR - chi-separation failed"
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
