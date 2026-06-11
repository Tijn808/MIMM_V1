#!/bin/bash
################################################################################
# Script properties
#
# Name        : 070_ROIstats.sh
# Description : Per-subject ROI statistics + QC figures (MIMM step 5, last per-subject).
# Arguments   : subjectName
# Exit code   : 0 if successfully executed, !0 if there was an error
#
# MIMM_FULL_FIGURES=1 produces the full per-subject figure set instead of QC-only.
################################################################################

################################################################################
# Version History
#
# 20260611 1.0 Tijn Saes   MIMM pipeline step 5 (ROI stats + QC figures)
################################################################################

################################################################################
# Includes
source $( dirname "$0" )/functions.source
source $( dirname "$0" )/project.config

################################################################################
# Constants
readonly SCRIPTNAME=$( basename $( readlink -f "$0" ) )
readonly SCRIPTVERSION=1.0

: "${MIMM_REPO:=$SCRIPTDIR/MIMM_V1}"
readonly MIMM_PIPE="$MIMM_REPO/MUMC_pipeline"
readonly ANALYSIS="$MIMM_PIPE/analysis"

################################################################################
# Variables
subjectName="$1"
logFile="$LOGDIR"/"$subjectName".log
subjectDir="$RESULTDIR"/"$subjectName"
workingDir=~/"$PROJECTNAME"/"$subjectName"

inputFileList=( mimm/MVF_basic.nii.gz atlas/theta_atlas.nii.gz )
outputList=( analysis figures )

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
    echo "ERROR - missing input files (run 060 first)"
    HeaderLog "END" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"
    exit 1
fi

# Create local working directory and copy mimm/ + atlas/ + chisep/ + qsm/ in
CreateWorkingDir "$workingDir"
echo "Copying input from subject directory to local working directory..."
cp -r "$subjectDir"/mimm  "$workingDir"/mimm
cp -r "$subjectDir"/atlas "$workingDir"/atlas
[ -d "$subjectDir"/chisep ] && cp -r "$subjectDir"/chisep "$workingDir"/chisep
[ -d "$subjectDir"/qsm ]    && cp -r "$subjectDir"/qsm    "$workingDir"/qsm

################################################################################
# Actual image processing part
echo "Extracting ROI statistics and generating figures..."
export MIMM_OUTPUT_DIR="$workingDir"
rc=0
python3 "$ANALYSIS/extract_roi_stats.py" || rc=1
if [ "$rc" -eq 0 ]; then
    if [ "${MIMM_FULL_FIGURES:-0}" = "1" ]; then
        python3 "$ANALYSIS/generate_figures.py"       || rc=1
        python3 "$ANALYSIS/plot_roi_stats.py"         || rc=1
        python3 "$ANALYSIS/analyse_overestimation.py" || rc=1
    else
        MIMM_QC_ONLY=1 python3 "$ANALYSIS/generate_figures.py" || rc=1
    fi
fi
unset MIMM_OUTPUT_DIR MIMM_QC_ONLY
if [ "$rc" -ne 0 ]; then
    echo "ERROR - ROI stats / figures failed"
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
