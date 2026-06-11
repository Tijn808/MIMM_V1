#!/bin/bash
################################################################################
# 070_ROIstats.sh - per-subject ROI statistics + QC figures.
#
# MIMM pipeline step 5 (last per-subject step). Output feeds 080 (cohort).
#   IN : results/<subj>/mimm/ + atlas/ + chisep/        (060, 040, 050)
#   OUT: results/<subj>/analysis/roi_stats.csv + figures/  (6 QC figures)
#
# MIMM_FULL_FIGURES=1 produces the full per-subject figure set instead of QC-only.
#
# Usage (from scripts/):  sh 070_ROIstats.sh <subjectName>
################################################################################
source "$( dirname "$0" )/functions.source" 2>/dev/null
source "$( dirname "$0" )/project.config"   2>/dev/null
readonly SCRIPTNAME=$( basename "$( readlink -f "$0" )" )
readonly SCRIPTVERSION=1.0

: "${MIMM_REPO:=$SCRIPTDIR/MIMM_V1}"
readonly MIMM_PIPE="$MIMM_REPO/MUMC_pipeline"
ANALYSIS="$MIMM_PIPE/analysis"

subjectName="$1"
subjectDir="$RESULTDIR/$subjectName"

HeaderLog "START" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"
fail() { echo "ERROR - $1"; unset MIMM_OUTPUT_DIR MIMM_QC_ONLY; HeaderLog "END" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"; exit 1; }

[ -n "$subjectName" ] || fail "no subjectName (usage: $SCRIPTNAME <subjectName>)"
[ -f "$subjectDir/mimm/MVF_basic.nii.gz" ] || fail "mimm/MVF_basic.nii.gz not found (run 060 first)"

echo "[070] ROI stats + figures: $subjectName"
export MIMM_OUTPUT_DIR="$subjectDir"
python3 "$ANALYSIS/extract_roi_stats.py" || fail "extract_roi_stats.py failed"
if [ "${MIMM_FULL_FIGURES:-0}" = "1" ]; then
    python3 "$ANALYSIS/generate_figures.py"       || fail "generate_figures.py failed"
    python3 "$ANALYSIS/plot_roi_stats.py"         || fail "plot_roi_stats.py failed"
    python3 "$ANALYSIS/analyse_overestimation.py" || fail "analyse_overestimation.py failed"
else
    MIMM_QC_ONLY=1 python3 "$ANALYSIS/generate_figures.py" || fail "QC figures failed"
fi
unset MIMM_OUTPUT_DIR MIMM_QC_ONLY

HeaderLog "END" "$SCRIPTDIR"/"$SCRIPTNAME" "$SCRIPTVERSION"
exit 0
