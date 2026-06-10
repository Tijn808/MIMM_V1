#!/bin/bash
################################################################################
# 040_cohort.sh - cohort-level analysis across all processed subjects.
#
# Run once after the cohort is processed (every results/<subj> has
# analysis/roi_stats.csv). Currently: myelin-per-ROI overview with outlier
# flagging (lesion candidates, pre-ground-truth). More cohort analyses can be
# added here later.
#
# Output: results/cohort_analysis/
#   cohort_mvf_zscore_heatmap.png   subjects x ROIs MVF z-scores (× = flagged)
#   cohort_mvf_flags.csv            flagged subject/ROI combos
#   cohort_mvf_flagcount.png        low-myelin flags per subject
#   cohort_mvf_matrix.csv           raw subject x ROI MVF table
#
# Usage (from scripts/):  sh 040_cohort.sh
################################################################################

source "$( dirname "$0" )/functions.source" 2>/dev/null
source "$( dirname "$0" )/project.config"   2>/dev/null
: "${MIMM_REPO:=$SCRIPTDIR/MIMM_V1}"

ANALYSIS="$MIMM_REPO/MUMC_pipeline/analysis"

echo "Cohort myelin-per-ROI flags (reading $RESULTDIR/*/analysis/roi_stats.csv)..."
python3 "$ANALYSIS/cohort_myelin_flags.py" "$RESULTDIR" || { echo "cohort analysis failed"; exit 1; }

echo ""
echo "Done. Results in: $RESULTDIR/cohort_analysis/"
