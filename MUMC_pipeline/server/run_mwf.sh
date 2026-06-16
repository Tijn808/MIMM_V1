#!/bin/bash
################################################################################
# run_mwf.sh - launch the MWF cohort (020) via Gerald's parallel_execute.sh,
#              with NO dash-flags to type (the thin client mangles '-').
#
# Usage (from scripts/):
#   sh run_mwf.sh                       # auto-detect subjects needing MWF, 4 screens
#   sh run_mwf.sh 8                     # same, 8 parallel screens
#   sh run_mwf.sh 2 IMPROMYMS_003 IMPROMYMS_004   # just these two, 2 screens (test)
#
# A subject is "needs MWF" if it has nifti/grase.nii but no MWF.nii yet.
# Each screen launches one MATLAB -> needs that many MATLAB licenses + RAM.
################################################################################
cd "$( dirname "$0" )" || exit 1
source ./project.config 2>/dev/null

# Limit each MATLAB to ONE thread so many can run in parallel without
# oversubscribing the CPU. MATLAB does its math via Intel MKL, whose threads are
# NOT limited by OMP_NUM_THREADS -- MKL_NUM_THREADS is needed too (per Gerald).
# These exports propagate through parallel_execute -> screen -> 020 -> matlab.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

N="${1:-4}"
shift 2>/dev/null

if [ "$#" -gt 0 ]; then
    subs="$@"
else
    subs=""
    for g in "$RESULTDIR"/*/nifti/grase.nii*; do
        [ -f "$g" ] || continue
        s=$( basename "$( dirname "$( dirname "$g" )" )" )
        [ -f "$RESULTDIR/$s/MWF.nii" ] && continue      # already computed
        subs="$subs $s"
    done
    subs=$( echo "$subs" | tr ' ' '\n' | sort -u | tr '\n' ' ' )
fi

if [ -z "$( echo "$subs" | tr -d ' ' )" ]; then
    echo "No subjects need MWF (all have MWF.nii, or no grase.nii found)."
    exit 0
fi

echo "Parallel screens : $N"
echo "MWF subjects     :$subs"
echo "Launching parallel_execute.sh ..."
./parallel_execute.sh -p 020_MWF_calc.sh -s $subs -n "$N"
