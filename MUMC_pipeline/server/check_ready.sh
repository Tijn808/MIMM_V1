#!/bin/bash
################################################################################
# check_ready.sh <subjectName>
#
# Pre-flight check for 030_MIMM.sh. Run it on the radstation BEFORE the real
# pipeline: it verifies the tools, the repo, the FSL atlases and the subject's
# input files are all in place, and that the results folder is writable.
# Reports PASS/FAIL per item and exits 1 if anything is missing, so the first
# real run does not fail 30 minutes in.
#
# Usage:  sh check_ready.sh <subjectName>      (run from the scripts/ folder)
################################################################################

source $( dirname "$0" )/functions.source 2>/dev/null
source $( dirname "$0" )/project.config   2>/dev/null

subjectName="$1"
: "${MIMM_REPO:=$SCRIPTDIR/MIMM_V1}"
subjectDir="${RESULTDIR:-.}/$subjectName"

fails=0
ok(){ echo "  [ ok ] $1"; }
no(){ echo "  [FAIL] $1"; fails=$((fails+1)); }

echo "=== MIMM readiness check: ${subjectName:-<no subject given>} ==="
[ -n "$subjectName" ] || no "no subjectName given (usage: sh check_ready.sh <subjectName>)"

echo "-- tools --"
command -v "${MATLAB_BIN:-matlab}" >/dev/null 2>&1 && ok "MATLAB (${MATLAB_BIN:-matlab})" || no "MATLAB not found - set MATLAB_BIN"
{ [ -n "$FSLDIR" ] && [ -x "$FSLDIR/bin/flirt" ]; } && ok "FSL ($FSLDIR)" || no "FSL not found - set FSLDIR"
command -v python3 >/dev/null 2>&1 && ok "python3" || no "python3 not found"
python3 -c "import nibabel,numpy,pandas,scipy,matplotlib" 2>/dev/null && ok "python modules" || no "python modules missing (nibabel numpy pandas scipy matplotlib)"

echo "-- chi-separation toolbox --"
{ [ -n "$CHISEP_DIR" ] && [ -d "$CHISEP_DIR" ]; } && ok "chi-sep toolbox dir ($CHISEP_DIR)" || no "CHISEP_DIR unset or not a directory"
# Find a ROMEO binary anywhere under it (handles nested copies) and check it runs
romeo_ok=no; romeo_found=no
for r in $(find "$CHISEP_DIR" -name romeo -path "*mritools*" -type f 2>/dev/null); do
    romeo_found=yes
    if "$r" --help >/dev/null 2>&1; then romeo_ok=yes; break; fi
done
[ "$romeo_found" = yes ] && ok "ROMEO binary present" || no "no ROMEO binary under CHISEP_DIR"
[ "$romeo_ok" = yes ]    && ok "a ROMEO binary runs on this machine" || no "ROMEO found but none ran (glibc mismatch?)"

echo "-- MIMM repo + dictionary --"
[ -f "$MIMM_REPO/MIMM_set_path.m" ] && ok "MIMM repo ($MIMM_REPO)" || no "MIMM repo not found at MIMM_REPO"
ls "$MIMM_REPO"/Dictionary/MIMM_dictionary_stochastic*.mat >/dev/null 2>&1 && ok "dictionary present" || no "dictionary .mat missing (did the clone include it?)"

echo "-- FSL atlases --"
[ -f "$FSLDIR/data/standard/FSL_HCP1065_FA_1mm.nii.gz" ] && ok "HCP1065 FA atlas" || no "HCP1065 atlas missing - FSL too old?"
[ -f "$FSLDIR/data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz" ] && ok "JHU label atlas" || no "JHU atlas missing in FSL"

echo "-- subject inputs (MUMC 010 layout: nifti/gremag.nii + grepha.nii) --"
[ -d "$subjectDir" ] && ok "subject dir: $subjectDir" || no "subject dir not found: $subjectDir"
{ ls "$subjectDir"/nifti/gremag.nii* >/dev/null 2>&1; } && ok "ME-GRE magnitude (nifti/gremag.nii)" || no "nifti/gremag.nii not found"
{ ls "$subjectDir"/nifti/grepha.nii* >/dev/null 2>&1; } && ok "ME-GRE phase (nifti/grepha.nii)" || no "nifti/grepha.nii not found"
echo "   (actual GRE files in nifti/:)"
ls "$subjectDir"/nifti 2>/dev/null | grep -iE "gremag|grepha|ME_GRE" | sed 's/^/     /' || echo "     (none matched)"

echo "-- write permission --"
if touch "$subjectDir/.mimm_writetest" 2>/dev/null; then ok "results dir writable"; rm -f "$subjectDir/.mimm_writetest"; else no "cannot write to $subjectDir"; fi

echo
if [ "$fails" -eq 0 ]; then
    echo "ALL CHECKS PASSED - ready to run:  sh 030_MIMM.sh $subjectName"
    exit 0
else
    echo "$fails CHECK(S) FAILED - fix the above before running 030_MIMM.sh"
    exit 1
fi
