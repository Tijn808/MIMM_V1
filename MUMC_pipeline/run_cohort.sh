#!/bin/bash
# ============================================================================
# MUMC MIMM cohort pipeline runner
# ============================================================================
# Runs the full per-subject pipeline across all sub-*/ directories, then
# runs the cohort-level aggregation and figures.
#
# Usage:
#   bash run_cohort.sh <cohort_dir> <mimm_root> <chisep_dir> [options]
#
# Required arguments:
#   cohort_dir   Parent directory containing sub-01/, sub-02/, ... folders
#   mimm_root    Path to this MIMM repository root
#   chisep_dir   Path to the Chi-separation toolbox
#
# Options:
#   --steps  STEPS   Comma-separated list of steps to run (default: all)
#                    matlab_pre  — prepare_mgre + QSM  (writes mag_e1; run first)
#                    register    — FSL atlas registration (needs mag_e1)
#                    dti         — topup/eddy/dtifit     (needs mag_e1, optional)
#                    matlab_post — chi-sep + MIMM        (needs register outputs)
#                    python      — all per-subject analysis scripts
#                    cohort      — aggregate + cohort figures
#                    matlab      — shorthand for matlab_pre + matlab_post
#   --subjects PATT  Glob pattern for subject dirs (default: sub-*)
#   --dry-run        Print commands without executing
#   --skip-done      Skip a subject's step if its expected output exists
#
# Examples:
#   bash run_cohort.sh /data/MUMC /home/user/MIMM /data/chisep
#   bash run_cohort.sh /data/MUMC /home/user/MIMM /data/chisep --steps matlab,python
#   bash run_cohort.sh /data/MUMC /home/user/MIMM /data/chisep --steps python,cohort
#   bash run_cohort.sh /data/MUMC /home/user/MIMM /data/chisep --subjects sub-0{1,2,3} --skip-done
# ============================================================================

set -e

# ── Parse required arguments ─────────────────────────────────────────────────
if [ $# -lt 3 ]; then
    echo "Usage: $0 <cohort_dir> <mimm_root> <chisep_dir> [options]"
    exit 1
fi

COHORT_DIR="$1"; MIMM_ROOT="$2"; CHISEP_DIR="$3"; shift 3

# ── Parse options ─────────────────────────────────────────────────────────────
# Correct order: matlab_pre (prepare+qsm) → register → dti → matlab_post (chisep+mimm)
# register and dti both need qsm/mag_e1.nii.gz which is written by the qsm step.
STEPS="matlab_pre,register,dti,matlab_post,python,cohort"
SUBJ_PATTERN="sub-*"
DRY_RUN=false
SKIP_DONE=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --steps)    STEPS="$2";       shift 2 ;;
        --subjects) SUBJ_PATTERN="$2"; shift 2 ;;
        --dry-run)  DRY_RUN=true;     shift   ;;
        --skip-done) SKIP_DONE=true;  shift   ;;
        *) echo "Unknown option: $1"; exit 1  ;;
    esac
done

# ── FSL path ──────────────────────────────────────────────────────────────────
FSL="${FSLDIR:-/home/tijn-saes/fsl}/bin"

# ── Python and analysis dir ───────────────────────────────────────────────────
ANALYSIS_DIR="$MIMM_ROOT/MUMC_pipeline/analysis"
PIPELINE_DIR="$MIMM_ROOT/MUMC_pipeline"

# ── Helper functions ──────────────────────────────────────────────────────────
run_cmd() {
    echo "  >> $*"
    if [ "$DRY_RUN" = false ]; then
        eval "$@"
    fi
}

skip_if_done() {
    local output="$1"; local label="$2"
    if [ "$SKIP_DONE" = true ] && [ -e "$output" ]; then
        echo "  [skip] $label (output exists: $output)"
        return 0
    fi
    return 1
}

# 'matlab' alone is a shorthand that enables both matlab_pre and matlab_post.
step_enabled() {
    local name="$1"
    if [[ "$name" == "matlab_pre" || "$name" == "matlab_post" ]]; then
        echo ",$STEPS," | grep -qE ",(${name}|matlab),"
    else
        echo ",$STEPS," | grep -q ",$name,"
    fi
}

log() { echo ""; echo "=== $* ==="; }

# ── Discover subjects ─────────────────────────────────────────────────────────
SUBJECTS=( "$COHORT_DIR"/$SUBJ_PATTERN )
if [ ${#SUBJECTS[@]} -eq 0 ] || [ ! -d "${SUBJECTS[0]}" ]; then
    echo "No subject directories found matching: $COHORT_DIR/$SUBJ_PATTERN"
    exit 1
fi
echo "Found ${#SUBJECTS[@]} subjects: ${SUBJECTS[*]##*/}"

# ── Per-subject loop ──────────────────────────────────────────────────────────
for SUBJ_DIR in "${SUBJECTS[@]}"; do
    SUBJ_ID=$(basename "$SUBJ_DIR")
    log "Subject: $SUBJ_ID"

    # ── Step 1: MATLAB pre — prepare_mgre + QSM ──────────────────────────────
    # MUST run first: writes qsm/mag_e1.nii.gz needed by register and dti.
    if step_enabled matlab_pre; then
        QSM_OUT="$SUBJ_DIR/qsm/mag_e1.nii.gz"
        if ! skip_if_done "$QSM_OUT" "matlab_pre (prepare+qsm)"; then
            echo "  [matlab_pre] $SUBJ_ID"
            run_cmd matlab -batch \
                "addpath('$PIPELINE_DIR'); \
                 run_subject('$SUBJ_DIR', '$MIMM_ROOT', '$CHISEP_DIR', \
                             {'prepare_mgre','qsm'})"
        fi
    fi

    # ── Step 2: Atlas registration (FSL) ─────────────────────────────────────
    # Needs qsm/mag_e1.nii.gz — must follow matlab_pre.
    if step_enabled register; then
        ATLAS_OUT="$SUBJ_DIR/atlas/theta_atlas.nii.gz"
        if ! skip_if_done "$ATLAS_OUT" "register_atlas"; then
            echo "  [register] $SUBJ_ID"
            run_cmd bash "$PIPELINE_DIR/atlas/register_atlas.sh" "$SUBJ_DIR"
        fi
    fi

    # ── Step 3: DTI preprocessing (optional — only if dti.nii.gz exists) ─────
    # Needs qsm/mag_e1.nii.gz — must follow matlab_pre.
    if step_enabled dti; then
        DTI_INPUT="$SUBJ_DIR/dti/dti.nii.gz"
        DTI_OUT="$SUBJ_DIR/dti/FA.nii.gz"
        if [ -f "$DTI_INPUT" ]; then
            if ! skip_if_done "$DTI_OUT" "preprocess_dti"; then
                echo "  [dti] $SUBJ_ID"
                REVERSE_B0="$SUBJ_DIR/dti/reverse_b0.nii.gz"
                if [ -f "$REVERSE_B0" ]; then
                    run_cmd bash "$PIPELINE_DIR/preprocess/preprocess_dti.sh" \
                        "$SUBJ_DIR" "$REVERSE_B0"
                else
                    run_cmd bash "$PIPELINE_DIR/preprocess/preprocess_dti.sh" \
                        "$SUBJ_DIR"
                fi
            fi
        else
            echo "  [skip dti] dti/dti.nii.gz not found — DTI step skipped"
        fi
    fi

    # ── Step 4: MATLAB post — chi-sep + MIMM ─────────────────────────────────
    # Needs atlas outputs from register and QSM outputs from matlab_pre.
    if step_enabled matlab_post; then
        MIMM_OUT="$SUBJ_DIR/mimm/MVF_basic.nii.gz"
        if ! skip_if_done "$MIMM_OUT" "matlab_post (chisep+mimm)"; then
            echo "  [matlab_post] $SUBJ_ID"
            run_cmd matlab -batch \
                "addpath('$PIPELINE_DIR'); \
                 run_subject('$SUBJ_DIR', '$MIMM_ROOT', '$CHISEP_DIR', \
                             {'chisep','mimm'})"
        fi
    fi

    # ── Step 4: Python analysis ───────────────────────────────────────────────
    if step_enabled python; then
        echo "  [python] $SUBJ_ID"
        export MIMM_OUTPUT_DIR="$SUBJ_DIR"

        ROI_CSV="$SUBJ_DIR/analysis/roi_stats.csv"
        if ! skip_if_done "$ROI_CSV" "extract_roi_stats"; then
            run_cmd python3 "$ANALYSIS_DIR/extract_roi_stats.py"
        fi

        run_cmd python3 "$ANALYSIS_DIR/generate_figures.py"
        run_cmd python3 "$ANALYSIS_DIR/plot_roi_stats.py"
        run_cmd python3 "$ANALYSIS_DIR/analyse_overestimation.py"
        run_cmd python3 "$ANALYSIS_DIR/plot_mwf.py"   # exits cleanly if no GRASE

        unset MIMM_OUTPUT_DIR
    fi

    echo "  [done] $SUBJ_ID"
done

# ── Cohort aggregation + figures ──────────────────────────────────────────────
if step_enabled cohort; then
    log "Cohort aggregation"

    # Write cohort_paths.py into a tmpdir and add it to PYTHONPATH.
    # This lets the cohort scripts import it without touching the user's
    # permanent cohort_paths.py (which may have custom settings).
    TMP_DIR=$(mktemp -d /tmp/mimm_cohort_XXXXXX)
    cat > "$TMP_DIR/cohort_paths.py" << PYEOF
import os, glob
COHORT_DIR      = '$COHORT_DIR'
SUBJECT_DIRS    = sorted(glob.glob(os.path.join(COHORT_DIR, 'sub-*')))
COHORT_OUT      = os.path.join(COHORT_DIR, 'cohort')
COHORT_ANALYSIS = os.path.join(COHORT_OUT, 'analysis')
COHORT_FIG      = os.path.join(COHORT_OUT, 'figures')
for d in [COHORT_OUT, COHORT_ANALYSIS, COHORT_FIG]:
    os.makedirs(d, exist_ok=True)
if not SUBJECT_DIRS:
    raise SystemExit(f'No sub-*/ directories found in {COHORT_DIR}')
print(f'Found {len(SUBJECT_DIRS)} subjects: {[os.path.basename(s) for s in SUBJECT_DIRS]}')
PYEOF

    export PYTHONPATH="$TMP_DIR:$ANALYSIS_DIR:${PYTHONPATH:-}"

    run_cmd python3 "$ANALYSIS_DIR/aggregate_cohort.py"
    run_cmd python3 "$ANALYSIS_DIR/plot_cohort_roi.py"
    run_cmd python3 "$ANALYSIS_DIR/plot_cohort_overestimation.py"

    unset PYTHONPATH
    rm -rf "$TMP_DIR"
    log "Cohort complete — outputs in $COHORT_DIR/cohort/"
fi

echo ""
echo "Pipeline complete."
echo "  Per-subject figures: $COHORT_DIR/sub-XX/figures/"
echo "  Cohort figures:      $COHORT_DIR/cohort/figures/"
