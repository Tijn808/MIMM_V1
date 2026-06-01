# =========================================================================
# MUMC Pipeline — Cohort Configuration Template (Python)
# =========================================================================
# 1. Copy this file to cohort_paths.py (same directory)
# 2. Fill in COHORT_DIR below
# 3. cohort_paths.py is gitignored — never committed
#
# Expected layout:
#   COHORT_DIR/
#     sub-01/{qsm,mimm,chisep,atlas,dti,grase,analysis}/
#     sub-02/...
#     cohort/                  ← created automatically
#       analysis/              ← cohort_roi_stats.csv, susceptibility_bias.csv
#       figures/               ← cohort figures
# =========================================================================

import os
import glob

COHORT_DIR = '/path/to/cohort_root'   # parent of all sub-XX/ folders

# Discover subjects: all sub-*/ directories that have run extract_roi_stats.py
# (i.e. have analysis/roi_stats.csv). Sorted for reproducibility.
SUBJECT_DIRS = sorted(glob.glob(os.path.join(COHORT_DIR, 'sub-*')))

# Derived cohort output dirs — do not edit
COHORT_OUT      = os.path.join(COHORT_DIR, 'cohort')
COHORT_ANALYSIS = os.path.join(COHORT_OUT, 'analysis')
COHORT_FIG      = os.path.join(COHORT_OUT, 'figures')

for d in [COHORT_OUT, COHORT_ANALYSIS, COHORT_FIG]:
    os.makedirs(d, exist_ok=True)

if not SUBJECT_DIRS:
    raise SystemExit(f'No sub-*/ directories found in {COHORT_DIR}')
print(f'Found {len(SUBJECT_DIRS)} subjects: {[os.path.basename(s) for s in SUBJECT_DIRS]}')
