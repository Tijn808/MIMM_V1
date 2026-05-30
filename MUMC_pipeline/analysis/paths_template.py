# =========================================================================
# MUMC Pipeline — Machine Configuration Template (Python)
# =========================================================================
# 1. Copy this file to paths.py (same directory)
# 2. Fill in the three paths below
# 3. paths.py is gitignored — never committed
# =========================================================================

import os

MIMM_ROOT  = '/path/to/MIMM'
INPUT_DIR  = '/path/to/raw/ME_GRE_data'
OUTPUT_DIR = '/path/to/output/subject_id'   # created automatically if missing

# Derived subdirs — do not edit
QSM_DIR      = os.path.join(OUTPUT_DIR, 'qsm')
MIMM_DIR     = os.path.join(OUTPUT_DIR, 'mimm')
CHISEP_DIR   = os.path.join(OUTPUT_DIR, 'chisep')
ATLAS_DIR    = os.path.join(OUTPUT_DIR, 'atlas')
FIG_DIR      = os.path.join(OUTPUT_DIR, 'figures')
ANALYSIS_DIR = os.path.join(OUTPUT_DIR, 'analysis')

for d in [QSM_DIR, MIMM_DIR, CHISEP_DIR, ATLAS_DIR, FIG_DIR, ANALYSIS_DIR]:
    os.makedirs(d, exist_ok=True)
