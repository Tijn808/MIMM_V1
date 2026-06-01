"""
Aggregate per-subject ROI statistics into cohort-level summaries.

Reads each subject's analysis/roi_stats.csv (written by extract_roi_stats.py),
stacks them, and computes between-subject mean, SD, and SEM per ROI per map.

Outputs (to COHORT_ANALYSIS/):
  cohort_roi_stats.csv       — one row per ROI, between-subject stats
  all_subjects_long.csv      — one row per (subject × ROI), for mixed-effects models

Between-subject SEM = SD_between / sqrt(n_subjects).  This is the correct error
bar for cohort figures where each point is a subject-mean (not a voxel).
"""

import numpy as np
import pandas as pd
import os

try:
    from cohort_paths import SUBJECT_DIRS, COHORT_ANALYSIS
except ImportError:
    raise SystemExit('Copy cohort_paths_template.py to cohort_paths.py and fill in COHORT_DIR.')

# ── Load all per-subject CSVs ─────────────────────────────────────────────────
print(f'Loading ROI stats from {len(SUBJECT_DIRS)} subjects...')
frames = []
missing = []

for subj_dir in SUBJECT_DIRS:
    csv = os.path.join(subj_dir, 'analysis', 'roi_stats.csv')
    if not os.path.exists(csv):
        missing.append(os.path.basename(subj_dir))
        continue
    df = pd.read_csv(csv)
    df.insert(0, 'subject', os.path.basename(subj_dir))
    frames.append(df)

if missing:
    print(f'  Warning: roi_stats.csv missing for: {missing}')
if not frames:
    raise SystemExit('No roi_stats.csv found across any subject.')

long_df = pd.concat(frames, ignore_index=True)
n_subjects = long_df['subject'].nunique()
print(f'  Loaded {len(long_df)} rows from {n_subjects} subjects.')

# ── Save long-format table (useful for mixed-effects models) ──────────────────
long_csv = os.path.join(COHORT_ANALYSIS, 'all_subjects_long.csv')
long_df.to_csv(long_csv, index=False, float_format='%.5f')
print(f'Saved: {long_csv}')

# ── Identify numeric columns that represent per-voxel stats ──────────────────
# These are the columns we want to average across subjects.
# Structure: {map}_{stat} where stat in (mean, sd, median, iqr, p05, p95, cv,
#            fa_weighted_mean) plus r_MVF_vs_chineg, p_MVF_vs_chineg, n_voxels.
# We aggregate only the _mean columns (the ROI mean per subject), and the
# r/p columns. n_voxels is kept as the average across subjects.
skip_cols = {'subject', 'ROI_index', 'ROI_name'}
mean_cols = [c for c in long_df.columns
             if c.endswith('_mean') and c not in skip_cols]
r_cols    = [c for c in long_df.columns
             if c.startswith('r_') or c.startswith('p_')]
vox_cols  = ['n_voxels'] if 'n_voxels' in long_df.columns else []

agg_cols  = mean_cols + r_cols + vox_cols

# ── Aggregate: mean, SD, SEM across subjects per ROI ─────────────────────────
groups = long_df.groupby(['ROI_index', 'ROI_name'], sort=False)

rows = []
for (roi_idx, roi_name), grp in groups:
    row = {'ROI_index': roi_idx, 'ROI_name': roi_name,
           'n_subjects': grp['subject'].nunique()}
    for col in agg_cols:
        if col not in grp.columns:
            continue
        vals = grp[col].dropna()
        n = len(vals)
        if n == 0:
            row[col]               = np.nan
            row[f'{col}_sd']       = np.nan
            row[f'{col}_sem']      = np.nan
        else:
            m = float(vals.mean())
            s = float(vals.std(ddof=1)) if n > 1 else np.nan
            row[col]               = m
            row[f'{col}_sd']       = s
            row[f'{col}_sem']      = s / np.sqrt(n) if n > 1 else np.nan
    rows.append(row)

cohort_df = pd.DataFrame(rows)

# Restore MVF-sorted order consistent with per-subject plots
if 'MVF_basic_mean' in cohort_df.columns:
    cohort_df = cohort_df.sort_values('MVF_basic_mean', ascending=True).reset_index(drop=True)

# Convenience: add chi_neg_abs (absolute value for plotting)
if 'chi_neg_chisep_mean' in cohort_df.columns:
    cohort_df['chi_neg_abs'] = cohort_df['chi_neg_chisep_mean'].abs()

out_csv = os.path.join(COHORT_ANALYSIS, 'cohort_roi_stats.csv')
cohort_df.to_csv(out_csv, index=False, float_format='%.5f')
print(f'Saved: {out_csv}  ({len(cohort_df)} ROIs, {n_subjects} subjects)')

# ── Terminal summary ──────────────────────────────────────────────────────────
if 'MVF_basic_mean' in cohort_df.columns:
    top5 = cohort_df.nlargest(5, 'MVF_basic_mean')[['ROI_name', 'MVF_basic_mean',
                                                       'MVF_basic_mean_sem', 'n_subjects']]
    print(f'\nTop 5 ROIs by cohort-mean MVF basic:')
    print(top5.to_string(index=False))
