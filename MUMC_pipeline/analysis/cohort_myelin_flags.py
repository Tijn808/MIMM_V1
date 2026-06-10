"""
cohort_myelin_flags.py - cohort myelin-per-ROI overview with outlier flagging.

No lesion masks needed (no ground truth yet): this pools every subject's
roi_stats.csv, looks at myelin (MVF, orientation-corrected) per ROI, and flags
subject/ROI combinations that deviate strongly from the cohort. A strongly LOW
value is a candidate for demyelination / lesion; a strongly HIGH value is more
likely a processing artefact. Confirm against lesion masks once they exist.

Reads:  RESULTDIR/<subject>/analysis/roi_stats.csv  (the MUMC layout)
Writes (into RESULTDIR/cohort_analysis/):
  cohort_mvf_zscore_heatmap.png  - subjects x ROIs heatmap of per-ROI MVF z-scores
  cohort_mvf_flags.csv           - every flagged subject/ROI combo, sorted by z
  cohort_mvf_flagcount.png       - number of flags per subject (lesion-burden proxy)
  cohort_mvf_matrix.csv          - the raw subject x ROI MVF table

Usage:
  python3 cohort_myelin_flags.py [RESULTDIR]      (default: $RESULTDIR or .)
  optional env: MIMM_METRIC (default MVF_atlas_mean), MIMM_ZTHR (default 2.0)
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

RESULTDIR = sys.argv[1] if len(sys.argv) > 1 else os.environ.get('RESULTDIR', '.')
METRIC    = os.environ.get('MIMM_METRIC', 'MVF_atlas_mean')
Z_THR     = float(os.environ.get('MIMM_ZTHR', '2.0'))
OUT       = os.path.join(RESULTDIR, 'cohort_analysis')
os.makedirs(OUT, exist_ok=True)

csvs = sorted(glob.glob(os.path.join(RESULTDIR, '*', 'analysis', 'roi_stats.csv')))
if not csvs:
    raise SystemExit(f'No roi_stats.csv found under {RESULTDIR}/*/analysis/')

# Subject x ROI matrix of the myelin metric.
cols = {}
for c in csvs:
    subj = os.path.basename(os.path.dirname(os.path.dirname(c)))
    df = pd.read_csv(c)
    if METRIC not in df.columns:
        print(f'  [skip] {subj}: no column {METRIC}')
        continue
    cols[subj] = df.set_index('ROI_name')[METRIC]
M = pd.DataFrame(cols).T                      # rows = subjects, cols = ROIs
M = M.loc[:, M.notna().any(axis=0)]           # drop all-NaN ROIs
M = M.sort_index()
print(f'Cohort: {M.shape[0]} subjects x {M.shape[1]} ROIs  (metric: {METRIC})')
M.to_csv(os.path.join(OUT, 'cohort_mvf_matrix.csv'))

# Per-ROI z-score across subjects.
mu = M.mean(axis=0)
sd = M.std(axis=0).replace(0, np.nan)
Z  = (M - mu) / sd

# Flag table.
flags = []
for s in Z.index:
    for r in Z.columns:
        z = Z.loc[s, r]
        if np.isfinite(z) and abs(z) >= Z_THR:
            flags.append({
                'subject': s, 'ROI': r,
                METRIC: round(float(M.loc[s, r]), 4),
                'cohort_mean': round(float(mu[r]), 4),
                'z': round(float(z), 2),
                'flag': 'LOW (demyelination?)' if z < 0 else 'HIGH (artefact?)',
            })
flags = pd.DataFrame(flags)
if len(flags):
    flags = flags.sort_values('z')
flags.to_csv(os.path.join(OUT, 'cohort_mvf_flags.csv'), index=False)
print(f'Flagged {len(flags)} subject/ROI combos at |z| >= {Z_THR}')
if len(flags):
    print('\nMost extreme LOW-myelin flags (lesion candidates):')
    print(flags[flags['z'] < 0].head(15).to_string(index=False))

# Per-subject flag count (low-myelin only = lesion-burden proxy).
low = flags[flags['z'] < 0] if len(flags) else flags
flagcount = low.groupby('subject').size().reindex(M.index, fill_value=0) if len(low) else pd.Series(0, index=M.index)

# ── Figure 1: subjects x ROIs z-score heatmap, flagged cells marked ──────────
fig, ax = plt.subplots(figsize=(max(12, M.shape[1] * 0.34), max(5, M.shape[0] * 0.45)),
                       facecolor='white')
vlim = max(3.0, np.nanmax(np.abs(Z.values)))
im = ax.imshow(Z.values, cmap='RdBu', vmin=-vlim, vmax=vlim, aspect='auto')
ax.set_xticks(range(M.shape[1])); ax.set_xticklabels(M.columns, rotation=90, fontsize=6)
ax.set_yticks(range(M.shape[0])); ax.set_yticklabels(M.index, fontsize=8)
# mark flagged cells
for j, r in enumerate(M.columns):
    for i, s in enumerate(M.index):
        z = Z.loc[s, r]
        if np.isfinite(z) and abs(z) >= Z_THR:
            ax.text(j, i, '×', ha='center', va='center', fontsize=7,
                    color='black', fontweight='bold')
ax.set_title(f'Cohort myelin (MVF atlas) per ROI — z-score vs cohort\n'
             f'blue = low (lesion candidate), red = high; × = |z| ≥ {Z_THR}',
             fontsize=11)
cb = fig.colorbar(im, ax=ax, fraction=0.015, pad=0.01)
cb.set_label('z-score (per ROI)', fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(OUT, 'cohort_mvf_zscore_heatmap.png'), dpi=150, bbox_inches='tight')
plt.close(); print('Saved: cohort_mvf_zscore_heatmap.png')

# ── Figure 2: per-subject low-myelin flag count ──────────────────────────────
fig, ax = plt.subplots(figsize=(10, 5), facecolor='white')
fc = flagcount.sort_values(ascending=False)
ax.bar(range(len(fc)), fc.values, color='#c0504d')
ax.set_xticks(range(len(fc))); ax.set_xticklabels(fc.index, rotation=90, fontsize=8)
ax.set_ylabel(f'# ROIs with low MVF (z ≤ -{Z_THR})')
ax.set_title('Low-myelin flags per subject (lesion-burden proxy, pre-ground-truth)', fontsize=11)
ax.grid(axis='y', ls='--', alpha=0.4)
plt.tight_layout()
plt.savefig(os.path.join(OUT, 'cohort_mvf_flagcount.png'), dpi=150, bbox_inches='tight')
plt.close(); print('Saved: cohort_mvf_flagcount.png')

print(f'\nOutputs in: {OUT}')
