#!/usr/bin/env python3
"""
Rank the JHU white-matter ROIs by how REPRESENTATIVE and RELIABLE they are, to
answer "which JHU regions should we focus on?". For each ROI three properties are
computed across the cohort:

  size       mean voxel count (bigger ROI -> more stable mean, less partial volume)
  stability  cross-subject coefficient of variation of MVF (lower = more reliable)
  represent  Pearson r between the ROI's per-subject MVF and that subject's
             whole-WM-mean MVF (higher = the ROI tracks global myelin, i.e. it is
             representative of the subject rather than idiosyncratic)

A combined score ranks ROIs; the top ones are a defensible "representative" subset.

Inputs:
  <results>/cohort_analysis/cohort_mvf_matrix.csv  (subjects x ROIs, from cohort_myelin_flags.py)
  <results>/*/analysis/roi_stats.csv               (for mean voxel counts)

Outputs to <results>/cohort_analysis/:
  representative_rois.csv          ranked table
  representative_rois.png          size vs stability, coloured by representativeness

Usage:  python3 cohort_representative_rois.py <results_dir>
"""
import sys, os, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    sys.exit('usage: cohort_representative_rois.py <results_dir>')
results_dir = sys.argv[1]
ca = os.path.join(results_dir, 'cohort_analysis')

M = pd.read_csv(os.path.join(ca, 'cohort_mvf_matrix.csv'), index_col=0)  # subjects x ROIs
M = M.loc[:, M.notna().sum(axis=0) >= 3]      # ROI must exist in >=3 subjects
wm_mean = M.mean(axis=1)                       # each subject's whole-WM-mean MVF

# mean voxel count per ROI from the per-subject roi_stats
vox = {}
for f in glob.glob(os.path.join(results_dir, '*', 'analysis', 'roi_stats.csv')):
    d = pd.read_csv(f)
    if 'ROI_name' in d and 'n_voxels' in d:
        for _, r in d.iterrows():
            vox.setdefault(r['ROI_name'], []).append(r['n_voxels'])
vox_mean = {k: float(np.mean(v)) for k, v in vox.items()}

rows = []
for roi in M.columns:
    col = M[roi].dropna()
    if len(col) < 3:
        continue
    mu, sd = col.mean(), col.std(ddof=1)
    cov = sd / mu if mu else np.nan
    common = col.index.intersection(wm_mean.index)
    represent = np.corrcoef(col.loc[common], wm_mean.loc[common])[0, 1] if len(common) >= 3 else np.nan
    rows.append({'ROI': roi, 'size_vox': vox_mean.get(roi, np.nan),
                 'mean_MVF': round(mu, 4), 'CoV': round(cov, 3),
                 'represent_r': round(represent, 3)})
T = pd.DataFrame(rows)

# Combined score: representativeness (high) + stability (low CoV) + size (large).
# Rank each property to [0,1], average. Robust to differing units.
T['rank_repr'] = T['represent_r'].rank(pct=True)
T['rank_stab'] = (-T['CoV']).rank(pct=True)
T['rank_size'] = T['size_vox'].rank(pct=True)
T['score'] = (T['rank_repr'] + T['rank_stab'] + T['rank_size']) / 3
T = T.sort_values('score', ascending=False).reset_index(drop=True)
T.to_csv(os.path.join(ca, 'representative_rois.csv'), index=False)

print(f'{len(T)} ROIs ranked. Top 12 most representative/reliable:')
print(T[['ROI', 'size_vox', 'CoV', 'represent_r', 'score']].head(12).to_string(index=False))
print('\nLeast representative (candidates to drop):')
print(T[['ROI', 'size_vox', 'CoV', 'represent_r', 'score']].tail(6).to_string(index=False))

# --- figure: size (log) vs CoV, colour = representativeness ---
fig, ax = plt.subplots(figsize=(9, 6))
sc = ax.scatter(T['size_vox'], T['CoV'], c=T['represent_r'], cmap='viridis',
                s=60, edgecolors='k', linewidths=0.4)
for _, r in T.head(12).iterrows():
    ax.annotate(r['ROI'], (r['size_vox'], r['CoV']), fontsize=6,
                xytext=(3, 3), textcoords='offset points')
ax.set_xscale('log')
ax.set_xlabel('ROI size (mean voxels, log scale)')
ax.set_ylabel('cross-subject CoV of MVF  (lower = more reliable)')
ax.set_title('JHU ROI representativeness\n(top-12 labelled; bottom-right = small & noisy)')
cb = fig.colorbar(sc); cb.set_label('representativeness r (ROI vs whole-WM MVF)')
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(os.path.join(ca, 'representative_rois.png'), dpi=150)
print(f'\nsaved: {os.path.join(ca, "representative_rois.csv")} and representative_rois.png')
