#!/usr/bin/env python3
"""
Correlate per-subject MIMM myelin (MVF) with demographics (age, sex).

Per subject, the myelin index is the mean of the per-ROI MVF (Atlas) across the
JHU white-matter ROIs in that subject's analysis/roi_stats.csv. This is joined
to a demographics table (Participant Id, age, sex) and plotted:
  - MVF vs age  (Pearson r + regression), healthy controls vs patients marked
  - MVF by sex  (group comparison, t-test)

NOTE: the cohort mixes MS patients and healthy controls, so a raw MVF-vs-age
correlation confounds normal aging with disease. Interpret accordingly.

Usage:
  python3 correlate_demographics.py <demographics.csv> <results_dir> [out_dir]
"""
import sys, os, glob, re
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv) < 3:
    sys.exit('usage: correlate_demographics.py <demographics.csv> <results_dir> [out_dir]')
demo_csv, results_dir = sys.argv[1], sys.argv[2]
out_dir = sys.argv[3] if len(sys.argv) > 3 else os.path.join(results_dir, 'cohort_analysis')
os.makedirs(out_dir, exist_ok=True)

HC = {'IMPROMYMS_001', 'IMPROMYMS_005', 'IMPROMYMS_006', 'IMPROMYMS_020'}

# --- demographics: detect id / age / sex columns robustly ---
demo = pd.read_csv(demo_csv)
cols = {c.lower().strip(): c for c in demo.columns}
id_col  = next(c for k, c in cols.items() if 'participant' in k or k == 'id' or 'subject' in k)
age_col = next(c for k, c in cols.items() if 'age' in k)
sex_col = next(c for k, c in cols.items() if 'sex' in k or 'gender' in k)
demo = demo.rename(columns={id_col: 'id', age_col: 'age', sex_col: 'sex'})
demo['id'] = demo['id'].astype(str).str.strip()
demo['sex'] = demo['sex'].astype(str).str.strip().str.capitalize()

# --- per-subject myelin index from roi_stats.csv ---
rows = []
for f in sorted(glob.glob(os.path.join(results_dir, '*', 'analysis', 'roi_stats.csv'))):
    sid = os.path.basename(os.path.dirname(os.path.dirname(f)))
    r = pd.read_csv(f)
    rows.append({'id': sid,
                 'MVF_atlas': float(np.nanmean(r['MVF_atlas_mean'])),
                 'MVF_basic': float(np.nanmean(r['MVF_basic_mean'])),
                 'chi_pos':   float(np.nanmean(r['chi_pos_chisep_mean']))})  # iron (chi+)
mvf = pd.DataFrame(rows)
if mvf.empty:
    sys.exit(f'No roi_stats.csv found under {results_dir}')

df = mvf.merge(demo[['id', 'age', 'sex']], on='id', how='inner')
df['group'] = np.where(df['id'].isin(HC), 'HC', 'patient')
df.to_csv(os.path.join(out_dir, 'demographics_mvf.csv'), index=False)
print(f'Merged {len(df)} subjects (of {len(mvf)} with MIMM, {len(demo)} in demographics)')

# --- figure: MVF vs age, MVF by sex, iron (chi+) vs age ---
fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(15, 5))

# MVF vs age
for grp, col, mk in [('patient', '#1f77b4', 'o'), ('HC', '#2ca02c', 's')]:
    g = df[df['group'] == grp]
    axA.scatter(g['age'], g['MVF_atlas'], c=col, marker=mk, s=40, label=grp, alpha=0.8)
x, y = df['age'].values, df['MVF_atlas'].values
m = np.isfinite(x) & np.isfinite(y)
r, p = stats.pearsonr(x[m], y[m])
b, a = np.polyfit(x[m], y[m], 1)
xs = np.linspace(x[m].min(), x[m].max(), 100)
axA.plot(xs, b * xs + a, 'k--', lw=1.3)
axA.text(0.04, 0.06, f'r = {r:.2f}, ' + ('p < 0.001' if p < 1e-3 else f'p = {p:.3f}') +
         f'  (n = {m.sum()})', transform=axA.transAxes, fontsize=10,
         bbox=dict(boxstyle='round', fc='white', ec='0.7'))
axA.set_xlabel('Age (years)'); axA.set_ylabel('MVF (Atlas), WM mean')
axA.set_title('Myelin vs age'); axA.legend(); axA.grid(alpha=0.25)

# MVF by sex
groups, labels = [], []
for s in sorted(df['sex'].dropna().unique()):
    vals = df[df['sex'] == s]['MVF_atlas'].dropna().values
    if len(vals):
        groups.append(vals); labels.append(f'{s}\n(n={len(vals)})')
axB.boxplot(groups, labels=labels, showmeans=True)
for i, vals in enumerate(groups, 1):
    axB.scatter(np.random.normal(i, 0.05, len(vals)), vals, c='#1f77b4', s=22, alpha=0.7)
if len(groups) == 2:
    t, pt = stats.ttest_ind(groups[0], groups[1])
    axB.text(0.04, 0.94, ('p < 0.001' if pt < 1e-3 else f't-test p = {pt:.3f}'),
             transform=axB.transAxes, va='top', fontsize=10,
             bbox=dict(boxstyle='round', fc='white', ec='0.7'))
axB.set_ylabel('MVF (Atlas), WM mean'); axB.set_title('Myelin by sex'); axB.grid(alpha=0.25)

# iron (chi+) vs age — iron is expected to accumulate with age
for grp, col, mk in [('patient', '#d62728', 'o'), ('HC', '#2ca02c', 's')]:
    g = df[df['group'] == grp]
    axC.scatter(g['age'], g['chi_pos'], c=col, marker=mk, s=40, label=grp, alpha=0.8)
xi, yi = df['age'].values, df['chi_pos'].values
mi = np.isfinite(xi) & np.isfinite(yi)
ri, pi = stats.pearsonr(xi[mi], yi[mi])
bi, ai = np.polyfit(xi[mi], yi[mi], 1)
xs2 = np.linspace(xi[mi].min(), xi[mi].max(), 100)
axC.plot(xs2, bi * xs2 + ai, 'k--', lw=1.3)
axC.text(0.04, 0.94, f'r = {ri:.2f}, ' + ('p < 0.001' if pi < 1e-3 else f'p = {pi:.3f}'),
         transform=axC.transAxes, va='top', fontsize=10,
         bbox=dict(boxstyle='round', fc='white', ec='0.7'))
axC.set_xlabel('Age (years)'); axC.set_ylabel('chi+ (iron), WM mean (ppm)')
axC.set_title('Iron vs age'); axC.legend(); axC.grid(alpha=0.25)

fig.suptitle('MIMM/chi-sep vs demographics', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.join(out_dir, 'mvf_vs_demographics.png')
fig.savefig(out, dpi=150)
print(f'MVF vs age:  r = {r:.3f}, p = {p:.4g}')
print(f'Iron vs age: r = {ri:.3f}, p = {pi:.4g}')
print(f'saved: {out}')
print(f'table: {os.path.join(out_dir, "demographics_mvf.csv")}')
