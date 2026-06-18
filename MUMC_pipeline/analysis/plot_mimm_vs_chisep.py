#!/usr/bin/env python3
"""
Depict MIMM myelin (MVF) against chi-separation myelin (|chi_neg|), per ROI.

Both are myelin-sensitive but independent: MIMM matches the mGRE magnitude decay
+ QSM to a fibre-model dictionary; chi-separation splits susceptibility into a
diamagnetic (myelin) component chi_neg. Agreement across white-matter ROIs is a
cross-method validation of MIMM.

Usage:
  python3 plot_mimm_vs_chisep.py [roi_stats.csv | RESULTDIR] [out.png]
- a single roi_stats.csv  -> one-subject scatter
- a results/ directory    -> pools every <subj>/analysis/roi_stats.csv and
                             plots the cohort (one point = ROI mean over subjects)
Defaults to the local test-subject CSV.
"""
import sys, os, glob
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

path = sys.argv[1] if len(sys.argv) > 1 else \
    '/home/tijn-saes/Documents/Internship/ME_GRE/analysis/roi_stats.csv'

if os.path.isdir(path):
    files = sorted(glob.glob(os.path.join(path, '*', 'analysis', 'roi_stats.csv')))
    if not files:
        sys.exit(f'No */analysis/roi_stats.csv under {path}')
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    df = df.groupby('ROI_name', as_index=False).mean(numeric_only=True)
    scope = f'cohort (n = {len(files)} subjects)'
    default_out = os.path.join(path, 'cohort_analysis', 'mimm_vs_chisep_cohort.png')
    os.makedirs(os.path.dirname(default_out), exist_ok=True)
else:
    df = pd.read_csv(path)
    scope = 'single subject'
    default_out = os.path.join(os.path.dirname(path), 'mimm_vs_chisep.png')

out = sys.argv[2] if len(sys.argv) > 2 else default_out
df['chi_neg_abs'] = df['chi_neg_chisep_mean'].abs()   # diamagnetic magnitude (ppm)

panels = [('MVF_basic_mean', 'MIMM Basic'),
          ('MVF_atlas_mean', 'MIMM Atlas')]

fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
for ax, (ycol, title) in zip(axes, panels):
    x = df['chi_neg_abs'].values
    y = df[ycol].values
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]

    ax.scatter(x, y, s=28, c='#1f77b4', alpha=0.75, edgecolors='none')

    # regression + Pearson r
    r, p = stats.pearsonr(x, y)
    b, a = np.polyfit(x, y, 1)                       # y = b*x + a
    xs = np.linspace(x.min(), x.max(), 100)
    ax.plot(xs, b * xs + a, 'k--', lw=1.4)
    ax.text(0.04, 0.94,
            f'r = {r:.2f}\n' + ('p < 0.001' if p < 1e-3 else f'p = {p:.3f}') +
            f'\nn = {len(x)} ROIs',
            transform=ax.transAxes, va='top', fontsize=10,
            bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.9))

    # label the 3 highest-myelin ROIs for orientation (placed to the left so
    # they don't run off the right edge)
    names = df.loc[m]['ROI_name'].values
    for i in np.argsort(y)[-3:]:
        ax.annotate(names[i][:22], (x[i], y[i]), fontsize=7, alpha=0.8,
                    ha='right', xytext=(-6, -2), textcoords='offset points')

    ax.set_xlabel('|χ⁻|  chi-separation  (ppm)')
    ax.set_title(title)
    ax.grid(alpha=0.25)
axes[0].set_ylabel('MVF  MIMM  (fraction)')
fig.suptitle(f'MIMM myelin vs chi-separation myelin, per JHU WM ROI  —  {scope}',
             fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(out, dpi=150)
print(f'r(Basic)  = {stats.pearsonr(df["chi_neg_abs"], df["MVF_basic_mean"])[0]:.3f}')
print(f'r(Atlas)  = {stats.pearsonr(df["chi_neg_abs"], df["MVF_atlas_mean"])[0]:.3f}')
print(f'saved: {out}')
