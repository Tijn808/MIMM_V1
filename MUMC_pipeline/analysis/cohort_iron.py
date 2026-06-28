#!/usr/bin/env python3
"""
Iron-channel analyses from the per-subject roi_stats.csv files:

1. Cross-method validation of MIMM's susceptibility channels against chi-sep,
   per ROI (cohort mean):
     - iron:   MIMM chi_iron  vs chi-sep chi_pos (chi+)
     - myelin: MIMM chi_myelin vs chi-sep chi_neg (chi-)
2. The systematic susceptibility offset between MIMM and chi-sep, via Bland-Altman
   of each channel (bias = the calibration offset).

Outputs to <results>/cohort_analysis/:
  cohort_iron_crossmethod.png    scatter + r, iron and myelin-chi
  cohort_chi_offset.png          Bland-Altman, iron and myelin-chi (offset = bias)

Usage:  python3 cohort_iron.py <results_dir | roi_stats.csv>
"""
import sys, os, glob
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory

path = sys.argv[1] if len(sys.argv) > 1 else \
    '/home/tijn-saes/Documents/Internship/ME_GRE/analysis/roi_stats.csv'
if os.path.isdir(path):
    files = sorted(glob.glob(os.path.join(path, '*', 'analysis', 'roi_stats.csv')))
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    df = df.groupby('ROI_name', as_index=False).mean(numeric_only=True)
    scope = f'cohort (n = {len(files)})'
    out_dir = os.path.join(path, 'cohort_analysis'); os.makedirs(out_dir, exist_ok=True)
else:
    df = pd.read_csv(path); scope = 'single subject'
    out_dir = os.path.dirname(path)

# (label, MIMM column, chi-sep column)
CHANNELS = [
    ('Iron  (chi+ )', 'chi_iron_atlas_mean',   'chi_pos_chisep_mean'),
    ('Myelin (chi-)', 'chi_myelin_atlas_mean', 'chi_neg_chisep_mean'),
]

# ---- Figure 1: cross-method scatter ----
fig, axes = plt.subplots(1, 2, figsize=(11, 5))
for ax, (label, mimm_c, cs_c) in zip(axes, CHANNELS):
    x = df[cs_c].values; y = df[mimm_c].values
    m = np.isfinite(x) & np.isfinite(y); x, y = x[m], y[m]
    ax.scatter(x, y, s=26, c='#d62728' if 'Iron' in label else '#1f77b4', alpha=0.75)
    r, p = stats.pearsonr(x, y)
    b, a = np.polyfit(x, y, 1); xs = np.linspace(x.min(), x.max(), 100)
    ax.plot(xs, b * xs + a, 'k--', lw=1.3)
    ax.text(0.04, 0.94, f'r = {r:.2f}\n' + ('p < 0.001' if p < 1e-3 else f'p = {p:.3f}'),
            transform=ax.transAxes, va='top', fontsize=10,
            bbox=dict(boxstyle='round', fc='white', ec='0.7'))
    ax.set_xlabel(f'chi-separation  {cs_c.split("_")[1]}  (ppm)')
    ax.set_ylabel(f'MIMM  {mimm_c.split("_")[1]}  (ppm)')
    ax.set_title(label); ax.grid(alpha=0.25)
fig.suptitle(f'MIMM vs chi-separation susceptibility, per JHU WM ROI  —  {scope}', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(os.path.join(out_dir, 'cohort_iron_crossmethod.png'), dpi=150)

# ---- Figure 2: Bland-Altman (offset) ----
fig, axes = plt.subplots(1, 2, figsize=(11, 5))
# paper-style (Sisman et al. 2025, Fig 4): blue points, red MEAN / +-1.96SD lines, each labelled
for ax, (label, mimm_c, cs_c) in zip(axes, CHANNELS):
    a1 = df[mimm_c].values; a2 = df[cs_c].values
    m = np.isfinite(a1) & np.isfinite(a2); a1, a2 = a1[m], a2[m]
    mean = (a1 + a2) / 2; diff = a1 - a2          # MIMM - chi-sep
    bias = diff.mean(); sd = diff.std(ddof=1)
    hi, lo = bias + 1.96 * sd, bias - 1.96 * sd
    ax.scatter(mean, diff, s=26, c='#1f77b4', alpha=0.75)
    ax.axhline(bias, color='#d62728', lw=1.8)
    ax.axhline(hi, color='#d62728', ls='--', lw=1.3)
    ax.axhline(lo, color='#d62728', ls='--', lw=1.3)
    tform = blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(0.99, bias, f'MEAN: {bias:+.4f}', transform=tform, va='bottom', ha='right',
            color='#d62728', fontsize=9)
    ax.text(0.99, hi, f'+1.96SD: {hi:+.4f}', transform=tform, va='bottom', ha='right',
            color='#d62728', fontsize=9)
    ax.text(0.99, lo, f'-1.96SD: {lo:+.4f}', transform=tform, va='top', ha='right',
            color='#d62728', fontsize=9)
    ax.set_xlabel('Mean of MIMM & chi-sep  (ppm)')
    ax.set_ylabel('MIMM - chi-sep  (ppm)')
    ax.set_title(label); ax.grid(alpha=0.25)
fig.suptitle(f'Systematic susceptibility offset, MIMM vs chi-separation  —  {scope}', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(os.path.join(out_dir, 'cohort_chi_offset.png'), dpi=150)

# ---- summary ----
for label, mimm_c, cs_c in CHANNELS:
    x = df[cs_c].values; y = df[mimm_c].values
    m = np.isfinite(x) & np.isfinite(y); x, y = x[m], y[m]
    r, _ = stats.pearsonr(x, y); diff = y - x; bias = diff.mean()
    mean_ba = (x + y) / 2
    ba_slope, _ = np.polyfit(mean_ba, diff, 1)
    print(f'{label}: r = {r:.3f},  bias = {bias:+.4f} ppm,  '
          f'BA slope = {ba_slope:+.3f} ppm/ppm  (n={len(x)} ROIs)')
print(f'saved figures to {out_dir}')
