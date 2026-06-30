#!/usr/bin/env python3
"""
Section 3.1 figure: MIMM agrees with chi-separation on BOTH source-separation
channels.

  Left  = myelin:  MIMM MVF (Basic + Atlas)  vs chi-separation |chi-|
  Right = iron:    MIMM chi_iron (Atlas)      vs chi-separation chi+

Each point is ONE SUBJECT: the whole-white-matter mean over all JHU ROIs,
voxel-weighted (n_voxels), so larger/more reliable tracts dominate and the small
noisy labels do not. Plotting one point per subject tests between-subject
agreement of the two methods within the cohort, rather than the tract-anatomy
correlation that a per-ROI plot is dominated by.

Usage:  python3 fig_crossmethod_channels.py <results_dir>
Output: <results_dir>/cohort_analysis/cohort_crossmethod_channels.png
"""
import sys, os, glob
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    sys.exit('usage: fig_crossmethod_channels.py <results_dir>')
results_dir = sys.argv[1]
files = sorted(glob.glob(os.path.join(results_dir, '*', 'analysis', 'roi_stats.csv')))
if not files:
    sys.exit('no roi_stats.csv found under <results_dir>/*/analysis/')
out_dir = os.path.join(results_dir, 'cohort_analysis'); os.makedirs(out_dir, exist_ok=True)


def wmean(vals, w):
    """Voxel-weighted mean over a subject's ROIs, ignoring NaN."""
    vals = np.asarray(vals, float); w = np.asarray(w, float)
    m = np.isfinite(vals) & np.isfinite(w) & (w > 0)
    return float(np.sum(vals[m] * w[m]) / np.sum(w[m])) if m.any() else np.nan


# one row per subject: whole-WM voxel-weighted mean of each channel
COLS = ['MVF_basic_mean', 'MVF_atlas_mean', 'chi_neg_chisep_mean',
        'chi_pos_chisep_mean', 'chi_iron_atlas_mean']
rows = []
for f in files:
    d = pd.read_csv(f)
    if 'n_voxels' not in d.columns:
        continue
    w = d['n_voxels'].values
    rows.append({c: wmean(d[c].values, w) for c in COLS if c in d.columns})
df = pd.DataFrame(rows)
n_sub = len(df)


def col(c):
    return df[c].values if c in df.columns else None


def scatter_r(ax, x, y, color, label):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y); x, y = x[m], y[m]
    ax.scatter(x, y, s=42, c=color, alpha=0.8, label=label)
    r, p = stats.pearsonr(x, y)
    b, a = np.polyfit(x, y, 1); xs = np.linspace(x.min(), x.max(), 100)
    ax.plot(xs, b * xs + a, '--', color=color, lw=1.2)
    return r, p


fig, (axM, axI) = plt.subplots(1, 2, figsize=(11.5, 5))

# --- myelin panel: MVF (Atlas) vs |chi-| ---
chineg = np.abs(col('chi_neg_chisep_mean'))
mvf = col('MVF_atlas_mean')
if mvf is None:
    sys.exit('MVF_atlas_mean column missing from roi_stats.csv')
r, p = scatter_r(axM, chineg, mvf, '#1f77b4', None)
axM.text(0.04, 0.96, f'r = {r:.2f}\n' + ('(p < 0.001)' if p < 1e-3 else f'(p = {p:.3f})'),
         transform=axM.transAxes, va='top',
         fontsize=10, bbox=dict(boxstyle='round', fc='white', ec='0.7'))
axM.set_xlabel('|chi-|  chi-separation  (ppm)')
axM.set_ylabel('MVF  MIMM  (fraction)')
axM.set_title('Myelin:  MIMM MVF vs chi-separation')
axM.grid(alpha=0.25)

# --- iron panel: MIMM chi_iron (Atlas) vs chi+ ---
chipos, chiiron = col('chi_pos_chisep_mean'), col('chi_iron_atlas_mean')
if chipos is None or chiiron is None:
    sys.exit('chi_pos_chisep_mean or chi_iron_atlas_mean column missing')
r, p = scatter_r(axI, chipos, chiiron, '#d62728', None)
axI.text(0.04, 0.96, f'r = {r:.2f}\n' + ('(p < 0.001)' if p < 1e-3 else f'(p = {p:.3f})'),
         transform=axI.transAxes, va='top',
         fontsize=10, bbox=dict(boxstyle='round', fc='white', ec='0.7'))
axI.set_xlabel('chi+  chi-separation  (ppm)')
axI.set_ylabel('chi_iron  MIMM  (ppm)')
axI.set_title('Iron:  MIMM chi_iron vs chi-separation'); axI.grid(alpha=0.25)

fig.suptitle(f'MIMM vs chi-separation, both channels, whole-WM mean per subject  (n={n_sub})', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = os.path.join(out_dir, 'cohort_crossmethod_channels.png')
fig.savefig(out, dpi=150)
print(f'saved: {out}  (n={n_sub} subjects, one point each)')
for lab, mc, cc in [('myelin MVF-atlas vs |chi-|', 'MVF_atlas_mean', 'chi_neg_chisep_mean'),
                    ('iron chi_iron vs chi+', 'chi_iron_atlas_mean', 'chi_pos_chisep_mean')]:
    x, y = col(cc), col(mc)
    if x is not None and y is not None:
        x = np.abs(x) if 'chi_neg' in cc else x
        m = np.isfinite(x) & np.isfinite(y)
        r, p = stats.pearsonr(x[m], y[m])
        print(f'  {lab}: r = {r:.3f}, p = {p:.3g}  (n={m.sum()})')
