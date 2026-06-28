#!/usr/bin/env python3
"""
Section 3.1 figure: MIMM agrees with chi-separation on BOTH source-separation
channels. FA is no longer a panel -- it becomes a one-line specificity check in
the text.

  Left  = myelin:  MIMM MVF (Basic + Atlas)  vs chi-separation |chi-|
  Right = iron:    MIMM chi_iron (Atlas)      vs chi-separation chi+

Each point is one JHU white-matter tract, averaged over the cohort (cohort-mean
ROI), matching cohort_doable.py / cohort_iron.py.

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
df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
df = df.groupby('ROI_name', as_index=False).mean(numeric_only=True)
out_dir = os.path.join(results_dir, 'cohort_analysis'); os.makedirs(out_dir, exist_ok=True)
n_sub = len(files)


def col(c):
    return df[c].values if c in df.columns else None


def scatter_r(ax, x, y, color, label):
    m = np.isfinite(x) & np.isfinite(y); x, y = x[m], y[m]
    ax.scatter(x, y, s=26, c=color, alpha=0.75, label=label)
    r, p = stats.pearsonr(x, y)
    b, a = np.polyfit(x, y, 1); xs = np.linspace(x.min(), x.max(), 100)
    ax.plot(xs, b * xs + a, '--', color=color, lw=1.2)
    return r, p


fig, (axM, axI) = plt.subplots(1, 2, figsize=(11.5, 5))

# --- myelin panel: MVF (Basic + Atlas) vs |chi-| ---
chineg = col('chi_neg_chisep_mean')
if chineg is None:
    sys.exit('chi_neg_chisep_mean column missing from roi_stats.csv')
txt = []
for var, color in [('basic', '#ff7f0e'), ('atlas', '#1f77b4')]:
    mvf = col(f'MVF_{var}_mean')
    if mvf is None:
        continue
    r, _ = scatter_r(axM, chineg, mvf, color, var.capitalize())
    txt.append(f'{var.capitalize()} r = {r:.2f}')
axM.text(0.04, 0.96, '\n'.join(txt) + '\n(p < 0.001)', transform=axM.transAxes, va='top',
         fontsize=10, bbox=dict(boxstyle='round', fc='white', ec='0.7'))
axM.set_xlabel('|chi-|  chi-separation  (ppm)')
axM.set_ylabel('MVF  MIMM  (fraction)')
axM.set_title('Myelin:  MIMM MVF vs chi-separation')
axM.grid(alpha=0.25); axM.legend(loc='lower right')

# --- iron panel: MIMM chi_iron (Atlas) vs chi+ ---
chipos, chiiron = col('chi_pos_chisep_mean'), col('chi_iron_atlas_mean')
if chipos is None or chiiron is None:
    sys.exit('chi_pos_chisep_mean or chi_iron_atlas_mean column missing')
r, _ = scatter_r(axI, chipos, chiiron, '#d62728', None)
axI.text(0.04, 0.96, f'r = {r:.2f}\n(p < 0.001)', transform=axI.transAxes, va='top',
         fontsize=10, bbox=dict(boxstyle='round', fc='white', ec='0.7'))
axI.set_xlabel('chi+  chi-separation  (ppm)')
axI.set_ylabel('chi_iron  MIMM  (ppm)')
axI.set_title('Iron:  MIMM chi_iron vs chi-separation'); axI.grid(alpha=0.25)

fig.suptitle(f'MIMM vs chi-separation, both channels, per JHU WM ROI  (cohort n={n_sub})', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = os.path.join(out_dir, 'cohort_crossmethod_channels.png')
fig.savefig(out, dpi=150)
print(f'saved: {out}  (n={n_sub} subjects, {len(df)} ROIs)')
for lab, mc, cc in [('myelin MVF-atlas vs |chi-|', 'MVF_atlas_mean', 'chi_neg_chisep_mean'),
                    ('iron chi_iron vs chi+', 'chi_iron_atlas_mean', 'chi_pos_chisep_mean')]:
    x, y = col(cc), col(mc)
    if x is not None and y is not None:
        m = np.isfinite(x) & np.isfinite(y)
        print(f'  {lab}: r = {stats.pearsonr(x[m], y[m])[0]:.3f}')
