#!/usr/bin/env python3
"""
All "doable now" cohort analyses from the per-subject roi_stats.csv files
(MIMM + chi-separation; no MWF / lesions / clinical scores needed).

Produces, in <results>/cohort_analysis/:
  cohort_crossmethod.png    - MVF vs |chi_neg| and MVF vs FA (cross-method validation)
  cohort_basic_vs_atlas.png - Basic-vs-Atlas MVF, and (Basic-Atlas) vs theta (orientation)
  cohort_patient_vs_hc.png  - per-subject WM-mean MVF, patients vs healthy controls

Usage:
  python3 cohort_doable.py <results_dir>
"""
import sys, os, glob
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    sys.exit('usage: cohort_doable.py <results_dir>')
results_dir = sys.argv[1]
out_dir = os.path.join(results_dir, 'cohort_analysis')
os.makedirs(out_dir, exist_ok=True)
HC = {'IMPROMYMS_001', 'IMPROMYMS_005', 'IMPROMYMS_006', 'IMPROMYMS_020'}

# ---- pool every subject's roi_stats.csv, keeping the subject id ----
frames = []
for f in sorted(glob.glob(os.path.join(results_dir, '*', 'analysis', 'roi_stats.csv'))):
    sid = os.path.basename(os.path.dirname(os.path.dirname(f)))
    d = pd.read_csv(f); d['subject'] = sid
    frames.append(d)
if not frames:
    sys.exit(f'No */analysis/roi_stats.csv under {results_dir}')
long = pd.concat(frames, ignore_index=True)
long['chi_neg_abs'] = long['chi_neg_chisep_mean'].abs()
nsub = long['subject'].nunique()
roi = long.groupby('ROI_name', as_index=False).mean(numeric_only=True)   # ROI means over cohort
print(f'{nsub} subjects, {roi.shape[0]} ROIs')


def scatter(ax, x, y, xlabel, ylabel, title):
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    ax.scatter(x, y, s=26, c='#1f77b4', alpha=0.75, edgecolors='none')
    r, p = stats.pearsonr(x, y)
    b, a = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 100)
    ax.plot(xs, b * xs + a, 'k--', lw=1.3)
    ax.text(0.04, 0.94, f'r = {r:.2f}\n' + ('p < 0.001' if p < 1e-3 else f'p = {p:.3f}'),
            transform=ax.transAxes, va='top', fontsize=10,
            bbox=dict(boxstyle='round', fc='white', ec='0.7'))
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_title(title); ax.grid(alpha=0.25)
    return r, p

# ============ 1. Cross-method validation ============
# Two correlation regimes, reported for an honest comparison with the MWF result:
#   - cohort-mean: average each ROI over subjects first (~50 points). Within-subject
#     noise is averaged out, so r is high but NOT comparable to a per-ROI r.
#   - per-ROI pooled: every subject x ROI is its own point (matches cohort_mwf.py),
#     so this r IS directly comparable to the MVF-vs-MWF r.
def pearson(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    return stats.pearsonr(x[m], y[m])[0], int(m.sum())

r_chi_atlas_pool, n_pool = pearson(long['chi_neg_abs'].values, long['MVF_atlas_mean'].values)
r_chi_basic_pool, _      = pearson(long['chi_neg_abs'].values, long['MVF_basic_mean'].values)
r_fa_pool, _             = pearson(long['FA_mean'].values,     long['MVF_atlas_mean'].values)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5))
# chi panel: both Basic and Atlas series on the cohort-mean ROIs
xb = roi['chi_neg_abs'].values
for yc, col, lab in [('MVF_basic_mean', '#ff7f0e', 'Basic'),
                     ('MVF_atlas_mean', '#1f77b4', 'Atlas')]:
    y = roi[yc].values; m = np.isfinite(xb) & np.isfinite(y)
    a1.scatter(xb[m], y[m], s=26, c=col, alpha=0.75, edgecolors='none', label=lab)
    b, a = np.polyfit(xb[m], y[m], 1); xs = np.linspace(xb[m].min(), xb[m].max(), 100)
    a1.plot(xs, b * xs + a, '--', color=col, lw=1.2)
r_chi, _       = pearson(xb, roi['MVF_atlas_mean'].values)   # cohort-mean ROIs
r_chi_basic, _ = pearson(xb, roi['MVF_basic_mean'].values)
a1.text(0.04, 0.94, f'Atlas r = {r_chi:.2f}\nBasic r = {r_chi_basic:.2f}\n'
        f'(cohort-mean ROIs)\nper-ROI pooled r = {r_chi_atlas_pool:.2f}',
        transform=a1.transAxes, va='top', fontsize=9,
        bbox=dict(boxstyle='round', fc='white', ec='0.7'))
a1.set_xlabel('|χ⁻| chi-separation (ppm)'); a1.set_ylabel('MVF')
a1.set_title('MIMM vs chi-separation'); a1.legend(loc='lower right'); a1.grid(alpha=0.25)
r_fa, _ = scatter(a2, roi['FA_mean'].values, roi['MVF_atlas_mean'].values,
                  'FA', 'MVF (Atlas)', 'MIMM vs FA')
fig.suptitle(f'Cross-method validation, per JHU WM ROI  (cohort n={nsub})', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(os.path.join(out_dir, 'cohort_crossmethod.png'), dpi=150)

# ============ 2. Basic vs Atlas + orientation ============
roi['overest'] = roi['MVF_basic_mean'] - roi['MVF_atlas_mean']
fig, (b1, b2) = plt.subplots(1, 2, figsize=(11, 5))
# Basic vs Atlas identity
b1.scatter(roi['MVF_atlas_mean'], roi['MVF_basic_mean'], s=26, c='#1f77b4', alpha=0.75)
lim = [min(roi['MVF_atlas_mean'].min(), roi['MVF_basic_mean'].min()),
       max(roi['MVF_atlas_mean'].max(), roi['MVF_basic_mean'].max())]
b1.plot(lim, lim, 'k--', lw=1.2, label='identity')
b1.set_xlabel('MVF Atlas'); b1.set_ylabel('MVF Basic')
b1.set_title('Basic vs Atlas'); b1.legend(); b1.grid(alpha=0.25)
md = roi['overest'].mean()
b1.text(0.04, 0.94, f'mean Basic−Atlas = {md:+.3f}\n({100*md/roi["MVF_atlas_mean"].mean():+.0f}%)',
        transform=b1.transAxes, va='top', fontsize=10,
        bbox=dict(boxstyle='round', fc='white', ec='0.7'))
# overestimation vs theta
r_theta, _ = scatter(b2, roi['theta_mean'].values, roi['overest'].values,
        'θ  fibre angle to B0 (deg)', 'MVF Basic − Atlas', 'Orientation-driven overestimation')
fig.suptitle(f'Orientation effect: Basic vs Atlas  (cohort n={nsub})', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(os.path.join(out_dir, 'cohort_basic_vs_atlas.png'), dpi=150)

# ============ 3. Patient vs HC ============
per_subj = long.groupby('subject', as_index=False)['MVF_atlas_mean'].mean()
per_subj['group'] = np.where(per_subj['subject'].isin(HC), 'HC', 'patient')
pat = per_subj[per_subj.group == 'patient']['MVF_atlas_mean'].values
hc  = per_subj[per_subj.group == 'HC']['MVF_atlas_mean'].values
fig, ax = plt.subplots(figsize=(5.5, 5))
ax.boxplot([pat, hc], labels=[f'patient\n(n={len(pat)})', f'HC\n(n={len(hc)})'], showmeans=True)
for i, vals in enumerate([pat, hc], 1):
    ax.scatter(np.random.normal(i, 0.05, len(vals)), vals, c='#1f77b4', s=24, alpha=0.7)
if len(hc) >= 2 and len(pat) >= 2:
    t, pt = stats.ttest_ind(pat, hc, equal_var=False)
    ax.text(0.04, 0.96, ('p < 0.001' if pt < 1e-3 else f't-test p = {pt:.3f}'),
            transform=ax.transAxes, va='top', fontsize=10,
            bbox=dict(boxstyle='round', fc='white', ec='0.7'))
ax.set_ylabel('MVF (Atlas), WM mean'); ax.set_title(f'Myelin: patients vs HC  (n={nsub})')
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'cohort_patient_vs_hc.png'), dpi=150)

# ---- summary ----
print(f'MVF(Atlas) vs |chi-neg|, cohort-mean ROIs : r = {r_chi:.3f}  (n={roi.shape[0]} ROIs)')
print(f'MVF(Basic) vs |chi-neg|, cohort-mean ROIs : r = {r_chi_basic:.3f}')
print(f'MVF(Atlas) vs |chi-neg|, per-ROI pooled   : r = {r_chi_atlas_pool:.3f}  (n={n_pool} subj*ROI)'
      f'   <-- compare with MVF-vs-MWF r')
print(f'MVF(Basic) vs |chi-neg|, per-ROI pooled   : r = {r_chi_basic_pool:.3f}')
print(f'MVF vs FA, cohort-mean ROIs : r = {r_fa:.3f}  (per-ROI pooled r = {r_fa_pool:.3f})')
print(f'Basic-Atlas mean overestimation: {md:+.4f} ({100*md/roi["MVF_atlas_mean"].mean():+.1f}%)')
print(f'(Basic-Atlas) vs theta:  r = {r_theta:.3f}   (orientation-driven)')
if len(hc) >= 2 and len(pat) >= 2:
    print(f'Patient {np.mean(pat):.3f} vs HC {np.mean(hc):.3f}  (t-test p={pt:.3g})')
print(f'figures saved to {out_dir}')
