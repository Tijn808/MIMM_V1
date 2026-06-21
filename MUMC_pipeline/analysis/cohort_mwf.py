#!/usr/bin/env python3
"""
MIMM MVF vs T2-GRASE MWF, per ROI, pooled across whoever has a registered MWF map
(grase/MWF.nii.gz, written by run_mwf_compare.sh). This is the only validation
against a reference that shares NO input data with MIMM (MWF is T2 relaxometry).

Runs on however many subjects are done -- no need for the full cohort.

Outputs to <results>/cohort_analysis/:
  cohort_mvf_vs_mwf.png    scatter (MVF vs MWF, r + regression) + Bland-Altman

Usage:  python3 cohort_mwf.py <results_dir>
"""
import sys, os, glob
import numpy as np
import nibabel as nib
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    sys.exit('usage: cohort_mwf.py <results_dir>')
results_dir = sys.argv[1]
out_dir = os.path.join(results_dir, 'cohort_analysis'); os.makedirs(out_dir, exist_ok=True)


def load(p):
    return np.asarray(nib.load(p).dataobj) if os.path.exists(p) else None

pts = []   # (subject_index, MVF, MWF) per ROI
subjects = []
for mp in sorted(glob.glob(os.path.join(results_dir, '*', 'grase', 'MWF.nii.gz'))):
    d = os.path.dirname(os.path.dirname(mp)); sid = os.path.basename(d)
    mwf = load(mp)
    mvf = load(os.path.join(d, 'mimm', 'MVF_Atlas.nii.gz'))
    labels = load(os.path.join(d, 'atlas', 'JHU_labels_subj.nii.gz'))
    brain = load(os.path.join(d, 'qsm', 'brain_mask.nii.gz'))
    if any(v is None for v in (mwf, mvf, labels, brain)):
        print(f'[skip] {sid} (missing input)'); continue
    mwf = np.clip(mwf.astype(float), 0, 0.5); mvf = mvf.astype(float)
    labels = labels.astype(int); brain = brain > 0
    si = len(subjects); subjects.append(sid); nroi = 0
    for roi in range(1, 51):
        m = (labels == roi) & brain & np.isfinite(mvf) & np.isfinite(mwf) & (mwf > 0)
        if m.sum() >= 20:
            pts.append((si, float(mvf[m].mean()), float(mwf[m].mean()))); nroi += 1
    print(f'{sid}: {nroi} ROIs')

if not pts:
    sys.exit('No subjects with grase/MWF.nii.gz found. Run run_mwf_compare.sh first.')
pts = np.array(pts); si = pts[:, 0].astype(int); mvf = pts[:, 1]; mwf = pts[:, 2]
n_sub = len(subjects)
print(f'\n{n_sub} subjects, {len(pts)} ROI points.')

fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 5.2))

# --- scatter: MVF vs MWF (do they track?) ---
sc = axA.scatter(mwf, mvf, c=si, cmap='tab20', s=18, alpha=0.7)
r, p = stats.pearsonr(mwf, mvf)
b, a = np.polyfit(mwf, mvf, 1); xs = np.linspace(mwf.min(), mwf.max(), 100)
axA.plot(xs, b * xs + a, 'k--', lw=1.4)
axA.text(0.04, 0.94, f'r = {r:.2f}\n' + ('p < 0.001' if p < 1e-3 else f'p = {p:.3f}') +
         f'\n{n_sub} subj, {len(pts)} ROIs', transform=axA.transAxes, va='top', fontsize=10,
         bbox=dict(boxstyle='round', fc='white', ec='0.7'))
axA.set_xlabel('MWF  T2-GRASE  (fraction)'); axA.set_ylabel('MVF  MIMM (Atlas)  (fraction)')
axA.set_title('MIMM MVF vs T2-GRASE MWF, per ROI'); axA.grid(alpha=0.25)

# --- Bland-Altman (note: different metrics/scales -> bias = metric offset) ---
mean = (mvf + mwf) / 2; diff = mvf - mwf
bias = diff.mean(); sd = diff.std(ddof=1)
axB.scatter(mean, diff, c=si, cmap='tab20', s=18, alpha=0.7)
axB.axhline(bias, color='k', lw=1.4)
axB.axhline(bias + 1.96 * sd, color='0.5', ls='--', lw=1)
axB.axhline(bias - 1.96 * sd, color='0.5', ls='--', lw=1)
axB.text(0.04, 0.94, f'bias = {bias:+.3f}\nLoA ± {1.96*sd:.3f}',
         transform=axB.transAxes, va='top', fontsize=10,
         bbox=dict(boxstyle='round', fc='white', ec='0.7'))
axB.set_xlabel('mean of MVF & MWF'); axB.set_ylabel('MVF − MWF')
axB.set_title('Bland–Altman (bias = MVF/MWF metric offset)'); axB.grid(alpha=0.25)

fig.suptitle(f'MIMM vs the independent gold standard (T2-GRASE MWF), n={n_sub}', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = os.path.join(out_dir, 'cohort_mvf_vs_mwf.png')
fig.savefig(out, dpi=150)
print(f'MVF vs MWF:  r = {r:.3f} (p={p:.2g}),  bias (MVF-MWF) = {bias:+.3f}')
print(f'saved: {out}')
