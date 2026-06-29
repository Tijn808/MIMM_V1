#!/usr/bin/env python3
"""
MIMM MVF vs T2-GRASE MWF agreement, done correctly for clustered data.

Agreement is quantified on ONE whole-NAWM mean per subject, so the subject -- not
the ROI -- is the unit of analysis and there is no pseudoreplication: Bland-Altman
bias + limits of agreement, plus the between-subject correlation. The regional
breakdown across the 50 JHU tracts is reported as a per-region BIAS TABLE
(mean MVF, mean MWF, MVF-MWF bias +/- SD across subjects), not a pooled scatter
that would count many correlated ROIs from one subject as independent.

MWF is an independent myelin reference (T2 relaxometry, shares no input with MIMM);
it is NOT a histological gold standard.

NAWM = brain, FA > 0.20, MWF > 0, excluding lesions where a lesion mask exists.

Outputs to <results>/cohort_analysis/:
  cohort_mvf_vs_mwf.png        NAWM scatter + Bland-Altman (per subject)
  cohort_mwf_region_bias.csv   per-tract mean MVF, mean MWF, bias +/- SD

Usage:  python3 cohort_mwf.py <results_dir>
"""
import sys, os, glob, csv
import numpy as np
import nibabel as nib
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory

# JHU ICBM-DTI-81 labels 1..50, in the pipeline's label order.
JHU_NAMES = [
    'Middle cerebellar peduncle', 'Pontine crossing tract', 'Genu of corpus callosum',
    'Body of corpus callosum', 'Splenium of corpus callosum', 'Fornix (column and body)',
    'Corticospinal tract R', 'Corticospinal tract L', 'Medial lemniscus R', 'Medial lemniscus L',
    'Inferior cerebellar peduncle R', 'Inferior cerebellar peduncle L',
    'Superior cerebellar peduncle R', 'Superior cerebellar peduncle L',
    'Cerebral peduncle R', 'Cerebral peduncle L',
    'Anterior limb of internal capsule R', 'Anterior limb of internal capsule L',
    'Posterior limb of internal capsule R', 'Posterior limb of internal capsule L',
    'Retrolenticular internal capsule R', 'Retrolenticular internal capsule L',
    'Anterior corona radiata R', 'Anterior corona radiata L',
    'Superior corona radiata R', 'Superior corona radiata L',
    'Posterior corona radiata R', 'Posterior corona radiata L',
    'Posterior thalamic radiation R', 'Posterior thalamic radiation L',
    'Sagittal stratum R', 'Sagittal stratum L', 'External capsule R', 'External capsule L',
    'Cingulum (cingulate gyrus) R', 'Cingulum (cingulate gyrus) L',
    'Cingulum (hippocampus) R', 'Cingulum (hippocampus) L',
    'Fornix / stria terminalis R', 'Fornix / stria terminalis L',
    'Superior longitudinal fasciculus R', 'Superior longitudinal fasciculus L',
    'Superior fronto-occipital fasciculus R', 'Superior fronto-occipital fasciculus L',
    'Uncinate fasciculus R', 'Uncinate fasciculus L', 'Tapetum R', 'Tapetum L',
    'Fornix (cres) R', 'Fornix (cres) L',
]
JHU_N = len(JHU_NAMES)

if len(sys.argv) < 2:
    sys.exit('usage: cohort_mwf.py <results_dir>')
results_dir = sys.argv[1]
out_dir = os.path.join(results_dir, 'cohort_analysis'); os.makedirs(out_dir, exist_ok=True)


def load(p):
    return np.asarray(nib.load(p).dataobj) if os.path.exists(p) else None


subjects, nawm_mvf, nawm_mwf = [], [], []
roi_mvf = {r: [] for r in range(1, JHU_N + 1)}   # per-ROI list of per-subject means
roi_mwf = {r: [] for r in range(1, JHU_N + 1)}

# the 10 white-matter ROIs of Sisman et al. (mapped to JHU labels; L+R combined)
PAPER_WM = {
    'Genu of corpus callosum': [3], 'Body of corpus callosum': [4],
    'Splenium of corpus callosum': [5], 'Corticospinal tract': [7, 8],
    'Anterior limb of internal capsule': [17, 18], 'Superior corona radiata': [25, 26],
    'Posterior thalamic radiation (optic radiation)': [29, 30], 'External capsule': [33, 34],
    'Cingulum (cingulate gyrus)': [35, 36], 'Superior longitudinal fasciculus': [41, 42],
}
paper_mvf = {k: [] for k in PAPER_WM}
paper_mwf = {k: [] for k in PAPER_WM}

for mp in sorted(glob.glob(os.path.join(results_dir, '*', 'grase', 'MWF.nii.gz'))):
    d = os.path.dirname(os.path.dirname(mp)); sid = os.path.basename(d)
    mwf = load(mp)
    mvf = load(os.path.join(d, 'mimm', 'MVF_Atlas.nii.gz'))
    labels = load(os.path.join(d, 'atlas', 'JHU_labels_subj.nii.gz'))
    brain = load(os.path.join(d, 'qsm', 'brain_mask.nii.gz'))
    fa = load(os.path.join(d, 'atlas', 'FA_atlas.nii.gz'))
    if any(v is None for v in (mwf, mvf, labels, brain, fa)):
        print(f'[skip] {sid} (missing input)'); continue
    mwf = np.clip(mwf.astype(float), 0, 0.5); mvf = mvf.astype(float)
    labels = labels.astype(int)
    nawm = (brain > 0) & (fa > 0.20) & np.isfinite(mvf) & np.isfinite(mwf) & (mwf > 0)
    lesion = load(os.path.join(d, 'lesion', 'lesion_mask.nii.gz'))
    if lesion is not None:
        nawm &= ~(lesion > 0.5)
    if nawm.sum() < 100:
        print(f'[skip] {sid} (NAWM too small)'); continue
    subjects.append(sid)
    nawm_mvf.append(float(mvf[nawm].mean()))
    nawm_mwf.append(float(mwf[nawm].mean()))
    for roi in range(1, JHU_N + 1):
        m = nawm & (labels == roi)
        if m.sum() >= 20:
            roi_mvf[roi].append(float(mvf[m].mean()))
            roi_mwf[roi].append(float(mwf[m].mean()))
    for name, labs in PAPER_WM.items():           # L+R combined paper ROIs
        m = nawm & np.isin(labels, labs)
        if m.sum() >= 20:
            paper_mvf[name].append(float(mvf[m].mean()))
            paper_mwf[name].append(float(mwf[m].mean()))
    print(f'{sid}: NAWM {int(nawm.sum())} vox')

n_sub = len(subjects)
if n_sub < 4:
    sys.exit('Too few subjects with a registered MWF map for a per-subject analysis.')
nawm_mvf = np.array(nawm_mvf); nawm_mwf = np.array(nawm_mwf)

# --- per-subject agreement on whole NAWM ---
diff = nawm_mvf - nawm_mwf
bias = diff.mean(); sd = diff.std(ddof=1)
r, p = stats.pearsonr(nawm_mwf, nawm_mvf)
print(f'\nNAWM agreement, subject as unit (n={n_sub}):')
print(f'  Bland-Altman bias (MVF-MWF) = {bias:+.4f}, LoA +/- {1.96*sd:.4f}')
print(f'  between-subject r = {r:+.3f} (p={p:.3f})')

fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 5.2))
# work in percent to match the reference (Sisman et al. 2025, Fig 4) BA style
mvf_pct, mwf_pct = nawm_mvf * 100, nawm_mwf * 100

# --- scatter ---
axA.scatter(mwf_pct, mvf_pct, s=45, c='#1f77b4', alpha=0.85)
b, a = np.polyfit(mwf_pct, mvf_pct, 1); xs = np.linspace(mwf_pct.min(), mwf_pct.max(), 100)
axA.plot(xs, b * xs + a, 'k--', lw=1.3)
axA.text(0.04, 0.94, f'r = {r:.2f}, p = {p:.3f}\nn = {n_sub} subjects (NAWM)',
         transform=axA.transAxes, va='top', fontsize=10,
         bbox=dict(boxstyle='round', fc='white', ec='0.7'))
axA.set_xlabel('MWF  T2-GRASE  (NAWM mean, %)')
axA.set_ylabel('MVF  MIMM Atlas  (NAWM mean, %)')
axA.set_title('MIMM MVF vs MWF, whole-NAWM per subject'); axA.grid(alpha=0.25)

# --- Bland-Altman, paper style: %, red MEAN / +-1.96SD lines each labelled ---
mean_pct = (mvf_pct + mwf_pct) / 2; diff_pct = mvf_pct - mwf_pct
bias_p = diff_pct.mean(); sd_p = diff_pct.std(ddof=1)
hi, lo = bias_p + 1.96 * sd_p, bias_p - 1.96 * sd_p
axB.scatter(mean_pct, diff_pct, s=45, c='#1f77b4', alpha=0.85)
axB.axhline(bias_p, color='#d62728', lw=1.8)
axB.axhline(hi, color='#d62728', ls='--', lw=1.3)
axB.axhline(lo, color='#d62728', ls='--', lw=1.3)
tform = blended_transform_factory(axB.transAxes, axB.transData)
axB.text(0.99, bias_p, f'MEAN: {bias_p:.2f}', transform=tform, va='bottom', ha='right',
         color='#d62728', fontsize=9)
axB.text(0.99, hi, f'+1.96SD: {hi:.2f}', transform=tform, va='bottom', ha='right',
         color='#d62728', fontsize=9)
axB.text(0.99, lo, f'-1.96SD: {lo:.2f}', transform=tform, va='top', ha='right',
         color='#d62728', fontsize=9)
axB.set_xlabel('Mean of MVF & MWF (%)'); axB.set_ylabel('MVF - MWF (%)')
axB.set_title('Bland-Altman (NAWM, per subject)'); axB.grid(alpha=0.25)

fig.suptitle(f'MIMM MVF vs the independent T2-GRASE MWF reference (NAWM, n={n_sub})', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = os.path.join(out_dir, 'cohort_mvf_vs_mwf.png')
fig.savefig(out, dpi=150)
print(f'saved: {out}')

def paired_stats(a, bb):
    """paired bias (a-b), SD, t-test p, Cohen's d for paired per-subject means."""
    dd = a - bb; bias = dd.mean(); sd = dd.std(ddof=1)
    p = stats.ttest_rel(a, bb).pvalue if len(a) >= 6 else float('nan')
    return bias, sd, p, (bias / sd if sd > 0 else float('nan'))


def write_bias_table(path, items):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['tract', 'n_subj', 'MVF_mean', 'MWF_mean', 'bias_MVF_minus_MWF',
                    'bias_SD', 'paired_p', 'cohens_d'])
        for nm, a, bb in items:
            if len(a) < 4:
                continue
            bias, sd, p, dv = paired_stats(a, bb)
            w.writerow([nm, len(a), f'{a.mean():.4f}', f'{bb.mean():.4f}',
                        f'{bias:+.4f}', f'{sd:.4f}', f'{p:.4g}', f'{dv:+.2f}'])


# --- full 50-label per-region table (appendix) ---
full = [(JHU_NAMES[r - 1], np.array(roi_mvf[r]), np.array(roi_mwf[r])) for r in range(1, JHU_N + 1)]
tbl = os.path.join(out_dir, 'cohort_mwf_region_bias.csv')
write_bias_table(tbl, full)
print(f'saved: {tbl}  (50 labels, with paired p + Cohen\'s d)')

# --- Sisman et al. 10 WM ROIs, L+R combined, with paired test (the body table) ---
paper = [(nm, np.array(paper_mvf[nm]), np.array(paper_mwf[nm])) for nm in PAPER_WM]
ptbl = os.path.join(out_dir, 'cohort_mwf_paper_wm.csv')
write_bias_table(ptbl, paper)
print(f'saved: {ptbl}  (10 WM ROIs, L+R combined)')
print('\nSisman WM ROIs (MVF vs MWF, paired):')
for nm, a, bb in paper:
    if len(a) < 4:
        continue
    bias, sd, p, dv = paired_stats(a, bb)
    star = ' *' if (np.isfinite(p) and p < 0.05) else ''
    print(f'  {nm:46s} MVF {a.mean():.3f}  MWF {bb.mean():.3f}  '
          f'bias {bias:+.4f} (p={p:.2g}, d={dv:+.2f}){star}')
