#!/usr/bin/env python3
"""
Cohort lesion analysis: MIMM myelin (and corroborating maps) in MS lesions vs
normal-appearing white matter (NAWM), pooled across patients.

For every subject with a registered lesion mask (lesion/lesion_mask.nii.gz,
produced by run_lesion.sh), recomputes the mean of each map in lesion voxels and
in NAWM (brain & FA>0.20 & not-lesion), then makes a paired cohort figure
(NAWM vs lesion per subject) with a paired t-test, for MVF, |chi_neg| and g-ratio.

Usage:
  python3 cohort_lesion.py <results_dir>
"""
import sys, os, glob
import numpy as np
import nibabel as nib
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    sys.exit('usage: cohort_lesion.py <results_dir>')
results_dir = sys.argv[1]
out_dir = os.path.join(results_dir, 'cohort_analysis')
os.makedirs(out_dir, exist_ok=True)


def load(p):
    return np.asarray(nib.load(p).dataobj) if os.path.exists(p) else None

# maps to compare (label -> (relative path, transform))
MAPS = [
    ('MVF (Atlas)',  'mimm/MVF_Atlas.nii.gz', lambda x: x),
    ('|chi_neg|',    'chisep/chi_neg.nii.gz', np.abs),
    ('g-ratio',      'mimm/g_ratio_Atlas.nii.gz', lambda x: x),
]

rows = []   # one dict per subject
for ld in sorted(glob.glob(os.path.join(results_dir, '*', 'lesion', 'lesion_mask.nii.gz'))):
    d = os.path.dirname(os.path.dirname(ld))
    sid = os.path.basename(d)
    lesion = load(ld)
    brain  = load(os.path.join(d, 'qsm', 'brain_mask.nii.gz'))
    fa     = load(os.path.join(d, 'atlas', 'FA_atlas.nii.gz'))
    if lesion is None or brain is None or fa is None:
        print(f'[skip] {sid} (missing lesion/brain/FA)'); continue
    lesion = (lesion > 0.5) & (brain > 0)
    nawm   = (brain > 0) & (fa > 0.20) & ~lesion
    if lesion.sum() < 50:
        print(f'[skip] {sid} ({int(lesion.sum())} lesion voxels)'); continue
    r = {'subject': sid, 'n_lesion': int(lesion.sum())}
    for label, rel, fn in MAPS:
        vol = load(os.path.join(d, rel))
        if vol is None:
            r[label] = (np.nan, np.nan); continue
        vol = fn(vol.astype(float))
        lv = vol[lesion]; nv = vol[nawm]
        lv = lv[np.isfinite(lv)]; nv = nv[np.isfinite(nv)]
        r[label] = (float(np.mean(nv)), float(np.mean(lv)))   # (NAWM, lesion)
    rows.append(r)
    print(f'{sid}: {r["n_lesion"]} lesion vox  '
          + '  '.join(f'{m}:{r[m][0]:.3f}->{r[m][1]:.3f}' for m, _, _ in MAPS if m in r))

if not rows:
    sys.exit('No subjects with a registered lesion mask found.')
n = len(rows)
print(f'\n{n} subjects with lesions.')

# ---- paired cohort figure ----
fig, axes = plt.subplots(1, len(MAPS), figsize=(4.2 * len(MAPS), 5))
for ax, (label, _, _) in zip(axes, MAPS):
    nawm = np.array([r[label][0] for r in rows])
    les  = np.array([r[label][1] for r in rows])
    m = np.isfinite(nawm) & np.isfinite(les)
    nawm, les = nawm[m], les[m]
    for a, b in zip(nawm, les):
        ax.plot([0, 1], [a, b], '-', color='0.6', lw=0.8, alpha=0.7)
    ax.plot([0]*len(nawm), nawm, 'o', color='#1f77b4', label='NAWM')
    ax.plot([1]*len(les),  les,  'o', color='#d62728', label='lesion')
    # mean +/- SE
    for xi, v in [(0, nawm), (1, les)]:
        ax.errorbar(xi, v.mean(), yerr=v.std(ddof=1)/np.sqrt(len(v)),
                    fmt='_', color='k', ms=22, capsize=5, lw=2)
    t, p = stats.ttest_rel(nawm, les)
    pct = 100 * (les.mean() - nawm.mean()) / nawm.mean()
    ax.set_xticks([0, 1]); ax.set_xticklabels(['NAWM', 'lesion'])
    ax.set_xlim(-0.4, 1.4)
    ax.set_title(f'{label}\n{pct:+.1f}%, ' +
                 ('p < 0.001' if p < 1e-3 else f'paired p = {p:.3f}') + f' (n={len(nawm)})')
    ax.grid(alpha=0.25, axis='y')
fig.suptitle(f'MS lesions vs normal-appearing WM, cohort (n={n} patients)', fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = os.path.join(out_dir, 'cohort_lesion_vs_nawm.png')
fig.savefig(out, dpi=150)

# ---- table ----
import csv
with open(os.path.join(out_dir, 'cohort_lesion_vs_nawm.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['subject', 'n_lesion'] + [f'{m}_{s}' for m, _, _ in MAPS for s in ('NAWM', 'lesion')])
    for r in rows:
        w.writerow([r['subject'], r['n_lesion']] +
                   [r[m][k] for m, _, _ in MAPS for k in (0, 1)])

print(f'\nPaired summary across {n} patients:')
for label, _, _ in MAPS:
    nawm = np.array([r[label][0] for r in rows]); les = np.array([r[label][1] for r in rows])
    m = np.isfinite(nawm) & np.isfinite(les); nawm, les = nawm[m], les[m]
    t, p = stats.ttest_rel(nawm, les)
    print(f'  {label:12s} NAWM {nawm.mean():.3f} -> lesion {les.mean():.3f}  '
          f'({100*(les.mean()-nawm.mean())/nawm.mean():+.1f}%, paired p={p:.2g})')
print(f'saved: {out}')
