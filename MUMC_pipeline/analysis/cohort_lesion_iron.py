#!/usr/bin/env python3
"""
Exploratory iron-rim analysis in MS lesions (paramagnetic-rim-lesion idea).

For every subject with a registered lesion mask, splits each subject's lesion
voxels into a RIM (outer 2-voxel shell) and a CORE (interior), and compares the
chi-separation iron map (chi_pos) in rim vs core vs NAWM. Chronic active
(paramagnetic-rim) lesions carry iron at the edge, so a rim > core difference is
the signal of interest.

NOTE: this is a crude cohort average over *all* lesion voxels, not per-lesion PRL
detection. Most lesions are not rim lesions, so a null/weak result is plausible
and is itself informative. Treat as exploratory.

Usage:  python3 cohort_lesion_iron.py <results_dir>
"""
import sys, os, glob
import numpy as np
import nibabel as nib
from scipy import ndimage, stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    sys.exit('usage: cohort_lesion_iron.py <results_dir>')
results_dir = sys.argv[1]
out_dir = os.path.join(results_dir, 'cohort_analysis'); os.makedirs(out_dir, exist_ok=True)


def load(p):
    return np.asarray(nib.load(p).dataobj) if os.path.exists(p) else None

rows = []
for ld in sorted(glob.glob(os.path.join(results_dir, '*', 'lesion', 'lesion_mask.nii.gz'))):
    d = os.path.dirname(os.path.dirname(ld)); sid = os.path.basename(d)
    lesion = load(ld); brain = load(os.path.join(d, 'qsm', 'brain_mask.nii.gz'))
    fa = load(os.path.join(d, 'atlas', 'FA_atlas.nii.gz'))
    chi_pos = load(os.path.join(d, 'chisep', 'chi_pos.nii.gz'))
    if any(v is None for v in (lesion, brain, fa, chi_pos)):
        print(f'[skip] {sid} (missing input)'); continue
    lesion = (lesion > 0.5) & (brain > 0)
    if lesion.sum() < 100:
        print(f'[skip] {sid} ({int(lesion.sum())} lesion vox)'); continue
    core = ndimage.binary_erosion(lesion, iterations=2)
    rim  = lesion & ~core
    nawm = (brain > 0) & (fa > 0.20) & ~lesion
    chi_pos = chi_pos.astype(float)

    def mean_in(mask):
        v = chi_pos[mask]; v = v[np.isfinite(v)]
        return float(np.mean(v)) if v.size else np.nan

    if core.sum() < 20:
        print(f'[skip] {sid} (core too small)'); continue
    r = {'subject': sid, 'rim': mean_in(rim), 'core': mean_in(core), 'nawm': mean_in(nawm),
         'n_rim': int(rim.sum()), 'n_core': int(core.sum())}
    rows.append(r)
    print(f'{sid}: chi+ rim {r["rim"]:.4f}  core {r["core"]:.4f}  NAWM {r["nawm"]:.4f}')

if not rows:
    sys.exit('No subjects with usable lesions.')
n = len(rows)

rim  = np.array([r['rim'] for r in rows]); core = np.array([r['core'] for r in rows])
nawm = np.array([r['nawm'] for r in rows])

fig, ax = plt.subplots(figsize=(6, 5))
for a, b, c in zip(nawm, core, rim):
    ax.plot([0, 1, 2], [a, b, c], '-', color='0.6', lw=0.8, alpha=0.7)
for xi, v, col, lab in [(0, nawm, '#1f77b4', 'NAWM'), (1, core, '#7f7f7f', 'lesion core'),
                        (2, rim, '#d62728', 'lesion rim')]:
    ax.plot([xi] * len(v), v, 'o', color=col)
    ax.errorbar(xi, np.nanmean(v), yerr=np.nanstd(v, ddof=1) / np.sqrt(np.isfinite(v).sum()),
                fmt='_', color='k', ms=22, capsize=5, lw=2)
t, p = stats.ttest_rel(rim, core, nan_policy='omit')
ax.set_xticks([0, 1, 2]); ax.set_xticklabels(['NAWM', 'core', 'rim'])
ax.set_ylabel('chi-sep iron (chi+)  (ppm)')
ax.set_title(f'Lesion iron: rim vs core (n={n})\nrim−core paired p = '
             + ('< 0.001' if p < 1e-3 else f'{p:.3f}'))
ax.grid(alpha=0.25, axis='y')
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'cohort_lesion_iron_rim.png'), dpi=150)
print(f'\nrim {np.nanmean(rim):.4f} vs core {np.nanmean(core):.4f}  '
      f'(paired p={p:.3g}, n={n})')
print(f'saved: {os.path.join(out_dir, "cohort_lesion_iron_rim.png")}')
