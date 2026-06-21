#!/usr/bin/env python3
"""
Two follow-up lesion analyses:

1. EROSION ROBUSTNESS — repeats the lesion-vs-NAWM comparison with the lesion
   mask eroded (CSF/NAWM-adjacent rim removed), for MVF and chi_neg. If the drop
   survives erosion, it is driven by the lesion core, not partial volume at the
   edge -> answers the partial-volume objection.

2. LESION BURDEN — per patient, total lesion volume vs whole-brain NAWM MVF.
   Tests whether diffuse myelin loss tracks lesion load (a disease-severity proxy).

Outputs to <results>/cohort_analysis/:
  cohort_lesion_erosion.png    NAWM vs full-lesion vs eroded-core (MVF, |chi_neg|)
  cohort_lesion_burden.png     lesion volume vs NAWM MVF, per patient

Usage:  python3 cohort_lesion_extra.py <results_dir>
"""
import sys, os, glob
import numpy as np
import nibabel as nib
from scipy import ndimage, stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    sys.exit('usage: cohort_lesion_extra.py <results_dir>')
results_dir = sys.argv[1]
out_dir = os.path.join(results_dir, 'cohort_analysis'); os.makedirs(out_dir, exist_ok=True)


def load(p):
    return np.asarray(nib.load(p).dataobj) if os.path.exists(p) else None

MAPS = [('MVF (Atlas)', 'mimm/MVF_Atlas.nii.gz', lambda x: x),
        ('|chi_neg|',    'chisep/chi_neg.nii.gz', np.abs)]

rows = []
for ld in sorted(glob.glob(os.path.join(results_dir, '*', 'lesion', 'lesion_mask.nii.gz'))):
    d = os.path.dirname(os.path.dirname(ld)); sid = os.path.basename(d)
    lesion = load(ld); brain = load(os.path.join(d, 'qsm', 'brain_mask.nii.gz'))
    fa = load(os.path.join(d, 'atlas', 'FA_atlas.nii.gz'))
    if any(v is None for v in (lesion, brain, fa)):
        print(f'[skip] {sid} (missing input)'); continue
    lesion = (lesion > 0.5) & (brain > 0)
    if lesion.sum() < 100:
        print(f'[skip] {sid} ({int(lesion.sum())} vox)'); continue
    eroded = ndimage.binary_erosion(lesion, iterations=1)   # drop 1-voxel rim
    nawm = (brain > 0) & (fa > 0.20) & ~lesion
    zooms = nib.load(ld).header.get_zooms()[:3]
    vox_mL = float(np.prod(zooms)) / 1000.0
    r = {'subject': sid, 'lesion_vol_mL': float(lesion.sum()) * vox_mL,
         'n_eroded': int(eroded.sum())}
    for label, rel, fn in MAPS:
        vol = load(os.path.join(d, rel))
        if vol is None:
            r[label] = (np.nan, np.nan, np.nan); continue
        vol = fn(vol.astype(float))
        def mean(m):
            v = vol[m]; v = v[np.isfinite(v)]; return float(np.mean(v)) if v.size else np.nan
        r[label] = (mean(nawm), mean(lesion), mean(eroded) if eroded.sum() >= 20 else np.nan)
    rows.append(r)
    print(f'{sid}: {r["lesion_vol_mL"]:.1f} mL  MVF NAWM/full/eroded '
          + '/'.join(f'{v:.3f}' for v in r['MVF (Atlas)']))

if not rows:
    sys.exit('No subjects with lesions.')
n = len(rows)

# ---- Figure 1: erosion robustness ----
fig, axes = plt.subplots(1, len(MAPS), figsize=(5 * len(MAPS), 5))
for ax, (label, _, _) in zip(axes, MAPS):
    nawm = np.array([r[label][0] for r in rows])
    full = np.array([r[label][1] for r in rows])
    erod = np.array([r[label][2] for r in rows])
    for a, b, c in zip(nawm, full, erod):
        ax.plot([0, 1, 2], [a, b, c], '-', color='0.6', lw=0.8, alpha=0.6)
    for xi, v, col in [(0, nawm, '#1f77b4'), (1, full, '#d62728'), (2, erod, '#9467bd')]:
        m = np.isfinite(v)
        ax.plot([xi]*m.sum(), v[m], 'o', color=col)
        ax.errorbar(xi, v[m].mean(), yerr=v[m].std(ddof=1)/np.sqrt(m.sum()),
                    fmt='_', color='k', ms=20, capsize=5, lw=2)
    m1 = np.isfinite(nawm) & np.isfinite(full); m2 = np.isfinite(nawm) & np.isfinite(erod)
    p_full = stats.ttest_rel(nawm[m1], full[m1]).pvalue
    p_erod = stats.ttest_rel(nawm[m2], erod[m2]).pvalue
    pc_full = 100*(full[m1].mean()-nawm[m1].mean())/nawm[m1].mean()
    pc_erod = 100*(erod[m2].mean()-nawm[m2].mean())/nawm[m2].mean()
    ax.set_xticks([0,1,2]); ax.set_xticklabels(['NAWM', 'lesion\n(full)', 'lesion\n(eroded)'])
    ax.set_title(f'{label}\nfull {pc_full:+.1f}% (p={p_full:.2g})\n'
                 f'eroded {pc_erod:+.1f}% (p={p_erod:.2g})', fontsize=10)
    ax.grid(alpha=0.25, axis='y')
fig.suptitle(f'Lesion effect survives mask erosion (n={n})', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(out_dir, 'cohort_lesion_erosion.png'), dpi=150)

# ---- Figure 2: lesion burden vs NAWM myelin ----
vol = np.array([r['lesion_vol_mL'] for r in rows])
nawm_mvf = np.array([r['MVF (Atlas)'][0] for r in rows])
m = np.isfinite(vol) & np.isfinite(nawm_mvf)
fig, ax = plt.subplots(figsize=(6, 5))
ax.scatter(vol[m], nawm_mvf[m], s=40, c='#1f77b4', alpha=0.8)
rr, pp = stats.pearsonr(vol[m], nawm_mvf[m])
if m.sum() > 2:
    b, a = np.polyfit(vol[m], nawm_mvf[m], 1)
    xs = np.linspace(vol[m].min(), vol[m].max(), 100); ax.plot(xs, b*xs+a, 'k--', lw=1.3)
ax.text(0.04, 0.06, f'r = {rr:.2f}, ' + ('p < 0.001' if pp < 1e-3 else f'p = {pp:.3f}')
        + f'  (n={m.sum()})', transform=ax.transAxes,
        bbox=dict(boxstyle='round', fc='white', ec='0.7'))
ax.set_xlabel('Lesion volume (mL)'); ax.set_ylabel('NAWM MVF (Atlas)')
ax.set_title('Diffuse myelin vs lesion burden'); ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'cohort_lesion_burden.png'), dpi=150)

# ---- summary ----
print(f'\nErosion robustness (n={n}):')
for label, _, _ in MAPS:
    nawm = np.array([r[label][0] for r in rows]); full = np.array([r[label][1] for r in rows])
    erod = np.array([r[label][2] for r in rows])
    m1 = np.isfinite(nawm)&np.isfinite(full); m2 = np.isfinite(nawm)&np.isfinite(erod)
    print(f'  {label:12s} full {100*(full[m1].mean()-nawm[m1].mean())/nawm[m1].mean():+.1f}% '
          f'(p={stats.ttest_rel(nawm[m1],full[m1]).pvalue:.2g}),  '
          f'eroded {100*(erod[m2].mean()-nawm[m2].mean())/nawm[m2].mean():+.1f}% '
          f'(p={stats.ttest_rel(nawm[m2],erod[m2]).pvalue:.2g})')
print(f'Lesion burden vs NAWM MVF: r = {rr:.3f}, p = {pp:.3g}')
print(f'saved figures to {out_dir}')
