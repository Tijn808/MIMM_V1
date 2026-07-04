#!/usr/bin/env python3
"""
fig_compartments.py  —  Fig 6: lesion compartment decomposition (myelin + axon).

Standalone version of make_thesis_figures.fig3(). For every patient with a lesion
mask, each lesion is compared against a LOCATION-MATCHED perilesional NAWM shell
(2-5 voxels outside the lesion). MIMM resolves the lesion into:
    MVF  = myelin volume fraction
    FVF  = fibre volume fraction (axon + myelin)
    AVF  = FVF - MVF  = axon volume fraction
Three paired panels show peri-NAWM vs lesion for each compartment, annotated with
the % change, paired-t p-value and Cohen's d. The point: lesions lose BOTH myelin
and axon (they are atrophic), which neither MWF nor chi-separation can show.

Usage:  python3 fig_compartments.py <results_dir>
Output: <results_dir>/cohort_analysis/thesis_figures/fig3_compartments.png
"""
import sys, os, glob
import numpy as np
import nibabel as nib
from scipy import ndimage, stats
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mimm_style import apply_style, C, COMPARTMENT   # shared deck palette (style only, numbers unchanged)
apply_style()
GAP, WIDTH = 2, 3          # perilesional shell: 2 voxels gap, 3 voxels wide
MIN_LES, MIN_PERI = 50, 50

if len(sys.argv) < 2:
    sys.exit('usage: fig_compartments.py <results_dir>')
RESULTS = sys.argv[1]
OUT = os.path.join(RESULTS, 'cohort_analysis', 'thesis_figures')
os.makedirs(OUT, exist_ok=True)


def load(p):
    return np.asarray(nib.load(p).dataobj) if os.path.exists(p) else None


def mean(vol, m):
    v = vol[m]; v = v[np.isfinite(v)]
    return float(np.mean(v)) if v.size else np.nan


rows = []
for ld in sorted(glob.glob(os.path.join(RESULTS, '*', 'lesion', 'lesion_mask.nii.gz'))):
    d = os.path.dirname(os.path.dirname(ld)); sid = os.path.basename(d)
    lesion = load(ld); brain = load(os.path.join(d, 'qsm', 'brain_mask.nii.gz'))
    fa = load(os.path.join(d, 'atlas', 'FA_atlas.nii.gz'))
    mvf = load(os.path.join(d, 'mimm', 'MVF_Atlas.nii.gz'))
    fvf = load(os.path.join(d, 'mimm', 'FVF_Atlas.nii.gz'))
    if any(v is None for v in (lesion, brain, fa, mvf, fvf)):
        print(f'[skip] {sid} (missing input)'); continue
    lesion = (lesion > 0.5) & (brain > 0)
    if lesion.sum() < MIN_LES:
        continue
    nawm = (brain > 0) & (fa > 0.20) & ~lesion
    inner = ndimage.binary_dilation(lesion, iterations=GAP)
    outer = ndimage.binary_dilation(lesion, iterations=GAP + WIDTH)
    peri = outer & ~inner & nawm
    if peri.sum() < MIN_PERI:
        continue
    mvf = mvf.astype(float); fvf = fvf.astype(float); avf = fvf - mvf
    # whole-brain NAWM, location-matched perilesional shell, lesion
    rows.append({'MVF (myelin)': (mean(mvf, nawm), mean(mvf, peri), mean(mvf, lesion)),
                 'AVF (axon)':   (mean(avf, nawm), mean(avf, peri), mean(avf, lesion)),
                 'FVF (fibre)':  (mean(fvf, nawm), mean(fvf, peri), mean(fvf, lesion))})

n = len(rows)
if n < 3:
    sys.exit(f'Only {n} usable patients - need >=3.')

fig, axes = plt.subplots(1, 3, figsize=(13, 5))
for ax, label in zip(axes, ['MVF (myelin)', 'AVF (axon)', 'FVF (fibre)']):
    whole = np.array([r[label][0] for r in rows])
    peri = np.array([r[label][1] for r in rows]); les = np.array([r[label][2] for r in rows])
    m = np.isfinite(whole) & np.isfinite(peri) & np.isfinite(les)
    whole, peri, les = whole[m], peri[m], les[m]
    for a, b, c in zip(whole, peri, les):
        ax.plot([0, 1, 2], [a, b, c], '-', color='0.65', lw=0.9, alpha=0.7)
    key = label.split()[0]                       # MVF / AVF / FVF -> compartment teal shade
    ax.plot([0]*len(whole), whole, 'o', color=C['secondary'], ms=6)   # NAWM (whole) grey
    ax.plot([1]*len(peri), peri, 'o', color=C['reference'], ms=6)     # NAWM (peri) grey
    ax.plot([2]*len(les), les, 'o', color=COMPARTMENT[key], ms=7)     # lesion = compartment colour
    for xi, v in [(0, whole), (1, peri), (2, les)]:
        ax.errorbar(xi, v.mean(), yerr=v.std(ddof=1)/np.sqrt(len(v)),
                    fmt='_', color='k', ms=26, capsize=6, lw=2.2)
    p = stats.ttest_rel(les, peri).pvalue
    pct = 100 * (les.mean() - peri.mean()) / peri.mean()
    d = (les - peri).mean() / (les - peri).std(ddof=1)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['NAWM\n(whole)', 'NAWM\n(peri)', 'lesion']); ax.set_xlim(-0.4, 2.4)
    ax.set_title(f'{label}   lesion vs peri-NAWM\n{pct:+.1f}%, '
                 + ('p < 0.001' if p < 1e-3 else f'p = {p:.3f}') + f', d = {d:+.2f}')
fig.suptitle(f'Lesion compartment decomposition vs perilesional NAWM (n = {n})', fontsize=15)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(OUT, 'fig3_compartments.png'))
print(f'saved: {os.path.join(OUT, "fig3_compartments.png")}  (n={n})')
