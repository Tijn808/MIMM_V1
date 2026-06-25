#!/usr/bin/env python3
"""
fig_resolution.py  —  Fig 7: lesion detection vs size, MIMM MVF vs MWF.

Standalone, publication-styled version. MWF is acquired at 1.5x1.5x4 mm, so a small
lesion is a fraction of one MWF voxel and is washed out by partial volume; MIMM at
1 mm isotropic should still resolve it. For every FLAIR-confirmed lesion the % drop
of MVF and of MWF (each vs its own location-matched perilesional NAWM shell) is
measured, then plotted against lesion volume (log axis) with tertile-binned means.

FALSE-POSITIVE GUARD: a connected component is kept only if it is genuinely
hyperintense on the registered FLAIR relative to its perilesional rim (>= +5%);
rejected components are counted and reported.

Usage:  python3 fig_resolution.py <results_dir>
Output: <results_dir>/cohort_analysis/thesis_figures/fig4_resolution.png
"""
import sys, os, glob
import numpy as np
import nibabel as nib
from scipy import ndimage, stats
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

mpl.rcParams.update({
    'font.size': 12, 'axes.titlesize': 14, 'axes.labelsize': 13,
    'xtick.labelsize': 11, 'ytick.labelsize': 11, 'legend.fontsize': 11,
    'savefig.dpi': 200, 'axes.grid': True, 'grid.alpha': 0.25, 'grid.linestyle': '--',
})
ATLAS, LESION = '#1f77b4', '#d62728'
GAP, WIDTH = 2, 3
MIN_LES, MIN_PERI = 3, 30          # min voxels per component / per perilesional shell
FLAIR_MARGIN = 0.05                # component must be >=5% brighter than its rim

if len(sys.argv) < 2:
    sys.exit('usage: fig_resolution.py <results_dir>')
RESULTS = sys.argv[1]
OUT = os.path.join(RESULTS, 'cohort_analysis', 'thesis_figures')
os.makedirs(OUT, exist_ok=True)


def load(p):
    return np.asarray(nib.load(p).dataobj) if os.path.exists(p) else None


def m(vol, msk):
    v = vol[msk]; v = v[np.isfinite(v)]
    return float(np.mean(v)) if v.size else np.nan


lesions = []
n_total = n_conf = n_rej_flair = n_rej_peri = 0
for ld in sorted(glob.glob(os.path.join(RESULTS, '*', 'lesion', 'lesion_mask.nii.gz'))):
    d = os.path.dirname(os.path.dirname(ld)); sid = os.path.basename(d)
    mask = load(ld); brain = load(os.path.join(d, 'qsm', 'brain_mask.nii.gz'))
    fa = load(os.path.join(d, 'atlas', 'FA_atlas.nii.gz'))
    mvf = load(os.path.join(d, 'mimm', 'MVF_Atlas.nii.gz'))
    mwf = load(os.path.join(d, 'grase', 'MWF.nii.gz'))
    flair = load(os.path.join(d, 'lesion', 'FLAIR_mgre.nii.gz'))
    if any(v is None for v in (mask, brain, fa, mvf, mwf, flair)):
        print(f'[skip] {sid} (missing input)'); continue
    mask = (mask > 0.5) & (brain > 0)
    mvf = mvf.astype(float); flair = flair.astype(float)
    mwf = np.clip(mwf.astype(float), 0, 0.5)
    nawm = (brain > 0) & (fa > 0.20) & ~mask
    vox_mL = float(np.prod(nib.load(ld).header.get_zooms()[:3])) / 1000.0

    labels, _ = ndimage.label(mask)
    pad = GAP + WIDTH + 1
    for i, sl in enumerate(ndimage.find_objects(labels), start=1):
        if sl is None:
            continue
        esl = tuple(slice(max(0, s.start - pad), min(dim, s.stop + pad))
                    for s, dim in zip(sl, mask.shape))
        comp = labels[esl] == i
        if comp.sum() < MIN_LES:
            continue
        n_total += 1
        inner = ndimage.binary_dilation(comp, iterations=GAP)
        outer = ndimage.binary_dilation(comp, iterations=GAP + WIDTH)
        peri = outer & ~inner & nawm[esl]
        if peri.sum() < MIN_PERI:
            n_rej_peri += 1; continue
        mvf_c, mwf_c, fl_c = mvf[esl], mwf[esl], flair[esl]
        # FLAIR confirmation: a real WMH is brighter than its surroundings
        fl_l, fl_p = m(fl_c, comp), m(fl_c, peri)
        if not (np.isfinite(fl_l) and np.isfinite(fl_p) and fl_p > 0
                and fl_l >= fl_p * (1 + FLAIR_MARGIN)):
            n_rej_flair += 1; continue
        mvf_p, mwf_p = m(mvf_c, peri), m(mwf_c, peri)
        if not (mvf_p > 0 and mwf_p > 0):
            continue
        n_conf += 1
        lesions.append((comp.sum() * vox_mL,
                        100 * (mvf_p - m(mvf_c, comp)) / mvf_p,
                        100 * (mwf_p - m(mwf_c, comp)) / mwf_p))

print(f'components >= {MIN_LES} vox:            {n_total}')
print(f'  rejected (no perilesional shell):  {n_rej_peri}')
print(f'  rejected (not FLAIR-hyperintense): {n_rej_flair}')
print(f'  FLAIR-confirmed lesions:           {n_conf}')
if n_conf < 6:
    sys.exit('Too few confirmed lesions for a size-stratified figure.')

vol = np.array([l[0] for l in lesions])
dmvf = np.array([l[1] for l in lesions])
dmwf = np.array([l[2] for l in lesions])
q1, q2 = np.quantile(vol, [1/3, 2/3])
bins = [(vol <= q1), (vol > q1) & (vol <= q2), (vol > q2)]
xb = [vol[s].mean() for s in bins]

fig, ax = plt.subplots(figsize=(8.5, 6))
ax.scatter(vol, dmvf, s=16, c=ATLAS, alpha=0.4)
ax.scatter(vol, dmwf, s=16, c=LESION, alpha=0.4)
ax.plot(xb, [dmvf[s].mean() for s in bins], '-o', color=ATLAS, lw=2.4, ms=9, label='MIMM MVF (binned)')
ax.plot(xb, [dmwf[s].mean() for s in bins], '-o', color=LESION, lw=2.4, ms=9, label='MWF (binned)')
ax.axhline(0, color='k', lw=0.8)
ax.set_xscale('log')
ax.set_xlabel('lesion volume (mL, log scale)')
ax.set_ylabel('drop vs perilesional NAWM (%)')
ax.set_title(f'Lesion detection vs size: MIMM resolves, MWF blurs\n'
             f'{n_conf} FLAIR-confirmed lesions ({n_rej_flair} rejected as not hyperintense)')
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'fig4_resolution.png'))
print(f'saved: {os.path.join(OUT, "fig4_resolution.png")}  ({n_conf} lesions)')
