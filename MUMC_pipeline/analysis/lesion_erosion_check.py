#!/usr/bin/env python3
"""
lesion_erosion_check.py -- is the slide-10 compartment drop driven by CSF in
periventricular lesion masks?

Concern: a periventricular lesion mask can include a rim of CSF partial-volume
voxels at the ventricle edge. Those read ~0 MVF/AVF and would pull the lesion
mean down, inflating the apparent myelin/axon loss.

Test: re-run the lesion vs perilesional-NAWM comparison (MVF, AVF, FVF) at
several mask EROSION levels. Eroding removes the outer rim (the voxels most
likely to be CSF). The perilesional shell is held FIXED (defined from the
original mask) so only the lesion-interior measurement changes. If the % drops
survive erosion, CSF is not driving the result; if they collapse, CSF partial
volume was inflating them.

The cohort is held fixed across erosion levels (only subjects whose mask still
has >= MIN_LES voxels after the LARGEST erosion are included), so the columns
are directly comparable, not measured on different subsets.

Also reports the fraction of lesion voxels that are CSF-like (very low FA AND
very low MVF) as a direct, if approximate, measure of contamination.

Usage:  python3 lesion_erosion_check.py <results_dir>
"""
import sys, os, glob
import numpy as np
import nibabel as nib
from scipy import ndimage, stats

GAP, WIDTH = 2, 3            # perilesional shell: 2-voxel gap, 3 voxels wide
MIN_LES, MIN_PERI = 50, 50
EROSIONS = [0, 1, 2]        # voxels eroded from the lesion mask

if len(sys.argv) < 2:
    sys.exit('usage: lesion_erosion_check.py <results_dir>')
RESULTS = sys.argv[1]


def load(p):
    return np.asarray(nib.load(p).dataobj, dtype=float) if os.path.exists(p) else None


def mean(vol, m):
    v = vol[m]; v = v[np.isfinite(v)]
    return float(np.mean(v)) if v.size else np.nan


rows = {e: [] for e in EROSIONS}   # per-subject (peri_mean, lesion_mean) per compartment
csf_fracs = []
emax = max(EROSIONS)

for ld in sorted(glob.glob(os.path.join(RESULTS, '*', 'lesion', 'lesion_mask.nii.gz'))):
    d = os.path.dirname(os.path.dirname(ld)); sid = os.path.basename(d)
    lesion0 = load(ld); brain = load(os.path.join(d, 'qsm', 'brain_mask.nii.gz'))
    fa = load(os.path.join(d, 'atlas', 'FA_atlas.nii.gz'))
    mvf = load(os.path.join(d, 'mimm', 'MVF_Atlas.nii.gz'))
    fvf = load(os.path.join(d, 'mimm', 'FVF_Atlas.nii.gz'))
    if any(v is None for v in (lesion0, brain, fa, mvf, fvf)):
        continue
    lesion0 = (lesion0 > 0.5) & (brain > 0)
    if lesion0.sum() < MIN_LES:
        continue
    # fix the cohort: mask must survive the largest erosion
    les_max = ndimage.binary_erosion(lesion0, iterations=emax) if emax else lesion0
    if les_max.sum() < MIN_LES:
        continue
    mvf = mvf.astype(float); fvf = fvf.astype(float); avf = fvf - mvf

    # perilesional NAWM shell, held FIXED from the original mask
    nawm = (brain > 0) & (fa > 0.20) & ~lesion0
    inner = ndimage.binary_dilation(lesion0, iterations=GAP)
    outer = ndimage.binary_dilation(lesion0, iterations=GAP + WIDTH)
    peri = outer & ~inner & nawm
    if peri.sum() < MIN_PERI:
        continue

    # direct CSF-contamination measure in the original mask
    csf_like = lesion0 & (fa < 0.10) & (mvf < 0.05)
    csf_fracs.append(csf_like.sum() / lesion0.sum())

    pm, pa, pf = mean(mvf, peri), mean(avf, peri), mean(fvf, peri)
    for e in EROSIONS:
        les_e = ndimage.binary_erosion(lesion0, iterations=e) if e else lesion0
        rows[e].append({'MVF': (pm, mean(mvf, les_e)),
                        'AVF': (pa, mean(avf, les_e)),
                        'FVF': (pf, mean(fvf, les_e))})

if not csf_fracs:
    sys.exit('no usable lesions found')

print(f'CSF-like fraction of lesion voxels (FA<0.10 AND MVF<0.05):  '
      f'mean {100*np.mean(csf_fracs):.1f}%   median {100*np.median(csf_fracs):.1f}%   '
      f'max {100*np.max(csf_fracs):.1f}%    (n={len(csf_fracs)} subjects)\n')

print(f'{"compartment":11s} {"mask":>7s} {"n":>3s} {"peri":>7s} {"lesion":>7s} '
      f'{"%change":>8s} {"p":>9s} {"d":>7s}')
for comp in ['MVF', 'AVF', 'FVF']:
    for e in EROSIONS:
        r = rows[e]
        peri = np.array([x[comp][0] for x in r]); les = np.array([x[comp][1] for x in r])
        m = np.isfinite(peri) & np.isfinite(les); peri, les = peri[m], les[m]
        if len(peri) < 3:
            continue
        pct = 100 * (les.mean() - peri.mean()) / peri.mean()
        dd = les - peri
        p = stats.ttest_rel(les, peri).pvalue
        dcohen = dd.mean() / dd.std(ddof=1)
        tag = 'orig' if e == 0 else f'-{e}vox'
        print(f'{comp:11s} {tag:>7s} {len(peri):>3d} {peri.mean():>7.3f} {les.mean():>7.3f} '
              f'{pct:>+7.1f}% {p:>9.2g} {dcohen:>+7.2f}')
    print()

print('Read: if %change and d stay ~stable from orig -> -2vox, the compartment')
print('loss is real tissue loss, not CSF partial volume in the mask rim.')
