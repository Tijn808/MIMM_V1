#!/usr/bin/env python3
"""
Does MIMM out-detect MWF on SMALL lesions? MWF is 1.5x1.5x4 mm, so a small lesion
is a fraction of one MWF voxel and its signal is washed out by partial volume;
MIMM at 1 mm iso should still resolve it. This is the resolution advantage tested
directly, per individual lesion.

GUARD AGAINST FALSE POSITIVES: small connected components of a SAMSEG mask can be
spurious. Each component is kept only if it is genuinely hyperintense on the
registered FLAIR relative to its own perilesional rim (a real WMH is bright on
FLAIR). The number rejected is reported — that is the "other small things" count.

For every FLAIR-confirmed lesion: measure the % drop of MVF and of MWF relative to
a perilesional NAWM shell (location-matched). Lesions are then binned by volume;
within each bin both drops are tested against zero. If, in the small-lesion bin,
the MVF drop is significant while the MWF drop is not, MIMM detects small lesions
that MWF misses.

Usage:  python3 cohort_lesion_size_mwf.py <results_dir>
"""
import sys, os, glob, csv
import numpy as np
import nibabel as nib
from scipy import ndimage, stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

GAP, WIDTH = 2, 3
MIN_LES, MIN_PERI = 3, 30      # min voxels per component / per perilesional shell
FLAIR_MARGIN = 0.05            # component must be >=5% brighter than its rim to count

if len(sys.argv) < 2:
    sys.exit('usage: cohort_lesion_size_mwf.py <results_dir>')
results_dir = sys.argv[1]
out_dir = os.path.join(results_dir, 'cohort_analysis'); os.makedirs(out_dir, exist_ok=True)


def load(p):
    return np.asarray(nib.load(p).dataobj) if os.path.exists(p) else None

lesions = []          # one dict per FLAIR-confirmed lesion
n_total = n_conf = n_rej_flair = n_rej_peri = 0
for ld in sorted(glob.glob(os.path.join(results_dir, '*', 'lesion', 'lesion_mask.nii.gz'))):
    d = os.path.dirname(os.path.dirname(ld)); sid = os.path.basename(d)
    mask  = load(ld); brain = load(os.path.join(d, 'qsm', 'brain_mask.nii.gz'))
    fa    = load(os.path.join(d, 'atlas', 'FA_atlas.nii.gz'))
    mvf   = load(os.path.join(d, 'mimm', 'MVF_Atlas.nii.gz'))
    mwf   = load(os.path.join(d, 'grase', 'MWF.nii.gz'))
    flair = load(os.path.join(d, 'lesion', 'FLAIR_mgre.nii.gz'))
    if any(v is None for v in (mask, brain, fa, mvf, mwf, flair)):
        miss = [n for n, v in (('lesion', mask), ('brain', brain), ('FA', fa),
                ('MVF', mvf), ('MWF', mwf), ('FLAIR', flair)) if v is None]
        print(f'[skip] {sid} (missing: {", ".join(miss)})'); continue
    mask  = (mask > 0.5) & (brain > 0)
    mvf   = mvf.astype(float); flair = flair.astype(float)
    mwf   = np.clip(mwf.astype(float), 0, 0.5)
    nawm  = (brain > 0) & (fa > 0.20) & ~mask     # excludes ALL lesions, not just this one
    zooms = nib.load(ld).header.get_zooms()[:3]
    vox_mL = float(np.prod(zooms)) / 1000.0

    labels, nlab = ndimage.label(mask)
    for i in range(1, nlab + 1):
        comp = labels == i
        if comp.sum() < MIN_LES:
            continue
        n_total += 1
        inner = ndimage.binary_dilation(comp, iterations=GAP)
        outer = ndimage.binary_dilation(comp, iterations=GAP + WIDTH)
        peri  = outer & ~inner & nawm
        if peri.sum() < MIN_PERI:
            n_rej_peri += 1; continue

        def m(vol, msk):
            v = vol[msk]; v = v[np.isfinite(v)]; return float(np.mean(v)) if v.size else np.nan

        # FLAIR confirmation: real WMH is brighter than its surroundings
        fl_les, fl_peri = m(flair, comp), m(flair, peri)
        if not (np.isfinite(fl_les) and np.isfinite(fl_peri) and fl_peri > 0
                and fl_les >= fl_peri * (1 + FLAIR_MARGIN)):
            n_rej_flair += 1; continue
        n_conf += 1

        mvf_l, mvf_p = m(mvf, comp), m(mvf, peri)
        mwf_l, mwf_p = m(mwf, comp), m(mwf, peri)
        if not all(np.isfinite(x) and x > 0 for x in (mvf_p, mwf_p)):
            continue
        lesions.append({'subject': sid, 'vol_mL': comp.sum() * vox_mL,
                        'mvf_drop_pct': 100 * (mvf_p - mvf_l) / mvf_p,
                        'mwf_drop_pct': 100 * (mwf_p - mwf_l) / mwf_p})

print(f'\nComponents >= {MIN_LES} vox: {n_total}')
print(f'  rejected (no perilesional shell):  {n_rej_peri}')
print(f'  rejected (not FLAIR-hyperintense): {n_rej_flair}   <- the "other small things"')
print(f'  FLAIR-confirmed lesions analysed:  {n_conf}')
if len(lesions) < 6:
    sys.exit('Too few confirmed lesions for a size-stratified test.')

vol  = np.array([l['vol_mL'] for l in lesions])
dmvf = np.array([l['mvf_drop_pct'] for l in lesions])
dmwf = np.array([l['mwf_drop_pct'] for l in lesions])

# size bins by tertiles of volume
q1, q2 = np.quantile(vol, [1/3, 2/3])
bins = [('small', vol <= q1), ('medium', (vol > q1) & (vol <= q2)), ('large', vol > q2)]
print(f'\nVolume tertiles: small <= {q1:.2f} mL < medium <= {q2:.2f} mL < large')
print(f'{"bin":7s} {"n":>4} {"MVF drop":>10} {"p(MVF)":>9}   {"MWF drop":>10} {"p(MWF)":>9}')
summary = []
for name, sel in bins:
    n = int(sel.sum())
    mv, mw = dmvf[sel], dmwf[sel]
    p_mv = stats.ttest_1samp(mv, 0, nan_policy='omit').pvalue
    p_mw = stats.ttest_1samp(mw, 0, nan_policy='omit').pvalue
    print(f'{name:7s} {n:>4} {np.nanmean(mv):>9.1f}% {p_mv:>9.3g}   {np.nanmean(mw):>9.1f}% {p_mw:>9.3g}')
    summary.append({'bin': name, 'n': n, 'mvf_drop_pct': round(float(np.nanmean(mv)), 1),
                    'mvf_p': float(p_mv), 'mwf_drop_pct': round(float(np.nanmean(mw)), 1),
                    'mwf_p': float(p_mw)})

with open(os.path.join(out_dir, 'cohort_lesion_size_mwf.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(summary[0].keys())); w.writeheader(); w.writerows(summary)

# verdict
sm = summary[0]
print('\n--- Small-lesion verdict ---')
if sm['mvf_p'] < 0.05 and sm['mwf_p'] >= 0.05:
    print(f'  Small lesions: MVF drops {sm["mvf_drop_pct"]:+.1f}% (p={sm["mvf_p"]:.3g}, significant) '
          f'but MWF {sm["mwf_drop_pct"]:+.1f}% (p={sm["mwf_p"]:.3g}, n.s.).')
    print('  -> MIMM detects small lesions that MWF misses: the resolution advantage.')
elif sm['mvf_p'] < 0.05 and sm['mwf_p'] < 0.05:
    print('  Both detect small lesions; no clear MIMM-only advantage at this size.')
else:
    print('  MVF does not significantly detect small lesions either; no advantage shown.')

# --- figure: drop vs volume, MVF vs MWF ---
fig, ax = plt.subplots(figsize=(8, 6))
order = np.argsort(vol)
ax.scatter(vol, dmvf, s=18, c='#1f77b4', alpha=0.5, label='MIMM MVF')
ax.scatter(vol, dmwf, s=18, c='#d62728', alpha=0.5, label='MWF')
xb = [np.nanmean(vol[sel]) for _, sel in bins]
ax.plot(xb, [s['mvf_drop_pct'] for s in summary], '-o', color='#1f77b4', lw=2, label='MVF (binned)')
ax.plot(xb, [s['mwf_drop_pct'] for s in summary], '-o', color='#d62728', lw=2, label='MWF (binned)')
ax.axhline(0, color='k', lw=0.8)
ax.set_xscale('log'); ax.set_xlabel('lesion volume (mL, log scale)')
ax.set_ylabel('drop vs perilesional NAWM (%)')
ax.set_title(f'Lesion detection vs size: MIMM MVF vs MWF\n'
             f'{n_conf} FLAIR-confirmed lesions ({n_rej_flair} rejected as not hyperintense)')
ax.legend(); ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'cohort_lesion_size_mwf.png'), dpi=150)
print(f'\nsaved: cohort_lesion_size_mwf.png / .csv')
