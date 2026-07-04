#!/usr/bin/env python3
"""
fig_brainmap.py  --  a real brain slice for the defence deck.

One representative lesion, shown as FLAIR | MIMM MVF | MIMM AVF | MWF, with the
lesion outlined in red on every panel. MVF and MWF share a colour scale so they
are directly comparable. The point (slide 10): inside the lesion BOTH the myelin
(MVF) and the axon (AVF) compartment darken -- the tissue is atrophic, not just
demyelinated -- and MIMM's MVF broadly tracks the independent MWF in the NAWM.

All maps are already co-registered in ME-GRE space by the pipeline, so this is a
pure rendering of the pipeline's own output (no numbers change).

Usage:
  python3 fig_brainmap.py <results_dir>                 # auto-pick best subject/slice
  python3 fig_brainmap.py <results_dir> <subject> <z>   # force a subject and axial slice
Output: <results_dir>/cohort_analysis/fig_brainmap.png (+ .svg)
"""
import sys, os, glob
import numpy as np
import nibabel as nib
from scipy import ndimage
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mimm_style import apply_style, C
apply_style()

if len(sys.argv) < 2:
    sys.exit('usage: fig_brainmap.py <results_dir> [subject] [z]')
RES = sys.argv[1]
FORCE_SUBJ = sys.argv[2] if len(sys.argv) > 2 else None
FORCE_Z    = int(sys.argv[3]) if len(sys.argv) > 3 else None
OUT = os.path.join(RES, 'cohort_analysis'); os.makedirs(OUT, exist_ok=True)


def load(p):
    return np.asarray(nib.load(p).dataobj, dtype=float) if os.path.exists(p) else None


def subj_paths(d):
    return dict(mvf=os.path.join(d, 'mimm', 'MVF_Atlas.nii.gz'),
                fvf=os.path.join(d, 'mimm', 'FVF_Atlas.nii.gz'),
                mwf=os.path.join(d, 'grase', 'MWF.nii.gz'),
                les=os.path.join(d, 'lesion', 'lesion_mask.nii.gz'),
                flair=os.path.join(d, 'lesion', 'FLAIR_mgre.nii.gz'),
                fa=os.path.join(d, 'atlas', 'FA_atlas.nii.gz'),
                brain=os.path.join(d, 'qsm', 'brain_mask.nii.gz'))


# We want a DISCRETE white-matter lesion surrounded by normal WM, not a big
# confluent periventricular blob hugging the ventricles (where the dark CSF
# reads as "myelin loss"). Score each connected lesion by area in a teaching-
# friendly range AND by how much of its surrounding ring is real white matter.
AREA_LO, AREA_HI = 25, 450        # voxels on the shown slice: discrete, not confluent


def score_candidates(les3d, fa3d):
    """yield (score, z, area, wm_ring_frac, comp2d) for each connected lesion."""
    lab, n = ndimage.label(les3d)
    out = []
    for k in range(1, n + 1):
        comp = lab == k
        if comp.sum() < 20:
            continue
        z = int(np.argmax(comp.sum(axis=(0, 1))))
        m2 = comp[:, :, z]
        area = float(m2.sum())
        if not (AREA_LO <= area <= AREA_HI):
            continue
        ring = ndimage.binary_dilation(m2, iterations=6) & ~ndimage.binary_dilation(m2, iterations=2)
        fa2 = fa3d[:, :, z] if fa3d is not None else None
        wm = float(np.mean(fa2[ring] > 0.20)) if (fa2 is not None and ring.sum()) else 0.0
        # reward discrete WM-surrounded lesions; wm_ring_frac dominates
        out.append((area * (wm ** 2), z, area, wm, m2.copy()))
    return out


cands = []
for ld in sorted(glob.glob(os.path.join(RES, '*', 'lesion', 'lesion_mask.nii.gz'))):
    d = os.path.dirname(os.path.dirname(ld)); sid = os.path.basename(d)
    if FORCE_SUBJ and sid != FORCE_SUBJ:
        continue
    p = subj_paths(d)
    if not all(os.path.exists(p[k]) for k in ('mvf', 'fvf', 'mwf', 'les', 'flair')):
        continue
    les = (load(p['les']) > 0.5)
    if les.sum() < 30:
        continue
    fa = load(p['fa'])
    for score, z, area, wm, comp2d in score_candidates(les, fa):
        cands.append((score, sid, d, z, area, wm, comp2d))

if not cands:
    sys.exit('no discrete WM lesion found (try relaxing AREA_LO/AREA_HI)')
cands.sort(key=lambda c: c[0], reverse=True)
print('[shortlist] top discrete WM lesions (subject, z, area_vox, wm_ring_frac):')
for c in cands[:8]:
    print(f'    {c[1]}  z={c[3]}  area={c[4]:.0f}  wm_ring={c[5]:.2f}')

# choose: if a z is forced, prefer the scored component on that exact slice
chosen = cands[0]
if FORCE_Z is not None:
    match = [c for c in cands if c[3] == FORCE_Z]
    chosen = match[0] if match else cands[0]
_, sid, d, z, _, _, foc = chosen
if FORCE_Z is not None:
    z = FORCE_Z
print(f'[pick] subject={sid}  z={z}  (one lesion isolated)')

p = subj_paths(d)
mvf = load(p['mvf']); fvf = load(p['fvf']); mwf = load(p['mwf'])
avf = fvf - mvf
les = (load(p['les']) > 0.5); flair = load(p['flair'])
brain = load(p['brain'])
if brain is not None:
    m = brain > 0
    for v in (mvf, avf, mwf):
        v[~m] = np.nan

# isolate ONE lesion component on this slice (not every lesion on the slice), so
# the outline and zoom show a single discrete lesion.
if foc is None or foc.shape != les[:, :, z].shape or not (les[:, :, z] & foc).any():
    lab2, n2 = ndimage.label(les[:, :, z])
    if n2 == 0:
        z = int(np.argmax(les.sum(axis=(0, 1)))); lab2, n2 = ndimage.label(les[:, :, z])
    foc = lab2 == (1 + int(np.argmax([(lab2 == k).sum() for k in range(1, n2 + 1)])))
    print(f'[info] isolated the largest lesion component on z={z}')

# --- tight bounding box around the single focused lesion, with margin ---
ys, xs = np.where(foc)
pad = 16
r0, r1 = max(ys.min() - pad, 0), min(ys.max() + pad, mvf.shape[0])
c0, c1 = max(xs.min() - pad, 0), min(xs.max() + pad, mvf.shape[1])


def sl(v):
    return np.rot90(v[r0:r1, c0:c1, z])   # radiological-ish display


les_s = sl(foc.astype(float))   # outline only the single focused lesion

# shared scale for the two myelin maps so MVF vs MWF is a fair visual comparison
mvf_s, mwf_s, avf_s, fl_s = sl(mvf), sl(mwf), sl(avf), sl(flair)
vmax_my = np.nanpercentile(np.concatenate([mvf_s[np.isfinite(mvf_s)],
                                           mwf_s[np.isfinite(mwf_s)]]), 98)
vmax_ax = np.nanpercentile(avf_s[np.isfinite(avf_s)], 98)

panels = [
    ('FLAIR',            fl_s,  'gray',  (np.nanpercentile(fl_s, 2), np.nanpercentile(fl_s, 98)), None),
    ('MIMM  MVF (myelin)', mvf_s, 'magma', (0, vmax_my), 'fraction'),
    ('MIMM  AVF (axon)',   avf_s, 'viridis', (0, vmax_ax), 'fraction'),
    ('MWF  (T2-GRASE)',    mwf_s, 'magma', (0, vmax_my), 'fraction'),
]

fig, axes = plt.subplots(1, 4, figsize=(14.5, 4.6))
for ax, (title, img, cmap, (lo, hi), cblab) in zip(axes, panels):
    im = ax.imshow(img, cmap=cmap, vmin=lo, vmax=hi, interpolation='nearest')
    ax.contour(les_s, levels=[0.5], colors=[C['highlight']], linewidths=1.6)
    ax.set_title(title, color=C['text'], fontsize=12, pad=6)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    if cblab:
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        cb.ax.tick_params(labelsize=8)

fig.suptitle('One lesion, four maps (representative patient) -- '
             'myelin and axon both drop inside the lesion (red)',
             color=C['text'], fontsize=13.5)
fig.tight_layout(rect=[0, 0, 1, 0.93])
for ext in ('png', 'svg'):
    fig.savefig(os.path.join(OUT, f'fig_brainmap.{ext}'), dpi=200, facecolor='white')
print('saved', os.path.join(OUT, 'fig_brainmap.png'), '(+ .svg)')

# quick numbers to quote: mean inside THIS lesion vs its perilesional ring
ring = (ndimage.binary_dilation(foc, iterations=5) & ~ndimage.binary_dilation(foc, iterations=2)
        & ~les[:, :, z])   # ring excludes any other lesion tissue
for name, vol in [('MVF', mvf[:, :, z]), ('AVF', avf[:, :, z]), ('MWF', mwf[:, :, z])]:
    li = np.nanmean(vol[foc]); ne = np.nanmean(vol[ring])
    print(f'  {name}: lesion={li:.3f}  peri-NAWM={ne:.3f}  ({100*(li-ne)/ne:+.1f}%)')
