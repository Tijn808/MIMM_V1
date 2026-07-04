#!/usr/bin/env python3
"""
fig_slicemap.py -- a clean whole-slice map of the pipeline's own output for the deck.

One axial slice: FLAIR | MIMM MVF (myelin) | MWF (T2-GRASE). MVF and MWF share a
colour scale so the eye sees they broadly agree across white matter (the visual
companion to the r=0.70 agreement). No lesion outline -- this simply shows what
the pipeline produces on a recognisable brain, the brain image the deck lacks.

Usage:
  python3 fig_slicemap.py <results_dir>                 # auto subject + WM-rich slice
  python3 fig_slicemap.py <results_dir> <subject> <z>
Output: <results_dir>/cohort_analysis/fig_slicemap.png (+ .svg)
"""
import sys, os, glob
import numpy as np
import nibabel as nib
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mimm_style import apply_style, C
apply_style()

if len(sys.argv) < 2:
    sys.exit('usage: fig_slicemap.py <results_dir> [subject] [z]')
RES = sys.argv[1]
FORCE_SUBJ = sys.argv[2] if len(sys.argv) > 2 else None
FORCE_Z = int(sys.argv[3]) if len(sys.argv) > 3 else None
OUT = os.path.join(RES, 'cohort_analysis'); os.makedirs(OUT, exist_ok=True)


def load(p):
    return np.asarray(nib.load(p).dataobj, dtype=float) if os.path.exists(p) else None


def paths(d):
    return dict(mvf=os.path.join(d, 'mimm', 'MVF_Atlas.nii.gz'),
                mwf=os.path.join(d, 'grase', 'MWF.nii.gz'),
                fa=os.path.join(d, 'atlas', 'FA_atlas.nii.gz'),
                flair=os.path.join(d, 'lesion', 'FLAIR_mgre.nii.gz'),
                brain=os.path.join(d, 'qsm', 'brain_mask.nii.gz'))


# pick subject: forced, else first with MVF+MWF+FLAIR+brain all present
subj = None
for md in sorted(glob.glob(os.path.join(RES, '*', 'mimm', 'MVF_Atlas.nii.gz'))):
    d = os.path.dirname(os.path.dirname(md)); sid = os.path.basename(d)
    if FORCE_SUBJ and sid != FORCE_SUBJ:
        continue
    p = paths(d)
    if all(os.path.exists(p[k]) for k in ('mvf', 'mwf', 'flair', 'brain')):
        subj = (sid, d, p); break
if subj is None:
    sys.exit('no subject with MVF+MWF+FLAIR+brain present')
sid, d, p = subj

mvf = load(p['mvf']); mwf = load(p['mwf']); flair = load(p['flair'])
fa = load(p['fa']); brain = load(p['brain'])
bm = brain > 0
mwf = np.clip(mwf, 0, 0.5)
for v in (mvf, mwf):
    v[~bm] = np.nan

# pick slice: WM-rich (most FA>0.2 voxels), else forced
if FORCE_Z is not None:
    z = FORCE_Z
else:
    wm = (fa > 0.20) & bm if fa is not None else bm
    z = int(np.argmax(wm.sum(axis=(0, 1))))
print(f'[pick] subject={sid}  z={z}')

# crop tightly to the brain on this slice
ys, xs = np.where(bm[:, :, z])
pad = 3
r0, r1 = max(ys.min() - pad, 0), min(ys.max() + pad, mvf.shape[0])
c0, c1 = max(xs.min() - pad, 0), min(xs.max() + pad, mvf.shape[1])


def sl(v):
    return np.rot90(v[r0:r1, c0:c1, z])


mvf_s, mwf_s, fl_s = sl(mvf), sl(mwf), sl(flair)
vmax = np.nanpercentile(np.concatenate([mvf_s[np.isfinite(mvf_s)],
                                        mwf_s[np.isfinite(mwf_s)]]), 98)

panels = [
    ('FLAIR',              fl_s,  'gray',  (np.nanpercentile(fl_s, 2), np.nanpercentile(fl_s, 98)), None),
    ('MIMM  MVF (myelin)', mvf_s, 'magma', (0, vmax), 'fraction'),
    ('MWF  (T2-GRASE)',    mwf_s, 'magma', (0, vmax), 'fraction'),
]

fig, axes = plt.subplots(1, 3, figsize=(12, 4.8))
for ax, (title, img, cmap, (lo, hi), cblab) in zip(axes, panels):
    im = ax.imshow(img, cmap=cmap, vmin=lo, vmax=hi, interpolation='nearest')
    ax.set_title(title, color=C['text'], fontsize=12, pad=6)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    if cblab:
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        cb.ax.tick_params(labelsize=8)

fig.suptitle('MIMM myelin map vs the independent MWF reference (one axial slice)',
             color=C['text'], fontsize=13.5)
fig.tight_layout(rect=[0, 0, 1, 0.94])
for ext in ('png', 'svg'):
    fig.savefig(os.path.join(OUT, f'fig_slicemap.{ext}'), dpi=200, facecolor='white')
print('saved', os.path.join(OUT, 'fig_slicemap.png'), '(+ .svg)')
