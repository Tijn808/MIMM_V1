#!/usr/bin/env python3
"""
fig_slicemap.py -- a clean whole-slice map of the pipeline's own output for the deck.

One axial slice: FLAIR | MIMM MVF (myelin) | MWF (T2-GRASE). MVF and MWF are
DIFFERENT quantities (volume vs water fraction), so each is scaled to its own
range: the eye compares the spatial PATTERN across white matter (the visual
companion to the r=0.70 agreement) without implying the absolute values match.

The picker PREFERS a subject/slice that contains ONE discrete, WM-surrounded
lesion (same connected-component + low-CSF-adjacency scoring as fig_brainmap.py,
NOT just "most lesion voxels on the slice" -- that would grab a whole confluent
periventricular disease burden and outline it as if it were one lesion, which
reads as tracing the ventricles). Only that single component is outlined in red,
on all three panels -- this is the intro slide ("here's what the pipeline
produces"), and the visible lesion is the hook into the compartment-decomposition
story (slide 11). Falls back to the plain WM-richest slice (no outline) if no
subject has a usable, discrete, non-confluent lesion.

Usage:
  python3 fig_slicemap.py <results_dir>                 # auto subject + slice with a lesion
  python3 fig_slicemap.py <results_dir> <subject> <z>
Output: <results_dir>/cohort_analysis/fig_slicemap.png (+ .svg)
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

# discrete-lesion scoring (same idea as fig_brainmap.py): reward a moderate-size,
# WM-surrounded lesion and heavily penalise ventricle/CSF adjacency, so we never
# outline a confluent periventricular disease burden as if it were "a lesion".
AREA_LO, AREA_HI = 25, 450


def best_discrete_lesion(les3d, fa3d, bm3d):
    """Return (score, z, area, comp2d) for the best single connected lesion
    component across the whole volume, or None if none qualifies."""
    lab, n = ndimage.label(les3d)
    best = None
    for k in range(1, n + 1):
        comp = lab == k
        if comp.sum() < 15:
            continue
        z = int(np.argmax(comp.sum(axis=(0, 1))))
        m2 = comp[:, :, z]
        area = float(m2.sum())
        if not (AREA_LO <= area <= AREA_HI):
            continue
        ring = ndimage.binary_dilation(m2, iterations=6) & ~ndimage.binary_dilation(m2, iterations=2)
        fa2 = fa3d[:, :, z] if fa3d is not None else None
        if fa2 is not None and ring.sum():
            wm = float(np.mean(fa2[ring] > 0.20))
            csf = float(np.mean(fa2[ring] < 0.10))
        else:
            wm, csf = 0.0, 1.0
        score = area * (wm ** 2) * ((1 - csf) ** 2)
        if best is None or score > best[0]:
            best = (score, z, area, m2.copy())
    return best

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
                brain=os.path.join(d, 'qsm', 'brain_mask.nii.gz'),
                les=os.path.join(d, 'lesion', 'lesion_mask.nii.gz'))


REQUIRED = ('mvf', 'mwf', 'flair', 'brain')

# candidate subjects: everything with the required maps present
cands = []
for md in sorted(glob.glob(os.path.join(RES, '*', 'mimm', 'MVF_Atlas.nii.gz'))):
    d = os.path.dirname(os.path.dirname(md)); sid = os.path.basename(d)
    if FORCE_SUBJ and sid != FORCE_SUBJ:
        continue
    p = paths(d)
    if all(os.path.exists(p[k]) for k in REQUIRED):
        cands.append((sid, d, p))
if not cands:
    sys.exit('no subject with MVF+MWF+FLAIR+brain present')

foc = None  # the single lesion component to outline (2D mask on slice z), if any
if FORCE_Z is None:
    # among candidates, prefer whichever subject has the best DISCRETE lesion
    # (moderate area, WM-surrounded, low CSF adjacency) -- never the whole raw
    # mask, which can be several lesions merged/confluent around the ventricles.
    best = None  # (score, sid, d, p, z, comp2d)
    for sid, d, p in cands:
        if not os.path.exists(p['les']):
            continue
        les = load(p['les']); brain = load(p['brain']); fa = load(p['fa'])
        bm_ = brain > 0
        les = (les > 0.5) & bm_
        if les.sum() < 15:
            continue
        r = best_discrete_lesion(les, fa, bm_)
        if r is None:
            continue
        score, z_, area, comp2d = r
        if best is None or score > best[0]:
            best = (score, sid, d, p, z_, comp2d, area)
    if best is not None:
        _, sid, d, p, z, foc, area = best
        print(f'[pick] subject={sid}  z={z}  (discrete lesion, area={area:.0f} vox)')
    else:
        sid, d, p = cands[0]
        fa = load(p['fa']); brain = load(p['brain']); bm = brain > 0
        wm = (fa > 0.20) & bm if fa is not None else bm
        z = int(np.argmax(wm.sum(axis=(0, 1))))
        print(f'[pick] subject={sid}  z={z}  (no subject had a usable discrete lesion; WM-richest slice, no outline)')
else:
    sid, d, p = cands[0]
    z = FORCE_Z
    print(f'[pick] subject={sid}  z={z}  (forced)')

mvf = load(p['mvf']); mwf = load(p['mwf']); flair = load(p['flair'])
fa = load(p['fa']); brain = load(p['brain'])
bm = brain > 0
mwf = np.clip(mwf, 0, 0.5)
for v in (mvf, mwf):
    v[~bm] = np.nan

les_s = foc  # already the single 2D component mask on slice z, or None

# crop tightly to the brain on this slice
ys, xs = np.where(bm[:, :, z])
pad = 3
r0, r1 = max(ys.min() - pad, 0), min(ys.max() + pad, mvf.shape[0])
c0, c1 = max(xs.min() - pad, 0), min(xs.max() + pad, mvf.shape[1])


def sl(v):
    return np.rot90(v[r0:r1, c0:c1, z])


mvf_s, mwf_s, fl_s = sl(mvf), sl(mwf), sl(flair)
# les_s is already a single-component 2D mask (full slice size) on slice z
les_out = np.rot90(les_s[r0:r1, c0:c1].astype(float)) if les_s is not None else None

# MVF and MWF are DIFFERENT quantities (volume vs water fraction), so each map is
# scaled to its OWN range (2nd-98th pct). This compares the spatial PATTERN without
# implying the absolute values are equivalent, and lets each use its full range.
vmax_mvf = np.nanpercentile(mvf_s[np.isfinite(mvf_s)], 98)
vmax_mwf = np.nanpercentile(mwf_s[np.isfinite(mwf_s)], 98)

panels = [
    ('FLAIR',              fl_s,  'gray',  (np.nanpercentile(fl_s, 2), np.nanpercentile(fl_s, 98)), None),
    ('MIMM  MVF (myelin)', mvf_s, 'magma', (0, vmax_mvf), 'fraction'),
    ('MWF  (T2-GRASE)',    mwf_s, 'magma', (0, vmax_mwf), 'fraction'),
]

fig, axes = plt.subplots(1, 3, figsize=(12, 4.8))
for ax, (title, img, cmap, (lo, hi), cblab) in zip(axes, panels):
    im = ax.imshow(img, cmap=cmap, vmin=lo, vmax=hi, interpolation='nearest')
    if les_out is not None:
        ax.contour(les_out, levels=[0.5], colors=[C['highlight']], linewidths=1.6)
    ax.set_title(title, color=C['text'], fontsize=12, pad=6)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    if cblab:
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        cb.ax.tick_params(labelsize=8)

suptitle = 'MIMM myelin map vs the independent MWF reference (one axial slice)'
if les_out is not None:
    suptitle += ' -- lesion outlined in red'
fig.suptitle(suptitle, color=C['text'], fontsize=13.5)
fig.tight_layout(rect=[0, 0, 1, 0.94])
for ext in ('png', 'svg'):
    fig.savefig(os.path.join(OUT, f'fig_slicemap.{ext}'), dpi=200, facecolor='white')
print('saved', os.path.join(OUT, 'fig_slicemap.png'), '(+ .svg)',
      '-- lesion outlined' if les_out is not None else '-- no lesion on this slice')
