#!/usr/bin/env python3
"""
fig_slicemap.py -- a clean whole-slice map of the pipeline's own output for the deck.

One axial slice: FLAIR | MIMM MVF (myelin) | MWF (T2-GRASE). MVF and MWF are
DIFFERENT quantities (volume vs water fraction), so each is scaled to its own
range: the eye compares the spatial PATTERN across white matter (the visual
companion to the r=0.70 agreement) without implying the absolute values match.

The picker PREFERS a subject/slice that actually contains a visible lesion (most
lesion voxels on a still WM-rich slice) and outlines it in red on all three
panels -- this is the intro slide ("here's what the pipeline produces"), and the
visible lesion is the hook into the compartment-decomposition story (slide 11).
Falls back to the plain WM-richest slice (no outline) if no subject has a usable
lesion mask.

Usage:
  python3 fig_slicemap.py <results_dir>                 # auto subject + slice with a lesion
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

les_available = {}
if FORCE_Z is None:
    # among candidates, prefer whichever subject has the biggest visible lesion
    # on a slice that is ALSO reasonably WM-rich (so it still reads as a normal
    # brain slice, not a lesion close-up).
    best = None  # (lesion_area, sid, d, p, z)
    for sid, d, p in cands:
        if not os.path.exists(p['les']):
            continue
        les = load(p['les']); brain = load(p['brain'])
        fa = load(p['fa'])
        les = (les > 0.5) & (brain > 0)
        if les.sum() < 20:
            continue
        wm = (fa > 0.20) if fa is not None else np.ones_like(les, bool)
        # score each slice: lesion voxels, but only among slices with decent WM
        wm_per_slice = wm.sum(axis=(0, 1))
        wm_ok = wm_per_slice > np.percentile(wm_per_slice[wm_per_slice > 0], 40)
        les_per_slice = les.sum(axis=(0, 1)).astype(float)
        les_per_slice[~wm_ok] = 0
        z = int(np.argmax(les_per_slice))
        area = les_per_slice[z]
        if area < 10:
            continue
        les_available[sid] = True
        if best is None or area > best[0]:
            best = (area, sid, d, p, z)
    if best is not None:
        _, sid, d, p, z = best
        print(f'[pick] subject={sid}  z={z}  (slice with a visible lesion, area={best[0]:.0f} vox)')
    else:
        sid, d, p = cands[0]
        fa = load(p['fa']); brain = load(p['brain']); bm = brain > 0
        wm = (fa > 0.20) & bm if fa is not None else bm
        z = int(np.argmax(wm.sum(axis=(0, 1))))
        print(f'[pick] subject={sid}  z={z}  (no subject had a usable lesion mask; WM-richest slice, no outline)')
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

# lesion outline for this exact slice, if this subject has one there
les_s = None
if os.path.exists(p['les']):
    les3d = load(p['les'])
    les3d = (les3d > 0.5) & bm
    if les3d[:, :, z].sum() >= 5:
        les_s = les3d

# crop tightly to the brain on this slice
ys, xs = np.where(bm[:, :, z])
pad = 3
r0, r1 = max(ys.min() - pad, 0), min(ys.max() + pad, mvf.shape[0])
c0, c1 = max(xs.min() - pad, 0), min(xs.max() + pad, mvf.shape[1])


def sl(v):
    return np.rot90(v[r0:r1, c0:c1, z])


mvf_s, mwf_s, fl_s = sl(mvf), sl(mwf), sl(flair)
les_out = np.rot90(les_s[r0:r1, c0:c1, z].astype(float)) if les_s is not None else None

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
