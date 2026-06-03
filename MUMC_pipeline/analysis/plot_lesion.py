"""
Lesion analysis for MIMM pipeline — dormant until lesion data arrives.

Requires lesion/lesion_mask.nii.gz (registered to ME-GRE by register_flair.sh).
If absent, exits cleanly.

Compares MS lesions against normal-appearing white matter (NAWM) and, where
T2-GRASE MWF is available, validates demyelination within lesions.

Figures:
  46_lesion_overlay.png       lesion locations on anatomy + MVF (3 views)
  47_lesion_vs_NAWM.png       voxel distributions, lesion vs NAWM, key maps
  48_lesion_MVF_vs_MWF.png    per-lesion scatter MVF vs MWF (if MWF present)
  49_lesion_summary.csv       per-lesion table (mean of each map per lesion)
"""

import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib.pyplot as plt
from scipy import ndimage
from scipy.stats import pearsonr
import os

_od = os.environ.get('MIMM_OUTPUT_DIR')
if _od:
    subj_dir = _od
    out_dir  = os.path.join(_od, 'figures')
    ANALYSIS_DIR = os.path.join(_od, 'analysis')
else:
    try:
        from paths import OUTPUT_DIR as subj_dir, FIG_DIR as out_dir, ANALYSIS_DIR
    except ImportError:
        raise SystemExit('Copy MUMC_pipeline/analysis/paths_template.py to paths.py.')
os.makedirs(out_dir, exist_ok=True)

# ── Dormant guard ─────────────────────────────────────────────────────────────
lesion_path = os.path.join(subj_dir, 'lesion', 'lesion_mask.nii.gz')
if not os.path.exists(lesion_path):
    print('lesion_mask.nii.gz not found — lesion data not yet available.')
    print(f'  Expected: {lesion_path}')
    print('  Run lesion/run_lst.m then lesion/register_flair.sh first.')
    raise SystemExit(0)

print(f'Lesion mask found: {lesion_path}')

def load(p):
    return np.array(nib.load(p).dataobj).astype(np.float32)

# ── Load masks and maps ───────────────────────────────────────────────────────
lesion = load(lesion_path) > 0.5
brain  = load(os.path.join(subj_dir, 'qsm', 'brain_mask.nii.gz')) > 0
fa     = load(os.path.join(subj_dir, 'atlas', 'FA_atlas.nii.gz'))
mag    = load(os.path.join(subj_dir, 'qsm', 'mag_e1.nii.gz'))

lesion = lesion & brain
wm     = brain & (fa > 0.20)
nawm   = wm & ~lesion          # normal-appearing WM
n_les  = int(lesion.sum())
print(f'  {n_les} lesion voxels, {int(nawm.sum())} NAWM voxels')

if n_les < 5:
    print('  Too few lesion voxels for analysis — exiting.')
    raise SystemExit(0)

# Key maps (name → (volume, unit)). chi_neg shown as |.| (diamagnetic).
maps = {
    'MVF_basic':  (load(os.path.join(subj_dir, 'mimm', 'MVF_basic.nii.gz')),       'fraction'),
    'g_ratio':    (load(os.path.join(subj_dir, 'mimm', 'g_ratio_basic.nii.gz')),   '-'),
    '|chi_neg|':  (np.abs(load(os.path.join(subj_dir, 'chisep', 'chi_neg.nii.gz'))), 'ppm'),
    'chi_pos':    (load(os.path.join(subj_dir, 'chisep', 'chi_pos.nii.gz')),       'ppm'),
    'R2s':        (load(os.path.join(subj_dir, 'qsm', 'R2star.nii.gz')),           's⁻¹'),
}
mwf_path = os.path.join(subj_dir, 'grase', 'MWF.nii.gz')
has_mwf = os.path.exists(mwf_path)
if has_mwf:
    maps['MWF'] = (load(mwf_path), 'fraction')

# ── Figure 46: lesion overlay on anatomy + MVF ───────────────────────────────
# Slices chosen at the most lesion-dense plane in each view, so lesions show.
def densest(axis):
    counts = lesion.sum(axis=tuple(i for i in range(3) if i != axis))
    return int(np.argmax(counts))
cx, cy, cz = densest(0), densest(1), densest(2)

def sl(vol, view):
    if view == 0:  return np.rot90(vol[:, :, cz])
    if view == 1:  return np.rot90(vol[:, cy, :])
    return np.rot90(vol[cx, :, :])

mvf_b = maps['MVF_basic'][0]
# Per-row display range: anatomy uses 1-99th percentile of the whole brain;
# MVF uses a fixed 0-0.5 scale.
mag_lo, mag_hi = np.percentile(mag[brain], [1, 99]) if brain.any() else (0, 1)
fig, axes = plt.subplots(2, 3, figsize=(15, 10), facecolor='black')
fig.suptitle(f'MS Lesions in ME-GRE space  ({n_les} voxels, {n_les/1000:.1f} mL)',
             color='white', fontsize=15, fontweight='bold', y=1.0)
view_names = ['Axial', 'Coronal', 'Sagittal']
rows = [('Anatomy + lesion outline', mag,   'gray', mag_lo, mag_hi),
        ('MVF + lesion outline',     mvf_b, 'hot',  0.0,    0.5)]
for r, (label, vol, cmap, vmin, vmax) in enumerate(rows):
    for v in range(3):
        a = axes[r, v]; a.set_facecolor('black')
        a.imshow(sl(vol, v), cmap=cmap, vmin=vmin, vmax=vmax,
                 aspect='auto', interpolation='nearest')
        les_sl = sl(lesion.astype(float), v)
        if les_sl.any():
            a.contour(les_sl, levels=[0.5], colors='cyan', linewidths=1.0)
        if r == 0:
            a.set_title(view_names[v], color='white', fontsize=11)
        a.axis('off')
    axes[r, 0].set_ylabel(label, color='white', fontsize=11)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, '46_lesion_overlay.png'),
            dpi=150, bbox_inches='tight', facecolor='black')
plt.close(); print('Saved: 46_lesion_overlay.png')

# ── Figure 47: lesion vs NAWM voxel distributions ────────────────────────────
# Subsample NAWM for fast violins; lesions usually small enough to use whole.
rng = np.random.default_rng(0)
def sample(mask, vol, n=20000):
    vals = vol[mask]; vals = vals[np.isfinite(vals)]
    if len(vals) > n:
        vals = rng.choice(vals, n, replace=False)
    return vals

n_maps = len(maps)
fig, axes = plt.subplots(1, n_maps, figsize=(3.2 * n_maps, 5), facecolor='#0d0d0d')
if n_maps == 1: axes = [axes]
fig.suptitle('Lesion vs NAWM — voxel distributions', color='white',
             fontsize=14, fontweight='bold')

summary = []
for ax, (name, (vol, unit)) in zip(axes, maps.items()):
    nawm_vals = sample(nawm, vol)
    les_vals  = vol[lesion]; les_vals = les_vals[np.isfinite(les_vals)]
    parts = ax.violinplot([nawm_vals, les_vals], showmeans=True, showextrema=False)
    for pc, col in zip(parts['bodies'], ['#5c9ee0', '#e05c5c']):
        pc.set_facecolor(col); pc.set_alpha(0.7)
    parts['cmeans'].set_color('white')
    ax.set_xticks([1, 2]); ax.set_xticklabels(['NAWM', 'Lesion'], color='white', fontsize=9)
    ax.set_title(f'{name}\n({unit})', color='white', fontsize=10)
    ax.tick_params(colors='white')
    for s in ax.spines.values(): s.set_edgecolor('#444444')
    ax.set_facecolor('#111111')
    summary.append({'map': name,
                    'NAWM_mean': float(np.mean(nawm_vals)),
                    'lesion_mean': float(np.mean(les_vals)),
                    'pct_change': float((np.mean(les_vals) - np.mean(nawm_vals)) /
                                        np.mean(nawm_vals) * 100) if np.mean(nawm_vals) else np.nan})
plt.tight_layout()
plt.savefig(os.path.join(out_dir, '47_lesion_vs_NAWM.png'),
            dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
plt.close(); print('Saved: 47_lesion_vs_NAWM.png')

print('\nLesion vs NAWM summary:')
print(pd.DataFrame(summary).to_string(index=False, float_format=lambda x: f'{x:.4f}'))

# ── Per-lesion table (connected components) ──────────────────────────────────
lesion_lbl, n_lesions = ndimage.label(lesion)
print(f'\n{n_lesions} discrete lesions found.')
rows = []
for li in range(1, n_lesions + 1):
    m = lesion_lbl == li
    nv = int(m.sum())
    if nv < 3:   # skip tiny specks
        continue
    row = {'lesion_id': li, 'n_voxels': nv}
    for name, (vol, _u) in maps.items():
        row[f'{name}_mean'] = float(np.nanmean(vol[m]))
    rows.append(row)
lesion_df = pd.DataFrame(rows)
lesion_csv = os.path.join(ANALYSIS_DIR, 'lesion_summary.csv')
os.makedirs(ANALYSIS_DIR, exist_ok=True)
lesion_df.to_csv(lesion_csv, index=False, float_format='%.5f')
print(f'Saved: {lesion_csv}  ({len(lesion_df)} lesions ≥3 voxels)')

# ── Figure 48: per-lesion MVF vs MWF (if MWF present) ────────────────────────
if has_mwf and len(lesion_df) >= 3:
    fig, ax = plt.subplots(figsize=(7, 7), facecolor='#0d0d0d')
    ax.set_facecolor('#111111')
    x = lesion_df['MWF_mean']; y = lesion_df['MVF_basic_mean']
    sizes = (lesion_df['n_voxels'] / lesion_df['n_voxels'].max() * 200 + 20)
    ax.scatter(x, y, s=sizes, color='#e05c5c', edgecolors='white',
               linewidths=0.4, alpha=0.8, zorder=3)
    if len(lesion_df) > 2:
        r, p = pearsonr(x, y)
        m, b = np.polyfit(x, y, 1)
        xl = np.linspace(x.min(), x.max(), 100)
        ax.plot(xl, m * xl + b, '--', color='cyan', lw=1.5,
                label=f'r = {r:.2f}, p = {p:.3f}')
        ax.legend(fontsize=9, facecolor='#1a1a1a', edgecolor='#444444', labelcolor='white')
    ax.set_xlabel('MWF — T2-GRASE (per-lesion mean)', color='white', fontsize=11)
    ax.set_ylabel('MVF — MIMM basic (per-lesion mean)', color='white', fontsize=11)
    ax.set_title('Demyelination within lesions: MVF vs MWF\n(dot size = lesion volume)',
                 color='white', fontsize=12, fontweight='bold')
    ax.tick_params(colors='white')
    for s in ax.spines.values(): s.set_edgecolor('#444444')
    ax.grid(color='#2a2a2a', lw=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '48_lesion_MVF_vs_MWF.png'),
                dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
    plt.close(); print('Saved: 48_lesion_MVF_vs_MWF.png')
else:
    print('Skipped: 48_lesion_MVF_vs_MWF.png (needs MWF + ≥3 lesions)')
