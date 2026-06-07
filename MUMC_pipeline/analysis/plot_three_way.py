"""
Three-way myelin-marker comparison — chi-separation |χ⁻|, MIMM atlas MVF, and
T2-GRASE MWF.

This is the headline cross-method figure set: it puts the three independent
myelin estimates side by side. Two are imaging-model based (MIMM atlas, chi-sep
diamagnetic susceptibility) and the third (MWF) is the closest thing to a myelin
gold standard available in vivo.

Dormant until T2-GRASE data arrives: requires grase/MWF.nii.gz (written by
ingest_grase.m). If absent, exits cleanly. MIMM and chi-sep maps are always
present once the MATLAB chain has run.

Outputs (figures 49–51):
  49_three_way_spatial.png   — MVF atlas | |χ⁻| | MWF spatial maps (own scales)
  50_three_way_scatter.png   — 3 pairwise ROI scatters, MWF as the reference axis
  51_three_way_ranking.png   — per-ROI grouped bars, each marker min-max normalised
                               (puts the three different units on one 0–1 axis so
                                tract ranking can be compared directly)
"""

import os
import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

_od = os.environ.get('MIMM_OUTPUT_DIR')
if _od:
    subj_dir     = _od
    out_dir      = os.path.join(_od, 'figures')
    ANALYSIS_DIR = os.path.join(_od, 'analysis')
else:
    try:
        from paths import OUTPUT_DIR as subj_dir, FIG_DIR as out_dir, ANALYSIS_DIR
    except ImportError:
        raise SystemExit('Copy MUMC_pipeline/analysis/paths_template.py to paths.py.')
os.makedirs(out_dir, exist_ok=True)

# ── Dormant guard ────────────────────────────────────────────────────────────
mwf_path = os.path.join(subj_dir, 'grase', 'MWF.nii.gz')
if not os.path.exists(mwf_path):
    print('MWF.nii.gz not found — T2-GRASE data not yet ingested.')
    print(f'  Expected: {mwf_path}')
    print('  Three-way comparison is dormant until then.')
    raise SystemExit(0)

print(f'T2-GRASE MWF found: {mwf_path}')


def load(path):
    return np.array(nib.load(path).dataobj).astype(np.float32)


brain = load(os.path.join(subj_dir, 'qsm', 'brain_mask.nii.gz')) > 0
mag   = load(os.path.join(subj_dir, 'qsm', 'mag_e1.nii.gz'))
fa    = load(os.path.join(subj_dir, 'atlas', 'FA_atlas.nii.gz'))
wm    = brain & (fa > 0.20)

mvf_a = load(os.path.join(subj_dir, 'mimm',   'MVF_Atlas.nii.gz'))
chineg = np.abs(load(os.path.join(subj_dir, 'chisep', 'chi_neg.nii.gz')))  # ppm, diamagnetic
mwf   = np.clip(load(mwf_path), 0, 0.5)
for v in (mvf_a, chineg, mwf):
    v[~brain] = 0

# Consistent colours per marker, reused across all three figures.
C_MVF, C_CHI, C_MWF = '#5c9ee0', '#c060d0', '#5cb85c'   # blue, purple, green


# ── Figure 49: spatial maps (each on its own physical scale) ─────────────────
def _slices(vol):
    coords = np.argwhere(vol > 0)
    cx, cy, cz = coords.mean(axis=0).astype(int)
    return np.rot90(vol[:, :, cz]), np.rot90(vol[:, cy, :]), np.rot90(vol[cx, :, :])


bg_slices = _slices(mag)
maps49 = [('MIMM atlas MVF', np.where(wm, mvf_a, 0),  'hot',     0, 0.50, 'fraction'),
          (r'chi-sep $|\chi^-|$', np.where(wm, chineg, 0), 'viridis', 0, 0.10, 'ppm'),
          ('T2-GRASE MWF',   np.where(brain, mwf, 0), 'hot',     0, 0.30, 'fraction')]
view_labels = ['Axial', 'Coronal', 'Sagittal']

fig, axes = plt.subplots(3, 3, figsize=(14, 12), facecolor='black')
fig.suptitle('Three myelin markers — MIMM atlas MVF | chi-sep |χ⁻| | T2-GRASE MWF',
             color='white', fontsize=15, fontweight='bold')
for r, (label, vol, cmap, vmin, vmax, unit) in enumerate(maps49):
    sl = _slices(vol)
    for c in range(3):
        a = axes[r, c]; a.set_facecolor('black')
        bg = bg_slices[c]
        a.imshow((bg - bg.min()) / (bg.max() - bg.min() + 1e-8),
                 cmap='gray', aspect='auto', interpolation='nearest')
        im = a.imshow(np.ma.masked_where(sl[c] == 0, sl[c]),
                      cmap=cmap, vmin=vmin, vmax=vmax,
                      aspect='auto', interpolation='nearest', alpha=0.85)
        if r == 0:
            a.set_title(view_labels[c], color='white', fontsize=11)
        a.axis('off')
    axes[r, 0].set_ylabel(label, color='white', fontsize=11, rotation=90, labelpad=8)
    cb = fig.colorbar(im, ax=axes[r, :], fraction=0.020, pad=0.01)
    cb.set_label(unit, color='white', fontsize=9)
    cb.ax.tick_params(colors='white', labelsize=8)

plt.savefig(os.path.join(out_dir, '49_three_way_spatial.png'),
            dpi=150, bbox_inches='tight', facecolor='black')
plt.close(); print('Saved: 49_three_way_spatial.png')


# ── ROI-level figures need roi_stats.csv with all three markers ──────────────
roi_csv = os.path.join(ANALYSIS_DIR, 'roi_stats.csv')
need = {'MVF_atlas_mean', 'chi_neg_chisep_mean', 'MWF_mean'}
if not (os.path.exists(roi_csv) and need.issubset(pd.read_csv(roi_csv, nrows=0).columns)):
    print('Skipped: 50/51 (roi_stats.csv lacks MWF — run extract_roi_stats.py with GRASE data)')
    raise SystemExit(0)

df = pd.read_csv(roi_csv).copy()
df['chineg'] = df['chi_neg_chisep_mean'].abs()
sem = lambda col: df[col] / np.sqrt(df['n_voxels'])
df['MVF_sem'] = sem('MVF_atlas_sd')
df['chineg_sem'] = sem('chi_neg_chisep_sd')
df['MWF_sem'] = sem('MWF_sd')

valid = df.dropna(subset=['MVF_atlas_mean', 'chineg', 'MWF_mean'])
valid = valid[valid['n_voxels'] > 0]


def _panel(ax, x, y, xerr, yerr, xlabel, ylabel, colour, identity=False):
    m = np.isfinite(x) & np.isfinite(y)
    r, p = pearsonr(x[m], y[m])
    sl, b = np.polyfit(x[m], y[m], 1)
    ax.set_facecolor('#111111')
    ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o', color=colour, markersize=5,
                alpha=0.85, ecolor='#777777', elinewidth=0.7, capsize=2)
    xf = np.linspace(np.nanmin(x), np.nanmax(x), 100)
    ax.plot(xf, sl * xf + b, '--', color='cyan', lw=1.5,
            label=f'r = {r:.2f}, p = {p:.1e}')
    if identity:
        lo = 0
        hi = max(np.nanmax(x), np.nanmax(y)) * 1.05
        ax.plot([lo, hi], [lo, hi], '-', color='#555555', lw=1.1, label='identity')
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel(xlabel, color='white', fontsize=10)
    ax.set_ylabel(ylabel, color='white', fontsize=10)
    ax.tick_params(colors='white')
    for s in ax.spines.values(): s.set_edgecolor('#444444')
    ax.grid(color='#2a2a2a', lw=0.5)
    ax.legend(fontsize=8.5, facecolor='#1a1a1a', edgecolor='#444444', labelcolor='white')
    return r, p


# ── Figure 50: pairwise scatters, MWF as the reference (gold standard) ────────
fig, axes = plt.subplots(1, 3, figsize=(18, 6), facecolor='#0d0d0d')
fig.suptitle('Three-way myelin agreement (per JHU ROI) — MWF as reference',
             color='white', fontsize=14, fontweight='bold')
x, xe = valid['MWF_mean'].values, valid['MWF_sem'].values
r1, _ = _panel(axes[0], x, valid['MVF_atlas_mean'].values, xe, valid['MVF_sem'].values,
               'T2-GRASE MWF (fraction)', 'MIMM atlas MVF (fraction)', C_MVF, identity=True)
r2, _ = _panel(axes[1], x, valid['chineg'].values, xe, valid['chineg_sem'].values,
               'T2-GRASE MWF (fraction)', r'chi-sep $|\chi^-|$ (ppm)', C_CHI)
r3, _ = _panel(axes[2], valid['MVF_atlas_mean'].values, valid['chineg'].values,
               valid['MVF_sem'].values, valid['chineg_sem'].values,
               'MIMM atlas MVF (fraction)', r'chi-sep $|\chi^-|$ (ppm)', C_CHI)
axes[0].set_title('MIMM vs MWF', color='white', fontsize=11)
axes[1].set_title('chi-sep vs MWF', color='white', fontsize=11)
axes[2].set_title('MIMM vs chi-sep', color='white', fontsize=11)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, '50_three_way_scatter.png'),
            dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
plt.close(); print(f'Saved: 50_three_way_scatter.png  (MIMM-MWF r={r1:.2f}, '
                   f'chi-MWF r={r2:.2f}, MIMM-chi r={r3:.2f})')


# ── Figure 51: per-ROI ranking, each marker min-max normalised to 0–1 ────────
# Different units (fraction vs ppm) can't share a raw axis, so each marker is
# rescaled to its own [0,1] range across ROIs. This isolates *ranking* agreement.
def _norm(s):
    s = s.astype(float)
    return (s - s.min()) / (s.max() - s.min() + 1e-12)


vr = valid.sort_values('MWF_mean', ascending=True).copy()
vr['MVF_n'] = _norm(vr['MVF_atlas_mean'])
vr['chi_n'] = _norm(vr['chineg'])
vr['MWF_n'] = _norm(vr['MWF_mean'])

names = vr['ROI_name'].tolist()
y = np.arange(len(vr)); h = 0.26
fig, ax = plt.subplots(figsize=(10, 14), facecolor='#0d0d0d')
ax.set_facecolor('#0d0d0d')
ax.barh(y + h, vr['MVF_n'], h, color=C_MVF, alpha=0.9, label='MIMM atlas MVF')
ax.barh(y,     vr['chi_n'], h, color=C_CHI, alpha=0.9, label='chi-sep |χ⁻|')
ax.barh(y - h, vr['MWF_n'], h, color=C_MWF, alpha=0.9, label='T2-GRASE MWF')
ax.set_yticks(y); ax.set_yticklabels(names, color='white', fontsize=7.5)
ax.set_xlabel('Relative myelin content (each marker min–max normalised 0–1)',
              color='white', fontsize=11)
ax.set_title('Tract ranking agreement across three myelin markers\n'
             '(ROIs sorted by MWF; aligned bars = consistent ranking)',
             color='white', fontsize=12, fontweight='bold', pad=12)
ax.set_xlim(0, 1.05)
ax.grid(axis='x', color='#333333', lw=0.5, ls='--')
ax.tick_params(colors='white')
for s in ax.spines.values(): s.set_edgecolor('#444444')
ax.legend(fontsize=10, facecolor='#1a1a1a', edgecolor='#555555',
          labelcolor='white', loc='lower right')
plt.tight_layout()
plt.savefig(os.path.join(out_dir, '51_three_way_ranking.png'),
            dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
plt.close(); print('Saved: 51_three_way_ranking.png')
