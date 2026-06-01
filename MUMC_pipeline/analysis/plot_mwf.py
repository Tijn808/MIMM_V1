"""
MWF validation figures — dormant until T2-GRASE data arrives.

Requires grase/MWF.nii.gz in the subject directory (written by ingest_grase.m).
If that file is absent this script exits cleanly with a message.

Per-subject outputs (figures 38–41):
  38_MWF_spatial.png        — MWF spatial map alongside MVF basic and atlas
  39_MVF_vs_MWF_ROI.png     — per-ROI scatter MVF vs MWF (identity line, both fractions)
  40_BA_MVF_vs_MWF.png      — Bland-Altman MVF vs MWF (voxel-level, WM mask)
  41_MVF_MWF_bar.png        — ROI bar chart: MVF basic, MVF atlas, MWF side by side

Cohort outputs produced by plot_cohort_roi.py when MWF is in roi_stats.csv:
  37_cohort_MVF_vs_MWF.png  — cohort scatter MVF vs MWF (between-subject SEM)
"""

import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import os

try:
    from paths import OUTPUT_DIR as subj_dir, FIG_DIR as out_dir, ANALYSIS_DIR
except ImportError:
    raise SystemExit('Copy MUMC_pipeline/analysis/paths_template.py to paths.py.')
os.makedirs(out_dir, exist_ok=True)

# Guard — exit cleanly if GRASE data is not yet available
mwf_path = os.path.join(subj_dir, 'grase', 'MWF.nii.gz')
if not os.path.exists(mwf_path):
    print('MWF.nii.gz not found — GRASE data not yet ingested.')
    print(f'  Expected: {mwf_path}')
    print('  Run MUMC_pipeline/grase/ingest_grase.m first.')
    raise SystemExit(0)

print(f'T2-GRASE MWF found: {mwf_path}')

def load(path):
    return np.array(nib.load(path).dataobj).astype(np.float32)

brain = load(os.path.join(subj_dir, 'qsm', 'brain_mask.nii.gz')) > 0
mag   = load(os.path.join(subj_dir, 'qsm', 'mag_e1.nii.gz'))
fa    = load(os.path.join(subj_dir, 'atlas', 'FA_atlas.nii.gz'))
wm    = brain & (fa > 0.20)

mvf_b = load(os.path.join(subj_dir, 'mimm', 'MVF_basic.nii.gz'))
mvf_a = load(os.path.join(subj_dir, 'mimm', 'MVF_Atlas.nii.gz'))
mwf   = load(mwf_path)

# Clip MWF to physiological range (0–0.5); registration artefacts can give
# values outside this window.
mwf = np.clip(mwf, 0, 0.5)
for v in [mvf_b, mvf_a, mwf]: v[~brain] = 0

def _spine(ax):
    ax.tick_params(colors='white')
    for s in ax.spines.values(): s.set_edgecolor('white')

# ── Figure 38: Spatial maps — MVF basic, MVF atlas, MWF ──────────────────────
def _slices(vol):
    coords = np.argwhere(vol > 0)
    cx, cy, cz = coords.mean(axis=0).astype(int)
    return np.rot90(vol[:,:,cz]), np.rot90(vol[:,cy,:]), np.rot90(vol[cx,:,:])

maps38 = [('MVF — Basic',  mvf_b, 'hot',    0, 0.50),
          ('MVF — Atlas',  mvf_a, 'hot',    0, 0.50),
          ('MWF (GRASE)',  mwf,   'hot',    0, 0.30)]

fig, axes = plt.subplots(3, 3, figsize=(14, 12), facecolor='black')
fig.suptitle('MVF Basic | MVF Atlas | T2-GRASE MWF',
             color='white', fontsize=15, fontweight='bold')
view_labels = ['Axial', 'Coronal', 'Sagittal']
for r, (label, vol, cmap, vmin, vmax) in enumerate(maps38):
    slices = _slices(np.where(wm, vol, 0) if r < 2 else np.where(brain, vol, 0))
    bg_slices = _slices(mag)
    for c in range(3):
        a = axes[r, c]
        a.set_facecolor('black')
        bg = bg_slices[c]
        a.imshow((bg - bg.min()) / (bg.max() - bg.min() + 1e-8),
                 cmap='gray', aspect='auto', interpolation='nearest')
        a.imshow(np.ma.masked_where(slices[c] == 0, slices[c]),
                 cmap=cmap, vmin=vmin, vmax=vmax,
                 aspect='auto', interpolation='nearest', alpha=0.85)
        if r == 0: a.set_title(view_labels[c], color='white', fontsize=11)
        a.axis('off')
    axes[r, 0].set_ylabel(label, color='white', fontsize=11,
                          rotation=90, labelpad=8)

plt.tight_layout()
plt.savefig(os.path.join(out_dir, '38_MWF_spatial.png'),
            dpi=150, bbox_inches='tight', facecolor='black')
plt.close(); print('Saved: 38_MWF_spatial.png')

# ── Figure 39: Per-ROI scatter MVF (basic) vs MWF ────────────────────────────
# Both are myelin-volume-like fractions — identity line is valid.
roi_csv = os.path.join(ANALYSIS_DIR, 'roi_stats.csv')
if os.path.exists(roi_csv) and 'MWF_mean' in pd.read_csv(roi_csv, nrows=0).columns:
    df = pd.read_csv(roi_csv)
    valid = df.dropna(subset=['MVF_basic_mean', 'MWF_mean'])
    if len(valid) > 3:
        r, p = pearsonr(valid['MWF_mean'], valid['MVF_basic_mean'])
        m, b = np.polyfit(valid['MWF_mean'], valid['MVF_basic_mean'], 1)
        lim = (0, max(valid['MWF_mean'].max(), valid['MVF_basic_mean'].max()) * 1.05)
        x_fit = np.linspace(lim[0], lim[1], 100)

        fig, ax = plt.subplots(figsize=(7, 7), facecolor='#0d0d0d')
        ax.set_facecolor('#111111')

        mvf_sem = df['MVF_basic_sd'] / np.sqrt(df['n_voxels'])
        mwf_sem = df['MWF_sd'] / np.sqrt(df['n_voxels'])
        ax.errorbar(valid['MWF_mean'], valid['MVF_basic_mean'],
                    xerr=mwf_sem[valid.index], yerr=mvf_sem[valid.index],
                    fmt='o', color='#e0a050', markersize=5, alpha=0.85,
                    ecolor='#888888', elinewidth=0.7, capsize=2)

        ax.plot(lim, lim, '-', color='#555555', lw=1.2, label='Identity (y = x)')
        ax.plot(x_fit, m * x_fit + b, '--', color='cyan', lw=1.5,
                label=f'r = {r:.2f}, p = {p:.3f}')

        for _, row in valid.iterrows():
            for name, short in [('Splenium of corpus callosum', 'Splenium CC'),
                                 ('Posterior limb of internal capsule R', 'Post. IC'),
                                 ('Cingulum (hippocampus) R', 'Cing. (hipp.)')]:
                if row['ROI_name'] == name:
                    ax.annotate(short, (row['MWF_mean'], row['MVF_basic_mean']),
                                textcoords='offset points', xytext=(6, 3),
                                fontsize=7.5, color='white', alpha=0.9)

        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel('MWF — T2-GRASE (fraction)', color='white', fontsize=11)
        ax.set_ylabel('MVF — MIMM Basic (fraction)', color='white', fontsize=11)
        ax.set_title('MIMM MVF vs T2-GRASE MWF per JHU ROI\n'
                     '(both myelin-volume fractions — identity line = perfect agreement)',
                     color='white', fontsize=12, fontweight='bold')
        ax.tick_params(colors='white')
        for s in ax.spines.values(): s.set_edgecolor('#444444')
        ax.grid(color='#2a2a2a', lw=0.5)
        ax.legend(fontsize=9, facecolor='#1a1a1a', edgecolor='#444444', labelcolor='white')
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, '39_MVF_vs_MWF_ROI.png'),
                    dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
        plt.close(); print('Saved: 39_MVF_vs_MWF_ROI.png')
else:
    print('Skipped: 39_MVF_vs_MWF_ROI.png (run extract_roi_stats.py first)')

# ── Figure 40: Bland-Altman MVF (basic) vs MWF (voxel, WM mask) ──────────────
# This is a valid Bland-Altman: both are myelin-volume fractions, same units.
v_mvf = mvf_b[wm].ravel().astype(float)
v_mwf = mwf[wm].ravel().astype(float)
finite = np.isfinite(v_mvf) & np.isfinite(v_mwf)
v_mvf, v_mwf = v_mvf[finite], v_mwf[finite]

mean_  = (v_mvf + v_mwf) / 2
diff_  = v_mvf - v_mwf
bias   = np.mean(diff_)
sd_    = np.std(diff_)
loa_u  = bias + 1.96 * sd_
loa_l  = bias - 1.96 * sd_
xlim   = (np.percentile(mean_, 1), np.percentile(mean_, 99))
ylim_  = (np.percentile(diff_, 1) - 0.02, np.percentile(diff_, 99) + 0.02)
keep   = ((mean_ >= xlim[0]) & (mean_ <= xlim[1]) &
          (diff_ >= ylim_[0]) & (diff_ <= ylim_[1]))

fig, ax = plt.subplots(figsize=(8, 6), facecolor='black')
ax.set_facecolor('#111111')
ax.hist2d(mean_[keep], diff_[keep], bins=150, cmap='hot', cmin=1)
ax.axhline(bias,  color='cyan',   lw=1.8, ls='--', label=f'Bias = {bias:.3f}')
ax.axhline(loa_u, color='yellow', lw=1.4, ls=':',  label=f'+1.96 SD = {loa_u:.3f}')
ax.axhline(loa_l, color='yellow', lw=1.4, ls=':',  label=f'−1.96 SD = {loa_l:.3f}')
ax.axhline(0, color='white', lw=0.8, alpha=0.3)
ax.set_xlim(xlim); ax.set_ylim(ylim_)
ax.set_xlabel('Mean (MVF, MWF)  [fraction]', color='white', fontsize=11)
ax.set_ylabel('MVF − MWF  [fraction]', color='white', fontsize=11)
ax.set_title('Bland-Altman: MIMM MVF vs T2-GRASE MWF\n(voxel-level, WM mask FA > 0.2)',
             color='white', fontsize=13, fontweight='bold')
ax.tick_params(colors='white')
for s in ax.spines.values(): s.set_edgecolor('white')
ax.legend(fontsize=10, facecolor='#222222', edgecolor='white', labelcolor='white')
plt.tight_layout()
plt.savefig(os.path.join(out_dir, '40_BA_MVF_vs_MWF.png'),
            dpi=150, bbox_inches='tight', facecolor='black')
plt.close(); print('Saved: 40_BA_MVF_vs_MWF.png')
print(f'  BA stats: bias={bias:.4f}  ±1.96SD=[{loa_l:.4f}, {loa_u:.4f}]')

# ── Figure 41: ROI bar chart — MVF basic, MVF atlas, MWF ─────────────────────
roi_csv = os.path.join(ANALYSIS_DIR, 'roi_stats.csv')
if os.path.exists(roi_csv) and 'MWF_mean' in pd.read_csv(roi_csv, nrows=0).columns:
    df = pd.read_csv(roi_csv).sort_values('MVF_basic_mean', ascending=True)
    names = df['ROI_name'].tolist()
    y = np.arange(len(df)); h = 0.28

    df['MVF_basic_sem'] = df['MVF_basic_sd'] / np.sqrt(df['n_voxels'])
    df['MVF_atlas_sem'] = df['MVF_atlas_sd'] / np.sqrt(df['n_voxels'])
    df['MWF_sem']       = df['MWF_sd']       / np.sqrt(df['n_voxels'])

    fig, ax = plt.subplots(figsize=(10, 14), facecolor='#0d0d0d')
    ax.set_facecolor('#0d0d0d')
    ax.barh(y + h, df['MVF_basic_mean'], h, xerr=df['MVF_basic_sem'],
            color='#e05c5c', alpha=0.85, error_kw=dict(ecolor='white', lw=0.8, capsize=2),
            label='MVF — Basic')
    ax.barh(y,     df['MVF_atlas_mean'], h, xerr=df['MVF_atlas_sem'],
            color='#5c9ee0', alpha=0.85, error_kw=dict(ecolor='white', lw=0.8, capsize=2),
            label='MVF — Atlas')
    ax.barh(y - h, df['MWF_mean'], h, xerr=df['MWF_sem'],
            color='#5cb85c', alpha=0.85, error_kw=dict(ecolor='white', lw=0.8, capsize=2),
            label='MWF — T2-GRASE')
    ax.set_yticks(y); ax.set_yticklabels(names, color='white', fontsize=7.5)
    ax.set_xlabel('Myelin fraction', color='white', fontsize=11)
    ax.set_title('MVF Basic | MVF Atlas | T2-GRASE MWF per JHU ROI\n(error bars = ±1 SEM)',
                 color='white', fontsize=12, fontweight='bold', pad=12)
    ax.set_xlim(0, 0.55); ax.axvline(0, color='#444444', lw=0.8)
    ax.grid(axis='x', color='#333333', lw=0.5, ls='--')
    ax.tick_params(colors='white')
    for s in ax.spines.values(): s.set_edgecolor('#444444')
    ax.legend(fontsize=10, facecolor='#1a1a1a', edgecolor='#555555',
              labelcolor='white', loc='lower right')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '41_MVF_MWF_bar.png'),
                dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
    plt.close(); print('Saved: 41_MVF_MWF_bar.png')
else:
    print('Skipped: 41_MVF_MWF_bar.png (run extract_roi_stats.py first)')
