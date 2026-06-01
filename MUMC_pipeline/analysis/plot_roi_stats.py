"""
Statistical figures for JHU ROI analysis of MIMM output maps.
"""

import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.stats import pearsonr
import os

try:
    from paths import OUTPUT_DIR as subj_dir, FIG_DIR as out_dir, ANALYSIS_DIR
except ImportError:
    raise SystemExit('Copy MUMC_pipeline/analysis/paths_template.py to paths.py and fill in your paths.')
os.makedirs(out_dir, exist_ok=True)

df = pd.read_csv(os.path.join(ANALYSIS_DIR, 'roi_stats.csv'))
df = df.sort_values('MVF_basic_mean', ascending=True).reset_index(drop=True)
df['chi_neg_abs'] = df['chi_neg_chisep_mean'].abs()   # chi_neg is diamagnetic (negative ppm)
names = df['ROI_name'].tolist()

# Error bars compare ROI MEANS, so use SEM (uncertainty on the mean), not SD
# (within-ROI voxel spread). SEM = SD / sqrt(n_voxels). With thousands of voxels
# per ROI this is far smaller than SD and does not swamp the plot.
# (Voxels are spatially correlated so this mildly overstates precision, but it is
# the correct quantity for comparing means.)
for _c in ['MVF_basic', 'MVF_atlas', 'chi_neg_chisep', 'chi_pos_chisep',
           'chi_iron_basic', 'chi_myelin_basic']:
    df[f'{_c}_sem'] = df[f'{_c}_sd'] / np.sqrt(df['n_voxels'])

# Per-ROI FA mean from atlas (needed for fig 29; not stored in roi_stats.csv)
_fa_path  = os.path.join(subj_dir, 'atlas', 'FA_atlas.nii.gz')
_lbl_path = os.path.join(subj_dir, 'atlas', 'JHU_labels_subj.nii.gz')
if os.path.exists(_fa_path) and os.path.exists(_lbl_path):
    _fa_vol  = np.array(nib.load(_fa_path).dataobj).astype(np.float32)
    _lbl_vol = np.array(nib.load(_lbl_path).dataobj).astype(int)
    df['fa_mean'] = [
        float(_fa_vol[_lbl_vol == idx].mean()) if (_lbl_vol == idx).any() else np.nan
        for idx in df['ROI_index']
    ]
else:
    df['fa_mean'] = np.nan

plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 9})

# ── Figure 1: MVF Basic vs Atlas — horizontal bar chart ──────────────────────
fig, ax = plt.subplots(figsize=(10, 14), facecolor='#0d0d0d')
ax.set_facecolor('#0d0d0d')

y  = np.arange(len(df))
h  = 0.35

bars_b = ax.barh(y + h/2, df['MVF_basic_mean'], h,
                 xerr=df['MVF_basic_sem'], color='#e05c5c', alpha=0.85,
                 error_kw=dict(ecolor='white', lw=0.8, capsize=2),
                 label='Basic (no prior)')
bars_a = ax.barh(y - h/2, df['MVF_atlas_mean'], h,
                 xerr=df['MVF_atlas_sem'], color='#5c9ee0', alpha=0.85,
                 error_kw=dict(ecolor='white', lw=0.8, capsize=2),
                 label='Atlas (HCP1065 prior)')

ax.set_yticks(y)
ax.set_yticklabels(names, color='white', fontsize=7.5)
ax.set_xlabel('MVF (fraction)', color='white', fontsize=11)
ax.set_title('Myelin Volume Fraction per JHU White Matter ROI\nBasic vs Atlas MIMM',
             color='white', fontsize=13, fontweight='bold', pad=12)
ax.tick_params(colors='white')
ax.xaxis.set_tick_params(color='white')
for spine in ax.spines.values():
    spine.set_edgecolor('#444444')
ax.set_xlim(0, 0.55)
ax.axvline(0, color='#444444', linewidth=0.8)
ax.grid(axis='x', color='#333333', linewidth=0.5, linestyle='--')
legend = ax.legend(fontsize=10, facecolor='#1a1a1a', edgecolor='#555555',
                   labelcolor='white', loc='lower right')

plt.tight_layout()
plt.savefig(os.path.join(out_dir, '18_ROI_MVF_basic_vs_atlas.png'),
            dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
plt.close()
print('Saved: 18_ROI_MVF_basic_vs_atlas.png')

# ── Figure 2: Heatmap — all ROIs × key parameters ───────────────────────────
params = {
    'MVF\nbasic':    ('MVF_basic_mean',     0, 0.45),
    'MVF\natlas':    ('MVF_atlas_mean',      0, 0.45),
    'FVF\nbasic':    ('FVF_basic_mean',      0, 0.80),
    'g-ratio\nbasic':('g_ratio_basic_mean', 0.5, 1.0),
    'R2*\nbasic':    ('R2s_basic_mean',      0, 40),
    '|χ⁻|\nchi-sep': ('chi_neg_abs',          0, 0.15),
    'χ⁺\nchi-sep':   ('chi_pos_chisep_mean',  0, 0.15),
}

# Normalise each column 0→1 for coloring
mat = np.zeros((len(df), len(params)))
for j, (col_label, (col, vmin, vmax)) in enumerate(params.items()):
    vals = df[col].values
    mat[:, j] = np.clip((vals - vmin) / (vmax - vmin), 0, 1)

fig, ax = plt.subplots(figsize=(10, 14), facecolor='#0d0d0d')
ax.set_facecolor('#0d0d0d')

im = ax.imshow(mat, aspect='auto', cmap='hot', vmin=0, vmax=1,
               interpolation='nearest')

ax.set_yticks(range(len(df)))
ax.set_yticklabels(names, color='white', fontsize=7.5)
ax.set_xticks(range(len(params)))
ax.set_xticklabels(list(params.keys()), color='white', fontsize=9)
ax.tick_params(colors='white', length=0)
for spine in ax.spines.values():
    spine.set_edgecolor('#333333')

# Annotate cells with actual values
for i in range(len(df)):
    for j, (col_label, (col, vmin, vmax)) in enumerate(params.items()):
        val = df[col].iloc[i]
        brightness = mat[i, j]
        txt_color = 'black' if brightness > 0.6 else 'white'
        ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                fontsize=5.5, color=txt_color)

ax.set_title('MIMM & χ-separation — ROI Parameter Heatmap\n(normalised per column)',
             color='white', fontsize=13, fontweight='bold', pad=12)

cbar = fig.colorbar(im, ax=ax, shrink=0.4, pad=0.02)
cbar.set_label('Normalised value', color='white', fontsize=9)
cbar.ax.yaxis.set_tick_params(color='white')
plt.setp(cbar.ax.yaxis.get_ticklabels(), color='white')
cbar.outline.set_edgecolor('white')

plt.tight_layout()
plt.savefig(os.path.join(out_dir, '19_ROI_heatmap.png'),
            dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
plt.close()
print('Saved: 19_ROI_heatmap.png')

# ── Figure 3: Scatter — MVF basic vs atlas per ROI ──────────────────────────
fig, ax = plt.subplots(figsize=(7, 7), facecolor='#0d0d0d')
ax.set_facecolor('#0d0d0d')

ax.errorbar(df['MVF_basic_mean'], df['MVF_atlas_mean'],
            xerr=df['MVF_basic_sem'], yerr=df['MVF_atlas_sem'],
            fmt='o', color='#e0a050', markersize=5, alpha=0.8,
            ecolor='#888888', elinewidth=0.7, capsize=2)

# Identity line
lim = (0, 0.45)
ax.plot(lim, lim, '--', color='#555555', linewidth=1, label='y = x')

# Label key tracts with custom offsets + leader lines to avoid collisions in
# the dense 0.18-0.25 cluster. (dx, dy in points, ha alignment.)
label_offsets = {
    'Posterior limb of internal capsule R': (8,   0, 'left'),
    'Splenium of corpus callosum':          (8,  -2, 'left'),
    'Corticospinal tract R':                (12, 14, 'left'),
    'Genu of corpus callosum':              (16, -16, 'left'),
    'Body of corpus callosum':              (20, -30, 'left'),
    'Cingulum (cingulate gyrus) R':         (-12, 16, 'right'),
}
for _, row in df.iterrows():
    if row['ROI_name'] in label_offsets:
        dx, dy, ha = label_offsets[row['ROI_name']]
        short = row['ROI_name'].replace(' R', '').replace(' L', '')
        ax.annotate(short, (row['MVF_basic_mean'], row['MVF_atlas_mean']),
                    textcoords='offset points', xytext=(dx, dy), ha=ha,
                    fontsize=7, color='white', alpha=0.9,
                    arrowprops=dict(arrowstyle='-', color='#888888', lw=0.5))

ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel('MVF — Basic (no prior)', color='white', fontsize=11)
ax.set_ylabel('MVF — Atlas (HCP1065 prior)', color='white', fontsize=11)
ax.set_title('MVF: Basic vs Atlas per JHU ROI\n(error bars = ±1 SEM)',
             color='white', fontsize=13, fontweight='bold')
ax.tick_params(colors='white')
for spine in ax.spines.values():
    spine.set_edgecolor('#444444')
ax.grid(color='#2a2a2a', linewidth=0.5)
ax.legend(fontsize=9, facecolor='#1a1a1a', edgecolor='#555555', labelcolor='white')

plt.tight_layout()
plt.savefig(os.path.join(out_dir, '20_ROI_MVF_scatter.png'),
            dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
plt.close()
print('Saved: 20_ROI_MVF_scatter.png')

# ── Figure 21: Tractography atlas — MVF Basic vs Atlas bar chart ─────────────
tracts_csv = os.path.join(subj_dir, 'analysis', 'roi_stats_tracts.csv')
if os.path.exists(tracts_csv):
    dft = pd.read_csv(tracts_csv).sort_values('MVF_basic_fa_weighted_mean', ascending=True)
    names_t = dft['ROI_name'].tolist()
    y_t = np.arange(len(dft))

    # SEM for error bars (see note at df load) — uncertainty on the ROI mean
    dft['MVF_basic_sem'] = dft['MVF_basic_sd'] / np.sqrt(dft['n_voxels'])
    dft['MVF_atlas_sem'] = dft['MVF_atlas_sd'] / np.sqrt(dft['n_voxels'])

    fig, ax = plt.subplots(figsize=(10, 7), facecolor='#0d0d0d')
    ax.set_facecolor('#0d0d0d')
    ax.barh(y_t + 0.2, dft['MVF_basic_mean'], 0.4,
            xerr=dft['MVF_basic_sem'], color='#e05c5c', alpha=0.85,
            error_kw=dict(ecolor='white', lw=0.8, capsize=2), label='Basic')
    ax.barh(y_t - 0.2, dft['MVF_atlas_mean'], 0.4,
            xerr=dft['MVF_atlas_sem'], color='#5c9ee0', alpha=0.85,
            error_kw=dict(ecolor='white', lw=0.8, capsize=2), label='Atlas')
    ax.set_yticks(y_t)
    ax.set_yticklabels(names_t, color='white', fontsize=8)
    ax.set_xlabel('MVF (fraction)', color='white', fontsize=11)
    ax.set_title('MVF per JHU Tractography ROI (thr25)\nBasic vs Atlas MIMM',
                 color='white', fontsize=13, fontweight='bold', pad=12)
    ax.tick_params(colors='white')
    for spine in ax.spines.values(): spine.set_edgecolor('#444444')
    ax.set_xlim(0, 0.55)
    ax.grid(axis='x', color='#333333', linewidth=0.5, linestyle='--')
    ax.legend(fontsize=10, facecolor='#1a1a1a', edgecolor='#555555',
              labelcolor='white', loc='lower right')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '21_ROI_tracts_MVF.png'),
                dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
    plt.close()
    print('Saved: 21_ROI_tracts_MVF.png')

    # ── Figure 22: Atlas comparison — DTI-81 vs Tractography CV + r ──────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), facecolor='#0d0d0d')
    fig.suptitle('JHU DTI-81 (48 ROIs) vs Tractography thr25 (20 ROIs)',
                 color='white', fontsize=13, fontweight='bold')

    for ax in axes: ax.set_facecolor('#111111')

    # Panel 1: within-ROI CV comparison
    cv_dti    = df['MVF_basic_cv'].dropna()
    cv_tracts = dft['MVF_basic_cv'].dropna()
    axes[0].hist(cv_dti,    bins=15, color='#e05c5c', alpha=0.75, label=f'DTI-81 (n=48, μ={cv_dti.mean():.2f})')
    axes[0].hist(cv_tracts, bins=10, color='#5c9ee0', alpha=0.75, label=f'Tracts (n=20, μ={cv_tracts.mean():.2f})')
    axes[0].set_xlabel('Within-ROI CV (MVF_basic)', color='white', fontsize=10)
    axes[0].set_ylabel('Number of ROIs', color='white', fontsize=10)
    axes[0].set_title('ROI Homogeneity\n(lower CV = more homogeneous)', color='white', fontsize=11)
    axes[0].tick_params(colors='white')
    for s in axes[0].spines.values(): s.set_edgecolor('#444444')
    axes[0].legend(fontsize=9, facecolor='#1a1a1a', edgecolor='#555555', labelcolor='white')

    # Panel 2: Pearson r MVF vs chi_neg
    r_dti    = df['r_MVF_vs_chineg'].dropna()
    r_tracts = dft['r_MVF_vs_chineg'].dropna()
    axes[1].hist(r_dti,    bins=15, color='#e05c5c', alpha=0.75, label=f'DTI-81 (μ={r_dti.mean():.2f})')
    axes[1].hist(r_tracts, bins=10, color='#5c9ee0', alpha=0.75, label=f'Tracts (μ={r_tracts.mean():.2f})')
    axes[1].set_xlabel('Pearson r (MVF_basic vs χ-sep χ_neg)', color='white', fontsize=10)
    axes[1].set_ylabel('Number of ROIs', color='white', fontsize=10)
    axes[1].set_title('MIMM–χ-sep Spatial Agreement\n(higher r = better agreement)', color='white', fontsize=11)
    axes[1].tick_params(colors='white')
    for s in axes[1].spines.values(): s.set_edgecolor('#444444')
    axes[1].legend(fontsize=9, facecolor='#1a1a1a', edgecolor='#555555', labelcolor='white')

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '22_atlas_comparison.png'),
                dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
    plt.close()
    print('Saved: 22_atlas_comparison.png')

# ── Tract family colour map (JHU label index → family) ───────────────────────
TRACT_FAMILIES = {
    'Corpus callosum':         ([3, 4, 5],                              '#f4c430'),
    'Internal capsule':        ([17, 18, 19, 20, 21, 22, 29, 30],      '#5c9ee0'),
    'Corticospinal/brainstem': ([1, 2, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], '#5cb85c'),
    'Corona radiata':          ([23, 24, 25, 26, 27, 28],              '#e8a838'),
    'Association fibres':      ([31, 32, 33, 34, 41, 42, 43, 44, 45, 46], '#9b59b6'),
    'Limbic/cingulum':         ([6, 35, 36, 37, 38, 39, 40, 47, 48, 49, 50], '#e05c9b'),
}

def family_color(roi_idx):
    for color_info in TRACT_FAMILIES.values():
        if roi_idx in color_info[0]:
            return color_info[1]
    return '#aaaaaa'

df_plot = df.copy()
df_plot['family_color'] = df_plot['ROI_index'].apply(family_color)

# ── Figure 26: MVF (MIMM) vs |χ⁻| (chi-sep) per ROI ─────────────────────────
# Both are myelin markers validated against histopathology (Şişman et al. 2024).
# Tests whether the two methods rank WM tracts by myelination consistently.
valid26 = df_plot.dropna(subset=['MVF_basic_mean', 'chi_neg_abs'])
r26, p26 = pearsonr(valid26['chi_neg_abs'], valid26['MVF_basic_mean'])
m26, b26 = np.polyfit(valid26['chi_neg_abs'], valid26['MVF_basic_mean'], 1)
x_line26 = np.linspace(valid26['chi_neg_abs'].min(), valid26['chi_neg_abs'].max(), 100)

fig, ax = plt.subplots(figsize=(8, 7), facecolor='#0d0d0d')
ax.set_facecolor('#111111')

for family, (indices, color) in TRACT_FAMILIES.items():
    sub = valid26[valid26['ROI_index'].isin(indices)]
    if sub.empty:
        continue
    ax.errorbar(sub['chi_neg_abs'], sub['MVF_basic_mean'],
                xerr=sub['chi_neg_chisep_sem'], yerr=sub['MVF_basic_sem'],
                fmt='o', color=color, markersize=6, alpha=0.85,
                ecolor='#555555', elinewidth=0.6, capsize=2, label=family, zorder=3)

ax.plot(x_line26, m26 * x_line26 + b26, '--', color='cyan', linewidth=1.5,
        label=f'r = {r26:.2f}, p = {p26:.3f}', zorder=4)

# Label key ROIs at the extremes
highlight26 = {
    'Splenium of corpus callosum': 'Splenium CC',
    'Genu of corpus callosum':     'Genu CC',
    'Body of corpus callosum':     'Body CC',
    'Posterior limb of internal capsule R': 'Post. IC',
    'Corticospinal tract R':       'CST',
    'Cingulum (hippocampus) R':    'Cing. (hipp.)',
    'Fornix (column and body)':    'Fornix',
}
for _, row in valid26.iterrows():
    if row['ROI_name'] in highlight26:
        ax.annotate(highlight26[row['ROI_name']],
                    (row['chi_neg_abs'], row['MVF_basic_mean']),
                    textcoords='offset points', xytext=(6, 3),
                    fontsize=7.5, color='white', alpha=0.9)

ax.set_xlabel('|χ⁻| chi-sep  (ppm)', color='white', fontsize=11)
ax.set_ylabel('MVF MIMM  (fraction)', color='white', fontsize=11)
ax.set_title('Myelin Marker Consistency: MIMM MVF vs χ-sep |χ⁻|\nper JHU DTI-81 ROI  (χ⁻ histopathologically validated, Şişman 2024)',
             color='white', fontsize=12, fontweight='bold')
ax.tick_params(colors='white')
for spine in ax.spines.values(): spine.set_edgecolor('#444444')
ax.grid(color='#2a2a2a', linewidth=0.5)
ax.legend(fontsize=8, facecolor='#1a1a1a', edgecolor='#444444',
          labelcolor='white', loc='upper left')

plt.tight_layout()
plt.savefig(os.path.join(out_dir, '26_MVF_vs_chineg.png'),
            dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
plt.close()
print('Saved: 26_MVF_vs_chineg.png')

# ── Figure 27: χ⁺ (chi-sep) vs chi_iron_est (MIMM) per ROI ──────────────────
# Both are paramagnetic susceptibility estimates from the same mGRE data
# via different pipelines. Tests iron estimate consistency.
valid27 = df_plot.dropna(subset=['chi_iron_basic_mean', 'chi_pos_chisep_mean'])
r27, p27 = pearsonr(valid27['chi_pos_chisep_mean'], valid27['chi_iron_basic_mean'])
m27, b27 = np.polyfit(valid27['chi_pos_chisep_mean'], valid27['chi_iron_basic_mean'], 1)
x_line27 = np.linspace(valid27['chi_pos_chisep_mean'].min(), valid27['chi_pos_chisep_mean'].max(), 100)

fig, ax = plt.subplots(figsize=(8, 7), facecolor='#0d0d0d')
ax.set_facecolor('#111111')

for family, (indices, color) in TRACT_FAMILIES.items():
    sub = valid27[valid27['ROI_index'].isin(indices)]
    if sub.empty:
        continue
    ax.errorbar(sub['chi_pos_chisep_mean'], sub['chi_iron_basic_mean'],
                xerr=sub['chi_pos_chisep_sem'], yerr=sub['chi_iron_basic_sem'],
                fmt='o', color=color, markersize=6, alpha=0.85,
                ecolor='#555555', elinewidth=0.6, capsize=2, label=family, zorder=3)

ax.plot(x_line27, m27 * x_line27 + b27, '--', color='cyan', linewidth=1.5,
        label=f'r = {r27:.2f}, p = {p27:.3f}', zorder=4)

highlight27 = {
    'Middle cerebellar peduncle':              'MCP',
    'Posterior limb of internal capsule R':    'Post. IC',
    'Corticospinal tract R':                   'CST',
    'Retrolenticular internal capsule R':      'Retrolent. IC',
    'Splenium of corpus callosum':             'Splenium CC',
    'Cingulum (hippocampus) R':                'Cing. (hipp.)',
}
for _, row in valid27.iterrows():
    if row['ROI_name'] in highlight27:
        ax.annotate(highlight27[row['ROI_name']],
                    (row['chi_pos_chisep_mean'], row['chi_iron_basic_mean']),
                    textcoords='offset points', xytext=(6, 3),
                    fontsize=7.5, color='white', alpha=0.9)

ax.set_xlabel('χ⁺ chi-sep  (ppm)', color='white', fontsize=11)
ax.set_ylabel('χ iron MIMM  (ppm)', color='white', fontsize=11)
ax.set_title('Iron Estimate Consistency: MIMM χ_iron vs χ-sep χ⁺\nper JHU DTI-81 ROI  (χ⁺ outperforms QSM for iron, Şişman 2024)',
             color='white', fontsize=12, fontweight='bold')
ax.tick_params(colors='white')
for spine in ax.spines.values(): spine.set_edgecolor('#444444')
ax.grid(color='#2a2a2a', linewidth=0.5)
ax.legend(fontsize=8, facecolor='#1a1a1a', edgecolor='#444444',
          labelcolor='white', loc='upper left')

plt.tight_layout()
plt.savefig(os.path.join(out_dir, '27_iron_chi_comparison.png'),
            dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
plt.close()
print('Saved: 27_iron_chi_comparison.png')

# ── Figure 28: |χ myelin| MIMM vs |χ⁻| chi-sep per ROI ──────────────────────
# Same physical quantity (diamagnetic susceptibility, ppm) from two different
# pipelines. Identity line shows perfect agreement.
# chi_myelin_basic_mean is already positive (np.abs applied in extract_roi_stats.py).
valid28 = df_plot.dropna(subset=['chi_myelin_basic_mean', 'chi_neg_abs'])
r28, p28 = pearsonr(valid28['chi_myelin_basic_mean'], valid28['chi_neg_abs'])
m28, b28 = np.polyfit(valid28['chi_myelin_basic_mean'], valid28['chi_neg_abs'], 1)
x_line28 = np.linspace(valid28['chi_myelin_basic_mean'].min(),
                       valid28['chi_myelin_basic_mean'].max(), 100)

# Axis range for identity line: use shared max
ax_max28 = max(valid28['chi_myelin_basic_mean'].max(), valid28['chi_neg_abs'].max()) * 1.05

fig, ax = plt.subplots(figsize=(8, 7), facecolor='#0d0d0d')
ax.set_facecolor('#111111')

for family, (indices, color) in TRACT_FAMILIES.items():
    sub = valid28[valid28['ROI_index'].isin(indices)]
    if sub.empty:
        continue
    ax.errorbar(sub['chi_myelin_basic_mean'], sub['chi_neg_abs'],
                xerr=sub['chi_myelin_basic_sem'], yerr=sub['chi_neg_chisep_sem'],
                fmt='o', color=color, markersize=6, alpha=0.85,
                ecolor='#555555', elinewidth=0.6, capsize=2, label=family, zorder=3)

ax.plot([0, ax_max28], [0, ax_max28], '-', color='#555555',
        linewidth=1.2, label='Identity (y = x)', zorder=2)
ax.plot(x_line28, m28 * x_line28 + b28, '--', color='cyan', linewidth=1.5,
        label=f'r = {r28:.2f}, p = {p28:.3f}', zorder=4)

highlight28 = {
    'Splenium of corpus callosum':           'Splenium CC',
    'Genu of corpus callosum':               'Genu CC',
    'Body of corpus callosum':               'Body CC',
    'Posterior limb of internal capsule R':  'Post. IC',
    'Corticospinal tract R':                 'CST',
    'Cingulum (hippocampus) R':              'Cing. (hipp.)',
    'Fornix (column and body)':              'Fornix',
}
for _, row in valid28.iterrows():
    if row['ROI_name'] in highlight28:
        ax.annotate(highlight28[row['ROI_name']],
                    (row['chi_myelin_basic_mean'], row['chi_neg_abs']),
                    textcoords='offset points', xytext=(6, 3),
                    fontsize=7.5, color='white', alpha=0.9)

ax.set_xlim(0, ax_max28); ax.set_ylim(0, ax_max28)
ax.set_xlabel('|χ myelin| MIMM  (ppm)', color='white', fontsize=11)
ax.set_ylabel('|χ⁻| chi-sep  (ppm)', color='white', fontsize=11)
ax.set_title('Diamagnetic Susceptibility Agreement: MIMM vs χ-sep\nper JHU DTI-81 ROI  (both in ppm, identity line = perfect agreement)',
             color='white', fontsize=12, fontweight='bold')
ax.tick_params(colors='white')
for spine in ax.spines.values(): spine.set_edgecolor('#444444')
ax.grid(color='#2a2a2a', linewidth=0.5)
ax.legend(fontsize=8, facecolor='#1a1a1a', edgecolor='#444444',
          labelcolor='white', loc='upper left')

plt.tight_layout()
plt.savefig(os.path.join(out_dir, '28_chi_myelin_vs_chineg_ROI.png'),
            dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
plt.close()
print('Saved: 28_chi_myelin_vs_chineg_ROI.png')

# ── Figure 29: MVF (MIMM basic) vs FA atlas per ROI ──────────────────────────
# FA is from DTI — fully independent of mGRE. Positive correlation expected
# since both track WM myelination/integrity. Tests MIMM against an independent
# modality without any shared signal.
if df['fa_mean'].notna().sum() > 3:
    valid29 = df_plot.copy()
    valid29['fa_mean'] = df['fa_mean']
    valid29 = valid29.dropna(subset=['MVF_basic_mean', 'fa_mean'])

    r29, p29 = pearsonr(valid29['fa_mean'], valid29['MVF_basic_mean'])
    m29, b29 = np.polyfit(valid29['fa_mean'], valid29['MVF_basic_mean'], 1)
    x_line29 = np.linspace(valid29['fa_mean'].min(), valid29['fa_mean'].max(), 100)

    fig, ax = plt.subplots(figsize=(8, 7), facecolor='#0d0d0d')
    ax.set_facecolor('#111111')

    for family, (indices, color) in TRACT_FAMILIES.items():
        sub = valid29[valid29['ROI_index'].isin(indices)]
        if sub.empty:
            continue
        ax.errorbar(sub['fa_mean'], sub['MVF_basic_mean'],
                    yerr=sub['MVF_basic_sem'],
                    fmt='o', color=color, markersize=6, alpha=0.85,
                    ecolor='#555555', elinewidth=0.6, capsize=2, label=family, zorder=3)

    ax.plot(x_line29, m29 * x_line29 + b29, '--', color='cyan', linewidth=1.5,
            label=f'r = {r29:.2f}, p = {p29:.3f}', zorder=4)

    highlight29 = {
        'Splenium of corpus callosum':           'Splenium CC',
        'Genu of corpus callosum':               'Genu CC',
        'Body of corpus callosum':               'Body CC',
        'Posterior limb of internal capsule R':  'Post. IC',
        'Corticospinal tract R':                 'CST',
        'Cingulum (hippocampus) R':              'Cing. (hipp.)',
        'Middle cerebellar peduncle':            'MCP',
    }
    for _, row in valid29.iterrows():
        if row['ROI_name'] in highlight29:
            ax.annotate(highlight29[row['ROI_name']],
                        (row['fa_mean'], row['MVF_basic_mean']),
                        textcoords='offset points', xytext=(6, 3),
                        fontsize=7.5, color='white', alpha=0.9)

    ax.set_xlabel('Mean FA (HCP1065 atlas)', color='white', fontsize=11)
    ax.set_ylabel('MVF MIMM basic  (fraction)', color='white', fontsize=11)
    ax.set_title('Independent Validation: MIMM MVF vs DTI FA\nper JHU DTI-81 ROI  (FA from independent DTI acquisition)',
                 color='white', fontsize=12, fontweight='bold')
    ax.tick_params(colors='white')
    for spine in ax.spines.values(): spine.set_edgecolor('#444444')
    ax.grid(color='#2a2a2a', linewidth=0.5)
    ax.legend(fontsize=8, facecolor='#1a1a1a', edgecolor='#444444',
              labelcolor='white', loc='upper left')

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, '29_MVF_vs_FA.png'),
                dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
    plt.close()
    print('Saved: 29_MVF_vs_FA.png')
else:
    print('Skipped: 29_MVF_vs_FA.png (FA atlas NIfTI not found)')
