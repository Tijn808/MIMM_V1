"""
Statistical figures for JHU ROI analysis of MIMM output maps.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os

subj_dir = '/home/tijn-saes/Documents/Internship/ME_GRE'
out_dir  = os.path.join(subj_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

df = pd.read_csv(os.path.join(subj_dir, 'analysis', 'roi_stats.csv'))
df = df.sort_values('MVF_basic_mean', ascending=True).reset_index(drop=True)
names = df['ROI_name'].tolist()

plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 9})

# ── Figure 1: MVF Basic vs Atlas — horizontal bar chart ──────────────────────
fig, ax = plt.subplots(figsize=(10, 14), facecolor='#0d0d0d')
ax.set_facecolor('#0d0d0d')

y  = np.arange(len(df))
h  = 0.35

bars_b = ax.barh(y + h/2, df['MVF_basic_mean'], h,
                 xerr=df['MVF_basic_sd'], color='#e05c5c', alpha=0.85,
                 error_kw=dict(ecolor='white', lw=0.8, capsize=2),
                 label='Basic (no prior)')
bars_a = ax.barh(y - h/2, df['MVF_atlas_mean'], h,
                 xerr=df['MVF_atlas_sd'], color='#5c9ee0', alpha=0.85,
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
    'χ neg\nchi-sep':('chi_neg_chisep_mean', 0, 0.15),
    'χ pos\nchi-sep':('chi_pos_chisep_mean', 0, 0.15),
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
            xerr=df['MVF_basic_sd'], yerr=df['MVF_atlas_sd'],
            fmt='o', color='#e0a050', markersize=5, alpha=0.8,
            ecolor='#666666', elinewidth=0.7, capsize=2)

# Identity line
lim = (0, 0.45)
ax.plot(lim, lim, '--', color='#555555', linewidth=1, label='y = x')

# Label a few key tracts
highlight = ['Splenium of corpus callosum', 'Genu of corpus callosum',
             'Body of corpus callosum', 'Posterior limb of internal capsule R',
             'Corticospinal tract R', 'Cingulum (cingulate gyrus) R']
for _, row in df.iterrows():
    if row['ROI_name'] in highlight:
        ax.annotate(row['ROI_name'].replace(' R', '').replace(' L', ''),
                    (row['MVF_basic_mean'], row['MVF_atlas_mean']),
                    textcoords='offset points', xytext=(6, 2),
                    fontsize=7, color='white', alpha=0.9)

ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel('MVF — Basic (no prior)', color='white', fontsize=11)
ax.set_ylabel('MVF — Atlas (HCP1065 prior)', color='white', fontsize=11)
ax.set_title('MVF: Basic vs Atlas per JHU ROI\n(error bars = ±1 SD within ROI)',
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

    fig, ax = plt.subplots(figsize=(10, 7), facecolor='#0d0d0d')
    ax.set_facecolor('#0d0d0d')
    ax.barh(y_t + 0.2, dft['MVF_basic_mean'], 0.4,
            xerr=dft['MVF_basic_sd'], color='#e05c5c', alpha=0.85,
            error_kw=dict(ecolor='white', lw=0.8, capsize=2), label='Basic')
    ax.barh(y_t - 0.2, dft['MVF_atlas_mean'], 0.4,
            xerr=dft['MVF_atlas_sd'], color='#5c9ee0', alpha=0.85,
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
