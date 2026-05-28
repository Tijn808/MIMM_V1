"""
Overestimation analysis: MIMM Basic vs Atlas MVF per JHU ROI.

Tests the core physics hypothesis:
  Basic MIMM overestimates MVF most in tracts running perpendicular to B0,
  because without an orientation prior it confuses angle-dependent signal
  decay with myelin content.

Produces:
  23_overestimation_ranked.png   -- all ROIs ranked by overestimation
  24_overestimation_vs_theta.png -- scatter: mean fibre angle vs overestimation
  25_overestimation_vs_FA.png    -- scatter: mean FA vs overestimation
"""

import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy import stats
import os

subj_dir = '/home/tijn-saes/Documents/Internship/ME_GRE'
out_dir  = os.path.join(subj_dir, 'figures')
os.makedirs(out_dir, exist_ok=True)

# ── Load ROI stats and compute overestimation ─────────────────────────────────
df = pd.read_csv(os.path.join(subj_dir, 'analysis', 'roi_stats.csv'))

df['overest_abs']     = df['MVF_basic_mean'] - df['MVF_atlas_mean']
df['overest_rel']     = (df['overest_abs'] / df['MVF_atlas_mean']) * 100
df['overest_abs_faw'] = df['MVF_basic_fa_weighted_mean'] - df['MVF_atlas_fa_weighted_mean']

# ── Add mean theta and FA per ROI from atlas maps ─────────────────────────────
def load(path):
    return np.array(nib.load(path).dataobj).astype(np.float32)

labels = load(f'{subj_dir}/atlas/JHU_labels_subj.nii.gz').astype(int)
theta  = load(f'{subj_dir}/atlas/theta_atlas.nii.gz')
fa     = load(f'{subj_dir}/atlas/FA_atlas.nii.gz')

theta_means, fa_means = [], []
for idx in df['ROI_index']:
    mask = labels == idx
    if mask.sum() == 0:
        theta_means.append(np.nan)
        fa_means.append(np.nan)
    else:
        theta_means.append(float(theta[mask].mean()))
        fa_means.append(float(fa[mask].mean()))

df['theta_mean'] = theta_means
df['fa_mean']    = fa_means

# Sort by overestimation for ranked plot
df_sorted = df.sort_values('overest_abs', ascending=True).reset_index(drop=True)

# ── Figure 23: Ranked overestimation bar chart ───────────────────────────────
fig, ax = plt.subplots(figsize=(10, 14), facecolor='#0d0d0d')
ax.set_facecolor('#0d0d0d')

colors = ['#e05c5c' if v > 0 else '#5c9ee0' for v in df_sorted['overest_abs']]
bars = ax.barh(range(len(df_sorted)), df_sorted['overest_abs'], color=colors, alpha=0.85)

# Annotate with relative %
for i, (abs_v, rel_v) in enumerate(zip(df_sorted['overest_abs'], df_sorted['overest_rel'])):
    x = abs_v + 0.002 if abs_v >= 0 else abs_v - 0.002
    ha = 'left' if abs_v >= 0 else 'right'
    if abs(abs_v) > 0.005:
        ax.text(x, i, f'{rel_v:+.0f}%', va='center', ha=ha,
                color='white', fontsize=7, alpha=0.85)

ax.set_yticks(range(len(df_sorted)))
ax.set_yticklabels(df_sorted['ROI_name'], color='white', fontsize=7.5)
ax.axvline(0, color='white', linewidth=0.8, alpha=0.5)
ax.set_xlabel('MVF_basic − MVF_atlas  (fraction)', color='white', fontsize=11)
ax.set_title('MVF Overestimation by Basic MIMM vs Atlas\nper JHU DTI-81 ROI  (red = basic overestimates)',
             color='white', fontsize=13, fontweight='bold', pad=12)
ax.tick_params(colors='white')
for spine in ax.spines.values(): spine.set_edgecolor('#444444')
ax.grid(axis='x', color='#333333', linewidth=0.5, linestyle='--')

plt.tight_layout()
plt.savefig(os.path.join(out_dir, '23_overestimation_ranked.png'),
            dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
plt.close()
print('Saved: 23_overestimation_ranked.png')

# ── Figure 24: Overestimation vs fibre angle θ ────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 7), facecolor='#0d0d0d')
ax.set_facecolor('#111111')

valid = df.dropna(subset=['theta_mean', 'overest_abs'])
sc = ax.scatter(valid['theta_mean'], valid['overest_abs'],
                c=valid['MVF_atlas_mean'], cmap='hot', s=60, alpha=0.85,
                vmin=0.05, vmax=0.40, zorder=3)

# Regression line
r, p = stats.pearsonr(valid['theta_mean'], valid['overest_abs'])
m, b = np.polyfit(valid['theta_mean'], valid['overest_abs'], 1)
x_line = np.linspace(valid['theta_mean'].min(), valid['theta_mean'].max(), 100)
ax.plot(x_line, m*x_line + b, '--', color='cyan', linewidth=1.5,
        label=f'r = {r:.2f}, p = {p:.3f}')

# Label notable points
highlight = ['Splenium of corpus callosum', 'Posterior thalamic radiation R',
             'Posterior thalamic radiation L', 'Posterior limb of internal capsule R',
             'Corticospinal tract R', 'Body of corpus callosum',
             'Superior longitudinal fasciculus R']
for _, row in valid.iterrows():
    if row['ROI_name'] in highlight:
        short = (row['ROI_name']
                 .replace('Posterior thalamic radiation', 'Post. thal. rad.')
                 .replace('Posterior limb of internal capsule', 'Post. IC')
                 .replace('Superior longitudinal fasciculus', 'SLF')
                 .replace(' R', '').replace(' L', ''))
        ax.annotate(short, (row['theta_mean'], row['overest_abs']),
                    textcoords='offset points', xytext=(6, 3),
                    fontsize=7.5, color='white', alpha=0.9,
                    path_effects=[pe.withStroke(linewidth=2, foreground='black')])

ax.axhline(0, color='white', linewidth=0.8, alpha=0.4)
cbar = fig.colorbar(sc, ax=ax, pad=0.02)
cbar.set_label('MVF_atlas (reference)', color='white', fontsize=9)
cbar.ax.yaxis.set_tick_params(color='white')
plt.setp(cbar.ax.yaxis.get_ticklabels(), color='white')
cbar.outline.set_edgecolor('white')

ax.set_xlabel('Mean fibre angle θ relative to B0  (degrees)', color='white', fontsize=11)
ax.set_ylabel('MVF_basic − MVF_atlas  (fraction)', color='white', fontsize=11)
ax.set_title('Overestimation vs Fibre Orientation\n(perpendicular fibres overestimated most)',
             color='white', fontsize=13, fontweight='bold')
ax.tick_params(colors='white')
for spine in ax.spines.values(): spine.set_edgecolor('#444444')
ax.grid(color='#2a2a2a', linewidth=0.5)
ax.legend(fontsize=10, facecolor='#1a1a1a', edgecolor='white', labelcolor='white')

plt.tight_layout()
plt.savefig(os.path.join(out_dir, '24_overestimation_vs_theta.png'),
            dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
plt.close()
print('Saved: 24_overestimation_vs_theta.png')

# ── Figure 25: Overestimation vs FA (FA = proxy for tract coherence) ──────────
fig, ax = plt.subplots(figsize=(8, 7), facecolor='#0d0d0d')
ax.set_facecolor('#111111')

sc = ax.scatter(valid['fa_mean'], valid['overest_abs'],
                c=valid['theta_mean'], cmap='RdBu_r', s=60, alpha=0.85,
                vmin=0, vmax=90, zorder=3)

r2, p2 = stats.pearsonr(valid['fa_mean'], valid['overest_abs'])
m2, b2 = np.polyfit(valid['fa_mean'], valid['overest_abs'], 1)
ax.plot(x_line := np.linspace(valid['fa_mean'].min(), valid['fa_mean'].max(), 100),
        m2*x_line + b2, '--', color='cyan', linewidth=1.5,
        label=f'r = {r2:.2f}, p = {p2:.3f}')

ax.axhline(0, color='white', linewidth=0.8, alpha=0.4)
cbar2 = fig.colorbar(sc, ax=ax, pad=0.02)
cbar2.set_label('Mean θ (degrees)', color='white', fontsize=9)
cbar2.ax.yaxis.set_tick_params(color='white')
plt.setp(cbar2.ax.yaxis.get_ticklabels(), color='white')
cbar2.outline.set_edgecolor('white')

ax.set_xlabel('Mean FA (HCP1065 atlas)', color='white', fontsize=11)
ax.set_ylabel('MVF_basic − MVF_atlas  (fraction)', color='white', fontsize=11)
ax.set_title('Overestimation vs Tract Coherence (FA)\n(colour = mean fibre angle)',
             color='white', fontsize=13, fontweight='bold')
ax.tick_params(colors='white')
for spine in ax.spines.values(): spine.set_edgecolor('#444444')
ax.grid(color='#2a2a2a', linewidth=0.5)
ax.legend(fontsize=10, facecolor='#1a1a1a', edgecolor='white', labelcolor='white')

plt.tight_layout()
plt.savefig(os.path.join(out_dir, '25_overestimation_vs_FA.png'),
            dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
plt.close()
print('Saved: 25_overestimation_vs_FA.png')

# ── Terminal summary ──────────────────────────────────────────────────────────
print(f'\nCorrelation: overestimation vs theta  r={r:.3f}  p={p:.4f}')
print(f'Correlation: overestimation vs FA     r={r2:.3f}  p={p2:.4f}')
print(f'\nTop 10 overestimated (basic - atlas):')
print(df.nlargest(10, 'overest_abs')[
    ['ROI_name', 'MVF_basic_mean', 'MVF_atlas_mean', 'overest_abs', 'overest_rel', 'theta_mean']
].to_string(index=False))
print(f'\nUnderestimated by basic (atlas > basic):')
print(df[df['overest_abs'] < -0.005][
    ['ROI_name', 'MVF_basic_mean', 'MVF_atlas_mean', 'overest_abs', 'overest_rel', 'theta_mean']
].to_string(index=False))
