"""
Overestimation analysis: MIMM Basic vs Atlas MVF per JHU ROI.

Tests the core physics hypothesis:
  Basic MIMM overestimates MVF most in tracts running perpendicular to B0,
  because without an orientation prior it confuses angle-dependent signal
  decay with myelin content.

Produces:
  23_overestimation_ranked.png        -- all ROIs ranked by overestimation
  24_overestimation_vs_theta.png      -- scatter: mean fibre angle vs overestimation
  25_overestimation_vs_FA.png         -- scatter: mean FA vs overestimation
  45_overestimation_spatial_JHU.png   -- spatial map with JHU boundaries overlaid
"""

import numpy as np
import nibabel as nib
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy import stats
import os

try:
    from paths import OUTPUT_DIR as subj_dir, FIG_DIR as out_dir, ANALYSIS_DIR
except ImportError:
    raise SystemExit('Copy MUMC_pipeline/analysis/paths_template.py to paths.py and fill in your paths.')
os.makedirs(out_dir, exist_ok=True)

# ── Load ROI stats and compute overestimation ─────────────────────────────────
df = pd.read_csv(os.path.join(ANALYSIS_DIR, 'roi_stats.csv'))

df['overest_abs']     = df['MVF_basic_mean'] - df['MVF_atlas_mean']
df['overest_rel']     = (df['overest_abs'] / df['MVF_atlas_mean']) * 100
df['overest_abs_faw'] = df['MVF_basic_fa_weighted_mean'] - df['MVF_atlas_fa_weighted_mean']

# ── Add mean theta and FA per ROI from atlas maps ─────────────────────────────
def load(path):
    return np.array(nib.load(path).dataobj).astype(np.float32)

labels = load(os.path.join(subj_dir, 'atlas', 'JHU_labels_subj.nii.gz')).astype(int)
theta  = load(os.path.join(subj_dir, 'atlas', 'theta_atlas.nii.gz'))
fa     = load(os.path.join(subj_dir, 'atlas', 'FA_atlas.nii.gz'))

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

# ── Figure 45: Spatial overestimation with JHU atlas overlay ─────────────────
# Two-row × three-column figure.
# Row 1: voxel-wise MVF difference (basic − atlas) — continuous spatial map.
# Row 2: ROI-mean overestimation projected back to each voxel — shows the
#         structured, atlas-level result in anatomical context.
# Both rows share the same diverging colormap and scale. JHU ROI boundaries
# are overlaid as thin white contour lines. Top-5 overestimated ROIs are
# labeled by name on the axial slice of row 2.

def _load(path):
    return np.array(nib.load(path).dataobj).astype(np.float32)

mvf_b = _load(os.path.join(subj_dir, 'mimm', 'MVF_basic.nii.gz'))
mvf_a = _load(os.path.join(subj_dir, 'mimm', 'MVF_Atlas.nii.gz'))
mag   = _load(os.path.join(subj_dir, 'qsm',  'mag_e1.nii.gz'))
brain = _load(os.path.join(subj_dir, 'qsm',  'brain_mask.nii.gz')) > 0

mvf_diff = mvf_b - mvf_a
mvf_diff[~brain] = 0.0

# Build ROI-mean overestimation volume: each voxel gets its ROI's mean overest.
overest_vol = np.full(mvf_diff.shape, np.nan, dtype=np.float32)
for _, row in df.iterrows():
    overest_vol[labels == row['ROI_index']] = row['overest_abs']
overest_vol[~brain] = np.nan   # transparent outside brain

# Symmetric colormap range: clip to 99th percentile of absolute values
_vals = mvf_diff[brain]
_clim = float(np.percentile(np.abs(_vals), 99))
_clim = max(_clim, 0.01)   # avoid zero range

# Slices centred on brain-mask centroid
_coords      = np.argwhere(brain)
cx, cy, cz   = _coords.mean(axis=0).astype(int)
view_labels  = ['Axial', 'Coronal', 'Sagittal']

def _s(vol):
    """Return (axial, coronal, sagittal) slices, rotated for display."""
    return (np.rot90(vol[:, :, cz]),
            np.rot90(vol[:, cy, :]),
            np.rot90(vol[cx, :, :]))

lbl_s     = _s(labels.astype(float))
diff_s    = _s(mvf_diff)
overest_s = _s(overest_vol)
mag_s     = _s(mag)

# Unique label levels for contour boundaries
_lev = np.unique(labels[brain])
_lev = _lev[_lev > 0] - 0.5   # boundary between each label and next

fig, axes = plt.subplots(2, 3, figsize=(16, 10), facecolor='black')
fig.suptitle('MVF Overestimation in Brain Space — MIMM Basic vs Atlas\n'
             'with JHU DTI-81 ROI Boundaries',
             color='white', fontsize=14, fontweight='bold', y=1.01)

row_titles = ['Voxel-wise MVF diff  (Basic − Atlas)',
              'ROI-mean overestimation  (colour = ROI mean Basic − Atlas)']
cmap       = 'RdBu_r'   # red = basic overestimates, blue = underestimates

for row, (slices, row_title) in enumerate(zip([diff_s, overest_s], row_titles)):
    for col in range(3):
        ax = axes[row, col]
        ax.set_facecolor('black')

        # Anatomy background
        bg = mag_s[col]
        bg_norm = (bg - bg.min()) / (bg.max() - bg.min() + 1e-8)
        ax.imshow(bg_norm, cmap='gray', aspect='auto', interpolation='nearest')

        # Overestimation overlay (masked where nan/zero-brain)
        sl = slices[col]
        masked = np.ma.masked_invalid(sl) if row == 1 else np.ma.masked_where(~np.isfinite(sl), sl)
        ax.imshow(masked, cmap=cmap, vmin=-_clim, vmax=_clim,
                  aspect='auto', interpolation='nearest', alpha=0.80)

        # JHU ROI boundaries
        if len(_lev) > 0:
            ax.contour(lbl_s[col], levels=_lev,
                       colors='white', linewidths=0.35, alpha=0.55)

        if row == 0:
            ax.set_title(view_labels[col], color='white', fontsize=11)
        ax.axis('off')

    axes[row, 0].set_ylabel(row_title, color='white', fontsize=10,
                            rotation=90, labelpad=8)

# Label top-5 overestimated ROIs on axial slice of row 2 (most informative view)
top5 = df.nlargest(5, 'overest_abs')
ax_lbl = axes[1, 0]   # axial slice of overestimation row
for _, row_r in top5.iterrows():
    roi_mask_ax = np.rot90(labels[:, :, cz] == row_r['ROI_index'])
    if roi_mask_ax.sum() == 0:
        continue
    ys, xs = np.where(roi_mask_ax)
    xc, yc = xs.mean(), ys.mean()
    short = (row_r['ROI_name']
             .replace('Posterior thalamic radiation', 'Post. thal. rad.')
             .replace('Posterior limb of internal capsule', 'Post. IC')
             .replace('Retrolenticular internal capsule', 'Retrolent. IC')
             .replace('Splenium of corpus callosum', 'Splenium CC')
             .replace('Body of corpus callosum', 'Body CC')
             .replace(' R', '').replace(' L', ''))
    ax_lbl.annotate(f'{short}\n{row_r["overest_rel"]:+.0f}%',
                    xy=(xc, yc), xytext=(xc + 12, yc - 12),
                    fontsize=6.5, color='white', ha='left',
                    path_effects=[pe.withStroke(linewidth=2, foreground='black')],
                    arrowprops=dict(arrowstyle='->', color='white', lw=0.8))

# Shared colorbar
sm = plt.cm.ScalarMappable(cmap=cmap,
                           norm=plt.Normalize(vmin=-_clim, vmax=_clim))
sm.set_array([])
cbar = fig.colorbar(sm, ax=axes, orientation='vertical',
                    fraction=0.015, pad=0.02, shrink=0.85)
cbar.set_label('MVF_basic − MVF_atlas  (fraction)', color='white', fontsize=10)
cbar.ax.yaxis.set_tick_params(color='white')
plt.setp(cbar.ax.yaxis.get_ticklabels(), color='white')
cbar.outline.set_edgecolor('white')

plt.tight_layout()
plt.savefig(os.path.join(out_dir, '45_overestimation_spatial_JHU.png'),
            dpi=150, bbox_inches='tight', facecolor='black', edgecolor='none')
plt.close()
print('Saved: 45_overestimation_spatial_JHU.png')

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
