"""
Subject-level Bland-Altman plots: MIMM MVF vs chi-separation chi_neg per JHU ROI.

Each dot = one subject's ROI mean.  No voxel-level pooling.

Two outputs:
  bland_altman_grid.png  — grid of subplots, one panel per selected ROI
  bland_altman_all.png   — single plot with all ROIs overlaid, color-coded by ROI

chi_neg (ppm) is scaled to MVF units via OLS before differencing,
matching the approach in bland_altman.m.

Usage: edit SUBJ_DIRS and OUT_DIR at the top, then run.
"""

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
from scipy import stats
import os

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SUBJ_DIRS = [
    # '/data/MUMC/sub-01',
    # '/data/MUMC/sub-02',
    # ...
]

OUT_DIR      = '/home/tijn-saes/Documents/Internship/ME_GRE/analysis'
MIMM_VARIANT = 'basic'   # 'basic' or 'Atlas'

# JHU ROIs to show in the grid (label index → display name).
# Selected as most relevant for WM research; edit freely.
ROIS_TO_PLOT = {
     3: 'Genu CC',
     4: 'Body CC',
     5: 'Splenium CC',
     7: 'CST R',
     8: 'CST L',
    19: 'Post. IC R',
    20: 'Post. IC L',
    25: 'Sup. CR R',
    26: 'Sup. CR L',
    41: 'SLF R',
    42: 'SLF L',
    35: 'Cingulum R',
}

# ---------------------------------------------------------------------------
# JHU full label dict (used for titles)
# ---------------------------------------------------------------------------
JHU_LABELS = {
     1: 'Middle cerebellar peduncle',   2: 'Pontine crossing tract',
     3: 'Genu of corpus callosum',      4: 'Body of corpus callosum',
     5: 'Splenium of corpus callosum',  6: 'Fornix (column and body)',
     7: 'Corticospinal tract R',        8: 'Corticospinal tract L',
     9: 'Medial lemniscus R',          10: 'Medial lemniscus L',
    11: 'Inf. cerebellar peduncle R',  12: 'Inf. cerebellar peduncle L',
    13: 'Sup. cerebellar peduncle R',  14: 'Sup. cerebellar peduncle L',
    15: 'Cerebral peduncle R',         16: 'Cerebral peduncle L',
    17: 'Ant. IC R',                   18: 'Ant. IC L',
    19: 'Post. IC R',                  20: 'Post. IC L',
    21: 'Retrolenticular IC R',        22: 'Retrolenticular IC L',
    23: 'Ant. corona radiata R',       24: 'Ant. corona radiata L',
    25: 'Sup. corona radiata R',       26: 'Sup. corona radiata L',
    27: 'Post. corona radiata R',      28: 'Post. corona radiata L',
    29: 'Post. thalamic radiation R',  30: 'Post. thalamic radiation L',
    31: 'Sagittal stratum R',          32: 'Sagittal stratum L',
    33: 'External capsule R',          34: 'External capsule L',
    35: 'Cingulum (cingulate) R',      36: 'Cingulum (cingulate) L',
    37: 'Cingulum (hippocampus) R',    38: 'Cingulum (hippocampus) L',
    39: 'Fornix/stria terminalis R',   40: 'Fornix/stria terminalis L',
    41: 'SLF R',                       42: 'SLF L',
    43: 'Sup. fronto-occ. fasc. R',    44: 'Sup. fronto-occ. fasc. L',
    45: 'Uncinate fasciculus R',       46: 'Uncinate fasciculus L',
    47: 'Tapetum R',                   48: 'Tapetum L',
    49: 'Fornix (cres) R',             50: 'Fornix (cres) L',
}
## there is another atlas 

def load(path):
    return np.array(nib.load(path).dataobj).astype(np.float32)


def subject_roi_means(subj_dir):
    """
    For one subject: return dict of {roi_idx: (mvf_mean, chineg_mean)}.
    Returns None if required files are missing.
    """
    mvf_file    = os.path.join(subj_dir, 'mimm',  f'MVF_{MIMM_VARIANT}.nii.gz')
    chineg_file = os.path.join(subj_dir, 'chisep', 'chi_neg.nii.gz')
    labels_file = os.path.join(subj_dir, 'atlas',  'JHU_labels_subj.nii.gz')

    for f in [mvf_file, chineg_file, labels_file]:
        if not os.path.exists(f):
            print(f'  Missing: {f} — skipping subject')
            return None

    mvf    = load(mvf_file)
    chineg = np.abs(load(chineg_file))   # diamagnetic → positive
    labels = load(labels_file).astype(int)

    roi_means = {}
    for idx in ROIS_TO_PLOT:
        mask = labels == idx
        if mask.sum() == 0:
            continue
        mvf_vals    = mvf[mask];    mvf_vals    = mvf_vals[np.isfinite(mvf_vals)]
        chineg_vals = chineg[mask]; chineg_vals = chineg_vals[np.isfinite(chineg_vals)]
        if len(mvf_vals) == 0 or len(chineg_vals) == 0:
            continue
        roi_means[idx] = (float(np.mean(mvf_vals)), float(np.mean(chineg_vals)))

    return roi_means


def ols_scale(x, y):
    """Scale y to x units via OLS: x = a + b*y.  Returns scaled y."""
    A = np.column_stack([np.ones_like(y), y])
    coef, _, _, _ = np.linalg.lstsq(A, x, rcond=None)
    return coef[0] + coef[1] * y, coef


def ba_stats(x, y_scaled):
    diff = x - y_scaled
    bias = float(np.mean(diff))
    sd   = float(np.std(diff, ddof=1))
    r, p = stats.pearsonr(x, y_scaled)
    return {
        'bias': bias, 'sd': sd,
        'loa_upper': bias + 1.96 * sd,
        'loa_lower': bias - 1.96 * sd,
        'r': float(r), 'p': float(p),
        'n': len(x),
    }


# ---------------------------------------------------------------------------
# Collect subject-level means
# ---------------------------------------------------------------------------
os.makedirs(OUT_DIR, exist_ok=True)

print(f'Processing {len(SUBJ_DIRS)} subjects...')
records = []   # list of {subj, roi_idx, mvf_mean, chineg_mean}

for subj_dir in SUBJ_DIRS:
    subj_id = os.path.basename(subj_dir)
    print(f'  {subj_id}')
    means = subject_roi_means(subj_dir)
    if means is None:
        continue
    for roi_idx, (mvf_mean, chineg_mean) in means.items():
        records.append({
            'subject':     subj_id,
            'roi_idx':     roi_idx,
            'roi_name':    ROIS_TO_PLOT[roi_idx],
            'MVF':         mvf_mean,
            'chi_neg_ppm': chineg_mean,
        })

df = pd.DataFrame(records)
df.to_csv(os.path.join(OUT_DIR, 'subject_roi_means.csv'), index=False, float_format='%.5f')
print(f'Collected {len(df)} subject × ROI observations from {df["subject"].nunique()} subjects.')

if df.empty or df['subject'].nunique() < 3:
    print('Need at least 3 subjects to generate Bland-Altman plots.')
    raise SystemExit(0)

# ---------------------------------------------------------------------------
# Grid plot: one Bland-Altman panel per ROI
# ---------------------------------------------------------------------------
n_rois = len(ROIS_TO_PLOT)
n_cols = 4
n_rows = int(np.ceil(n_rois / n_cols))

fig, axes = plt.subplots(n_rows, n_cols,
                         figsize=(n_cols * 4.5, n_rows * 4),
                         facecolor='black')
fig.suptitle(f'Bland-Altman: MIMM MVF ({MIMM_VARIANT}) vs |χ_neg| (chi-sep)\n'
             f'One dot = one subject  |  n = {df["subject"].nunique()} subjects',
             color='white', fontsize=13, fontweight='bold', y=1.01)

axes_flat = axes.ravel() if n_rows > 1 else [axes] if n_cols == 1 else axes.ravel()
summary_rows = []

for ax_i, (roi_idx, roi_short) in enumerate(ROIS_TO_PLOT.items()):
    ax = axes_flat[ax_i]
    ax.set_facecolor('#111111')

    roi_df = df[df['roi_idx'] == roi_idx].dropna()
    if len(roi_df) < 3:
        ax.set_title(roi_short, color='gray', fontsize=9)
        ax.text(0.5, 0.5, 'n < 3', color='gray', ha='center', va='center',
                transform=ax.transAxes)
        ax.axis('off')
        continue

    mvf    = roi_df['MVF'].values
    chineg = roi_df['chi_neg_ppm'].values

    chineg_scaled, coef = ols_scale(mvf, chineg)
    s = ba_stats(mvf, chineg_scaled)
    summary_rows.append({'roi': JHU_LABELS.get(roi_idx, roi_short), **s})

    mean_xy = (mvf + chineg_scaled) / 2
    diff_xy = mvf - chineg_scaled

    ax.scatter(mean_xy, diff_xy,
               s=60, color='#4a9eda', edgecolors='white', linewidths=0.4, zorder=3)

    ax.axhline(s['bias'],      color='cyan',   linewidth=1.5, linestyle='--', zorder=2)
    ax.axhline(s['loa_upper'], color='#ffcc44', linewidth=1.2, linestyle=':',  zorder=2)
    ax.axhline(s['loa_lower'], color='#ffcc44', linewidth=1.2, linestyle=':',  zorder=2)
    ax.axhline(0,              color='white',   linewidth=0.6, alpha=0.3,      zorder=1)

    # Annotate stats in corner
    ax.text(0.97, 0.97,
            f'r={s["r"]:.2f}  bias={s["bias"]:.3f}\nLoA [{s["loa_lower"]:.3f}, {s["loa_upper"]:.3f}]',
            transform=ax.transAxes, ha='right', va='top',
            fontsize=7, color='white',
            bbox=dict(facecolor='#222222', edgecolor='none', alpha=0.7))

    ax.set_title(f'{roi_short}  (n={s["n"]})', color='white', fontsize=9)
    ax.set_xlabel('Mean (MVF, scaled χ_neg)', color='white', fontsize=7)
    ax.set_ylabel('MVF − scaled χ_neg', color='white', fontsize=7)
    ax.tick_params(colors='white', labelsize=7)
    for spine in ax.spines.values():
        spine.set_edgecolor('#444444')

# Hide unused panels
for ax_i in range(n_rois, len(axes_flat)):
    axes_flat[ax_i].axis('off')

plt.tight_layout()
grid_path = os.path.join(OUT_DIR, 'bland_altman_grid.png')
plt.savefig(grid_path, dpi=150, bbox_inches='tight', facecolor='black', edgecolor='none')
plt.close()
print(f'Saved: bland_altman_grid.png')

# ---------------------------------------------------------------------------
# Summary plot: all ROIs overlaid, color-coded by ROI
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 7), facecolor='black')
ax.set_facecolor('#111111')

roi_list   = sorted(ROIS_TO_PLOT.keys())
colors     = cm.tab20(np.linspace(0, 1, len(roi_list)))
legend_handles = []

for color, roi_idx in zip(colors, roi_list):
    roi_df = df[df['roi_idx'] == roi_idx].dropna()
    if len(roi_df) < 3:
        continue
    mvf    = roi_df['MVF'].values
    chineg = roi_df['chi_neg_ppm'].values
    chineg_scaled, _ = ols_scale(mvf, chineg)

    mean_xy = (mvf + chineg_scaled) / 2
    diff_xy = mvf - chineg_scaled

    h = ax.scatter(mean_xy, diff_xy, s=50, color=color,
                   edgecolors='white', linewidths=0.3, alpha=0.85, zorder=3,
                   label=ROIS_TO_PLOT[roi_idx])
    legend_handles.append(h)

ax.axhline(0, color='white', linewidth=0.8, alpha=0.4, zorder=1)

ax.set_xlabel('Mean of MIMM MVF and scaled |χ_neg|', color='white', fontsize=11)
ax.set_ylabel('MIMM MVF − scaled |χ_neg|', color='white', fontsize=11)
ax.set_title(f'Bland-Altman: MIMM MVF vs |χ_neg|  —  all ROIs\n'
             f'One dot = one subject  |  n = {df["subject"].nunique()} subjects  |  variant: {MIMM_VARIANT}',
             color='white', fontsize=12, fontweight='bold')
ax.tick_params(colors='white')
for spine in ax.spines.values():
    spine.set_edgecolor('#444444')
ax.legend(handles=legend_handles, fontsize=7, ncol=2,
          facecolor='#1a1a1a', edgecolor='#444444', labelcolor='white',
          loc='upper right')

plt.tight_layout()
all_path = os.path.join(OUT_DIR, 'bland_altman_all.png')
plt.savefig(all_path, dpi=150, bbox_inches='tight', facecolor='black', edgecolor='none')
plt.close()
print(f'Saved: bland_altman_all.png')

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
if summary_rows:
    summary_df = pd.DataFrame(summary_rows)
    print(f'\nBland-Altman summary (n={df["subject"].nunique()} subjects):')
    print(summary_df[['roi', 'n', 'bias', 'sd', 'loa_lower', 'loa_upper', 'r', 'p']
                     ].to_string(index=False, float_format=lambda x: f'{x:.4f}'))
    summary_df.to_csv(os.path.join(OUT_DIR, 'bland_altman_summary.csv'),
                      index=False, float_format='%.5f')
    print(f'\nSaved: bland_altman_summary.csv')
