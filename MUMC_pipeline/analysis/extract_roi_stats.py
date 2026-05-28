"""
Extract per-ROI statistics from MIMM and chi-separation maps using the JHU atlas.

Per map, per ROI:
  mean, SD, median, IQR, 5th/95th percentiles, CV, FA-weighted mean

Per ROI (map-independent):
  Pearson r between MVF_basic and chi_neg_chisep (voxel-level validation metric)

Post-loop:
  Laterality index for all paired R/L tracts (labels 7–50)

Saves results to subj_dir/analysis/roi_stats.csv
             and subj_dir/analysis/laterality.csv
"""

import numpy as np
import nibabel as nib
import pandas as pd
from scipy import stats
import os

subj_dir = '/home/tijn-saes/Documents/Internship/ME_GRE'
out_dir  = os.path.join(subj_dir, 'analysis')
os.makedirs(out_dir, exist_ok=True)

# --- JHU label names (index 1–50) ---
# Labels 1–6 are midline structures; 7–50 are paired R (odd) / L (even)
JHU_LABELS = {
     1: 'Middle cerebellar peduncle',
     2: 'Pontine crossing tract',
     3: 'Genu of corpus callosum',
     4: 'Body of corpus callosum',
     5: 'Splenium of corpus callosum',
     6: 'Fornix (column and body)',
     7: 'Corticospinal tract R',
     8: 'Corticospinal tract L',
     9: 'Medial lemniscus R',
    10: 'Medial lemniscus L',
    11: 'Inferior cerebellar peduncle R',
    12: 'Inferior cerebellar peduncle L',
    13: 'Superior cerebellar peduncle R',
    14: 'Superior cerebellar peduncle L',
    15: 'Cerebral peduncle R',
    16: 'Cerebral peduncle L',
    17: 'Anterior limb of internal capsule R',
    18: 'Anterior limb of internal capsule L',
    19: 'Posterior limb of internal capsule R',
    20: 'Posterior limb of internal capsule L',
    21: 'Retrolenticular internal capsule R',
    22: 'Retrolenticular internal capsule L',
    23: 'Anterior corona radiata R',
    24: 'Anterior corona radiata L',
    25: 'Superior corona radiata R',
    26: 'Superior corona radiata L',
    27: 'Posterior corona radiata R',
    28: 'Posterior corona radiata L',
    29: 'Posterior thalamic radiation R',
    30: 'Posterior thalamic radiation L',
    31: 'Sagittal stratum R',
    32: 'Sagittal stratum L',
    33: 'External capsule R',
    34: 'External capsule L',
    35: 'Cingulum (cingulate gyrus) R',
    36: 'Cingulum (cingulate gyrus) L',
    37: 'Cingulum (hippocampus) R',
    38: 'Cingulum (hippocampus) L',
    39: 'Fornix / stria terminalis R',
    40: 'Fornix / stria terminalis L',
    41: 'Superior longitudinal fasciculus R',
    42: 'Superior longitudinal fasciculus L',
    43: 'Superior fronto-occipital fasciculus R',
    44: 'Superior fronto-occipital fasciculus L',
    45: 'Uncinate fasciculus R',
    46: 'Uncinate fasciculus L',
    47: 'Tapetum R',
    48: 'Tapetum L',
    49: 'Fornix (cres) R',
    50: 'Fornix (cres) L',
}

# Paired R/L label indices for DTI-81: (R_label, L_label, tract_name)
LATERAL_PAIRS = [
    (idx, idx + 1, JHU_LABELS[idx].replace(' R', ''))
    for idx in range(7, 50, 2)
]

# --- JHU tractography atlas label names (NIfTI values 1–20, thr25) ---
JHU_TRACTS = {
     1: 'Anterior thalamic radiation L',
     2: 'Anterior thalamic radiation R',
     3: 'Corticospinal tract L',
     4: 'Corticospinal tract R',
     5: 'Cingulum (cingulate gyrus) L',
     6: 'Cingulum (cingulate gyrus) R',
     7: 'Cingulum (hippocampus) L',
     8: 'Cingulum (hippocampus) R',
     9: 'Forceps major',
    10: 'Forceps minor',
    11: 'Inferior fronto-occipital fasciculus L',
    12: 'Inferior fronto-occipital fasciculus R',
    13: 'Inferior longitudinal fasciculus L',
    14: 'Inferior longitudinal fasciculus R',
    15: 'Superior longitudinal fasciculus L',
    16: 'Superior longitudinal fasciculus R',
    17: 'Uncinate fasciculus L',
    18: 'Uncinate fasciculus R',
    19: 'Superior longitudinal fasciculus (temporal) L',
    20: 'Superior longitudinal fasciculus (temporal) R',
}


def load(path):
    return np.array(nib.load(path).dataobj).astype(np.float32)


def roi_stats(vals, fa_weights):
    """Compute all per-map statistics for a 1D array of voxel values."""
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return {k: np.nan for k in
                ['mean', 'sd', 'median', 'iqr', 'p05', 'p95', 'cv', 'fa_weighted_mean']}
    mean = float(np.mean(vals))
    sd   = float(np.std(vals))
    return {
        'mean':           mean,
        'sd':             sd,
        'median':         float(np.median(vals)),
        'iqr':            float(np.percentile(vals, 75) - np.percentile(vals, 25)),
        'p05':            float(np.percentile(vals, 5)),
        'p95':            float(np.percentile(vals, 95)),
        'cv':             float(sd / mean) if mean != 0 else np.nan,
        'fa_weighted_mean': float(np.average(vals, weights=fa_weights[np.isfinite(vals)])),
    }


print('Loading maps...')
labels = load(f'{subj_dir}/atlas/JHU_labels_subj.nii.gz').astype(int)
fa_dti   = f'{subj_dir}/dti/FA.nii.gz'
fa_atlas = f'{subj_dir}/atlas/FA_atlas.nii.gz'
fa = load(fa_dti if os.path.exists(fa_dti) else fa_atlas)

maps = {
    'MVF_basic':        load(f'{subj_dir}/mimm/MVF_basic.nii.gz'),
    'MVF_atlas':        load(f'{subj_dir}/mimm/MVF_Atlas.nii.gz'),
    'FVF_basic':        load(f'{subj_dir}/mimm/FVF_basic.nii.gz'),
    'FVF_atlas':        load(f'{subj_dir}/mimm/FVF_Atlas.nii.gz'),
    'g_ratio_basic':    load(f'{subj_dir}/mimm/g_ratio_basic.nii.gz'),
    'g_ratio_atlas':    load(f'{subj_dir}/mimm/g_ratio_Atlas.nii.gz'),
    'R2s_basic':        load(f'{subj_dir}/mimm/R2s_basic.nii.gz'),
    'R2s_atlas':        load(f'{subj_dir}/mimm/R2s_Atlas.nii.gz'),
    'chi_myelin_basic': np.abs(load(f'{subj_dir}/mimm/chi_myelin_basic.nii.gz')),
    'chi_myelin_atlas': np.abs(load(f'{subj_dir}/mimm/chi_myelin_Atlas.nii.gz')),
    'chi_iron_basic':   load(f'{subj_dir}/mimm/chi_iron_est_basic.nii.gz'),
    'chi_iron_atlas':   load(f'{subj_dir}/mimm/chi_iron_est_Atlas.nii.gz'),
    'chi_neg_chisep':   load(f'{subj_dir}/chisep/chi_neg.nii.gz'),
    'chi_pos_chisep':   load(f'{subj_dir}/chisep/chi_pos.nii.gz'),
}

# -------------------------------------------------------------------------
# Main loop: per-ROI statistics
# -------------------------------------------------------------------------
print('Extracting ROI statistics...')
rows = []

for idx, name in sorted(JHU_LABELS.items()):
    roi_mask = labels == idx
    n_vox = int(roi_mask.sum())
    if n_vox == 0:
        continue

    fa_weights = fa[roi_mask].clip(min=0)   # negative FA values are artefacts

    row = {'ROI_index': idx, 'ROI_name': name, 'n_voxels': n_vox}

    for map_name, vol in maps.items():
        s = roi_stats(vol[roi_mask], fa_weights)
        for stat_name, val in s.items():
            row[f'{map_name}_{stat_name}'] = val

    # Within-ROI Pearson r: MVF_basic vs chi_neg_chisep
    mvf  = maps['MVF_basic'][roi_mask]
    cneg = maps['chi_neg_chisep'][roi_mask]
    finite = np.isfinite(mvf) & np.isfinite(cneg)
    if finite.sum() > 2:
        r, p = stats.pearsonr(mvf[finite], cneg[finite])
        row['r_MVF_vs_chineg'] = float(r)
        row['p_MVF_vs_chineg'] = float(p)
    else:
        row['r_MVF_vs_chineg'] = np.nan
        row['p_MVF_vs_chineg'] = np.nan

    rows.append(row)

df = pd.DataFrame(rows)
out_csv = os.path.join(out_dir, 'roi_stats.csv')
df.to_csv(out_csv, index=False, float_format='%.5f')
print(f'Saved: {out_csv}  ({len(rows)} ROIs)')

# -------------------------------------------------------------------------
# Laterality index: LI = (R - L) / (R + L) for paired tracts
# -------------------------------------------------------------------------
print('Computing laterality indices...')
lat_rows = []

for r_idx, l_idx, tract_name in LATERAL_PAIRS:
    r_row = df[df['ROI_index'] == r_idx]
    l_row = df[df['ROI_index'] == l_idx]
    if r_row.empty or l_row.empty:
        continue

    lat_row = {'tract': tract_name}
    for map_name in maps:
        r_mean = r_row[f'{map_name}_mean'].values[0]
        l_mean = l_row[f'{map_name}_mean'].values[0]
        denom  = r_mean + l_mean
        lat_row[f'{map_name}_LI'] = float((r_mean - l_mean) / denom) if denom != 0 else np.nan
    lat_rows.append(lat_row)

lat_df = pd.DataFrame(lat_rows)
lat_csv = os.path.join(out_dir, 'laterality.csv')
lat_df.to_csv(lat_csv, index=False, float_format='%.5f')
print(f'Saved: {lat_csv}  ({len(lat_rows)} tract pairs)')

# -------------------------------------------------------------------------
# Terminal summary
# -------------------------------------------------------------------------
print(f'\nTop 10 ROIs by FA-weighted MVF (basic):')
print(df.nlargest(10, 'MVF_basic_fa_weighted_mean')[
    ['ROI_name', 'MVF_basic_fa_weighted_mean', 'MVF_atlas_fa_weighted_mean',
     'g_ratio_basic_mean', 'r_MVF_vs_chineg']
].to_string(index=False))

print(f'\nMost asymmetric tracts (|LI| for MVF_basic):')
lat_df['MVF_basic_LI_abs'] = lat_df['MVF_basic_LI'].abs()
print(lat_df.nlargest(5, 'MVF_basic_LI_abs')[['tract', 'MVF_basic_LI']].to_string(index=False))

# -------------------------------------------------------------------------
# JHU Tractography atlas (thr25) — 20 tracts
# -------------------------------------------------------------------------
tracts_path = f'{subj_dir}/atlas/JHU_tracts_subj.nii.gz'
if os.path.exists(tracts_path):
    print('\n--- JHU Tractography atlas (thr25, 20 tracts) ---')
    tract_labels = load(tracts_path).astype(int)
    tract_rows = []
    for idx, name in sorted(JHU_TRACTS.items()):
        roi_mask = tract_labels == idx
        n_vox = int(roi_mask.sum())
        if n_vox == 0:
            continue
        fa_weights = fa[roi_mask].clip(min=0)
        row = {'ROI_index': idx, 'ROI_name': name, 'n_voxels': n_vox}
        for map_name, vol in maps.items():
            s = roi_stats(vol[roi_mask], fa_weights)
            for stat_name, val in s.items():
                row[f'{map_name}_{stat_name}'] = val
        mvf  = maps['MVF_basic'][roi_mask]
        cneg = maps['chi_neg_chisep'][roi_mask]
        finite = np.isfinite(mvf) & np.isfinite(cneg)
        if finite.sum() > 2:
            r, p = stats.pearsonr(mvf[finite], cneg[finite])
            row['r_MVF_vs_chineg'] = float(r)
            row['p_MVF_vs_chineg'] = float(p)
        else:
            row['r_MVF_vs_chineg'] = np.nan
            row['p_MVF_vs_chineg'] = np.nan
        tract_rows.append(row)

    tract_df = pd.DataFrame(tract_rows)
    tract_csv = os.path.join(out_dir, 'roi_stats_tracts.csv')
    tract_df.to_csv(tract_csv, index=False, float_format='%.5f')
    print(f'Saved: {tract_csv}  ({len(tract_rows)} tracts)')
    print(tract_df.nlargest(10, 'MVF_basic_fa_weighted_mean')[
        ['ROI_name', 'MVF_basic_fa_weighted_mean', 'MVF_atlas_fa_weighted_mean',
         'g_ratio_basic_mean', 'r_MVF_vs_chineg', 'n_voxels']
    ].to_string(index=False))

    # ---- Atlas comparison: within-ROI CV (lower = more homogeneous) ----
    print('\n--- Atlas comparison: mean within-ROI CV for MVF_basic ---')
    print(f'  DTI-81  (48 ROIs, {int(df["n_voxels"].mean())} vox/ROI avg): '
          f'CV = {df["MVF_basic_cv"].mean():.3f}')
    print(f'  Tracts  (20 ROIs, {int(tract_df["n_voxels"].mean())} vox/ROI avg): '
          f'CV = {tract_df["MVF_basic_cv"].mean():.3f}')
    print('  (Lower CV = more homogeneous ROIs = better atlas for microstructure)')
else:
    print('\nJHU tractography atlas not found — skipping (run register_atlas.sh first)')