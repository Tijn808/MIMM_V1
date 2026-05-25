"""
Extract mean ± SD of MIMM and chi-separation maps per JHU white matter ROI.
Saves results to subj_dir/analysis/roi_stats.csv.
"""

import numpy as np
import nibabel as nib
import pandas as pd
import os

subj_dir  = '/home/tijn-saes/Documents/Internship/ME_GRE'
out_dir   = os.path.join(subj_dir, 'analysis')
os.makedirs(out_dir, exist_ok=True)

# --- JHU label names (index 1–50) ---
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

def load(path):
    return np.array(nib.load(path).dataobj).astype(np.float32)

print('Loading maps...')
labels = load(f'{subj_dir}/atlas/JHU_labels_subj.nii.gz').astype(int)

maps = {
    'MVF_basic':            load(f'{subj_dir}/mimm/MVF_basic.nii.gz'),
    'MVF_atlas':            load(f'{subj_dir}/mimm/MVF_Atlas.nii.gz'),
    'FVF_basic':            load(f'{subj_dir}/mimm/FVF_basic.nii.gz'),
    'FVF_atlas':            load(f'{subj_dir}/mimm/FVF_Atlas.nii.gz'),
    'g_ratio_basic':        load(f'{subj_dir}/mimm/g_ratio_basic.nii.gz'),
    'g_ratio_atlas':        load(f'{subj_dir}/mimm/g_ratio_Atlas.nii.gz'),
    'R2s_basic':            load(f'{subj_dir}/mimm/R2s_basic.nii.gz'),
    'R2s_atlas':            load(f'{subj_dir}/mimm/R2s_Atlas.nii.gz'),
    'chi_myelin_basic':     np.abs(load(f'{subj_dir}/mimm/chi_myelin_basic.nii.gz')),
    'chi_myelin_atlas':     np.abs(load(f'{subj_dir}/mimm/chi_myelin_Atlas.nii.gz')),
    'chi_iron_basic':       load(f'{subj_dir}/mimm/chi_iron_est_basic.nii.gz'),
    'chi_iron_atlas':       load(f'{subj_dir}/mimm/chi_iron_est_Atlas.nii.gz'),
    'chi_neg_chisep':       load(f'{subj_dir}/chisep/chi_neg.nii.gz'),
    'chi_pos_chisep':       load(f'{subj_dir}/chisep/chi_pos.nii.gz'),
}

print('Extracting ROI statistics...')
rows = []
for idx, name in sorted(JHU_LABELS.items()):
    roi_mask = labels == idx
    n_vox = roi_mask.sum()
    if n_vox == 0:
        continue
    row = {'ROI_index': idx, 'ROI_name': name, 'n_voxels': n_vox}
    for map_name, vol in maps.items():
        vals = vol[roi_mask]
        row[f'{map_name}_mean'] = float(np.mean(vals))
        row[f'{map_name}_sd']   = float(np.std(vals))
    rows.append(row)

df = pd.DataFrame(rows)
out_csv = os.path.join(out_dir, 'roi_stats.csv')
df.to_csv(out_csv, index=False, float_format='%.5f')
print(f'Saved: {out_csv}')
print(f'\n{len(rows)} ROIs extracted')
print(f'\nTop 10 ROIs by MVF_basic:')
print(df.nlargest(10, 'MVF_basic_mean')[['ROI_name', 'MVF_basic_mean', 'MVF_atlas_mean', 'g_ratio_basic_mean', 'FVF_basic_mean']].to_string(index=False))
