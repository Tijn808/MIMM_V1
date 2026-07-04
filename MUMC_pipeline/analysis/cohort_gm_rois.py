#!/usr/bin/env python3
"""
Subcortical grey-matter ROIs, to match the 15-ROI set of Sisman et al. (2024):
thalamus, caudate, putamen, globus pallidus (pallidum) and hippocampus. The JHU
atlas used for the white-matter ROIs is WM-only, so these five GM structures are
taken from the FSL Harvard-Oxford subcortical atlas (same MNI space) and warped
into each subject's ME-GRE space with the *existing* mni2subj_warp.nii.gz -- the
same transform that placed the JHU labels. No new registration is computed.

Reports cohort-mean MVF, MWF and chi+ (iron) per GM structure. The globus
pallidus is the key one: it is iron-rich, so a valid myelin method should read
LOW myelin there despite the high iron (the region where 3PCF overestimates).

Output: <results>/cohort_analysis/cohort_gm_rois.csv

Usage: python3 cohort_gm_rois.py <results_dir>
Note: FSL must be available (uses applywarp). Source FSL first if needed.
"""
import sys, os, glob, subprocess
import csv
import numpy as np
import nibabel as nib
from scipy import stats

FSLDIR = os.environ.get('FSLDIR', os.path.expanduser('~/fsl'))
HO = os.path.join(FSLDIR, 'data', 'atlases', 'HarvardOxford',
                  'HarvardOxford-sub-maxprob-thr25-1mm.nii.gz')
APPLYWARP = os.path.join(FSLDIR, 'bin', 'applywarp')
if not os.path.exists(APPLYWARP):
    APPLYWARP = 'applywarp'   # hope it is on PATH

# Harvard-Oxford subcortical maxprob labels (Left, Right) for the paper's 5 GM ROIs
GM = {'Thalamus': (4, 15), 'Caudate': (5, 16), 'Putamen': (6, 17),
      'Globus pallidus': (7, 18), 'Hippocampus': (9, 19)}

if len(sys.argv) < 2:
    sys.exit('usage: cohort_gm_rois.py <results_dir>')
results_dir = sys.argv[1]
out_dir = os.path.join(results_dir, 'cohort_analysis'); os.makedirs(out_dir, exist_ok=True)
if not os.path.exists(HO):
    sys.exit(f'Harvard-Oxford atlas not found at:\n  {HO}\n'
             'Set FSLDIR or edit the HO path at the top of this script.')


def load(p):
    return np.asarray(nib.load(p).dataobj) if os.path.exists(p) else None


def first_existing(d, *rels):
    for r in rels:
        p = os.path.join(d, r)
        if os.path.exists(p):
            return p
    return None


agg = {k: [] for k in GM}   # per subject: (MVF, MWF, chi_pos)
n_done = 0
for ad in sorted(glob.glob(os.path.join(results_dir, '*', 'atlas', 'mni2subj_warp.nii.gz'))):
    d = os.path.dirname(os.path.dirname(ad)); sid = os.path.basename(d)
    ref = first_existing(d, 'atlas/mag_e1_brain.nii.gz', 'qsm/brain_mask.nii.gz')
    if ref is None:
        print(f'[skip] {sid} (no reference image)'); continue
    ho_subj = os.path.join(d, 'atlas', 'HO_sub_subj.nii.gz')
    if not os.path.exists(ho_subj):
        cmd = [APPLYWARP, '-i', HO, '-r', ref, '-w', ad, '-o', ho_subj, '--interp=nn']
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except Exception as e:
            print(f'[skip] {sid}: applywarp failed -> {e}\n  manual: {" ".join(cmd)}')
            continue
    ho = load(ho_subj)
    mvf = load(os.path.join(d, 'mimm', 'MVF_Atlas.nii.gz'))
    mwf = load(first_existing(d, 'grase/MWF.nii.gz') or 'none')
    chip = load(first_existing(d, 'chisep/chi_pos.nii.gz', 'chisep/chi_positive.nii.gz',
                               'qsm/chi_pos.nii.gz') or 'none')
    if ho is None or mvf is None:
        print(f'[skip] {sid} (missing HO warp or MVF)'); continue
    ho = ho.astype(int); mvf = mvf.astype(float)
    if mwf is not None:
        mwf = np.clip(mwf.astype(float), 0, 0.5)
    if chip is not None:
        chip = chip.astype(float)
    for name, (lL, lR) in GM.items():
        m = (ho == lL) | (ho == lR)
        if m.sum() < 20:
            agg[name].append((np.nan, np.nan, np.nan)); continue
        agg[name].append((
            float(np.nanmean(mvf[m])),
            float(np.nanmean(mwf[m])) if mwf is not None else np.nan,
            float(np.nanmean(chip[m])) if chip is not None else np.nan))
    n_done += 1
    print(f'{sid}: GM ROIs extracted')

if n_done == 0:
    sys.exit('No subjects processed (check FSL / atlas paths).')

out = os.path.join(out_dir, 'cohort_gm_rois.csv')
with open(out, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['structure', 'n_subj', 'MVF_mean', 'MVF_sem', 'MWF_mean', 'MWF_sem',
                'chi_pos_mean', 'bias_MVF_minus_MWF', 'bias_SD', 'paired_p', 'cohens_d'])
    for name in GM:
        arr = np.array(agg[name], dtype=float)
        if arr.size == 0 or not np.isfinite(arr[:, 0]).any():
            continue
        mvf, mwf, chi = arr[:, 0], arr[:, 1], arr[:, 2]
        pair = np.isfinite(mvf) & np.isfinite(mwf)
        if pair.sum() >= 6:
            diff = mvf[pair] - mwf[pair]
            bias = diff.mean(); sd = diff.std(ddof=1)
            p = stats.ttest_rel(mvf[pair], mwf[pair]).pvalue
            dval = bias / sd
        else:
            bias = sd = p = dval = np.nan
        mvf_f = mvf[np.isfinite(mvf)]; mwf_f = mwf[np.isfinite(mwf)]
        mvf_sem = mvf_f.std(ddof=1) / np.sqrt(mvf_f.size) if mvf_f.size > 1 else np.nan
        mwf_sem = mwf_f.std(ddof=1) / np.sqrt(mwf_f.size) if mwf_f.size > 1 else np.nan
        w.writerow([name, int(np.isfinite(mvf).sum()),
                    f'{np.nanmean(mvf):.4f}', f'{mvf_sem:.4f}',
                    f'{np.nanmean(mwf):.4f}', f'{mwf_sem:.4f}',
                    f'{np.nanmean(chi):.4f}', f'{bias:+.4f}', f'{sd:.4f}',
                    f'{p:.4g}', f'{dval:+.2f}'])
print(f'\nsaved: {out}  ({n_done} subjects)')
print('\nGM ROI: MIMM MVF vs MWF (paired t-test) and iron (chi+):')
for name in GM:
    arr = np.array(agg[name], dtype=float)
    if arr.size == 0 or not np.isfinite(arr[:, 0]).any():
        continue
    mvf, mwf, chi = arr[:, 0], arr[:, 1], arr[:, 2]
    pair = np.isfinite(mvf) & np.isfinite(mwf)
    if pair.sum() >= 6:
        diff = mvf[pair] - mwf[pair]
        p = stats.ttest_rel(mvf[pair], mwf[pair]).pvalue
        d = diff.mean() / diff.std(ddof=1)
        star = ' *' if p < 0.05 else ''
        print(f'  {name:16s} MVF {np.nanmean(mvf):.3f}  MWF {np.nanmean(mwf):.3f}  '
              f'bias {diff.mean():+.3f} (p={p:.2g}, d={d:+.2f}){star}  chi+ {np.nanmean(chi):.4f}')
    else:
        print(f'  {name:16s} MVF {np.nanmean(mvf):.3f}  (no paired MWF)  chi+ {np.nanmean(chi):.4f}')
