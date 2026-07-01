#!/usr/bin/env python3
"""
check_r2prime.py  --  Does swapping the R2* proxy for a true R2' (from T2-GRASE)
change the Figure 1 per-subject cross-method correlations?

We do NOT re-run chi-separation. Using the identities
    chi+ = R2' / alpha ,   chi- = QSM - chi+
the switch from the R2* proxy to a true R2' = R2* - R2 shifts the whole-WM means by
exactly R2/alpha:
    chi-_true = chi-_proxy + R2bar/alpha
    chi+_true = chi+_proxy - R2bar/alpha        (alpha = 124 s^-1/ppm, Shin 2021)

So per subject we only need R2bar = whole-WM mean R2 from the T2-GRASE decay.

Per subject:
  1. monoexponential R2 fit from nifti/grase.nii (first NFIT echoes),
  2. register R2 to ME-GRE space with the EXISTING grase/grase2mgre.mat (flirt applyxfm),
  3. whole-WM mean R2 (brain_mask & FA_atlas > 0.20),
  4. correct the proxy chi means and recompute the between-subject Pearson r.

Proxy whole-WM means (MVF, MIMM iron, chi-, chi+) are read from analysis/roi_stats.csv
exactly as fig_crossmethod does (voxel-weighted over JHU ROIs), so r_proxy reproduces
the paper's ~0.90 and r_true is the apples-to-apples comparison.

Usage:  python3 check_r2prime.py <results_dir>
"""
import sys, os, glob, subprocess, tempfile
import numpy as np, pandas as pd, nibabel as nib
from scipy import stats

ALPHA = 124.0                       # s^-1 / ppm  (Shin et al. 2021)
DTE   = 0.010                       # echo spacing (s)  -- EDIT if not 10 ms
NECHO = 32                          # echoes in grase.nii
NFIT  = 12                          # echoes used for the monoexp fit (10-120 ms; avoids noise floor)
TE_ALL = np.arange(1, NECHO + 1) * DTE
TE_FIT = TE_ALL[:NFIT]
FSL   = os.environ.get('FSLDIR', '')
FLIRT = os.path.join(FSL, 'bin', 'flirt') if FSL else 'flirt'


def load(p):
    return np.asarray(nib.load(p).dataobj).astype(float) if os.path.exists(p) else None


def wmean(v, w):
    v = np.asarray(v, float); w = np.asarray(w, float)
    m = np.isfinite(v) & np.isfinite(w) & (w > 0)
    return float(np.sum(v[m] * w[m]) / np.sum(w[m])) if m.any() else np.nan


def fit_r2(grase4d):
    """monoexponential log-linear R2 per voxel over the first NFIT echoes."""
    S = grase4d[..., :NFIT].reshape(-1, NFIT).T          # (NFIT, nvox)
    R2 = np.zeros(S.shape[1])
    good = np.all(S > 0, axis=0)
    if good.any():
        A = np.vstack([TE_FIT, np.ones_like(TE_FIT)]).T
        coef, *_ = np.linalg.lstsq(A, np.log(S[:, good]), rcond=None)
        R2[good] = -coef[0]
    return R2.reshape(grase4d.shape[:-1])


if len(sys.argv) < 2:
    sys.exit('usage: check_r2prime.py <results_dir>')
results = sys.argv[1]
rows = []
for rs in sorted(glob.glob(os.path.join(results, '*', 'analysis', 'roi_stats.csv'))):
    d = os.path.dirname(os.path.dirname(rs)); sid = os.path.basename(d)
    df = pd.read_csv(rs)
    if 'n_voxels' not in df.columns:
        continue
    w = df['n_voxels'].values

    def m(c):
        return wmean(df[c].values, w) if c in df.columns else np.nan

    mvf, iron = m('MVF_atlas_mean'), m('chi_iron_atlas_mean')
    chineg, chipos = m('chi_neg_chisep_mean'), m('chi_pos_chisep_mean')

    grase = os.path.join(d, 'nifti', 'grase.nii')
    xfm   = os.path.join(d, 'grase', 'grase2mgre.mat')
    ref   = os.path.join(d, 'qsm', 'mag_e1.nii.gz')
    brain = load(os.path.join(d, 'qsm', 'brain_mask.nii.gz'))
    fa    = load(os.path.join(d, 'atlas', 'FA_atlas.nii.gz'))
    if not all(os.path.exists(p) for p in (grase, xfm, ref)) or brain is None or fa is None:
        print(f'[skip] {sid} (missing grase/xfm/ref/mask/fa)'); continue

    gimg = nib.load(grase)
    g = np.asarray(gimg.dataobj).astype(float)
    if g.ndim != 4 or g.shape[-1] != NECHO:
        print(f'[skip] {sid} (grase 4th dim = {g.shape[-1] if g.ndim == 4 else "n/a"}, expected {NECHO})')
        continue
    r2 = fit_r2(g)

    with tempfile.TemporaryDirectory() as tmp:
        r2_native = os.path.join(tmp, 'r2.nii.gz')
        nib.save(nib.Nifti1Image(r2.astype(np.float32), gimg.affine), r2_native)
        r2_reg = os.path.join(tmp, 'r2_mgre.nii.gz')
        cmd = [FLIRT, '-in', r2_native, '-ref', ref, '-applyxfm', '-init', xfm,
               '-out', r2_reg, '-interp', 'trilinear']
        if subprocess.run(cmd, capture_output=True).returncode != 0:
            print(f'[skip] {sid} (flirt failed)'); continue
        r2m = load(r2_reg)

    wm = (brain > 0) & (fa > 0.20) & np.isfinite(r2m) & (r2m > 0) & (r2m < 100)
    r2_wm = float(np.mean(r2m[wm])) if wm.sum() else np.nan
    corr = r2_wm / ALPHA
    rows.append(dict(sid=sid, mvf=mvf, iron=iron, chineg=chineg, chipos=chipos,
                     chineg_true=chineg + corr, chipos_true=chipos - corr, r2_wm=r2_wm))
    print(f'{sid}: R2_WM = {r2_wm:5.1f} s^-1   correction = {corr:+.3f} ppm')

df = pd.DataFrame(rows)
if len(df) < 4:
    sys.exit(f'\nOnly {len(df)} usable subjects.')


def r(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    mm = np.isfinite(x) & np.isfinite(y)
    return stats.pearsonr(x[mm], y[mm])[0], int(mm.sum())


print(f'\nn = {len(df)} subjects   mean R2_WM = {df.r2_wm.mean():.1f} +- {df.r2_wm.std():.1f} s^-1'
      f'  (CV = {100 * df.r2_wm.std() / df.r2_wm.mean():.1f}%)')
rp, n = r(np.abs(df.chineg), df.mvf); rt, _ = r(np.abs(df.chineg_true), df.mvf)
print(f'MYELIN  |chi-| vs MIMM MVF :  r_proxy = {rp:.3f}   r_true = {rt:.3f}   (n={n})')
rp, n = r(df.chipos, df.iron);        rt, _ = r(df.chipos_true, df.iron)
print(f'IRON    chi+   vs MIMM iron:  r_proxy = {rp:.3f}   r_true = {rt:.3f}   (n={n})')
print('\nIf r_true is within ~0.02 of r_proxy, the R2* proxy does not change Figure 1 '
      'and you keep it with a verified sentence.')
