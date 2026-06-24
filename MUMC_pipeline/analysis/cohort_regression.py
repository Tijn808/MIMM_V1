#!/usr/bin/env python3
"""
Multiple regression of whole-WM-mean MVF on demographics (age, sex) and, in a
second model, lesion burden as a covariate. Answers two supervisor questions:
  (1) age/sex correlation+regression (not just univariate)
  (3) does adding lesion load as a covariate change the picture

Model A:  MVF ~ age + sex                  (all subjects with demographics)
Model B:  MVF ~ age + sex + lesion_volume  (subjects with known lesion volume:
          HCs = 0 mL; patients with a lesion mask = measured; patients without a
          mask are dropped, since their lesion load is unknown)

Inputs:
  <results>/cohort_analysis/demographics_mvf.csv   (id, MVF_atlas, age, sex, group)
  <results>/cohort_analysis/cohort_lesion_vs_nawm.csv  (subject, n_lesion)  [optional]

Usage:  python3 cohort_regression.py <results_dir>
"""
import sys, os
import numpy as np
import pandas as pd
from scipy import stats

if len(sys.argv) < 2:
    sys.exit('usage: cohort_regression.py <results_dir>')
results_dir = sys.argv[1]
ca = os.path.join(results_dir, 'cohort_analysis')

demo = pd.read_csv(os.path.join(ca, 'demographics_mvf.csv'))
demo['id'] = demo['id'].astype(str).str.strip()
demo['sex_male'] = (demo['sex'].astype(str).str.strip().str.lower().str[0] == 'm').astype(float)

# Lesion volume (mL) per subject. 1 mm iso -> n_voxels / 1000 mL.
les_csv = os.path.join(ca, 'cohort_lesion_vs_nawm.csv')
les_vol = {}
if os.path.exists(les_csv):
    ldf = pd.read_csv(les_csv)
    for _, r in ldf.iterrows():
        les_vol[str(r['subject']).strip()] = float(r['n_lesion']) / 1000.0
    print(f'Lesion volumes loaded for {len(les_vol)} patients with masks.')
else:
    print('No cohort_lesion_vs_nawm.csv — run cohort_lesion.py first for Model B.')

HC = {'IMPROMYMS_001', 'IMPROMYMS_005', 'IMPROMYMS_006', 'IMPROMYMS_020'}
def lesion_for(sid):
    if sid in les_vol:
        return les_vol[sid]
    if sid in HC:
        return 0.0          # controls have no MS lesions
    return np.nan           # patient without a mask: unknown


def ols(y, X, names):
    """Plain OLS with t-tests. X already includes an intercept column."""
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    n, p = X.shape
    dof = n - p
    sigma2 = (resid @ resid) / dof
    XtX_inv = np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(XtX_inv) * sigma2)
    t = beta / se
    pvals = 2 * stats.t.sf(np.abs(t), dof)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - (resid @ resid) / ss_tot
    r2_adj = 1 - (1 - r2) * (n - 1) / dof
    print(f'  n = {n}, R2 = {r2:.3f}, adj-R2 = {r2_adj:.3f}')
    for nm, b, s, tt, pp in zip(names, beta, se, t, pvals):
        star = '*' if pp < 0.05 else ' '
        print(f'    {nm:14s} beta = {b:+.5f}  SE = {s:.5f}  t = {tt:+.2f}  p = {pp:.3f} {star}')
    return beta, pvals


print('\n=== Model A:  MVF ~ age + sex ===')
dA = demo.dropna(subset=['MVF_atlas', 'age', 'sex_male'])
yA = dA['MVF_atlas'].values
XA = np.column_stack([np.ones(len(dA)), dA['age'].values, dA['sex_male'].values])
ols(yA, XA, ['intercept', 'age', 'sex(male)'])

print('\n=== Model B:  MVF ~ age + sex + lesion_volume ===')
demo['lesion_mL'] = demo['id'].map(lesion_for)
dB = demo.dropna(subset=['MVF_atlas', 'age', 'sex_male', 'lesion_mL'])
if len(dB) >= 5:
    yB = dB['MVF_atlas'].values
    XB = np.column_stack([np.ones(len(dB)), dB['age'].values,
                          dB['sex_male'].values, dB['lesion_mL'].values])
    ols(yB, XB, ['intercept', 'age', 'sex(male)', 'lesion(mL)'])
    print(f'  (Model B drops patients with unknown lesion status; '
          f'{len(dB)} of {len(demo)} subjects used.)')
else:
    print('  Too few subjects with known lesion volume for Model B.')
