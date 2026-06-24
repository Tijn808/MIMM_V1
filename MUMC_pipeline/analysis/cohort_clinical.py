#!/usr/bin/env python3
"""
Correlate MIMM/chi-sep metrics with clinical severity. Looks for a relationship
between the imaging and disease burden — the clinically valuable question not yet
asked.

Per patient it builds:
  WM_MVF      whole-white-matter mean MVF (from roi_stats.csv)
  WM_iron     whole-white-matter mean chi_pos (chi-sep iron)
  WM_R2s      whole-white-matter mean R2*  (for R2*-vs-age, an iron proxy)
  lesion_MVF  mean MVF inside lesions      (from cohort_lesion_vs_nawm.csv)
  lesion_mL   total lesion volume          (n_lesion / 1000)
and correlates them with clinical variables: EDSS and disease duration (years).

Clinical CSV (--, any columns named loosely): participant id, edss,
and EITHER disease_duration_years OR date_diagnosis (duration computed to 2026).

Usage:
  python3 cohort_clinical.py <results_dir> [clinical.csv]
Without a clinical.csv it still runs R2*-vs-age and iron-vs-age from
demographics_mvf.csv.
"""
import sys, os, glob, datetime
import numpy as np
import pandas as pd
import nibabel as nib
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def _load(p):
    return np.asarray(nib.load(p).dataobj) if os.path.exists(p) else None

if len(sys.argv) < 2:
    sys.exit('usage: cohort_clinical.py <results_dir> [clinical.csv]')
results_dir = sys.argv[1]
clin_csv = sys.argv[2] if len(sys.argv) > 2 else None
ca = os.path.join(results_dir, 'cohort_analysis')

# --- per-subject imaging summaries from roi_stats ---
# WM_AVF = WM_FVF - WM_MVF is the axon volume fraction: the measure that, in
# theory, should track disability (EDSS) better than MWF, since EDSS is driven by
# axonal loss and MWF (myelin water) is blind to axons.
rows = []
for f in sorted(glob.glob(os.path.join(results_dir, '*', 'analysis', 'roi_stats.csv'))):
    sid = os.path.basename(os.path.dirname(os.path.dirname(f)))
    df = pd.read_csv(f)
    def col(c):
        return float(np.nanmean(df[c])) if c in df.columns else np.nan
    wm_mvf, wm_fvf = col('MVF_atlas_mean'), col('FVF_atlas_mean')
    rows.append({'id': sid, 'WM_MVF': wm_mvf, 'WM_FVF': wm_fvf,
                 'WM_AVF': wm_fvf - wm_mvf if np.isfinite(wm_fvf) and np.isfinite(wm_mvf) else np.nan,
                 'WM_iron': col('chi_pos_chisep_mean'), 'WM_R2s': col('R2s_atlas_mean')})
img = pd.DataFrame(rows)

# whole-WM mean MWF from the registered T2-GRASE map (for the MWF-vs-EDSS head-to-head)
mwf_rows = []
for mp in sorted(glob.glob(os.path.join(results_dir, '*', 'grase', 'MWF.nii.gz'))):
    d = os.path.dirname(os.path.dirname(mp)); sid = os.path.basename(d)
    mwf = _load(mp); brain = _load(os.path.join(d, 'qsm', 'brain_mask.nii.gz'))
    fa = _load(os.path.join(d, 'atlas', 'FA_atlas.nii.gz'))
    if any(v is None for v in (mwf, brain, fa)):
        continue
    mwf = np.clip(mwf.astype(float), 0, 0.5)
    wm = (brain > 0) & (fa > 0.20) & np.isfinite(mwf) & (mwf > 0)
    mwf_rows.append({'id': sid, 'WM_MWF': float(mwf[wm].mean()) if wm.sum() else np.nan})
if mwf_rows:
    img = img.merge(pd.DataFrame(mwf_rows), on='id', how='left')

# lesion MVF + volume
les_csv = os.path.join(ca, 'cohort_lesion_vs_nawm.csv')
if os.path.exists(les_csv):
    l = pd.read_csv(les_csv)
    lcol = 'MVF (Atlas)_lesion' if 'MVF (Atlas)_lesion' in l.columns else None
    img = img.merge(pd.DataFrame({'id': l['subject'],
                                  'lesion_MVF': l[lcol] if lcol else np.nan,
                                  'lesion_mL': l['n_lesion'] / 1000.0}), on='id', how='left')

# age from demographics_mvf.csv (for R2*/iron-vs-age, always available)
demo_csv = os.path.join(ca, 'demographics_mvf.csv')
if os.path.exists(demo_csv):
    dm = pd.read_csv(demo_csv)[['id', 'age']]
    img = img.merge(dm, on='id', how='left')


def corr(x, y, xlabel, ylabel):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4:
        print(f'  {ylabel:11s} vs {xlabel:18s}: n={m.sum()} too few'); return
    r, p = stats.pearsonr(x[m], y[m])
    star = ' *' if p < 0.05 else ''
    print(f'  {ylabel:11s} vs {xlabel:18s}: r = {r:+.3f}, p = {p:.3f}  (n={m.sum()}){star}')

# --- age relationships (no clinical file needed) ---
print('=== imaging vs age (iron accumulation proxy) ===')
if 'age' in img.columns:
    corr(img['age'].values, img['WM_R2s'].values, 'age', 'WM R2*')
    corr(img['age'].values, img['WM_iron'].values, 'age', 'WM iron')

# --- clinical relationships ---
if clin_csv and os.path.exists(clin_csv):
    c = pd.read_csv(clin_csv)
    cols = {k.lower().strip(): k for k in c.columns}
    idc = next((v for k, v in cols.items() if 'participant' in k or k == 'id' or 'subject' in k), None)
    edssc = next((v for k, v in cols.items() if 'edss' in k), None)
    durc = next((v for k, v in cols.items() if 'duration' in k), None)
    diagc = next((v for k, v in cols.items() if 'diagnos' in k and 'date' in k), None)
    c = c.rename(columns={idc: 'id'}); c['id'] = c['id'].astype(str).str.strip()
    if edssc: c = c.rename(columns={edssc: 'edss'})
    if durc:
        c = c.rename(columns={durc: 'dur_years'})
    elif diagc:
        def yrs(s):
            try:
                return (datetime.date(2026, 6, 24) - pd.to_datetime(s).date()).days / 365.25
            except Exception:
                return np.nan
        c['dur_years'] = c[diagc].map(yrs)
    img = img.merge(c[[col for col in ['id', 'edss', 'dur_years'] if col in c.columns]], on='id', how='left')

    print('\n=== imaging vs clinical severity ===')
    METRICS = ['WM_MVF', 'WM_FVF', 'WM_AVF', 'WM_MWF', 'WM_iron', 'lesion_MVF', 'lesion_mL']
    for clinvar, lab in [('edss', 'EDSS'), ('dur_years', 'disease duration')]:
        if clinvar not in img.columns:
            print(f'  ({lab} not found in clinical file)'); continue
        for met in METRICS:
            if met in img.columns:
                corr(img[clinvar].values, img[met].values, lab, met)

    # --- head-to-head: do MIMM's axon measures track EDSS better than MWF? ---
    if 'edss' in img.columns:
        print('\n--- EDSS head-to-head: MIMM axon (AVF/FVF) vs MWF ---')
        print('  EDSS is axon-driven; MWF is blind to axons. If AVF/FVF track EDSS')
        print('  and MWF does not, that is the MIMM-specific advantage.')
        def r_of(met):
            if met not in img.columns:
                return None
            x, y = img['edss'].values, img[met].values
            m = np.isfinite(x) & np.isfinite(y)
            return (stats.pearsonr(x[m], y[m]) + (int(m.sum()),)) if m.sum() >= 4 else None
        for met, tag in [('WM_AVF', 'MIMM axon'), ('WM_FVF', 'MIMM fibre'),
                         ('WM_MVF', 'MIMM myelin'), ('WM_MWF', 'MWF (reference)')]:
            res = r_of(met)
            if res:
                r, p, nn = res
                print(f'    {tag:18s} ({met:7s}) vs EDSS: r = {r:+.3f}, p = {p:.3f} (n={nn})'
                      + ('  *' if p < 0.05 else ''))

    # figure: the two most clinically interesting scatters
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, (clinvar, met, xl, yl) in zip(axes, [('edss', 'WM_MVF', 'EDSS', 'WM-mean MVF'),
                                                 ('dur_years', 'lesion_mL', 'disease duration (y)', 'lesion volume (mL)')]):
        if clinvar in img.columns and met in img.columns:
            x, y = img[clinvar].values, img[met].values
            m = np.isfinite(x) & np.isfinite(y)
            if m.sum() >= 4:
                ax.scatter(x[m], y[m], s=40, c='#1f77b4', alpha=0.8)
                r, p = stats.pearsonr(x[m], y[m])
                b, a = np.polyfit(x[m], y[m], 1); xs = np.linspace(x[m].min(), x[m].max(), 100)
                ax.plot(xs, b*xs+a, 'k--', lw=1.3)
                ax.text(0.04, 0.94, f'r={r:.2f}, p={p:.3f} (n={m.sum()})', transform=ax.transAxes,
                        va='top', bbox=dict(boxstyle='round', fc='white', ec='0.7'))
            ax.set_xlabel(xl); ax.set_ylabel(yl); ax.grid(alpha=0.25)
    fig.suptitle('MIMM imaging vs clinical severity', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(ca, 'cohort_clinical.png'), dpi=150)
    print(f'\nsaved: cohort_clinical.png')
else:
    print('\nNo clinical CSV given — EDSS / disease-duration correlations skipped.')
    print('Provide one:  python3 cohort_clinical.py <results_dir> clinical.csv')
