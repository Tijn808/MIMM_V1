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
and correlates them (Spearman) with the clinical outcomes that are available:
  edss      EDSS disability score (floored: most patients <= 2.5 in this cohort)
  sdmt      SDMT processing speed (best variance; WM-myelin sensitive) <- primary
  t25ftwt   timed 25-ft walk (s)
  hpt_dom   9-hole peg, dominant hand mean (s)
The key test is the head-to-head: do MIMM's axon measures (AVF/FVF) track an
outcome where MWF (myelin water, axon-blind) does not?

Clinical CSV: produce it from the Castor export with make_clinical_csv.py
(participant_id, edss, sdmt, t25ftwt, hpt_dom, ...). Disease duration is NOT
computed — the export carries no diagnosis date.

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

# --- clinical relationships (Spearman: small n, floored EDSS, skewed tests) ---
def scorr(x, y, metlabel, outlabel):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 6:
        print(f'  {outlabel:9s} vs {metlabel:11s}: n={int(m.sum())} too few'); return
    rho, p = stats.spearmanr(x[m], y[m])
    print(f'  {outlabel:9s} vs {metlabel:11s}: rho = {rho:+.3f}, p = {p:.3f}  (n={int(m.sum())})'
          + (' *' if p < 0.05 else ''))

if clin_csv and os.path.exists(clin_csv):
    c = pd.read_csv(clin_csv)
    c.columns = [str(k).strip() for k in c.columns]
    low = {k.lower(): k for k in c.columns}
    idc = next((low[k] for k in low if 'participant' in k or k == 'id' or 'subject' in k), None)
    c = c.rename(columns={idc: 'id'}); c['id'] = c['id'].astype(str).str.strip()
    # canonicalise outcome columns (decoder already names them, this just tolerates variants)
    for canon, keys in [('edss', ['edss']), ('sdmt', ['sdmt']),
                        ('t25ftwt', ['t25', 'walk']), ('hpt_dom', ['hpt_dom', 'peg_dom']),
                        ('hpt_nondom', ['hpt_nondom'])]:
        if canon not in c.columns:
            hit = next((low[k] for k in low if any(s in k for s in keys)), None)
            if hit:
                c = c.rename(columns={hit: canon})
    OUTCOMES = [('edss', 'EDSS'), ('sdmt', 'SDMT'), ('t25ftwt', 'T25ftWT'), ('hpt_dom', '9HPT')]
    keep = ['id'] + [o for o, _ in OUTCOMES if o in c.columns]
    img = img.merge(c[keep], on='id', how='left')

    METRICS = ['WM_MVF', 'WM_FVF', 'WM_AVF', 'WM_MWF', 'WM_iron', 'lesion_MVF', 'lesion_mL']
    print('\n=== imaging vs clinical outcomes (Spearman) ===')
    print('NB exploratory: small n, EDSS is floored (most patients <=2.5), tests are skewed.')
    for ov, olab in OUTCOMES:
        if ov not in img.columns:
            continue
        print(f'-- {olab} --')
        for met in METRICS:
            if met in img.columns:
                scorr(img[met].values, img[ov].values, met, olab)

    # --- head-to-head on the two key outcomes: MIMM axon vs the MWF reference ---
    # SDMT (processing speed) is the WM-myelin-sensitive cognitive measure with real
    # variance here; EDSS is the conventional disability scale (axon-driven but floored).
    # If AVF/FVF track an outcome where MWF does not, that is the MIMM-specific payoff.
    for ov, olab, note in [('sdmt', 'SDMT', 'cognition / WM-myelin sensitive'),
                           ('edss', 'EDSS', 'disability / axon-driven, floored')]:
        if ov not in img.columns:
            continue
        print(f'\n--- {olab} head-to-head: MIMM axon (AVF/FVF) vs MWF  [{note}] ---')
        for met, tag in [('WM_AVF', 'MIMM axon'), ('WM_FVF', 'MIMM fibre'),
                         ('WM_MVF', 'MIMM myelin'), ('WM_MWF', 'MWF (reference)')]:
            if met not in img.columns:
                continue
            x, y = img[ov].values, img[met].values
            m = np.isfinite(x) & np.isfinite(y)
            if m.sum() >= 6:
                rho, p = stats.spearmanr(x[m], y[m])
                print(f'    {tag:18s} ({met:7s}) vs {olab}: rho = {rho:+.3f}, p = {p:.3f} (n={int(m.sum())})'
                      + ('  *' if p < 0.05 else ''))

    # figure: the SDMT head-to-head — does the axonal compartment track cognition
    # where the myelin-water reference cannot?
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    panels = [('sdmt', 'WM_AVF', 'SDMT (processing speed)', 'WM-mean AVF (MIMM axon)'),
              ('sdmt', 'WM_MWF', 'SDMT (processing speed)', 'WM-mean MWF (reference)')]
    for ax, (ov, met, xl, yl) in zip(axes, panels):
        if ov in img.columns and met in img.columns:
            x, y = img[ov].values, img[met].values
            m = np.isfinite(x) & np.isfinite(y)
            if m.sum() >= 6:
                ax.scatter(x[m], y[m], s=40, c='#1f77b4', alpha=0.8)
                rho, p = stats.spearmanr(x[m], y[m])
                b, a = np.polyfit(x[m], y[m], 1); xs = np.linspace(x[m].min(), x[m].max(), 100)
                ax.plot(xs, b*xs+a, 'k--', lw=1.3)
                ax.text(0.04, 0.94, f'rho={rho:.2f}, p={p:.3f} (n={int(m.sum())})', transform=ax.transAxes,
                        va='top', bbox=dict(boxstyle='round', fc='white', ec='0.7'))
            ax.set_xlabel(xl); ax.set_ylabel(yl); ax.grid(alpha=0.25)
    fig.suptitle('Does the axonal compartment track cognition where MWF cannot?', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(ca, 'cohort_clinical.png'), dpi=150)
    print(f'\nsaved: cohort_clinical.png')
else:
    print('\nNo clinical CSV given — clinical correlations skipped.')
    print('Build one:  python3 make_clinical_csv.py <IMPROMYMS_export.csv> clinical.csv')
    print('Then run :  python3 cohort_clinical.py <results_dir> clinical.csv')
