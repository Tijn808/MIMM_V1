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
from mimm_style import apply_style, C   # shared deck palette (style only, numbers unchanged)
apply_style()


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
    # chi-separation diamagnetic (myelin) channel; try the likely column names
    wm_chineg = next((col(c) for c in ('chi_neg_chisep_mean', 'chi_neg_mean',
                                       'abs_chi_neg_chisep_mean') if c in df.columns), np.nan)
    # MIMM g-ratio (axon/fibre geometry) -- a MIMM-unique output
    wm_gratio = next((col(c) for c in ('g_ratio_atlas_mean', 'g_ratio_Atlas_mean',
                                       'gratio_atlas_mean') if c in df.columns), np.nan)
    rows.append({'id': sid, 'WM_MVF': wm_mvf, 'WM_FVF': wm_fvf,
                 'WM_AVF': wm_fvf - wm_mvf if np.isfinite(wm_fvf) and np.isfinite(wm_mvf) else np.nan,
                 'WM_gratio': wm_gratio, 'WM_chineg': wm_chineg,
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


def williams_test(r1, r2, r12, n):
    """Williams's t (1959): is corr(y,x1) [=r1] significantly different from
    corr(y,x2) [=r2] when the two correlations share the outcome y? r12 is the
    correlation between the two predictors x1 and x2; these are DEPENDENT
    (overlapping) correlations measured on the same n subjects. Returns (t, p).
    Applied to Spearman correlations it is an approximation (exact for Pearson)."""
    if not all(np.isfinite([r1, r2, r12])) or n < 5:
        return np.nan, np.nan
    detR = 1 - r1**2 - r2**2 - r12**2 + 2 * r1 * r2 * r12
    num = (r1 - r2) * np.sqrt((n - 1) * (1 + r12))
    den = np.sqrt(2 * detR * (n - 1) / (n - 3) + ((r1 + r2)**2 / 4) * (1 - r12)**3)
    if den == 0:
        return np.nan, np.nan
    t = num / den
    return t, 2 * stats.t.sf(abs(t), n - 3)


def compare_rhos(ov, m1, m2):
    """Williams test comparing the |Spearman rho| of m1 vs m2 against outcome ov,
    on the subjects with all three present. Returns a result string or None."""
    if any(k not in img.columns for k in (ov, m1, m2)):
        return None
    x, a, b = img[ov].values, img[m1].values, img[m2].values
    mask = np.isfinite(x) & np.isfinite(a) & np.isfinite(b)
    n = int(mask.sum())
    if n < 8:
        return None
    r1 = stats.spearmanr(x[mask], a[mask]).correlation   # outcome vs m1
    r2 = stats.spearmanr(x[mask], b[mask]).correlation   # outcome vs m2
    r12 = stats.spearmanr(a[mask], b[mask]).correlation  # m1 vs m2
    _, p = williams_test(r1, r2, r12, n)
    sig = ' *' if (np.isfinite(p) and p < 0.05) else ''
    return (f'    {m1:8s} (rho={r1:+.2f}) vs {m2:8s} (rho={r2:+.2f}): '
            f'Williams p = {p:.3f} (n={n}){sig}')

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

    METRICS = ['WM_MVF', 'WM_FVF', 'WM_AVF', 'WM_gratio', 'WM_chineg', 'WM_MWF', 'WM_iron', 'lesion_MVF', 'lesion_mL']
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
    for ov, olab, note in [('hpt_dom', '9HPT', 'dexterity / strongest signal'),
                           ('sdmt', 'SDMT', 'cognition / WM-myelin sensitive'),
                           ('edss', 'EDSS', 'disability / axon-driven, floored')]:
        if ov not in img.columns:
            continue
        print(f'\n--- {olab} head-to-head: every technique\'s myelin metric  [{note}] ---')
        for met, tag in [('WM_AVF', 'MIMM axon'), ('WM_FVF', 'MIMM fibre'),
                         ('WM_MVF', 'MIMM myelin'), ('WM_chineg', 'chi-sep myelin'),
                         ('WM_MWF', 'MWF (reference)')]:
            if met not in img.columns:
                continue
            x, y = img[ov].values, img[met].values
            m = np.isfinite(x) & np.isfinite(y)
            if m.sum() >= 6:
                rho, p = stats.spearmanr(x[m], y[m])
                print(f'    {tag:18s} ({met:8s}) vs {olab}: rho = {rho:+.3f}, p = {p:.3f} (n={int(m.sum())})'
                      + ('  *' if p < 0.05 else ''))
        # are the MIMM-fibre rho and the comparator rhos significantly DIFFERENT?
        # (dependent correlations sharing the outcome -> Williams's test)
        print(f'    -- rho-vs-rho (is MIMM fibre a stronger correlate?) --')
        for comp in ('WM_MVF', 'WM_chineg', 'WM_MWF'):
            line = compare_rhos(ov, 'WM_FVF', comp)
            if line:
                print(line)

    # --- partial correlation: does the axon compartment track 9HPT once myelin
    #     is held constant? (Spearman partial: rank-transform, then partial Pearson) ---
    def partial_spearman(ov, target, control):
        cols = [ov, target, control]
        if any(c not in img.columns for c in cols):
            return None
        sub = img[cols].apply(pd.to_numeric, errors='coerce').dropna()
        if len(sub) < 8:
            return None
        rk = sub.rank()
        r = rk.corr().values  # 0=ov, 1=target, 2=control
        r_ot, r_oc, r_tc = r[0, 1], r[0, 2], r[1, 2]
        denom = np.sqrt((1 - r_oc**2) * (1 - r_tc**2))
        if denom == 0:
            return None
        pr = (r_ot - r_oc * r_tc) / denom
        n = len(sub); dfree = n - 3
        t = pr * np.sqrt(dfree / (1 - pr**2)) if abs(pr) < 1 else np.inf
        return pr, 2 * stats.t.sf(abs(t), dfree), n

    print('\n=== partial correlation: axon vs 9HPT, controlling for myelin ===')
    for tgt, lab in [('WM_AVF', 'AVF'), ('WM_FVF', 'FVF')]:
        pc = partial_spearman('hpt_dom', tgt, 'WM_MVF')
        if pc:
            print(f'  {lab} vs 9HPT | MVF held constant: partial rho = {pc[0]:+.3f}, '
                  f'p = {pc[1]:.3f} (n={pc[2]})' + ('  *' if pc[1] < 0.05 else ''))

    # ---- comprehensive clinical figure ----
    # (A) headline regplot FVF vs 9HPT with bootstrap CI band
    # (B,C) forest plots: every metric's Spearman rho vs 9HPT / SDMT, with
    #       bootstrap 95% CIs, oriented so positive = better function.
    METR = [('WM_FVF', 'FVF (fibre)', True), ('WM_AVF', 'AVF (axon)', True),
            ('WM_gratio', 'g-ratio', True),
            ('WM_MVF', 'MVF (myelin)', False), ('WM_chineg', 'chi-  (chi-sep)', False),
            ('WM_MWF', 'MWF (reference)', False)]
    ORIENT = {'hpt_dom': -1, 'sdmt': +1}   # 9HPT: lower time = better -> flip sign

    def boot_ci(x, y, n_boot=2000):
        rng = np.random.default_rng(0); n = len(x); out = []
        for _ in range(n_boot):
            i = rng.integers(0, n, n)
            if np.std(x[i]) > 0 and np.std(y[i]) > 0:
                out.append(stats.spearmanr(x[i], y[i]).correlation)
        out = np.array(out)
        return np.nanpercentile(out, 2.5), np.nanpercentile(out, 97.5)

    def metric_stat(ov, met):
        x, y = img[ov].values, img[met].values
        m = np.isfinite(x) & np.isfinite(y)
        if m.sum() < 6:
            return None
        rho, p = stats.spearmanr(x[m], y[m])
        lo, hi = boot_ci(x[m], y[m])
        return rho, p, lo, hi

    def forest(ax, ov, title):
        sgn = ORIENT[ov]; data = []
        for met, lab, is_c in METR:
            s = metric_stat(ov, met)
            if s is None:
                continue
            rho, p, lo, hi = s
            orho = sgn * rho; olo, ohi = sorted([sgn * lo, sgn * hi])
            data.append((lab, is_c, orho, olo, ohi, p))
        data.sort(key=lambda d: d[2])                  # weakest at bottom, strongest at top
        for yi, (lab, is_c, orho, olo, ohi, p) in enumerate(data):
            col = C['MIMM'] if is_c else C['reference']
            ax.errorbar(orho, yi, xerr=[[max(0, orho - olo)], [max(0, ohi - orho)]],
                        fmt='o', color=col, ecolor=col, capsize=3, ms=8, lw=1.8, zorder=3)
            if p < 0.05:
                ax.text(ohi + 0.03, yi, '*', va='center', color=col, fontsize=14)
        ax.axvline(0, color='k', lw=0.9, zorder=1)
        ax.set_yticks(range(len(data))); ax.set_yticklabels([d[0] for d in data], fontsize=9)
        ax.set_xlim(-0.55, 0.95); ax.set_xlabel('Spearman rho  (positive = better function)')
        ax.set_title(title, fontsize=11); ax.grid(alpha=0.2, axis='x')

    fig, axd = plt.subplot_mosaic([['A', 'B'], ['A', 'C']], figsize=(13, 6.5))
    # (A) headline scatter with bootstrap regression band
    ax = axd['A']
    x, y = img['hpt_dom'].values, img['WM_FVF'].values
    m = np.isfinite(x) & np.isfinite(y); x, y = x[m], y[m]
    rho, p = stats.spearmanr(x, y)
    xs = np.linspace(x.min(), x.max(), 100)
    rng = np.random.default_rng(0); preds = []
    for _ in range(1000):
        i = rng.integers(0, len(x), len(x))
        bb, aa = np.polyfit(x[i], y[i], 1); preds.append(bb * xs + aa)
    preds = np.array(preds)
    ax.fill_between(xs, np.percentile(preds, 2.5, 0), np.percentile(preds, 97.5, 0),
                    color=C['MIMM'], alpha=0.18, zorder=1)
    b, a = np.polyfit(x, y, 1); ax.plot(xs, b * xs + a, 'k--', lw=1.5, zorder=2)
    ax.scatter(x, y, s=44, c=C['MIMM'], alpha=0.9, zorder=3, edgecolor='white', lw=0.5)
    ax.text(0.04, 0.06, f'FVF vs 9HPT\nrho = {rho:.2f}, p = {p:.3f} (n={len(x)})',
            transform=ax.transAxes, va='bottom',
            bbox=dict(boxstyle='round', fc='white', ec='0.7'))
    ax.set_xlabel('9-hole peg test, dominant hand (s)')
    ax.set_ylabel('WM-mean FVF (MIMM fibre)'); ax.grid(alpha=0.25)
    forest(axd['B'], 'hpt_dom', '9-hole peg (dexterity)')
    forest(axd['C'], 'sdmt', 'SDMT (processing speed)')
    from matplotlib.patches import Patch
    axd['C'].legend(handles=[Patch(color=C['MIMM'], label='MIMM-specific (FVF, AVF, g-ratio)'),
                             Patch(color=C['reference'], label='single-number myelin (MVF, chi-, MWF)')],
                    loc='lower right', fontsize=8, framealpha=0.9)
    fig.suptitle('Clinical correlates: the fibre/axon compartment vs the single-number measures '
                 '(bootstrap 95% CI)', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(os.path.join(ca, 'cohort_clinical.png'), dpi=150)
    print('\nsaved: cohort_clinical.png  (regplot + bootstrap-CI forest plots)')

    # --- 9HPT-only variant for the deck (panels A + B only, no SDMT panel) ---
    fig2, axd2 = plt.subplot_mosaic([['A', 'B']], figsize=(11, 5))
    ax2 = axd2['A']
    ax2.fill_between(xs, np.percentile(preds, 2.5, 0), np.percentile(preds, 97.5, 0),
                     color=C['MIMM'], alpha=0.18, zorder=1)
    ax2.plot(xs, b * xs + a, 'k--', lw=1.5, zorder=2)
    ax2.scatter(x, y, s=44, c=C['MIMM'], alpha=0.9, zorder=3, edgecolor='white', lw=0.5)
    ax2.text(0.04, 0.06, f'FVF vs 9HPT\nrho = {rho:.2f}, p = {p:.3f} (n={len(x)})',
             transform=ax2.transAxes, va='bottom',
             bbox=dict(boxstyle='round', fc='white', ec='0.7'))
    ax2.set_xlabel('9-hole peg test, dominant hand (s)')
    ax2.set_ylabel('WM-mean FVF (MIMM fibre)'); ax2.grid(alpha=0.25)
    forest(axd2['B'], 'hpt_dom', '9-hole peg (dexterity)')
    axd2['B'].legend(handles=[Patch(color=C['MIMM'], label='MIMM-specific (FVF, AVF, g-ratio)'),
                              Patch(color=C['reference'], label='single-number myelin (MVF, chi-, MWF)')],
                     loc='lower right', fontsize=8, framealpha=0.9)
    fig2.suptitle('Clinical correlate: the fibre/axon compartment vs the single-number measures, '
                  '9-hole peg test (bootstrap 95% CI)', fontsize=12)
    fig2.tight_layout(rect=[0, 0, 1, 0.94])
    fig2.savefig(os.path.join(ca, 'cohort_clinical_9hpt.png'), dpi=150)
    print('saved: cohort_clinical_9hpt.png  (9HPT-only: regplot + forest plot, no SDMT panel)')
else:
    print('\nNo clinical CSV given — clinical correlations skipped.')
    print('Build one:  python3 make_clinical_csv.py <IMPROMYMS_export.csv> clinical.csv')
    print('Then run :  python3 cohort_clinical.py <results_dir> clinical.csv')
