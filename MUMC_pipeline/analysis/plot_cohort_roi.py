"""
Cohort-level ROI figures for MIMM pipeline.

Reads cohort_roi_stats.csv (written by aggregate_cohort.py). Each dot = one
ROI mean averaged across all subjects. Error bars = between-subject SEM
(SD_between / sqrt(n_subjects)) — the correct uncertainty for comparing
cohort ROI means.

Produces figures 31–36 (cohort equivalents of per-subject 18, 20, 26, 27, 28, 29).
Saved to COHORT_FIG/.
"""

import numpy as np
import pandas as pd
import warnings; warnings.simplefilter('ignore', pd.errors.PerformanceWarning)  # ~50-row ROI tables; fragmentation is irrelevant
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import os

try:
    from cohort_paths import COHORT_ANALYSIS, COHORT_FIG, SUBJECT_DIRS
except ImportError:
    raise SystemExit('Copy cohort_paths_template.py to cohort_paths.py and fill in COHORT_DIR.')

csv = os.path.join(COHORT_ANALYSIS, 'cohort_roi_stats.csv')
if not os.path.exists(csv):
    raise SystemExit(f'cohort_roi_stats.csv not found — run aggregate_cohort.py first.\n{csv}')

df   = pd.read_csv(csv)
N    = int(df['n_subjects'].max())   # number of subjects contributing to most ROIs
names = df['ROI_name'].tolist()

os.makedirs(COHORT_FIG, exist_ok=True)
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 9})

# Helper: format the subtitle with subject count
def sub_title(base):
    return f'{base}\n(cohort N = {N} subjects, error bars = ±1 SEM between subjects)'

# ── Tract family colour map (same as plot_roi_stats.py) ──────────────────────
TRACT_FAMILIES = {
    'Corpus callosum':         ([3, 4, 5],                                    '#f4c430'),
    'Internal capsule':        ([17, 18, 19, 20, 21, 22, 29, 30],            '#5c9ee0'),
    'Corticospinal/brainstem': ([1, 2, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], '#5cb85c'),
    'Corona radiata':          ([23, 24, 25, 26, 27, 28],                    '#e8a838'),
    'Association fibres':      ([31, 32, 33, 34, 41, 42, 43, 44, 45, 46],    '#9b59b6'),
    'Limbic/cingulum':         ([6, 35, 36, 37, 38, 39, 40, 47, 48, 49, 50], '#e05c9b'),
}

def family_color(roi_idx):
    for info in TRACT_FAMILIES.values():
        if roi_idx in info[0]: return info[1]
    return '#aaaaaa'

df['family_color'] = df['ROI_index'].apply(family_color)

# Cohort chi_neg abs column (may already exist from aggregation)
if 'chi_neg_abs' not in df.columns and 'chi_neg_chisep_mean' in df.columns:
    df['chi_neg_abs'] = df['chi_neg_chisep_mean'].abs()

def _spine_style(ax):
    ax.tick_params(colors='white')
    for s in ax.spines.values(): s.set_edgecolor('#444444')
    ax.grid(color='#2a2a2a', linewidth=0.5)

# ── Figure 31: Cohort MVF Basic vs Atlas — bar chart ─────────────────────────
if {'MVF_basic_mean', 'MVF_basic_mean_sem', 'MVF_atlas_mean', 'MVF_atlas_mean_sem'}.issubset(df.columns):
    y = np.arange(len(df)); h = 0.35
    fig, ax = plt.subplots(figsize=(10, 14), facecolor='#0d0d0d')
    ax.set_facecolor('#0d0d0d')
    ax.barh(y + h/2, df['MVF_basic_mean'], h,
            xerr=df['MVF_basic_mean_sem'], color='#e05c5c', alpha=0.85,
            error_kw=dict(ecolor='white', lw=0.8, capsize=2), label='Basic (no prior)')
    ax.barh(y - h/2, df['MVF_atlas_mean'], h,
            xerr=df['MVF_atlas_mean_sem'], color='#5c9ee0', alpha=0.85,
            error_kw=dict(ecolor='white', lw=0.8, capsize=2), label='Atlas (HCP1065 prior)')
    ax.set_yticks(y); ax.set_yticklabels(names, color='white', fontsize=7.5)
    ax.set_xlabel('MVF (fraction)', color='white', fontsize=11)
    ax.set_title(sub_title('Myelin Volume Fraction per JHU WM ROI — Basic vs Atlas MIMM'),
                 color='white', fontsize=12, fontweight='bold', pad=12)
    ax.set_xlim(0, 0.55)
    ax.axvline(0, color='#444444', lw=0.8)
    ax.grid(axis='x', color='#333333', lw=0.5, linestyle='--')
    ax.legend(fontsize=10, facecolor='#1a1a1a', edgecolor='#555555',
              labelcolor='white', loc='lower right')
    _spine_style(ax)
    plt.tight_layout()
    plt.savefig(os.path.join(COHORT_FIG, '31_cohort_MVF_bar.png'),
                dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
    plt.close(); print('Saved: 31_cohort_MVF_bar.png')

# ── Figure 32: Cohort MVF Basic vs Atlas — scatter with identity ──────────────
if {'MVF_basic_mean', 'MVF_atlas_mean'}.issubset(df.columns):
    lim = (0, 0.45)
    fig, ax = plt.subplots(figsize=(7, 7), facecolor='#0d0d0d')
    ax.set_facecolor('#0d0d0d')
    ax.errorbar(df['MVF_basic_mean'], df['MVF_atlas_mean'],
                xerr=df.get('MVF_basic_mean_sem'), yerr=df.get('MVF_atlas_mean_sem'),
                fmt='o', color='#e0a050', markersize=5, alpha=0.8,
                ecolor='#888888', elinewidth=0.7, capsize=2)
    ax.plot(lim, lim, '--', color='#555555', lw=1, label='y = x')
    label_offsets = {
        'Posterior limb of internal capsule R': (8,   0, 'left'),
        'Splenium of corpus callosum':          (8,  -2, 'left'),
        'Corticospinal tract R':                (12, 14, 'left'),
        'Genu of corpus callosum':              (16, -16, 'left'),
        'Body of corpus callosum':              (20, -30, 'left'),
        'Cingulum (cingulate gyrus) R':         (-12, 16, 'right'),
    }
    for _, row in df.iterrows():
        if row['ROI_name'] in label_offsets:
            dx, dy, ha = label_offsets[row['ROI_name']]
            short = row['ROI_name'].replace(' R', '').replace(' L', '')
            ax.annotate(short, (row['MVF_basic_mean'], row['MVF_atlas_mean']),
                        textcoords='offset points', xytext=(dx, dy), ha=ha,
                        fontsize=7, color='white', alpha=0.9,
                        arrowprops=dict(arrowstyle='-', color='#888888', lw=0.5))
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel('MVF — Basic (no prior)', color='white', fontsize=11)
    ax.set_ylabel('MVF — Atlas (HCP1065 prior)', color='white', fontsize=11)
    ax.set_title(sub_title('Cohort MVF: Basic vs Atlas per JHU ROI'),
                 color='white', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9, facecolor='#1a1a1a', edgecolor='#555555', labelcolor='white')
    _spine_style(ax)
    plt.tight_layout()
    plt.savefig(os.path.join(COHORT_FIG, '32_cohort_MVF_scatter.png'),
                dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
    plt.close(); print('Saved: 32_cohort_MVF_scatter.png')

# ── Shared scatter helper ─────────────────────────────────────────────────────
def roi_scatter(df, xcol, ycol, xerr_col, yerr_col, xlabel, ylabel, title,
                fname, highlight, identity=False):
    valid = df.dropna(subset=[xcol, ycol])
    if len(valid) < 3:
        print(f'Skipped: {fname} (insufficient data)')
        return
    r, p = pearsonr(valid[xcol], valid[ycol])
    m, b = np.polyfit(valid[xcol], valid[ycol], 1)
    xline = np.linspace(valid[xcol].min(), valid[xcol].max(), 100)

    fig, ax = plt.subplots(figsize=(8, 7), facecolor='#0d0d0d')
    ax.set_facecolor('#111111')
    for family, (indices, color) in TRACT_FAMILIES.items():
        sub = valid[valid['ROI_index'].isin(indices)]
        if sub.empty: continue
        xe = sub[xerr_col] if xerr_col and xerr_col in sub.columns else None
        ye = sub[yerr_col] if yerr_col and yerr_col in sub.columns else None
        ax.errorbar(sub[xcol], sub[ycol], xerr=xe, yerr=ye,
                    fmt='o', color=color, markersize=6, alpha=0.85,
                    ecolor='#555555', elinewidth=0.6, capsize=2,
                    label=family, zorder=3)

    if identity:
        amax = max(valid[xcol].max(), valid[ycol].max()) * 1.05
        ax.plot([0, amax], [0, amax], '-', color='#555555', lw=1.2,
                label='Identity (y = x)', zorder=2)
        ax.set_xlim(0, amax); ax.set_ylim(0, amax)

    ax.plot(xline, m * xline + b, '--', color='cyan', lw=1.5,
            label=f'r = {r:.2f}, p = {p:.3f}', zorder=4)

    for _, row in valid.iterrows():
        if row['ROI_name'] in highlight:
            ax.annotate(highlight[row['ROI_name']],
                        (row[xcol], row[ycol]),
                        textcoords='offset points', xytext=(6, 3),
                        fontsize=7.5, color='white', alpha=0.9)

    ax.set_xlabel(xlabel, color='white', fontsize=11)
    ax.set_ylabel(ylabel, color='white', fontsize=11)
    ax.set_title(sub_title(title), color='white', fontsize=11, fontweight='bold')
    _spine_style(ax)
    ax.legend(fontsize=8, facecolor='#1a1a1a', edgecolor='#444444',
              labelcolor='white', loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(COHORT_FIG, fname),
                dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
    plt.close(); print(f'Saved: {fname}')

HL_MYELIN = {
    'Splenium of corpus callosum':          'Splenium CC',
    'Genu of corpus callosum':              'Genu CC',
    'Body of corpus callosum':              'Body CC',
    'Posterior limb of internal capsule R': 'Post. IC',
    'Corticospinal tract R':                'CST',
    'Cingulum (hippocampus) R':             'Cing. (hipp.)',
    'Fornix (column and body)':             'Fornix',
}
HL_IRON = {
    'Middle cerebellar peduncle':           'MCP',
    'Posterior limb of internal capsule R': 'Post. IC',
    'Corticospinal tract R':                'CST',
    'Retrolenticular internal capsule R':   'Retrolent. IC',
    'Splenium of corpus callosum':          'Splenium CC',
    'Cingulum (hippocampus) R':             'Cing. (hipp.)',
}

# ── Figure 33: Cohort MVF vs |χ⁻| ────────────────────────────────────────────
if 'chi_neg_abs' in df.columns:
    roi_scatter(df,
                xcol='chi_neg_abs',       ycol='MVF_basic_mean',
                xerr_col='chi_neg_chisep_mean_sem', yerr_col='MVF_basic_mean_sem',
                xlabel='|χ⁻| chi-sep  (ppm)', ylabel='MVF MIMM  (fraction)',
                title='Cohort: MIMM MVF vs χ-sep |χ⁻| per JHU ROI',
                fname='33_cohort_MVF_vs_chineg.png',
                highlight=HL_MYELIN)

# ── Figure 34: Cohort chi_iron vs χ⁺ ─────────────────────────────────────────
if {'chi_iron_basic_mean', 'chi_pos_chisep_mean'}.issubset(df.columns):
    roi_scatter(df,
                xcol='chi_pos_chisep_mean',    ycol='chi_iron_basic_mean',
                xerr_col='chi_pos_chisep_mean_sem', yerr_col='chi_iron_basic_mean_sem',
                xlabel='χ⁺ chi-sep  (ppm)', ylabel='χ iron MIMM  (ppm)',
                title='Cohort: MIMM χ_iron vs χ-sep χ⁺ per JHU ROI',
                fname='34_cohort_iron_chi.png',
                highlight=HL_IRON)

# ── Figure 35: Cohort |χ_myelin| vs |χ⁻| (same units — identity line) ────────
if {'chi_myelin_basic_mean', 'chi_neg_abs'}.issubset(df.columns):
    roi_scatter(df,
                xcol='chi_myelin_basic_mean',     ycol='chi_neg_abs',
                xerr_col='chi_myelin_basic_mean_sem', yerr_col='chi_neg_chisep_mean_sem',
                xlabel='|χ myelin| MIMM  (ppm)', ylabel='|χ⁻| chi-sep  (ppm)',
                title='Cohort: Diamagnetic susceptibility agreement MIMM vs χ-sep',
                fname='35_cohort_chi_myelin_vs_chineg.png',
                highlight=HL_MYELIN, identity=True)

# ── Figure 36: Cohort MVF vs FA ───────────────────────────────────────────────
# FA_mean is now in cohort_roi_stats.csv (extract_roi_stats.py maps dict).
# MVF basic is used (not atlas) — atlas uses the FA atlas for orientation,
# which would artificially inflate the correlation.
if {'FA_mean', 'MVF_basic_mean'}.issubset(df.columns):
    HL_FA = {
        'Splenium of corpus callosum':           'Splenium CC',
        'Genu of corpus callosum':               'Genu CC',
        'Posterior limb of internal capsule R':  'Post. IC',
        'Corticospinal tract R':                 'CST',
        'Cingulum (hippocampus) R':              'Cing. (hipp.)',
        'Middle cerebellar peduncle':            'MCP',
        'Fornix (column and body)':              'Fornix',
    }
    roi_scatter(df,
                xcol='FA_mean',             ycol='MVF_basic_mean',
                xerr_col='FA_mean_sem',     yerr_col='MVF_basic_mean_sem',
                xlabel='Mean FA (HCP1065 atlas)', ylabel='MVF MIMM basic  (fraction)',
                title='Cohort: MIMM MVF vs DTI FA per JHU ROI\n(FA independent of mGRE)',
                fname='36_cohort_MVF_vs_FA.png',
                highlight=HL_FA)
else:
    print('Skipped: 36_cohort_MVF_vs_FA.png (FA_mean not in cohort CSV — re-run extract_roi_stats.py)')

# ── Figure 37: Cohort MVF vs MWF (T2-GRASE) — dormant ───────────────────────
# Built in plot_mwf.py. Dormant until GRASE data arrives.
print('Skipped: 37_cohort_MVF_vs_MWF.png (see plot_mwf.py — dormant until GRASE data)')
