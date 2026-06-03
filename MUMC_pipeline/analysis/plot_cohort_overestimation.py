"""
Cohort-level overestimation analysis: MIMM Basic vs Atlas MVF.

Cohort equivalents of per-subject figures 23–25. Each point = one ROI,
values averaged across subjects. Error bars = between-subject SEM.

Requires cohort_roi_stats.csv (run aggregate_cohort.py first).
FA_mean and theta_mean must be in the CSV (re-run extract_roi_stats.py
if the CSV was generated before these were added to the maps dict).

Produces:
  42_cohort_overestimation_ranked.png  — all ROIs ranked by cohort-mean overestimation
  43_cohort_overestimation_vs_theta.png — overestimation vs fibre angle
  44_cohort_overestimation_vs_FA.png   — overestimation vs FA
"""

import numpy as np
import pandas as pd
import warnings; warnings.simplefilter('ignore', pd.errors.PerformanceWarning)  # ~50-row ROI tables; fragmentation is irrelevant
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy.stats import pearsonr
import os

try:
    from cohort_paths import COHORT_ANALYSIS, COHORT_FIG, SUBJECT_DIRS
except ImportError:
    raise SystemExit('Copy cohort_paths_template.py to cohort_paths.py and fill in COHORT_DIR.')

csv = os.path.join(COHORT_ANALYSIS, 'cohort_roi_stats.csv')
if not os.path.exists(csv):
    raise SystemExit(f'cohort_roi_stats.csv not found — run aggregate_cohort.py first.\n{csv}')

df = pd.read_csv(csv)
N  = int(df['n_subjects'].max())
os.makedirs(COHORT_FIG, exist_ok=True)

# ── Check required columns ────────────────────────────────────────────────────
missing = [c for c in ['MVF_basic_mean', 'MVF_atlas_mean', 'FA_mean', 'theta_mean']
           if c not in df.columns]
if missing:
    raise SystemExit(
        f'Missing columns: {missing}\n'
        'Re-run extract_roi_stats.py (FA and theta now in maps dict) '
        'then re-run aggregate_cohort.py.')

# ── Compute cohort-mean overestimation ────────────────────────────────────────
df['overest_abs'] = df['MVF_basic_mean'] - df['MVF_atlas_mean']
df['overest_rel'] = (df['overest_abs'] / df['MVF_atlas_mean']) * 100

# Propagated SEM for the difference (independent estimates, add in quadrature)
if {'MVF_basic_mean_sem', 'MVF_atlas_mean_sem'}.issubset(df.columns):
    df['overest_sem'] = np.sqrt(df['MVF_basic_mean_sem']**2 + df['MVF_atlas_mean_sem']**2)
else:
    df['overest_sem'] = np.nan

df_sorted = df.sort_values('overest_abs', ascending=True).reset_index(drop=True)

def _spine(ax):
    ax.tick_params(colors='white')
    for s in ax.spines.values(): s.set_edgecolor('#444444')
    ax.grid(axis='x', color='#333333', lw=0.5, ls='--')

# ── Figure 42: Cohort ranked overestimation bar chart ─────────────────────────
fig, ax = plt.subplots(figsize=(10, 14), facecolor='#0d0d0d')
ax.set_facecolor('#0d0d0d')

colors = ['#e05c5c' if v > 0 else '#5c9ee0' for v in df_sorted['overest_abs']]
ax.barh(range(len(df_sorted)), df_sorted['overest_abs'], color=colors, alpha=0.85,
        xerr=df_sorted['overest_sem'] if df_sorted['overest_sem'].notna().any() else None,
        error_kw=dict(ecolor='#aaaaaa', lw=0.7, capsize=2))

for i, (abs_v, rel_v) in enumerate(zip(df_sorted['overest_abs'], df_sorted['overest_rel'])):
    x  = abs_v + 0.002 if abs_v >= 0 else abs_v - 0.002
    ha = 'left' if abs_v >= 0 else 'right'
    if abs(abs_v) > 0.003:
        ax.text(x, i, f'{rel_v:+.0f}%', va='center', ha=ha,
                color='white', fontsize=7, alpha=0.85)

ax.set_yticks(range(len(df_sorted)))
ax.set_yticklabels(df_sorted['ROI_name'], color='white', fontsize=7.5)
ax.axvline(0, color='white', lw=0.8, alpha=0.5)
ax.set_xlabel('MVF_basic − MVF_atlas  (fraction)', color='white', fontsize=11)
ax.set_title(f'Cohort MVF Overestimation — Basic vs Atlas MIMM\n'
             f'(N = {N} subjects, error bars = ±1 SEM propagated, red = basic overestimates)',
             color='white', fontsize=12, fontweight='bold', pad=12)
_spine(ax)
plt.tight_layout()
plt.savefig(os.path.join(COHORT_FIG, '42_cohort_overestimation_ranked.png'),
            dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
plt.close(); print('Saved: 42_cohort_overestimation_ranked.png')

# ── Shared scatter helper ─────────────────────────────────────────────────────
def overest_scatter(xcol, xlabel, colour_col, clabel, cmap, vmin, vmax,
                    highlight, fname):
    valid = df.dropna(subset=[xcol, 'overest_abs'])
    if len(valid) < 5:
        print(f'Skipped: {fname} (insufficient data)')
        return
    r, p = pearsonr(valid[xcol], valid['overest_abs'])
    m, b = np.polyfit(valid[xcol], valid['overest_abs'], 1)
    x_line = np.linspace(valid[xcol].min(), valid[xcol].max(), 100)

    fig, ax = plt.subplots(figsize=(8, 7), facecolor='#0d0d0d')
    ax.set_facecolor('#111111')

    sc = ax.scatter(valid[xcol], valid['overest_abs'],
                    c=valid[colour_col], cmap=cmap, s=60, alpha=0.85,
                    vmin=vmin, vmax=vmax, zorder=3)

    # Error bars on overestimation (between-subject SEM)
    if 'overest_sem' in valid.columns and valid['overest_sem'].notna().any():
        ax.errorbar(valid[xcol], valid['overest_abs'],
                    yerr=valid['overest_sem'],
                    fmt='none', ecolor='#555555', elinewidth=0.6, capsize=2, zorder=2)

    ax.plot(x_line, m * x_line + b, '--', color='cyan', lw=1.5,
            label=f'r = {r:.2f}, p = {p:.3f}', zorder=4)
    ax.axhline(0, color='white', lw=0.8, alpha=0.4)

    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label(clabel, color='white', fontsize=9)
    cbar.ax.yaxis.set_tick_params(color='white')
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color='white')
    cbar.outline.set_edgecolor('white')

    for _, row in valid.iterrows():
        if row['ROI_name'] in highlight:
            short = highlight[row['ROI_name']]
            ax.annotate(short, (row[xcol], row['overest_abs']),
                        textcoords='offset points', xytext=(6, 3),
                        fontsize=7.5, color='white', alpha=0.9,
                        path_effects=[pe.withStroke(linewidth=2, foreground='black')])

    ax.set_xlabel(xlabel, color='white', fontsize=11)
    ax.set_ylabel('MVF_basic − MVF_atlas  (fraction)', color='white', fontsize=11)
    ax.set_title(f'Cohort Overestimation vs {xlabel.split("(")[0].strip()}\n'
                 f'(N = {N} subjects)',
                 color='white', fontsize=13, fontweight='bold')
    ax.tick_params(colors='white')
    for s in ax.spines.values(): s.set_edgecolor('#444444')
    ax.grid(color='#2a2a2a', lw=0.5)
    ax.legend(fontsize=10, facecolor='#1a1a1a', edgecolor='white', labelcolor='white')
    plt.tight_layout()
    plt.savefig(os.path.join(COHORT_FIG, fname),
                dpi=150, bbox_inches='tight', facecolor='#0d0d0d')
    plt.close(); print(f'Saved: {fname}  (r={r:.2f}, p={p:.4f})')

HL_THETA = {
    'Splenium of corpus callosum':           'Splenium CC',
    'Posterior thalamic radiation R':        'Post. thal. rad.',
    'Posterior thalamic radiation L':        'Post. thal. rad.',
    'Posterior limb of internal capsule R':  'Post. IC',
    'Corticospinal tract R':                 'CST',
    'Body of corpus callosum':               'Body CC',
    'Superior longitudinal fasciculus R':    'SLF',
}

# ── Figure 43: Cohort overestimation vs fibre angle θ ────────────────────────
overest_scatter('theta_mean',
                'Mean fibre angle θ relative to B0  (degrees)',
                colour_col='MVF_atlas_mean', clabel='MVF_atlas (reference)',
                cmap='hot', vmin=0.05, vmax=0.40,
                highlight=HL_THETA,
                fname='43_cohort_overestimation_vs_theta.png')

# ── Figure 44: Cohort overestimation vs FA ────────────────────────────────────
overest_scatter('FA_mean',
                'Mean FA (HCP1065 atlas)',
                colour_col='theta_mean', clabel='Mean θ (degrees)',
                cmap='RdBu_r', vmin=0, vmax=90,
                highlight=HL_THETA,
                fname='44_cohort_overestimation_vs_FA.png')

# ── Terminal summary ──────────────────────────────────────────────────────────
print(f'\nTop 5 cohort-overestimated ROIs (basic − atlas):')
print(df.nlargest(5, 'overest_abs')[
    ['ROI_name', 'overest_abs', 'overest_rel', 'theta_mean', 'n_subjects']
].to_string(index=False))
