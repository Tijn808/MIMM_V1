#!/usr/bin/env python3
"""
fig_gp_robustness.py -- iron-robustness figure (defence slide 8).

MIMM MVF vs the T2-GRASE MWF: whole-NAWM baseline, the 5 subcortical grey-matter
structures, and a Lesions group, in that order. The NAWM group shows the two
methods agree in normal white matter; in the iron-rich globus pallidus, iron
shortens T2 and the MWF over-reads myelin, MIMM, which models iron explicitly,
stays low. The GP gap (Cohen's d ~ -3.96) is the finding, highlighted -- NAWM is
what makes it visible as the EXCEPTION, not the rule. Lesions sit last, separated
by a divider, because they are a different phenomenon (tissue loss, not an iron
confound) -- the detailed myelin/axon/fibre lesion decomposition is its own slide.

Reads <results_dir>/cohort_analysis/cohort_gm_rois.csv (run cohort_gm_rois.py first,
it emits MVF_sem / MWF_sem for the error bars) and, if present:
  cohort_nawm_reference.csv    (written by cohort_mwf.py)               -> NAWM group
  cohort_lesion_reference.csv  (written by cohort_lesion_mwf_sensitivity.py) -> Lesions group

Usage: python3 fig_gp_robustness.py <results_dir>
"""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mimm_style import apply_style, C   # shared deck palette / rcParams

apply_style()

if len(sys.argv) < 2:
    sys.exit('usage: fig_gp_robustness.py <results_dir>')
ca = os.path.join(sys.argv[1], 'cohort_analysis')
csv = os.path.join(ca, 'cohort_gm_rois.csv')
if not os.path.exists(csv):
    sys.exit(f'missing {csv}\nRun cohort_gm_rois.py first (it writes MVF_sem/MWF_sem).')
df = pd.read_csv(csv)

# consistent left->right order: whole-NAWM baseline first, then the 5 GM structures
# (most iron-rich, globus pallidus, in the middle of that group), then Lesions last
# (a pathological tissue state, not a GM structure -- kept visually separate).
ORDER = ['NAWM', 'Thalamus', 'Caudate', 'Putamen', 'Globus pallidus', 'Hippocampus', 'Lesions']
df['ord'] = df['structure'].apply(lambda s: ORDER.index(s) if s in ORDER else 99)
df = df.sort_values('ord').reset_index(drop=True)

# prepend the whole-NAWM reference group if available (cohort_mwf.py writes it) --
# shows the two methods AGREE in normal white matter, so the GP gap reads as the
# exception, not business as usual.
nawm_csv = os.path.join(ca, 'cohort_nawm_reference.csv')
if os.path.exists(nawm_csv):
    nawm = pd.read_csv(nawm_csv)
    nawm['ord'] = 0
    df = pd.concat([nawm, df], ignore_index=True)
else:
    print(f'[note] {nawm_csv} not found -- run cohort_mwf.py first to include the NAWM baseline group.')

# append the Lesions group if available (cohort_lesion_mwf_sensitivity.py writes it)
# -- per supervisor request: lesion tissue MVF vs MWF alongside NAWM/GM structures.
lesion_csv = os.path.join(ca, 'cohort_lesion_reference.csv')
if os.path.exists(lesion_csv):
    lesion = pd.read_csv(lesion_csv)
    lesion['ord'] = ORDER.index('Lesions')
    df = pd.concat([df, lesion], ignore_index=True)
else:
    print(f'[note] {lesion_csv} not found -- run cohort_lesion_mwf_sensitivity.py first to include the Lesions group.')

df = df.sort_values('ord').reset_index(drop=True)

x = np.arange(len(df)); w = 0.38
fig, ax = plt.subplots()
b1 = ax.bar(x - w / 2, df['MVF_mean'], w, yerr=df.get('MVF_sem'), capsize=3,
            color=C['MIMM'], edgecolor=C['text'], linewidth=0.6, label='MIMM MVF')
b2 = ax.bar(x + w / 2, df['MWF_mean'], w, yerr=df.get('MWF_sem'), capsize=3,
            color=C['reference'], edgecolor=C['text'], linewidth=0.6, label='T2-GRASE MWF')

# highlight the globus-pallidus gap: the one finding on this figure
gp = df.index[df['structure'] == 'Globus pallidus']
if len(gp):
    i = int(gp[0]); d = df.loc[i, 'cohens_d']
    for b in (b1[i], b2[i]):
        b.set_edgecolor(C['highlight']); b.set_linewidth(2.4)
    top = max(df.loc[i, 'MVF_mean'], df.loc[i, 'MWF_mean'])
    top += (df.loc[i, 'MWF_sem'] if 'MWF_sem' in df else 0) + 0.03
    ax.annotate(f"d = {d:+.2f}", xy=(i, top), ha='center', va='bottom',
                color=C['highlight'], fontweight='bold', fontsize=11)

ax.set_xticks(x); ax.set_xticklabels(df['structure'], rotation=20, ha='right')
ax.set_ylabel('Myelin fraction')

# separators between the three groups (NAWM | GM structures | Lesions), so the
# grouping reads clearly: "agreement here" vs "the iron-rich structures" vs
# "pathological tissue" (a different phenomenon, kept visually distinct)
if 'NAWM' in df['structure'].values:
    sep_x = df.index[df['structure'] == 'NAWM'][0] + 0.5
    ax.axvline(sep_x, color=C['reference'], lw=1, ls=':', alpha=0.6)
if 'Lesions' in df['structure'].values:
    sep_x = df.index[df['structure'] == 'Lesions'][0] - 0.5
    ax.axvline(sep_x, color=C['reference'], lw=1, ls=':', alpha=0.6)
# no baked-in title: the slide headline carries the (qualified) interpretive claim,
# and the axes alone cannot show MIMM is "correct" — only that the two diverge.
ax.legend(loc='upper left')
ymax = float(np.nanmax([df['MVF_mean'].max(), df['MWF_mean'].max()]))
ax.set_ylim(0, ymax * 1.28)

out = os.path.join(ca, 'cohort_gp_robustness.png')
fig.savefig(out)
fig.savefig(out.replace('.png', '.svg'))
print('saved', out, 'and .svg')
