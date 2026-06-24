#!/usr/bin/env python3
"""
Is MIMM more sensitive than MWF to demyelination?  Head-to-head: for every patient
with a lesion mask AND a registered T2-GRASE MWF map, measure both MIMM MVF and
MWF in lesion vs NAWM, paired within subject, and compare the effect sizes.

"More sensitive" = a larger, more significant lesion-vs-NAWM drop. We report, for
each metric, the % change, the paired t-test p, and the paired Cohen's d
(mean within-subject drop / SD of that drop) — d is the scale-free sensitivity
measure that lets MVF and MWF be compared directly despite different units.

Inputs (per subject):
  lesion/lesion_mask.nii.gz, qsm/brain_mask.nii.gz, atlas/FA_atlas.nii.gz,
  mimm/MVF_Atlas.nii.gz, grase/MWF.nii.gz

Outputs to <results>/cohort_analysis/:
  cohort_lesion_mwf_sensitivity.png
  cohort_lesion_mwf_sensitivity.csv

Usage:  python3 cohort_lesion_mwf_sensitivity.py <results_dir>
"""
import sys, os, glob, csv
import numpy as np
import nibabel as nib
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    sys.exit('usage: cohort_lesion_mwf_sensitivity.py <results_dir>')
results_dir = sys.argv[1]
out_dir = os.path.join(results_dir, 'cohort_analysis'); os.makedirs(out_dir, exist_ok=True)


def load(p):
    return np.asarray(nib.load(p).dataobj) if os.path.exists(p) else None

# (label, relative path, transform, clip)
METRICS = [('MIMM MVF', 'mimm/MVF_Atlas.nii.gz', None),
           ('MWF',      'grase/MWF.nii.gz',      (0, 0.5))]

rows = []
for ld in sorted(glob.glob(os.path.join(results_dir, '*', 'lesion', 'lesion_mask.nii.gz'))):
    d = os.path.dirname(os.path.dirname(ld)); sid = os.path.basename(d)
    lesion = load(ld); brain = load(os.path.join(d, 'qsm', 'brain_mask.nii.gz'))
    fa = load(os.path.join(d, 'atlas', 'FA_atlas.nii.gz'))
    vols = {lab: load(os.path.join(d, rel)) for lab, rel, _ in METRICS}
    missing = [lab for lab, v in vols.items() if v is None] + \
              [n for n, v in (('lesion', lesion), ('brain', brain), ('FA', fa)) if v is None]
    if missing:
        print(f'[skip] {sid} (missing: {", ".join(missing)})'); continue
    lesion = (lesion > 0.5) & (brain > 0)
    if lesion.sum() < 50:
        print(f'[skip] {sid} ({int(lesion.sum())} lesion vox < 50)'); continue
    nawm = (brain > 0) & (fa > 0.20) & ~lesion
    r = {'subject': sid}
    for lab, rel, clip in METRICS:
        v = vols[lab].astype(float)
        if clip:
            v = np.clip(v, *clip)
        lv = v[lesion]; nv = v[nawm]
        lv = lv[np.isfinite(lv) & (lv != 0)]; nv = nv[np.isfinite(nv) & (nv != 0)]
        r[lab] = (float(np.mean(nv)), float(np.mean(lv)))   # (NAWM, lesion)
    rows.append(r)
    print(f'{sid}: ' + '  '.join(f'{lab} {r[lab][0]:.3f}->{r[lab][1]:.3f}' for lab, _, _ in METRICS))

if not rows:
    sys.exit('No subjects with both a lesion mask and an MWF map.')
n = len(rows)
print(f'\n{n} patients with both lesion mask and MWF.')

# --- paired stats per metric ---
fig, axes = plt.subplots(1, len(METRICS), figsize=(5 * len(METRICS), 5))
summary = []
for ax, (lab, _, _) in zip(axes, METRICS):
    nawm = np.array([r[lab][0] for r in rows]); les = np.array([r[lab][1] for r in rows])
    m = np.isfinite(nawm) & np.isfinite(les); nawm, les = nawm[m], les[m]
    diff = les - nawm
    t, p = stats.ttest_rel(nawm, les)
    pct = 100 * diff.mean() / nawm.mean()
    d = diff.mean() / diff.std(ddof=1)               # paired Cohen's d
    summary.append({'metric': lab, 'n': int(m.sum()),
                    'NAWM': round(nawm.mean(), 4), 'lesion': round(les.mean(), 4),
                    'pct_change': round(pct, 1), 'paired_p': float(p),
                    'cohens_d': round(float(d), 2)})
    for a, b in zip(nawm, les):
        ax.plot([0, 1], [a, b], '-', color='0.6', lw=0.8, alpha=0.7)
    ax.plot([0]*len(nawm), nawm, 'o', color='#1f77b4'); ax.plot([1]*len(les), les, 'o', color='#d62728')
    for xi, v in [(0, nawm), (1, les)]:
        ax.errorbar(xi, v.mean(), yerr=v.std(ddof=1)/np.sqrt(len(v)), fmt='_', color='k', ms=22, capsize=5, lw=2)
    ax.set_xticks([0, 1]); ax.set_xticklabels(['NAWM', 'lesion']); ax.set_xlim(-0.4, 1.4)
    ax.set_title(f'{lab}\n{pct:+.1f}%, ' + ('p<0.001' if p < 1e-3 else f'p={p:.3f}') +
                 f'\nCohen d = {d:+.2f}  (n={int(m.sum())})')
    ax.grid(alpha=0.25, axis='y')
fig.suptitle(f'Sensitivity to demyelination: MIMM MVF vs MWF  (same {n} patients)', fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(out_dir, 'cohort_lesion_mwf_sensitivity.png'), dpi=150)

with open(os.path.join(out_dir, 'cohort_lesion_mwf_sensitivity.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(summary[0].keys())); w.writeheader(); w.writerows(summary)

print('\n=== Sensitivity summary (|Cohen d| larger = more sensitive) ===')
for s in summary:
    print(f'  {s["metric"]:10s} {s["pct_change"]:+6.1f}%  p={s["paired_p"]:.3g}  d={s["cohens_d"]:+.2f}')
mvf_d = abs(next(s['cohens_d'] for s in summary if s['metric'] == 'MIMM MVF'))
mwf_d = abs(next(s['cohens_d'] for s in summary if s['metric'] == 'MWF'))
verdict = 'MIMM MVF' if mvf_d > mwf_d else 'MWF'
print(f'\n-> {verdict} shows the larger lesion effect size on these {n} patients '
      f'(|d| {max(mvf_d, mwf_d):.2f} vs {min(mvf_d, mwf_d):.2f}).')
print(f'saved: cohort_lesion_mwf_sensitivity.png / .csv')
