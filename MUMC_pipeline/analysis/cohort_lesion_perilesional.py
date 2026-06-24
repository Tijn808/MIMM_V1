#!/usr/bin/env python3
"""
Location-matched lesion analysis: compare lesion tissue to a PERILESIONAL NAWM
shell (normal-appearing white matter just outside the lesion) instead of global
NAWM. This controls for lesion LOCATION — MS lesions cluster periventricularly,
and baseline g-ratio/MVF vary by region, so a global-NAWM comparison can mix a
"where lesions sit" effect into the "demyelination" effect. A perilesional rim in
the same anatomical neighbourhood removes that confound.

For each patient and each map (MVF, |chi_neg|, g-ratio):
  lesion        = mask & brain
  perilesional  = a shell `GAP`..`GAP+WIDTH` voxels outside the lesion, kept to
                  NAWM (FA>0.20, not lesion). The GAP skips the lesion-edge
                  partial-volume voxels so the control is clean.
  global NAWM   = brain & FA>0.20 & not-lesion   (the original §3.4 reference)

Paired t-tests across patients, with % change, p and Cohen's d, for BOTH the
perilesional and the global comparison, so the two can be read side by side.

Usage:  python3 cohort_lesion_perilesional.py <results_dir>
"""
import sys, os, glob, csv
import numpy as np
import nibabel as nib
from scipy import ndimage, stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

GAP, WIDTH, MIN_VOX = 2, 3, 50      # shell = 2..5 voxels out; min voxels per region

if len(sys.argv) < 2:
    sys.exit('usage: cohort_lesion_perilesional.py <results_dir>')
results_dir = sys.argv[1]
out_dir = os.path.join(results_dir, 'cohort_analysis'); os.makedirs(out_dir, exist_ok=True)


def load(p):
    return np.asarray(nib.load(p).dataobj) if os.path.exists(p) else None

MAPS = [('MVF (Atlas)', 'mimm/MVF_Atlas.nii.gz', lambda x: x),
        ('|chi_neg|',   'chisep/chi_neg.nii.gz', np.abs),
        ('g-ratio',     'mimm/g_ratio_Atlas.nii.gz', lambda x: x)]

rows = []
for ld in sorted(glob.glob(os.path.join(results_dir, '*', 'lesion', 'lesion_mask.nii.gz'))):
    d = os.path.dirname(os.path.dirname(ld)); sid = os.path.basename(d)
    lesion = load(ld); brain = load(os.path.join(d, 'qsm', 'brain_mask.nii.gz'))
    fa = load(os.path.join(d, 'atlas', 'FA_atlas.nii.gz'))
    if any(v is None for v in (lesion, brain, fa)):
        print(f'[skip] {sid} (missing lesion/brain/FA)'); continue
    lesion = (lesion > 0.5) & (brain > 0)
    if lesion.sum() < MIN_VOX:
        print(f'[skip] {sid} ({int(lesion.sum())} lesion vox < {MIN_VOX})'); continue
    nawm   = (brain > 0) & (fa > 0.20) & ~lesion
    inner  = ndimage.binary_dilation(lesion, iterations=GAP)
    outer  = ndimage.binary_dilation(lesion, iterations=GAP + WIDTH)
    peri   = outer & ~inner & nawm
    if peri.sum() < MIN_VOX:
        print(f'[skip] {sid} (perilesional shell {int(peri.sum())} vox < {MIN_VOX})'); continue
    r = {'subject': sid, 'n_lesion': int(lesion.sum()), 'n_peri': int(peri.sum())}
    for label, rel, fn in MAPS:
        vol = load(os.path.join(d, rel))
        if vol is None:
            r[label] = (np.nan, np.nan, np.nan); continue
        vol = fn(vol.astype(float))
        def mean(m):
            v = vol[m]; v = v[np.isfinite(v)]; return float(np.mean(v)) if v.size else np.nan
        r[label] = (mean(nawm), mean(lesion), mean(peri))   # (global NAWM, lesion, perilesional)
    rows.append(r)
    print(f'{sid}: lesion {r["n_lesion"]} vox, peri {r["n_peri"]} vox  '
          f'g-ratio glob/les/peri ' + '/'.join(f'{v:.3f}' for v in r['g-ratio']))

if not rows:
    sys.exit('No subjects with a usable lesion + perilesional shell.')
n = len(rows)
print(f'\n{n} patients with a lesion and a perilesional NAWM shell.')


def paired(les, ref):
    m = np.isfinite(les) & np.isfinite(ref); les, ref = les[m], ref[m]
    diff = les - ref
    t, p = stats.ttest_rel(les, ref)
    pct = 100 * diff.mean() / ref.mean()
    d = diff.mean() / diff.std(ddof=1)
    return pct, float(p), float(d), int(m.sum())

# --- figure: perilesional (location-matched) paired comparison, 3 maps ---
fig, axes = plt.subplots(1, len(MAPS), figsize=(4.6 * len(MAPS), 5))
summary = []
for ax, (label, _, _) in zip(axes, MAPS):
    glob = np.array([r[label][0] for r in rows])
    les  = np.array([r[label][1] for r in rows])
    peri = np.array([r[label][2] for r in rows])
    pct_g, p_g, d_g, _ = paired(les, glob)     # original: lesion vs global NAWM
    pct_p, p_p, d_p, n_ok = paired(les, peri)   # location-matched: lesion vs perilesional
    summary.append({'map': label,
                    'global_pct': round(pct_g, 1), 'global_p': p_g, 'global_d': round(d_g, 2),
                    'peri_pct': round(pct_p, 1), 'peri_p': p_p, 'peri_d': round(d_p, 2), 'n': n_ok})
    # plot the location-matched comparison
    mok = np.isfinite(les) & np.isfinite(peri)
    for a, b in zip(peri[mok], les[mok]):
        ax.plot([0, 1], [a, b], '-', color='0.6', lw=0.8, alpha=0.7)
    ax.plot([0]*mok.sum(), peri[mok], 'o', color='#2ca02c')
    ax.plot([1]*mok.sum(), les[mok], 'o', color='#d62728')
    for xi, v in [(0, peri[mok]), (1, les[mok])]:
        ax.errorbar(xi, v.mean(), yerr=v.std(ddof=1)/np.sqrt(len(v)), fmt='_', color='k', ms=20, capsize=5, lw=2)
    ax.set_xticks([0, 1]); ax.set_xticklabels(['peri-\nNAWM', 'lesion']); ax.set_xlim(-0.4, 1.4)
    ax.set_title(f'{label}\nlocation-matched: {pct_p:+.1f}%, ' +
                 ('p<0.001' if p_p < 1e-3 else f'p={p_p:.3f}') + f', d={d_p:+.2f}\n'
                 f'(global NAWM was {pct_g:+.1f}%, p={p_g:.3f})', fontsize=9)
    ax.grid(alpha=0.25, axis='y')
fig.suptitle(f'Lesion vs PERILESIONAL NAWM (location-matched), n={n}', fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(os.path.join(out_dir, 'cohort_lesion_perilesional.png'), dpi=150)

with open(os.path.join(out_dir, 'cohort_lesion_perilesional.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(summary[0].keys())); w.writeheader(); w.writerows(summary)

print('\n=== lesion vs GLOBAL NAWM  vs  lesion vs PERILESIONAL NAWM ===')
print(f'{"map":12s} {"global %":>9} {"glob p":>8} {"glob d":>7}   {"peri %":>8} {"peri p":>8} {"peri d":>7}')
for s in summary:
    print(f'{s["map"]:12s} {s["global_pct"]:>8.1f}% {s["global_p"]:>8.3g} {s["global_d"]:>7.2f}   '
          f'{s["peri_pct"]:>7.1f}% {s["peri_p"]:>8.3g} {s["peri_d"]:>7.2f}')
print('\nIf an effect shrinks or loses significance from global -> perilesional, '
      'it was partly a location/regional baseline effect, not pure demyelination.')
print(f'saved: cohort_lesion_perilesional.png / .csv')
