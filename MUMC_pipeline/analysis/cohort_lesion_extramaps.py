#!/usr/bin/env python3
"""
Extra lesion-vs-NAWM contrasts not covered by the main analysis, to look for new
effects. For each map, lesion is compared to BOTH global NAWM and a location-
matched perilesional NAWM shell (paired across patients), with % change, p and
Cohen's d.

Maps:
  MVF, FVF   -> demyelination vs axon loss. If MVF drops but FVF holds, the lesion
                is demyelination-dominant (axons preserved) — something MWF cannot
                show. If FVF drops too, there is axonal loss.
  MIMM iron, chi-sep iron -> does iron change in lesions (not just at the rim)?
  R2*        -> a second, independent iron/relaxation marker in lesions.

Usage:  python3 cohort_lesion_extramaps.py <results_dir>
"""
import sys, os, glob, csv
import numpy as np
import nibabel as nib
from scipy import ndimage, stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

GAP, WIDTH, MIN_VOX = 2, 3, 50

if len(sys.argv) < 2:
    sys.exit('usage: cohort_lesion_extramaps.py <results_dir>')
results_dir = sys.argv[1]
out_dir = os.path.join(results_dir, 'cohort_analysis'); os.makedirs(out_dir, exist_ok=True)


def load(p):
    return np.asarray(nib.load(p).dataobj) if os.path.exists(p) else None

MAPS = [('MVF',          'mimm/MVF_Atlas.nii.gz'),
        ('FVF',          'mimm/FVF_Atlas.nii.gz'),
        ('MIMM iron',    'mimm/chi_iron_est_Atlas.nii.gz'),
        ('chi-sep iron', 'chisep/chi_pos.nii.gz'),
        ('R2*',          'mimm/R2s_Atlas.nii.gz')]

rows = []
for ld in sorted(glob.glob(os.path.join(results_dir, '*', 'lesion', 'lesion_mask.nii.gz'))):
    d = os.path.dirname(os.path.dirname(ld)); sid = os.path.basename(d)
    lesion = load(ld); brain = load(os.path.join(d, 'qsm', 'brain_mask.nii.gz'))
    fa = load(os.path.join(d, 'atlas', 'FA_atlas.nii.gz'))
    if any(v is None for v in (lesion, brain, fa)):
        print(f'[skip] {sid} (missing lesion/brain/FA)'); continue
    lesion = (lesion > 0.5) & (brain > 0)
    if lesion.sum() < MIN_VOX:
        print(f'[skip] {sid} ({int(lesion.sum())} lesion vox)'); continue
    nawm  = (brain > 0) & (fa > 0.20) & ~lesion
    inner = ndimage.binary_dilation(lesion, iterations=GAP)
    outer = ndimage.binary_dilation(lesion, iterations=GAP + WIDTH)
    peri  = outer & ~inner & nawm
    if peri.sum() < MIN_VOX:
        print(f'[skip] {sid} (peri shell {int(peri.sum())} vox)'); continue
    r = {'subject': sid}
    for label, rel in MAPS:
        vol = load(os.path.join(d, rel))
        if vol is None:
            r[label] = (np.nan, np.nan, np.nan); continue
        vol = vol.astype(float)
        def mean(m):
            v = vol[m]; v = v[np.isfinite(v)]; return float(np.mean(v)) if v.size else np.nan
        r[label] = (mean(nawm), mean(lesion), mean(peri))   # (global NAWM, lesion, peri)
    rows.append(r)
    print(f'{sid}: MVF {r["MVF"][1]:.3f} (NAWM {r["MVF"][0]:.3f}), '
          f'FVF {r["FVF"][1]:.3f} (NAWM {r["FVF"][0]:.3f})')

if not rows:
    sys.exit('No usable subjects.')
n = len(rows)
print(f'\n{n} patients.')


def paired(les, ref):
    m = np.isfinite(les) & np.isfinite(ref); les, ref = les[m], ref[m]
    if m.sum() < 3:
        return np.nan, np.nan, np.nan, int(m.sum())
    diff = les - ref
    p = stats.ttest_rel(les, ref).pvalue
    pct = 100 * diff.mean() / ref.mean()
    d = diff.mean() / diff.std(ddof=1)
    return pct, float(p), float(d), int(m.sum())

summary = []
for label, _ in MAPS:
    glob = np.array([r[label][0] for r in rows])
    les  = np.array([r[label][1] for r in rows])
    peri = np.array([r[label][2] for r in rows])
    pg, pg_p, dg, _ = paired(les, glob)
    pp, pp_p, dp, n_ok = paired(les, peri)
    summary.append({'map': label, 'global_pct': round(pg, 1), 'global_p': pg_p, 'global_d': round(dg, 2),
                    'peri_pct': round(pp, 1), 'peri_p': pp_p, 'peri_d': round(dp, 2), 'n': n_ok})

with open(os.path.join(out_dir, 'cohort_lesion_extramaps.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(summary[0].keys())); w.writeheader(); w.writerows(summary)

print('\n=== lesion vs NAWM — extra maps (global | perilesional) ===')
print(f'{"map":14s} {"glob %":>7} {"glob p":>9} {"glob d":>7}   {"peri %":>7} {"peri p":>9} {"peri d":>7}')
for s in summary:
    star_g = '*' if s['global_p'] < 0.05 else ' '
    star_p = '*' if s['peri_p'] < 0.05 else ' '
    print(f'{s["map"]:14s} {s["global_pct"]:>6.1f}% {s["global_p"]:>9.3g}{star_g}{s["global_d"]:>6.2f}   '
          f'{s["peri_pct"]:>6.1f}% {s["peri_p"]:>9.3g}{star_p}{s["peri_d"]:>6.2f}')

# --- demyelination vs axon-loss readout ---
mvf = next(s for s in summary if s['map'] == 'MVF')
fvf = next(s for s in summary if s['map'] == 'FVF')
print('\n--- Demyelination vs axon loss (perilesional, location-matched) ---')
print(f'  MVF {mvf["peri_pct"]:+.1f}% (p={mvf["peri_p"]:.3g}),  FVF {fvf["peri_pct"]:+.1f}% (p={fvf["peri_p"]:.3g})')
if fvf['peri_p'] >= 0.05 or abs(fvf['peri_pct']) < abs(mvf['peri_pct']) / 2:
    print('  -> MVF drops but FVF is preserved/weaker: demyelination-dominant '
          '(axons relatively spared). This is a MIMM-specific readout.')
else:
    print('  -> FVF drops alongside MVF: evidence of axonal loss, not pure demyelination.')

# --- figure: perilesional comparison for the extra maps ---
plot_maps = ['FVF', 'MIMM iron', 'chi-sep iron', 'R2*']
fig, axes = plt.subplots(1, len(plot_maps), figsize=(4.3 * len(plot_maps), 5))
for ax, label in zip(axes, plot_maps):
    les  = np.array([r[label][1] for r in rows])
    peri = np.array([r[label][2] for r in rows])
    m = np.isfinite(les) & np.isfinite(peri); les, peri = les[m], peri[m]
    for a, b in zip(peri, les):
        ax.plot([0, 1], [a, b], '-', color='0.6', lw=0.8, alpha=0.7)
    ax.plot([0]*len(peri), peri, 'o', color='#2ca02c'); ax.plot([1]*len(les), les, 'o', color='#d62728')
    for xi, v in [(0, peri), (1, les)]:
        ax.errorbar(xi, v.mean(), yerr=v.std(ddof=1)/np.sqrt(len(v)), fmt='_', color='k', ms=20, capsize=5, lw=2)
    s = next(x for x in summary if x['map'] == label)
    ax.set_xticks([0, 1]); ax.set_xticklabels(['peri-\nNAWM', 'lesion']); ax.set_xlim(-0.4, 1.4)
    ax.set_title(f'{label}\n{s["peri_pct"]:+.1f}%, ' +
                 ('p<0.001' if s['peri_p'] < 1e-3 else f'p={s["peri_p"]:.3f}') + f', d={s["peri_d"]:+.2f}', fontsize=10)
    ax.grid(alpha=0.25, axis='y')
fig.suptitle(f'Extra lesion maps vs perilesional NAWM (n={n})', fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(os.path.join(out_dir, 'cohort_lesion_extramaps.png'), dpi=150)
print(f'\nsaved: cohort_lesion_extramaps.png / .csv')
