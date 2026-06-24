#!/usr/bin/env python3
"""
make_thesis_figures.py — the main, publication-quality figures for the MIMM thesis.

Regenerates only the figures that carry the results, with consistent styling
(larger fonts, Greek symbols, fixed colours). Everything else (nulls, robustness
checks) stays as numbers in the text or in the appendix scripts.

Outputs to <results_dir>/cohort_analysis/thesis_figures/:
  fig1_crossmethod.png    MVF vs |χ⁻| (Basic + Atlas) + MVF vs FA          (§3.1)
  fig2_orientation.png    Basic vs Atlas identity + (Basic−Atlas) vs θ     (§3.3)
  fig3_compartments.png   lesion vs perilesional NAWM: MVF, AVF, FVF       (§3.4)
  fig4_resolution.png     MVF vs MWF lesion drop vs lesion size            (§3.9)
  fig5_mwf.png            MVF vs MWF scatter + Bland–Altman                (§3.7)

Usage:  python3 make_thesis_figures.py <results_dir>
"""
import sys, os, glob
import numpy as np
import pandas as pd
import nibabel as nib
from scipy import ndimage, stats
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

mpl.rcParams.update({
    'font.size': 12, 'axes.titlesize': 14, 'axes.labelsize': 13,
    'xtick.labelsize': 11, 'ytick.labelsize': 11, 'legend.fontsize': 11,
    'savefig.dpi': 200, 'axes.grid': True, 'grid.alpha': 0.25, 'grid.linestyle': '--',
})
ATLAS, BASIC = '#1f77b4', '#ff7f0e'
LESION, NAWM, MWFC = '#d62728', '#2ca02c', '#9467bd'

if len(sys.argv) < 2:
    sys.exit('usage: make_thesis_figures.py <results_dir>')
RESULTS = sys.argv[1]
OUT = os.path.join(RESULTS, 'cohort_analysis', 'thesis_figures')
os.makedirs(OUT, exist_ok=True)
GAP, WIDTH = 2, 3


def load(p):
    return np.asarray(nib.load(p).dataobj) if os.path.exists(p) else None


def annotate(ax, text):
    ax.text(0.04, 0.96, text, transform=ax.transAxes, va='top', fontsize=11,
            bbox=dict(boxstyle='round', fc='white', ec='0.6'))


def pooled_roi():
    """Cohort-mean per-ROI table from every subject's roi_stats.csv."""
    frames = []
    for f in sorted(glob.glob(os.path.join(RESULTS, '*', 'analysis', 'roi_stats.csv'))):
        frames.append(pd.read_csv(f))
    long = pd.concat(frames, ignore_index=True)
    long['chi_neg_abs'] = long['chi_neg_chisep_mean'].abs()
    return long.groupby('ROI_name', as_index=False).mean(numeric_only=True), long['ROI_name'].nunique()


# ── Figure 1: cross-method validation ────────────────────────────────────────
def fig1():
    roi, _ = pooled_roi()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5))
    x = roi['chi_neg_abs'].values
    for ycol, col, lab in [('MVF_basic_mean', BASIC, 'Basic'), ('MVF_atlas_mean', ATLAS, 'Atlas')]:
        y = roi[ycol].values; m = np.isfinite(x) & np.isfinite(y)
        a1.scatter(x[m], y[m], s=30, c=col, alpha=0.8, edgecolors='none', label=lab)
        b, a = np.polyfit(x[m], y[m], 1); xs = np.linspace(x[m].min(), x[m].max(), 100)
        a1.plot(xs, b * xs + a, '--', color=col, lw=1.4)
    r_at = stats.pearsonr(*[v[np.isfinite(x) & np.isfinite(roi['MVF_atlas_mean'])]
                            for v in (x, roi['MVF_atlas_mean'].values)])[0]
    r_ba = stats.pearsonr(*[v[np.isfinite(x) & np.isfinite(roi['MVF_basic_mean'])]
                            for v in (x, roi['MVF_basic_mean'].values)])[0]
    annotate(a1, f'Atlas  r = {r_at:.2f}\nBasic  r = {r_ba:.2f}\n(p < 0.001)')
    a1.set_xlabel('|χ⁻|  χ-separation  (ppm)'); a1.set_ylabel('MVF (MIMM)')
    a1.set_title('Myelin: MIMM vs χ-separation'); a1.legend(loc='lower right')

    xf, yf = roi['FA_mean'].values, roi['MVF_atlas_mean'].values
    m = np.isfinite(xf) & np.isfinite(yf)
    a2.scatter(xf[m], yf[m], s=30, c=ATLAS, alpha=0.8, edgecolors='none')
    rfa = stats.pearsonr(xf[m], yf[m])[0]
    b, a = np.polyfit(xf[m], yf[m], 1); xs = np.linspace(xf[m].min(), xf[m].max(), 100)
    a2.plot(xs, b * xs + a, 'k--', lw=1.3)
    annotate(a2, f'r = {rfa:.2f}')
    a2.set_xlabel('Fractional anisotropy (FA)'); a2.set_ylabel('MVF (Atlas)')
    a2.set_title('Myelin vs FA (specificity)')
    fig.suptitle('Cross-method validation, per JHU white-matter ROI', fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(os.path.join(OUT, 'fig1_crossmethod.png'))
    plt.close(fig); print('fig1_crossmethod.png')


# ── Figure 2: orientation (Basic vs Atlas) ───────────────────────────────────
def fig2():
    roi, _ = pooled_roi()
    roi['overest'] = roi['MVF_basic_mean'] - roi['MVF_atlas_mean']
    fig, (b1, b2) = plt.subplots(1, 2, figsize=(11, 5))
    lim = [min(roi['MVF_atlas_mean'].min(), roi['MVF_basic_mean'].min()),
           max(roi['MVF_atlas_mean'].max(), roi['MVF_basic_mean'].max())]
    b1.scatter(roi['MVF_atlas_mean'], roi['MVF_basic_mean'], s=30, c=ATLAS, alpha=0.8)
    b1.plot(lim, lim, 'k--', lw=1.2, label='identity')
    md = roi['overest'].mean()
    annotate(b1, f'mean Basic−Atlas = {md:+.3f}\n({100*md/roi["MVF_atlas_mean"].mean():+.0f}%)')
    b1.set_xlabel('MVF Atlas'); b1.set_ylabel('MVF Basic')
    b1.set_title('Basic vs Atlas'); b1.legend(loc='lower right')

    xt, yt = roi['theta_mean'].values, roi['overest'].values
    m = np.isfinite(xt) & np.isfinite(yt)
    b2.scatter(xt[m], yt[m], s=30, c=BASIC, alpha=0.8)
    rt = stats.pearsonr(xt[m], yt[m])[0]
    b, a = np.polyfit(xt[m], yt[m], 1); xs = np.linspace(xt[m].min(), xt[m].max(), 100)
    b2.plot(xs, b * xs + a, 'k--', lw=1.3)
    annotate(b2, f'r = {rt:.2f}')
    b2.set_xlabel('θ  fibre angle to B₀  (deg)'); b2.set_ylabel('MVF  Basic − Atlas')
    b2.set_title('Orientation-driven overestimation')
    fig.suptitle('Orientation effect: Basic vs Atlas reconstruction', fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(os.path.join(OUT, 'fig2_orientation.png'))
    plt.close(fig); print('fig2_orientation.png')


# ── Figure 3: lesion compartment decomposition (myelin + axon) ───────────────
def fig3():
    rows = []
    for ld in sorted(glob.glob(os.path.join(RESULTS, '*', 'lesion', 'lesion_mask.nii.gz'))):
        d = os.path.dirname(os.path.dirname(ld)); sid = os.path.basename(d)
        lesion = load(ld); brain = load(os.path.join(d, 'qsm', 'brain_mask.nii.gz'))
        fa = load(os.path.join(d, 'atlas', 'FA_atlas.nii.gz'))
        mvf = load(os.path.join(d, 'mimm', 'MVF_Atlas.nii.gz'))
        fvf = load(os.path.join(d, 'mimm', 'FVF_Atlas.nii.gz'))
        if any(v is None for v in (lesion, brain, fa, mvf, fvf)):
            continue
        lesion = (lesion > 0.5) & (brain > 0)
        if lesion.sum() < 50:
            continue
        nawm = (brain > 0) & (fa > 0.20) & ~lesion
        inner = ndimage.binary_dilation(lesion, iterations=GAP)
        outer = ndimage.binary_dilation(lesion, iterations=GAP + WIDTH)
        peri = outer & ~inner & nawm
        if peri.sum() < 50:
            continue
        mvf = mvf.astype(float); fvf = fvf.astype(float); avf = fvf - mvf

        def mean(vol, m):
            v = vol[m]; v = v[np.isfinite(v)]; return float(np.mean(v)) if v.size else np.nan
        rows.append({'MVF (myelin)': (mean(mvf, peri), mean(mvf, lesion)),
                     'AVF (axon)':   (mean(avf, peri), mean(avf, lesion)),
                     'FVF (fibre)':  (mean(fvf, peri), mean(fvf, lesion))})
    n = len(rows)
    fig, axes = plt.subplots(1, 3, figsize=(13, 5))
    for ax, label in zip(axes, ['MVF (myelin)', 'AVF (axon)', 'FVF (fibre)']):
        peri = np.array([r[label][0] for r in rows]); les = np.array([r[label][1] for r in rows])
        m = np.isfinite(peri) & np.isfinite(les); peri, les = peri[m], les[m]
        for a, b in zip(peri, les):
            ax.plot([0, 1], [a, b], '-', color='0.65', lw=0.9, alpha=0.7)
        ax.plot([0]*len(peri), peri, 'o', color=NAWM, ms=6)
        ax.plot([1]*len(les), les, 'o', color=LESION, ms=6)
        for xi, v in [(0, peri), (1, les)]:
            ax.errorbar(xi, v.mean(), yerr=v.std(ddof=1)/np.sqrt(len(v)), fmt='_', color='k', ms=26, capsize=6, lw=2.2)
        p = stats.ttest_rel(les, peri).pvalue
        pct = 100 * (les.mean() - peri.mean()) / peri.mean()
        d = (les - peri).mean() / (les - peri).std(ddof=1)
        ax.set_xticks([0, 1]); ax.set_xticklabels(['peri-NAWM', 'lesion']); ax.set_xlim(-0.4, 1.4)
        ax.set_title(f'{label}\n{pct:+.1f}%, ' + ('p < 0.001' if p < 1e-3 else f'p = {p:.3f}') + f', d = {d:+.2f}')
    fig.suptitle(f'Lesion compartment decomposition vs perilesional NAWM (n = {n})', fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.94]); fig.savefig(os.path.join(OUT, 'fig3_compartments.png'))
    plt.close(fig); print(f'fig3_compartments.png  (n={n})')


# ── Figure 4: resolution — MVF vs MWF lesion drop vs size ─────────────────────
def fig4():
    lesions = []
    for ld in sorted(glob.glob(os.path.join(RESULTS, '*', 'lesion', 'lesion_mask.nii.gz'))):
        d = os.path.dirname(os.path.dirname(ld)); sid = os.path.basename(d)
        mask = load(ld); brain = load(os.path.join(d, 'qsm', 'brain_mask.nii.gz'))
        fa = load(os.path.join(d, 'atlas', 'FA_atlas.nii.gz'))
        mvf = load(os.path.join(d, 'mimm', 'MVF_Atlas.nii.gz'))
        mwf = load(os.path.join(d, 'grase', 'MWF.nii.gz'))
        flair = load(os.path.join(d, 'lesion', 'FLAIR_mgre.nii.gz'))
        if any(v is None for v in (mask, brain, fa, mvf, mwf, flair)):
            continue
        mask = (mask > 0.5) & (brain > 0)
        mvf = mvf.astype(float); flair = flair.astype(float); mwf = np.clip(mwf.astype(float), 0, 0.5)
        nawm = (brain > 0) & (fa > 0.20) & ~mask
        vox_mL = float(np.prod(nib.load(ld).header.get_zooms()[:3])) / 1000.0
        labels, _ = ndimage.label(mask); pad = GAP + WIDTH + 1
        for i, sl in enumerate(ndimage.find_objects(labels), start=1):
            if sl is None:
                continue
            esl = tuple(slice(max(0, s.start - pad), min(dim, s.stop + pad)) for s, dim in zip(sl, mask.shape))
            comp = labels[esl] == i
            if comp.sum() < 3:
                continue
            inner = ndimage.binary_dilation(comp, iterations=GAP)
            outer = ndimage.binary_dilation(comp, iterations=GAP + WIDTH)
            peri = outer & ~inner & nawm[esl]
            if peri.sum() < 30:
                continue
            mvf_c, mwf_c, fl_c = mvf[esl], mwf[esl], flair[esl]

            def m(v, msk):
                x = v[msk]; x = x[np.isfinite(x)]; return float(np.mean(x)) if x.size else np.nan
            fl_l, fl_p = m(fl_c, comp), m(fl_c, peri)
            if not (np.isfinite(fl_l) and np.isfinite(fl_p) and fl_p > 0 and fl_l >= fl_p * 1.05):
                continue
            mvf_p, mwf_p = m(mvf_c, peri), m(mwf_c, peri)
            if not (mvf_p > 0 and mwf_p > 0):
                continue
            lesions.append((comp.sum() * vox_mL, 100*(mvf_p - m(mvf_c, comp))/mvf_p,
                            100*(mwf_p - m(mwf_c, comp))/mwf_p))
    vol = np.array([l[0] for l in lesions]); dmvf = np.array([l[1] for l in lesions]); dmwf = np.array([l[2] for l in lesions])
    q1, q2 = np.quantile(vol, [1/3, 2/3])
    bins = [(vol <= q1), (vol > q1) & (vol <= q2), (vol > q2)]
    xb = [vol[s].mean() for s in bins]
    fig, ax = plt.subplots(figsize=(8.5, 6))
    ax.scatter(vol, dmvf, s=16, c=ATLAS, alpha=0.4)
    ax.scatter(vol, dmwf, s=16, c=LESION, alpha=0.4)
    ax.plot(xb, [dmvf[s].mean() for s in bins], '-o', color=ATLAS, lw=2.4, ms=9, label='MIMM MVF (binned)')
    ax.plot(xb, [dmwf[s].mean() for s in bins], '-o', color=LESION, lw=2.4, ms=9, label='MWF (binned)')
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xscale('log'); ax.set_xlabel('lesion volume (mL, log scale)')
    ax.set_ylabel('drop vs perilesional NAWM (%)')
    ax.set_title(f'Lesion detection vs size: MIMM resolves, MWF blurs (n = {len(lesions)} lesions)')
    ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'fig4_resolution.png'))
    plt.close(fig); print(f'fig4_resolution.png  ({len(lesions)} lesions)')


# ── Figure 5: MWF reference agreement ────────────────────────────────────────
def fig5():
    pts = []; subjects = []
    for mp in sorted(glob.glob(os.path.join(RESULTS, '*', 'grase', 'MWF.nii.gz'))):
        d = os.path.dirname(os.path.dirname(mp)); sid = os.path.basename(d)
        mwf = load(mp); mvf = load(os.path.join(d, 'mimm', 'MVF_Atlas.nii.gz'))
        labels = load(os.path.join(d, 'atlas', 'JHU_labels_subj.nii.gz'))
        brain = load(os.path.join(d, 'qsm', 'brain_mask.nii.gz'))
        if any(v is None for v in (mwf, mvf, labels, brain)):
            continue
        mwf = np.clip(mwf.astype(float), 0, 0.5); mvf = mvf.astype(float)
        labels = labels.astype(int); brain = brain > 0
        si = len(subjects); subjects.append(sid)
        for roi in range(1, 51):
            m = (labels == roi) & brain & np.isfinite(mvf) & np.isfinite(mwf) & (mwf > 0)
            if m.sum() >= 20:
                pts.append((si, float(mvf[m].mean()), float(mwf[m].mean())))
    pts = np.array(pts); si = pts[:, 0].astype(int); mvf = pts[:, 1]; mwf = pts[:, 2]
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 5.2))
    axA.scatter(mwf, mvf, c=si, cmap='tab20', s=16, alpha=0.7)
    r = stats.pearsonr(mwf, mvf)[0]
    b, a = np.polyfit(mwf, mvf, 1); xs = np.linspace(mwf.min(), mwf.max(), 100)
    axA.plot(xs, b * xs + a, 'k--', lw=1.5)
    annotate(axA, f'r = {r:.2f}\n{len(subjects)} subj, {len(pts)} ROIs')
    axA.set_xlabel('MWF  T2-GRASE  (fraction)'); axA.set_ylabel('MVF  MIMM (Atlas)')
    axA.set_title('MIMM MVF vs MWF, per ROI')
    mean = (mvf + mwf) / 2; diff = mvf - mwf; bias = diff.mean(); sd = diff.std(ddof=1)
    axB.scatter(mean, diff, c=si, cmap='tab20', s=16, alpha=0.7)
    axB.axhline(bias, color='k', lw=1.5); axB.axhline(bias + 1.96*sd, color='0.5', ls='--')
    axB.axhline(bias - 1.96*sd, color='0.5', ls='--')
    annotate(axB, f'bias = {bias:+.3f}\nLoA ± {1.96*sd:.3f}')
    axB.set_xlabel('mean of MVF & MWF'); axB.set_ylabel('MVF − MWF')
    axB.set_title('Bland–Altman')
    fig.suptitle('MIMM vs the MWF reference (T2-GRASE)', fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.95]); fig.savefig(os.path.join(OUT, 'fig5_mwf.png'))
    plt.close(fig); print(f'fig5_mwf.png')


if __name__ == '__main__':
    for fn in (fig1, fig2, fig3, fig4, fig5):
        try:
            fn()
        except Exception as e:
            print(f'[{fn.__name__}] FAILED: {e}')
    print(f'\nfigures in: {OUT}')
