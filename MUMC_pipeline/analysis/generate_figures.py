"""
Generate presentation-ready figures for all MIMM pipeline outputs.
Saves PNG files to subj_dir/figures/.
"""

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.stats import pearsonr
import os

_od = os.environ.get('MIMM_OUTPUT_DIR')
if _od:
    subj_dir = _od
    out_dir  = os.path.join(_od, 'figures')
else:
    try:
        from paths import OUTPUT_DIR as subj_dir, FIG_DIR as out_dir
    except ImportError:
        raise SystemExit('Copy MUMC_pipeline/analysis/paths_template.py to paths.py and fill in your paths.')
os.makedirs(out_dir, exist_ok=True)

# QC-only mode (set MIMM_QC_ONLY=1): produce just the per-subject sanity-check
# figures for a cohort run, and skip the full analysis set. The statistical /
# cohort figures are built later from each subject's roi_stats.csv.
QC_FIGS = {'01_anatomy.png', '02_QSM.png', '04_MVF_basic.png',
           '05_MVF_atlas.png', '12_theta_atlas.png', '13_error.png'}
QC_ONLY = os.environ.get('MIMM_QC_ONLY') == '1'

def load(path):
    return np.array(nib.load(path).dataobj).astype(np.float32)

def get_slices(vol, mask=None):
    """Return central axial, coronal, sagittal slices, masked."""
    if mask is not None:
        coords = np.argwhere(mask > 0)
        cx, cy, cz = coords.mean(axis=0).astype(int)
    else:
        cx, cy, cz = [s // 2 for s in vol.shape[:3]]
    ax  = np.rot90(vol[:, :, cz])
    cor = np.rot90(vol[:, cy, :])
    sag = np.rot90(vol[cx, :, :])
    return ax, cor, sag

def make_bland_altman(vol1, label1, vol2, label2, title, unit, mask, fname):
    """Bland-Altman plot restricted to white matter (FA > 0.2), axes clipped to 99th percentile."""
    wm_mask = mask & (fa > 0.2)
    v1 = vol1[wm_mask].ravel()
    v2 = vol2[wm_mask].ravel()

    mean = (v1 + v2) / 2
    diff = v1 - v2
    bias = np.mean(diff)
    sd   = np.std(diff)
    loa_upper = bias + 1.96 * sd
    loa_lower = bias - 1.96 * sd

    # Clip axes to 1st–99th percentile to remove outlier tail
    xlim = (0, np.percentile(mean, 99))
    ylim = (np.percentile(diff, 1), np.percentile(diff, 99))
    # Add 20% padding around y
    ypad = (ylim[1] - ylim[0]) * 0.2
    ylim = (ylim[0] - ypad, ylim[1] + ypad)

    # Only plot voxels within the clipped range
    keep = (mean >= xlim[0]) & (mean <= xlim[1]) & (diff >= ylim[0]) & (diff <= ylim[1])

    fig, ax = plt.subplots(figsize=(8, 6), facecolor='black')
    ax.set_facecolor('#111111')

    ax.hist2d(mean[keep], diff[keep], bins=150, cmap='hot', cmin=1)

    ax.axhline(bias,      color='cyan',   linewidth=1.8, linestyle='--', label=f'Bias = {bias:.3f} {unit}')
    ax.axhline(loa_upper, color='yellow', linewidth=1.4, linestyle=':',  label=f'+1.96 SD = {loa_upper:.3f} {unit}')
    ax.axhline(loa_lower, color='yellow', linewidth=1.4, linestyle=':',  label=f'−1.96 SD = {loa_lower:.3f} {unit}')
    ax.axhline(0,         color='white',  linewidth=0.8, linestyle='-',  alpha=0.3)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(f'Mean ({label1}, {label2})  [{unit}]', color='white', fontsize=11)
    ax.set_ylabel(f'Difference ({label1} − {label2})  [{unit}]', color='white', fontsize=11)
    ax.set_title(title, color='white', fontsize=13, fontweight='bold')
    ax.tick_params(colors='white')
    for spine in ax.spines.values():
        spine.set_edgecolor('white')
    ax.legend(fontsize=10, facecolor='#222222', edgecolor='white', labelcolor='white')

    plt.tight_layout()
    path = os.path.join(out_dir, fname)
    plt.savefig(path, dpi=150, bbox_inches='tight', facecolor='black', edgecolor='none')
    plt.close()
    print(f'Saved: {fname}')

def make_comparison_figure(vol1, label1, vol2, label2, title, cmap, vmin, vmax, unit, mask, fname, bg=None):
    """Two maps side by side, each shown as axial/coronal/sagittal."""
    ax1, cor1, sag1 = get_slices(vol1, mask)
    ax2, cor2, sag2 = get_slices(vol2, mask)

    fig, axes = plt.subplots(2, 3, figsize=(14, 9), facecolor='black')
    fig.suptitle(title, color='white', fontsize=15, fontweight='bold', y=1.01)

    view_labels = ['Axial', 'Coronal', 'Sagittal']
    rows = [(label1, [ax1, cor1, sag1]), (label2, [ax2, cor2, sag2])]

    for r, (row_label, slices) in enumerate(rows):
        for c, (sl, vlab) in enumerate(zip(slices, view_labels)):
            a = axes[r, c]
            a.set_facecolor('black')
            if bg is not None:
                bg_slices = get_slices(bg, mask)
                bg_sl = bg_slices[c]
                bg_norm = (bg_sl - bg_sl.min()) / (bg_sl.max() - bg_sl.min() + 1e-8)
                a.imshow(bg_norm, cmap='gray', aspect='auto', interpolation='nearest')
                a.imshow(np.ma.masked_where(sl == 0, sl),
                         cmap=cmap, vmin=vmin, vmax=vmax,
                         aspect='auto', interpolation='nearest', alpha=0.85)
            else:
                a.imshow(sl, cmap=cmap, vmin=vmin, vmax=vmax,
                         aspect='auto', interpolation='nearest')
            if r == 0:
                a.set_title(vlab, color='white', fontsize=11)
            a.axis('off')
        axes[r, 0].set_ylabel(row_label, color='white', fontsize=12, rotation=90, labelpad=8)
        axes[r, 0].yaxis.set_label_position('left')
        axes[r, 0].axis('off')

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, orientation='vertical',
                        fraction=0.02, pad=0.02, shrink=0.85)
    cbar.set_label(unit, color='white', fontsize=11)
    cbar.ax.yaxis.set_tick_params(color='white')
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color='white')
    cbar.outline.set_edgecolor('white')

    plt.tight_layout()
    path = os.path.join(out_dir, fname)
    plt.savefig(path, dpi=150, bbox_inches='tight',
                facecolor='black', edgecolor='none')
    plt.close()
    print(f'Saved: {fname}')

def make_figure(vol, title, cmap, vmin, vmax, unit, mask, fname, bg=None):
    if QC_ONLY and fname not in QC_FIGS:
        return
    ax, cor, sag = get_slices(vol, mask)

    fig, axes = plt.subplots(1, 3, figsize=(14, 5), facecolor='black')
    fig.suptitle(title, color='white', fontsize=16, fontweight='bold', y=1.01)

    labels = ['Axial', 'Coronal', 'Sagittal']
    slices = [ax, cor, sag]

    for i, (a, sl, lab) in enumerate(zip(axes, slices, labels)):
        a.set_facecolor('black')
        if bg is not None:
            bg_ax, bg_cor, bg_sag = get_slices(bg, mask)
            bg_sl = [bg_ax, bg_cor, bg_sag][i]
            bg_norm = (bg_sl - bg_sl.min()) / (bg_sl.max() - bg_sl.min() + 1e-8)
            a.imshow(bg_norm, cmap='gray', aspect='auto', interpolation='nearest')
            a.imshow(np.ma.masked_where(sl == 0, sl),
                     cmap=cmap, vmin=vmin, vmax=vmax,
                     aspect='auto', interpolation='nearest', alpha=0.85)
        else:
            im = a.imshow(sl, cmap=cmap, vmin=vmin, vmax=vmax,
                          aspect='auto', interpolation='nearest')
        a.set_title(lab, color='white', fontsize=11)
        a.axis('off')

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, orientation='vertical',
                        fraction=0.02, pad=0.02, shrink=0.85)
    cbar.set_label(unit, color='white', fontsize=11)
    cbar.ax.yaxis.set_tick_params(color='white')
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color='white')
    cbar.outline.set_edgecolor('white')

    plt.tight_layout()
    path = os.path.join(out_dir, fname)
    plt.savefig(path, dpi=150, bbox_inches='tight',
                facecolor='black', edgecolor='none')
    plt.close()
    print(f'Saved: {fname}')

# --- Load data ---
print('Loading data...')
brain  = load(os.path.join(subj_dir, 'qsm',   'brain_mask.nii.gz')) > 0
mag    = load(os.path.join(subj_dir, 'qsm',   'mag_e1.nii.gz'))
qsm    = load(os.path.join(subj_dir, 'qsm',   'QSM.nii.gz'))
r2s    = load(os.path.join(subj_dir, 'qsm',   'R2star.nii.gz'))
mvf_b  = load(os.path.join(subj_dir, 'mimm',  'MVF_basic.nii.gz'))
mvf_a  = load(os.path.join(subj_dir, 'mimm',  'MVF_Atlas.nii.gz'))
fvf    = load(os.path.join(subj_dir, 'mimm',  'FVF_basic.nii.gz'))
gratio = load(os.path.join(subj_dir, 'mimm',  'g_ratio_basic.nii.gz'))
chim   = load(os.path.join(subj_dir, 'mimm',  'chi_myelin_basic.nii.gz'))
chii   = load(os.path.join(subj_dir, 'mimm',  'chi_iron_est_basic.nii.gz'))
err    = load(os.path.join(subj_dir, 'mimm',  'error_basic.nii.gz'))
fa     = load(os.path.join(subj_dir, 'atlas', 'FA_atlas.nii.gz'))
theta  = load(os.path.join(subj_dir, 'atlas', 'theta_atlas.nii.gz'))

# Chi-separation outputs
chisep_neg = load(os.path.join(subj_dir, 'chisep', 'chi_neg.nii.gz'))   # diamagnetic (myelin)
chisep_pos = load(os.path.join(subj_dir, 'chisep', 'chi_pos.nii.gz'))   # paramagnetic (iron)

# Mask values outside brain
for v in [qsm, r2s, mvf_b, mvf_a, fvf, gratio, chim, chii, err, fa, theta,
          chisep_neg, chisep_pos]:
    v[~brain] = 0

# White matter mask: brain AND FA_atlas > 0.20 (excludes GM for microstructure maps)
# Population-average FA is low; 0.20 gives ~180k voxels, primarily WM
wm = brain & (fa > 0.20)

# Both MIMM chi_myelin and chi-sep chi_neg are diamagnetic (negative ppm); take abs throughout
chim_abs       = np.abs(chim)
chisep_neg_abs = np.abs(chisep_neg)

print('Generating figures...')

# 1. Anatomy
make_figure(mag, 'Anatomy — Echo 1 Magnitude', 'gray',
            np.percentile(mag[brain], 1), np.percentile(mag[brain], 99),
            'a.u.', brain, '01_anatomy.png')

# 2. QSM
make_figure(qsm, 'QSM — Susceptibility Map', 'RdBu_r',
            -0.10, 0.10, 'ppm', brain, '02_QSM.png', bg=mag)

# 3. R2*
make_figure(r2s, 'R2* Map', 'hot',
            0, 35, 's⁻¹', brain, '03_R2star.png', bg=mag)

# 4. MVF Basic — WM mask to exclude GM false positives
make_figure(mvf_b, 'MVF — Basic (no orientation prior)', 'hot',
            0, 0.50, 'fraction', wm, '04_MVF_basic.png', bg=mag)

# 5. MVF Atlas — WM mask
make_figure(mvf_a, 'MVF — Atlas (HCP1065 orientation prior)', 'hot',
            0, 0.50, 'fraction', wm, '05_MVF_atlas.png', bg=mag)

# 6. MVF Basic vs Atlas difference — WM mask
diff = mvf_b - mvf_a
diff[~wm] = 0
make_figure(diff, 'MVF Difference — Basic minus Atlas', 'RdBu_r',
            -0.15, 0.15, 'fraction', wm, '06_MVF_diff.png', bg=mag)

# 7. FVF — WM mask
make_figure(fvf, 'FVF — Fibre Volume Fraction', 'viridis',
            0, 0.80, 'fraction', wm, '07_FVF.png', bg=mag)

# 8. g-ratio — WM mask
make_figure(gratio, 'g-ratio', 'plasma',
            0.5, 1.0, '-', wm, '08_g_ratio.png', bg=mag)

# 9. chi_myelin — WM mask
make_figure(chim, 'χ myelin (diamagnetic)', 'Blues_r',
            -0.055, 0, 'ppm', wm, '09_chi_myelin.png', bg=mag)

# 10. chi_iron — WM mask
make_figure(chii, 'χ iron (paramagnetic)', 'Reds',
            0, 0.10, 'ppm', wm, '10_chi_iron.png', bg=mag)

# 11. FA atlas — tighten scale to show actual WM FA range
make_figure(fa, 'FA Atlas (HCP1065) in Subject Space', 'hot',
            0, 0.5, '-', brain, '11_FA_atlas.png', bg=mag)

# 12. Theta atlas — WM mask
make_figure(theta, 'Fibre Angle θ relative to B0', 'hsv',
            0, 90, 'degrees', wm, '12_theta_atlas.png', bg=mag)

# 13. Error map
make_figure(err, 'MIMM Dictionary Matching Error', 'hot',
            0, 0.01, '-', brain, '13_error.png', bg=mag)

# In QC-only mode, stop here: the comparison / Bland-Altman / consistency figures
# (14-17, 30) and the ROI analyses are built later, at cohort level.
if QC_ONLY:
    print('QC-only mode: 6 sanity figures done (01, 02, 04, 05, 12, 13).')
    raise SystemExit(0)

# 14. Chi myelin comparison: MIMM vs chi-separation
make_comparison_figure(
    chim_abs,       'MIMM |χ myelin|',
    chisep_neg_abs, 'χ-sep |χ⁻|',
    'Diamagnetic Susceptibility — MIMM vs χ-separation',
    'Blues_r', 0, 0.15, 'ppm', brain, '14_chi_myelin_comparison.png', bg=mag)

# 15. Chi iron comparison: MIMM vs chi-separation
make_comparison_figure(
    chii,       'MIMM χ iron (paramagnetic)',
    chisep_pos, 'χ-sep χ_pos (paramagnetic)',
    'Paramagnetic Susceptibility — MIMM vs χ-separation',
    'Reds', 0, 0.15, 'ppm', brain, '15_chi_iron_comparison.png', bg=mag)

# 16. Bland-Altman: chi iron — MIMM vs chi-sep
make_bland_altman(chii, 'MIMM χ iron', chisep_pos, 'χ-sep χ_pos',
                  'Bland-Altman — Paramagnetic Susceptibility (iron)',
                  'ppm', brain, '16_BA_chi_iron.png')

# 17. Bland-Altman: chi myelin — MIMM vs chi-sep
make_bland_altman(chim_abs, 'MIMM |χ myelin|', chisep_neg_abs, 'χ-sep |χ⁻|',
                  'Bland-Altman — Diamagnetic Susceptibility (myelin)',
                  'ppm', brain, '17_BA_chi_myelin.png')

# 30. Internal consistency check: MIMM chi_total (= chi_myelin + chi_iron) vs QSM
# MIMM softly matches the sum of its two susceptibility components to the measured
# QSM (chi_error term in MIMM.m, weighted by lambda_chi=0.015). If chi_total
# systematically underestimates QSM, the model under-accounts for total tissue
# susceptibility — the likely root of the diamagnetic/paramagnetic offsets in
# figures 16-17. Whole-brain mask (matching is brain-wide, not WM-only).
# chim is stored negative (diamagnetic), chii positive (paramagnetic): signed sum.
chi_total = (chii + chim)

v_qsm = qsm[brain].ravel()
v_tot = chi_total[brain].ravel()
finite = np.isfinite(v_qsm) & np.isfinite(v_tot)
v_qsm, v_tot = v_qsm[finite], v_tot[finite]

diff = v_tot - v_qsm
bias = np.mean(diff)
sd   = np.std(diff)
loa_u, loa_l = bias + 1.96 * sd, bias - 1.96 * sd
r_tot, _ = pearsonr(v_qsm, v_tot)
slope, intercept = np.polyfit(v_qsm, v_tot, 1)

fig, axes = plt.subplots(1, 2, figsize=(15, 6), facecolor='black')
fig.suptitle('MIMM χ_total (χ_myelin + χ_iron) vs Measured QSM — Internal Consistency',
             color='white', fontsize=15, fontweight='bold', y=1.02)

# Panel A: Bland-Altman
axA = axes[0]
axA.set_facecolor('#111111')
m_xy = (v_tot + v_qsm) / 2
xlimA = (np.percentile(m_xy, 1), np.percentile(m_xy, 99))
ylimA = (np.percentile(diff, 1), np.percentile(diff, 99))
ypad = (ylimA[1] - ylimA[0]) * 0.2
ylimA = (ylimA[0] - ypad, ylimA[1] + ypad)
keep = ((m_xy >= xlimA[0]) & (m_xy <= xlimA[1]) &
        (diff >= ylimA[0]) & (diff <= ylimA[1]))
axA.hist2d(m_xy[keep], diff[keep], bins=150, cmap='hot', cmin=1)
axA.axhline(bias,  color='cyan',   lw=1.8, ls='--', label=f'Bias = {bias:.3f} ppm')
axA.axhline(loa_u, color='yellow', lw=1.4, ls=':',  label=f'+1.96 SD = {loa_u:.3f} ppm')
axA.axhline(loa_l, color='yellow', lw=1.4, ls=':',  label=f'−1.96 SD = {loa_l:.3f} ppm')
axA.axhline(0,     color='white',  lw=0.8, alpha=0.3)
axA.set_xlim(xlimA); axA.set_ylim(ylimA)
axA.set_xlabel('Mean (χ_total, QSM)  [ppm]', color='white', fontsize=11)
axA.set_ylabel('χ_total − QSM  [ppm]', color='white', fontsize=11)
axA.set_title('Bland-Altman', color='white', fontsize=12)
axA.tick_params(colors='white')
for s in axA.spines.values(): s.set_edgecolor('white')
axA.legend(fontsize=9, facecolor='#222222', edgecolor='white', labelcolor='white')

# Panel B: density scatter with identity + regression
axB = axes[1]
axB.set_facecolor('#111111')
lim = (np.percentile(v_qsm, 1), np.percentile(v_qsm, 99))
keepB = (v_qsm >= lim[0]) & (v_qsm <= lim[1]) & (v_tot >= lim[0]) & (v_tot <= lim[1])
axB.hist2d(v_qsm[keepB], v_tot[keepB], bins=150, cmap='hot', cmin=1)
axB.plot(lim, lim, '-', color='#888888', lw=1.2, label='Identity (y = x)')
x_fit = np.linspace(lim[0], lim[1], 100)
axB.plot(x_fit, slope * x_fit + intercept, '--', color='cyan', lw=1.5,
         label=f'Fit: slope={slope:.2f}, r={r_tot:.2f}')
axB.set_xlim(lim); axB.set_ylim(lim)
axB.set_xlabel('Measured QSM  [ppm]', color='white', fontsize=11)
axB.set_ylabel('MIMM χ_total  [ppm]', color='white', fontsize=11)
axB.set_title('Scatter (slope < 1 = χ_total compresses QSM range)', color='white', fontsize=12)
axB.tick_params(colors='white')
for s in axB.spines.values(): s.set_edgecolor('white')
axB.legend(fontsize=9, facecolor='#222222', edgecolor='white', labelcolor='white')

plt.tight_layout()
plt.savefig(os.path.join(out_dir, '30_chi_total_vs_QSM.png'),
            dpi=150, bbox_inches='tight', facecolor='black', edgecolor='none')
plt.close()
print('Saved: 30_chi_total_vs_QSM.png')
print(f'  chi_total vs QSM: bias={bias:.4f} ppm, slope={slope:.3f}, r={r_tot:.3f}')

print(f'\nAll figures saved to: {out_dir}')
