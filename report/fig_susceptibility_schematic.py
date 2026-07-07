#!/usr/bin/env python3
"""
fig_susceptibility_schematic.py -- conceptual schematic: myelin and iron shift
the local field in opposite directions, chi-separation splits them apart.

This is illustrative (not tied to one real voxel's numbers -- that's
fig_qsm_match.py, the "QSM half of the match" backup slide). Labels are placed
OUTSIDE the bar ends so they can never be clipped by the bar itself, regardless
of font size or bar length.

Output: report/fig_susceptibility_schematic.png (+ .svg)
"""
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

mpl.rcParams.update({'font.size': 12, 'font.family': 'DejaVu Sans', 'savefig.dpi': 200})
TEAL = '#12A4A9'; NAVY = '#001F3E'; AMBER = '#E9A63C'
OUT = '/home/tijn-saes/Documents/Internship/MIMM/report'

# illustrative susceptibility ranges (ppm) -- schematic, not one subject's values
myelin_lo, myelin_hi = -0.055, -0.005   # diamagnetic: negative, weakens the field
iron_lo, iron_hi = 0.005, 0.032         # paramagnetic: positive, strengthens the field

fig, ax = plt.subplots(figsize=(6.6, 3.6))

# bars as horizontal spans (barh with left offset), NOT anchored at 0, so each
# bar's true extent is visible and the label sits just past its right edge
ax.barh(1.0, myelin_hi - myelin_lo, left=myelin_lo, height=0.55,
        color=TEAL, edgecolor=NAVY, lw=0.8, alpha=0.92, zorder=3)
ax.barh(0.0, iron_hi - iron_lo, left=iron_lo, height=0.55,
        color=AMBER, edgecolor=NAVY, lw=0.8, alpha=0.92, zorder=3)

# labels OUTSIDE the bar (never clipped, regardless of text length/font size)
ax.text(myelin_lo - 0.004, 1.0, 'myelin  $\\chi^{-}$  (diamagnetic)',
        ha='right', va='center', color=NAVY, fontsize=11, fontweight='bold')
ax.text(iron_hi + 0.004, 0.0, 'iron  $\\chi^{+}$  (paramagnetic)',
        ha='left', va='center', color=NAVY, fontsize=11, fontweight='bold')

ax.axvline(0, color=NAVY, lw=1.1, zorder=2)

# weakens/strengthens annotation above the axis -- arrows start well clear of 0
# and the two labels are centred over their own arrow, not the midpoint, so
# they never collide. No separate "0" label: the axis line + tick already mark it.
ytop = 1.62
ax.annotate('', xy=(-0.070, ytop), xytext=(-0.010, ytop),
            arrowprops=dict(arrowstyle='-|>', color='0.4', lw=1.2))
ax.annotate('', xy=(0.070, ytop), xytext=(0.010, ytop),
            arrowprops=dict(arrowstyle='-|>', color='0.4', lw=1.2))
ax.text(-0.040, ytop + 0.10, 'weakens the field', ha='center', va='bottom',
        color='0.35', fontsize=9.5)
ax.text(0.040, ytop + 0.10, 'strengthens the field', ha='center', va='bottom',
        color='0.35', fontsize=9.5)

ax.set_xlim(-0.075, 0.075)
ax.set_ylim(-0.6, 2.15)
ax.set_yticks([])
ax.set_xlabel('susceptibility (ppm)')
for s in ('top', 'right', 'left'):
    ax.spines[s].set_visible(False)
ax.tick_params(left=False)
ax.set_title('Myelin and iron shift the field in opposite directions,\n'
              'and χ-separation splits them', fontsize=12.5, color=NAVY, pad=10)

fig.tight_layout()
fig.savefig(f'{OUT}/fig_susceptibility_schematic.png', bbox_inches='tight', facecolor='white')
fig.savefig(f'{OUT}/fig_susceptibility_schematic.svg', bbox_inches='tight', facecolor='white')
print('saved report/fig_susceptibility_schematic.png + .svg')
