"""
mimm_style.py  —  one shared visual language for the MIMM defence deck.

Goal: every figure (results re-renders on the radstation AND the local
schematics) looks like it came from one hand. You do that two ways:
  1) call apply_style() once at the top of each figure script, and
  2) pull colours by MEANING from C[...], never by raw hex.

Usage
-----
    from mimm_style import apply_style, C, C_CHI, COMPARTMENT, FIGSIZE
    apply_style()
    ...
    ax.scatter(mwf, mvf, color=C["MIMM"])                 # MIMM data = teal
    ax.plot([lo, hi], [lo, hi], color=C["reference"], ls="--")  # reference series = grey

The golden rule: one hue = one meaning across the WHOLE deck. Teal is always
MIMM/myelin, amber is always iron, grey is always the reference-method data,
navy is always axes/text/measured. Never reuse a hue for a second meaning.
"""

import matplotlib as mpl

# ---------------------------------------------------------------------------
# Semantic palette — reference these by name, do not scatter hex in scripts.
# ---------------------------------------------------------------------------
C = {
    "MIMM":      "#12A4A9",  # teal  — MIMM / myelin (MVF, MIMM compartments)
    "reference": "#8A94A0",  # grey  — reference METHODS as data (MWF, chi-sep-as-method)
    "iron":      "#E9A63C",  # amber — iron (chi+)
    "secondary": "#C2C8CF",  # light grey — a THIRD, de-emphasised also-ran series
    "highlight": "#E15241",  # coral — the ONE highlighted finding per figure
    "text":      "#001F3E",  # navy  — text, axes, tick labels, measured values
    "grid":      "#E5E7EB",  # very light — gridlines
}
# NOTE: reference is GREY, not navy. navy (#001F3E) is text/axes ONLY, if a
# reference series were navy it would vanish into the axes. Do not "fix" this
# back to navy on the radstation; tell Claude and it changes in one place.

# chi-separation is context-dependent:
#   * as a WHOLE METHOD compared to MIMM  -> C["reference"] (grey)
#   * shown as its DECOMPOSED CHANNELS    -> use C_CHI below
# Lock this so the agreement figure and the offset figure never disagree.
C_CHI = {
    "chi_myelin": "#12A4A9",  # chi- -> teal  (myelin family)
    "chi_iron":   "#E9A63C",  # chi+ -> amber (iron)
}

# MIMM compartments are ALL MIMM, so they stay in the teal family and are
# separated by SHADE, never by borrowing amber/grey. Highlight the key bar
# in-script with C["highlight"] (see helper below).
COMPARTMENT = {
    "FVF": "#0E7E82",  # dark teal   — fibre (the whole; FVF = MVF + AVF)
    "MVF": "#12A4A9",  # base teal   — myelin
    "AVF": "#7FCBCE",  # light teal  — axon
}

# Fixed geometry so every single-panel figure renders at one scale.
FIGSIZE = (6.4, 4.2)   # inches; keep constant across all result figures
DPI = 200


def apply_style():
    """Apply shared rcParams. Call once at the top of every figure script."""
    mpl.rcParams.update({
        # geometry
        "figure.figsize": FIGSIZE,
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        # type
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 12,
        "axes.titlesize": 14,
        "axes.titleweight": "regular",
        "axes.labelsize": 12.5,
        # axes / spines
        "axes.edgecolor": C["text"],
        "axes.labelcolor": C["text"],
        "axes.linewidth": 1.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        # grid
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": C["grid"],
        "grid.linewidth": 0.6,
        # ticks
        "xtick.color": C["text"],
        "ytick.color": C["text"],
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        # lines / markers
        "lines.linewidth": 2.2,
        "lines.markersize": 8,
        # legend
        "legend.frameon": False,
        "legend.fontsize": 11,
        # misc
        "text.color": C["text"],
    })


# ---------------------------------------------------------------------------
# Small helpers for the two conventions that are easy to get wrong.
# ---------------------------------------------------------------------------
def highlight_bar(ax, bars, index, label=None):
    """Recolour one bar to the highlight colour and (optionally) annotate it.

    Use for 'the one finding per figure' — e.g. the globus-pallidus bar on the
    GP chart, or the significant compartment on the clinical figure.
    """
    b = bars[index]
    b.set_edgecolor(C["highlight"])
    b.set_linewidth(2.5)
    if label:
        ax.annotate(label,
                    xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                    xytext=(0, 10), textcoords="offset points",
                    ha="center", va="bottom",
                    color=C["highlight"], fontweight="bold", fontsize=11)


def dim(artist):
    """De-emphasise a non-key series (the 'also-rans' go light grey)."""
    try:
        artist.set_color(C["secondary"])
    except Exception:
        for a in artist:
            a.set_color(C["secondary"])


if __name__ == "__main__":
    # tiny self-check render so you can eyeball the palette
    import matplotlib.pyplot as plt
    apply_style()
    fig, ax = plt.subplots()
    names = list(C.keys())
    for i, k in enumerate(names):
        ax.barh(i, 1, color=C[k])
        ax.text(1.02, i, k, va="center", color=C["text"])
    ax.set_yticks([]); ax.set_xlim(0, 1.6); ax.set_title("mimm_style palette")
    fig.savefig("mimm_style_swatches.png")
    print("wrote mimm_style_swatches.png")
