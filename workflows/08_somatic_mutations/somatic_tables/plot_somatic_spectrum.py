#!/usr/bin/env python3
"""
plot_somatic_spectrum.py
Render the somatic 96SBS spectrum for all species as a single combined figure:
one row per species, two column blocks (absolute count on the left, fraction on
the right). Each cell shows all 96 categories as bars in canonical SBS order,
colored by the 6 mutation types (COSMIC / SigProfiler palette).

Reads somatic_spectrum_counts.tsv produced by somatic_spectrum_table.py.
"""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

MUT_TYPES = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
MUT_COLORS = {
    'C>A': '#03BCEE',
    'C>G': '#010101',
    'C>T': '#E32926',
    'T>A': '#CAC9C9',
    'T>C': '#A1CE63',
    'T>G': '#EBC6C4',
}
BASES = ['A', 'C', 'G', 'T']
SPECIES = ['AFR', 'BIC', 'DUM', 'LIN', 'MIM', 'SAR', 'TEN']

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
IN_TSV = os.path.join(OUT_DIR, 'somatic_spectrum_counts.tsv')
OUT_PDF = os.path.join(OUT_DIR, 'somatic_spectrum_96SBS.pdf')


def build_sbs96_order():
    """Canonical SBS96 order: list of (mutation_type, category, context_label)."""
    cats = []
    for mut in MUT_TYPES:
        ref = mut[0]
        for five in BASES:
            for three in BASES:
                cat = f'{five}[{mut}]{three}'
                cats.append((mut, cat, f'{five}{ref}{three}'))
    return cats


def read_table(path):
    """Return {species: {sbs96_category: count}} and {species: total}."""
    counts = {sp: {} for sp in SPECIES}
    totals = {sp: 0 for sp in SPECIES}
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        i_sp = header.index('species')
        i_cat = header.index('sbs96_category')
        i_cnt = header.index('count')
        for line in f:
            fields = line.rstrip('\n').split('\t')
            sp = fields[i_sp]
            cat = fields[i_cat]
            if cat == 'overall':
                totals[sp] = int(fields[i_cnt])
                continue
            counts[sp][cat] = int(fields[i_cnt])
    return counts, totals


def main():
    sbs96 = build_sbs96_order()
    cats = [c for _, c, _ in sbs96]
    labels = [lab for _, _, lab in sbs96]
    bar_colors = [MUT_COLORS[m] for m, _, _ in sbs96]
    x = range(len(cats))

    counts, totals = read_table(IN_TSV)

    n_rows = len(SPECIES)
    fig, axes = plt.subplots(
        n_rows, 1, figsize=(16, 2.2 * n_rows + 1),
        sharex=True, squeeze=False,
        gridspec_kw={'hspace': 0.25},
    )

    for r, sp in enumerate(SPECIES):
        total = totals[sp]
        raw = [counts[sp].get(c, 0) for c in cats]

        # Left axis: absolute count. Bars drawn once.
        ax = axes[r][0]
        ax.bar(x, raw, color=bar_colors, edgecolor='none', width=0.85)
        ax.spines['top'].set_visible(False)
        ax.margins(x=0.005)
        ax.set_ylabel(f'{sp}\n(n={total})', fontsize=10,
                      rotation=0, ha='right', va='center', labelpad=28)
        if r == 0:
            ax.set_title('Count (left axis)  /  Fraction (right axis)',
                         fontsize=12)

        # Right axis: fraction = count / total, same bars at a rescaled y-axis.
        ax2 = ax.twinx()
        ax2.spines['top'].set_visible(False)
        lo, hi = ax.get_ylim()
        if total > 0:
            ax2.set_ylim(lo / total, hi / total)
        ax2.set_ylabel('Fraction', fontsize=8, rotation=-90, labelpad=12)
        ax2.tick_params(labelsize=7)

        # Context labels only on the bottom row.
        if r == n_rows - 1:
            ax.set_xticks(list(x))
            ax.set_xticklabels(labels, rotation=90, fontsize=4.5,
                               family='monospace')
        else:
            ax.set_xticks([])

    legend_handles = [Patch(facecolor=MUT_COLORS[m], label=m) for m in MUT_TYPES]
    fig.legend(handles=legend_handles, loc='upper center', ncol=6,
               frameon=False, fontsize=11, bbox_to_anchor=(0.5, 1.0))
    fig.suptitle('Somatic 96SBS mutation spectrum per species', y=1.02,
                 fontsize=14)

    fig.savefig(OUT_PDF, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'Wrote {OUT_PDF}')


if __name__ == '__main__':
    main()
