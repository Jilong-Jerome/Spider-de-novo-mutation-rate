#!/usr/bin/env python3
"""
plot_spectrum.py
Render a standard SBS96 bar plot: 6 mutation-type panels, 16 trinucleotide
contexts each, canonical colors.
"""
import argparse
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

MUT_TYPES = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
# Canonical SBS96 colors (as used by COSMIC / SigProfiler)
MUT_COLORS = {
    'C>A': '#03BCEE',
    'C>G': '#010101',
    'C>T': '#E32926',
    'T>A': '#CAC9C9',
    'T>C': '#A1CE63',
    'T>G': '#EBC6C4',
}
BASES = ['A', 'C', 'G', 'T']


def read_spectrum(path):
    data = {}
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            row = dict(zip(header, fields))
            data[row['sbs96_category']] = int(row['count'])
    return data


def main():
    parser = argparse.ArgumentParser(description='Plot SBS96 spectrum')
    parser.add_argument('--input', required=True, help='Spectrum TSV (from compute_spectrum.py)')
    parser.add_argument('--output', required=True, help='Output PDF/PNG')
    parser.add_argument('--title', default='', help='Plot title')
    parser.add_argument('--normalize', action='store_true',
                        help='Plot proportions instead of raw counts')
    args = parser.parse_args()

    counts = read_spectrum(args.input)
    total = sum(counts.values())
    if args.normalize and total > 0:
        values = {k: v / total for k, v in counts.items()}
        ylabel = 'Proportion'
    else:
        values = counts
        ylabel = 'Count'

    # Build canonical context order: for each mut_type, 16 contexts ordered by (5', 3')
    fig, axes = plt.subplots(1, 6, figsize=(18, 4), sharey=True)
    for i, mut in enumerate(MUT_TYPES):
        ax = axes[i]
        cats = []
        vals = []
        labels = []
        for five in BASES:
            for three in BASES:
                cat = f'{five}[{mut}]{three}'
                cats.append(cat)
                vals.append(values.get(cat, 0))
                labels.append(f'{five}{mut[0]}{three}')
        ax.bar(range(16), vals, color=MUT_COLORS[mut], edgecolor='none', width=0.8)
        ax.set_title(mut, fontsize=12, color=MUT_COLORS[mut] if mut != 'C>G' else 'black')
        ax.set_xticks(range(16))
        ax.set_xticklabels(labels, rotation=90, fontsize=6, family='monospace')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if i == 0:
            ax.set_ylabel(ylabel)

    title = args.title if args.title else os.path.basename(args.input).replace('.tsv', '')
    fig.suptitle(f'{title}  (n={total})', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    fig.savefig(args.output, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'Plotted {total} DNMs -> {args.output}')


if __name__ == '__main__':
    main()
