#!/usr/bin/env python3
"""
plot_dnm_genome.py

Genome-wide DNM tick-mark visualisation.
One subplot per autosome; tick marks at each DNM position.
"""
import argparse
import re
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def load_fai(fai_path):
    """Return {chrom: length} dict."""
    chrom_lengths = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.strip().split('\t')
            chrom_lengths[parts[0]] = int(parts[1])
    return chrom_lengths


def chrom_sort_key(chrom):
    """Sort chromosomes numerically by extracting trailing integer."""
    m = re.search(r'(\d+)$', chrom)
    return int(m.group(1)) if m else 0


def main():
    parser = argparse.ArgumentParser(
        description='Plot genome-wide DNM tick-mark visualisation.'
    )
    parser.add_argument('--dnm_file',   required=True)
    parser.add_argument('--genome_fai', required=True)
    parser.add_argument('--species',    required=True,
                        help='Species prefix (e.g. afr)')
    parser.add_argument('--x_chroms',   nargs='+', default=[])
    parser.add_argument('--output',     required=True,
                        help='Output PDF path')
    args = parser.parse_args()

    x_chroms_set  = set(args.x_chroms)
    chrom_lengths = load_fai(args.genome_fai)

    # Load DNMs, filter autosomes, deduplicate shared sibling mutations
    # (same chrom/pos/ref/alt/father/mother = one parental event, show one tick)
    dnm_df = pd.read_csv(args.dnm_file, sep='\t')
    dnm_df = dnm_df[~dnm_df['chrom'].isin(x_chroms_set)].copy()
    dnm_df = dnm_df.drop_duplicates(
        subset=['chrom', 'pos', 'ref', 'alt', 'father', 'mother']
    )

    # Autosomal chromosomes present in FAI (sorted numerically)
    autosomes = sorted(
        [c for c in chrom_lengths if c not in x_chroms_set],
        key=chrom_sort_key,
    )
    n_chroms  = len(autosomes)
    total_dnm = len(dnm_df)

    # Uniform x-axis scale: all panels share the length of the longest autosome
    max_chrom_len = max(chrom_lengths[c] for c in autosomes)

    fig, axes = plt.subplots(
        n_chroms, 1,
        figsize=(12, max(3, n_chroms * 0.7)),
        sharex=True,
    )
    if n_chroms == 1:
        axes = [axes]

    fig.suptitle(
        f"{args.species.upper()} — {total_dnm} de novo mutations (autosomes)",
        fontsize=13, fontweight='bold', y=1.0,
    )

    for ax, chrom in zip(axes, autosomes):
        chrom_len = chrom_lengths[chrom]
        # Grey baseline drawn only to actual chromosome end, not full axis width
        ax.hlines(0, 0, chrom_len, colors='lightgrey', linewidth=1.5, zorder=1)

        # DNM ticks
        chrom_dnms = dnm_df[dnm_df['chrom'] == chrom]['pos'].values
        if len(chrom_dnms) > 0:
            ax.vlines(
                chrom_dnms, -0.4, 0.4,
                colors='steelblue', linewidth=0.8, alpha=0.7, zorder=2,
            )

        ax.set_xlim(0, max_chrom_len)
        ax.set_ylim(-0.6, 0.6)
        ax.set_yticks([])
        ax.set_ylabel(chrom, rotation=0, ha='right', va='center', fontsize=8)
        ax.tick_params(axis='x', labelsize=7)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)

    # Format x-axis in Mb on the shared axis (bottom panel only shows ticks)
    axes[-1].xaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _: f"{x/1e6:.0f} Mb")
    )

    plt.tight_layout()
    plt.savefig(args.output, bbox_inches='tight')
    plt.close()
    print(f"Saved: {args.output}")


if __name__ == '__main__':
    main()
