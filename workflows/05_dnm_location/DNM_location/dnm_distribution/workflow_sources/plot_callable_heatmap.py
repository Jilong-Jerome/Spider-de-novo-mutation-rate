#!/usr/bin/env python3
"""
plot_callable_heatmap.py

Callable fraction heatmap: offspring (rows) × chromosome (columns).
"""
import argparse
import re
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def load_fai(fai_path):
    """Return {chrom: length} dict."""
    chrom_lengths = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.strip().split('\t')
            chrom_lengths[parts[0]] = int(parts[1])
    return chrom_lengths


def chrom_sort_key(chrom):
    """Sort chromosomes numerically by trailing integer."""
    m = re.search(r'(\d+)$', chrom)
    return int(m.group(1)) if m else 0


def extract_offspring(offspring_chrom_id):
    """
    offspring_chrom_id format: {offspring}_{chrom}
    where chrom itself may contain underscores (e.g. afr_1).
    Strip the last two underscore-separated tokens to recover offspring name.

    Example:
      "AFR_family1_S1_offspring_afr_1" → "AFR_family1_S1_offspring"
    """
    parts = offspring_chrom_id.rsplit('_', 2)
    # parts[-2] and parts[-1] form the chrom (e.g. "afr", "1")
    return parts[0]


def main():
    parser = argparse.ArgumentParser(
        description='Plot callable-fraction heatmap per individual per chromosome.'
    )
    parser.add_argument('--callable_summary', required=True,
                        help='Callable summary TSV (no header): chrom, callable_sites, offspring_chrom_id')
    parser.add_argument('--genome_fai',       required=True)
    parser.add_argument('--species',          required=True,
                        help='Species prefix (e.g. afr)')
    parser.add_argument('--x_chroms',         nargs='+', default=[])
    parser.add_argument('--output',           required=True,
                        help='Output PDF path')
    args = parser.parse_args()

    x_chroms_set  = set(args.x_chroms)
    chrom_lengths = load_fai(args.genome_fai)

    # Load callable summary (no header)
    df = pd.read_csv(
        args.callable_summary, sep='\t', header=None,
        names=['chrom', 'callable_sites', 'offspring_chrom_id'],
    )

    # Filter autosomes
    df = df[~df['chrom'].isin(x_chroms_set)].copy()

    # Extract offspring name
    df['offspring'] = df['offspring_chrom_id'].apply(extract_offspring)

    # Callable fraction = callable_sites / chrom_length
    df['chrom_length']       = df['chrom'].map(chrom_lengths)
    df['callable_fraction']  = df['callable_sites'] / df['chrom_length']

    # Pivot: rows = offspring, cols = chrom
    pivot = df.pivot_table(
        index='offspring',
        columns='chrom',
        values='callable_fraction',
        aggfunc='mean',   # each (offspring, chrom) should be unique
    )

    # Sort columns numerically, rows lexicographically (family → individual)
    autosomes = sorted(
        [c for c in pivot.columns if c not in x_chroms_set],
        key=chrom_sort_key,
    )
    pivot = pivot[autosomes]
    pivot = pivot.sort_index()

    # Plot
    fig_height = max(4, len(pivot) * 0.35 + 2)
    fig_width  = max(8, len(autosomes) * 0.6 + 2)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    sns.heatmap(
        pivot,
        ax=ax,
        cmap='Blues',
        vmin=0, vmax=1,
        annot=False,
        linewidths=0.2,
        linecolor='lightgrey',
        cbar_kws={'label': 'Callable fraction'},
    )

    ax.set_title(
        f"{args.species.upper()} — callable genome fraction per individual per chromosome",
        fontsize=11, fontweight='bold', pad=12,
    )
    ax.set_xlabel('Chromosome', fontsize=9)
    ax.set_ylabel('Individual', fontsize=9)
    ax.tick_params(axis='x', labelsize=7, rotation=45)
    ax.tick_params(axis='y', labelsize=7)

    plt.tight_layout()
    plt.savefig(args.output, bbox_inches='tight')
    plt.close()
    print(f"Saved: {args.output}")


if __name__ == '__main__':
    main()
