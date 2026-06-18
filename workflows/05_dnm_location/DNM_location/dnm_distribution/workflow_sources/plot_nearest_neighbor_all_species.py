#!/usr/bin/env python3
"""
plot_nearest_neighbor_all_species.py

Multi-panel nearest-neighbour permutation result plot across species.

Inputs are the per-species outputs from nearest_neighbor_test.py:
  - {sp}_nearest_neighbor_simulations.tsv
  - {sp}_nearest_neighbor_summary.txt
"""
import argparse
import math
import os
import re

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


LAYOUT = [
    ('DUM', 'dum', 'TEN', 'ten'),
    ('SAR', 'sar', 'BIC', 'bic'),
    ('MIM', 'mim', 'AFR', 'afr'),
    (None, None, 'LIN', 'lin'),
]


def parse_number(value):
    """Parse a summary value that may contain commas, units, or text."""
    match = re.search(r'[-+]?(?:\d+(?:,\d{3})*|\d+)(?:\.\d+)?', value)
    if not match:
        raise ValueError(f"Could not parse numeric value from: {value!r}")
    return float(match.group(0).replace(',', ''))


def parse_summary(summary_path):
    """Return key nearest-neighbour statistics from a summary text file."""
    fields = {}
    with open(summary_path) as fh:
        for line in fh:
            if ':' not in line:
                continue
            key, value = line.split(':', 1)
            fields[key.strip()] = value.strip()

    return {
        'observed_bp': parse_number(fields['Observed ANND']),
        'sim_mean_bp': parse_number(fields['Simulated ANND mean']),
        'z_score': parse_number(fields['Z-score']),
        'p_two_tailed': parse_number(fields['p (two-tailed)']),
        'n_dnms': int(parse_number(fields['Unique autosomal DNM events'])),
        'n_simulations': int(parse_number(fields['N simulations'])),
    }


def result_paths(input_dir, SP, sp):
    species_dir = os.path.join(input_dir, SP)
    return (
        os.path.join(species_dir, f"{sp}_nearest_neighbor_simulations.tsv"),
        os.path.join(species_dir, f"{sp}_nearest_neighbor_summary.txt"),
    )


def plot_species(ax, input_dir, SP, sp, column_label):
    sim_path, summary_path = result_paths(input_dir, SP, sp)
    if not os.path.exists(sim_path):
        raise FileNotFoundError(sim_path)
    if not os.path.exists(summary_path):
        raise FileNotFoundError(summary_path)

    sim_df = pd.read_csv(sim_path, sep='\t')
    stats = parse_summary(summary_path)
    sim_mb = sim_df['avg_nearest_neighbour_dist_bp'] / 1e6

    ax.hist(
        sim_mb, bins=50,
        color='steelblue', alpha=0.75, edgecolor='white', linewidth=0.4,
        label=f"Simulated (n={stats['n_simulations']})",
    )
    ax.axvline(
        stats['observed_bp'] / 1e6,
        color='crimson', linewidth=2,
        label=f"Observed ({stats['observed_bp'] / 1e6:.2f} Mb)",
    )
    ax.axvline(
        stats['sim_mean_bp'] / 1e6,
        color='steelblue', linewidth=1.5, linestyle='--',
        label=f"Sim. mean ({stats['sim_mean_bp'] / 1e6:.2f} Mb)",
    )

    ax.set_title(
        f"{SP} ({column_label})\n"
        f"p = {stats['p_two_tailed']:.4f} | Z = {stats['z_score']:.2f} | "
        f"DNMs = {stats['n_dnms']}",
        fontsize=10, fontweight='bold',
    )
    ax.tick_params(axis='both', labelsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(fontsize=7, frameon=False)


def main():
    parser = argparse.ArgumentParser(
        description='Plot all species nearest-neighbour test results.'
    )
    parser.add_argument('--input_dir', required=True,
                        help='Directory containing nearest_neighbor_test/{SP} outputs')
    parser.add_argument('--output', required=True,
                        help='Output PDF path')
    args = parser.parse_args()

    n_rows = len(LAYOUT)
    fig, axes = plt.subplots(
        n_rows, 2,
        figsize=(12, max(10, 2.6 * n_rows)),
        sharex=False,
        sharey=False,
    )

    for row_i, (social_SP, social_sp, subsocial_SP, subsocial_sp) in enumerate(LAYOUT):
        left_ax, right_ax = axes[row_i]

        if social_SP is None:
            left_ax.axis('off')
        else:
            plot_species(left_ax, args.input_dir, social_SP, social_sp, 'social')

        plot_species(right_ax, args.input_dir, subsocial_SP, subsocial_sp, 'subsocial')

    fig.supxlabel('Average nearest-neighbour distance (Mb)', fontsize=11)
    fig.supylabel('Number of simulations', fontsize=11)
    fig.suptitle(
        'DNM spatial distribution test across species',
        fontsize=13, fontweight='bold',
    )

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    plt.tight_layout(rect=[0.03, 0.03, 1, 0.96])
    plt.savefig(args.output, bbox_inches='tight')
    plt.close()
    print(f"Saved: {args.output}")


if __name__ == '__main__':
    main()
