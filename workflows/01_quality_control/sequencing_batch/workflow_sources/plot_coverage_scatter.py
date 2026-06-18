#!/usr/bin/env python3
"""
Scatter plot of total sequencing coverage per individual, grouped by species.

Each point is one individual.  Points are spread deterministically within each
species column (sorted by family then individual) to avoid overplotting.
Marker shape encodes sociality (circle = social, triangle = subsocial).
Colour encodes species (same palette as the heatmap).

Output: sequencing_coverage_scatter.pdf
"""

import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ------------------------------------------------------------------
# Constants (must match plot_heatmap.py)
# ------------------------------------------------------------------
SPECIES_ORDER = [
    'S_dumicola', 'S_mimosarum', 'S_sarasinorum',   # social first
    'S_africanus', 'S_bicolor', 'S_lineatus', 'S_tentoriicola',
]

SPECIES_LABELS = {
    'S_africanus':    'S. africanus',
    'S_bicolor':      'S. bicolor',
    'S_dumicola':     'S. dumicola',
    'S_lineatus':     'S. lineatus',
    'S_mimosarum':    'S. mimosarum',
    'S_sarasinorum':  'S. sarasinorum',
    'S_tentoriicola': 'S. tentoriicola',
}

SPECIES_COLORS = {
    'S_africanus':    '#E41A1C',
    'S_bicolor':      '#FF7F00',
    'S_dumicola':     '#4DAF4A',
    'S_lineatus':     '#377EB8',
    'S_mimosarum':    '#984EA3',
    'S_sarasinorum':  '#A65628',
    'S_tentoriicola': '#F781BF',
}

SOCIAL_SPECIES  = {'S_dumicola', 'S_mimosarum', 'S_sarasinorum'}
SOCIALITY_MARKER = {'Social': 'o', 'Subsocial': '^'}
SOCIALITY_COLORS = {'Social': '#D73027', 'Subsocial': '#4575B4'}

JITTER_WIDTH = 0.35   # half-width of the spread within each species column


def main():
    parser = argparse.ArgumentParser(
        description='Plot total coverage per individual as a scatter / strip chart.')
    parser.add_argument('--input',  required=True,
                        help='Path to sequencing_summary_with_coverage.tsv')
    parser.add_argument('--output', required=True,
                        help='Output PDF path')
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Load TOTAL rows only
    # ------------------------------------------------------------------
    df = pd.read_csv(args.input, sep='\t')
    df = df[df['round_id'] == 'TOTAL'].copy()
    df['coverage'] = pd.to_numeric(df['coverage'], errors='coerce')
    df['sociality'] = df['species'].apply(
        lambda s: 'Social' if s in SOCIAL_SPECIES else 'Subsocial')

    # Sort within each species: family → individual (deterministic x positions)
    df = df.sort_values(['species', 'family', 'individual']).reset_index(drop=True)

    # ------------------------------------------------------------------
    # Assign x positions: species column centre ± spread
    # ------------------------------------------------------------------
    species_present = [s for s in SPECIES_ORDER if s in df['species'].values]
    sp_x_center = {sp: i + 1 for i, sp in enumerate(species_present)}

    x_coords = []
    for sp, grp in df.groupby('species', sort=False):
        n = len(grp)
        if n == 1:
            offsets = [0.0]
        else:
            offsets = np.linspace(-JITTER_WIDTH, JITTER_WIDTH, n)
        x_coords.extend([sp_x_center[sp] + o for o in offsets])

    df['x'] = x_coords

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(14, 7))
    # Reserve the right 25 % of the canvas for the two outside legends so they
    # are fully within the figure boundary and never clipped.
    fig.subplots_adjust(right=0.75)

    for sp in species_present:
        sub = df[df['species'] == sp]
        for _, row in sub.iterrows():
            marker = SOCIALITY_MARKER[row['sociality']]
            ax.scatter(row['x'], row['coverage'],
                       color=SPECIES_COLORS[sp],
                       marker=marker,
                       s=55, alpha=0.85, linewidths=0.4,
                       edgecolors='white', zorder=3)

    # ------------------------------------------------------------------
    # Species column separators and x-axis labels
    # ------------------------------------------------------------------
    n_sp = len(species_present)
    for i in range(n_sp - 1):
        ax.axvline(i + 1.5, color='#cccccc', linewidth=0.8, linestyle='--', zorder=1)

    # Sociality background shading
    social_xs   = [sp_x_center[s] for s in species_present if s in SOCIAL_SPECIES]
    subsocial_xs = [sp_x_center[s] for s in species_present if s not in SOCIAL_SPECIES]
    for xs, col in [(social_xs, '#fdecea'), (subsocial_xs, '#eaf0fb')]:
        for x in xs:
            ax.axvspan(x - 0.5, x + 0.5, color=col, alpha=0.4, zorder=0)

    ax.set_xticks(range(1, n_sp + 1))
    ax.set_xticklabels(
        [SPECIES_LABELS[s] for s in species_present],
        rotation=30, ha='right', fontsize=11, style='italic')
    ax.set_xlim(0.4, n_sp + 0.6)

    ax.set_ylabel('Total sequencing coverage (×)', fontsize=12)
    ax.set_xlabel('')
    ax.tick_params(axis='y', labelsize=10)

    # Light horizontal guide lines
    ymax = df['coverage'].max()
    for ref in [10, 20, 30, 40, 50]:
        if ref < ymax * 1.15:
            ax.axhline(ref, color='#dddddd', linewidth=0.8, linestyle='-', zorder=1)
            # Label inside the axes (left edge) so it never touches the right legend
            ax.text(0.55, ref, f'{ref}×',
                    va='center', ha='left', fontsize=8, color='#aaaaaa', zorder=2)

    ax.set_ylim(0, ymax * 1.1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # ------------------------------------------------------------------
    # Legend
    # ------------------------------------------------------------------
    species_handles = [
        mpatches.Patch(color=SPECIES_COLORS[s],
                       label=SPECIES_LABELS[s])
        for s in species_present
    ]
    sociality_handles = [
        mlines.Line2D([], [], color='#555555',
                      marker=SOCIALITY_MARKER[cat], linestyle='None',
                      markersize=8, label=cat)
        for cat in ['Social', 'Subsocial']
    ]

    # Place both legends outside the axes to the right.
    # Store both handles and pass them to bbox_extra_artists in savefig so
    # that tight-layout correctly expands the canvas to include full text.
    leg1 = ax.legend(handles=species_handles,
                     title='Species', title_fontsize=10,
                     fontsize=9,
                     bbox_to_anchor=(1.02, 1.0), loc='upper left',
                     borderaxespad=0, framealpha=0.9)
    ax.add_artist(leg1)
    leg2 = ax.legend(handles=sociality_handles,
                     title='Sociality', title_fontsize=10,
                     fontsize=9,
                     bbox_to_anchor=(1.02, 0.22), loc='upper left',
                     borderaxespad=0, framealpha=0.9)

    # ------------------------------------------------------------------
    # Title
    # ------------------------------------------------------------------
    ax.set_title('Total sequencing coverage per individual',
                 fontsize=14, fontweight='bold', pad=12)
    n_ind = len(df)
    ax.text(0.5, 1.01,
            f'{n_ind} individuals across {n_sp} species  ·  '
            f'circle = social, triangle = subsocial',
            transform=ax.transAxes, ha='center', va='bottom',
            fontsize=9, color='#555555')

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    fig.savefig(args.output, bbox_inches='tight',
                bbox_extra_artists=(leg1, leg2), dpi=150)
    plt.close(fig)
    print(f'Scatter plot saved to {args.output}')


if __name__ == '__main__':
    main()
