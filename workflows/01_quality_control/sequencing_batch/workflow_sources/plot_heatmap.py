#!/usr/bin/env python3
"""
Visualise the sequencing scheme as a coverage heatmap.

All individuals were sequenced on a single platform (DNBSEQ-G400); the columns
are flowcell x lane combinations on that platform (flowcell IDs are run
identifiers, not different platforms).

Rows    : individuals  (labelled "individual  [species / family]")
Columns : unique (flowcell, lane) combinations, labelled "flowcell·lane"

Cell colour encodes the approximate per-lane sequencing coverage (from the
'coverage' column produced by add_coverage.py).  Cells where the individual
was not sequenced in that lane are shown in light grey.

Left annotation strips (left to right):
  1. Role bar       – parent vs. proband (offspring)
  2. Sociality bar  – social (DUM, MIM, SAR) vs. subsocial (AFR, BIC, LIN, TEN)
  3. Species bar    – one colour per species

A continuous colour bar for coverage is placed to the right of the heatmap.

Output: sequencing_heatmap.pdf  (vector PDF, suitable for publication)
"""

import argparse
import os
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.colorbar as mcolorbar
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ------------------------------------------------------------------
# Colour dictionaries
# ------------------------------------------------------------------
SPECIES_COLORS = {
    'S_africanus':    '#E41A1C',
    'S_bicolor':      '#FF7F00',
    'S_dumicola':     '#4DAF4A',
    'S_lineatus':     '#377EB8',
    'S_mimosarum':    '#984EA3',
    'S_sarasinorum':  '#A65628',
    'S_tentoriicola': '#F781BF',
}

SOCIAL_SPECIES = {'S_dumicola', 'S_mimosarum', 'S_sarasinorum'}

SOCIALITY_COLORS = {
    'Social':    '#D73027',   # warm red
    'Subsocial': '#4575B4',   # cool blue
}

ROLE_COLORS = {
    'parent':  '#1b9e77',   # teal-green
    'proband': '#7570b3',   # muted purple
}

NO_DATA_COLOR = '#e8e8e8'   # light grey for unsequenced cells


def species_to_sociality(species):
    return 'Social' if species in SOCIAL_SPECIES else 'Subsocial'


def role_of(individual):
    """Parent vs. proband.  Probands are offspring labelled '...S<number>' or,
    in S. bicolor families 3-5, with a bare numeric suffix '..._<number>'
    (e.g. '150.9_1'); parents end in a sex token (Ma/Fema/M/F/male/female)."""
    return 'proband' if re.search(r'(?:S\d+|_\d+)$', individual) else 'parent'


def main():
    parser = argparse.ArgumentParser(
        description='Plot sequencing-scheme coverage heatmap.')
    parser.add_argument('--input',  required=True,
                        help='Path to sequencing_summary_with_coverage.tsv')
    parser.add_argument('--output', required=True,
                        help='Output PDF path')
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Load data – keep only per-lane rows (drop TOTAL)
    # ------------------------------------------------------------------
    df = pd.read_csv(args.input, sep='\t')
    df = df[df['round_id'] != 'TOTAL'].copy()
    df['coverage'] = pd.to_numeric(df['coverage'], errors='coerce').fillna(0.0)

    # Composite column key
    df['lane_key'] = df['round_id'] + '·' + df['lane']

    # Row label
    df['ind_label'] = (df['individual'] + '  ['
                       + df['species'].str.replace('S_', '', regex=False)
                       + ' / ' + df['family'] + ']')

    # Sort rows: sociality → species → family → individual
    df['sociality'] = df['species'].map(species_to_sociality)
    df['role'] = df['individual'].map(role_of)
    df_sorted = df.sort_values(['sociality', 'species', 'family', 'individual'])

    row_order = (df_sorted[['ind_label', 'species', 'family', 'individual',
                            'sociality', 'role']]
                 .drop_duplicates())
    col_order = sorted(df['lane_key'].unique())

    # Build coverage pivot table (0 = not sequenced in that lane)
    cov_pivot = (
        df.pivot_table(index='ind_label', columns='lane_key',
                       values='coverage', aggfunc='sum')
        .reindex(index=row_order['ind_label'], columns=col_order, fill_value=0.0)
    )

    n_rows, n_cols = cov_pivot.shape

    # ------------------------------------------------------------------
    # Colormap: mask 0-coverage cells so they render as NO_DATA_COLOR
    # ------------------------------------------------------------------
    data = cov_pivot.values.astype(float)
    masked = np.ma.masked_where(data == 0, data)

    cmap = plt.cm.YlOrRd.copy()
    cmap.set_bad(NO_DATA_COLOR)

    vmax = np.nanpercentile(data[data > 0], 98) if (data > 0).any() else 1.0

    # ------------------------------------------------------------------
    # Figure layout
    # ------------------------------------------------------------------
    col_width        = 0.18   # inches per column
    row_height       = 0.18   # inches per row
    top_margin       = 3.5    # room for column labels
    bottom_margin    = 0.8
    right_pad        = 0.8

    # Distance from figure left edge to the first annotation strip.
    # Labels extend LEFTWARD from the tick spine, so this just needs to be
    # wide enough to fit the longest label (~40 chars × ~0.036 in ≈ 1.5 in).
    role_left        = 3.0    # inches: figure left → role strip (leftmost) left edge
    label_strip_gap  = 0.15   # gap between label text and the strip
    ax_labels_w      = 0.01   # thin invisible reference axes
    bar_width        = 0.25   # width of each annotation strip
    strip_gap        = 0.06   # gap between adjacent strips
    strip_heat_gap   = 0.1    # gap between last strip and heatmap
    cbar_gap         = 0.15
    cbar_w           = 0.18

    soc_left          = role_left + bar_width + strip_gap
    sp_left           = soc_left + bar_width + strip_gap
    heatmap_left_inch = sp_left + bar_width + strip_heat_gap
    heatmap_w_inch    = n_cols * col_width
    heatmap_h_inch    = n_rows * row_height

    fig_w = heatmap_left_inch + heatmap_w_inch + cbar_gap + cbar_w + right_pad
    fig_h = top_margin + heatmap_h_inch + bottom_margin

    fig = plt.figure(figsize=(fig_w, fig_h))

    def to_frac(x_inch, y_inch, w_inch, h_inch):
        return [x_inch / fig_w, y_inch / fig_h,
                w_inch / fig_w, h_inch / fig_h]

    # Row-label axes: a near-zero-width spine positioned label_strip_gap inches
    # to the LEFT of ax_soc.  With default LEFT-side ticks, labels extend
    # LEFTWARD from this spine into the white margin – never touching the strips.
    ax_labels_x = role_left - label_strip_gap - ax_labels_w
    ax_labels = fig.add_axes(to_frac(ax_labels_x, bottom_margin,
                                      ax_labels_w, heatmap_h_inch))
    ax_labels.set_xlim(0, 1)
    ax_labels.set_ylim(n_rows - 0.5, -0.5)
    ax_labels.set_xticks([])
    ax_labels.set_yticks(range(n_rows))
    ax_labels.set_yticklabels(row_order['ind_label'], fontsize=6.5)
    # DEFAULT left-side ticks: labels extend leftward (away from strips)
    ax_labels.tick_params(axis='y', which='major', length=0)
    for spine in ax_labels.spines.values():
        spine.set_visible(False)
    ax_labels.set_facecolor('none')

    # Role strip (leftmost)
    ax_role = fig.add_axes(to_frac(role_left, bottom_margin,
                                    bar_width, heatmap_h_inch))

    # Sociality strip
    ax_soc = fig.add_axes(to_frac(soc_left, bottom_margin,
                                   bar_width, heatmap_h_inch))

    # Species strip
    ax_sp = fig.add_axes(to_frac(sp_left, bottom_margin,
                                  bar_width, heatmap_h_inch))

    # Heatmap axes
    ax = fig.add_axes(to_frac(heatmap_left_inch, bottom_margin,
                               heatmap_w_inch, heatmap_h_inch))

    # Coverage colorbar
    cb_left = heatmap_left_inch + heatmap_w_inch + cbar_gap
    ax_cbar = fig.add_axes(to_frac(cb_left, bottom_margin + heatmap_h_inch * 0.1,
                                    cbar_w, heatmap_h_inch * 0.8))

    # ------------------------------------------------------------------
    # Draw heatmap with pcolormesh (clean polygon edges, no PDF bleed)
    # ------------------------------------------------------------------
    X_edges = np.arange(n_cols + 1) - 0.5   # cell edges at -0.5, 0.5, …, n_cols-0.5
    Y_edges = np.arange(n_rows + 1) - 0.5

    ax.pcolormesh(X_edges, Y_edges, masked,
                  cmap=cmap, vmin=0, vmax=vmax,
                  edgecolors='#aaaaaa', linewidth=0.35,
                  antialiased=False)

    # Column labels on top
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(col_order, rotation=90, fontsize=6, ha='right')
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    ax.tick_params(axis='x', which='major', length=0)

    # No row labels on the heatmap axes — handled by ax_labels
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([])
    ax.tick_params(axis='y', which='major', length=0)

    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)   # top-down orientation

    # ------------------------------------------------------------------
    # Species annotation strip  –  Rectangle patches: one solid polygon
    # per row, no interpolation possible, guaranteed crisp boundaries.
    # ------------------------------------------------------------------
    for i, sp in enumerate(row_order['species'].values):
        ax_sp.add_patch(Rectangle(
            xy=(0, i - 0.5), width=1, height=1,
            facecolor=SPECIES_COLORS.get(sp, '#999999'),
            edgecolor='none', lw=0))
    ax_sp.set_xlim(0, 1)
    ax_sp.set_ylim(n_rows - 0.5, -0.5)
    ax_sp.set_xticks([])
    ax_sp.set_yticks([])
    ax_sp.set_title('Sp.', fontsize=8, pad=2)
    for spine in ax_sp.spines.values():
        spine.set_visible(False)

    # ------------------------------------------------------------------
    # Sociality annotation strip  –  same Rectangle approach
    # ------------------------------------------------------------------
    for i, soc in enumerate(row_order['sociality'].values):
        ax_soc.add_patch(Rectangle(
            xy=(0, i - 0.5), width=1, height=1,
            facecolor=SOCIALITY_COLORS[soc],
            edgecolor='none', lw=0))
    ax_soc.set_xlim(0, 1)
    ax_soc.set_ylim(n_rows - 0.5, -0.5)
    ax_soc.set_xticks([])
    ax_soc.set_yticks([])
    ax_soc.set_title('Soc.', fontsize=8, pad=2)
    for spine in ax_soc.spines.values():
        spine.set_visible(False)

    # ------------------------------------------------------------------
    # Role annotation strip  –  parent vs. proband, same Rectangle approach
    # ------------------------------------------------------------------
    for i, ro in enumerate(row_order['role'].values):
        ax_role.add_patch(Rectangle(
            xy=(0, i - 0.5), width=1, height=1,
            facecolor=ROLE_COLORS[ro],
            edgecolor='none', lw=0))
    ax_role.set_xlim(0, 1)
    ax_role.set_ylim(n_rows - 0.5, -0.5)
    ax_role.set_xticks([])
    ax_role.set_yticks([])
    ax_role.set_title('Role', fontsize=8, pad=2)
    for spine in ax_role.spines.values():
        spine.set_visible(False)

    # ------------------------------------------------------------------
    # Coverage colorbar
    # ------------------------------------------------------------------
    norm = mcolors.Normalize(vmin=0, vmax=vmax)
    cb = mcolorbar.ColorbarBase(ax_cbar, cmap=cmap, norm=norm,
                                orientation='vertical')
    cb.set_label('Coverage (×)', fontsize=10)
    cb.ax.tick_params(labelsize=9)

    # Patch for no-data in colorbar legend
    no_data_patch = mpatches.Patch(color=NO_DATA_COLOR, label='Not sequenced')

    # ------------------------------------------------------------------
    # Legends
    # ------------------------------------------------------------------
    species_handles = [
        mpatches.Patch(color=col, label=sp.replace('S_', 'S. '))
        for sp, col in SPECIES_COLORS.items()
    ]
    sociality_handles = [
        mpatches.Patch(color=col, label=cat)
        for cat, col in SOCIALITY_COLORS.items()
    ] + [no_data_patch]
    role_handles = [
        mpatches.Patch(color=col, label=cat.capitalize())
        for cat, col in ROLE_COLORS.items()
    ]

    # Single combined legend placed BELOW the figure body (y=0 is the figure
    # bottom edge).  loc='upper center' means the top of the legend box sits
    # at that anchor, so the entire box is outside the figure.  bbox_inches='tight'
    # in savefig will expand the canvas to include it.
    spacer = mpatches.Patch(visible=False, label='')   # blank divider
    all_handles = (species_handles + [spacer] + sociality_handles
                   + [spacer] + role_handles)
    fig.legend(handles=all_handles,
               loc='upper center',
               bbox_to_anchor=(0.5, 0.0),
               ncol=len(all_handles),
               fontsize=9, framealpha=0.9,
               handlelength=1.2, handletextpad=0.5, columnspacing=1.0)

    # ------------------------------------------------------------------
    # Title
    # ------------------------------------------------------------------
    fig.text(0.5, 0.995,
             'Sequencing coverage per individual × lane',
             ha='center', va='top', fontsize=15, fontweight='bold')
    fig.text(0.5, 0.978,
             f'single platform: DNBSEQ-G400  ·  {n_rows} individuals  ·  '
             f'{n_cols} unique (flowcell × lane) combinations  ·  '
             f'colour = approx. coverage (×)',
             ha='center', va='top', fontsize=10, color='#555555')

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    fig.savefig(args.output, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f'Heatmap saved to {args.output}')


if __name__ == '__main__':
    main()
