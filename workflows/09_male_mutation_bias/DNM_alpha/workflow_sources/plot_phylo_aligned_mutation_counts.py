#!/usr/bin/env python3
"""
Plot per-species unique and sibling-shared DNM locus counts aligned to a
species cladogram.
"""
import argparse
import csv
import os

import yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import numpy as np

plt.rcParams.update({
    'font.size': 13,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 11,
    'legend.title_fontsize': 12,
})

UNIQUE_COLOR = '#9E9E9E'
SHARED_COLOR = '#111111'
PATERNAL_COLOR = '#1F4E79'
MATERNAL_COLOR = '#5B2C6F'
UNPHASED_COLOR = '#BDBDBD'
SOCIAL_COLOR = '#8B2E2E'
SUBSOCIAL_COLOR = '#2F4F7F'
PANEL_D_SPECIES_ORDER = ['AFR', 'BIC', 'LIN', 'TEN']


def canonical_species(label):
    return label.split(':', 1)[0].strip().upper()


def parse_newick(s):
    s = s.strip().rstrip(';').strip()
    pos = [0]

    def _clade():
        children = []
        if pos[0] < len(s) and s[pos[0]] == '(':
            pos[0] += 1
            children.append(_clade())
            while pos[0] < len(s) and s[pos[0]] == ',':
                pos[0] += 1
                children.append(_clade())
            if pos[0] >= len(s) or s[pos[0]] != ')':
                raise ValueError(f'Expected ) at pos {pos[0]}: {s}')
            pos[0] += 1
        label_start = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in ',()':
            pos[0] += 1
        label = canonical_species(s[label_start:pos[0]].strip())
        return (children, label)

    return _clade()


def assign_layout(tree):
    leaf_counter = [0]
    leaf_order = []

    def walk(clade):
        children, label = clade
        if not children:
            y = leaf_counter[0]
            leaf_counter[0] += 1
            leaf_order.append(label)
            return {'height': 0, 'y': y, 'label': label, 'children': []}
        kids = [walk(c) for c in children]
        y = sum(k['y'] for k in kids) / len(kids)
        height = 1 + max(k['height'] for k in kids)
        return {'height': height, 'y': y, 'label': label, 'children': kids}

    root = walk(tree)
    max_height = root['height']

    def assign_x(node):
        node['x'] = max_height - node['height']
        for child in node['children']:
            assign_x(child)
    assign_x(root)

    segments = []

    def emit(node):
        if not node['children']:
            return
        ys = [child['y'] for child in node['children']]
        segments.append(((node['x'], min(ys)), (node['x'], max(ys))))
        for child in node['children']:
            segments.append(((node['x'], child['y']), (child['x'], child['y'])))
            emit(child)
    emit(root)

    return segments, leaf_order, max_height


def read_summary(path):
    rows = {}
    with open(path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sp = row['species']
            rows[sp] = {
                'species': sp,
                'unique': int(row['n_unique_mutation_loci']),
                'shared': int(row['n_sibling_shared_mutation_loci']),
            }
    return rows


def read_phasing(path):
    rows = {}
    with open(path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sp = row['species']
            sharing_class = row['sharing_class']
            rows[(sp, sharing_class)] = {
                'P': int(row['n_P_loci']),
                'M': int(row['n_M_loci']),
                'U': int(row['n_U_loci']),
            }
    return rows


def read_parental_rates(path):
    parental_rows = []
    social_rows = []
    with open(path, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            record_type = row.get('record_type', 'subsocial_parental')
            if record_type == 'social_reference':
                social_rows.append({
                    'species': row['species'],
                    'phylo_order': int(row['phylo_order']),
                    'social_germline_rate': float(row['social_germline_rate']),
                })
            else:
                parental_rows.append({
                    'species': row['species'],
                    'phylo_order': int(row['phylo_order']),
                    'paternal_rate_observed': float(row['paternal_rate_observed']),
                    'paternal_rate_ci_low': float(row['paternal_rate_ci_low']),
                    'paternal_rate_ci_high': float(row['paternal_rate_ci_high']),
                    'maternal_rate_observed': float(row['maternal_rate_observed']),
                    'maternal_rate_ci_low': float(row['maternal_rate_ci_low']),
                    'maternal_rate_ci_high': float(row['maternal_rate_ci_high']),
                })
    order = {sp: i for i, sp in enumerate(PANEL_D_SPECIES_ORDER)}
    return (
        sorted(parental_rows, key=lambda r: order.get(r['species'], len(order))),
        sorted(social_rows, key=lambda r: r['phylo_order']),
    )


def render(
    sharing_rows, phasing_rows, parental_rate_rows, social_rate_rows,
    segments, leaf_order, max_depth,
    social_set, subsocial_set, output_pdf,
):
    ordered = [
        sharing_rows.get(sp, {'species': sp, 'unique': 0, 'shared': 0})
        for sp in leaf_order
    ]
    n = len(ordered)
    y_pos = np.arange(n)

    fig = plt.figure(figsize=(9, 9))
    outer = GridSpec(
        2, 3, figure=fig,
        width_ratios=[1, 1.5, 1.5],
        height_ratios=[1.0, 0.8],
        wspace=0.30,
        hspace=0.72,
    )
    ax_tree = fig.add_subplot(outer[0, 0])
    ax_counts = fig.add_subplot(outer[0, 1], sharey=ax_tree)
    phase_grid = GridSpecFromSubplotSpec(
        1, 2, subplot_spec=outer[0, 2], width_ratios=[2, 1], wspace=0.24,
    )
    ax_unique_phase = fig.add_subplot(phase_grid[0, 0], sharey=ax_tree)
    ax_shared_phase = fig.add_subplot(phase_grid[0, 1], sharey=ax_tree)
    parental_grid = GridSpecFromSubplotSpec(
        1, 2, subplot_spec=outer[1, :], width_ratios=[1, 1], wspace=0.18,
    )
    ax_maternal_rate = fig.add_subplot(parental_grid[0, 0])
    ax_paternal_rate = fig.add_subplot(parental_grid[0, 1], sharey=ax_maternal_rate)

    for (x1, y1), (x2, y2) in segments:
        ax_tree.plot([x1, x2], [y1, y2], color='black', linewidth=1.2)

    tip_x = max_depth + 0.15
    for i, row in enumerate(ordered):
        sp = row['species']
        if sp in social_set:
            color = SOCIAL_COLOR
        elif sp in subsocial_set:
            color = SUBSOCIAL_COLOR
        else:
            color = '#888888'
        ax_tree.plot(
            tip_x, i, marker='s', markersize=8, color=color,
            markeredgecolor='black', markeredgewidth=0.4,
        )
        ax_tree.text(tip_x + 0.28, i, sp, ha='left', va='center')

    ax_tree.set_xlim(-0.2, max_depth + 1.6)
    for spine in ax_tree.spines.values():
        spine.set_visible(False)
    ax_tree.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    ax_tree.text(
        0.0, 1.05, 'A', transform=ax_tree.transAxes,
        ha='left', va='bottom', fontweight='bold',
    )

    lifestyle_handles = [
        plt.Line2D(
            [0], [0], marker='s', linestyle='', color=SOCIAL_COLOR,
            markeredgecolor='black', markeredgewidth=0.4, markersize=8,
            label='social',
        ),
        plt.Line2D(
            [0], [0], marker='s', linestyle='', color=SUBSOCIAL_COLOR,
            markeredgecolor='black', markeredgewidth=0.4, markersize=8,
            label='subsocial',
        ),
    ]
    unique = np.array([row['unique'] for row in ordered], dtype=float)
    shared = np.array([row['shared'] for row in ordered], dtype=float)
    bar_h = 0.48
    unique_bars = ax_counts.barh(
        y_pos, unique, height=bar_h, color=UNIQUE_COLOR,
        alpha=1.0, label='unique',
    )
    shared_bars = ax_counts.barh(
        y_pos, shared, left=unique, height=bar_h, color=SHARED_COLOR,
        alpha=1.0, label='sibling-shared',
    )
    ax_counts.grid(axis='x', linestyle=':', alpha=0.4)
    ax_counts.tick_params(left=False, labelleft=False)
    ax_counts.text(
        0.0, 1.05, 'B', transform=ax_counts.transAxes,
        ha='left', va='bottom', fontweight='bold',
    )

    phase_specs = [
        ('P', PATERNAL_COLOR, 'P'),
        ('M', MATERNAL_COLOR, 'M'),
        ('U', UNPHASED_COLOR, 'U'),
    ]
    for ax, sharing_class, title in [
        (ax_unique_phase, 'unique', 'unique'),
        (ax_shared_phase, 'sibling_shared', 'sibling-shared'),
    ]:
        left = np.zeros(n, dtype=float)
        for phase, color, label in phase_specs:
            vals = np.array([
                phasing_rows.get((row['species'], sharing_class), {}).get(phase, 0)
                for row in ordered
            ], dtype=float)
            ax.barh(
                y_pos, vals, left=left, height=bar_h,
                color=color, alpha=1.0, label=label,
            )
            left += vals
        ax.set_title(title)
        ax.grid(axis='x', linestyle=':', alpha=0.4)
        ax.tick_params(left=False, labelleft=False)

    ax_unique_phase.text(
        0.0, 1.05, 'C', transform=ax_unique_phase.transAxes,
        ha='left', va='bottom', fontweight='bold',
    )
    handles, labels = ax_unique_phase.get_legend_handles_labels()
    phase_handles = handles
    phase_labels = labels

    ax_tree.set_ylim(n - 0.4, -0.6)

    x = np.arange(len(parental_rate_rows))
    paternal = np.array([r['paternal_rate_observed'] for r in parental_rate_rows])
    paternal_lo = np.array([r['paternal_rate_ci_low'] for r in parental_rate_rows])
    paternal_hi = np.array([r['paternal_rate_ci_high'] for r in parental_rate_rows])
    maternal = np.array([r['maternal_rate_observed'] for r in parental_rate_rows])
    maternal_lo = np.array([r['maternal_rate_ci_low'] for r in parental_rate_rows])
    maternal_hi = np.array([r['maternal_rate_ci_high'] for r in parental_rate_rows])
    ax_paternal_rate.vlines(
        x, paternal_lo, paternal_hi, color=PATERNAL_COLOR,
        linewidth=1.2, alpha=0.9,
    )
    ax_paternal_rate.plot(
        x, paternal, 'o', color=PATERNAL_COLOR, markersize=5,
    )
    ax_maternal_rate.vlines(
        x, maternal_lo, maternal_hi, color=MATERNAL_COLOR,
        linewidth=1.2, alpha=0.9,
    )
    ax_maternal_rate.plot(
        x, maternal, 'o', color=MATERNAL_COLOR, markersize=5,
    )

    ymax = max(
        np.nanmax(paternal_hi) if paternal_hi.size else 0,
        np.nanmax(maternal_hi) if maternal_hi.size else 0,
        max([r['social_germline_rate'] for r in social_rate_rows] or [0]),
    )
    social_line_styles = [':', '--', '-.']
    for ax, title in [
        (ax_maternal_rate, 'maternal'),
        (ax_paternal_rate, 'paternal'),
    ]:
        for i, row in enumerate(social_rate_rows):
            ax.axhline(
                row['social_germline_rate'], color=SOCIAL_COLOR,
                linestyle=social_line_styles[i % len(social_line_styles)],
                linewidth=0.9, alpha=0.8,
                label=f"{row['species']} social" if ax is ax_maternal_rate else None,
            )
        ax.set_xticks(x)
        ax.set_xticklabels([r['species'] for r in parental_rate_rows])
        ax.set_title(title)
        ax.grid(axis='y', linestyle=':', alpha=0.4)
        if ymax > 0:
            ax.set_ylim(0, ymax * 1.12)
    ax_maternal_rate.set_ylabel('mutation rate')
    handles, labels = ax_maternal_rate.get_legend_handles_labels()
    ax_paternal_rate.legend(
        handles, labels, loc='center left', bbox_to_anchor=(1.03, 0.5),
        frameon=False, fontsize=9,
    )
    ax_paternal_rate.tick_params(labelleft=False)
    ax_maternal_rate.text(
        0.0, 1.05, 'D', transform=ax_maternal_rate.transAxes,
        ha='left', va='bottom', fontweight='bold',
    )

    fig.subplots_adjust(top=0.91, bottom=0.09, left=0.07, right=0.84)
    counts_pos = ax_counts.get_position()
    tree_pos = ax_tree.get_position()
    unique_phase_pos = ax_unique_phase.get_position()
    shared_phase_pos = ax_shared_phase.get_position()
    top_row_label_y = counts_pos.y0 - 0.050
    top_row_legend_y = counts_pos.y0 - 0.090
    compact_legend_kwargs = {
        'frameon': False,
        'borderaxespad': 0.0,
        'handletextpad': 0.35,
        'columnspacing': 0.8,
        'labelspacing': 0.25,
        'handlelength': 1.2,
    }
    fig.text(
        (counts_pos.x0 + counts_pos.x1) / 2, top_row_label_y,
        'number of mutations', ha='center', va='center',
    )
    fig.text(
        (unique_phase_pos.x0 + shared_phase_pos.x1) / 2, top_row_label_y,
        'number of mutations', ha='center', va='center',
    )
    fig.legend(
        handles=lifestyle_handles,
        loc='upper center',
        bbox_to_anchor=((tree_pos.x0 + tree_pos.x1) / 2, top_row_legend_y),
        ncol=1, **compact_legend_kwargs,
    )
    fig.legend(
        handles=[unique_bars[0], shared_bars[0]],
        labels=['unique', 'sibling-shared'],
        loc='upper center',
        bbox_to_anchor=((counts_pos.x0 + counts_pos.x1) / 2, top_row_legend_y),
        ncol=1, **compact_legend_kwargs,
    )
    fig.legend(
        handles=phase_handles, labels=phase_labels,
        loc='upper center',
        bbox_to_anchor=((unique_phase_pos.x0 + shared_phase_pos.x1) / 2, top_row_legend_y),
        ncol=3, **compact_legend_kwargs,
    )
    os.makedirs(os.path.dirname(output_pdf), exist_ok=True)
    fig.savefig(output_pdf)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', required=True)
    ap.add_argument('--summary_tsv', required=True)
    ap.add_argument('--phasing_tsv', required=True)
    ap.add_argument('--parental_rates_tsv', required=True)
    ap.add_argument('--output_pdf', required=True)
    args = ap.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    groups = cfg.get('spectrum_groups', {})
    social_set = set(groups.get('social', {}).get('species', []))
    subsocial_set = set(groups.get('subsocial', {}).get('species', []))

    tree = parse_newick(cfg['species_tree_newick'])
    segments, leaf_order, max_depth = assign_layout(tree)
    sharing_rows = read_summary(args.summary_tsv)
    phasing_rows = read_phasing(args.phasing_tsv)
    parental_rate_rows, social_rate_rows = read_parental_rates(args.parental_rates_tsv)

    render(
        sharing_rows, phasing_rows, parental_rate_rows, social_rate_rows,
        segments, leaf_order, max_depth,
        social_set, subsocial_set, args.output_pdf,
    )


if __name__ == '__main__':
    main()
