#!/usr/bin/env python3
"""
plot_rates.py
Render rate plots from the master rate TSV and the somatic-vs-germline
fold-change TSV produced by aggregate_rates.py.

Outputs three PDFs:
  - per_species_rates.pdf : faceted by species, x = 8 classes, two series.
  - group_rates.pdf       : faceted by group (all_species/social/subsocial).
  - fold_change.pdf       : log2(somatic/germline) per group x class.
Optionally also renders three folded SBS96 PDFs when 96-context TSVs and
output paths are provided.
"""
import argparse
import os
import math

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

CLASSES = ['overall', 'C>A', 'C>G', 'C>T_nonCpG', 'C>T_CpG', 'T>A', 'T>C', 'T>G']
MUT_TYPES = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
BASES = ['A', 'C', 'G', 'T']
LEVEL_COLORS = {'germline': '#2D7DD2', 'somatic': '#E15554'}
GROUP_COLORS = {'all_species': '#2B2B2B', 'social': '#8B2E2E', 'subsocial': '#2F4F7F'}
MUT_COLORS = {
    'C>A': '#03BCEE',
    'C>G': '#010101',
    'C>T': '#E32926',
    'T>A': '#CAC9C9',
    'T>C': '#A1CE63',
    'T>G': '#EBC6C4',
}


def load_master(path):
    rows = []
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            rows.append(dict(zip(header, fields)))
    return rows


def load_fc(path):
    rows = []
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            row = dict(zip(header, fields))
            rows.append(row)
    return rows


def load_fold_change_tests(path):
    tests = {}
    if not path:
        return tests
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        for line in f:
            if line.startswith('#'):
                continue
            row = dict(zip(header, line.rstrip('\n').split('\t')))
            tests[row['mutation_class']] = row
    return tests


def to_float(s):
    if s is None or s == '':
        return float('nan')
    try:
        return float(s)
    except ValueError:
        return float('nan')


def panel_grouped_bars(ax, classes, germ_data, som_data, title):
    """germ_data/som_data: dict class -> (rate, ci_low, ci_high)."""
    x = np.arange(len(classes))
    w = 0.4
    g_rate = [germ_data.get(c, (np.nan,))[0] for c in classes]
    s_rate = [som_data.get(c, (np.nan,))[0] for c in classes]
    def err(data, c):
        if c not in data:
            return 0.0, 0.0
        rate, lo, hi = data[c]
        if math.isnan(rate):
            return 0.0, 0.0
        low = max(rate - lo, 0.0) if not math.isnan(lo) else 0.0
        high = max(hi - rate, 0.0) if not math.isnan(hi) else 0.0
        return low, high
    g_lo, g_hi = np.array([err(germ_data, c) for c in classes]).T
    s_lo, s_hi = np.array([err(som_data, c) for c in classes]).T
    ax.bar(x - w/2, g_rate, w, yerr=[g_lo, g_hi], capsize=2,
           color=LEVEL_COLORS['germline'], label='germline')
    ax.bar(x + w/2, s_rate, w, yerr=[s_lo, s_hi], capsize=2,
           color=LEVEL_COLORS['somatic'], label='somatic')
    ax.set_yscale('log')
    ax.set_xticks(x)
    ax.set_xticklabels(classes, rotation=45, ha='right', fontsize=7)
    ax.set_ylabel('rate (log)')
    ax.set_title(title, fontsize=9)
    ax.grid(axis='y', linestyle=':', alpha=0.4)


def index_master(rows):
    """Return idx[(group, species, level)] = {class: (rate, ci_low, ci_high)}."""
    idx = {}
    for r in rows:
        key = (r['group'], r['species'], r['level'])
        idx.setdefault(key, {})[r['mutation_class']] = (
            to_float(r['rate']),
            to_float(r['ci_low']),
            to_float(r['ci_high']),
        )
    return idx


def sbs96_order():
    cats = []
    for mut in MUT_TYPES:
        ref = mut[0]
        for five in BASES:
            for three in BASES:
                trinuc = f'{five}{ref}{three}'
                cats.append((mut, trinuc, f'{five}[{mut}]{three}'))
    return cats


def index_master96(rows):
    idx = {}
    for r in rows:
        key = (r['group'], r['species'], r['level'])
        cat = (r['mutation_type'], r['trinuc_context'], r['sbs96_category'])
        idx.setdefault(key, {})[cat] = (
            to_float(r['rate']),
            to_float(r['ci_low']),
            to_float(r['ci_high']),
        )
    return idx


def rate_yerr(data, cat):
    if cat not in data:
        return 0.0, 0.0
    rate, lo, hi = data[cat]
    if math.isnan(rate):
        return 0.0, 0.0
    low = max(rate - lo, 0.0) if not math.isnan(lo) else 0.0
    high = max(hi - rate, 0.0) if not math.isnan(hi) else 0.0
    return low, high


def plot_sbs96_rate_page(fig_title, germ_data, som_data):
    fig, axes = plt.subplots(1, 6, figsize=(22, 4.5), sharey=True)
    for i, mut in enumerate(MUT_TYPES):
        ax = axes[i]
        cats = []
        labels = []
        ref = mut[0]
        for five in BASES:
            for three in BASES:
                trinuc = f'{five}{ref}{three}'
                cats.append((mut, trinuc, f'{five}[{mut}]{three}'))
                labels.append(trinuc)
        x = np.arange(len(cats))
        w = 0.38
        g_rate = [germ_data.get(c, (np.nan,))[0] for c in cats]
        s_rate = [som_data.get(c, (np.nan,))[0] for c in cats]
        g_yerr = np.array([rate_yerr(germ_data, c) for c in cats]).T
        s_yerr = np.array([rate_yerr(som_data, c) for c in cats]).T
        ax.bar(x - w / 2, g_rate, w, yerr=g_yerr, capsize=1.5,
               color=LEVEL_COLORS['germline'], label='germline')
        ax.bar(x + w / 2, s_rate, w, yerr=s_yerr, capsize=1.5,
               color=LEVEL_COLORS['somatic'], label='somatic')
        ax.set_yscale('log')
        ax.set_title(mut, fontsize=11, color=MUT_COLORS[mut] if mut != 'C>G' else 'black')
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=90, fontsize=6, family='monospace')
        ax.grid(axis='y', linestyle=':', alpha=0.4)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if i == 0:
            ax.set_ylabel('rate (log)')
            ax.legend(fontsize=8)
    fig.suptitle(fig_title, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    return fig


def plot_sbs96_rate_pdf(rows96, output, per_species):
    idx = index_master96(rows96)
    if per_species:
        entities = [('per_species', sp, sp)
                    for sp in sorted({r['species'] for r in rows96 if r['group'] == 'per_species'})]
    else:
        entities = [(g, 'merged', g)
                    for g in sorted({r['group'] for r in rows96 if r['group'] != 'per_species'})]
    os.makedirs(os.path.dirname(output), exist_ok=True)
    with PdfPages(output) as pdf:
        for group, species, title in entities:
            fig = plot_sbs96_rate_page(
                title,
                idx.get((group, species, 'germline'), {}),
                idx.get((group, species, 'somatic'), {}),
            )
            pdf.savefig(fig)
            plt.close(fig)


def index_fc96(rows):
    idx = {}
    for r in rows:
        cat = (r['mutation_type'], r['trinuc_context'], r['sbs96_category'])
        idx.setdefault(r['group'], {})[cat] = (
            to_float(r['log2_fold_change']),
            to_float(r['fc_ci_low']),
            to_float(r['fc_ci_high']),
        )
    return idx


def fc_yerr(data, cat):
    if cat not in data:
        return 0.0, 0.0
    l2, lo, hi = data[cat]
    if math.isnan(l2):
        return 0.0, 0.0
    low = l2 - math.log2(lo) if lo and lo > 0 else 0.0
    high = math.log2(hi) - l2 if hi and hi > 0 else 0.0
    return max(low, 0.0), max(high, 0.0)


def plot_sbs96_fold_change_pdf(rows96, output):
    by_group = index_fc96([r for r in rows96 if r['species'] == 'merged'])
    groups = sorted(by_group)
    os.makedirs(os.path.dirname(output), exist_ok=True)
    fig, axes = plt.subplots(1, 6, figsize=(22, 4.5), sharey=True)
    for i, mut in enumerate(MUT_TYPES):
        ax = axes[i]
        cats = []
        labels = []
        ref = mut[0]
        for five in BASES:
            for three in BASES:
                trinuc = f'{five}{ref}{three}'
                cats.append((mut, trinuc, f'{five}[{mut}]{three}'))
                labels.append(trinuc)
        x = np.arange(len(cats))
        w = 0.8 / len(groups) if groups else 0.8
        for j, group in enumerate(groups):
            data = by_group[group]
            vals = [data.get(c, (np.nan,))[0] for c in cats]
            yerr = np.array([fc_yerr(data, c) for c in cats]).T
            ax.bar(x + (j - (len(groups) - 1) / 2) * w, vals, w,
                   yerr=yerr, capsize=1.5, label=group)
        ax.axhline(0, color='k', linewidth=0.6)
        ax.set_title(mut, fontsize=11, color=MUT_COLORS[mut] if mut != 'C>G' else 'black')
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=90, fontsize=6, family='monospace')
        ax.grid(axis='y', linestyle=':', alpha=0.4)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if i == 0:
            ax.set_ylabel('log2(somatic / germline)')
            ax.legend(fontsize=8)
    fig.suptitle('Somatic vs germline rate fold change by SBS96 context', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(output)
    plt.close(fig)


def group_order(groups):
    preferred = ['all_species', 'social', 'subsocial']
    return [g for g in preferred if g in groups] + sorted(g for g in groups if g not in preferred)


def significance_label(p):
    if math.isnan(p):
        return ''
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    return ''


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--master', required=True)
    ap.add_argument('--fold_change', required=True)
    ap.add_argument('--fold_change_test')
    ap.add_argument('--output_per_species', required=True)
    ap.add_argument('--output_group', required=True)
    ap.add_argument('--output_fold_change', required=True)
    ap.add_argument('--master_96')
    ap.add_argument('--fold_change_96')
    ap.add_argument('--output_per_species_96')
    ap.add_argument('--output_group_96')
    ap.add_argument('--output_fold_change_96')
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.output_per_species), exist_ok=True)

    rows = load_master(args.master)
    idx = index_master(rows)

    # Per-species: one panel per species
    species = sorted({r['species'] for r in rows
                      if r['group'] == 'per_species'})
    n = len(species)
    cols = min(3, n)
    rows_n = math.ceil(n / cols) if n else 1
    fig, axes = plt.subplots(rows_n, cols, figsize=(4.5 * cols, 3.5 * rows_n),
                             squeeze=False)
    for i, sp in enumerate(species):
        ax = axes[i // cols][i % cols]
        germ = idx.get(('per_species', sp, 'germline'), {})
        som = idx.get(('per_species', sp, 'somatic'), {})
        panel_grouped_bars(ax, CLASSES, germ, som, sp)
        if i == 0:
            ax.legend(fontsize=8)
    for j in range(n, rows_n * cols):
        axes[j // cols][j % cols].axis('off')
    fig.tight_layout()
    fig.savefig(args.output_per_species)
    plt.close(fig)

    # Group: one panel per group
    groups = sorted({r['group'] for r in rows if r['group'] != 'per_species'})
    cols = min(3, len(groups)) or 1
    rows_n = math.ceil(len(groups) / cols) if groups else 1
    fig, axes = plt.subplots(rows_n, cols, figsize=(4.5 * cols, 3.5 * rows_n),
                             squeeze=False)
    for i, g in enumerate(groups):
        ax = axes[i // cols][i % cols]
        germ = idx.get((g, 'merged', 'germline'), {})
        som = idx.get((g, 'merged', 'somatic'), {})
        panel_grouped_bars(ax, CLASSES, germ, som, g)
        if i == 0:
            ax.legend(fontsize=8)
    for j in range(len(groups), rows_n * cols):
        axes[j // cols][j % cols].axis('off')
    fig.tight_layout()
    fig.savefig(args.output_group)
    plt.close(fig)

    # Fold change: horizontal raw fold-change plot per group x class (group level only)
    fc_rows = [r for r in load_fc(args.fold_change) if r['species'] == 'merged']
    fc_tests = load_fold_change_tests(args.fold_change_test)
    by_group = {}
    for r in fc_rows:
        by_group.setdefault(r['group'], {})[r['mutation_class']] = (
            to_float(r['fold_change']),
            to_float(r['fc_ci_low']),
            to_float(r['fc_ci_high']),
        )
    n_g = len(by_group)
    fig, ax = plt.subplots(figsize=(7, 5))
    if n_g:
        y = np.arange(len(CLASSES))
        h = 0.8 / n_g
        for i, g in enumerate(group_order(by_group)):
            data = by_group[g]
            fc = [data.get(c, (np.nan,))[0] for c in CLASSES]
            lo = [data.get(c, (np.nan, np.nan, np.nan))[1] for c in CLASSES]
            hi = [data.get(c, (np.nan, np.nan, np.nan))[2] for c in CLASSES]
            yerr = [
                [fc[k] - lo[k] if lo[k] and lo[k] > 0 and not math.isnan(fc[k]) else 0
                 for k in range(len(CLASSES))],
                [hi[k] - fc[k] if hi[k] and hi[k] > 0 and not math.isnan(fc[k]) else 0
                 for k in range(len(CLASSES))],
            ]
            ax.barh(y + (i - (n_g - 1) / 2) * h, fc, h, xerr=yerr, capsize=2,
                    label=g, color=GROUP_COLORS.get(g))
        star_x = []
        for k, cls in enumerate(CLASSES):
            max_x = 1.0
            for data in by_group.values():
                fc, _, hi = data.get(cls, (np.nan, np.nan, np.nan))
                if not math.isnan(hi) and hi > 0:
                    max_x = max(max_x, hi)
                elif not math.isnan(fc) and fc > 0:
                    max_x = max(max_x, fc)
            star_x.append(max_x)
            row = fc_tests.get(cls)
            label = significance_label(to_float(row.get('p_adjusted')) if row else float('nan'))
            if label:
                ax.text(max_x * 1.08, k, label, va='center', ha='left',
                        fontsize=10, fontweight='bold')
        if any(fc_tests.get(cls) and significance_label(to_float(fc_tests[cls].get('p_adjusted')))
               for cls in CLASSES):
            ax.set_xlim(right=max(star_x) * 1.25)
        ax.axvline(1, color='0.25', linewidth=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(CLASSES, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel('fold change (somatic / germline)')
        ax.legend(fontsize=8)
        ax.grid(axis='x', linestyle=':', alpha=0.4)
        ax.set_title('Somatic vs germline rate fold change')
    fig.tight_layout()
    fig.savefig(args.output_fold_change)
    plt.close(fig)

    print(f'Wrote {args.output_per_species}')
    print(f'Wrote {args.output_group}')
    print(f'Wrote {args.output_fold_change}')

    args96 = [
        args.master_96,
        args.fold_change_96,
        args.output_per_species_96,
        args.output_group_96,
        args.output_fold_change_96,
    ]
    if any(args96):
        if not all(args96):
            raise SystemExit('All 96-context plotting arguments must be provided together.')
        rows96 = load_master(args.master_96)
        plot_sbs96_rate_pdf(rows96, args.output_per_species_96, per_species=True)
        plot_sbs96_rate_pdf(rows96, args.output_group_96, per_species=False)
        plot_sbs96_fold_change_pdf(load_fc(args.fold_change_96), args.output_fold_change_96)
        print(f'Wrote {args.output_per_species_96}')
        print(f'Wrote {args.output_group_96}')
        print(f'Wrote {args.output_fold_change_96}')


if __name__ == '__main__':
    main()
