#!/usr/bin/env python3
"""
Generate the all-species combined germline/somatic mutation-rate visual and
companion TSV datasets from existing dnm_rate master outputs.
"""
import argparse
import csv
import math
import os
import random
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
from scipy.stats import chi2


CLASSES = ['overall', 'C>A', 'C>G', 'C>T_nonCpG', 'C>T_CpG', 'T>A', 'T>C', 'T>G']
SPECTRUM_CLASSES = ['C>A', 'C>G', 'C>T_nonCpG', 'C>T_CpG', 'T>A', 'T>C', 'T>G']
GERMLINE_COLOR = '#8A6A00'
SOMATIC_COLOR = '#4A4A4A'
FOLD_CHANGE_COLOR = '#222222'
MUTATION_CLASS_COLORS = {
    'C>A': '#03BCEE',
    'C>G': '#010101',
    'C>T_nonCpG': '#E32926',
    'C>T_CpG': '#E32926',
    'T>A': '#CAC9C9',
    'T>C': '#A1CE63',
    'T>G': '#EBC6C4',
}
SPECTRUM_LABELS = {
    'C>A': 'C>A',
    'C>G': 'C>G',
    'C>T_nonCpG': 'C>T\n(non-CpG)',
    'C>T_CpG': 'C>T\n(CpG)',
    'T>A': 'T>A',
    'T>C': 'T>C',
    'T>G': 'T>G',
}


def poisson_ci(n, denom):
    if denom <= 0:
        return float('nan'), float('nan')
    if n == 0:
        ci_low = 0.0
    else:
        ci_low = chi2.ppf(0.025, 2 * n) / 2.0 / denom
    ci_high = chi2.ppf(0.975, 2 * (n + 1)) / 2.0 / denom
    return ci_low, ci_high


def read_tsv(path):
    with open(path, newline='') as f:
        reader = csv.DictReader((line for line in f if not line.startswith('#')), delimiter='\t')
        return list(reader)


def parse_int(value):
    if value is None or value == '':
        return None
    return int(value)


def fmt(value):
    if value is None:
        return ''
    if isinstance(value, float):
        if math.isnan(value):
            return 'nan'
        return f'{value:.10g}'
    return str(value)


def percentile(values, q):
    if not values:
        return float('nan')
    return float(np.percentile(np.array(values, dtype=float), q))


def load_germline(per_trio_path):
    raw_rows = read_tsv(per_trio_path)
    input_rows = []
    by_species_offspring = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [0, 0])))
    totals = {cls: [0, 0] for cls in CLASSES}

    for row in raw_rows:
        if row['mutation_class'] not in CLASSES:
            continue
        if row['chrom_group'] != 'autosome':
            continue

        denom_bp = parse_int(row['callable_denom_bp'])
        n_mut = parse_int(row['n_germline_mutations'])
        if denom_bp is None or n_mut is None:
            continue

        input_rows.append({
            'level': 'germline',
            'species': row['species'],
            'offspring_id': row['offspring_id'],
            'chrom': row['chrom'],
            'chrom_group': row['chrom_group'],
            'mutation_class': row['mutation_class'],
            'n_mutations': n_mut,
            'denominator': denom_bp,
            'included': row['included'],
        })

        if row['included'] != 'yes':
            continue

        cls = row['mutation_class']
        species = row['species']
        offspring = row['offspring_id']
        by_species_offspring[species][offspring][cls][0] += n_mut
        by_species_offspring[species][offspring][cls][1] += denom_bp
        totals[cls][0] += n_mut
        totals[cls][1] += denom_bp

    return input_rows, by_species_offspring, totals


def load_somatic(per_species_path):
    raw_rows = read_tsv(per_species_path)
    input_rows = []
    totals = {cls: [0, 0] for cls in CLASSES}

    for row in raw_rows:
        cls = row['mutation_class']
        if cls not in CLASSES:
            continue
        n_mut = int(row['n_somatic_mutations'])
        denom = int(row['callable_trinuc_count'])
        input_rows.append({
            'level': 'somatic',
            'species': row['species'],
            'offspring_id': '',
            'chrom': '',
            'chrom_group': '',
            'mutation_class': cls,
            'n_mutations': n_mut,
            'denominator': denom,
            'included': 'yes',
        })
        totals[cls][0] += n_mut
        totals[cls][1] += denom

    return input_rows, totals


def load_spectrum_test(spectrum_test_path):
    rows = read_tsv(spectrum_test_path)
    out = {}
    for row in rows:
        if row['level'] != '7':
            continue
        cls = row['category']
        if cls not in SPECTRUM_CLASSES:
            continue
        out[cls] = {
            'germline_count': int(row['germline_count']),
            'somatic_count': int(row['somatic_count']),
            'germline_prop': float(row['germline_prop']),
            'somatic_prop': float(row['somatic_prop']),
            'p_adjusted': float(row['p_adjusted']),
            'significant': row['significant'],
        }
    return out


def bootstrap_spectrum_ci(counts, n_reps, seed):
    rng = np.random.default_rng(seed)
    counts_vec = np.array(counts, dtype=float)
    total = counts_vec.sum()
    if total <= 0:
        return [(float('nan'), float('nan')) for _ in counts]
    probs = counts_vec / total
    draws = rng.multinomial(int(total), probs, size=n_reps).astype(float) / total
    lo = np.quantile(draws, 0.025, axis=0)
    hi = np.quantile(draws, 0.975, axis=0)
    return list(zip(lo, hi))


def bootstrap_germline_ci(by_species_offspring, cls, n_reps, seed):
    rng = random.Random(seed)
    species_offspring = {
        species: list(offspring_map.keys())
        for species, offspring_map in by_species_offspring.items()
    }
    rates = []

    for _ in range(n_reps):
        n_total = 0
        denom_bp_total = 0
        for species, offspring_ids in species_offspring.items():
            if not offspring_ids:
                continue
            offspring_map = by_species_offspring[species]
            for _ in range(len(offspring_ids)):
                offspring = rng.choice(offspring_ids)
                n, denom_bp = offspring_map[offspring].get(cls, [0, 0])
                n_total += n
                denom_bp_total += denom_bp
        denom_hap = 2 * denom_bp_total
        if denom_hap > 0:
            rates.append(n_total / denom_hap)

    return percentile(rates, 2.5), percentile(rates, 97.5)


def build_plot_rows(germline_totals, somatic_totals, by_species_offspring,
                    spectrum_rows, n_reps, seed):
    rows = []

    for i, cls in enumerate(CLASSES):
        n_g, denom_bp_g = germline_totals[cls]
        denom_g = 2 * denom_bp_g
        rate_g = n_g / denom_g if denom_g > 0 else float('nan')
        lo_g, hi_g = bootstrap_germline_ci(by_species_offspring, cls, n_reps, seed + i)
        rows.append({
            'row_type': 'rate',
            'level': 'germline',
            'mutation_class': cls,
            'n_mutations': n_g,
            'denominator': denom_g,
            'rate': rate_g,
            'ci_low': lo_g,
            'ci_high': hi_g,
            'ci_method': 'stratified_trio_bootstrap_percentile',
            'n_bootstrap_replicates': n_reps,
        })

        n_s, denom_s = somatic_totals[cls]
        rate_s = n_s / denom_s if denom_s > 0 else float('nan')
        lo_s, hi_s = poisson_ci(n_s, denom_s)
        rows.append({
            'row_type': 'rate',
            'level': 'somatic',
            'mutation_class': cls,
            'n_mutations': n_s,
            'denominator': denom_s,
            'rate': rate_s,
            'ci_low': lo_s,
            'ci_high': hi_s,
            'ci_method': 'exact_poisson',
            'n_bootstrap_replicates': '',
        })

        if denom_g <= 0 or denom_s <= 0 or n_g == 0 or n_s == 0:
            fold_change = rate_s / rate_g if rate_g and rate_g > 0 else float('nan')
            fc_low = float('nan')
            fc_high = float('nan')
        else:
            fold_change = rate_s / rate_g
            log_fc = math.log(fold_change)
            se = math.sqrt(1.0 / n_g + 1.0 / n_s)
            fc_low = math.exp(log_fc - 1.96 * se)
            fc_high = math.exp(log_fc + 1.96 * se)

        rows.append({
            'row_type': 'fold_change',
            'level': 'somatic_vs_germline',
            'mutation_class': cls,
            'n_mutations': '',
            'denominator': '',
            'rate': '',
            'ci_low': '',
            'ci_high': '',
            'ci_method': 'poisson_log_rate_ratio',
            'n_bootstrap_replicates': '',
            'germline_rate': rate_g,
            'somatic_rate': rate_s,
            'fold_change': fold_change,
            'fc_ci_low': fc_low,
            'fc_ci_high': fc_high,
        })

    total_g = sum(r['germline_count'] for r in spectrum_rows.values())
    total_s = sum(r['somatic_count'] for r in spectrum_rows.values())
    g_counts = [spectrum_rows[cls]['germline_count'] for cls in SPECTRUM_CLASSES]
    s_counts = [spectrum_rows[cls]['somatic_count'] for cls in SPECTRUM_CLASSES]
    g_ci = bootstrap_spectrum_ci(g_counts, n_reps, seed + 1000)
    s_ci = bootstrap_spectrum_ci(s_counts, n_reps, seed + 2000)

    for idx, cls in enumerate(SPECTRUM_CLASSES):
        spectrum = spectrum_rows[cls]
        rows.append({
            'row_type': 'spectrum',
            'level': 'germline',
            'mutation_class': cls,
            'ci_method': 'fisher_exact_bh_from_spectrum_workflow',
            'spectrum_count': spectrum['germline_count'],
            'spectrum_total': total_g,
            'spectrum_proportion': spectrum['germline_prop'],
            'spectrum_ci_low': g_ci[idx][0],
            'spectrum_ci_high': g_ci[idx][1],
            'p_adjusted': spectrum['p_adjusted'],
            'significant': spectrum['significant'],
        })
        rows.append({
            'row_type': 'spectrum',
            'level': 'somatic',
            'mutation_class': cls,
            'ci_method': 'fisher_exact_bh_from_spectrum_workflow',
            'spectrum_count': spectrum['somatic_count'],
            'spectrum_total': total_s,
            'spectrum_proportion': spectrum['somatic_prop'],
            'spectrum_ci_low': s_ci[idx][0],
            'spectrum_ci_high': s_ci[idx][1],
            'p_adjusted': spectrum['p_adjusted'],
            'significant': spectrum['significant'],
        })

    return rows


def write_input_rows(rows, output):
    os.makedirs(os.path.dirname(output), exist_ok=True)
    fields = [
        'level', 'species', 'offspring_id', 'chrom', 'chrom_group',
        'mutation_class', 'n_mutations', 'denominator', 'included',
    ]
    with open(output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def write_plot_rows(rows, output):
    fields = [
        'row_type', 'level', 'mutation_class', 'n_mutations', 'denominator', 'rate',
        'ci_low', 'ci_high', 'ci_method', 'n_bootstrap_replicates',
        'germline_rate', 'somatic_rate', 'fold_change', 'fc_ci_low', 'fc_ci_high',
        'spectrum_count', 'spectrum_total', 'spectrum_proportion',
        'spectrum_ci_low', 'spectrum_ci_high', 'p_adjusted', 'significant',
    ]
    with open(output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, '')) for field in fields})


def plot_combined_rates(rows, output_pdf):
    rate_rows = [row for row in rows if row['row_type'] == 'rate']
    fc_rows = [row for row in rows if row['row_type'] == 'fold_change']
    spectrum_rows = [row for row in rows if row['row_type'] == 'spectrum']
    by_key = {(row['level'], row['mutation_class']): row for row in rate_rows}
    fc_by_class = {row['mutation_class']: row for row in fc_rows}
    sig_classes = {
        row['mutation_class'] for row in spectrum_rows
        if row.get('significant') == 'yes'
    }
    y = np.arange(len(CLASSES))
    offsets = {'germline': -0.16, 'somatic': 0.16}
    colors = {'germline': GERMLINE_COLOR, 'somatic': SOMATIC_COLOR}

    fig = plt.figure(figsize=(10, 6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.8, 1, 3], wspace=0.18)
    ax_rate = fig.add_subplot(gs[0, 0])
    ax_fc = fig.add_subplot(gs[0, 1], sharey=ax_rate)
    ax_spec = fig.add_subplot(gs[0, 2])
    for level in ['germline', 'somatic']:
        x = []
        xerr_low = []
        xerr_high = []
        ypos = []
        for idx, cls in enumerate(CLASSES):
            row = by_key[(level, cls)]
            rate = row['rate']
            lo = row['ci_low']
            hi = row['ci_high']
            x.append(rate)
            xerr_low.append(max(rate - lo, 0.0) if not math.isnan(lo) else 0.0)
            xerr_high.append(max(hi - rate, 0.0) if not math.isnan(hi) else 0.0)
            ypos.append(idx + offsets[level])
        ax_rate.errorbar(
            x,
            ypos,
            xerr=[xerr_low, xerr_high],
            fmt='o',
            markersize=5,
            linewidth=1.2,
            elinewidth=1.1,
            capsize=2.5,
            color=colors[level],
            markerfacecolor=colors[level],
            markeredgecolor='white',
            markeredgewidth=0.5,
        )

    fc_x = []
    fc_low = []
    fc_high = []
    for idx, cls in enumerate(CLASSES):
        row = fc_by_class[cls]
        fc = row['fold_change']
        lo = row['fc_ci_low']
        hi = row['fc_ci_high']
        fc_x.append(fc)
        fc_low.append(max(fc - lo, 0.0) if not math.isnan(lo) else 0.0)
        fc_high.append(max(hi - fc, 0.0) if not math.isnan(hi) else 0.0)

    ax_fc.errorbar(
        fc_x,
        y,
        xerr=[fc_low, fc_high],
        fmt='o',
        markersize=5,
        linewidth=1.2,
        elinewidth=1.1,
        capsize=2.5,
        color=FOLD_CHANGE_COLOR,
        markerfacecolor=FOLD_CHANGE_COLOR,
        markeredgecolor='white',
        markeredgewidth=0.5,
    )
    ax_fc.axvline(1.0, color='#8C8C8C', linestyle='--', linewidth=0.9, zorder=0)

    ax_rate.set_yticks(y)
    ax_rate.set_yticklabels(CLASSES)
    ax_rate.invert_yaxis()
    ax_rate.set_xscale('log')
    ax_rate.set_xlabel('Mutation rate per callable site')
    ax_rate.set_ylabel('Mutation class')
    ax_rate.set_title('Combined rate', fontsize=10)
    ax_rate.grid(axis='x', linestyle=':', linewidth=0.6, alpha=0.45)
    rate_handles = [
        Patch(facecolor=GERMLINE_COLOR, edgecolor='none', label='Germline'),
        Patch(facecolor=SOMATIC_COLOR, edgecolor='none', label='Somatic'),
    ]
    ax_rate.legend(
        handles=rate_handles,
        frameon=False,
        loc='lower right',
        fontsize=8,
        handlelength=1.0,
        borderpad=0.2,
        labelspacing=0.3,
    )

    ax_fc.set_xlabel('Somatic / germline fold change')
    ax_fc.set_title('Fold change', fontsize=10)
    ax_fc.grid(axis='x', linestyle=':', linewidth=0.6, alpha=0.45)
    ax_fc.tick_params(axis='y', left=False, labelleft=False)

    x_spec = np.arange(len(SPECTRUM_CLASSES))
    width = 0.38
    spectrum_by_key = {
        (row['level'], row['mutation_class']): row
        for row in spectrum_rows
    }
    class_colors = [MUTATION_CLASS_COLORS[cls] for cls in SPECTRUM_CLASSES]
    for level, x_offset, alpha, hatch in [
        ('germline', -width / 2, 1.0, None),
        ('somatic', width / 2, 0.55, '///'),
    ]:
        vals = []
        err_low = []
        err_high = []
        for cls in SPECTRUM_CLASSES:
            row = spectrum_by_key[(level, cls)]
            prop = row['spectrum_proportion']
            vals.append(prop)
            err_low.append(max(prop - row['spectrum_ci_low'], 0.0))
            err_high.append(max(row['spectrum_ci_high'] - prop, 0.0))
        ax_spec.bar(
            x_spec + x_offset,
            vals,
            width=width,
            color=class_colors,
            edgecolor='black',
            linewidth=0.35,
            alpha=alpha,
            hatch=hatch,
            yerr=[err_low, err_high],
            error_kw=dict(ecolor='black', lw=0.75, capsize=2.5),
        )

    spectrum_tick_labels = [
        f'{SPECTRUM_LABELS[cls]}*' if cls in sig_classes else SPECTRUM_LABELS[cls]
        for cls in SPECTRUM_CLASSES
    ]
    tick_text = ax_spec.set_xticks(x_spec)
    tick_text = ax_spec.set_xticklabels(spectrum_tick_labels, fontsize=8)
    for text, cls in zip(tick_text, SPECTRUM_CLASSES):
        if cls in sig_classes:
            text.set_fontweight('bold')
    ax_spec.set_ylabel('Fraction')
    ax_spec.set_xlabel('Mutation class')
    ax_spec.set_title('Spectrum proportion', fontsize=10)
    ax_spec.grid(axis='y', linestyle=':', linewidth=0.6, alpha=0.45)
    ax_spec.set_ylim(bottom=0.0)
    spectrum_handles = [
        Patch(facecolor='#BDBDBD', edgecolor='black', linewidth=0.35,
              label='Germline'),
        Patch(facecolor='#BDBDBD', edgecolor='black', linewidth=0.35,
              hatch='///', alpha=0.55, label='Somatic'),
    ]
    ax_spec.legend(
        handles=spectrum_handles,
        frameon=False,
        loc='upper right',
        fontsize=8,
        handlelength=1.0,
        borderpad=0.2,
        labelspacing=0.3,
    )

    for ax in (ax_rate, ax_fc, ax_spec):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.subplots_adjust(left=0.12, right=0.985, bottom=0.14, top=0.92, wspace=0.18)
    fig.savefig(output_pdf)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--per_trio_germline', required=True)
    ap.add_argument('--per_species_somatic', required=True)
    ap.add_argument('--spectrum_test', required=True)
    ap.add_argument('--output_input_data', required=True)
    ap.add_argument('--output_plot_data', required=True)
    ap.add_argument('--output_pdf', required=True)
    ap.add_argument('--bootstrap_replicates', type=int, default=10000)
    ap.add_argument('--bootstrap_seed', type=int, default=20260430)
    args = ap.parse_args()

    germline_input, by_species_offspring, germline_totals = load_germline(args.per_trio_germline)
    somatic_input, somatic_totals = load_somatic(args.per_species_somatic)
    spectrum_rows = load_spectrum_test(args.spectrum_test)
    input_rows = germline_input + somatic_input

    plot_rows = build_plot_rows(
        germline_totals,
        somatic_totals,
        by_species_offspring,
        spectrum_rows,
        args.bootstrap_replicates,
        args.bootstrap_seed,
    )

    write_input_rows(input_rows, args.output_input_data)
    write_plot_rows(plot_rows, args.output_plot_data)
    plot_combined_rates(plot_rows, args.output_pdf)


if __name__ == '__main__':
    main()
