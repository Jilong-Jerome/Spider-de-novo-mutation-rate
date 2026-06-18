#!/usr/bin/env python3
"""
Final individual mutation-rate variation outputs.

This script creates per-species observed variance/CV tables, pairwise CV
simulation tests, and a tree-aligned PDF plot from the current workflow
configuration files and raw DNM/callable-site inputs.
"""

import argparse
import os
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DNM_WORKFLOW_SOURCES = os.path.join(REPO_ROOT, 'dnm_variance', 'workflow_sources')
sys.path.insert(0, DNM_WORKFLOW_SOURCES)

from dnm_individual_variance_simulation import (  # noqa: E402
    DIPLOID_FACTOR,
    count_dnms_per_individual,
    compute_individual_stats,
    parse_callable_sites,
    parse_dnm_file,
)


SPECIES_ORDER = ['lin', 'mim', 'afr', 'sar', 'bic', 'dum', 'ten']
SPECIES_LABELS = {sp: sp.upper() for sp in SPECIES_ORDER}
TREE_TEXT = '(LIN:subsocial,((MIM:social,AFR:subsocial),((SAR:social,BIC:subsocial),(DUM:social,TEN:subsocial))))'


def load_species_stats(config_path, sociality_by_species):
    with open(config_path) as f:
        config = yaml.safe_load(f)

    species = config['species_prefix']
    x_chroms_set = set(config['x_chromosomes'])
    excluded_trios_set = set(config.get('excluded_trios', []))

    callable_per_ind = parse_callable_sites(config['callable_file'], x_chroms_set)
    dnm_df = parse_dnm_file(config['dnm_file'], x_chroms_set, excluded_trios_set)
    dnm_counts = count_dnms_per_individual(dnm_df)
    individual_stats = compute_individual_stats(
        callable_per_ind,
        dnm_counts,
        excluded_trios_set,
    )

    if not individual_stats:
        raise ValueError(f'No individuals available after filtering for {species}')

    total_dnms = sum(s['unique_dnms'] for s in individual_stats)
    total_callable = sum(s['callable_sites'] for s in individual_stats)
    mean_rate = total_dnms / (DIPLOID_FACTOR * total_callable)
    rates = np.array([s['rate'] for s in individual_stats], dtype=float)
    variance = float(np.var(rates, ddof=0))
    mean_individual_rate = float(np.mean(rates))
    cv = float(np.std(rates, ddof=0) / mean_individual_rate) if mean_individual_rate > 0 else 0.0
    callable_arr = np.array([s['callable_sites'] for s in individual_stats], dtype=float)

    return {
        'species': species,
        'sociality': sociality_by_species.get(species, 'unclassified'),
        'n_trios': len(individual_stats),
        'total_dnms': total_dnms,
        'total_callable_sites': total_callable,
        'mean_mutation_rate': mean_rate,
        'individual_variance_ddof0': variance,
        'individual_cv_ddof0': cv,
        'callable_arr': callable_arr,
        'lam_arr': DIPLOID_FACTOR * callable_arr * mean_rate,
    }


def simulate_cv(species_record, n_simulations, rng):
    counts = rng.poisson(
        lam=species_record['lam_arr'],
        size=(n_simulations, len(species_record['lam_arr'])),
    )
    rates = counts / species_record['callable_arr'][np.newaxis, :]
    means = rates.mean(axis=1)
    stds = rates.std(axis=1, ddof=0)
    return np.where(means > 0, stds / means, 0.0)


def ordered_species_records(records_by_species):
    missing = [sp for sp in SPECIES_ORDER if sp not in records_by_species]
    if missing:
        raise ValueError(f'Missing species required for final tree plot: {", ".join(missing)}')
    return [records_by_species[sp] for sp in SPECIES_ORDER]


def write_species_tables(records, output_variance_tsv, output_cv_tsv):
    table = pd.DataFrame(records)
    shared_cols = [
        'species',
        'sociality',
        'n_trios',
        'total_dnms',
        'total_callable_sites',
        'mean_mutation_rate',
    ]

    variance_table = table[shared_cols + ['individual_variance_ddof0']]
    cv_table = table[shared_cols + ['individual_cv_ddof0']]

    variance_table.to_csv(output_variance_tsv, sep='\t', index=False)
    cv_table.to_csv(output_cv_tsv, sep='\t', index=False)


def run_pairwise_tests(records_by_species, sister_pairs, n_simulations, random_seed):
    rng = np.random.default_rng(random_seed)
    sim_cvs = {}
    for species in sorted(records_by_species):
        sim_cvs[species] = simulate_cv(records_by_species[species], n_simulations, rng)

    pair_rows = []
    pair_simulated_deltas = []
    for subsocial_species, social_species in sister_pairs:
        if subsocial_species not in records_by_species:
            raise ValueError(f'Missing subsocial species in sister pair: {subsocial_species}')
        if social_species not in records_by_species:
            raise ValueError(f'Missing social species in sister pair: {social_species}')

        cv_subsocial = records_by_species[subsocial_species]['individual_cv_ddof0']
        cv_social = records_by_species[social_species]['individual_cv_ddof0']
        observed_delta = cv_subsocial - cv_social
        simulated_delta = sim_cvs[subsocial_species] - sim_cvs[social_species]
        p_value = float(np.mean(simulated_delta >= observed_delta))

        pair_rows.append({
            'subsocial_species': subsocial_species,
            'social_species': social_species,
            'cv_subsocial': cv_subsocial,
            'cv_social': cv_social,
            'delta_cv_subsocial_minus_social': observed_delta,
            'p_value_upper_tail': p_value,
            'n_simulations': n_simulations,
            'random_seed': random_seed,
        })
        pair_simulated_deltas.append({
            'subsocial_species': subsocial_species,
            'social_species': social_species,
            'simulated_delta_cv': simulated_delta,
            'observed_delta_cv': observed_delta,
            'p_value_upper_tail': p_value,
        })

    if pair_rows:
        all_positive = np.ones(n_simulations, dtype=bool)
        for row in pair_rows:
            all_positive &= sim_cvs[row['subsocial_species']] > sim_cvs[row['social_species']]
        combined_positive_p = float(np.mean(all_positive))
        observed_positive_pairs = sum(
            row['delta_cv_subsocial_minus_social'] > 0 for row in pair_rows
        )
    else:
        combined_positive_p = np.nan
        observed_positive_pairs = 0

    return pair_rows, pair_simulated_deltas, sim_cvs, combined_positive_p, observed_positive_pairs


def write_cv_simulations_table(records, sim_cvs, output_cv_simulations_tsv):
    rows = []
    for record in records:
        species = record['species']
        for simulation_id, simulated_cv in enumerate(sim_cvs[species]):
            rows.append({
                'species': species,
                'sociality': record['sociality'],
                'simulation_id': simulation_id,
                'simulated_cv': simulated_cv,
                'observed_cv': record['individual_cv_ddof0'],
            })

    pd.DataFrame(rows).to_csv(output_cv_simulations_tsv, sep='\t', index=False)


def write_pairwise_outputs(
    pair_rows,
    records,
    output_pairwise_tsv,
    output_summary,
    combined_positive_p,
    observed_positive_pairs,
    n_simulations,
    random_seed,
):
    pd.DataFrame(pair_rows).to_csv(output_pairwise_tsv, sep='\t', index=False)

    subsocial_cvs = [
        r['individual_cv_ddof0'] for r in records if r['sociality'] == 'subsocial'
    ]
    social_cvs = [
        r['individual_cv_ddof0'] for r in records if r['sociality'] == 'social'
    ]

    with open(output_summary, 'w') as out:
        out.write('Final Individual-Level Mutation-Rate CV Pairwise Tests\n')
        out.write('=' * 58 + '\n\n')
        out.write(f'N simulations: {n_simulations}\n')
        out.write(f'Random seed: {random_seed}\n')
        out.write('Null model: uniform species mutation rate with observed trio count and callable-site counts.\n\n')
        out.write('Per-species observed individual-level CV:\n')
        out.write(f"{'Species':<8} {'Sociality':<12} {'N_trios':>8} {'CV':>12}\n")
        out.write('-' * 44 + '\n')
        for record in records:
            out.write(
                f"{record['species']:<8} {record['sociality']:<12} "
                f"{record['n_trios']:>8} {record['individual_cv_ddof0']:>12.6f}\n"
            )
        out.write('\n')

        if subsocial_cvs and social_cvs:
            out.write(f"Mean CV, subsocial: {np.mean(subsocial_cvs):.6f}\n")
            out.write(f"Mean CV, social:    {np.mean(social_cvs):.6f}\n\n")

        out.write('Sister-pair tests, subsocial - social:\n')
        out.write(
            f"{'Pair':<16} {'CV_sub':>12} {'CV_soc':>12} "
            f"{'Delta_CV':>12} {'P_upper':>10}\n"
        )
        out.write('-' * 68 + '\n')
        for row in pair_rows:
            pair = f"{row['subsocial_species']}:{row['social_species']}"
            out.write(
                f"{pair:<16} {row['cv_subsocial']:>12.6f} "
                f"{row['cv_social']:>12.6f} "
                f"{row['delta_cv_subsocial_minus_social']:>+12.6f} "
                f"{row['p_value_upper_tail']:>10.4f}\n"
            )
        out.write('\n')
        out.write(
            f'Observed positive sister pairs (subsocial > social): '
            f'{observed_positive_pairs}/{len(pair_rows)}\n'
        )
        out.write(
            f'P-value for all {len(pair_rows)} simulated pairs showing '
            f'subsocial > social: {combined_positive_p:.4f}\n'
        )


def draw_node(ax, x, y, tree, y_lookup):
    if isinstance(tree, str):
        label = SPECIES_LABELS[tree]
        ax.plot([x, 4.0], [y_lookup[tree], y_lookup[tree]], color='black', lw=1.2)
        ax.text(4.08, y_lookup[tree], label, va='center', ha='left', fontsize=9)
        return y_lookup[tree]

    child_ys = []
    child_x = x + 0.75
    for child in tree:
        child_y = draw_node(ax, child_x, None, child, y_lookup)
        ax.plot([x, child_x], [child_y, child_y], color='black', lw=1.2)
        child_ys.append(child_y)
    ax.plot([x, x], [min(child_ys), max(child_ys)], color='black', lw=1.2)
    return float(np.mean(child_ys))


def draw_tree_panel(ax, y_lookup):
    tree = (
        'lin',
        (
            ('mim', 'afr'),
            (
                ('sar', 'bic'),
                ('dum', 'ten'),
            ),
        ),
    )
    draw_node(ax, 0.1, None, tree, y_lookup)
    ax.set_xlim(-0.05, 4.9)
    ax.set_ylim(-0.6, len(SPECIES_ORDER) - 0.4)
    ax.set_title('a. Species tree', loc='left', fontsize=11, fontweight='bold')
    ax.text(0.0, -0.45, TREE_TEXT, fontsize=6.5, ha='left', va='bottom')
    ax.axis('off')


def draw_metric_panel(ax, records, y_lookup, metric, title, xlabel, color_by_sociality):
    values = [r[metric] for r in records]
    xmax = max(values) * 1.18 if values and max(values) > 0 else 1.0

    for record in records:
        y = y_lookup[record['species']]
        value = record[metric]
        color = color_by_sociality.get(record['sociality'], '#666666')
        ax.barh(y, value, height=0.42, color=color, edgecolor='black', linewidth=0.4)
        ax.text(value + xmax * 0.02, y, format_metric_value(value), va='center', fontsize=8)

    ax.set_ylim(-0.6, len(SPECIES_ORDER) - 0.4)
    ax.set_xlim(0, xmax)
    ax.set_yticks([y_lookup[sp] for sp in SPECIES_ORDER])
    ax.set_yticklabels([])
    ax.tick_params(axis='y', length=0)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_title(title, loc='left', fontsize=11, fontweight='bold')
    ax.grid(axis='x', color='#d9d9d9', linewidth=0.6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)


def format_metric_value(value):
    if value == 0:
        return '0'
    if abs(value) < 0.001:
        return f'{value:.2e}'
    return f'{value:.3f}'


def write_plot(records, output_plot_pdf):
    y_lookup = {sp: len(SPECIES_ORDER) - 1 - idx for idx, sp in enumerate(SPECIES_ORDER)}
    color_by_sociality = {
        'social': '#3B78B8',
        'subsocial': '#D0782A',
        'unclassified': '#777777',
    }

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(11.5, 5.2),
        gridspec_kw={'width_ratios': [1.45, 1.05, 1.0], 'wspace': 0.22},
        sharey=True,
    )

    draw_tree_panel(axes[0], y_lookup)
    draw_metric_panel(
        axes[1],
        records,
        y_lookup,
        'individual_variance_ddof0',
        'b. Individual-rate variance',
        'Variance',
        color_by_sociality,
    )
    draw_metric_panel(
        axes[2],
        records,
        y_lookup,
        'individual_cv_ddof0',
        'c. Individual-rate CV',
        'Coefficient of variation',
        color_by_sociality,
    )

    handles = [
        plt.Line2D([0], [0], color=color_by_sociality['social'], lw=6, label='Social'),
        plt.Line2D([0], [0], color=color_by_sociality['subsocial'], lw=6, label='Subsocial'),
    ]
    fig.legend(handles=handles, loc='lower center', ncol=2, frameon=False, bbox_to_anchor=(0.62, 0.01))
    fig.subplots_adjust(bottom=0.15)
    fig.savefig(output_plot_pdf, bbox_inches='tight')
    plt.close(fig)


def write_pairwise_distribution_plot(pair_simulated_deltas, output_pairwise_plot_pdf):
    n_pairs = len(pair_simulated_deltas)
    if n_pairs == 0:
        raise ValueError('No pairwise simulated CV differences available for plotting')

    fig, axes = plt.subplots(
        n_pairs,
        1,
        figsize=(7.0, 2.6 * n_pairs),
        sharex=False,
        constrained_layout=True,
    )
    if n_pairs == 1:
        axes = [axes]

    for ax, pair_data in zip(axes, pair_simulated_deltas):
        simulated_delta = pair_data['simulated_delta_cv']
        observed_delta = pair_data['observed_delta_cv']
        pair_label = f"{pair_data['subsocial_species'].upper()} - {pair_data['social_species'].upper()}"

        ax.hist(
            simulated_delta,
            bins=40,
            density=True,
            color='#9DB8D5',
            edgecolor='white',
            linewidth=0.5,
        )
        ax.axvline(
            observed_delta,
            color='#B53A2D',
            linestyle='--',
            linewidth=1.8,
            label='Observed difference',
        )
        ax.axvline(0, color='#555555', linestyle=':', linewidth=1.0)

        ax.set_title(pair_label, loc='left', fontsize=11, fontweight='bold')
        ax.set_ylabel('Density', fontsize=9)
        ax.grid(axis='y', color='#d9d9d9', linewidth=0.6)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.text(
            0.98,
            0.88,
            f"Observed delta = {observed_delta:+.3f}\n"
            f"P_upper = {pair_data['p_value_upper_tail']:.4f}",
            transform=ax.transAxes,
            ha='right',
            va='top',
            fontsize=9,
            bbox={'boxstyle': 'round,pad=0.25', 'facecolor': 'white', 'edgecolor': '#bbbbbb'},
        )

    axes[-1].set_xlabel('Simulated CV difference (subsocial - social)', fontsize=9)
    fig.suptitle(
        'Pairwise null distributions of individual-rate CV differences',
        fontsize=12,
        fontweight='bold',
    )
    fig.savefig(output_pairwise_plot_pdf, bbox_inches='tight')
    plt.close(fig)


def write_cv_simulations_plot(records, sim_cvs, output_cv_simulations_plot_pdf):
    n_species = len(records)
    n_cols = 2
    n_rows = int(np.ceil(n_species / n_cols))
    color_by_sociality = {
        'social': '#3B78B8',
        'subsocial': '#D0782A',
        'unclassified': '#777777',
    }

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(9.0, 2.35 * n_rows),
        constrained_layout=True,
    )
    axes = np.asarray(axes).reshape(-1)

    for ax, record in zip(axes, records):
        species = record['species']
        simulated_cv = sim_cvs[species]
        observed_cv = record['individual_cv_ddof0']
        p_value = float(np.mean(simulated_cv >= observed_cv))
        color = color_by_sociality.get(record['sociality'], '#777777')

        ax.hist(
            simulated_cv,
            bins=40,
            density=True,
            color=color,
            alpha=0.65,
            edgecolor='white',
            linewidth=0.5,
        )
        ax.axvline(
            observed_cv,
            color='#B53A2D',
            linestyle='--',
            linewidth=1.8,
        )
        ax.set_title(
            f"{species.upper()} ({record['sociality']})",
            loc='left',
            fontsize=10,
            fontweight='bold',
        )
        ax.set_ylabel('Density', fontsize=8)
        ax.grid(axis='y', color='#d9d9d9', linewidth=0.6)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.text(
            0.98,
            0.88,
            f"Observed CV = {observed_cv:.3f}\nP_upper = {p_value:.4f}",
            transform=ax.transAxes,
            ha='right',
            va='top',
            fontsize=8,
            bbox={'boxstyle': 'round,pad=0.25', 'facecolor': 'white', 'edgecolor': '#bbbbbb'},
        )

    for ax in axes[n_species:]:
        ax.axis('off')

    for ax in axes[-n_cols:]:
        ax.set_xlabel('Simulated individual-rate CV', fontsize=8)

    fig.suptitle(
        'Per-species null distributions of individual-rate CV',
        fontsize=12,
        fontweight='bold',
    )
    fig.savefig(output_cv_simulations_plot_pdf, bbox_inches='tight')
    plt.close(fig)


def parse_sister_pairs(pair_strings):
    pairs = []
    for pair_string in pair_strings:
        parts = pair_string.split(':')
        if len(parts) != 2:
            raise ValueError(f'Invalid sister pair: {pair_string}')
        pairs.append((parts[0], parts[1]))
    return pairs


def main():
    parser = argparse.ArgumentParser(
        description='Create final individual mutation-rate variation outputs.'
    )
    parser.add_argument('--config_files', nargs='+', required=True)
    parser.add_argument('--cv_config_file', required=True)
    parser.add_argument('--output_variance_tsv', required=True)
    parser.add_argument('--output_cv_tsv', required=True)
    parser.add_argument('--output_cv_simulations_tsv', required=True)
    parser.add_argument('--output_pairwise_tsv', required=True)
    parser.add_argument('--output_summary', required=True)
    parser.add_argument('--output_plot_pdf', required=True)
    parser.add_argument('--output_pairwise_plot_pdf', required=True)
    parser.add_argument('--output_cv_simulations_plot_pdf', required=True)
    args = parser.parse_args()

    with open(args.cv_config_file) as f:
        cv_config = yaml.safe_load(f)

    sociality_by_species = {}
    for species in cv_config['subsocial']:
        sociality_by_species[species] = 'subsocial'
    for species in cv_config['social']:
        sociality_by_species[species] = 'social'

    os.makedirs(os.path.dirname(args.output_variance_tsv), exist_ok=True)
    os.makedirs(os.path.dirname(args.output_cv_tsv), exist_ok=True)
    os.makedirs(os.path.dirname(args.output_cv_simulations_tsv), exist_ok=True)
    os.makedirs(os.path.dirname(args.output_pairwise_tsv), exist_ok=True)
    os.makedirs(os.path.dirname(args.output_summary), exist_ok=True)
    os.makedirs(os.path.dirname(args.output_plot_pdf), exist_ok=True)
    os.makedirs(os.path.dirname(args.output_pairwise_plot_pdf), exist_ok=True)
    os.makedirs(os.path.dirname(args.output_cv_simulations_plot_pdf), exist_ok=True)

    records_by_species = {}
    for config_file in sorted(args.config_files):
        record = load_species_stats(config_file, sociality_by_species)
        records_by_species[record['species']] = record

    records = ordered_species_records(records_by_species)
    write_species_tables(records, args.output_variance_tsv, args.output_cv_tsv)

    pair_rows, pair_simulated_deltas, sim_cvs, combined_positive_p, observed_positive_pairs = run_pairwise_tests(
        records_by_species=records_by_species,
        sister_pairs=parse_sister_pairs(cv_config['sister_pairs']),
        n_simulations=cv_config.get('n_simulations', 10000),
        random_seed=cv_config.get('random_seed', 42),
    )
    write_cv_simulations_table(records, sim_cvs, args.output_cv_simulations_tsv)
    write_pairwise_outputs(
        pair_rows=pair_rows,
        records=records,
        output_pairwise_tsv=args.output_pairwise_tsv,
        output_summary=args.output_summary,
        combined_positive_p=combined_positive_p,
        observed_positive_pairs=observed_positive_pairs,
        n_simulations=cv_config.get('n_simulations', 10000),
        random_seed=cv_config.get('random_seed', 42),
    )
    write_plot(records, args.output_plot_pdf)
    write_pairwise_distribution_plot(pair_simulated_deltas, args.output_pairwise_plot_pdf)
    write_cv_simulations_plot(records, sim_cvs, args.output_cv_simulations_plot_pdf)


if __name__ == '__main__':
    main()
