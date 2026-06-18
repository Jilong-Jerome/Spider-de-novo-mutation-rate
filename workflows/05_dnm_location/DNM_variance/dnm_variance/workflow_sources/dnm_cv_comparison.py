#!/usr/bin/env python3
"""
DNM CV Comparison
=================
Tests whether the coefficient of variation (CV = std/mean) in family-level
autosomal mutation rates is higher in subsocial species than social species,
using a simulation-based null model.

Null model: for each species, the observed mean rate is taken as the true rate.
Per-family DNM counts are simulated under Poisson (same offspring structure as
observed), yielding a null distribution of CV differences across species.

Two tests:
  1. Global: mean_CV(subsocial) - mean_CV(social) vs simulation null
  2. Pairwise: ΔCV per sister-species pair vs simulation null
"""

import argparse
import os
import sys

import numpy as np
import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dnm_variance_simulation import (
    DIPLOID_FACTOR,
    parse_callable_sites,
    parse_dnm_file,
    assign_and_deduplicate,
    build_family_to_offspring,
    compute_family_stats,
)


def compute_species_cv_and_structure(config_path):
    """
    Load config, parse data, compute observed CV and the per-family structure
    needed for simulation.

    Returns dict with keys:
        species_prefix, sociality (set by caller), n_families, mean_rate, cv,
        lam_arr (np array), family_sizes (list), fam_callable (np array)
    """
    with open(config_path) as f:
        config = yaml.safe_load(f)

    sp = config['species_prefix']
    x_chroms_set = set(config['x_chromosomes'])
    excluded_trios_set = set(config.get('excluded_trios', []))

    callable_per_ind = parse_callable_sites(config['callable_file'], x_chroms_set)
    dnm_df = parse_dnm_file(config['dnm_file'], x_chroms_set, excluded_trios_set)
    family_dnm_counts = assign_and_deduplicate(dnm_df)
    family_to_offspring = build_family_to_offspring(callable_per_ind, excluded_trios_set)
    family_stats = compute_family_stats(family_to_offspring, family_dnm_counts, callable_per_ind)

    total_dnms = sum(fs['unique_dnms'] for fs in family_stats)
    total_callable = sum(fs['callable_sites'] for fs in family_stats)
    empirical_rate = total_dnms / (DIPLOID_FACTOR * total_callable)

    family_rates = np.array([fs['family_rate'] for fs in family_stats])
    obs_cv = np.std(family_rates, ddof=0) / np.mean(family_rates)

    # Build simulation structure
    lam_arr = []
    family_sizes = []
    fam_callable = []
    for fs in family_stats:
        family_sizes.append(len(fs['offspring']))
        fam_callable.append(fs['callable_sites'])
        for ind_id in fs['offspring']:
            lam_arr.append(DIPLOID_FACTOR * callable_per_ind[ind_id] * empirical_rate)

    return {
        'species_prefix': sp,
        'n_families': len(family_stats),
        'mean_rate': empirical_rate,
        'cv': obs_cv,
        'lam_arr': np.array(lam_arr),
        'family_sizes': family_sizes,
        'fam_callable': np.array(fam_callable, dtype=float),
    }


def simulate_cv(species_data, n_simulations, rng):
    """
    Run vectorized Poisson simulations for one species.
    Returns array of simulated CVs, shape (n_simulations,).
    """
    lam_arr = species_data['lam_arr']
    family_sizes = species_data['family_sizes']
    fam_callable = species_data['fam_callable']

    counts_matrix = rng.poisson(lam=lam_arr, size=(n_simulations, len(lam_arr)))

    cumulative = np.cumsum(family_sizes)[:-1]
    blocks = np.split(counts_matrix, cumulative, axis=1)
    fam_counts = np.column_stack([b.sum(axis=1) for b in blocks])

    sim_rates = fam_counts / fam_callable[np.newaxis, :]
    sim_means = sim_rates.mean(axis=1)
    sim_stds = sim_rates.std(axis=1, ddof=0)

    # Guard against zero mean (shouldn't happen with real data, but be safe)
    sim_cv = np.where(sim_means > 0, sim_stds / sim_means, 0.0)
    return sim_cv


def main():
    parser = argparse.ArgumentParser(
        description='Compare CV of family-level mutation rates: subsocial vs social.'
    )
    parser.add_argument('--config_files', nargs='+', required=True,
                        help='All species YAML config paths')
    parser.add_argument('--subsocial', nargs='+', required=True,
                        help='Species prefixes classified as subsocial')
    parser.add_argument('--social', nargs='+', required=True,
                        help='Species prefixes classified as social')
    parser.add_argument('--sister_pairs', nargs='+', required=True,
                        help='Sister-species pairs as subsocial:social (e.g. bic:sar afr:mim ten:dum)')
    parser.add_argument('--n_simulations', type=int, default=10000)
    parser.add_argument('--random_seed', type=int, default=42)
    parser.add_argument('--output_tsv', required=True)
    parser.add_argument('--output_summary', required=True)
    args = parser.parse_args()

    subsocial_set = set(args.subsocial)
    social_set = set(args.social)
    rng = np.random.default_rng(args.random_seed)

    # Parse sister pairs
    sister_pairs = []
    for pair_str in args.sister_pairs:
        sub_sp, soc_sp = pair_str.split(':')
        sister_pairs.append((sub_sp, soc_sp))

    # --- Load data for all species ---
    print("Loading data for all species...", file=sys.stderr)
    all_data = {}
    for config_path in sorted(args.config_files):
        data = compute_species_cv_and_structure(config_path)
        sp = data['species_prefix']
        if sp in subsocial_set:
            data['sociality'] = 'subsocial'
        elif sp in social_set:
            data['sociality'] = 'social'
        else:
            print(f"  WARNING: {sp} not in subsocial or social list — skipping.", file=sys.stderr)
            continue
        all_data[sp] = data
        print(f"  {sp}: CV={data['cv']:.4f}, n_families={data['n_families']}, "
              f"mean_rate={data['mean_rate']:.4e}, sociality={data['sociality']}", file=sys.stderr)

    # Verify all requested species are present
    for sp in list(subsocial_set) + list(social_set):
        if sp not in all_data:
            print(f"  WARNING: {sp} config not found among provided config files.", file=sys.stderr)

    subsocial_sps = [sp for sp in args.subsocial if sp in all_data]
    social_sps = [sp for sp in args.social if sp in all_data]

    # Observed CVs
    obs_cv_sub = np.array([all_data[sp]['cv'] for sp in subsocial_sps])
    obs_cv_soc = np.array([all_data[sp]['cv'] for sp in social_sps])
    obs_global_diff = obs_cv_sub.mean() - obs_cv_soc.mean()

    # --- Run simulations per species ---
    print(f"\nRunning {args.n_simulations} simulations per species...", file=sys.stderr)
    sim_cvs = {}
    for sp, data in all_data.items():
        sim_cvs[sp] = simulate_cv(data, args.n_simulations, rng)
        print(f"  {sp}: sim CV mean={sim_cvs[sp].mean():.4f}", file=sys.stderr)

    # --- Test 1: Global ---
    sim_cv_sub_mat = np.column_stack([sim_cvs[sp] for sp in subsocial_sps])  # (N, 4)
    sim_cv_soc_mat = np.column_stack([sim_cvs[sp] for sp in social_sps])     # (N, 3)
    sim_global_diff = sim_cv_sub_mat.mean(axis=1) - sim_cv_soc_mat.mean(axis=1)
    p_global = np.mean(sim_global_diff >= obs_global_diff)

    # --- Test 2: Pairwise ---
    pair_results = []
    for sub_sp, soc_sp in sister_pairs:
        if sub_sp not in all_data or soc_sp not in all_data:
            print(f"  WARNING: pair {sub_sp}:{soc_sp} missing data — skipping.", file=sys.stderr)
            continue
        obs_delta = all_data[sub_sp]['cv'] - all_data[soc_sp]['cv']
        sim_delta = sim_cvs[sub_sp] - sim_cvs[soc_sp]
        p_pair = np.mean(sim_delta >= obs_delta)
        pair_results.append({
            'subsocial': sub_sp,
            'social': soc_sp,
            'cv_subsocial': all_data[sub_sp]['cv'],
            'cv_social': all_data[soc_sp]['cv'],
            'delta_cv': obs_delta,
            'p_value': p_pair,
        })

    # Combined sign: fraction of sims where all pairs show subsocial > social
    if pair_results:
        all_positive = np.ones(args.n_simulations, dtype=bool)
        for pr in pair_results:
            all_positive &= (sim_cvs[pr['subsocial']] > sim_cvs[pr['social']])
        p_combined_sign = np.mean(all_positive)
        obs_n_positive = sum(1 for pr in pair_results if pr['delta_cv'] > 0)
    else:
        p_combined_sign = np.nan
        obs_n_positive = 0

    # --- Write TSV ---
    with open(args.output_tsv, 'w') as out:
        out.write('species\tsociality\tn_families\tmean_rate\tcv\n')
        for sp in subsocial_sps + social_sps:
            d = all_data[sp]
            out.write(f"{sp}\t{d['sociality']}\t{d['n_families']}\t"
                      f"{d['mean_rate']:.6e}\t{d['cv']:.6f}\n")
    print(f"\nWrote CV table to {args.output_tsv}", file=sys.stderr)

    # --- Write summary ---
    with open(args.output_summary, 'w') as out:
        out.write("DNM Family-Level CV Comparison: Subsocial vs Social\n")
        out.write("=" * 60 + "\n\n")
        out.write(f"N simulations: {args.n_simulations}  |  Random seed: {args.random_seed}\n\n")

        out.write("Per-species CV (observed):\n")
        out.write(f"{'Species':<10} {'Sociality':<12} {'N_fam':>6} {'Mean_rate':>14} {'CV':>10}\n")
        out.write("-" * 56 + "\n")
        for sp in subsocial_sps + social_sps:
            d = all_data[sp]
            out.write(f"{sp:<10} {d['sociality']:<12} {d['n_families']:>6} "
                      f"{d['mean_rate']:>14.4e} {d['cv']:>10.6f}\n")
        out.write("\n")

        out.write(f"Mean CV — subsocial ({', '.join(subsocial_sps)}): "
                  f"{obs_cv_sub.mean():.6f}\n")
        out.write(f"Mean CV — social    ({', '.join(social_sps)}): "
                  f"{obs_cv_soc.mean():.6f}\n\n")

        out.write("Test 1 — Global (subsocial vs social, simulation-based null)\n")
        out.write("-" * 60 + "\n")
        out.write(f"Observed difference (mean_CV_subsocial - mean_CV_social): "
                  f"{obs_global_diff:+.6f}\n")
        out.write(f"Null distribution: mean={sim_global_diff.mean():.6f}, "
                  f"std={sim_global_diff.std():.6f}\n")
        out.write(f"P-value (fraction sim >= obs, one-tailed upper): {p_global:.4f}\n")
        if p_global < 0.05:
            out.write("Interpretation: SIGNIFICANT — subsocial species have higher "
                      "family-rate CV than social species (p < 0.05).\n\n")
        else:
            out.write("Interpretation: NOT SIGNIFICANT — difference in CV is consistent "
                      "with Poisson sampling noise (p >= 0.05).\n\n")

        out.write("Test 2 — Pairwise sister-species comparisons (phylogenetically controlled)\n")
        out.write("-" * 60 + "\n")
        out.write(f"{'Pair':<20} {'CV_sub':>10} {'CV_soc':>10} "
                  f"{'ΔCV':>10} {'p_value':>10}\n")
        out.write("-" * 62 + "\n")
        for pr in pair_results:
            out.write(f"{pr['subsocial']+':'+pr['social']:<20} "
                      f"{pr['cv_subsocial']:>10.6f} {pr['cv_social']:>10.6f} "
                      f"{pr['delta_cv']:>+10.6f} {pr['p_value']:>10.4f}\n")
        out.write("\n")
        out.write(f"Observed positive pairs (subsocial > social): "
                  f"{obs_n_positive}/{len(pair_results)}\n")
        out.write(f"P-value (fraction of sims where ALL {len(pair_results)} pairs "
                  f"show subsocial > social): {p_combined_sign:.4f}\n")

    print(f"Wrote summary to {args.output_summary}", file=sys.stderr)
    print("Done.", file=sys.stderr)


if __name__ == '__main__':
    main()
