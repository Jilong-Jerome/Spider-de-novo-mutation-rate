#!/usr/bin/env python3
"""
DNM Variance Simulation
=======================
Tests whether inter-family variance in autosomal DNM rates exceeds
what is expected under a Poisson null model (each offspring accumulates
DNMs independently at the species-level empirical mutation rate).
"""

import argparse
import sys
import numpy as np
import pandas as pd


# Genomes are diploid: each callable position carries two alleles that could
# acquire a de novo mutation, so per-bp rates use 2 * callable_sites.
DIPLOID_FACTOR = 2


def parse_callable_sites(callable_file, x_chroms_set):
    """
    Parse callable sites file (no header, 3 cols: chrom, count, individual_chrom).
    Sum autosomal callable sites per individual (skip x_chroms).
    individual_id = col3[:-(len(chrom)+1)]
    """
    callable_per_ind = {}
    with open(callable_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            chrom = parts[0]
            count = int(parts[1])
            ind_chrom = parts[2]
            if chrom in x_chroms_set:
                continue
            individual_id = ind_chrom[:-(len(chrom) + 1)]
            callable_per_ind[individual_id] = callable_per_ind.get(individual_id, 0) + count
    return callable_per_ind


def parse_dnm_file(dnm_file, x_chroms_set, excluded_trios_set):
    """
    Load DNM TSV (with header). Filter out x chromosomes and excluded trios.
    Returns filtered DataFrame.
    """
    dnm_df = pd.read_csv(dnm_file, sep='\t')
    dnm_df = dnm_df[~dnm_df['chrom'].isin(x_chroms_set)]
    if excluded_trios_set:
        dnm_df = dnm_df[~dnm_df['child'].isin(excluded_trios_set)]
    return dnm_df


def derive_family_key(child_id):
    """
    Derive family key from child ID.
    e.g. 'AFR_family1_S3_offspring' -> 'AFR_family1'
    """
    return child_id.rsplit('_S', 1)[0]


def assign_and_deduplicate(dnm_df):
    """
    Assign DNMs to families via child ID parsing.
    Deduplicate by (chrom, pos, ref, alt) within each family.
    Returns dict: {family_key: unique_dnm_count}
    """
    dnm_df = dnm_df.copy()
    dnm_df['family_key'] = dnm_df['child'].apply(derive_family_key)

    family_dnm_counts = {}
    for family_key, group in dnm_df.groupby('family_key'):
        unique_events = group.drop_duplicates(subset=['chrom', 'pos', 'ref', 'alt'])
        family_dnm_counts[family_key] = len(unique_events)

    return family_dnm_counts


def build_family_to_offspring(callable_per_ind, excluded_trios_set):
    """
    Build family -> [offspring_ids] mapping from callable_per_ind keys.
    This includes families with 0 DNMs (not in DNM file).
    """
    family_to_offspring = {}
    for ind_id in callable_per_ind:
        if ind_id in excluded_trios_set:
            continue
        family_key = derive_family_key(ind_id)
        if family_key not in family_to_offspring:
            family_to_offspring[family_key] = []
        family_to_offspring[family_key].append(ind_id)
    return family_to_offspring


def compute_family_stats(family_to_offspring, family_dnm_counts, callable_per_ind):
    """
    Compute per-family stats: unique_dnms, callable_sites, family_rate.
    Warns and skips families with 0 callable sites.
    Returns list of dicts with keys: family, unique_dnms, callable_sites, family_rate.
    """
    stats = []
    for family_key, offspring_ids in sorted(family_to_offspring.items()):
        callable_sites_f = sum(callable_per_ind.get(i, 0) for i in offspring_ids)
        if callable_sites_f == 0:
            print(f"WARNING: Family {family_key} has 0 callable sites — skipping.", file=sys.stderr)
            continue
        unique_dnms_f = family_dnm_counts.get(family_key, 0)
        family_rate_f = unique_dnms_f / (DIPLOID_FACTOR * callable_sites_f)
        stats.append({
            'family': family_key,
            'offspring': offspring_ids,
            'unique_dnms': unique_dnms_f,
            'callable_sites': callable_sites_f,
            'family_rate': family_rate_f,
        })
    return stats


def run_simulation(family_stats, callable_per_ind, empirical_rate, n_simulations, rng):
    """
    Vectorized Poisson simulation.
    Returns array of simulated inter-family variances, shape (n_simulations,).
    """
    # Build ordered arrays
    family_sizes = []
    lam_arr = []
    fam_callable = []

    for fs in family_stats:
        offspring_ids = fs['offspring']
        family_sizes.append(len(offspring_ids))
        fam_callable.append(fs['callable_sites'])
        for ind_id in offspring_ids:
            lam = DIPLOID_FACTOR * callable_per_ind[ind_id] * empirical_rate
            lam_arr.append(lam)

    lam_arr = np.array(lam_arr)
    fam_callable = np.array(fam_callable, dtype=float)
    family_sizes = np.array(family_sizes, dtype=int)

    N_individuals = len(lam_arr)
    N_families = len(family_stats)

    # Draw all counts at once: shape (n_simulations, N_individuals)
    counts_matrix = rng.poisson(lam=lam_arr, size=(n_simulations, N_individuals))

    # Aggregate per family: split along axis=1 by family sizes
    cumulative = np.cumsum(family_sizes)[:-1]
    family_blocks = np.split(counts_matrix, cumulative, axis=1)

    # Sum within each family block -> shape (n_simulations, N_families)
    fam_counts = np.column_stack([block.sum(axis=1) for block in family_blocks])

    # Per-family simulated rates: (n_simulations, N_families)
    sim_rates = fam_counts / (DIPLOID_FACTOR * fam_callable[np.newaxis, :])

    # Inter-family variance per simulation (ddof=0)
    sim_variances = np.var(sim_rates, axis=1, ddof=0)

    return sim_variances


def main():
    parser = argparse.ArgumentParser(
        description='Test inter-family DNM rate variance against Poisson null.'
    )
    parser.add_argument('--dnm_file', required=True, help='Path to DNM TSV file')
    parser.add_argument('--callable_file', required=True, help='Path to callable sites TSV')
    parser.add_argument('--species_prefix', required=True, help='Species prefix (lowercase)')
    parser.add_argument('--x_chroms', nargs='+', required=True, help='X chromosome names to exclude')
    parser.add_argument('--n_simulations', type=int, default=1000, help='Number of simulations')
    parser.add_argument('--random_seed', type=int, default=42, help='Random seed')
    parser.add_argument('--excluded_trios', nargs='*', default=[], help='Individual IDs to exclude')
    parser.add_argument('--output_simulations', required=True, help='Output TSV for simulation results')
    parser.add_argument('--output_summary', required=True, help='Output summary text file')
    args = parser.parse_args()

    x_chroms_set = set(args.x_chroms)
    excluded_trios_set = set(args.excluded_trios) if args.excluded_trios else set()
    rng = np.random.default_rng(args.random_seed)

    print(f"[{args.species_prefix}] Parsing callable sites...", file=sys.stderr)
    callable_per_ind = parse_callable_sites(args.callable_file, x_chroms_set)
    print(f"  Found {len(callable_per_ind)} individuals in callable file.", file=sys.stderr)

    print(f"[{args.species_prefix}] Parsing DNM file...", file=sys.stderr)
    dnm_df = parse_dnm_file(args.dnm_file, x_chroms_set, excluded_trios_set)
    print(f"  {len(dnm_df)} autosomal DNM rows after filtering.", file=sys.stderr)

    print(f"[{args.species_prefix}] Assigning DNMs to families...", file=sys.stderr)
    family_dnm_counts = assign_and_deduplicate(dnm_df)

    print(f"[{args.species_prefix}] Building family-to-offspring map...", file=sys.stderr)
    family_to_offspring = build_family_to_offspring(callable_per_ind, excluded_trios_set)
    print(f"  Found {len(family_to_offspring)} families.", file=sys.stderr)

    print(f"[{args.species_prefix}] Computing per-family stats...", file=sys.stderr)
    family_stats = compute_family_stats(family_to_offspring, family_dnm_counts, callable_per_ind)

    # Empirical rate
    total_dnms = sum(fs['unique_dnms'] for fs in family_stats)
    total_callable = sum(fs['callable_sites'] for fs in family_stats)
    empirical_rate = total_dnms / (DIPLOID_FACTOR * total_callable)
    print(f"  Empirical rate: {empirical_rate:.6e}", file=sys.stderr)

    # Observed inter-family variance
    family_rates = np.array([fs['family_rate'] for fs in family_stats])
    obs_variance = np.var(family_rates, ddof=0)
    print(f"  Observed inter-family variance: {obs_variance:.6e}", file=sys.stderr)

    print(f"[{args.species_prefix}] Running {args.n_simulations} simulations...", file=sys.stderr)
    sim_variances = run_simulation(
        family_stats, callable_per_ind, empirical_rate, args.n_simulations, rng
    )

    # P-value (one-tailed upper)
    p_value = np.mean(sim_variances >= obs_variance)
    sim_mean = np.mean(sim_variances)
    sim_std = np.std(sim_variances)

    print(f"  P-value: {p_value:.4f}", file=sys.stderr)

    # Write simulations TSV
    sim_df = pd.DataFrame({
        'simulation_id': np.arange(args.n_simulations),
        'simulated_variance': sim_variances,
    })
    sim_df.to_csv(args.output_simulations, sep='\t', index=False)
    print(f"  Wrote simulations to {args.output_simulations}", file=sys.stderr)

    # Write summary
    with open(args.output_summary, 'w') as out:
        out.write(f"DNM Variance Summary — {args.species_prefix.upper()}\n")
        out.write("=" * 60 + "\n\n")

        out.write(f"Empirical species-level mutation rate: {empirical_rate:.6e}\n")
        out.write(f"X chromosomes excluded: {', '.join(sorted(x_chroms_set))}\n")
        if excluded_trios_set:
            out.write(f"Excluded trios: {', '.join(sorted(excluded_trios_set))}\n")
        else:
            out.write("Excluded trios: none\n")
        out.write(f"N simulations: {args.n_simulations}\n")
        out.write(f"Random seed: {args.random_seed}\n\n")

        out.write("Per-family breakdown:\n")
        out.write(f"{'Family':<30} {'N_offspring':>12} {'Unique_DNMs':>12} {'Callable_sites':>15} {'Rate':>14}\n")
        out.write("-" * 90 + "\n")
        for fs in family_stats:
            out.write(
                f"{fs['family']:<30} {len(fs['offspring']):>12} {fs['unique_dnms']:>12} "
                f"{fs['callable_sites']:>15} {fs['family_rate']:>14.6e}\n"
            )
        out.write("\n")

        out.write(f"Observed inter-family variance (ddof=0): {obs_variance:.6e}\n\n")
        out.write(f"Null distribution (from {args.n_simulations} Poisson simulations):\n")
        out.write(f"  Mean simulated variance:  {sim_mean:.6e}\n")
        out.write(f"  Std simulated variance:   {sim_std:.6e}\n\n")
        out.write(f"P-value (fraction sim >= obs): {p_value:.4f}\n\n")

        if p_value < 0.05:
            interpretation = (
                "SIGNIFICANT: Observed inter-family variance exceeds the Poisson null "
                "expectation (p < 0.05), suggesting families differ in their true mutation rates."
            )
        else:
            interpretation = (
                "NOT SIGNIFICANT: Observed inter-family variance is consistent with the "
                "Poisson null expectation (p >= 0.05)."
            )
        out.write(f"Interpretation: {interpretation}\n")

    print(f"  Wrote summary to {args.output_summary}", file=sys.stderr)
    print(f"[{args.species_prefix}] Done.", file=sys.stderr)


if __name__ == '__main__':
    main()
