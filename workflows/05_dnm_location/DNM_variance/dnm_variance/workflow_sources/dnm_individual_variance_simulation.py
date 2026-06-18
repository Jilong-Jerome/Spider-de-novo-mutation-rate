#!/usr/bin/env python3
"""
DNM Individual-Level Variance Simulation
=========================================
Tests whether inter-individual variance in autosomal DNM rates exceeds
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


def count_dnms_per_individual(dnm_df):
    """
    Count unique (chrom, pos, ref, alt) DNMs per individual (child).
    Returns dict: {individual_id: unique_dnm_count}
    """
    ind_dnm_counts = {}
    for child, group in dnm_df.groupby('child'):
        unique_events = group.drop_duplicates(subset=['chrom', 'pos', 'ref', 'alt'])
        ind_dnm_counts[child] = len(unique_events)
    return ind_dnm_counts


def compute_individual_stats(callable_per_ind, ind_dnm_counts, excluded_trios_set):
    """
    Compute per-individual stats: unique_dnms, callable_sites, rate.
    Uses callable_per_ind keys as the authoritative list (includes 0-DNM individuals).
    Warns and skips individuals with 0 callable sites.
    """
    stats = []
    for ind_id in sorted(callable_per_ind.keys()):
        if ind_id in excluded_trios_set:
            continue
        callable_sites = callable_per_ind[ind_id]
        if callable_sites == 0:
            print(f"WARNING: Individual {ind_id} has 0 callable sites — skipping.", file=sys.stderr)
            continue
        unique_dnms = ind_dnm_counts.get(ind_id, 0)
        rate = unique_dnms / (DIPLOID_FACTOR * callable_sites)
        stats.append({
            'individual': ind_id,
            'unique_dnms': unique_dnms,
            'callable_sites': callable_sites,
            'rate': rate,
        })
    return stats


def run_simulation(individual_stats, empirical_rate, n_simulations, rng):
    """
    Vectorized Poisson simulation at individual level.
    Returns array of simulated inter-individual variances, shape (n_simulations,).
    """
    lam_arr = np.array([DIPLOID_FACTOR * s['callable_sites'] * empirical_rate for s in individual_stats])
    callable_arr = np.array([s['callable_sites'] for s in individual_stats], dtype=float)

    # Draw all counts: shape (n_simulations, N_individuals)
    counts_matrix = rng.poisson(lam=lam_arr, size=(n_simulations, len(lam_arr)))

    # Per-individual simulated rates
    sim_rates = counts_matrix / (DIPLOID_FACTOR * callable_arr[np.newaxis, :])

    # Inter-individual variance per simulation (ddof=0)
    sim_variances = np.var(sim_rates, axis=1, ddof=0)

    return sim_variances


def main():
    parser = argparse.ArgumentParser(
        description='Test inter-individual DNM rate variance against Poisson null.'
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

    print(f"[{args.species_prefix}] Counting DNMs per individual...", file=sys.stderr)
    ind_dnm_counts = count_dnms_per_individual(dnm_df)

    print(f"[{args.species_prefix}] Computing per-individual stats...", file=sys.stderr)
    individual_stats = compute_individual_stats(callable_per_ind, ind_dnm_counts, excluded_trios_set)
    print(f"  {len(individual_stats)} individuals included.", file=sys.stderr)

    # Empirical rate
    total_dnms = sum(s['unique_dnms'] for s in individual_stats)
    total_callable = sum(s['callable_sites'] for s in individual_stats)
    empirical_rate = total_dnms / (DIPLOID_FACTOR * total_callable)
    print(f"  Empirical rate: {empirical_rate:.6e}", file=sys.stderr)

    # Observed inter-individual variance
    ind_rates = np.array([s['rate'] for s in individual_stats])
    obs_variance = np.var(ind_rates, ddof=0)
    print(f"  Observed inter-individual variance: {obs_variance:.6e}", file=sys.stderr)

    print(f"[{args.species_prefix}] Running {args.n_simulations} simulations...", file=sys.stderr)
    sim_variances = run_simulation(individual_stats, empirical_rate, args.n_simulations, rng)

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
        out.write(f"DNM Individual-Level Variance Summary — {args.species_prefix.upper()}\n")
        out.write("=" * 60 + "\n\n")

        out.write(f"Empirical species-level mutation rate: {empirical_rate:.6e}\n")
        out.write(f"X chromosomes excluded: {', '.join(sorted(x_chroms_set))}\n")
        if excluded_trios_set:
            out.write(f"Excluded individuals: {', '.join(sorted(excluded_trios_set))}\n")
        else:
            out.write("Excluded individuals: none\n")
        out.write(f"N individuals: {len(individual_stats)}\n")
        out.write(f"N simulations: {args.n_simulations}\n")
        out.write(f"Random seed: {args.random_seed}\n\n")

        out.write("Per-individual breakdown:\n")
        out.write(f"{'Individual':<40} {'Unique_DNMs':>12} {'Callable_sites':>15} {'Rate':>14}\n")
        out.write("-" * 85 + "\n")
        for s in individual_stats:
            out.write(
                f"{s['individual']:<40} {s['unique_dnms']:>12} "
                f"{s['callable_sites']:>15} {s['rate']:>14.6e}\n"
            )
        out.write("\n")

        out.write(f"Observed inter-individual variance (ddof=0): {obs_variance:.6e}\n\n")
        out.write(f"Null distribution (from {args.n_simulations} Poisson simulations):\n")
        out.write(f"  Mean simulated variance:  {sim_mean:.6e}\n")
        out.write(f"  Std simulated variance:   {sim_std:.6e}\n\n")
        out.write(f"P-value (fraction sim >= obs): {p_value:.4f}\n\n")

        if p_value < 0.05:
            interpretation = (
                "SIGNIFICANT: Observed inter-individual variance exceeds the Poisson null "
                "expectation (p < 0.05), suggesting individuals differ in their true mutation rates."
            )
        else:
            interpretation = (
                "NOT SIGNIFICANT: Observed inter-individual variance is consistent with the "
                "Poisson null expectation (p >= 0.05)."
            )
        out.write(f"Interpretation: {interpretation}\n")

    print(f"  Wrote summary to {args.output_summary}", file=sys.stderr)
    print(f"[{args.species_prefix}] Done.", file=sys.stderr)


if __name__ == '__main__':
    main()
