#!/usr/bin/env python3
"""
nearest_neighbor_test.py

Permutation test for spatial clustering of de novo mutations in the
callable genome.

Observed metric
---------------
Average nearest-neighbour distance (ANND) among deduplicated autosomal DNMs.
For each DNM, the nearest other DNM on the *same chromosome* is found; the
average of those minimum distances is the observed ANND.  Positions on
chromosomes that carry only one DNM contribute no within-chromosome neighbour
and are excluded from the average (same rule applied to simulations).

Null distribution
-----------------
N positions (= number of unique DNM events) are drawn uniformly at random
from the union of all individuals' callable genome segments.  The ANND of
the draw is recorded.  This is repeated n_simulations times.

P-values (two-tailed)
---------------------
  p_cluster      : fraction of simulations with ANND ≤ observed
                   (DNMs more clustered than chance)
  p_dispersed    : fraction of simulations with ANND ≥ observed
                   (DNMs more evenly spaced than chance)
  p_two_tailed   : 2 × min(p_cluster, p_dispersed)

Outputs
-------
  --output_tsv     per-simulation ANND values
  --output_summary plain-text report
  --output_plot    histogram PDF with observed value marked
"""

import argparse
import glob
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def load_callable_union(bed_dir, species_prefix, x_chroms_set):
    """
    Load all per-offspring callable BED TSVs, compute the per-chromosome union
    of segments, and return a list of (chrom, start, end).
    """
    pattern   = os.path.join(bed_dir, f"{species_prefix}_*_callable.bed")
    bed_files = glob.glob(pattern)
    if not bed_files:
        raise FileNotFoundError(f"No callable BED files found matching: {pattern}")

    df = pd.concat(
        [pd.read_csv(f, sep='\t') for f in bed_files],
        ignore_index=True,
    )
    df = df[~df['chrom'].isin(x_chroms_set)].copy()

    union_segments = []
    for chrom, grp in df.groupby('chrom'):
        intervals = sorted(zip(grp['start'].values, grp['end'].values))
        seg_s, seg_e = intervals[0]
        for s, e in intervals[1:]:
            if s <= seg_e:
                seg_e = max(seg_e, e)
            else:
                union_segments.append((chrom, int(seg_s), int(seg_e)))
                seg_s, seg_e = s, e
        union_segments.append((chrom, int(seg_s), int(seg_e)))

    return union_segments


def build_sampling_structure(segments):
    """
    Convert segment list into numpy arrays for fast weighted sampling.
    Returns (chroms, starts, lengths, total_bp).
    """
    chroms  = np.array([s[0] for s in segments])
    starts  = np.array([s[1] for s in segments], dtype=np.int64)
    lengths = np.array([s[2] - s[1] for s in segments], dtype=np.int64)
    total_bp = int(lengths.sum())
    return chroms, starts, lengths, total_bp


def sample_positions(chroms, starts, lengths, total_bp, n, rng):
    """
    Sample n positions uniformly from callable segments.
    Returns dict {chrom: sorted list of int positions}.
    """
    weights  = lengths / total_bp
    seg_idxs = rng.choice(len(chroms), size=n, p=weights)
    offsets  = (rng.random(n) * lengths[seg_idxs]).astype(np.int64)
    positions = starts[seg_idxs] + offsets

    result = defaultdict(list)
    for chrom, pos in zip(chroms[seg_idxs], positions):
        result[chrom].append(int(pos))
    return {c: sorted(ps) for c, ps in result.items()}


def avg_nearest_neighbour_dist(positions_by_chrom):
    """
    Compute average nearest-neighbour distance across all positions that
    have at least one same-chromosome neighbour.

    Returns np.nan if no chromosome carries ≥ 2 positions.
    """
    nn_distances = []
    for positions in positions_by_chrom.values():
        if len(positions) < 2:
            continue
        pos   = np.array(sorted(positions), dtype=np.int64)
        diffs = np.diff(pos)                              # gaps between consecutive
        nn_left  = np.concatenate([[np.inf], diffs])      # distance to left neighbour
        nn_right = np.concatenate([diffs, [np.inf]])      # distance to right neighbour
        nn_dist  = np.minimum(nn_left, nn_right)
        nn_distances.extend(nn_dist[np.isfinite(nn_dist)].tolist())

    return float(np.mean(nn_distances)) if nn_distances else np.nan


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Permutation test: are DNMs spatially clustered in the '
                    'callable genome?'
    )
    parser.add_argument('--bed_dir',        required=True,
                        help='Directory containing per-offspring callable BED TSVs')
    parser.add_argument('--species_prefix', required=True,
                        help='Species prefix for BED file glob (e.g. afr)')
    parser.add_argument('--dnm_file',       required=True,
                        help='DNM TSV (chrom, pos, ref, alt, child, father, mother, ...)')
    parser.add_argument('--x_chroms',       nargs='+', default=[])
    parser.add_argument('--n_simulations',  type=int, default=1000,
                        help='Number of permutation draws (default: 1000)')
    parser.add_argument('--seed',           type=int, default=42,
                        help='Random seed for reproducibility (default: 42)')
    parser.add_argument('--output_tsv',     required=True,
                        help='Per-simulation ANND TSV output')
    parser.add_argument('--output_summary', required=True,
                        help='Plain-text summary output')
    parser.add_argument('--output_plot',    required=True,
                        help='Histogram PDF output')
    args = parser.parse_args()

    x_chroms_set = set(args.x_chroms)

    # ------------------------------------------------------------------ #
    # 1. Build callable genome sampling structure
    # ------------------------------------------------------------------ #
    print("Loading callable genome union ...", flush=True)
    segments = load_callable_union(args.bed_dir, args.species_prefix, x_chroms_set)
    chroms_arr, starts_arr, lengths_arr, total_bp = build_sampling_structure(segments)
    print(f"  {len(segments):,} union segments, {total_bp:,} callable bp", flush=True)

    # ------------------------------------------------------------------ #
    # 2. Load and deduplicate DNMs
    # ------------------------------------------------------------------ #
    dnm_df = pd.read_csv(args.dnm_file, sep='\t')
    dnm_df = dnm_df[~dnm_df['chrom'].isin(x_chroms_set)].copy()
    dnm_dedup = dnm_df.drop_duplicates(
        subset=['chrom', 'pos', 'ref', 'alt', 'father', 'mother']
    )
    n_dnm = len(dnm_dedup)
    print(f"  {n_dnm} unique autosomal DNM events", flush=True)

    if n_dnm < 2:
        sys.exit("Fewer than 2 DNMs — cannot compute nearest-neighbour distance.")

    # ------------------------------------------------------------------ #
    # 3. Observed ANND
    # ------------------------------------------------------------------ #
    real_by_chrom = defaultdict(list)
    for _, row in dnm_dedup.iterrows():
        real_by_chrom[row['chrom']].append(int(row['pos']))
    real_by_chrom = {c: sorted(ps) for c, ps in real_by_chrom.items()}

    observed_annd = avg_nearest_neighbour_dist(real_by_chrom)
    print(f"  Observed ANND = {observed_annd:,.1f} bp", flush=True)

    # ------------------------------------------------------------------ #
    # 4. Permutation simulations
    # ------------------------------------------------------------------ #
    print(f"Running {args.n_simulations} simulations ...", flush=True)
    rng = np.random.default_rng(args.seed)
    sim_annds = []

    for i in range(args.n_simulations):
        sim_pos   = sample_positions(
            chroms_arr, starts_arr, lengths_arr, total_bp, n_dnm, rng
        )
        sim_annds.append(avg_nearest_neighbour_dist(sim_pos))
        if (i + 1) % 100 == 0:
            print(f"  {i + 1}/{args.n_simulations}", flush=True)

    sim_annds = np.array(sim_annds)

    # ------------------------------------------------------------------ #
    # 5. P-values
    # ------------------------------------------------------------------ #
    p_cluster    = float(np.mean(sim_annds <= observed_annd))
    p_dispersed  = float(np.mean(sim_annds >= observed_annd))
    p_two_tailed = 2 * min(p_cluster, p_dispersed)

    sim_mean = float(np.mean(sim_annds))
    sim_std  = float(np.std(sim_annds))
    z_score  = (observed_annd - sim_mean) / sim_std if sim_std > 0 else np.nan

    # ------------------------------------------------------------------ #
    # 6. Write per-simulation TSV
    # ------------------------------------------------------------------ #
    sim_df = pd.DataFrame({
        'simulation': np.arange(1, args.n_simulations + 1),
        'avg_nearest_neighbour_dist_bp': sim_annds,
    })
    sim_df.to_csv(args.output_tsv, sep='\t', index=False)

    # ------------------------------------------------------------------ #
    # 7. Write summary
    # ------------------------------------------------------------------ #
    with open(args.output_summary, 'w') as f:
        f.write(f"Species:                        {args.species_prefix}\n")
        f.write(f"Unique autosomal DNM events:    {n_dnm}\n")
        f.write(f"Callable genome (union):        {total_bp:,} bp\n")
        f.write(f"N simulations:                  {args.n_simulations}\n")
        f.write(f"Random seed:                    {args.seed}\n")
        f.write(f"\nObserved ANND:                  {observed_annd:,.1f} bp\n")
        f.write(f"Simulated ANND mean:            {sim_mean:,.1f} bp\n")
        f.write(f"Simulated ANND std:             {sim_std:,.1f} bp\n")
        f.write(f"Z-score:                        {z_score:.3f}\n")
        f.write(f"\np (clustering,   ANND ≤ obs):  {p_cluster:.4f}\n")
        f.write(f"p (dispersion,   ANND ≥ obs):  {p_dispersed:.4f}\n")
        f.write(f"p (two-tailed):                 {p_two_tailed:.4f}\n")
        if p_two_tailed < 0.05:
            direction = "clustered" if observed_annd < sim_mean else "overdispersed"
            f.write(f"\nConclusion: DNMs are significantly {direction} relative to "
                    f"a uniform callable-genome null (p = {p_two_tailed:.4f}).\n")
        else:
            f.write(f"\nConclusion: No significant deviation from uniform distribution "
                    f"in the callable genome (p = {p_two_tailed:.4f}).\n")

    # ------------------------------------------------------------------ #
    # 8. Plot
    # ------------------------------------------------------------------ #
    fig, ax = plt.subplots(figsize=(8, 4))

    ax.hist(
        sim_annds / 1e6, bins=50,
        color='steelblue', alpha=0.75, edgecolor='white', linewidth=0.4,
        label=f'Simulated (n={args.n_simulations})',
    )
    ax.axvline(
        observed_annd / 1e6,
        color='crimson', linewidth=2,
        label=f'Observed ({observed_annd/1e6:.2f} Mb)',
    )
    ax.axvline(
        sim_mean / 1e6,
        color='steelblue', linewidth=1.5, linestyle='--',
        label=f'Sim. mean ({sim_mean/1e6:.2f} Mb)',
    )

    ax.set_xlabel('Average nearest-neighbour distance (Mb)', fontsize=10)
    ax.set_ylabel('Number of simulations', fontsize=10)
    ax.set_title(
        f"{args.species_prefix.upper()} — DNM spatial distribution test\n"
        f"p (two-tailed) = {p_two_tailed:.4f}   |   Z = {z_score:.2f}",
        fontsize=10, fontweight='bold',
    )
    ax.legend(fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(args.output_plot, bbox_inches='tight')
    plt.close()

    print(f"Observed ANND = {observed_annd/1e6:.3f} Mb  |  "
          f"Sim mean = {sim_mean/1e6:.3f} Mb  |  "
          f"p (two-tailed) = {p_two_tailed:.4f}")
    print(f"Summary: {args.output_summary}")
    print(f"Plot:    {args.output_plot}")


if __name__ == '__main__':
    main()
