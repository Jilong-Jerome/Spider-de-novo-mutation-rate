#!/usr/bin/env python3
"""
dnm_distribution_test.py

Chi-squared goodness-of-fit test: are DNMs proportional to callable coverage?

Outputs:
  --output_tsv       per-window results table
  --output_per_ind   per-individual per-window DNM counts (for future exclusion use)
  --output_summary   text summary of the test
"""
import argparse
import sys
import pandas as pd
from scipy.stats import chisquare


def load_window_callable(path):
    df = pd.read_csv(path, sep='\t')
    # Expected columns: chrom, window_start, window_end, total_callable_sites
    return df


def load_dnms(path, x_chroms_set):
    df = pd.read_csv(path, sep='\t')
    # Filter autosomes
    df = df[~df['chrom'].isin(x_chroms_set)].copy()
    return df


def assign_window(pos, window_size):
    idx = (pos - 1) // window_size
    return idx * window_size + 1


def main():
    parser = argparse.ArgumentParser(
        description='Chi-squared GOF test for DNM distribution across windows.'
    )
    parser.add_argument('--window_callable', required=True,
                        help='Species-level window callable TSV')
    parser.add_argument('--dnm_file',   required=True,
                        help='DNM TSV (with chrom, pos, child columns)')
    parser.add_argument('--species',    required=True,
                        help='Species prefix (e.g. afr)')
    parser.add_argument('--window_size', type=int, default=100_000_000)
    parser.add_argument('--x_chroms',   nargs='+', default=[])
    parser.add_argument('--output_tsv',     required=True)
    parser.add_argument('--output_per_ind', required=True)
    parser.add_argument('--output_summary', required=True)
    args = parser.parse_args()

    x_chroms_set = set(args.x_chroms)

    # ------------------------------------------------------------------ #
    # Load data
    # ------------------------------------------------------------------ #
    callable_df = load_window_callable(args.window_callable)
    dnm_df      = load_dnms(args.dnm_file, x_chroms_set)

    if dnm_df.empty:
        sys.exit("No autosomal DNMs found — aborting test.")

    # Assign each DNM to its window
    dnm_df = dnm_df.copy()
    dnm_df['window_start'] = dnm_df['pos'].apply(
        lambda p: assign_window(p, args.window_size)
    )

    # ------------------------------------------------------------------ #
    # Per-individual per-window counts (stored for future use, un-deduplicated)
    # Each child's raw count is preserved here regardless of sharing.
    # ------------------------------------------------------------------ #
    per_ind = (
        dnm_df.groupby(['child', 'chrom', 'window_start'])
        .size()
        .reset_index(name='dnm_count')
    )
    # Attach window_end from callable table
    callable_ends = callable_df[['chrom', 'window_start', 'window_end']].drop_duplicates()
    per_ind = per_ind.merge(callable_ends, on=['chrom', 'window_start'], how='left')
    per_ind = per_ind.rename(columns={'child': 'individual'})
    per_ind = per_ind[['individual', 'chrom', 'window_start', 'window_end', 'dnm_count']]
    per_ind.to_csv(args.output_per_ind, sep='\t', index=False)

    # ------------------------------------------------------------------ #
    # Deduplicate shared DNMs for the species-level test.
    #
    # A DNM shared by multiple siblings (same chrom, pos, ref, alt,
    # father, mother) is one parental mutation event.  Counting it once
    # per child would inflate the observed count in that window and mimic
    # local enrichment.  We collapse to unique mutation events by keeping
    # one row per (chrom, pos, ref, alt, father, mother).
    #
    # Recurrent mutations at the same site in different families are kept
    # as separate events because they represent independent occurrences.
    # ------------------------------------------------------------------ #
    dnm_dedup = dnm_df.drop_duplicates(
        subset=['chrom', 'pos', 'ref', 'alt', 'father', 'mother']
    )
    n_raw   = len(dnm_df)
    n_dedup = len(dnm_dedup)
    n_shared_removed = n_raw - n_dedup

    # ------------------------------------------------------------------ #
    # Species-level observed DNMs per window (deduplicated)
    # ------------------------------------------------------------------ #
    obs = (
        dnm_dedup.groupby(['chrom', 'window_start'])
        .size()
        .reset_index(name='observed_dnms')
    )

    # Left-join callable table with observed DNMs (fill missing with 0)
    merged = callable_df.merge(obs, on=['chrom', 'window_start'], how='left')
    merged['observed_dnms'] = merged['observed_dnms'].fillna(0).astype(int)

    # ------------------------------------------------------------------ #
    # Expected DNMs proportional to callable fraction
    # ------------------------------------------------------------------ #
    total_callable = merged['total_callable_sites'].sum()
    total_dnms     = merged['observed_dnms'].sum()

    merged['callable_fraction'] = merged['total_callable_sites'] / total_callable
    merged['expected_dnms']     = total_dnms * merged['callable_fraction']
    merged['enrichment_ratio']  = merged.apply(
        lambda r: r['observed_dnms'] / r['expected_dnms']
        if r['expected_dnms'] > 0 else float('nan'),
        axis=1,
    )

    # ------------------------------------------------------------------ #
    # Exclude only windows with zero callable sites (expected = 0).
    # These windows are genuinely inaccessible — no DNM can be called
    # there — so observed is also 0 and they contribute nothing to the
    # test.  All other windows, including those with small expected
    # counts, are retained: a window with low expected but high observed
    # is precisely the enrichment signal we want to detect.
    #
    # With this approach sum(expected) == sum(observed) == total_dnms
    # is guaranteed by construction, so scipy.stats.chisquare will not
    # raise a precision error.
    # ------------------------------------------------------------------ #
    testable   = merged[merged['expected_dnms'] > 0].copy()
    n_zero     = (merged['expected_dnms'] == 0).sum()
    n_used     = len(testable)

    if n_used < 2:
        sys.exit(
            f"Only {n_used} windows with callable sites — "
            "cannot perform chi-squared test."
        )

    # ------------------------------------------------------------------ #
    # Chi-squared GOF test
    # ------------------------------------------------------------------ #
    observed_arr = testable['observed_dnms'].values.astype(float)
    expected_arr = testable['expected_dnms'].values.astype(float)

    chi2_stat, p_value = chisquare(f_obs=observed_arr, f_exp=expected_arr)
    dof = n_used - 1

    # ------------------------------------------------------------------ #
    # Write outputs
    # ------------------------------------------------------------------ #
    out_cols = [
        'chrom', 'window_start', 'window_end',
        'total_callable_sites', 'callable_fraction',
        'observed_dnms', 'expected_dnms', 'enrichment_ratio',
    ]
    merged[out_cols].to_csv(args.output_tsv, sep='\t', index=False)

    with open(args.output_summary, 'w') as sumf:
        sumf.write(f"Species:                              {args.species}\n")
        sumf.write(f"Total DNM rows (autosomal):           {n_raw}\n")
        sumf.write(f"Shared DNM rows removed:              {n_shared_removed}\n")
        sumf.write(f"Unique mutation events tested:        {total_dnms}\n")
        sumf.write(f"Windows tested:                       {n_used}\n")
        sumf.write(f"Windows excluded (zero callable):     {n_zero}\n")
        sumf.write(f"Chi-squared statistic:                {chi2_stat:.4f}\n")
        sumf.write(f"Degrees of freedom:                   {dof}\n")
        sumf.write(f"p-value:                              {p_value:.6g}\n")
        sumf.write(
            f"\nInterpretation: "
            f"{'DNMs are uniformly distributed (p > 0.05)' if p_value > 0.05 else 'Significant deviation from uniform distribution (p <= 0.05)'}\n"
        )

    print(f"Chi2={chi2_stat:.4f}  df={dof}  p={p_value:.6g}")
    print(f"Results: {args.output_tsv}")
    print(f"Summary: {args.output_summary}")


if __name__ == '__main__':
    main()
