#!/usr/bin/env python3
"""
merge_window_callable.py

Sum per-offspring window callable counts into a species-level table.

Input:  all {sp}_{offspring}_window_callable.tsv files in input_dir
Output: {sp}_window_callable.tsv
        columns: chrom  window_start  window_end  total_callable_sites
"""
import argparse
import glob
import os
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(
        description='Merge per-offspring window callable counts.'
    )
    parser.add_argument('--input_dir', required=True,
                        help='Directory containing per-offspring window TSV files')
    parser.add_argument('--species',   required=True,
                        help='Species prefix (e.g. afr)')
    parser.add_argument('--output',    required=True,
                        help='Output TSV path')
    args = parser.parse_args()

    pattern = os.path.join(args.input_dir, f"{args.species}_*_window_callable.tsv")
    input_files = glob.glob(pattern)
    if not input_files:
        raise FileNotFoundError(
            f"No per-offspring callable TSV files found matching: {pattern}"
        )

    # Accumulate totals: (chrom, window_start) → total_callable_sites
    totals = defaultdict(int)
    # Keep window_end for output
    window_ends = {}

    for fpath in input_files:
        with open(fpath) as fh:
            header = fh.readline()   # skip header
            for line in fh:
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue
                chrom        = parts[1]
                window_start = int(parts[2])
                window_end   = int(parts[3])
                callable_n   = int(parts[4])
                key = (chrom, window_start)
                totals[key]      += callable_n
                window_ends[key]  = window_end   # same across files

    # Sort by chrom then window_start for readable output
    sorted_keys = sorted(totals.keys(), key=lambda x: (x[0], x[1]))

    with open(args.output, 'w') as out:
        out.write('chrom\twindow_start\twindow_end\ttotal_callable_sites\n')
        for chrom, ws in sorted_keys:
            we  = window_ends[(chrom, ws)]
            cnt = totals[(chrom, ws)]
            out.write(f"{chrom}\t{ws}\t{we}\t{cnt}\n")

    print(f"Merged {len(input_files)} offspring files → {args.output}")
    print(f"Total windows: {len(sorted_keys)}")


if __name__ == '__main__':
    main()
