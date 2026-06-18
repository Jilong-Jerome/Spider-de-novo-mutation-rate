#!/usr/bin/env python3
"""
compute_window_callable.py

For one offspring, parse all autosomal *_callable_sites.vcf files in
chrom_split/, count callable positions per 100 Mb window, and write a TSV.

Output columns: offspring  chrom  window_start  window_end  callable_sites
"""
import argparse
import os
import glob
from collections import defaultdict


def load_fai(fai_path):
    """Return {chrom: length} dict from a .fai file."""
    chrom_lengths = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.strip().split('\t')
            chrom_lengths[parts[0]] = int(parts[1])
    return chrom_lengths


def build_windows(chrom_lengths, x_chroms_set, window_size):
    """
    Return list of (chrom, window_start, window_end) for all autosomal windows.
    Positions are 1-based; window_end = min(window_start + window_size - 1, chrom_len).
    """
    windows = []
    for chrom, length in chrom_lengths.items():
        if chrom in x_chroms_set:
            continue
        start = 1
        while start <= length:
            end = min(start + window_size - 1, length)
            windows.append((chrom, start, end))
            start += window_size
    return windows


def extract_chrom_from_filename(vcf_filename, offspring):
    """
    VCF filename pattern: {offspring}_{chrom}_callable_sites.vcf
    Strip the offspring prefix and '_callable_sites' suffix to get chrom.
    Example: AFR_family1_S1_offspring_afr_1_callable_sites.vcf
             offspring = AFR_family1_S1_offspring
             → chrom = afr_1
    """
    basename = os.path.basename(vcf_filename)
    # Remove offspring prefix (includes trailing underscore)
    after_offspring = basename[len(offspring) + 1:]   # e.g. afr_1_callable_sites.vcf
    chrom = after_offspring.replace('_callable_sites.vcf', '')
    return chrom


def count_callable_sites(vcf_path, window_size):
    """
    Stream a VCF file, count positions per window.
    Returns {window_start: count} for one chromosome.
    """
    counts = defaultdict(int)
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.split('\t', 2)
            if len(fields) < 2:
                continue
            pos = int(fields[1])
            window_idx   = (pos - 1) // window_size
            window_start = window_idx * window_size + 1
            counts[window_start] += 1
    return counts


def main():
    parser = argparse.ArgumentParser(
        description='Count callable sites per window for one offspring.'
    )
    parser.add_argument('--vcf_dir',    required=True,
                        help='Directory containing *_callable_sites.vcf files')
    parser.add_argument('--offspring',  required=True,
                        help='Offspring identifier (e.g. AFR_family1_S1_offspring)')
    parser.add_argument('--genome_fai', required=True,
                        help='Path to genome .fai file')
    parser.add_argument('--window_size', type=int, default=100_000_000,
                        help='Window size in bp (default: 100 Mb)')
    parser.add_argument('--x_chroms',  nargs='+', default=[],
                        help='Sex chromosome names to exclude')
    parser.add_argument('--output',    required=True,
                        help='Output TSV path')
    args = parser.parse_args()

    chrom_lengths = load_fai(args.genome_fai)
    x_chroms_set  = set(args.x_chroms)

    # Build zero-initialised window table keyed by (chrom, window_start)
    windows = build_windows(chrom_lengths, x_chroms_set, args.window_size)
    window_ends = {
        (chrom, ws): we
        for chrom, ws, we in windows
    }
    counts = defaultdict(int)   # (chrom, window_start) → callable_sites

    vcf_files = glob.glob(os.path.join(args.vcf_dir, '*_callable_sites.vcf'))
    for vcf_path in vcf_files:
        chrom = extract_chrom_from_filename(vcf_path, args.offspring)
        if chrom in x_chroms_set:
            continue
        if chrom not in chrom_lengths:
            # Chromosome not in FAI — skip (e.g. unplaced scaffolds)
            continue
        chrom_counts = count_callable_sites(vcf_path, args.window_size)
        for ws, cnt in chrom_counts.items():
            counts[(chrom, ws)] += cnt

    # Write output: all windows, including those with 0 callable sites
    with open(args.output, 'w') as out:
        out.write('offspring\tchrom\twindow_start\twindow_end\tcallable_sites\n')
        for chrom, ws, we in sorted(windows, key=lambda x: (x[0], x[1])):
            cnt = counts.get((chrom, ws), 0)
            out.write(f"{args.offspring}\t{chrom}\t{ws}\t{we}\t{cnt}\n")

    print(f"Done: {args.output}")


if __name__ == '__main__':
    main()
