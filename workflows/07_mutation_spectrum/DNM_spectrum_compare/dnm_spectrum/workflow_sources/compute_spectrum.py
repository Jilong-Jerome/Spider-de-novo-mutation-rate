#!/usr/bin/env python3
"""
compute_spectrum.py
Merge annotated DNMs from multiple species, filter by chromosome group,
deduplicate shared DNMs, and count SBS96 categories.
"""
import argparse
import os
from itertools import product

PYRIMIDINE_TYPES = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
BASES = ['A', 'C', 'G', 'T']


def build_sbs96_order():
    """Return the canonical SBS96 category list in standard order."""
    categories = []
    for mut_type in PYRIMIDINE_TYPES:
        for five_prime in BASES:
            for three_prime in BASES:
                cat = f'{five_prime}[{mut_type}]{three_prime}'
                categories.append(cat)
    return categories


def main():
    parser = argparse.ArgumentParser(description='Compute SBS96 spectrum from annotated DNMs')
    parser.add_argument('--input_files', nargs='+', required=True,
                        help='Annotated DNM TSV files (from extract_context.py)')
    parser.add_argument('--chrom_group', required=True, choices=['autosome', 'x_chromosome'],
                        help='Chromosome group to include')
    parser.add_argument('--output', required=True, help='Output spectrum TSV')
    args = parser.parse_args()

    # Read all annotated DNMs
    all_rows = []
    for fpath in args.input_files:
        if not os.path.exists(fpath):
            continue
        with open(fpath) as f:
            header = f.readline().strip().split('\t')
            for line in f:
                fields = line.strip().split('\t')
                row = dict(zip(header, fields))
                all_rows.append(row)

    # Filter by chromosome group
    rows = [r for r in all_rows if r['chrom_group'] == args.chrom_group]

    # Deduplicate shared DNMs: same (chrom, pos, ref, alt, father, mother)
    seen = set()
    unique_rows = []
    for r in rows:
        key = (r['chrom'], r['pos'], r['ref'], r['alt'], r['father'], r['mother'])
        if key not in seen:
            seen.add(key)
            unique_rows.append(r)

    # Count SBS96 categories
    counts = {}
    for r in unique_rows:
        cat = r['sbs96_category']
        counts[cat] = counts.get(cat, 0) + 1

    # Write output with all 96 categories (including zeros)
    sbs96_order = build_sbs96_order()
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, 'w') as out:
        out.write('mutation_type\ttrinuc_context\tsbs96_category\tcount\n')
        for cat in sbs96_order:
            mut_type = cat[2:5]  # e.g. C>A from A[C>A]T
            trinuc = cat[0] + cat[2] + cat[6]  # e.g. ACT from A[C>A]T
            out.write(f'{mut_type}\t{trinuc}\t{cat}\t{counts.get(cat, 0)}\n')

    total = sum(counts.get(cat, 0) for cat in sbs96_order)
    print(f"Spectrum: {total} unique DNMs across {len(sbs96_order)} categories -> {args.output}")


if __name__ == '__main__':
    main()
