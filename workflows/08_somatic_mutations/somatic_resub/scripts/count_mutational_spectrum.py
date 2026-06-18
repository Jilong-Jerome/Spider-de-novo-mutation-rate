#!/usr/bin/env python3
"""Count somatic mutations per (mutation_type, trinucleotide context) combination.

Reads a filtered somatic mutation TSV and outputs a 96-category mutational
spectrum table (6 mutation types x 16 trinucleotide contexts each).
"""

import argparse
import csv
import itertools
from collections import Counter


# The 6 canonical mutation types (pyrimidine reference)
MUTATION_TYPES = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']

BASES = ['A', 'C', 'G', 'T']


def get_all_contexts():
    """Generate all 96 (mutation_type, context) categories, sorted."""
    categories = []
    for mut_type in MUTATION_TYPES:
        ref = mut_type[0]
        contexts = sorted(
            f"{b5}{ref}{b3}" for b5, b3 in itertools.product(BASES, BASES)
        )
        for ctx in contexts:
            categories.append((mut_type, ctx))
    return categories


def main():
    parser = argparse.ArgumentParser(description='Count mutational spectrum from somatic TSV')
    parser.add_argument('--input', required=True, help='Input TSV file')
    parser.add_argument('--species', required=True, help='Species code')
    parser.add_argument('--output', required=True, help='Output TSV file')
    args = parser.parse_args()

    counts = Counter()
    total = 0

    with open(args.input, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ref = row['oriented_ref']
            alt = row['oriented_alt']
            context = row['oriented_context']
            mut_type = f"{ref}>{alt}"
            counts[(mut_type, context)] += 1
            total += 1

    all_categories = get_all_contexts()

    with open(args.output, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['species', 'mutation_type', 'context', 'count', 'fraction'])
        for mut_type, context in all_categories:
            c = counts.get((mut_type, context), 0)
            fraction = c / total if total > 0 else 0.0
            writer.writerow([args.species, mut_type, context, c, f"{fraction:.6g}"])


if __name__ == '__main__':
    main()
