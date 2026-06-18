#!/usr/bin/env python3
"""
compute_somatic_spectrum.py
Aggregate the per-species somatic mutation table
(*_final_mutation_table_no_clusters.tsv) into a 96-category SBS spectrum,
matching the schema produced by the legacy somatic project
(species, mutation_type, context, count, fraction).
The somatic table is already strand-collapsed via oriented_ref/oriented_alt/
oriented_context, so no further collapsing is needed.
"""
import argparse
import os

PYRIMIDINE_TYPES = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
BASES = ['A', 'C', 'G', 'T']


def build_sbs96_order():
    cats = []
    for mut_type in PYRIMIDINE_TYPES:
        for five in BASES:
            for three in BASES:
                cats.append((mut_type, f'{five}{mut_type[0]}{three}'))
    return cats


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', required=True,
                   help='Raw somatic mutation table (oriented_ref/alt/context).')
    p.add_argument('--species', required=True)
    p.add_argument('--output', required=True)
    args = p.parse_args()

    counts = {}
    with open(args.input) as f:
        header = f.readline().rstrip('\n').split('\t')
        idx_ref = header.index('oriented_ref')
        idx_alt = header.index('oriented_alt')
        idx_ctx = header.index('oriented_context')
        for line in f:
            fields = line.rstrip('\n').split('\t')
            ref = fields[idx_ref]
            alt = fields[idx_alt]
            ctx = fields[idx_ctx]
            if ref not in ('C', 'T') or len(ctx) != 3 or ctx[1] != ref:
                continue
            mut = f'{ref}>{alt}'
            if mut not in PYRIMIDINE_TYPES:
                continue
            key = (mut, ctx)
            counts[key] = counts.get(key, 0) + 1

    total = sum(counts.values())
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, 'w') as out:
        out.write('species\tmutation_type\tcontext\tcount\tfraction\n')
        for mut, ctx in build_sbs96_order():
            n = counts.get((mut, ctx), 0)
            frac = (n / total) if total > 0 else 0.0
            out.write(f'{args.species}\t{mut}\t{ctx}\t{n}\t{frac:g}\n')

    print(f'{args.species}: {total} somatic SNVs across 96 categories -> {args.output}')


if __name__ == '__main__':
    main()
