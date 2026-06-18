#!/usr/bin/env python3
"""
aggregate_callable_trinuc.py
Sum the per-(offspring, chrom) trinuc-count TSVs for one species into a
single species-level autosome callable trinuc table (32 rows).
"""
import argparse
import os

PYRIMIDINE_CANONICAL_TRINUCS = []
for first in 'ACGT':
    for mid in 'CT':
        for last in 'ACGT':
            PYRIMIDINE_CANONICAL_TRINUCS.append(first + mid + last)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--species', required=True)
    ap.add_argument('--input_files', nargs='+', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    totals = {t: 0 for t in PYRIMIDINE_CANONICAL_TRINUCS}
    n_files = 0
    n_offspring = set()
    for path in args.input_files:
        if not os.path.exists(path):
            raise FileNotFoundError(f'Missing per-VCF trinuc TSV: {path}')
        n_files += 1
        n_offspring.add(os.path.basename(os.path.dirname(path)))
        with open(path) as f:
            saw_header = False
            for line in f:
                if line.startswith('#'):
                    continue
                if not saw_header:
                    saw_header = True  # the trinuc_canonical\tcount header line
                    continue
                t, c = line.rstrip('\n').split('\t')
                totals[t] = totals.get(t, 0) + int(c)

    total_bp = sum(totals.values())
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, 'w') as out:
        out.write(
            f'# species={args.species} n_offspring={len(n_offspring)} '
            f'n_offspring_chrom_files={n_files} total_callable_bp={total_bp}\n'
        )
        out.write('trinuc_canonical\tcount\n')
        for t in PYRIMIDINE_CANONICAL_TRINUCS:
            out.write(f'{t}\t{totals[t]}\n')

    print(
        f'{args.species}: aggregated {n_files} per-(offspring, chrom) tables '
        f'across {len(n_offspring)} offspring; total callable bp = {total_bp} '
        f'-> {args.output}'
    )


if __name__ == '__main__':
    main()
