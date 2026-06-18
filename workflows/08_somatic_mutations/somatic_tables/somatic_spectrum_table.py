#!/usr/bin/env python3
"""
somatic_spectrum_table.py
Build a per-species table of somatic mutation counts (overall + all 96 SBS
categories, oriented trinucleotide context, zeros included) together with the
callable sites for each category.

The somatic mutation tables are already strand-collapsed to the pyrimidine
reference via oriented_ref/oriented_alt/oriented_context, so counting is done
directly off those columns (same approach as DNM_spectrum_compare/.../
compute_somatic_spectrum.py). The callable count for an SBS96 category is the
callable opportunity count of its trinucleotide context; the three types that
share a context (e.g. C>A/C>G/C>T on ACA) therefore share the same callable
value. The overall callable count is the sum over all 32 trinucleotides.
"""
import os
import sys

PYRIMIDINE_TYPES = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
BASES = ['A', 'C', 'G', 'T']
SPECIES = ['AFR', 'BIC', 'DUM', 'LIN', 'MIM', 'SAR', 'TEN']

DATA_ROOT = ('/faststorage/project/spider2/social_spiders_2020/people/jilong/'
             'chapter3_dnm/science_resub/DNM_spectrum_compare/data')
SOMATIC_DIR = os.path.join(DATA_ROOT, 'dnm', 'somatic')
CALLABLE_DIR = os.path.join(DATA_ROOT, 'callable_sites', 'somatic')

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_TSV = os.path.join(OUT_DIR, 'somatic_spectrum_counts.tsv')


def build_sbs96_order():
    """Canonical SBS96 order: (mutation_type, context, sbs96_category)."""
    cats = []
    for mut in PYRIMIDINE_TYPES:
        ref = mut[0]
        for five in BASES:
            for three in BASES:
                ctx = f'{five}{ref}{three}'
                cat = f'{five}[{mut}]{three}'
                cats.append((mut, ctx, cat))
    return cats


def read_callable(species):
    """Return {trinuc_canonical: count} for a species' callable trinuc file."""
    path = os.path.join(CALLABLE_DIR, f'{species}_trinuc_counts.tsv')
    callable_counts = {}
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        i_tri = header.index('trinuc_canonical')
        i_cnt = header.index('count')
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) <= max(i_tri, i_cnt):
                continue
            callable_counts[fields[i_tri]] = int(fields[i_cnt])
    return callable_counts


def count_somatic(species):
    """Return ({(mut, ctx): count}, total_kept, dropped) for a species."""
    path = os.path.join(SOMATIC_DIR, f'{species}_final_mutation_table_no_clusters.tsv')
    counts = {}
    total = 0
    dropped = 0
    with open(path) as f:
        header = f.readline().rstrip('\n').split('\t')
        i_ref = header.index('oriented_ref')
        i_alt = header.index('oriented_alt')
        i_ctx = header.index('oriented_context')
        for line in f:
            fields = line.rstrip('\n').split('\t')
            ref = fields[i_ref]
            alt = fields[i_alt]
            ctx = fields[i_ctx]
            mut = f'{ref}>{alt}'
            if (ref not in ('C', 'T') or len(ctx) != 3 or ctx[1] != ref
                    or mut not in PYRIMIDINE_TYPES):
                dropped += 1
                continue
            key = (mut, ctx)
            counts[key] = counts.get(key, 0) + 1
            total += 1
    return counts, total, dropped


def main():
    sbs96 = build_sbs96_order()
    with open(OUT_TSV, 'w') as out:
        out.write('species\tmutation_type\tcontext\tsbs96_category\t'
                  'count\tcallable_sites\n')
        for sp in SPECIES:
            counts, total, dropped = count_somatic(sp)
            callable_counts = read_callable(sp)
            total_callable = sum(callable_counts.values())

            # Overall row first.
            out.write(f'{sp}\toverall\tNA\toverall\t{total}\t{total_callable}\n')

            # 96 SBS rows in canonical order (zeros included).
            for mut, ctx, cat in sbs96:
                n = counts.get((mut, ctx), 0)
                callable_n = callable_counts.get(ctx, 0)
                out.write(f'{sp}\t{mut}\t{ctx}\t{cat}\t{n}\t{callable_n}\n')

            print(f'{sp}: {total} somatic SNVs kept, {dropped} rows dropped; '
                  f'total callable={total_callable}', file=sys.stderr)

    print(f'Wrote {OUT_TSV}', file=sys.stderr)


if __name__ == '__main__':
    main()
