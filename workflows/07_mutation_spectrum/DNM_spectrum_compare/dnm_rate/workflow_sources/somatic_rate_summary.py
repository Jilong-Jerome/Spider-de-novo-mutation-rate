#!/usr/bin/env python3
"""
somatic_rate_summary.py
Per-species somatic mutation rate summary across 8 mutation classes
(overall + 7 collapsed: C>A, C>G, C>T_nonCpG, C>T_CpG, T>A, T>C, T>G).

- Numerator: rows from data/dnm/somatic/{SP}_final_mutation_table_no_clusters.tsv
  using (oriented_ref, oriented_alt, oriented_context). CpG split: oriented_context[2] == 'G'.
- Denominator: aggregated trinuc counts from data/callable_sites/somatic/{SP}_trinuc_counts.tsv.
- Rate: per-bp-per-read (per-molecule); no factor-of-2 for diploid.
- 95% CI: exact Poisson via chi2 quantiles divided by the trinuc denominator.
"""
import argparse
import os
from scipy.stats import chi2

CLASSES = ['overall', 'C>A', 'C>G', 'C>T_nonCpG', 'C>T_CpG', 'T>A', 'T>C', 'T>G']


def class_denominators(trinuc_counts):
    overall = sum(trinuc_counts.values())
    c_total = sum(c for t, c in trinuc_counts.items() if t[1] == 'C')
    cpg_total = sum(c for t, c in trinuc_counts.items() if t[1] == 'C' and t[2] == 'G')
    non_cpg_c = c_total - cpg_total
    t_total = sum(c for t, c in trinuc_counts.items() if t[1] == 'T')
    return {
        'overall': overall,
        'C>A': c_total,
        'C>G': c_total,
        'C>T_nonCpG': non_cpg_c,
        'C>T_CpG': cpg_total,
        'T>A': t_total,
        'T>C': t_total,
        'T>G': t_total,
    }


def classify(ref, alt, ctx):
    mut = f'{ref}>{alt}'
    if mut == 'C>T':
        return 'C>T_CpG' if ctx[2] == 'G' else 'C>T_nonCpG'
    return mut


def poisson_ci(n, denom):
    if denom <= 0:
        return float('nan'), float('nan')
    if n == 0:
        ci_low = 0.0
    else:
        ci_low = chi2.ppf(0.025, 2 * n) / 2.0 / denom
    ci_high = chi2.ppf(0.975, 2 * (n + 1)) / 2.0 / denom
    return ci_low, ci_high


def load_trinuc_counts(path):
    counts = {}
    with open(path) as f:
        header = f.readline()
        if not header.startswith('trinuc'):
            f.seek(0)
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('trinuc'):
                continue
            t, c = line.rstrip('\n').split('\t')
            counts[t] = int(c)
    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--species', required=True)
    ap.add_argument('--somatic_dnm', required=True)
    ap.add_argument('--trinuc_counts', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    num = {c: 0 for c in CLASSES}
    n_total = 0
    n_skipped = 0
    with open(args.somatic_dnm) as f:
        header = f.readline().rstrip('\n').split('\t')
        idx_ref = header.index('oriented_ref')
        idx_alt = header.index('oriented_alt')
        idx_ctx = header.index('oriented_context')
        for line in f:
            n_total += 1
            fields = line.rstrip('\n').split('\t')
            ref = fields[idx_ref]
            alt = fields[idx_alt]
            ctx = fields[idx_ctx]
            if ref not in ('C', 'T') or len(ctx) != 3 or ctx[1] != ref or alt not in 'ACGT' or alt == ref:
                n_skipped += 1
                continue
            cls = classify(ref, alt, ctx)
            if cls not in num:
                n_skipped += 1
                continue
            num[cls] += 1
            num['overall'] += 1

    trinuc_counts = load_trinuc_counts(args.trinuc_counts)
    denom = class_denominators(trinuc_counts)

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, 'w') as out:
        out.write(f'# species={args.species} n_total_rows={n_total} '
                  f'n_skipped={n_skipped}\n')
        out.write('species\tmutation_class\tn_mutations\tdenom_trinuc\t'
                  'rate_per_bp_per_read\tci_low\tci_high\n')
        for cls in CLASSES:
            n = num[cls]
            d = denom[cls]
            rate = (n / d) if d > 0 else float('nan')
            ci_low, ci_high = poisson_ci(n, d)
            out.write(
                f'{args.species}\t{cls}\t{n}\t{d}\t{rate:g}\t{ci_low:g}\t{ci_high:g}\n'
            )

    print(f'{args.species}: {num["overall"]} somatic SNVs; '
          f'overall rate={num["overall"]/denom["overall"]:.3e} -> {args.output}')


if __name__ == '__main__':
    main()
