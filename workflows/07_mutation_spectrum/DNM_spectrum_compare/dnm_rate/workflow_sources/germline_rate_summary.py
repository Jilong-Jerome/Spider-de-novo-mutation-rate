#!/usr/bin/env python3
"""
germline_rate_summary.py
Per-species germline mutation rate summary across 8 mutation classes
(overall + 7 collapsed: C>A, C>G, C>T_nonCpG, C>T_CpG, T>A, T>C, T>G).

- Numerator: autosomal germline DNMs after offspring exclusions and dedup.
- Denominator: class-conditioned bp from the species autosome callable trinuc
  table (built exactly from the per-(offspring, chrom) callable VCFs).
- Rate: per-haploid-bp-per-generation = N / (2 * denom).
- 95% CI: exact Poisson via chi2 quantiles, divided by (2 * denom).
"""
import argparse
import os
from scipy.stats import chi2

CLASSES = ['overall', 'C>A', 'C>G', 'C>T_nonCpG', 'C>T_CpG', 'T>A', 'T>C', 'T>G']
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def class_denominators(trinuc_counts):
    """Map 32 pyrimidine-canonical trinuc counts -> per-class denominators."""
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


def reverse_complement(seq):
    return ''.join(COMPLEMENT[b] for b in reversed(seq.upper()))


def classify(ref, alt, trinuc):
    """Classify a DNM from raw alleles and local context using SBS strand collapse."""
    ref = ref.upper()
    alt = alt.upper()
    trinuc = trinuc.upper()
    if ref not in COMPLEMENT or alt not in COMPLEMENT or len(trinuc) != 3:
        return None
    if trinuc[1] != ref:
        return None
    if ref in ('C', 'T'):
        collapsed_trinuc = trinuc
        mut_type = f'{ref}>{alt}'
    else:
        collapsed_trinuc = reverse_complement(trinuc)
        mut_type = f'{COMPLEMENT[ref]}>{COMPLEMENT[alt]}'

    if mut_type == 'C>T':
        return 'C>T_CpG' if collapsed_trinuc[2] == 'G' else 'C>T_nonCpG'
    return mut_type


def poisson_ci(n, denom_haploid):
    """Exact Poisson 95% CI on rate = n / denom_haploid (denom already 2*bp)."""
    if denom_haploid <= 0:
        return float('nan'), float('nan')
    if n == 0:
        ci_low = 0.0
    else:
        ci_low = chi2.ppf(0.025, 2 * n) / 2.0 / denom_haploid
    ci_high = chi2.ppf(0.975, 2 * (n + 1)) / 2.0 / denom_haploid
    return ci_low, ci_high


def load_callable_offspring(callable_file):
    """Return set of offspring_id parsed from column 3 of the callable TSV."""
    out = set()
    with open(callable_file) as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 3:
                continue
            chrom = fields[0]
            tag = fields[2]
            suffix = f'_{chrom}'
            if not tag.endswith(suffix):
                continue
            out.add(tag[:-len(suffix)])
    return out


def load_callable_trinuc(path):
    counts = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('trinuc_canonical'):
                continue
            t, c = line.rstrip('\n').split('\t')
            counts[t] = int(c)
    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--species', required=True)
    ap.add_argument('--dnm_file', required=True)
    ap.add_argument('--callable_file', required=True)
    ap.add_argument('--callable_trinuc', required=True,
                    help='Species autosome 32-trinuc TSV (aggregate output).')
    ap.add_argument('--x_chromosomes', nargs='*', default=[])
    ap.add_argument('--exclude_trios', nargs='*', default=[])
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    x_set = set(args.x_chromosomes)
    excl_set = set(args.exclude_trios)

    callable_off = load_callable_offspring(args.callable_file)
    kept_off = callable_off - excl_set

    # Load DNMs and apply filters: drop X-chrom, drop offspring not kept, dedup.
    dnms = []
    n_total_rows = 0
    n_dropped_x = 0
    n_dropped_off = 0
    n_skipped_unclassified = 0
    seen = set()
    with open(args.dnm_file) as f:
        header = f.readline().rstrip('\n').split('\t')
        idx = {col: i for i, col in enumerate(header)}
        for line in f:
            n_total_rows += 1
            fields = line.rstrip('\n').split('\t')
            chrom = fields[idx['chrom']]
            if chrom in x_set:
                n_dropped_x += 1
                continue
            child = fields[idx['child']]
            if child not in kept_off:
                n_dropped_off += 1
                continue
            father = fields[idx['father']]
            mother = fields[idx['mother']]
            pos = fields[idx['pos']]
            ref = fields[idx['ref']]
            alt = fields[idx['alt']]
            trinuc = fields[idx['n_before']] + fields[idx['nuc_ref']] + fields[idx['nuc_after']]
            cls = classify(ref, alt, trinuc)
            if cls not in CLASSES:
                n_skipped_unclassified += 1
                continue
            key = (chrom, pos, ref, alt, father, mother)
            if key in seen:
                continue
            seen.add(key)
            dnms.append((chrom, child, cls))

    n_unique = len(dnms)

    # Numerator per class
    num = {c: 0 for c in CLASSES}
    for chrom, child, cls in dnms:
        num[cls] += 1
        num['overall'] += 1

    # Denominator per class (haploid: multiply bp by 2, then divide for rate)
    trinuc_counts = load_callable_trinuc(args.callable_trinuc)
    bp_denom = class_denominators(trinuc_counts)  # diploid bp totals
    haploid_denom = {c: 2 * bp_denom[c] for c in CLASSES}

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, 'w') as out:
        out.write(f'# species={args.species} n_offspring_kept={len(kept_off)} '
                  f'n_dnm_total_rows={n_total_rows} n_dropped_x={n_dropped_x} '
                  f'n_dropped_offspring_not_in_callable={n_dropped_off} '
                  f'n_skipped_unclassified={n_skipped_unclassified} '
                  f'n_unique_autosomal_dnm={n_unique}\n')
        out.write('species\tmutation_class\tn_mutations\tdenom_bp\t'
                  'rate_per_haploid_bp_per_gen\tci_low\tci_high\tn_offspring_kept\n')
        for cls in CLASSES:
            n = num[cls]
            denom_hap = haploid_denom[cls]
            rate = (n / denom_hap) if denom_hap > 0 else float('nan')
            ci_low, ci_high = poisson_ci(n, denom_hap)
            out.write(
                f'{args.species}\t{cls}\t{n}\t{bp_denom[cls]}\t'
                f'{rate:g}\t{ci_low:g}\t{ci_high:g}\t{len(kept_off)}\n'
            )

    print(f'{args.species}: {n_unique} unique autosomal germline DNMs; '
          f'kept_offspring={len(kept_off)}; '
          f'overall rate={num["overall"]/haploid_denom["overall"]:.3e} -> {args.output}')


if __name__ == '__main__':
    main()
