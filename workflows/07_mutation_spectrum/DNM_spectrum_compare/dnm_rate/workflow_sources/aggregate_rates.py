#!/usr/bin/env python3
"""
aggregate_rates.py
Combine the per-species germline + somatic rate summaries plus the
per-(offspring, chrom) callable trinuc tables into:
  - mutation_rate_master.tsv               (per-species + per-group rates)
  - per_trio_germline_mutations_callable.tsv  (long: species x offspring x chrom x class)
  - per_species_somatic_mutations_callable.tsv (long: species x class)
  - somatic_vs_germline_fold_change.tsv    (per-species + per-group fold change)
  - social_vs_subsocial_fold_change_test.tsv
  - mutation_rate_96_context_master.tsv    (per-species + per-group 96 SBS rates)
  - somatic_vs_germline_96_context_fold_change.tsv
"""
import argparse
import math
import os
from collections import defaultdict
from scipy.stats import chi2
from scipy.optimize import minimize
import yaml

CLASSES = ['overall', 'C>A', 'C>G', 'C>T_nonCpG', 'C>T_CpG', 'T>A', 'T>C', 'T>G']
MUT_TYPES = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
BASES = ['A', 'C', 'G', 'T']
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


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


def reverse_complement(seq):
    return ''.join(COMPLEMENT[b] for b in reversed(seq.upper()))


def classify_germline(ref, alt, trinuc):
    ref = ref.upper()
    alt = alt.upper()
    trinuc = trinuc.upper()
    if ref not in COMPLEMENT or alt not in COMPLEMENT or len(trinuc) != 3:
        return None
    if trinuc[1] != ref:
        return None
    if ref in ('C', 'T'):
        collapsed_trinuc = trinuc
        mut = f'{ref}>{alt}'
    else:
        collapsed_trinuc = reverse_complement(trinuc)
        mut = f'{COMPLEMENT[ref]}>{COMPLEMENT[alt]}'
    if mut == 'C>T':
        return 'C>T_CpG' if collapsed_trinuc[2] == 'G' else 'C>T_nonCpG'
    return mut


def classify_germline_sbs96(ref, alt, trinuc):
    ref = ref.upper()
    alt = alt.upper()
    trinuc = trinuc.upper()
    if ref not in COMPLEMENT or alt not in COMPLEMENT or len(trinuc) != 3:
        return None
    if trinuc[1] != ref:
        return None
    if ref in ('C', 'T'):
        collapsed_trinuc = trinuc
        mut = f'{ref}>{alt}'
    else:
        collapsed_trinuc = reverse_complement(trinuc)
        mut = f'{COMPLEMENT[ref]}>{COMPLEMENT[alt]}'
    if mut not in MUT_TYPES:
        return None
    return mut, collapsed_trinuc, f'{collapsed_trinuc[0]}[{mut}]{collapsed_trinuc[2]}'


def classify_somatic(ref, alt, ctx):
    mut = f'{ref}>{alt}'
    if mut == 'C>T':
        return 'C>T_CpG' if ctx[2] == 'G' else 'C>T_nonCpG'
    return mut


def classify_somatic_sbs96(ref, alt, ctx):
    ref = ref.upper()
    alt = alt.upper()
    ctx = ctx.upper()
    if ref not in ('C', 'T') or alt not in COMPLEMENT or len(ctx) != 3 or ctx[1] != ref:
        return None
    mut = f'{ref}>{alt}'
    if mut not in MUT_TYPES:
        return None
    return mut, ctx, f'{ctx[0]}[{mut}]{ctx[2]}'


def build_sbs96_order():
    cats = []
    for mut in MUT_TYPES:
        ref = mut[0]
        for five in BASES:
            for three in BASES:
                trinuc = f'{five}{ref}{three}'
                cats.append((mut, trinuc, f'{five}[{mut}]{three}'))
    return cats


def poisson_ci(n, denom):
    if denom <= 0:
        return float('nan'), float('nan')
    if n == 0:
        ci_low = 0.0
    else:
        ci_low = chi2.ppf(0.025, 2 * n) / 2.0 / denom
    ci_high = chi2.ppf(0.975, 2 * (n + 1)) / 2.0 / denom
    return ci_low, ci_high


def bh_adjust(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    adj = [0.0] * n
    running = 1.0
    for rank_from_end, i in enumerate(reversed(order), start=1):
        rank = n - rank_from_end + 1
        val = min(pvals[i] * n / rank, running)
        running = val
        adj[i] = min(max(val, 0.0), 1.0)
    return adj


def poisson_loglik(counts, means):
    ll = 0.0
    for y, mu in zip(counts, means):
        if mu <= 0:
            return float('-inf')
        ll += y * math.log(mu) - mu
    return ll


def poisson_interaction_lrt(counts, exposures):
    """
    Test group:level interaction in a 2x2 Poisson rate table with offsets.
    Cell order: social germline, social somatic, subsocial germline,
    subsocial somatic.
    """
    if any(e <= 0 for e in exposures):
        return float('nan'), float('nan')

    full_means = [y if y > 0 else 1e-300 for y in counts]
    ll_full = poisson_loglik(counts, full_means)

    def neg_ll(beta):
        b0, b_group, b_level = beta
        eta = [
            b0,
            b0 + b_level,
            b0 + b_group,
            b0 + b_group + b_level,
        ]
        means = [exposures[i] * math.exp(eta[i]) for i in range(4)]
        return -poisson_loglik(counts, means)

    total_rate = (sum(counts) + 0.5) / sum(exposures)
    init = [math.log(total_rate), 0.0, 0.0]
    res = minimize(neg_ll, init, method='BFGS')
    if not res.success:
        res = minimize(neg_ll, init, method='Nelder-Mead')
    ll_null = -res.fun
    stat = max(2 * (ll_full - ll_null), 0.0)
    p = 1 - chi2.cdf(stat, 1)
    return stat, p


def load_trinuc_table(path):
    counts = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('trinuc'):
                continue
            t, c = line.rstrip('\n').split('\t')
            counts[t] = int(c)
    return counts


def load_callable_offspring(callable_file):
    out = set()
    with open(callable_file) as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 3:
                continue
            chrom = fields[0]
            tag = fields[2]
            suffix = f'_{chrom}'
            if tag.endswith(suffix):
                out.add(tag[:-len(suffix)])
    return out


def count_germline_sbs96(dnm_file, callable_file, x_chromosomes, exclude_trios):
    """Count unique autosomal germline DNMs by folded SBS96 category."""
    counts = defaultdict(int)
    x_set = set(x_chromosomes or [])
    kept_off = load_callable_offspring(callable_file) - set(exclude_trios or [])
    seen = set()
    with open(dnm_file) as f:
        header = f.readline().rstrip('\n').split('\t')
        idx = {col: i for i, col in enumerate(header)}
        for line in f:
            fields = line.rstrip('\n').split('\t')
            chrom = fields[idx['chrom']]
            child = fields[idx['child']]
            if chrom in x_set or child not in kept_off:
                continue
            ref = fields[idx['ref']]
            alt = fields[idx['alt']]
            father = fields[idx['father']]
            mother = fields[idx['mother']]
            pos = fields[idx['pos']]
            trinuc = fields[idx['n_before']] + fields[idx['nuc_ref']] + fields[idx['nuc_after']]
            cat = classify_germline_sbs96(ref, alt, trinuc)
            if cat is None:
                continue
            key = (chrom, pos, ref, alt, father, mother)
            if key in seen:
                continue
            seen.add(key)
            counts[cat] += 1
    return counts


def count_somatic_sbs96(somatic_dnm_file):
    counts = defaultdict(int)
    with open(somatic_dnm_file) as f:
        header = f.readline().rstrip('\n').split('\t')
        idx_ref = header.index('oriented_ref')
        idx_alt = header.index('oriented_alt')
        idx_ctx = header.index('oriented_context')
        for line in f:
            fields = line.rstrip('\n').split('\t')
            cat = classify_somatic_sbs96(fields[idx_ref], fields[idx_alt], fields[idx_ctx])
            if cat is not None:
                counts[cat] += 1
    return counts


def load_callable_offspring_chrom_bp(callable_file):
    """Map (offspring_id, chrom) -> callable_bp from the species callable TSV."""
    out = {}
    with open(callable_file) as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 3:
                continue
            chrom = fields[0]
            try:
                bp = int(fields[1])
            except ValueError:
                continue
            tag = fields[2]
            suffix = f'_{chrom}'
            if not tag.endswith(suffix):
                continue
            offspring = tag[:-len(suffix)]
            out[(offspring, chrom)] = bp
    return out


def load_germline_summary(path):
    """Return list of dict rows from a per-species germline rate summary."""
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('species\t'):
                header = line.rstrip('\n').split('\t')
                continue
            fields = line.rstrip('\n').split('\t')
            rows.append(dict(zip(header, fields)))
    return rows


def load_somatic_summary(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('species\t'):
                header = line.rstrip('\n').split('\t')
                continue
            fields = line.rstrip('\n').split('\t')
            rows.append(dict(zip(header, fields)))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--spectrum_config', required=True)
    ap.add_argument('--output_path', required=True)
    ap.add_argument('--output_master', required=True)
    ap.add_argument('--output_per_trio', required=True)
    ap.add_argument('--output_per_species_somatic', required=True)
    ap.add_argument('--output_fold_change', required=True)
    ap.add_argument('--output_fold_change_test', required=True)
    ap.add_argument('--output_master_96', required=True)
    ap.add_argument('--output_fold_change_96', required=True)
    args = ap.parse_args()

    with open(args.spectrum_config) as f:
        cfg = yaml.safe_load(f)
    species_configs = cfg['species']
    spectrum_groups = cfg['spectrum_groups']
    species_list = list(species_configs.keys())

    # ------------------------------------------------------------------
    # Per-species summaries: load germline + somatic
    # ------------------------------------------------------------------
    germline_rows = {}
    somatic_rows = {}
    for SP in species_list:
        gpath = f'{args.output_path}/rate/germline/{SP}/{SP}_germline_rate_summary.tsv'
        spath = f'{args.output_path}/rate/somatic/{SP}/{SP}_somatic_rate_summary.tsv'
        if os.path.exists(gpath):
            germline_rows[SP] = {r['mutation_class']: r for r in load_germline_summary(gpath)}
        if os.path.exists(spath):
            somatic_rows[SP] = {r['mutation_class']: r for r in load_somatic_summary(spath)}

    # ------------------------------------------------------------------
    # Per-species folded SBS96 numerators + denominators
    # ------------------------------------------------------------------
    sbs96_order = build_sbs96_order()
    germline_96 = {}
    somatic_96 = {}
    germline_trinuc = {}
    somatic_trinuc = {}
    for SP in species_list:
        sp_conf = species_configs[SP]
        dnm_file = sp_conf.get('dnm_file')
        callable_file = sp_conf.get('germline_callable_file')
        gtri_path = f'{args.output_path}/rate/germline_callable_trinuc/{SP}/{SP}_autosome_callable_trinuc.tsv'
        if dnm_file and callable_file and os.path.exists(dnm_file) and os.path.exists(callable_file):
            germline_96[SP] = count_germline_sbs96(
                dnm_file,
                callable_file,
                sp_conf.get('x_chromosomes', []),
                sp_conf.get('exclude_trios', []),
            )
        if os.path.exists(gtri_path):
            germline_trinuc[SP] = load_trinuc_table(gtri_path)

        somatic_dnm_file = sp_conf.get('somatic_dnm_file')
        somatic_trinuc_file = sp_conf.get('somatic_trinuc_counts_file')
        if somatic_dnm_file and os.path.exists(somatic_dnm_file):
            somatic_96[SP] = count_somatic_sbs96(somatic_dnm_file)
        if somatic_trinuc_file and os.path.exists(somatic_trinuc_file):
            somatic_trinuc[SP] = load_trinuc_table(somatic_trinuc_file)

    # ------------------------------------------------------------------
    # MASTER rate TSV: per-species + per-group rows for both levels
    # ------------------------------------------------------------------
    master_lines = ['group\tspecies\tlevel\tmutation_class\tn_mutations\t'
                    'denom\trate\tci_low\tci_high']

    for SP in species_list:
        if SP in germline_rows:
            for cls in CLASSES:
                r = germline_rows[SP][cls]
                n = int(r['n_mutations'])
                bp = int(r['denom_bp'])
                denom_hap = 2 * bp
                rate = (n / denom_hap) if denom_hap > 0 else float('nan')
                lo, hi = poisson_ci(n, denom_hap)
                master_lines.append(
                    f'per_species\t{SP}\tgermline\t{cls}\t{n}\t{bp}\t{rate:g}\t{lo:g}\t{hi:g}'
                )
        if SP in somatic_rows:
            for cls in CLASSES:
                r = somatic_rows[SP][cls]
                n = int(r['n_mutations'])
                d = int(r['denom_trinuc'])
                rate = (n / d) if d > 0 else float('nan')
                lo, hi = poisson_ci(n, d)
                master_lines.append(
                    f'per_species\t{SP}\tsomatic\t{cls}\t{n}\t{d}\t{rate:g}\t{lo:g}\t{hi:g}'
                )

    # Group rows (sum N + sum denom across species)
    for group_name, group_conf in spectrum_groups.items():
        for cls in CLASSES:
            # Germline aggregate across group
            n_g = 0
            bp_g = 0
            n_in_g = 0
            for SP in group_conf['species']:
                if SP not in germline_rows:
                    continue
                r = germline_rows[SP][cls]
                n_g += int(r['n_mutations'])
                bp_g += int(r['denom_bp'])
                n_in_g += 1
            if n_in_g > 0:
                denom_hap = 2 * bp_g
                rate = (n_g / denom_hap) if denom_hap > 0 else float('nan')
                lo, hi = poisson_ci(n_g, denom_hap)
                master_lines.append(
                    f'{group_name}\tmerged\tgermline\t{cls}\t{n_g}\t{bp_g}\t{rate:g}\t{lo:g}\t{hi:g}'
                )
            # Somatic aggregate
            n_s = 0
            d_s = 0
            n_in_s = 0
            for SP in group_conf['species']:
                if SP not in somatic_rows:
                    continue
                r = somatic_rows[SP][cls]
                n_s += int(r['n_mutations'])
                d_s += int(r['denom_trinuc'])
                n_in_s += 1
            if n_in_s > 0:
                rate = (n_s / d_s) if d_s > 0 else float('nan')
                lo, hi = poisson_ci(n_s, d_s)
                master_lines.append(
                    f'{group_name}\tmerged\tsomatic\t{cls}\t{n_s}\t{d_s}\t{rate:g}\t{lo:g}\t{hi:g}'
                )

    os.makedirs(os.path.dirname(args.output_master), exist_ok=True)
    with open(args.output_master, 'w') as f:
        f.write('\n'.join(master_lines) + '\n')

    # ------------------------------------------------------------------
    # PER-SPECIES SOMATIC TSV
    # ------------------------------------------------------------------
    pss_lines = ['species\tmutation_class\tn_somatic_mutations\tcallable_trinuc_count']
    for SP in species_list:
        if SP not in somatic_rows:
            continue
        for cls in CLASSES:
            r = somatic_rows[SP][cls]
            pss_lines.append(
                f'{SP}\t{cls}\t{r["n_mutations"]}\t{r["denom_trinuc"]}'
            )
    with open(args.output_per_species_somatic, 'w') as f:
        f.write('\n'.join(pss_lines) + '\n')

    # ------------------------------------------------------------------
    # PER-TRIO PER-CHROM GERMLINE TSV
    # ------------------------------------------------------------------
    per_trio_lines = ['# Long-format per-(species, offspring, chrom, mutation_class)']
    per_trio_lines.append(
        'species\toffspring_id\tchrom\tchrom_group\tmutation_class\t'
        'n_germline_mutations\tcallable_denom_bp\tn_callable_total_bp\tincluded'
    )

    for SP in species_list:
        sp_conf = species_configs[SP]
        x_set = set(sp_conf.get('x_chromosomes', []))
        excl_set = set(sp_conf.get('exclude_trios', []))
        callable_file = sp_conf.get('germline_callable_file')
        dnm_file = sp_conf.get('dnm_file')
        if not callable_file or not os.path.exists(callable_file):
            continue
        if not dnm_file or not os.path.exists(dnm_file):
            continue

        # Build the universe of (offspring, chrom) pairs from the callable TSV
        offchrom_bp = load_callable_offspring_chrom_bp(callable_file)
        # Per-trio per-class DNM counts (no dedup)
        dnm_counts = defaultdict(lambda: defaultdict(int))  # (off, chrom) -> {class: n}
        with open(dnm_file) as f:
            header = f.readline().rstrip('\n').split('\t')
            idx = {col: i for i, col in enumerate(header)}
            for line in f:
                fields = line.rstrip('\n').split('\t')
                chrom = fields[idx['chrom']]
                child = fields[idx['child']]
                ref = fields[idx['ref']]
                alt = fields[idx['alt']]
                trinuc = fields[idx['n_before']] + fields[idx['nuc_ref']] + fields[idx['nuc_after']]
                cls = classify_germline(ref, alt, trinuc)
                if cls in CLASSES:
                    dnm_counts[(child, chrom)][cls] += 1
                    dnm_counts[(child, chrom)]['overall'] += 1

        # Class-conditioned denominators per (offspring, chrom) — autosomes only
        per_trio_chrom_denom = {}
        for (offspring, chrom), bp in offchrom_bp.items():
            if chrom in x_set:
                continue
            tri_path = (
                f'{args.output_path}/rate/germline_callable_trinuc/{SP}/'
                f'{offspring}/{offspring}_{chrom}_trinuc.tsv'
            )
            if not os.path.exists(tri_path):
                continue
            tcounts = load_trinuc_table(tri_path)
            per_trio_chrom_denom[(offspring, chrom)] = class_denominators(tcounts)

        # Emit one row per (offspring, chrom, class). Include excluded trios
        # as well so the user has full provenance.
        for (offspring, chrom), bp in sorted(offchrom_bp.items()):
            chrom_group = 'x_chromosome' if chrom in x_set else 'autosome'
            included = 'no' if offspring in excl_set else 'yes'
            cls_denom = per_trio_chrom_denom.get((offspring, chrom))
            counts = dnm_counts.get((offspring, chrom), {})
            for cls in CLASSES:
                n = counts.get(cls, 0)
                if cls_denom is not None:
                    denom_bp = cls_denom[cls]
                elif cls == 'overall' and chrom_group == 'x_chromosome':
                    # X chrom: report total bp under "overall"; class-specific NA
                    denom_bp = bp
                else:
                    denom_bp = ''
                per_trio_lines.append(
                    f'{SP}\t{offspring}\t{chrom}\t{chrom_group}\t{cls}\t{n}\t'
                    f'{denom_bp}\t{bp}\t{included}'
                )

    with open(args.output_per_trio, 'w') as f:
        f.write('\n'.join(per_trio_lines) + '\n')

    # ------------------------------------------------------------------
    # FOLD CHANGE TSV (per-species + per-group)
    # ------------------------------------------------------------------
    fc_lines = ['group\tspecies\tmutation_class\tgermline_rate\tsomatic_rate\t'
                'fold_change\tlog2_fold_change\tfc_ci_low\tfc_ci_high']

    def emit_fc(group_label, species_label, cls, n_g, bp_g, n_s, d_s):
        denom_hap = 2 * bp_g
        if denom_hap <= 0 or d_s <= 0 or n_g == 0 or n_s == 0:
            mu_g = (n_g / denom_hap) if denom_hap > 0 else float('nan')
            mu_s = (n_s / d_s) if d_s > 0 else float('nan')
            fc = (mu_s / mu_g) if (mu_g and mu_g > 0) else float('nan')
            l2 = math.log2(fc) if fc and fc > 0 else float('nan')
            fc_lines.append(
                f'{group_label}\t{species_label}\t{cls}\t{mu_g:g}\t{mu_s:g}\t'
                f'{fc:g}\t{l2:g}\t\t'
            )
            return
        mu_g = n_g / denom_hap
        mu_s = n_s / d_s
        fc = mu_s / mu_g
        log_fc = math.log(fc)
        se = math.sqrt(1.0 / n_g + 1.0 / n_s)
        lo = math.exp(log_fc - 1.96 * se)
        hi = math.exp(log_fc + 1.96 * se)
        fc_lines.append(
            f'{group_label}\t{species_label}\t{cls}\t{mu_g:g}\t{mu_s:g}\t'
            f'{fc:g}\t{log_fc / math.log(2):g}\t{lo:g}\t{hi:g}'
        )

    for SP in species_list:
        if SP not in germline_rows or SP not in somatic_rows:
            continue
        for cls in CLASSES:
            g = germline_rows[SP][cls]
            s = somatic_rows[SP][cls]
            emit_fc('per_species', SP, cls,
                    int(g['n_mutations']), int(g['denom_bp']),
                    int(s['n_mutations']), int(s['denom_trinuc']))

    for group_name, group_conf in spectrum_groups.items():
        for cls in CLASSES:
            n_g = bp_g = n_s = d_s = 0
            n_in = 0
            for SP in group_conf['species']:
                if SP not in germline_rows or SP not in somatic_rows:
                    continue
                n_g += int(germline_rows[SP][cls]['n_mutations'])
                bp_g += int(germline_rows[SP][cls]['denom_bp'])
                n_s += int(somatic_rows[SP][cls]['n_mutations'])
                d_s += int(somatic_rows[SP][cls]['denom_trinuc'])
                n_in += 1
            if n_in > 0:
                emit_fc(group_name, 'merged', cls, n_g, bp_g, n_s, d_s)

    with open(args.output_fold_change, 'w') as f:
        f.write('\n'.join(fc_lines) + '\n')

    # ------------------------------------------------------------------
    # SOCIAL VS SUBSOCIAL FOLD-CHANGE INTERACTION TEST
    # ------------------------------------------------------------------
    fc_test_rows = []
    for cls in CLASSES:
        group_cells = {}
        for group_name in ('social', 'subsocial'):
            group_conf = spectrum_groups.get(group_name)
            if not group_conf:
                continue
            n_g = bp_g = n_s = d_s = 0
            for SP in group_conf['species']:
                if SP not in germline_rows or SP not in somatic_rows:
                    continue
                n_g += int(germline_rows[SP][cls]['n_mutations'])
                bp_g += int(germline_rows[SP][cls]['denom_bp'])
                n_s += int(somatic_rows[SP][cls]['n_mutations'])
                d_s += int(somatic_rows[SP][cls]['denom_trinuc'])
            group_cells[group_name] = {
                'g_n': n_g,
                'g_d': 2 * bp_g,
                's_n': n_s,
                's_d': d_s,
            }
        if 'social' not in group_cells or 'subsocial' not in group_cells:
            continue
        soc = group_cells['social']
        sub = group_cells['subsocial']
        fc_soc = ((soc['s_n'] / soc['s_d']) / (soc['g_n'] / soc['g_d'])
                  if soc['g_n'] > 0 and soc['g_d'] > 0 and soc['s_d'] > 0 else float('nan'))
        fc_sub = ((sub['s_n'] / sub['s_d']) / (sub['g_n'] / sub['g_d'])
                  if sub['g_n'] > 0 and sub['g_d'] > 0 and sub['s_d'] > 0 else float('nan'))
        ratio = (fc_soc / fc_sub) if fc_soc and fc_sub and fc_soc > 0 and fc_sub > 0 else float('nan')
        stat, p = poisson_interaction_lrt(
            [soc['g_n'], soc['s_n'], sub['g_n'], sub['s_n']],
            [soc['g_d'], soc['s_d'], sub['g_d'], sub['s_d']],
        )
        fc_test_rows.append({
            'mutation_class': cls,
            'germline_social_n': soc['g_n'],
            'germline_social_denom': soc['g_d'],
            'somatic_social_n': soc['s_n'],
            'somatic_social_denom': soc['s_d'],
            'germline_subsocial_n': sub['g_n'],
            'germline_subsocial_denom': sub['g_d'],
            'somatic_subsocial_n': sub['s_n'],
            'somatic_subsocial_denom': sub['s_d'],
            'fold_change_social': fc_soc,
            'fold_change_subsocial': fc_sub,
            'ratio_of_fold_changes': ratio,
            'lrt_chi2': stat,
            'p_value': p,
        })

    adj = bh_adjust([r['p_value'] if not math.isnan(r['p_value']) else 1.0 for r in fc_test_rows])
    for r, p_adj in zip(fc_test_rows, adj):
        r['p_adjusted'] = p_adj
        r['significant'] = 'yes' if p_adj < 0.05 else 'no'

    fc_test_cols = [
        'mutation_class',
        'germline_social_n', 'germline_social_denom',
        'somatic_social_n', 'somatic_social_denom',
        'germline_subsocial_n', 'germline_subsocial_denom',
        'somatic_subsocial_n', 'somatic_subsocial_denom',
        'fold_change_social', 'fold_change_subsocial', 'ratio_of_fold_changes',
        'lrt_chi2', 'p_value', 'p_adjusted', 'significant',
    ]
    with open(args.output_fold_change_test, 'w') as f:
        f.write('\t'.join(fc_test_cols) + '\n')
        for r in fc_test_rows:
            f.write('\t'.join(f'{r[c]:g}' if isinstance(r[c], float) else str(r[c])
                              for c in fc_test_cols) + '\n')

    # ------------------------------------------------------------------
    # FOLDED SBS96 RATE TSV + FOLD CHANGE TSV
    # ------------------------------------------------------------------
    master96_lines = [
        'group\tspecies\tlevel\tmutation_type\ttrinuc_context\tsbs96_category\t'
        'n_mutations\tdenom\trate\tci_low\tci_high'
    ]

    def emit_master96(group_label, species_label, level, mut, trinuc, sbs96, n, denom_bp, rate_denom):
        rate = (n / rate_denom) if rate_denom > 0 else float('nan')
        lo, hi = poisson_ci(n, rate_denom)
        master96_lines.append(
            f'{group_label}\t{species_label}\t{level}\t{mut}\t{trinuc}\t{sbs96}\t'
            f'{n}\t{denom_bp}\t{rate:g}\t{lo:g}\t{hi:g}'
        )

    for SP in species_list:
        if SP in germline_96 and SP in germline_trinuc:
            counts = germline_96[SP]
            denoms = germline_trinuc[SP]
            for mut, trinuc, sbs96 in sbs96_order:
                bp = denoms.get(trinuc, 0)
                emit_master96(
                    'per_species', SP, 'germline', mut, trinuc, sbs96,
                    counts.get((mut, trinuc, sbs96), 0), bp, 2 * bp
                )
        if SP in somatic_96 and SP in somatic_trinuc:
            counts = somatic_96[SP]
            denoms = somatic_trinuc[SP]
            for mut, trinuc, sbs96 in sbs96_order:
                d = denoms.get(trinuc, 0)
                emit_master96(
                    'per_species', SP, 'somatic', mut, trinuc, sbs96,
                    counts.get((mut, trinuc, sbs96), 0), d, d
                )

    for group_name, group_conf in spectrum_groups.items():
        for mut, trinuc, sbs96 in sbs96_order:
            n_g = bp_g = n_s = d_s = 0
            n_in_g = n_in_s = 0
            for SP in group_conf['species']:
                if SP in germline_96 and SP in germline_trinuc:
                    n_g += germline_96[SP].get((mut, trinuc, sbs96), 0)
                    bp_g += germline_trinuc[SP].get(trinuc, 0)
                    n_in_g += 1
                if SP in somatic_96 and SP in somatic_trinuc:
                    n_s += somatic_96[SP].get((mut, trinuc, sbs96), 0)
                    d_s += somatic_trinuc[SP].get(trinuc, 0)
                    n_in_s += 1
            if n_in_g > 0:
                emit_master96(group_name, 'merged', 'germline', mut, trinuc, sbs96, n_g, bp_g, 2 * bp_g)
            if n_in_s > 0:
                emit_master96(group_name, 'merged', 'somatic', mut, trinuc, sbs96, n_s, d_s, d_s)

    with open(args.output_master_96, 'w') as f:
        f.write('\n'.join(master96_lines) + '\n')

    fc96_lines = [
        'group\tspecies\tmutation_type\ttrinuc_context\tsbs96_category\t'
        'germline_rate\tsomatic_rate\tfold_change\tlog2_fold_change\tfc_ci_low\tfc_ci_high'
    ]

    def emit_fc96(group_label, species_label, mut, trinuc, sbs96, n_g, bp_g, n_s, d_s):
        denom_hap = 2 * bp_g
        if denom_hap <= 0 or d_s <= 0 or n_g == 0 or n_s == 0:
            mu_g = (n_g / denom_hap) if denom_hap > 0 else float('nan')
            mu_s = (n_s / d_s) if d_s > 0 else float('nan')
            fc = (mu_s / mu_g) if (mu_g and mu_g > 0) else float('nan')
            l2 = math.log2(fc) if fc and fc > 0 else float('nan')
            fc96_lines.append(
                f'{group_label}\t{species_label}\t{mut}\t{trinuc}\t{sbs96}\t'
                f'{mu_g:g}\t{mu_s:g}\t{fc:g}\t{l2:g}\t\t'
            )
            return
        mu_g = n_g / denom_hap
        mu_s = n_s / d_s
        fc = mu_s / mu_g
        log_fc = math.log(fc)
        se = math.sqrt(1.0 / n_g + 1.0 / n_s)
        lo = math.exp(log_fc - 1.96 * se)
        hi = math.exp(log_fc + 1.96 * se)
        fc96_lines.append(
            f'{group_label}\t{species_label}\t{mut}\t{trinuc}\t{sbs96}\t'
            f'{mu_g:g}\t{mu_s:g}\t{fc:g}\t{log_fc / math.log(2):g}\t{lo:g}\t{hi:g}'
        )

    for SP in species_list:
        if SP not in germline_96 or SP not in germline_trinuc:
            continue
        if SP not in somatic_96 or SP not in somatic_trinuc:
            continue
        for mut, trinuc, sbs96 in sbs96_order:
            emit_fc96(
                'per_species', SP, mut, trinuc, sbs96,
                germline_96[SP].get((mut, trinuc, sbs96), 0),
                germline_trinuc[SP].get(trinuc, 0),
                somatic_96[SP].get((mut, trinuc, sbs96), 0),
                somatic_trinuc[SP].get(trinuc, 0),
            )

    for group_name, group_conf in spectrum_groups.items():
        for mut, trinuc, sbs96 in sbs96_order:
            n_g = bp_g = n_s = d_s = 0
            n_in = 0
            for SP in group_conf['species']:
                if SP not in germline_96 or SP not in germline_trinuc:
                    continue
                if SP not in somatic_96 or SP not in somatic_trinuc:
                    continue
                n_g += germline_96[SP].get((mut, trinuc, sbs96), 0)
                bp_g += germline_trinuc[SP].get(trinuc, 0)
                n_s += somatic_96[SP].get((mut, trinuc, sbs96), 0)
                d_s += somatic_trinuc[SP].get(trinuc, 0)
                n_in += 1
            if n_in > 0:
                emit_fc96(group_name, 'merged', mut, trinuc, sbs96, n_g, bp_g, n_s, d_s)

    with open(args.output_fold_change_96, 'w') as f:
        f.write('\n'.join(fc96_lines) + '\n')

    print(
        f'Wrote master TSV ({len(master_lines) - 1} data rows), '
        f'per-trio germline ({len(per_trio_lines) - 2} data rows), '
        f'per-species somatic ({len(pss_lines) - 1} data rows), '
        f'fold change ({len(fc_lines) - 1} data rows), '
        f'fold-change test ({len(fc_test_rows)} data rows), '
        f'96-context master ({len(master96_lines) - 1} data rows), '
        f'96-context fold change ({len(fc96_lines) - 1} data rows).'
    )


if __name__ == '__main__':
    main()
