#!/usr/bin/env python3
"""
Estimate parent-specific mutation rates for subsocial species.

The bootstrap samples families and offspring hierarchically. In each replicate,
unphased mutations are assigned by sampling observed phased P/M mutations from
the bootstrapped family. Families with fewer phased mutations naturally
contribute more assignment uncertainty through the resampled phased pool.
"""
import argparse
import csv
import math
import os
import random
from collections import Counter, defaultdict

import yaml


REQUIRED_COLUMNS = {'chrom', 'pos', 'ref', 'alt', 'child', 'Phased'}


def canonical_species(label):
    return label.split(':', 1)[0].strip().upper()


def parse_child(child):
    parts = child.split('_')
    if len(parts) < 4:
        raise ValueError(f'Cannot parse child ID: {child}')
    return parts[0].upper(), parts[1]


def read_leaf_order(newick):
    labels = []
    token = []
    for char in newick.strip().rstrip(';'):
        if char in '(),':
            if token:
                labels.append(canonical_species(''.join(token)))
                token = []
        else:
            token.append(char)
    if token:
        labels.append(canonical_species(''.join(token)))
    return labels


def percentile(values, pct):
    if not values:
        return float('nan')
    xs = sorted(values)
    if len(xs) == 1:
        return xs[0]
    pos = (len(xs) - 1) * pct / 100.0
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return xs[lo]
    frac = pos - lo
    return xs[lo] * (1.0 - frac) + xs[hi] * frac


def mean(values):
    return sum(values) / len(values) if values else float('nan')


def read_config(path):
    with open(path) as f:
        cfg = yaml.safe_load(f)
    species_cfg = cfg.get('species', {})
    exclude = {
        sp: set(values.get('exclude_trios', []) or [])
        for sp, values in species_cfg.items()
    }
    subsocial = set(
        cfg.get('spectrum_groups', {}).get('subsocial', {}).get('species', [])
    )
    social = set(
        cfg.get('spectrum_groups', {}).get('social', {}).get('species', [])
    )
    if not subsocial:
        subsocial = {
            sp for sp, values in species_cfg.items()
            if values.get('lifestyle') == 'subsocial'
        }
    if not social:
        social = {
            sp for sp, values in species_cfg.items()
            if values.get('lifestyle') == 'social'
        }
    leaf_order = read_leaf_order(cfg['species_tree_newick'])
    ordered_subsocial = [sp for sp in leaf_order if sp in subsocial]
    ordered_subsocial.extend(sp for sp in species_cfg if sp in subsocial and sp not in ordered_subsocial)
    ordered_social = [sp for sp in leaf_order if sp in social]
    ordered_social.extend(sp for sp in species_cfg if sp in social and sp not in ordered_social)
    all_species = list(dict.fromkeys(ordered_subsocial + ordered_social))
    return cfg, ordered_subsocial, ordered_social, all_species, exclude


def read_callable_bp(callable_dir, species):
    out = {}
    for sp in species:
        path = os.path.join(callable_dir, f'{sp}_callable_minDP26.tsv')
        per_child = Counter()
        with open(path) as f:
            for line in f:
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 3:
                    continue
                chrom, bp, tag = fields[:3]
                suffix = f'_{chrom}'
                if not tag.endswith(suffix):
                    raise ValueError(f'Cannot parse callable tag in {path}: {tag}')
                child = tag[:-len(suffix)]
                per_child[child] += int(bp)
        out[sp] = dict(per_child)
    return out


def read_mutations(dnm_tsv, species, exclude, callable_bp):
    per_species = {
        sp: defaultdict(lambda: {'offspring': set(), 'loci': {}})
        for sp in species
    }
    for sp in species:
        for child in callable_bp[sp]:
            if child in exclude.get(sp, set()):
                continue
            child_sp, family = parse_child(child)
            if child_sp != sp:
                raise ValueError(f'Callable offspring species mismatch: {child}')
            per_species[sp][family]['offspring'].add(child)

    with open(dnm_tsv, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f'{dnm_tsv} is missing required columns: {", ".join(sorted(missing))}'
            )
        for row in reader:
            child = row['child']
            sp, family = parse_child(child)
            if sp not in per_species:
                continue
            if child in exclude.get(sp, set()):
                continue
            phase = row['Phased']
            if phase == 'F':
                continue
            if phase not in {'P', 'M', 'U'}:
                raise ValueError(f'Unexpected Phased value for {child}: {phase}')
            if child not in callable_bp[sp]:
                raise ValueError(f'Missing callable bp for included offspring: {child}')
            key = (family, row['chrom'], row['pos'], row['ref'], row['alt'])
            locus = per_species[sp][family]['loci'].setdefault(
                key, {'children': set(), 'phases': []}
            )
            locus['children'].add(child)
            locus['phases'].append(phase)

    for sp in species:
        for family, payload in per_species[sp].items():
            resolved = {}
            for key, locus in payload['loci'].items():
                resolved[key] = {
                    'children': locus['children'],
                    'phase': resolve_locus_phase(locus['phases'], sp, key),
                }
            payload['loci'] = resolved

    return per_species


def resolve_locus_phase(phases, sp, key):
    values = set(phases)
    parents = values & {'P', 'M'}
    if len(parents) > 1:
        raise ValueError(
            f'Conflicting paternal and maternal phasing at {sp} locus {key}: '
            f'{",".join(sorted(values))}'
        )
    if parents:
        return next(iter(parents))
    return 'U'


def family_totals(family_payload):
    totals = Counter()
    for locus in family_payload['loci'].values():
        totals[locus['phase']] += 1
    return totals


def phased_pool_from_payload(family_payload):
    pool = []
    for locus in family_payload['loci'].values():
        if locus['phase'] in {'P', 'M'}:
            pool.append(locus['phase'])
    return pool


def phased_pool_from_families(families):
    pool = []
    for family_payload in families.values():
        pool.extend(phased_pool_from_payload(family_payload))
    return pool


def sampled_loci_for_offspring(family_payload, sampled_children):
    selected = {}
    sampled = set(sampled_children)
    for key, locus in family_payload['loci'].items():
        if locus['children'] & sampled:
            selected[key] = locus
    return selected


def bootstrap_species(families, callable_bp, n_bootstrap, seed):
    rng = random.Random(seed)
    family_names = sorted(families)
    observed_p = 0
    observed_m = 0
    unphased = 0
    n_offspring = 0
    observed_callable_bp = 0
    for family in family_names:
        totals = family_totals(families[family])
        observed_p += totals['P']
        observed_m += totals['M']
        unphased += totals['U']
        for child in families[family]['offspring']:
            n_offspring += 1
            observed_callable_bp += callable_bp[child]
    species_phase_pool = phased_pool_from_families(families)
    if not species_phase_pool and unphased:
        raise ValueError(
            'Cannot assign unphased mutations because this species has no '
            'observed phased P/M mutations.'
        )

    paternal_rates = []
    maternal_rates = []
    for _ in range(n_bootstrap):
        boot_p = 0
        boot_m = 0
        boot_callable = 0
        for _family_i in family_names:
            family = rng.choice(family_names)
            offspring = sorted(families[family]['offspring'])
            sampled_children = []
            for _offspring_i in offspring:
                sampled_children.append(rng.choice(offspring))
            unique_children = sorted(set(sampled_children))
            for child in unique_children:
                boot_callable += callable_bp[child]
            sampled_loci = sampled_loci_for_offspring(
                families[family], unique_children
            )
            phase_pool = [
                locus['phase'] for locus in sampled_loci.values()
                if locus['phase'] in {'P', 'M'}
            ] or species_phase_pool
            for locus in sampled_loci.values():
                if locus['phase'] == 'P':
                    boot_p += 1
                elif locus['phase'] == 'M':
                    boot_m += 1
                else:
                    if not phase_pool:
                        raise ValueError(
                            'Cannot assign unphased mutations because no phased '
                            f'P/M mutations were sampled for family {family}.'
                        )
                    if rng.choice(phase_pool) == 'P':
                        boot_p += 1
                    else:
                        boot_m += 1
        denom = 2 * boot_callable
        paternal_rates.append(boot_p / denom if denom else float('nan'))
        maternal_rates.append(boot_m / denom if denom else float('nan'))

    return {
        'n_families': len(family_names),
        'n_offspring': n_offspring,
        'observed_P_mutations': observed_p,
        'observed_M_mutations': observed_m,
        'unphased_mutations': unphased,
        'callable_bp': observed_callable_bp,
        'paternal_rate_observed': observed_parental_rates(families, callable_bp)[0],
        'maternal_rate_observed': observed_parental_rates(families, callable_bp)[1],
        'paternal_rate_mean': mean(paternal_rates),
        'paternal_rate_ci_low': percentile(paternal_rates, 2.5),
        'paternal_rate_ci_high': percentile(paternal_rates, 97.5),
        'maternal_rate_mean': mean(maternal_rates),
        'maternal_rate_ci_low': percentile(maternal_rates, 2.5),
        'maternal_rate_ci_high': percentile(maternal_rates, 97.5),
    }


def observed_parental_rates(families, callable_bp):
    family_names = sorted(families)
    species_p = 0
    species_m = 0
    for family in family_names:
        totals = family_totals(families[family])
        species_p += totals['P']
        species_m += totals['M']
    if species_p + species_m == 0:
        raise ValueError('Cannot assign unphased mutations because species has no observed P/M mutations.')
    species_p_frac = species_p / (species_p + species_m)

    assigned_p_total = 0.0
    assigned_m_total = 0.0
    callable_total = 0
    for family in family_names:
        totals = family_totals(families[family])
        if totals['P'] + totals['M'] > 0:
            p_frac = totals['P'] / (totals['P'] + totals['M'])
        else:
            p_frac = species_p_frac
        assigned_p_total += totals['P'] + totals['U'] * p_frac
        assigned_m_total += totals['M'] + totals['U'] * (1.0 - p_frac)
        for child in families[family]['offspring']:
            callable_total += callable_bp[child]
    denom = 2 * callable_total
    return assigned_p_total / denom, assigned_m_total / denom


def social_germline_summary(families, callable_bp):
    n_mut = 0
    callable_total = 0
    n_offspring = 0
    for family in sorted(families):
        n_mut += len(families[family]['loci'])
        for child in families[family]['offspring']:
            callable_total += callable_bp[child]
            n_offspring += 1
    denom = 2 * callable_total
    return {
        'n_offspring': n_offspring,
        'n_mutations': n_mut,
        'callable_bp': callable_total,
        'social_germline_rate': n_mut / denom if denom else float('nan'),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', required=True)
    ap.add_argument('--dnm_tsv', required=True)
    ap.add_argument('--callable_dir', required=True)
    ap.add_argument('--output_tsv', required=True)
    args = ap.parse_args()

    cfg, subsocial_species, social_species, all_species, exclude = read_config(args.config)
    callable_bp = read_callable_bp(args.callable_dir, all_species)
    mutations = read_mutations(args.dnm_tsv, all_species, exclude, callable_bp)
    n_bootstrap = int(cfg.get('n_bootstrap_parental_rates', 10000))
    seed = int(cfg.get('bootstrap_seed', 42))
    order = {sp: i for i, sp in enumerate(subsocial_species)}

    fieldnames = [
        'record_type',
        'species',
        'phylo_order',
        'n_families',
        'n_offspring',
        'n_mutations',
        'observed_P_mutations',
        'observed_M_mutations',
        'unphased_mutations',
        'callable_bp',
        'paternal_rate_observed',
        'paternal_rate_mean',
        'paternal_rate_ci_low',
        'paternal_rate_ci_high',
        'maternal_rate_observed',
        'maternal_rate_mean',
        'maternal_rate_ci_low',
        'maternal_rate_ci_high',
        'social_germline_rate',
    ]
    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)
    with open(args.output_tsv, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for i, sp in enumerate(subsocial_species):
            if not mutations[sp]:
                continue
            row = bootstrap_species(
                mutations[sp], callable_bp[sp], n_bootstrap, seed + order[sp],
            )
            row['record_type'] = 'subsocial_parental'
            row['species'] = sp
            row['phylo_order'] = i
            row['n_mutations'] = ''
            row['social_germline_rate'] = ''
            writer.writerow(row)
        for i, sp in enumerate(social_species):
            if not mutations[sp]:
                continue
            row = social_germline_summary(mutations[sp], callable_bp[sp])
            row.update({
                'record_type': 'social_reference',
                'species': sp,
                'phylo_order': i,
                'n_families': len(mutations[sp]),
                'observed_P_mutations': '',
                'observed_M_mutations': '',
                'unphased_mutations': '',
                'paternal_rate_observed': '',
                'paternal_rate_mean': '',
                'paternal_rate_ci_low': '',
                'paternal_rate_ci_high': '',
                'maternal_rate_observed': '',
                'maternal_rate_mean': '',
                'maternal_rate_ci_low': '',
                'maternal_rate_ci_high': '',
            })
            writer.writerow(row)


if __name__ == '__main__':
    main()
