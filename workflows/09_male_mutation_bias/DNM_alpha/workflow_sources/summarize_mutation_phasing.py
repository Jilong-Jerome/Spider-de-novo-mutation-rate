#!/usr/bin/env python3
"""
Summarize parental phasing for unique and sibling-shared DNM loci.

Rows with Phased == F and configured excluded trios are removed before loci are
classified. Shared loci with one parent label plus U are assigned to the parent
label; loci containing both P and M fail explicitly.
"""
import argparse
import csv
import os
from collections import Counter, defaultdict

import yaml


REQUIRED_COLUMNS = {'chrom', 'pos', 'ref', 'alt', 'child', 'Phased'}
PHASES = ('P', 'M', 'U')
SHARING_CLASSES = ('unique', 'sibling_shared')


def canonical_species(label):
    return label.split(':', 1)[0].strip().upper()


def parse_child(child):
    parts = child.split('_')
    if len(parts) < 4:
        raise ValueError(f'Cannot parse child ID: {child}')
    return parts[0].upper(), parts[1]


def read_config(path):
    with open(path) as f:
        cfg = yaml.safe_load(f)
    species_cfg = cfg.get('species', {})
    species = list(species_cfg.keys())
    exclude = {
        sp: set(values.get('exclude_trios', []) or [])
        for sp, values in species_cfg.items()
    }
    return cfg, species, exclude


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


def assign_locus_phase(phases, sp, key):
    values = set(phases)
    parents = values & {'P', 'M'}
    if len(parents) > 1:
        raise ValueError(
            f'Conflicting paternal and maternal phasing at {sp} locus {key}: '
            f'{",".join(sorted(values))}'
        )
    if parents:
        return next(iter(parents))
    if values == {'U'}:
        return 'U'
    unexpected = values - {'P', 'M', 'U'}
    if unexpected:
        raise ValueError(
            f'Unexpected phasing labels at {sp} locus {key}: '
            f'{",".join(sorted(values))}'
        )
    return 'U'


def summarize(dnm_tsv, species, exclude):
    loci = defaultdict(lambda: defaultdict(lambda: {'children': set(), 'phases': []}))

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
            if sp not in species:
                continue
            if child in exclude.get(sp, set()):
                continue
            phase = row['Phased']
            if phase == 'F':
                continue
            if phase not in PHASES:
                raise ValueError(f'Unexpected Phased value for {child}: {phase}')

            key = (family, row['chrom'], row['pos'], row['ref'], row['alt'])
            loci[sp][key]['children'].add(child)
            loci[sp][key]['phases'].append(phase)

    rows = []
    for sp in species:
        counts = {sharing_class: Counter() for sharing_class in SHARING_CLASSES}
        for key, payload in loci.get(sp, {}).items():
            sharing_class = (
                'sibling_shared' if len(payload['children']) >= 2 else 'unique'
            )
            phase = assign_locus_phase(payload['phases'], sp, key)
            counts[sharing_class][phase] += 1

        for sharing_class in SHARING_CLASSES:
            phase_counts = counts[sharing_class]
            rows.append({
                'species': sp,
                'sharing_class': sharing_class,
                'n_P_loci': phase_counts['P'],
                'n_M_loci': phase_counts['M'],
                'n_U_loci': phase_counts['U'],
                'n_total_loci': sum(phase_counts[p] for p in PHASES),
            })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', required=True)
    ap.add_argument('--dnm_tsv', required=True)
    ap.add_argument('--output_tsv', required=True)
    args = ap.parse_args()

    cfg, species, exclude = read_config(args.config)
    leaf_order = read_leaf_order(cfg['species_tree_newick'])
    if leaf_order:
        ordered_species = [sp for sp in leaf_order if sp in species]
        ordered_species.extend(sp for sp in species if sp not in ordered_species)
    else:
        ordered_species = species

    rows = summarize(args.dnm_tsv, ordered_species, exclude)

    fieldnames = [
        'species',
        'phylo_order',
        'sharing_class',
        'n_P_loci',
        'n_M_loci',
        'n_U_loci',
        'n_total_loci',
    ]
    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)
    with open(args.output_tsv, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        order = {sp: i for i, sp in enumerate(ordered_species)}
        for row in rows:
            row = dict(row)
            row['phylo_order'] = order[row['species']]
            writer.writerow(row)


if __name__ == '__main__':
    main()
