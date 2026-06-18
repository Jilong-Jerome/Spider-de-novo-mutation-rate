#!/usr/bin/env python3
"""
Summarize unique and sibling-shared germline DNM loci per species.

Sibling sharing is evaluated within species and family after excluding trios
listed in the workflow configuration. A locus is defined by family, chrom,
pos, ref, and alt; child rows carrying the same locus are collapsed before
classification.
"""
import argparse
import csv
import os
from collections import defaultdict

import yaml


REQUIRED_COLUMNS = {'chrom', 'pos', 'ref', 'alt', 'child', 'Phased'}


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


def summarize(dnm_tsv, species, exclude):
    loci = defaultdict(lambda: defaultdict(set))
    included_rows = defaultdict(int)
    excluded_rows = defaultdict(int)
    families = defaultdict(set)
    offspring = defaultdict(set)

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
                excluded_rows[sp] += 1
                continue
            if row['Phased'] == 'F':
                excluded_rows[sp] += 1
                continue

            key = (family, row['chrom'], row['pos'], row['ref'], row['alt'])
            loci[sp][key].add(child)
            included_rows[sp] += 1
            families[sp].add(family)
            offspring[sp].add(child)

    rows = []
    for sp in species:
        unique = 0
        shared = 0
        for children in loci.get(sp, {}).values():
            if len(children) >= 2:
                shared += 1
            else:
                unique += 1
        rows.append({
            'species': sp,
            'n_unique_mutation_loci': unique,
            'n_sibling_shared_mutation_loci': shared,
            'n_total_mutation_loci': unique + shared,
            'n_child_dnm_rows_included': included_rows[sp],
            'n_child_dnm_rows_excluded': excluded_rows[sp],
            'n_families_included': len(families[sp]),
            'n_offspring_included': len(offspring[sp]),
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
        'n_unique_mutation_loci',
        'n_sibling_shared_mutation_loci',
        'n_total_mutation_loci',
        'n_child_dnm_rows_included',
        'n_child_dnm_rows_excluded',
        'n_families_included',
        'n_offspring_included',
    ]
    os.makedirs(os.path.dirname(args.output_tsv), exist_ok=True)
    with open(args.output_tsv, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for i, row in enumerate(rows):
            row = dict(row)
            row['phylo_order'] = i
            writer.writerow(row)


if __name__ == '__main__':
    main()
