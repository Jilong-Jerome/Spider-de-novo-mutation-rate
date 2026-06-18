#!/usr/bin/env python3
"""
Family-structure sanity check for the sequencing inventory.

Counts how many parents and offspring were successfully sequenced (i.e. have
FASTQ / appear in the sequencing inventory), as UNIQUE individuals, broken down
per family, per species and overall; and, against the expected full family
structure (a complete family = 2 parents + 6 offspring = 8 individuals), how
many individuals that should be there are MISSING, again per family / species /
overall.

Naming conventions (verified across the inventory):
  * offspring : end in 'S<number>' (most species) OR a bare '_<number>'
                (S. bicolor families 3-5, e.g. '150.9_1')
  * mother    : end in a female token  Fema / female / Female / F
  * father    : everything else ending in a male token  Ma / male / Male / M
Mother is tested before father because 'Fema' ends in '...ma'.

Two males (S. mimosarum '7.M', '8.M') are each the sire of two families, so they
fill a father slot in two families while being ONE sequenced sample. Counts are
therefore reported both as unique individuals and, where they differ, as
family-slot records (a parent shared by two families counts once as a unique
individual but fills two family slots).

Input : sequencing_summary_with_coverage.tsv
        columns: species, family, individual, round_id, lane, size_gb, coverage
        (rows with round_id == 'TOTAL' carry per-individual totals and are skipped)
Outputs:
        family_structure_per_family.tsv
        family_structure_per_species.tsv
        family_structure_overall.tsv
"""

import argparse
import os
import re
from collections import defaultdict

import pandas as pd

EXPECTED_PARENTS = 2     # one father + one mother
EXPECTED_OFFSPRING = 6
EXPECTED_TOTAL = EXPECTED_PARENTS + EXPECTED_OFFSPRING   # 8

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(THIS_DIR)
DEFAULT_INPUT = os.path.join(ROOT, 'sequencing_summary_with_coverage.tsv')
OUT_FAMILY = os.path.join(ROOT, 'family_structure_per_family.tsv')
OUT_SPECIES = os.path.join(ROOT, 'family_structure_per_species.tsv')
OUT_OVERALL = os.path.join(ROOT, 'family_structure_overall.tsv')


def classify(individual):
    """Return 'offspring', 'mother' or 'father' for an individual name."""
    if re.search(r'(?:S\d+|_\d+)$', individual):
        return 'offspring'
    if re.search(r'(?:Fema|female|Female|F)$', individual):
        return 'mother'
    return 'father'


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', default=DEFAULT_INPUT)
    parser.add_argument('--out-family', default=OUT_FAMILY)
    parser.add_argument('--out-species', default=OUT_SPECIES)
    parser.add_argument('--out-overall', default=OUT_OVERALL)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t')
    lanes = df[df['round_id'] != 'TOTAL']

    # Unique individuals, each tagged with its (species, family) context(s).
    # An individual appearing under >1 family (a shared sire) yields one record
    # per family context here, so per-family slot counts are correct; we de-dup
    # to unique individuals separately for the species/overall tables.
    seen_slot = set()
    fam = defaultdict(lambda: {'father': [], 'mother': [], 'offspring': []})
    ind_role = {}            # individual -> role (unique individual)
    ind_species = {}         # individual -> species
    for _, r in lanes.iterrows():
        key = (r['species'], r['family'], r['individual'])
        if key in seen_slot:
            continue
        seen_slot.add(key)
        role = classify(r['individual'])
        fam[(r['species'], r['family'])][role].append(r['individual'])
        ind_role[r['individual']] = role
        ind_species[r['individual']] = r['species']

    # ------------------------------------------------------------------
    # Per-family table
    # ------------------------------------------------------------------
    fam_cols = ['species', 'family', 'n_father', 'n_mother', 'n_parents',
                'n_offspring', 'n_sequenced', 'expected_parents',
                'expected_offspring', 'expected_total', 'missing_father',
                'missing_mother', 'missing_offspring', 'n_missing', 'n_surplus',
                'complete']
    fam_rows = []
    for (sp, fm), d in sorted(fam.items()):
        n_father = len(d['father'])
        n_mother = len(d['mother'])
        n_off = len(d['offspring'])
        n_parents = n_father + n_mother
        n_seq = n_parents + n_off
        miss_father = max(0, 1 - n_father)
        miss_mother = max(0, 1 - n_mother)
        miss_off = max(0, EXPECTED_OFFSPRING - n_off)
        n_missing = miss_father + miss_mother + miss_off
        n_surplus = max(0, n_off - EXPECTED_OFFSPRING)
        fam_rows.append({
            'species': sp, 'family': fm,
            'n_father': n_father, 'n_mother': n_mother, 'n_parents': n_parents,
            'n_offspring': n_off, 'n_sequenced': n_seq,
            'expected_parents': EXPECTED_PARENTS,
            'expected_offspring': EXPECTED_OFFSPRING,
            'expected_total': EXPECTED_TOTAL,
            'missing_father': miss_father, 'missing_mother': miss_mother,
            'missing_offspring': miss_off, 'n_missing': n_missing,
            'n_surplus': n_surplus,
            'complete': 'yes' if n_missing == 0 else 'no',
        })
    fam_df = pd.DataFrame(fam_rows, columns=fam_cols)
    fam_df.to_csv(args.out_family, sep='\t', index=False)

    # ------------------------------------------------------------------
    # Per-species table  (parents/offspring as UNIQUE individuals)
    # ------------------------------------------------------------------
    sp_unique = defaultdict(lambda: {'father': set(), 'mother': set(),
                                     'offspring': set(), 'families': set()})
    for ind, role in ind_role.items():
        sp = ind_species[ind]
        sp_unique[sp][role].add(ind)
    for (sp, fm) in fam:
        sp_unique[sp]['families'].add(fm)

    sp_cols = ['species', 'n_families', 'n_fathers_unique', 'n_mothers_unique',
               'n_parents_unique', 'n_offspring_unique', 'n_sequenced_unique',
               'expected_total', 'n_missing', 'n_surplus']
    sp_rows = []
    fam_by_species = defaultdict(list)
    for row in fam_rows:
        fam_by_species[row['species']].append(row)
    for sp in sorted(sp_unique):
        u = sp_unique[sp]
        n_fam = len(u['families'])
        n_fa = len(u['father'])
        n_mo = len(u['mother'])
        n_off = len(u['offspring'])
        n_par = n_fa + n_mo
        miss = sum(r['n_missing'] for r in fam_by_species[sp])
        surp = sum(r['n_surplus'] for r in fam_by_species[sp])
        sp_rows.append({
            'species': sp, 'n_families': n_fam,
            'n_fathers_unique': n_fa, 'n_mothers_unique': n_mo,
            'n_parents_unique': n_par, 'n_offspring_unique': n_off,
            'n_sequenced_unique': n_par + n_off,
            'expected_total': n_fam * EXPECTED_TOTAL,
            'n_missing': miss, 'n_surplus': surp,
        })
    sp_df = pd.DataFrame(sp_rows, columns=sp_cols)
    sp_df.to_csv(args.out_species, sep='\t', index=False)

    # ------------------------------------------------------------------
    # Overall table  (key / value)
    # ------------------------------------------------------------------
    n_fam_total = len(fam)
    n_species = len({sp for sp, _ in fam})
    uniq_fathers = {i for i, r in ind_role.items() if r == 'father'}
    uniq_mothers = {i for i, r in ind_role.items() if r == 'mother'}
    uniq_offspring = {i for i, r in ind_role.items() if r == 'offspring'}
    n_uniq_parents = len(uniq_fathers) + len(uniq_mothers)
    n_uniq_off = len(uniq_offspring)
    n_uniq_total = n_uniq_parents + n_uniq_off
    # family-slot records: an individual counted once per family it belongs to
    n_slot_records = sum(r['n_sequenced'] for r in fam_rows)
    shared = sorted(i for i in ind_role
                    if sum(i in d_role for (s, f), d in fam.items()
                           for d_role in (d['father'], d['mother'],
                                          d['offspring'])) > 1)
    total_missing = sum(r['n_missing'] for r in fam_rows)
    missing_parents = sum(r['missing_father'] + r['missing_mother']
                          for r in fam_rows)
    missing_offspring = sum(r['missing_offspring'] for r in fam_rows)
    total_surplus = sum(r['n_surplus'] for r in fam_rows)

    overall = [
        ('n_species', n_species),
        ('n_families', n_fam_total),
        ('n_sequenced_unique_individuals', n_uniq_total),
        ('n_parents_unique', n_uniq_parents),
        ('n_fathers_unique', len(uniq_fathers)),
        ('n_mothers_unique', len(uniq_mothers)),
        ('n_offspring_unique', n_uniq_off),
        ('n_family_slot_records', n_slot_records),
        ('n_shared_parents', len(shared)),
        ('shared_parents', ';'.join(shared) if shared else ''),
        ('expected_total_slots', n_fam_total * EXPECTED_TOTAL),
        ('n_missing', total_missing),
        ('n_missing_parents', missing_parents),
        ('n_missing_offspring', missing_offspring),
        ('n_surplus_offspring', total_surplus),
    ]
    pd.DataFrame(overall, columns=['metric', 'value']).to_csv(
        args.out_overall, sep='\t', index=False)

    # ------------------------------------------------------------------
    # Console summary
    # ------------------------------------------------------------------
    print(f'Per-family  table -> {args.out_family}   ({len(fam_rows)} families)')
    print(f'Per-species table -> {args.out_species}   ({len(sp_rows)} species)')
    print(f'Overall     table -> {args.out_overall}')
    print()
    print(f'Sequenced unique individuals: {n_uniq_total}  '
          f'({n_uniq_parents} parents = {len(uniq_fathers)} fathers + '
          f'{len(uniq_mothers)} mothers; {n_uniq_off} offspring)')
    print(f'Family-slot records: {n_slot_records}  '
          f'(shared parents filling two families: {", ".join(shared) or "none"})')
    print(f'Expected slots (2 parents + 6 offspring x {n_fam_total} families): '
          f'{n_fam_total * EXPECTED_TOTAL}')
    print(f'Missing vs template: {total_missing} '
          f'({missing_parents} parent(s) + {missing_offspring} offspring); '
          f'surplus offspring: {total_surplus}')
    incomplete = [r for r in fam_rows if r['n_missing'] > 0]
    for r in incomplete:
        bits = []
        if r['missing_father']:
            bits.append('no father')
        if r['missing_mother']:
            bits.append('no mother')
        if r['missing_offspring']:
            bits.append(f"{r['missing_offspring']} offspring short")
        print(f"   {r['species']} {r['family']}: " + ', '.join(bits))


if __name__ == '__main__':
    main()
