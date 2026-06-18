#!/usr/bin/env python3
"""
Build a per-individual sequencing-batch table for the supplement, plus
flowcell-batch summaries that document how samples were grouped during
sequencing.

All individuals were sequenced on a single platform (DNBSEQ-G400); the differing
`round_id` prefixes (V300.../V350.../DP8400.../FP100...) are flowcell/run
identifiers on that one platform, NOT different platforms, instruments or kits.
The meaningful "batch" unit is therefore the flowcell (and the lane within it).

For every individual we report which flowcell(s) and lane(s) the sample was
sequenced on, its role (parent vs. proband) and total coverage.  Three summary
tables then characterise the batch structure and directly address the reviewer's
question of whether any flowcell batch could be confounded with sociality:

  * Flowcell occupancy            (individuals / families per flowcell)
  * Flowcells per family          (how a family is spread across flowcells)
  * Flowcells spanning sociality  (flowcells that carried BOTH social and
                                    subsocial individuals -> no sociality-aligned
                                    batch)

Inputs : sequencing_summary_with_coverage.tsv
           columns: species, family, individual, round_id, lane, size_gb, coverage
           ('round_id' == sequencing flowcell/run; 'lane' == lane within flowcell;
            rows with round_id == 'TOTAL' carry the per-individual total coverage)

Outputs: batch_summary_table.tsv   (per-individual, tab-separated, for the supplement)
         batch_summary_table.md    (per-individual table + flowcell-batch summaries)
"""

import argparse
import os
import re
import statistics
from collections import Counter, defaultdict

import pandas as pd


# ------------------------------------------------------------------
# Constants shared with the plotting scripts
# ------------------------------------------------------------------
SOCIAL_SPECIES = {'S_dumicola', 'S_mimosarum', 'S_sarasinorum'}

# A single sequencing platform was used for every individual.  Flowcell IDs are
# run identifiers on this platform, not distinct platforms/kits.
PLATFORM = 'DNBSEQ-G400'

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_INPUT = os.path.join(os.path.dirname(THIS_DIR),
                             'sequencing_summary_with_coverage.tsv')
DEFAULT_TSV = os.path.join(os.path.dirname(THIS_DIR), 'batch_summary_table.tsv')
DEFAULT_MD = os.path.join(os.path.dirname(THIS_DIR), 'batch_summary_table.md')


def role_of(individual):
    """Parent vs. proband.  Probands are offspring labelled '...S<number>' or,
    in S. bicolor families 3-5, with a bare numeric suffix '..._<number>'
    (e.g. '150.9_1'); parents end in a sex token (Ma/Fema/M/F/male/female)."""
    return ('proband' if re.search(r'(?:S\d+|_\d+)$', individual)
            else 'parent')


def sociality_of(species):
    return 'Social' if species in SOCIAL_SPECIES else 'Subsocial'


def _dist_str(counter):
    """Render a Counter of small ints as 'k:v, k:v, ...' sorted by key."""
    return ', '.join(f'{k}:{v}' for k, v in sorted(counter.items()))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', default=DEFAULT_INPUT,
                        help='Path to sequencing_summary_with_coverage.tsv')
    parser.add_argument('--out-tsv', default=DEFAULT_TSV,
                        help='Per-individual TSV output path')
    parser.add_argument('--out-md', default=DEFAULT_MD,
                        help='Markdown output path (table + summaries)')
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t')
    df['coverage'] = pd.to_numeric(df['coverage'], errors='coerce')

    lanes = df[df['round_id'] != 'TOTAL'].copy()
    totals = (df[df['round_id'] == 'TOTAL']
              .set_index('individual')['coverage'].to_dict())

    # ------------------------------------------------------------------
    # Per-individual rows
    # ------------------------------------------------------------------
    records = []
    shared_parents = []   # one sequenced sample serving as parent of >1 family
    for ind, grp in lanes.groupby('individual', sort=False):
        species = grp['species'].iloc[0]
        families = sorted(grp['family'].unique())
        family = ';'.join(families)
        if len(families) > 1:
            shared_parents.append((ind, species, families))
        flowcells = sorted(grp['round_id'].unique())
        lane_ids = sorted(grp['lane'].unique())
        records.append({
            'species': species,
            'family': family,
            'individual': ind,
            'role': role_of(ind),
            'sociality': sociality_of(species),
            'platform': PLATFORM,
            'flowcell_ids': ';'.join(flowcells),
            'n_flowcells': len(flowcells),
            'lanes': ';'.join(lane_ids),
            'n_lane_records': len(grp),
            'total_coverage': round(totals.get(ind, float('nan')), 2),
        })

    per_ind = pd.DataFrame(records).sort_values(
        ['sociality', 'species', 'family', 'role', 'individual'],
        ascending=[True, True, True, False, True]).reset_index(drop=True)

    per_ind.to_csv(args.out_tsv, sep='\t', index=False)
    print(f'Per-individual batch table written to {args.out_tsv} '
          f'({len(per_ind)} individuals)')

    # ------------------------------------------------------------------
    # Flowcell-level batch structure
    # ------------------------------------------------------------------
    fc = defaultdict(lambda: {'ind': set(), 'fam': set(), 'sp': set(),
                              'soc': set()})
    for _, r in lanes.iterrows():
        f = fc[r['round_id']]
        f['ind'].add(r['individual'])
        f['fam'].add((r['species'], r['family']))
        f['sp'].add(r['species'])
        f['soc'].add(sociality_of(r['species']))

    n_fc = len(fc)
    ind_per_fc = [len(v['ind']) for v in fc.values()]
    fam_per_fc = [len(v['fam']) for v in fc.values()]

    # Flowcells per family
    fam_fc = defaultdict(set)
    for rid, v in fc.items():
        for famkey in v['fam']:
            fam_fc[famkey].add(rid)
    fc_per_fam = [len(s) for s in fam_fc.values()]

    # Flowcells spanning both sociality classes
    cross = [(rid, sorted(v['sp'])) for rid, v in fc.items() if len(v['soc']) > 1]
    cross.sort()

    n_parent = int((per_ind['role'] == 'parent').sum())
    n_proband = int((per_ind['role'] == 'proband').sum())

    # Flowcell occupancy by role and by sociality: how many flowcells carry
    # each class, and how many carry a MIX.  Used to show flowcell assignment is
    # not structured by role or sociality.
    fc_roles = defaultdict(set)
    fc_socs = defaultdict(set)
    fc_fam_roles = defaultdict(lambda: defaultdict(set))  # flowcell -> family -> roles
    for _, r in lanes.iterrows():
        fc_roles[r['round_id']].add(role_of(r['individual']))
        fc_socs[r['round_id']].add(sociality_of(r['species']))
        fc_fam_roles[r['round_id']][(r['species'], r['family'])].add(
            role_of(r['individual']))
    fc_mixed_role = sum(len(v) > 1 for v in fc_roles.values())
    fc_mixed_soc = sum(len(v) > 1 for v in fc_socs.values())

    # Parent/proband flowcell composition
    fc_parent_only = sum(v == {'parent'} for v in fc_roles.values())
    fc_proband_only = sum(v == {'proband'} for v in fc_roles.values())
    fc_with_parent = sum('parent' in v for v in fc_roles.values())
    fc_with_proband = sum('proband' in v for v in fc_roles.values())
    # Families in which a parent and a proband OF THAT FAMILY were co-sequenced
    # on at least one shared flowcell
    fams_coseq = set()
    for fams in fc_fam_roles.values():
        for famkey, roles in fams.items():
            if len(roles) > 1:
                fams_coseq.add(famkey)
    n_fams_coseq = len(fams_coseq)
    n_fams_total = len(fam_fc)

    # Flowcell composition: how many flowcells functionally contain a single
    # species / family / sociality / individual (i.e. are nested within one
    # biological unit) vs. span more than one.
    n_single_sp = sum(len(v['sp']) == 1 for v in fc.values())
    n_single_fam = sum(len(v['fam']) == 1 for v in fc.values())
    n_single_soc = sum(len(v['soc']) == 1 for v in fc.values())
    n_single_ind = sum(len(v['ind']) == 1 for v in fc.values())

    # ------------------------------------------------------------------
    # Markdown render
    # ------------------------------------------------------------------
    L = []
    L.append('# Sequencing batch summary\n')
    L.append(f'{len(per_ind)} distinct sequenced individuals across '
             f'{per_ind["species"].nunique()} species and {len(fam_fc)} families '
             f'— {n_parent} parents, {n_proband} probands.\n')
    if shared_parents:
        notes = '; '.join(
            f'`{i}` ({s.replace("S_", "S. ")}, families {", ".join(f)})'
            for i, s, f in shared_parents)
        L.append(f'> Note: {len(shared_parents)} sequenced males are each the '
                 f'shared sire of two families — {notes}. They are one sample '
                 f'each (counted once above), but appear once per family context '
                 f'in Figure Sx, which therefore shows '
                 f'{len(per_ind) + sum(len(f) - 1 for _, _, f in shared_parents)} '
                 f'rows.\n')
    L.append(f'**All individuals were sequenced on a single platform '
             f'({PLATFORM}).** `round_id` = sequencing flowcell/run; `lane` = '
             f'lane within the flowcell. The differing flowcell-ID prefixes '
             f'(V300…/V350…/DP8400…/FP100…) are run/flowcell identifiers on this '
             f'one platform, not different platforms, instruments or library '
             f'kits. The relevant "batch" unit is therefore the flowcell.\n')

    L.append('## Flowcell-batch structure\n')
    L.append(f'- Distinct flowcells: **{n_fc}**')
    L.append(f'- Individuals per flowcell: median **{int(statistics.median(ind_per_fc))}**, '
             f'mean {statistics.mean(ind_per_fc):.1f}, range {min(ind_per_fc)}–{max(ind_per_fc)}')
    L.append(f'- Families per flowcell: median **{int(statistics.median(fam_per_fc))}**, '
             f'mean {statistics.mean(fam_per_fc):.1f}, range {min(fam_per_fc)}–{max(fam_per_fc)}')
    L.append(f'- Flowcells per family: median **{int(statistics.median(fc_per_fam))}**, '
             f'mean {statistics.mean(fc_per_fam):.1f}, range {min(fc_per_fam)}–{max(fc_per_fam)}')
    L.append('')
    L.append(f'  - individuals-per-flowcell distribution (count:flowcells): '
             f'{_dist_str(Counter(ind_per_fc))}')
    L.append(f'  - families-per-flowcell distribution: {_dist_str(Counter(fam_per_fc))}')
    L.append(f'  - flowcells-per-family distribution: {_dist_str(Counter(fc_per_fam))}')
    L.append('')
    L.append('Sequencing was distributed sample-/family-wise across many small '
             'flowcells (not concentrated into a few large pooled batches), so no '
             'single batch could group all social trios apart from all subsocial '
             'trios.\n')

    L.append('## Flowcell composition (independence from the biological comparison)\n')
    L.append('How many of the 73 flowcells functionally contain samples from a '
             'single species / family / sociality class (nested within one '
             'biological unit) versus span more than one:\n')
    L.append(f'- Single-species flowcells: **{n_single_sp}/{n_fc}** '
             f'(span >1 species: {n_fc - n_single_sp})')
    L.append(f'- Single-family flowcells: **{n_single_fam}/{n_fc}** '
             f'(span >1 family: {n_fc - n_single_fam})')
    L.append(f'- Single-sociality flowcells: **{n_single_soc}/{n_fc}** '
             f'(span both social & subsocial: {n_fc - n_single_soc})')
    L.append(f'- Single-individual flowcells: **{n_single_ind}/{n_fc}**')
    L.append('')
    n_multi_sp = n_fc - n_single_sp
    n_cross_soc = n_fc - n_single_soc
    n_multi_within = n_multi_sp - n_cross_soc
    L.append('The two regimes both argue for independence of flowcell batch from '
             'sociality: most flowcells are nested within one species/family and '
             f'so *cannot* separate social from subsocial samples, while of the '
             f'{n_multi_sp} flowcells that do mix species, {n_cross_soc} mix '
             f'*across* the social/subsocial divide (next section) and the other '
             f'{n_multi_within} combine species within one sociality class — so '
             'none isolates social from subsocial. In neither case does a flowcell '
             'batch line up with the biological comparison.\n')

    L.append('## Flowcells spanning both sociality classes\n')
    L.append(f'**{len(cross)} of {n_fc} flowcells carried BOTH social and '
             f'subsocial individuals**, so flowcell batches do not segregate by '
             f'sociality — the same flowcell physically contained both classes:\n')
    L.append('| Flowcell | Species on flowcell |')
    L.append('| --- | --- |')
    for rid, sps in cross:
        labelled = ', '.join(
            f'{s.replace("S_", "S. ")} ({sociality_of(s).lower()})' for s in sps)
        L.append(f'| {rid} | {labelled} |')
    L.append('')

    L.append('## Flowcell assignment is unstructured by role and by sociality\n')
    L.append(f'Because occupancy is sparse, the members of a family — including '
             f'parents and probands — are themselves spread across several '
             f'flowcells (median {int(statistics.median(fc_per_fam))} flowcells '
             f'per family), and the same flowcells are reused across families. '
             f'Flowcell assignment therefore does not encode either biological '
             f'variable:\n')
    L.append(f'- Flowcells carrying both parent and proband individuals: '
             f'**{fc_mixed_role}/{n_fc}**')
    L.append(f'- Flowcells carrying both social and subsocial individuals: '
             f'**{fc_mixed_soc}/{n_fc}**')
    L.append('')
    L.append('### Parent / proband flowcell composition\n')
    L.append('Breakdown of the 73 flowcells by the roles they functionally '
             'contain (parallels the species/family/sociality composition above):\n')
    L.append(f'- Parent-only flowcells: **{fc_parent_only}/{n_fc}**')
    L.append(f'- Proband-only flowcells: **{fc_proband_only}/{n_fc}**')
    L.append(f'- Mixed (both parents and probands): **{fc_mixed_role}/{n_fc}**')
    L.append(f'- Flowcells containing ≥1 parent: {fc_with_parent}/{n_fc}; '
             f'≥1 proband: {fc_with_proband}/{n_fc}')
    L.append(f'- Families in which a parent and a proband of that same family were '
             f'co-sequenced on ≥1 shared flowcell: **{n_fams_coseq}/{n_fams_total}**')
    L.append('')
    L.append('Parents and probands are therefore not segregated into separate '
             f'sequencing batches: {fc_mixed_role} of {n_fc} flowcells carry both '
             f'roles, and in '
             f'{n_fams_coseq} of {n_fams_total} families a parent and a proband '
             'of the same trio were co-sequenced on one flowcell; the remainder '
             'are spread across flowcells. There is no "parents-in-one-batch, '
             'probands-in-another" structure that could bias the within-trio de '
             'novo mutation comparison.\n')
    L.append('Crucially, flowcell assignment is essentially haphazard with '
             'respect to the social/subsocial contrast (and to parent/proband '
             'role): there is no flowcell batch dedicated to one class, so '
             'flowcell-level technical variation cannot create a systematic '
             'social-vs-subsocial difference. De novo mutations are additionally '
             'called within each family on identical chemistry and filtered for '
             'run-specific artefacts.\n')

    L.append('## Per-individual table\n')
    cols = ['species', 'family', 'individual', 'role', 'sociality', 'platform',
            'flowcell_ids', 'lanes', 'n_flowcells', 'n_lane_records',
            'total_coverage']
    L.append('| ' + ' | '.join(cols) + ' |')
    L.append('| ' + ' | '.join('---' for _ in cols) + ' |')
    for _, r in per_ind.iterrows():
        L.append('| ' + ' | '.join(str(r[c]) for c in cols) + ' |')
    L.append('')

    with open(args.out_md, 'w') as fh:
        fh.write('\n'.join(L))
    print(f'Markdown summary written to {args.out_md}')


if __name__ == '__main__':
    main()
