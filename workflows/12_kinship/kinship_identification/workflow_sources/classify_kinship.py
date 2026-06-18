#!/usr/bin/env python3
"""
classify_kinship.py

Merge KING (kinship coefficient + IBS0) and PLINK (--genome Z0/Z1/Z2/PI_HAT)
pairwise relatedness estimates, classify every individual pair into a
relationship category, and compare the inferred class against the known truth
encoded in the sample names.

Relationship inference
-----------------------
Degree is assigned from the KING robust kinship coefficient:

    kinship > kinship_dup            -> duplicate/MZ
    kinship_first  .. kinship_dup    -> 1st degree
    kinship_second .. kinship_first  -> 2nd degree
    kinship_third  .. kinship_second -> 3rd degree
    kinship <= kinship_third         -> unrelated

1st-degree pairs are split into parent-offspring vs full-sibling using IBS0
(parent-offspring share one allele at every locus, so IBS0 ~ 0; full-sibs have
IBS0 clearly > 0). PLINK Z0 is reported alongside as an independent cross-check
(parent-offspring Z0 ~ 0, Z1 ~ 1; full-sib Z0 ~ 0.25, Z1 ~ 0.5, Z2 ~ 0.25).

Ground truth from sample names  ``{SP}_family{N}_{role}_{...}``
--------------------------------------------------------------
    different family                        -> unrelated
    two parents (F/M) of the same family    -> unrelated (mates / founders)
    parent vs offspring, same family        -> parent-offspring
    offspring vs offspring, same family      -> full-sibling

Usage:
    python classify_kinship.py \
        --king_kin0 {SP}_king.kin0 \
        --king_related {SP}_king_rel.kin0 \
        --plink_genome {SP}_plink.genome \
        --output_tsv {SP}_kinship_classified.tsv \
        --kinship_dup 0.354 --kinship_first 0.177 \
        --kinship_second 0.0884 --kinship_third 0.0442 \
        --ibs0_po_max 0.005
"""
import argparse
import csv
import os


def pair_key(a, b):
    """Order-independent key for an unordered individual pair."""
    return tuple(sorted((a, b)))


def parse_individual(name):
    """Return (family_id, role) where role is one of {'parent', 'offspring'}.

    Sample names look like ``MIM_family1_S1_offspring`` or
    ``MIM_family1_F_female`` / ``MIM_family1_M_male``.
    """
    parts = name.split('_')
    family = parts[1] if len(parts) > 1 else name
    role_token = parts[2] if len(parts) > 2 else ''
    if role_token[:1].upper() == 'S':
        role = 'offspring'
    elif role_token[:1].upper() in ('F', 'M'):
        role = 'parent'
    else:
        role = 'unknown'
    return family, role


def known_relationship(a, b):
    fam_a, role_a = parse_individual(a)
    fam_b, role_b = parse_individual(b)
    if fam_a != fam_b:
        return 'unrelated'
    # same family
    if role_a == 'parent' and role_b == 'parent':
        return 'unrelated'          # the two founders are mates
    if 'unknown' in (role_a, role_b):
        return 'unknown'
    if role_a != role_b:            # one parent + one offspring
        return 'parent-offspring'
    if role_a == 'offspring':       # both offspring
        return 'full-sibling'
    return 'unknown'


def read_king_kin0(path):
    """KING .kin0 -> {pair_key: {'kinship':..,'ibs0':..,'n_snp':..}}."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        header = fh.readline().split()
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.split()
            if not f:
                continue
            id1, id2 = f[idx['ID1']], f[idx['ID2']]
            rec = {
                'kinship': float(f[idx['Kinship']]),
                'ibs0': float(f[idx['IBS0']]) if 'IBS0' in idx else float('nan'),
                'n_snp': f[idx['N_SNP']] if 'N_SNP' in idx else 'NA',
            }
            out[pair_key(id1, id2)] = rec
    return out


def read_king_related(path):
    """KING --related .kin0 -> {pair_key: InfType}."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        header = fh.readline().split()
        idx = {name: i for i, name in enumerate(header)}
        if 'InfType' not in idx:
            return out
        for line in fh:
            f = line.split()
            if not f:
                continue
            out[pair_key(f[idx['ID1']], f[idx['ID2']])] = f[idx['InfType']]
    return out


def read_plink_genome(path):
    """PLINK .genome -> {pair_key: {'z0','z1','z2','pi_hat'}} (whitespace-delimited)."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        header = fh.readline().split()
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.split()
            if not f:
                continue
            id1, id2 = f[idx['IID1']], f[idx['IID2']]
            out[pair_key(id1, id2)] = {
                'z0': float(f[idx['Z0']]),
                'z1': float(f[idx['Z1']]),
                'z2': float(f[idx['Z2']]),
                'pi_hat': float(f[idx['PI_HAT']]),
            }
    return out


def classify(kinship, ibs0, args):
    """Assign a relationship category from kinship coefficient and IBS0."""
    if kinship > args.kinship_dup:
        return 'duplicate/MZ'
    if kinship > args.kinship_first:
        # 1st degree: split parent-offspring vs full-sibling using IBS0
        if ibs0 == ibs0 and ibs0 <= args.ibs0_po_max:   # ibs0==ibs0 filters NaN
            return 'parent-offspring'
        return 'full-sibling'
    if kinship > args.kinship_second:
        return '2nd-degree'
    if kinship > args.kinship_third:
        return '3rd-degree'
    return 'unrelated'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--king_kin0', required=True)
    ap.add_argument('--king_related', default=None)
    ap.add_argument('--plink_genome', required=True)
    ap.add_argument('--output_tsv', required=True)
    ap.add_argument('--kinship_dup', type=float, default=0.354)
    ap.add_argument('--kinship_first', type=float, default=0.177)
    ap.add_argument('--kinship_second', type=float, default=0.0884)
    ap.add_argument('--kinship_third', type=float, default=0.0442)
    ap.add_argument('--ibs0_po_max', type=float, default=0.005)
    args = ap.parse_args()

    king = read_king_kin0(args.king_kin0)
    related = read_king_related(args.king_related)
    plink = read_plink_genome(args.plink_genome)

    all_pairs = set(king) | set(plink)

    rows = []
    for key in sorted(all_pairs):
        id1, id2 = key
        k = king.get(key, {})
        p = plink.get(key, {})
        kinship = k.get('kinship', float('nan'))
        ibs0 = k.get('ibs0', float('nan'))
        inferred = classify(kinship, ibs0, args) if kinship == kinship else 'NA'
        truth = known_relationship(id1, id2)
        rows.append({
            'ID1': id1,
            'ID2': id2,
            'n_snp': k.get('n_snp', 'NA'),
            'kinship': kinship,
            'IBS0': ibs0,
            'Z0': p.get('z0', float('nan')),
            'Z1': p.get('z1', float('nan')),
            'Z2': p.get('z2', float('nan')),
            'PI_HAT': p.get('pi_hat', float('nan')),
            'inferred_class': inferred,
            'king_inftype': related.get(key, 'NA'),
            'known_class': truth,
            'agree': (inferred == truth),
        })

    fieldnames = ['ID1', 'ID2', 'n_snp', 'kinship', 'IBS0', 'Z0', 'Z1', 'Z2',
                  'PI_HAT', 'inferred_class', 'king_inftype', 'known_class', 'agree']
    with open(args.output_tsv, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    # Console summary
    n = len(rows)
    n_known = sum(1 for r in rows if r['known_class'] != 'unknown')
    n_agree = sum(1 for r in rows if r['agree'] and r['known_class'] != 'unknown')
    print("Wrote %d pairs to %s" % (n, args.output_tsv))
    if n_known:
        print("Agreement vs known truth: %d/%d (%.1f%%)"
              % (n_agree, n_known, 100.0 * n_agree / n_known))
        mismatches = [r for r in rows
                      if not r['agree'] and r['known_class'] != 'unknown']
        for r in mismatches:
            print("  MISMATCH: %s  %s  inferred=%s  known=%s  "
                  "kinship=%.4f IBS0=%.4f"
                  % (r['ID1'], r['ID2'], r['inferred_class'], r['known_class'],
                     r['kinship'], r['IBS0']))


if __name__ == '__main__':
    main()
