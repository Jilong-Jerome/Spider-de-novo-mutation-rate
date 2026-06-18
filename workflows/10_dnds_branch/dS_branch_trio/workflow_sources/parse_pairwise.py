#!/usr/bin/env python3
"""
parse_pairwise.py - Parse pairwise codeml output (runmode = -2)

Reads 2ML.dN and 2ML.dS lower-triangular distance matrices and outputs
a tab-separated table of pairwise dN, dS, dN/dS for all species pairs.

Usage:
    python3 parse_pairwise.py <work_dir> <output_tab> <replicate_id>
"""

import sys
import os


def parse_distance_matrix(filepath):
    """Parse a PAML lower-triangular distance matrix (2ML.dN or 2ML.dS).

    Format:
         3
        BIC
        SAR            0.0123
        PAC            0.0456  0.0789

    Returns:
        species: list of species names
        dist: dict of {(sp1, sp2): distance} for all pairs (sp1 < sp2 alphabetically)
    """
    with open(filepath) as f:
        lines = [l.rstrip() for l in f if l.strip()]

    n = int(lines[0].strip())
    species = []
    dist = {}

    for i in range(1, n + 1):
        parts = lines[i].split()
        sp = parts[0]
        species.append(sp)
        values = [float(x) for x in parts[1:]]
        for j, val in enumerate(values):
            pair = tuple(sorted([species[j], sp]))
            dist[pair] = val

    return species, dist


def main():
    work_dir = sys.argv[1]
    output_tab = sys.argv[2]
    rep_id = sys.argv[3]

    dn_file = os.path.join(work_dir, '2ML.dN')
    ds_file = os.path.join(work_dir, '2ML.dS')

    if not os.path.isfile(dn_file) or not os.path.isfile(ds_file):
        print(f"ERROR: missing {dn_file} or {ds_file}")
        sys.exit(1)

    species, dn_dist = parse_distance_matrix(dn_file)
    _, ds_dist = parse_distance_matrix(ds_file)

    with open(output_tab, 'w') as out:
        out.write("id\tsp1\tsp2\tdN\tdS\tdNdS\n")
        for pair in sorted(dn_dist.keys()):
            sp1, sp2 = pair
            dN = dn_dist[pair]
            dS = ds_dist.get(pair, 0)
            dNdS = dN / dS if dS > 0 else 'NA'
            if dNdS != 'NA':
                out.write(f"{rep_id}\t{sp1}\t{sp2}\t{dN:.6f}\t{dS:.6f}\t{dNdS:.6f}\n")
            else:
                out.write(f"{rep_id}\t{sp1}\t{sp2}\t{dN:.6f}\t{dS:.6f}\tNA\n")

    print(f"Pairwise results written to {output_tab}")


if __name__ == '__main__':
    main()
