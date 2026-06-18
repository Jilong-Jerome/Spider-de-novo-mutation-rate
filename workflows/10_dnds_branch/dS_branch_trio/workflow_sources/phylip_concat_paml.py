#!/usr/bin/env python3
"""
phylip_concat_paml.py - Concatenate multi-block phylip into PAML-compatible format

Reads a phylip file with multiple alignment blocks (output from align_filter_region_fa.py),
concatenates all blocks by matching sequence names, and writes a single PAML-compatible
sequential phylip file.

Usage:
    python3 phylip_concat_paml.py <input.phy> <output.paml.phy>
"""

import sys
from Bio import AlignIO
from collections import OrderedDict


def main():
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Parse all alignment blocks
    alignments = list(AlignIO.parse(open(input_file), "phylip"))

    if not alignments:
        print("ERROR: No alignment blocks found in input file")
        sys.exit(1)

    print(f"Read {len(alignments)} alignment blocks")

    # Concatenate sequences per species
    species_order = [rec.id for rec in alignments[0]]
    concat_seqs = OrderedDict((sp, []) for sp in species_order)

    for aln in alignments:
        for rec in aln:
            concat_seqs[rec.id].append(str(rec.seq))

    # Join all blocks
    for sp in concat_seqs:
        concat_seqs[sp] = ''.join(concat_seqs[sp])

    n_species = len(species_order)
    seq_len = len(concat_seqs[species_order[0]])

    print(f"Concatenated: {n_species} species x {seq_len} bp")

    # Write PAML-compatible sequential phylip
    with open(output_file, 'w') as out:
        out.write(f" {n_species} {seq_len}\n")
        for sp in species_order:
            # PAML phylip: name followed by two spaces then sequence
            out.write(f"{sp}  {concat_seqs[sp]}\n")

    print(f"Written: {output_file}")


if __name__ == '__main__':
    main()
