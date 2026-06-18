#!/usr/bin/env python3
"""
concat_genes.py - Concatenate per-gene FASTA files into a supermatrix alignment

Reads a gene ID list and concatenates CDS sequences per species across all genes.
Each gene FASTA has 3 sequences (BIC, SAR, PAC) of the same length.
Output is a single FASTA with 3 sequences where each is the concatenation of all genes.

Usage:
    python3 concat_genes.py <id_list_file> <fasta_dir> <output_fasta>
"""

import sys
import os


def main():
    id_list_file = sys.argv[1]
    fasta_dir = sys.argv[2]
    output_fasta = sys.argv[3]

    # Read gene IDs
    with open(id_list_file) as f:
        gene_ids = [line.strip() for line in f if line.strip()]

    print(f"Concatenating {len(gene_ids)} genes...")

    # Accumulate sequences per species
    species_order = []
    species_seqs = {}
    n_missing = 0

    for gene_id in gene_ids:
        fa_path = os.path.join(fasta_dir, f"{gene_id}.fa")
        if not os.path.isfile(fa_path):
            print(f"WARNING: missing {fa_path}")
            n_missing += 1
            continue

        # Parse simple FASTA (3 sequences, single-line each)
        current_sp = None
        with open(fa_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    current_sp = line[1:]
                    if current_sp not in species_seqs:
                        species_order.append(current_sp)
                        species_seqs[current_sp] = []
                else:
                    species_seqs[current_sp].append(line)

    if n_missing > 0:
        print(f"WARNING: {n_missing} genes missing")

    # Write concatenated FASTA
    with open(output_fasta, 'w') as out:
        for sp in species_order:
            out.write(f">{sp}\n")
            out.write(''.join(species_seqs[sp]))
            out.write('\n')

    # Report total length
    total_len = len(''.join(species_seqs[species_order[0]])) if species_order else 0
    print(f"Concatenated alignment: {len(species_order)} species x {total_len} bp")


if __name__ == '__main__':
    main()
