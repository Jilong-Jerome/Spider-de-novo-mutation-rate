#!/usr/bin/env python3
"""
Extract human testis-specific gene protein sequences from UniProt proteome.

Filters the proteome FASTA by matching gene names (GN= field) against
the list of testis-specific genes.
"""

import argparse
import re
import sys
from pathlib import Path


def parse_gene_name(header: str) -> str | None:
    """Extract gene name from UniProt FASTA header.

    Example header:
    >sp|A0A087X1C5|CP2D7_HUMAN Cytochrome P450 2D7 OS=Homo sapiens OX=9606 GN=CYP2D7 PE=1 SV=1
    """
    match = re.search(r'\sGN=(\S+)', header)
    if match:
        return match.group(1)
    return None


def load_gene_list(gene_file: Path) -> set:
    """Load testis-specific gene names from file."""
    genes = set()
    with open(gene_file) as f:
        for line in f:
            gene = line.strip()
            if gene:
                genes.add(gene)
    return genes


def extract_proteins(proteome_file: Path, gene_list: set, output_file: Path) -> dict:
    """Extract proteins matching gene list from proteome.

    Returns statistics about extraction.
    """
    stats = {
        'total_proteins': 0,
        'matched_proteins': 0,
        'genes_found': set(),
        'genes_not_found': set()
    }

    current_header = None
    current_seq = []
    matched = False

    with open(proteome_file) as fin, open(output_file, 'w') as fout:
        for line in fin:
            line = line.rstrip()

            if line.startswith('>'):
                # Write previous protein if matched
                if matched and current_header:
                    fout.write(f"{current_header}\n")
                    fout.write('\n'.join(current_seq) + '\n')
                    stats['matched_proteins'] += 1

                # Start new protein
                stats['total_proteins'] += 1
                current_header = line
                current_seq = []

                # Check if gene matches
                gene_name = parse_gene_name(line)
                if gene_name and gene_name in gene_list:
                    matched = True
                    stats['genes_found'].add(gene_name)
                else:
                    matched = False
            else:
                current_seq.append(line)

        # Don't forget last protein
        if matched and current_header:
            fout.write(f"{current_header}\n")
            fout.write('\n'.join(current_seq) + '\n')
            stats['matched_proteins'] += 1

    stats['genes_not_found'] = gene_list - stats['genes_found']
    return stats


def main():
    parser = argparse.ArgumentParser(
        description='Extract testis-specific gene proteins from UniProt proteome'
    )
    parser.add_argument('--proteome', required=True, type=Path,
                        help='Input UniProt proteome FASTA file')
    parser.add_argument('--genes', required=True, type=Path,
                        help='File with testis-specific gene names (one per line)')
    parser.add_argument('--output', required=True, type=Path,
                        help='Output FASTA file with matched proteins')
    parser.add_argument('--stats', type=Path,
                        help='Optional: output statistics file')

    args = parser.parse_args()

    # Load gene list
    print(f"Loading gene list from {args.genes}...")
    gene_list = load_gene_list(args.genes)
    print(f"  Loaded {len(gene_list)} gene names")

    # Extract proteins
    print(f"Extracting proteins from {args.proteome}...")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    stats = extract_proteins(args.proteome, gene_list, args.output)

    # Print summary
    print(f"\nResults:")
    print(f"  Total proteins in proteome: {stats['total_proteins']}")
    print(f"  Proteins matched: {stats['matched_proteins']}")
    print(f"  Genes found: {len(stats['genes_found'])} / {len(gene_list)}")
    print(f"  Genes not found: {len(stats['genes_not_found'])}")

    if stats['genes_not_found']:
        print(f"\nFirst 10 genes not found:")
        for gene in sorted(stats['genes_not_found'])[:10]:
            print(f"    {gene}")

    # Write stats file if requested
    if args.stats:
        with open(args.stats, 'w') as f:
            f.write("Testis Protein Extraction Summary\n")
            f.write("=" * 40 + "\n\n")
            f.write(f"Total proteins in proteome: {stats['total_proteins']}\n")
            f.write(f"Proteins matched: {stats['matched_proteins']}\n")
            f.write(f"Genes in input list: {len(gene_list)}\n")
            f.write(f"Genes found: {len(stats['genes_found'])}\n")
            f.write(f"Genes not found: {len(stats['genes_not_found'])}\n\n")

            f.write("Genes found:\n")
            for gene in sorted(stats['genes_found']):
                f.write(f"  {gene}\n")

            f.write("\nGenes not found:\n")
            for gene in sorted(stats['genes_not_found']):
                f.write(f"  {gene}\n")

    print(f"\nOutput written to {args.output}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
