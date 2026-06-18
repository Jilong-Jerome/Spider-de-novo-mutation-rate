#!/usr/bin/env python3
"""
prepare_gene_lists.py - Split passing genes into autosomal vs X-linked

Parses GFF3 to build gene_id -> chromosome mapping, then filters
passing_genes.txt for genes on configured autosomal chromosomes.

Usage:
    python3 prepare_gene_lists.py \
        --gff3 BIC.gff3 \
        --passing passing_genes.txt \
        --output_dir steps/07_gene_lists/ \
        --autosomal_chromosomes bic_1,bic_2,bic_3,bic_4,bic_5,bic_6,bic_7,bic_8,bic_9,bic_10,bic_11,bic_12,bic_13,bic_14 \
        --x_chromosome_prefixes bic_X
"""

import argparse
import os
import re


DEFAULT_AUTOSOMAL_CHROMS = ','.join(f'bic_{i}' for i in range(1, 15))
DEFAULT_X_CHROM_PREFIXES = 'bic_X'


def parse_csv_arg(value):
    """Parse a comma-separated command-line list, ignoring empty items."""
    return [item.strip() for item in value.split(',') if item.strip()]


def parse_gff3_genes(gff3_path):
    """Parse GFF3 file and return {gene_id: chromosome} for gene features."""
    gene_chrom = {}
    with open(gff3_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            if fields[2] != 'gene':
                continue
            chrom = fields[0]
            attrs = fields[8]
            # Extract gene ID from attributes
            match = re.search(r'ID=([^;]+)', attrs)
            if match:
                gene_id = match.group(1)
                gene_chrom[gene_id] = chrom
    return gene_chrom


def main():
    parser = argparse.ArgumentParser(description='Split passing genes into autosomal vs X-linked')
    parser.add_argument('--gff3', required=True, help='Path to reference GFF3 annotation')
    parser.add_argument('--passing', required=True, help='Path to passing_genes.txt')
    parser.add_argument('--output_dir', required=True, help='Output directory for gene lists')
    parser.add_argument(
        '--autosomal_chromosomes',
        default=DEFAULT_AUTOSOMAL_CHROMS,
        help='Comma-separated autosomal chromosome IDs to keep'
    )
    parser.add_argument(
        '--x_chromosome_prefixes',
        default=DEFAULT_X_CHROM_PREFIXES,
        help='Comma-separated chromosome prefixes to count as X-linked'
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Parse GFF3
    gene_chrom = parse_gff3_genes(args.gff3)
    print(f"Parsed {len(gene_chrom)} genes from GFF3")

    # Read passing genes
    with open(args.passing) as f:
        passing_genes = [line.strip() for line in f if line.strip()]
    print(f"Read {len(passing_genes)} passing genes")

    auto_chroms = set(parse_csv_arg(args.autosomal_chromosomes))
    x_chrom_prefixes = tuple(parse_csv_arg(args.x_chromosome_prefixes))
    print(f"Autosomal chromosomes: {','.join(sorted(auto_chroms))}")
    print(f"X chromosome prefixes: {','.join(x_chrom_prefixes)}")

    # Write gene_chrom_map.tsv and split into auto/X
    auto_genes = []
    x_genes = []
    map_path = os.path.join(args.output_dir, 'gene_chrom_map.tsv')
    with open(map_path, 'w') as fmap:
        fmap.write('gene_id\tchrom\n')
        for gene_id in passing_genes:
            chrom = gene_chrom.get(gene_id, 'unknown')
            fmap.write(f'{gene_id}\t{chrom}\n')
            if chrom in auto_chroms:
                auto_genes.append(gene_id)
            elif x_chrom_prefixes and chrom.startswith(x_chrom_prefixes):
                x_genes.append(gene_id)
            else:
                print(f"WARNING: gene {gene_id} on unknown chrom {chrom}")

    # Write autosomal gene list
    auto_path = os.path.join(args.output_dir, 'auto_passing_genes.txt')
    with open(auto_path, 'w') as f:
        for gene_id in auto_genes:
            f.write(f'{gene_id}\n')

    print(f"Autosomal genes: {len(auto_genes)}")
    print(f"X-linked genes: {len(x_genes)}")
    print(f"Written: {auto_path}")
    print(f"Written: {map_path}")


if __name__ == '__main__':
    main()
