#!/usr/bin/env python3
"""
calculate_heterozygosity.py - Calculate heterozygosity metrics for an individual

This script reads genotype data and calculates:
1. Heterozygosity for whole gene, exon-only, and intron-only regions
2. List of all heterozygous sites
"""

import argparse
import os
import sys


def calculate_heterozygosity(genotypes_file, individual, gene, output_het, output_sites):
    """
    Calculate heterozygosity metrics from genotypes file.

    Args:
        genotypes_file: Path to genotypes TSV from extract_genotypes.py
        individual: Individual identifier
        gene: Gene name
        output_het: Output file for heterozygosity metrics
        output_sites: Output file for heterozygous site list
    """
    # Initialize counters
    counts = {
        'whole_gene': {'callable': 0, 'het': 0},
        'exon_only': {'callable': 0, 'het': 0},
        'intron_only': {'callable': 0, 'het': 0}
    }

    het_sites = []

    # Read genotypes file
    with open(genotypes_file, 'r') as f:
        header = f.readline().strip().split('\t')

        # Get column indices
        col_idx = {name: i for i, name in enumerate(header)}

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < len(header):
                continue

            chrom = parts[col_idx['chrom']]
            pos = parts[col_idx['pos']]
            ref = parts[col_idx['ref']]
            alt = parts[col_idx['alt']]
            genotype = parts[col_idx['genotype']]
            is_callable = parts[col_idx['is_callable']] == 'True'
            is_het = parts[col_idx['is_het']] == 'True'
            region_type = parts[col_idx['region_type']]

            if not is_callable:
                continue

            # Count for whole gene
            counts['whole_gene']['callable'] += 1
            if is_het:
                counts['whole_gene']['het'] += 1
                het_sites.append({
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'genotype': genotype,
                    'region_type': region_type
                })

            # Count for specific regions
            if region_type == 'exon':
                counts['exon_only']['callable'] += 1
                if is_het:
                    counts['exon_only']['het'] += 1
            elif region_type == 'intron':
                counts['intron_only']['callable'] += 1
                if is_het:
                    counts['intron_only']['het'] += 1

    # Calculate heterozygosity rates
    results = []
    for region in ['whole_gene', 'exon_only', 'intron_only']:
        callable_sites = counts[region]['callable']
        het_count = counts[region]['het']

        if callable_sites > 0:
            het_rate = het_count / callable_sites
        else:
            het_rate = 'NA'

        results.append({
            'individual': individual,
            'gene': gene,
            'region': region,
            'het_sites': het_count,
            'callable_sites': callable_sites,
            'heterozygosity': het_rate
        })

    # Write heterozygosity output
    with open(output_het, 'w') as f:
        f.write("individual\tgene\tregion\thet_sites\tcallable_sites\theterozygosity\n")
        for r in results:
            het_str = f"{r['heterozygosity']:.6f}" if isinstance(r['heterozygosity'], float) else r['heterozygosity']
            f.write(f"{r['individual']}\t{r['gene']}\t{r['region']}\t"
                    f"{r['het_sites']}\t{r['callable_sites']}\t{het_str}\n")

    # Write het sites list
    with open(output_sites, 'w') as f:
        f.write("individual\tgene\tchrom\tpos\tref\talt\tgenotype\tregion_type\n")
        for site in het_sites:
            f.write(f"{individual}\t{gene}\t{site['chrom']}\t{site['pos']}\t"
                    f"{site['ref']}\t{site['alt']}\t{site['genotype']}\t{site['region_type']}\n")

    # Print summary
    print(f"Heterozygosity Summary for {individual} - {gene}:")
    print("-" * 60)
    for r in results:
        het_str = f"{r['heterozygosity']:.6f}" if isinstance(r['heterozygosity'], float) else r['heterozygosity']
        print(f"  {r['region']}: {r['het_sites']}/{r['callable_sites']} = {het_str}")
    print(f"\nTotal heterozygous sites: {len(het_sites)}")
    print(f"Output files:")
    print(f"  Heterozygosity: {output_het}")
    print(f"  Het sites: {output_sites}")


def main():
    parser = argparse.ArgumentParser(
        description='Calculate heterozygosity metrics for an individual'
    )
    parser.add_argument('--genotypes', required=True, help='Genotypes TSV file')
    parser.add_argument('--individual', required=True, help='Individual identifier')
    parser.add_argument('--gene', required=True, help='Gene name')
    parser.add_argument('--output_het', required=True, help='Output heterozygosity TSV')
    parser.add_argument('--output_sites', required=True, help='Output het sites TSV')

    args = parser.parse_args()

    # Create output directories if needed
    for output in [args.output_het, args.output_sites]:
        output_dir = os.path.dirname(output)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

    calculate_heterozygosity(
        genotypes_file=args.genotypes,
        individual=args.individual,
        gene=args.gene,
        output_het=args.output_het,
        output_sites=args.output_sites
    )


if __name__ == "__main__":
    main()
