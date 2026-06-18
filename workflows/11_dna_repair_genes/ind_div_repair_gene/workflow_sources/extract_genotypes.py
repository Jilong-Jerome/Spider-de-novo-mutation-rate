#!/usr/bin/env python3
"""
extract_genotypes.py - Extract genotypes from gVCF for a gene region

This script extracts diploid genotypes from a gVCF file for positions within
a specified gene region, classifying each position as exon, intron, or other.
"""

import argparse
import os
import sys
from collections import defaultdict


def load_bed_intervals(bed_file):
    """Load BED file and return list of (chrom, start, end) tuples"""
    intervals = []
    if not os.path.exists(bed_file):
        return intervals

    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                intervals.append((chrom, start, end))
    return intervals


def build_interval_lookup(intervals):
    """
    Build a dictionary for fast interval lookup.
    Returns dict: chrom -> sorted list of (start, end) tuples
    """
    lookup = defaultdict(list)
    for chrom, start, end in intervals:
        lookup[chrom].append((start, end))

    # Sort intervals by start position
    for chrom in lookup:
        lookup[chrom].sort()

    return lookup


def position_in_intervals(chrom, pos, interval_lookup):
    """Check if a position falls within any interval"""
    if chrom not in interval_lookup:
        return False

    intervals = interval_lookup[chrom]
    # Binary search could be used for optimization, but linear is fine for gene-sized regions
    for start, end in intervals:
        if start <= pos <= end:
            return True
        if start > pos:
            break
    return False


def classify_region(chrom, pos, exon_lookup, intron_lookup):
    """Classify a position as exon, intron, or other"""
    if position_in_intervals(chrom, pos, exon_lookup):
        return "exon"
    elif position_in_intervals(chrom, pos, intron_lookup):
        return "intron"
    else:
        return "other"


def parse_vcf_genotype(gt_string):
    """Parse genotype string like '0/0' or '0/1' or '1/1'"""
    if gt_string in ['.', './.', '.|.']:
        return None
    # Handle both / and | separators
    gt_string = gt_string.replace('|', '/')
    return gt_string


def get_alleles_from_genotype(ref, alt, genotype):
    """
    Get actual alleles from genotype.
    Returns tuple of (allele1, allele2) or None if uncallable.
    """
    if genotype is None:
        return None

    # Handle <NON_REF> and missing data
    alt_alleles = alt.split(',') if alt else []
    # Filter out <NON_REF> and empty
    alt_alleles = [a for a in alt_alleles if a and a != '<NON_REF>' and a != '.']

    allele_list = [ref] + alt_alleles

    try:
        parts = genotype.split('/')
        if len(parts) != 2:
            return None

        idx1 = int(parts[0])
        idx2 = int(parts[1])

        # Check if indices are valid
        if idx1 >= len(allele_list) or idx2 >= len(allele_list):
            # If genotype refers to <NON_REF>, treat as ref
            if idx1 == 0 and idx2 == 0:
                return (ref, ref)
            elif idx1 == 0:
                return (ref, ref)  # Treat as hom ref if alt is <NON_REF>
            elif idx2 == 0:
                return (ref, ref)
            return None

        return (allele_list[idx1], allele_list[idx2])

    except (ValueError, IndexError):
        return None


def extract_genotypes(gvcf_file, gene_bed, exon_bed, intron_bed, output_file,
                      min_dp=10, min_gq=20):
    """
    Extract genotypes from gVCF for positions within gene region.
    """
    # Load interval data
    gene_intervals = load_bed_intervals(gene_bed)
    exon_intervals = load_bed_intervals(exon_bed)
    intron_intervals = load_bed_intervals(intron_bed)

    if not gene_intervals:
        raise ValueError(f"No gene intervals found in {gene_bed}")

    # Get gene boundaries
    gene_chrom = gene_intervals[0][0]
    gene_start = min(s for c, s, e in gene_intervals)
    gene_end = max(e for c, s, e in gene_intervals)

    # Build lookup structures
    exon_lookup = build_interval_lookup(exon_intervals)
    intron_lookup = build_interval_lookup(intron_intervals)

    # Process gVCF
    results = []
    in_gene_region = False

    with open(gvcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue

            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            info = parts[7]
            format_fields = parts[8].split(':')
            sample_data = parts[9].split(':')

            # Check if we're in the gene region
            if chrom != gene_chrom:
                if in_gene_region:
                    break  # Moved past gene region
                continue

            if pos < gene_start:
                continue
            if pos > gene_end:
                break  # Past gene region

            in_gene_region = True

            # Parse format fields
            format_dict = dict(zip(format_fields, sample_data))

            # Extract genotype
            gt = format_dict.get('GT', '.')
            genotype = parse_vcf_genotype(gt)

            # Extract depth and genotype quality
            try:
                dp = int(format_dict.get('DP', '0'))
            except ValueError:
                dp = 0

            try:
                gq = int(format_dict.get('GQ', '0'))
            except ValueError:
                gq = 0

            # Determine callability
            is_callable = (dp >= min_dp and gq >= min_gq and genotype is not None)

            # Classify region
            region_type = classify_region(chrom, pos, exon_lookup, intron_lookup)

            # Get alleles
            alleles = get_alleles_from_genotype(ref, alt, genotype)

            # Determine if heterozygous
            is_het = False
            if alleles and alleles[0] != alleles[1]:
                is_het = True

            # Format output
            alt_out = alt if alt and alt != '<NON_REF>' else '.'
            gt_out = genotype if genotype else '.'

            results.append({
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt_out,
                'genotype': gt_out,
                'depth': dp,
                'gq': gq,
                'is_callable': is_callable,
                'is_het': is_het,
                'region_type': region_type
            })

    # Write output
    with open(output_file, 'w') as f:
        # Header
        f.write("chrom\tpos\tref\talt\tgenotype\tdepth\tgq\tis_callable\tis_het\tregion_type\n")

        for r in results:
            f.write(f"{r['chrom']}\t{r['pos']}\t{r['ref']}\t{r['alt']}\t"
                    f"{r['genotype']}\t{r['depth']}\t{r['gq']}\t"
                    f"{r['is_callable']}\t{r['is_het']}\t{r['region_type']}\n")

    # Print summary
    total = len(results)
    callable_count = sum(1 for r in results if r['is_callable'])
    het_count = sum(1 for r in results if r['is_callable'] and r['is_het'])

    exon_sites = sum(1 for r in results if r['region_type'] == 'exon')
    intron_sites = sum(1 for r in results if r['region_type'] == 'intron')

    print(f"Total sites in gene: {total}")
    print(f"Callable sites: {callable_count}")
    print(f"Heterozygous sites: {het_count}")
    print(f"Exon sites: {exon_sites}")
    print(f"Intron sites: {intron_sites}")
    print(f"Output: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Extract genotypes from gVCF for gene region'
    )
    parser.add_argument('--gvcf', required=True, help='Path to gVCF file')
    parser.add_argument('--gene_bed', required=True, help='Gene coordinates BED file')
    parser.add_argument('--exon_bed', required=True, help='Exon coordinates BED file')
    parser.add_argument('--intron_bed', required=True, help='Intron coordinates BED file')
    parser.add_argument('--output', required=True, help='Output TSV file')
    parser.add_argument('--min_dp', type=int, default=10, help='Minimum depth for callable site')
    parser.add_argument('--min_gq', type=int, default=20, help='Minimum GQ for callable site')

    args = parser.parse_args()

    # Create output directory if needed
    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)

    extract_genotypes(
        gvcf_file=args.gvcf,
        gene_bed=args.gene_bed,
        exon_bed=args.exon_bed,
        intron_bed=args.intron_bed,
        output_file=args.output,
        min_dp=args.min_dp,
        min_gq=args.min_gq
    )


if __name__ == "__main__":
    main()
