#!/usr/bin/env python3
"""
parse_gff3.py - Extract gene, exon, and intron coordinates from GFF3 file

This script parses a GFF3 annotation file to extract coordinates for a specific
gene/transcript, outputting BED files for the gene boundaries, exons, and introns.

Introns are defined strictly as regions between consecutive exons, excluding UTRs.
"""

import argparse
import os
import sys
from collections import defaultdict


def parse_gff3_attributes(attr_string):
    """Parse GFF3 attribute string into dictionary"""
    attrs = {}
    for item in attr_string.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            attrs[key.strip()] = value.strip()
    return attrs


def parse_gff3(gff3_file, gene_id):
    """
    Parse GFF3 file and extract coordinates for a specific gene/transcript.

    Args:
        gff3_file: Path to GFF3 file
        gene_id: Gene/transcript ID to search for (e.g., AFR.HiC_4_file_1_file_1_jg528.t1)

    Returns:
        dict with 'gene', 'exons', 'introns' coordinates
    """
    gene_info = None
    exons = []

    # The gene_id might be a transcript ID (ending with .t1, .t2, etc.)
    # We need to find exons with Parent matching this transcript ID

    with open(gff3_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            chrom, source, feature_type, start, end, score, strand, phase, attributes = parts
            attrs = parse_gff3_attributes(attributes)

            start = int(start)
            end = int(end)

            # Check if this is the mRNA/transcript we're looking for
            if feature_type == 'mRNA' and attrs.get('ID') == gene_id:
                gene_info = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'id': gene_id
                }

            # Check if this is an exon belonging to our transcript
            if feature_type == 'exon' and attrs.get('Parent') == gene_id:
                exons.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand
                })

    # If we didn't find mRNA entry, try to infer from exons
    if gene_info is None and exons:
        gene_info = {
            'chrom': exons[0]['chrom'],
            'start': min(e['start'] for e in exons),
            'end': max(e['end'] for e in exons),
            'strand': exons[0]['strand'],
            'id': gene_id
        }

    if gene_info is None:
        raise ValueError(f"Gene/transcript '{gene_id}' not found in GFF3 file")

    if not exons:
        raise ValueError(f"No exons found for gene/transcript '{gene_id}'")

    # Sort exons by start position
    exons.sort(key=lambda x: x['start'])

    # Calculate introns (strictly between consecutive exons)
    introns = []
    for i in range(len(exons) - 1):
        intron_start = exons[i]['end'] + 1
        intron_end = exons[i + 1]['start'] - 1

        if intron_start <= intron_end:
            introns.append({
                'chrom': exons[i]['chrom'],
                'start': intron_start,
                'end': intron_end,
                'strand': exons[i]['strand']
            })

    return {
        'gene': gene_info,
        'exons': exons,
        'introns': introns
    }


def write_bed(features, output_file, feature_name="feature"):
    """Write features to BED format file"""
    with open(output_file, 'w') as f:
        if isinstance(features, dict):
            # Single feature (gene)
            f.write(f"{features['chrom']}\t{features['start']}\t{features['end']}\t{features.get('id', feature_name)}\t0\t{features.get('strand', '+')}\n")
        else:
            # List of features (exons/introns)
            for i, feat in enumerate(features, 1):
                f.write(f"{feat['chrom']}\t{feat['start']}\t{feat['end']}\t{feature_name}_{i}\t0\t{feat.get('strand', '+')}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Extract gene/exon/intron coordinates from GFF3 file'
    )
    parser.add_argument('--gff3', required=True, help='Path to GFF3 annotation file')
    parser.add_argument('--gene_id', required=True, help='Gene/transcript ID to extract')
    parser.add_argument('--output_dir', required=True, help='Output directory for BED files')

    args = parser.parse_args()

    # Create output directory if needed
    os.makedirs(args.output_dir, exist_ok=True)

    # Parse GFF3
    try:
        coords = parse_gff3(args.gff3, args.gene_id)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Write output files
    gene_bed = os.path.join(args.output_dir, 'gene_coords.bed')
    exon_bed = os.path.join(args.output_dir, 'exon_coords.bed')
    intron_bed = os.path.join(args.output_dir, 'intron_coords.bed')

    write_bed(coords['gene'], gene_bed, "gene")
    write_bed(coords['exons'], exon_bed, "exon")
    write_bed(coords['introns'], intron_bed, "intron")

    # Print summary
    print(f"Gene: {coords['gene']['chrom']}:{coords['gene']['start']}-{coords['gene']['end']}")
    print(f"  Length: {coords['gene']['end'] - coords['gene']['start'] + 1} bp")
    print(f"  Exons: {len(coords['exons'])}")
    total_exon_length = sum(e['end'] - e['start'] + 1 for e in coords['exons'])
    print(f"  Total exon length: {total_exon_length} bp")
    print(f"  Introns: {len(coords['introns'])}")
    total_intron_length = sum(i['end'] - i['start'] + 1 for i in coords['introns'])
    print(f"  Total intron length: {total_intron_length} bp")

    print(f"\nOutput files:")
    print(f"  {gene_bed}")
    print(f"  {exon_bed}")
    print(f"  {intron_bed}")


if __name__ == "__main__":
    main()
