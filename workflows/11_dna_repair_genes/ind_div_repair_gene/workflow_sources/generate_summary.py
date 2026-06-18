#!/usr/bin/env python3
"""
generate_summary.py - Compile summary reports from all analysis results

This script concatenates individual heterozygosity, Dxy, and het site files
into species-wide summary tables.
"""

import argparse
import os
import sys
import glob


def concatenate_files(input_files, output_file, add_source=False):
    """
    Concatenate multiple TSV files into one, keeping header from first file only.

    Args:
        input_files: List of input file paths
        output_file: Output file path
        add_source: If True, add 'source_file' column
    """
    if not input_files:
        print(f"Warning: No input files provided for {output_file}")
        # Create empty file with placeholder
        with open(output_file, 'w') as f:
            f.write("# No data files found\n")
        return 0

    header_written = False
    total_lines = 0

    with open(output_file, 'w') as out:
        for input_file in input_files:
            if not os.path.exists(input_file):
                print(f"Warning: File not found: {input_file}")
                continue

            with open(input_file, 'r') as inp:
                header = inp.readline()

                if not header_written:
                    if add_source:
                        out.write(header.strip() + "\tsource_file\n")
                    else:
                        out.write(header)
                    header_written = True

                for line in inp:
                    if line.strip():
                        if add_source:
                            source = os.path.basename(input_file)
                            out.write(line.strip() + f"\t{source}\n")
                        else:
                            out.write(line)
                        total_lines += 1

    return total_lines


def calculate_total_heterozygosity(het_files, output_file):
    """
    Calculate cross-gene aggregated heterozygosity per individual.

    For each individual, for each region (whole_gene, exon_only, intron_only):
    - total_het_sites = sum(het_sites across all genes)
    - total_callable_sites = sum(callable_sites across all genes)
    - total_heterozygosity = total_het_sites / total_callable_sites

    Args:
        het_files: List of per-gene heterozygosity TSV files
        output_file: Output file for total heterozygosity

    Input format (per-gene):
        individual  gene  region  het_sites  callable_sites  heterozygosity

    Output format:
        individual  region  total_het_sites  total_callable_sites  total_heterozygosity
    """
    # Aggregate data: {individual: {region: {'het': 0, 'callable': 0}}}
    aggregated = {}

    for het_file in het_files:
        if not os.path.exists(het_file):
            print(f"Warning: File not found: {het_file}")
            continue

        with open(het_file, 'r') as f:
            header = f.readline().strip().split('\t')
            col_idx = {name: i for i, name in enumerate(header)}

            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < len(header):
                    continue

                individual = parts[col_idx['individual']]
                region = parts[col_idx['region']]
                het_sites = int(parts[col_idx['het_sites']])
                callable_sites = int(parts[col_idx['callable_sites']])

                if individual not in aggregated:
                    aggregated[individual] = {}
                if region not in aggregated[individual]:
                    aggregated[individual][region] = {'het': 0, 'callable': 0}

                aggregated[individual][region]['het'] += het_sites
                aggregated[individual][region]['callable'] += callable_sites

    # Write output
    with open(output_file, 'w') as f:
        f.write("individual\tregion\ttotal_het_sites\ttotal_callable_sites\ttotal_heterozygosity\n")

        for individual in sorted(aggregated.keys()):
            for region in ['whole_gene', 'exon_only', 'intron_only']:
                if region not in aggregated[individual]:
                    continue

                total_het = aggregated[individual][region]['het']
                total_callable = aggregated[individual][region]['callable']

                if total_callable > 0:
                    total_heterozygosity = total_het / total_callable
                    het_str = f"{total_heterozygosity:.6f}"
                else:
                    het_str = "NA"

                f.write(f"{individual}\t{region}\t{total_het}\t{total_callable}\t{het_str}\n")

    return len(aggregated)


def calculate_total_dxy(dxy_files, output_file):
    """
    Calculate cross-gene aggregated Dxy per family pair.

    For each family, for each region (whole_gene, exon_only, intron_only):
    - total_callable_sites_both = sum(callable_sites_both across all genes)
    - total_dxy_sum = sum(dxy_sum across all genes)
    - total_dxy = total_dxy_sum / total_callable_sites_both

    Args:
        dxy_files: List of per-gene Dxy TSV files
        output_file: Output file for total Dxy

    Input format (per-gene):
        gene  family  region  callable_sites_both  dxy_sum  dxy

    Output format:
        family  region  total_callable_sites_both  total_dxy_sum  total_dxy
    """
    # Aggregate data: {family: {region: {'callable': 0, 'dxy_sum': 0.0}}}
    aggregated = {}

    for dxy_file in dxy_files:
        if not os.path.exists(dxy_file):
            print(f"Warning: File not found: {dxy_file}")
            continue

        with open(dxy_file, 'r') as f:
            header = f.readline().strip().split('\t')
            col_idx = {name: i for i, name in enumerate(header)}

            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < len(header):
                    continue

                family = parts[col_idx['family']]
                region = parts[col_idx['region']]
                callable_sites_both = int(parts[col_idx['callable_sites_both']])
                dxy_sum = float(parts[col_idx['dxy_sum']])

                if family not in aggregated:
                    aggregated[family] = {}
                if region not in aggregated[family]:
                    aggregated[family][region] = {'callable': 0, 'dxy_sum': 0.0}

                aggregated[family][region]['callable'] += callable_sites_both
                aggregated[family][region]['dxy_sum'] += dxy_sum

    # Write output
    with open(output_file, 'w') as f:
        f.write("family\tregion\ttotal_callable_sites_both\ttotal_dxy_sum\ttotal_dxy\n")

        for family in sorted(aggregated.keys()):
            for region in ['whole_gene', 'exon_only', 'intron_only']:
                if region not in aggregated[family]:
                    continue

                total_callable = aggregated[family][region]['callable']
                total_dxy_sum = aggregated[family][region]['dxy_sum']

                if total_callable > 0:
                    total_dxy = total_dxy_sum / total_callable
                    dxy_str = f"{total_dxy:.6f}"
                else:
                    dxy_str = "NA"

                dxy_sum_str = f"{total_dxy_sum:.4f}"
                f.write(f"{family}\t{region}\t{total_callable}\t{dxy_sum_str}\t{dxy_str}\n")

    return len(aggregated)


def main():
    parser = argparse.ArgumentParser(
        description='Generate summary reports from analysis results'
    )
    parser.add_argument('--species', required=True, help='Species identifier')
    parser.add_argument('--het_files', nargs='+', required=True,
                        help='Heterozygosity TSV files')
    parser.add_argument('--dxy_files', nargs='+', required=True,
                        help='Dxy TSV files')
    parser.add_argument('--sites_files', nargs='+', required=True,
                        help='Het sites TSV files')
    parser.add_argument('--output_het', required=True,
                        help='Output heterozygosity summary TSV')
    parser.add_argument('--output_dxy', required=True,
                        help='Output Dxy summary TSV')
    parser.add_argument('--output_sites', required=True,
                        help='Output het sites catalog TSV')
    parser.add_argument('--output_total_het', required=True,
                        help='Output individual total heterozygosity TSV')
    parser.add_argument('--output_total_dxy', required=True,
                        help='Output family total Dxy TSV')

    args = parser.parse_args()

    # Create output directories
    for output in [args.output_het, args.output_dxy, args.output_sites,
                   args.output_total_het, args.output_total_dxy]:
        output_dir = os.path.dirname(output)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

    print(f"Generating summary reports for species: {args.species}")
    print("=" * 60)

    # Heterozygosity summary
    print(f"\nProcessing heterozygosity files ({len(args.het_files)} files)...")
    het_lines = concatenate_files(args.het_files, args.output_het)
    print(f"  Total records: {het_lines}")
    print(f"  Output: {args.output_het}")

    # Dxy summary
    print(f"\nProcessing Dxy files ({len(args.dxy_files)} files)...")
    dxy_lines = concatenate_files(args.dxy_files, args.output_dxy)
    print(f"  Total records: {dxy_lines}")
    print(f"  Output: {args.output_dxy}")

    # Het sites catalog
    print(f"\nProcessing het sites files ({len(args.sites_files)} files)...")
    sites_lines = concatenate_files(args.sites_files, args.output_sites)
    print(f"  Total records: {sites_lines}")
    print(f"  Output: {args.output_sites}")

    # Cross-gene aggregated heterozygosity per individual
    print(f"\nCalculating total heterozygosity per individual...")
    num_individuals = calculate_total_heterozygosity(args.het_files, args.output_total_het)
    print(f"  Total individuals: {num_individuals}")
    print(f"  Output: {args.output_total_het}")

    # Cross-gene aggregated Dxy per family
    print(f"\nCalculating total Dxy per family...")
    num_families = calculate_total_dxy(args.dxy_files, args.output_total_dxy)
    print(f"  Total families: {num_families}")
    print(f"  Output: {args.output_total_dxy}")

    print("\n" + "=" * 60)
    print("Summary generation complete!")


if __name__ == "__main__":
    main()
