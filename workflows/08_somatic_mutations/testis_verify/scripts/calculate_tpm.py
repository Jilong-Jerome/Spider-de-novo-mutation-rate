#!/usr/bin/env python3
"""
Calculate TPM (Transcripts Per Million) from featureCounts output.

Merges count files from multiple samples and normalizes to TPM.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import numpy as np


def parse_featurecounts(count_file: Path) -> pd.DataFrame:
    """Parse featureCounts output file.

    Returns DataFrame with columns: gene_id, length, count
    """
    # featureCounts output has header lines starting with #
    # and a header row with column names
    df = pd.read_csv(count_file, sep='\t', comment='#')

    # Columns are: Geneid, Chr, Start, End, Strand, Length, <sample>
    # Rename for clarity
    sample_name = df.columns[-1]
    result = pd.DataFrame({
        'gene_id': df['Geneid'],
        'length': df['Length'],
        'count': df[sample_name]
    })

    return result


def calculate_tpm(counts: pd.Series, lengths: pd.Series) -> pd.Series:
    """Calculate TPM from counts and gene lengths.

    TPM = (count / length) * 1e6 / sum(count / length)
    """
    # Reads per kilobase
    rpk = counts / (lengths / 1000)

    # Scale to million
    scaling_factor = rpk.sum() / 1e6

    # TPM
    tpm = rpk / scaling_factor

    return tpm


def main():
    parser = argparse.ArgumentParser(
        description='Calculate TPM from featureCounts output and merge samples'
    )
    parser.add_argument('--counts', required=True, nargs='+', type=Path,
                        help='featureCounts output files')
    parser.add_argument('--sample-names', required=True, nargs='+',
                        help='Sample names (same order as count files)')
    parser.add_argument('--output', required=True, type=Path,
                        help='Output TPM matrix file')

    args = parser.parse_args()

    if len(args.counts) != len(args.sample_names):
        print("Error: Number of count files must match number of sample names")
        return 1

    print(f"Processing {len(args.counts)} samples...")

    # Parse all count files
    all_data = {}
    gene_lengths = None

    for count_file, sample_name in zip(args.counts, args.sample_names):
        print(f"  Reading {sample_name} from {count_file}...")
        df = parse_featurecounts(count_file)

        # Store gene lengths (should be same across samples)
        if gene_lengths is None:
            gene_lengths = df.set_index('gene_id')['length']

        all_data[sample_name] = df.set_index('gene_id')['count']

    # Create count matrix
    count_matrix = pd.DataFrame(all_data)
    count_matrix['length'] = gene_lengths

    print(f"\nTotal genes: {len(count_matrix)}")

    # Calculate TPM for each sample
    tpm_data = {}
    for sample_name in args.sample_names:
        counts = count_matrix[sample_name]
        lengths = count_matrix['length']

        # Handle zero lengths
        valid_mask = lengths > 0
        tpm = pd.Series(0.0, index=count_matrix.index)
        tpm[valid_mask] = calculate_tpm(counts[valid_mask], lengths[valid_mask])

        tpm_data[sample_name] = tpm

        # Print summary stats
        expressed = (tpm > 1).sum()
        print(f"  {sample_name}: {expressed} genes with TPM > 1")

    # Create TPM matrix
    tpm_matrix = pd.DataFrame(tpm_data)
    tpm_matrix['length'] = gene_lengths

    # Reorder columns
    cols = ['length'] + list(args.sample_names)
    tpm_matrix = tpm_matrix[cols]

    # Save
    args.output.parent.mkdir(parents=True, exist_ok=True)
    tpm_matrix.to_csv(args.output, sep='\t', float_format='%.2f')

    print(f"\nTPM matrix written to {args.output}")

    # Also save raw counts
    count_output = args.output.parent / 'count_matrix.tsv'
    count_matrix.to_csv(count_output, sep='\t')
    print(f"Count matrix written to {count_output}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
