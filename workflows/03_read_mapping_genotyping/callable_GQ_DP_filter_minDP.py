import pandas as pd
import gzip
import os
import sys

import pandas as pd
import gzip
import os

# Load the mean depth information with per-chromosome detail
def load_mean_depth_info(filepath):
    mean_depth_info = pd.read_csv(filepath, sep='\t')
    mean_depth_info.set_index(['sample_id', 'chromosome'], inplace=True)
    return mean_depth_info

# Function to determine if a file should be opened as gzipped based on its extension
def open_file_autodetect(filepath, mode='rt'):
    if filepath.endswith('.gz'):
        return gzip.open(filepath, mode)
    else:
        return open(filepath, mode)

# Function to filter VCF lines based on depth, FILTER column, per-chromosome depth, and minimal GQ or RGQ
def filter_vcf_line(line, header, mean_depth_info, minDP):
    fields = line.split('\t')
    chromosome = fields[0]
    if fields[6] != "PASS":
        return False

    format_fields = fields[8].split(':')
    dp_index = format_fields.index('DP') if 'DP' in format_fields else None
    gq_index = format_fields.index('GQ') if 'GQ' in format_fields else None
    rgq_index = format_fields.index('RGQ') if 'RGQ' in format_fields else None

    if dp_index is None:
        return False  # Skip this site if DP is not present

    sample_fields = fields[9:]
    min_quality = float('inf')  # Initialize with infinity for comparison

    for sample_field, sample_header in zip(sample_fields, header[9:]):
        sample_info = sample_field.split(':')
        dp = int(sample_info[dp_index].strip("\n")) if sample_info[dp_index].strip("\n") != '.' else 0
        quality = 0  # Initialize quality score

        # Use GQ or RGQ if available, preferring GQ when both are present
        if gq_index is not None:
            quality = int(sample_info[gq_index].strip("\n")) if sample_info[gq_index].strip("\n") != '.' else 0
        elif rgq_index is not None:
            quality = int(sample_info[rgq_index].strip("\n")) if sample_info[rgq_index].strip("\n") != '.' else 0

        min_quality = min(min_quality, quality)  # Update min_quality if current quality is lower

        try:
            mean_depth = mean_depth_info.loc[(sample_header, chromosome), 'mean_depth']
        except KeyError:
            print("Mean depth info missing for sample {sample_header}, chromosome {chromosome}. Skipping...".format(sample_header=sample_header,chromosome=chromosome))
            return False

        if not ((0.5 * mean_depth <= dp <= 2 * mean_depth) and  (dp >= minDP)):
            return False  # Filter out this site if DP is not within range or smaller than minimum threshold

    if min_quality < 60:
        return False  # Filter out this site if any individual's quality (GQ or RGQ) is below 60

    return True  # Keep this site if all samples pass the filters

# Main function to filter the VCF file
def filter_vcf(vcf_filepath, mean_depth_filepath, output_filepath, minDP):
    mean_depth_info = load_mean_depth_info(mean_depth_filepath)

    with open_file_autodetect(vcf_filepath, 'rt') as vcf_file:
        header_lines = []
        for line in vcf_file:
            if line.startswith('##'):
                header_lines.append(line)
            elif line.startswith('#CHROM'):
                header_lines.append('##FILTER=<ID=PASS,Description="Passed all filters including DP range check based on individual chromosome mean depth, minimal quality (GQ or RGQ) >= 60, and FILTER column is PASS">\n')
                header_lines.append(line)
                break

        with open_file_autodetect(output_filepath, 'wt') as output_file:
            output_file.writelines(header_lines)

            for line in vcf_file:
                if filter_vcf_line(line, header_lines[-1].strip().split('\t'), mean_depth_info, minDP):
                    output_file.write(line)

# Example usage

vcf_filepath = sys.argv[1]#'path/to/your/input.vcf'  # This can be a .vcf or .vcf.gz file
mean_depth_filepath = sys.argv[2]#'path/to/your/mean_depth_info_per_chromosome.txt'
output_filepath = sys.argv[3]#'path/to/your/filtered_output.vcf'  # Change to .vcf.gz if desired
minDP_threshold = int(sys.argv[4])
filter_vcf(vcf_filepath, mean_depth_filepath, output_filepath, minDP)

