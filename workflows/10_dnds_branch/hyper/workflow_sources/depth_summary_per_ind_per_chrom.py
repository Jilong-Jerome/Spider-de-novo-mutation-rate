#!/usr/bin/env python3

import pysam
import numpy as np
import sys
from collections import defaultdict

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} input.vcf.gz")
    sys.exit(1)

vcf_file = sys.argv[1]
out_file = sys.argv[2]
# Open gzipped VCF
vcf_in = pysam.VariantFile(vcf_file)

# Extract sample names
samples = list(vcf_in.header.samples)

# Get all contigs from header
chromosomes = list(vcf_in.header.contigs)

# Output file
output_file = out_file #vcf_file + ".depth_summary.tsv"

with open(output_file, "w") as out:
    out.write("individual\tchromosome\tmean_coverage\tstd_coverage\n")

    for chrom in chromosomes:
        print(f"Processing chromosome: {chrom}", file=sys.stderr)
        
        # Reset per chromosome
        depths = defaultdict(list)
        record_counter = 0

        # Fetch records only from current chromosome (indexed access)
        for record in vcf_in.fetch(chrom):
            for sample in samples:
                sample_data = record.samples[sample]
                dp = sample_data.get("DP", None)
                if dp is not None:
                    depths[sample].append(dp)
            record_counter += 1
            if record_counter % 100000 == 0:
                print(f"  Processed {record_counter} records on {chrom}...", file=sys.stderr)

        # Write results for this chromosome
        for sample in samples:
            dp_list = depths[sample]
            if dp_list:
                mean_dp = np.mean(dp_list)
                std_dp = np.std(dp_list)
            else:
                mean_dp = 0
                std_dp = 0
            out.write(f"{sample}\t{chrom}\t{mean_dp:.2f}\t{std_dp:.2f}\n")

vcf_in.close()

print(f"All chromosomes processed. Output written to {output_file}")

