import sys
import pysam

def extract_allele_depths(vcf_file, output_file, target_individuals, sp, minDP):
    # Open VCF file using pysam
    vcf = pysam.VariantFile(vcf_file)

    # Prepare output file
    with open(output_file, 'w') as out_file:

        # Iterate through each variant in the VCF file
        for record in vcf.fetch():
            chrom = record.chrom
            pos = record.pos

            # Process each specified individual
            for ind_id in target_individuals:
                if ind_id in record.samples:
                    sample = record.samples[ind_id]
                    gt = sample['GT']  # Genotype tuple (e.g., (0, 1) for heterozygous)

                    # Check if site is heterozygous
                    if gt == (0, 1) or gt == (1, 0):
                        # Extract allele depth information
                        if 'AD' in sample:
                            ref_depth, alt_depth = sample['AD']
                            # Write to output file
                            out_file.write(f"{ind_id}\t{chrom}\t{pos}\t{ref_depth}\t{alt_depth}\t{sp}\t{minDP}\n")

    print(f"Results saved to {output_file}")

# Example usage
vcf_file = sys.argv[1] # "path_to_your_vcf_file.vcf"
output_file = sys.argv[2] # "allele_depths.tsv"
female = sys.argv[3]
male = sys.argv[4]
offspring = sys.argv[5]
sp = sys.argv[6]
minDP = sys.argv[7]
target_individuals = [female,male,offspring]  # Replace with your specific individuals

extract_allele_depths(vcf_file, output_file, target_individuals,sp,minDP)

