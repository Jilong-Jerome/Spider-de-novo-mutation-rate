import pysam
import os
import sys

def split_vcf_by_chromosome(input_vcf_path, output_dir,sp):
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Open the input VCF file (pysam automatically handles gzipped files with .gz extension)
    vcf_in = pysam.VariantFile(input_vcf_path, "r")
    
    # Create a dictionary to hold file handles for each chromosome
    chromosome_files = {}
    
    # Iterate through each record in the input VCF
    for record in vcf_in:
        chrom = record.chrom
        if chrom not in chromosome_files:
            # Create a new gzipped VCF file for this chromosome with mode 'w.gz'
            output_path = os.path.join(output_dir, "{sp}_{chrom}.vcf".format(sp=sp,chrom=chrom))
            chromosome_files[chrom] = pysam.VariantFile(output_path, "w", header=vcf_in.header)
        
        # Write the record to the appropriate chromosome VCF file
        chromosome_files[chrom].write(record)
    
    # Close all file handles
    for chrom, file in chromosome_files.items():
        file.close()
    # Close the input VCF file
    vcf_in.close()

input_vcf_path = sys.argv[1]#"path/to/your/input.vcf.gz"
output_dir = sys.argv[2]#"path/to/output/directory"
sp = sys.argv[3]
split_vcf_by_chromosome(input_vcf_path, output_dir, sp)

