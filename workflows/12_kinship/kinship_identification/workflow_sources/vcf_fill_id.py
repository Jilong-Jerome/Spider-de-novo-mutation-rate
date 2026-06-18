import sys
def update_vcf_ids(input_vcf, output_vcf):
    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith("#"):
                # Write header lines as is
                outfile.write(line)
            else:
                # Split the line by tab to access fields
                fields = line.strip().split('\t')
                
                chrom = fields[0]  # Chromosome
                pos = fields[1]    # Position
                var_id = fields[2] # Variant ID

                # If the variant ID is ".", update it to chrom:pos
                if var_id == ".":
                    fields[2] = f"{chrom}:{pos}"

                # Write the updated line
                outfile.write("\t".join(fields) + "\n")


# Example usage
input_vcf = sys.argv[1]#"input_file.vcf"  # Replace with your input VCF file path
output_vcf = sys.argv[2]#"output_file.vcf"  # Replace with desired output VCF file path

update_vcf_ids(input_vcf, output_vcf)

