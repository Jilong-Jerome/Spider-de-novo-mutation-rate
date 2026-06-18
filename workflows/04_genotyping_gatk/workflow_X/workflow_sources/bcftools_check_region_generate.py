# Import the necessary library
import sys

# Check for proper usage and command line arguments
if len(sys.argv) != 3:
    print("Usage: python script.py input.vcf output.tsv")
    sys.exit(1)

input_vcf = sys.argv[1]
output_tsv = sys.argv[2]

try:
    with open(input_vcf, 'r') as vcf, open(output_tsv, 'w') as output:
        # Write the header to the output file
        #output.write("CHROM\tBEG\tEND\n")

        for line in vcf:
            # Skip header lines
            if line.startswith("#"):
                continue
            
            # Split the line by tab to extract fields
            fields = line.strip().split('\t')
            
            # Extract the chromosome and position
            chrom = fields[0]
            pos = int(fields[1])
            
            # Calculate position-150 and position+150
            pos_minus_150 = pos - 150
            pos_plus_150 = pos + 150
            
            # Write the output line
            output_line = f"{chrom}\t{pos_minus_150}\t{pos_plus_150}\n"
            output.write(output_line)

    print("Processing complete. Output saved to:", output_tsv)
except Exception as e:
    print("An error occurred:", str(e))

