import sys
# Define the path to your input and output files
input_file = sys.argv[1]#'/path/to/input.tsv'  # Update this to your input file path
output_file = sys.argv[2]#'/path/to/output.tsv'  # Update this to your desired output file path

# Function to check the criteria for filtering
def passes_criteria(child_AD, father_AD, mother_AD):
    # Criteria (a): checked before this function is called

    # Criteria (b): Check if there's support for the reference allele in the child
    child_ref, child_alt = map(int, child_AD.split(","))
    if child_ref == 0:
        return False
    
    # Criteria (c): Check if there's no support for the alternative allele in the parents
    father_ref, father_alt = map(int, father_AD.split(","))
    mother_ref, mother_alt = map(int, mother_AD.split(","))
    if father_alt > 0 or mother_alt > 0:
        return False
    
    return True

# Open the input file for reading and the output file for writing
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        # Split the line into columns
        cols = line.strip().split('\t')
        
        # Extract relevant columns
        chrom, pos, ref, alt, child_DP, child_AD, father_DP, father_AD, mother_DP, mother_AD = cols
        
        # Criteria (a): Filter out if the site is not a SNP (i.e., if the alternative allele is ".")
        if alt == ".":
            continue
        
        # Apply the other criteria
        if passes_criteria(child_AD, father_AD, mother_AD):
            # Write the passing sites to the output file
            outfile.write(line)

