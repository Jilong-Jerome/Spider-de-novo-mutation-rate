import sys
input_vcf = sys.argv[1]#'/path/to/input.vcf'  # Update this with your VCF file path
output_vcf = sys.argv[2]#'/path/to/output_filtered.vcf'  # Update this with your desired output VCF file path
# Function to read variants and headers from the VCF file
def read_variants_and_headers_from_vcf(vcf_file):
    headers = []
    sites = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                headers.append(line)  # Store header lines
            elif line.startswith('#CHROM'):
                headers.append(line)  # Store the column header line
            else:
                parts = line.strip().split('\t')
                chrom, pos = parts[0], parts[1]
                sites.append((chrom, pos, line))  # Keep the entire line for writing back to file
    return headers, sites

# Check if a site is within 150bp of any other site
def is_within_distance(site, sites):
    chrom1, pos1 = site[0], int(site[1])
    for other_site in sites:
        chrom2, pos2 = other_site[0], int(other_site[1])
        if chrom1 == chrom2 and abs(pos1 - pos2) <= 150 and (chrom1, pos1) != (chrom2, pos2):
            return True
    return False

# Read the headers and sites from the input VCF
headers, sites = read_variants_and_headers_from_vcf(input_vcf)

# Filter sites based on distance
filtered_sites = [site for site in sites if not is_within_distance(site, sites)]

# Write the headers and filtered sites to a new VCF file
with open(output_vcf, 'w') as f:
    for header in headers:
        f.write(header)  # Write the header lines to the output VCF
    for site in filtered_sites:
        f.write(site[2])  # Write the original line from the VCF
