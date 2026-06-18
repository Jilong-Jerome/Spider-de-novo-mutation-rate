import vcf
import sys
def read_vcf(file_path):
    vcf_reader = vcf.Reader(open(file_path, 'r'))
    return vcf_reader

def filter_and_write_vcf(query_vcf_path, database_vcf_path, output_vcf_path, max_distance=150):
    # Read both VCF files
    query_vcf = read_vcf(query_vcf_path)
    database_sites = list(read_vcf(database_vcf_path))  # Convert to list for easier processing

    # Create a VCF writer for the output file
    output_vcf = vcf.Writer(open(output_vcf_path, 'w'), query_vcf)

    for query_site in query_vcf:
        keep_site = True
        for db_site in database_sites:
            distance = abs(query_site.POS - db_site.POS)
            if distance <= max_distance and distance != 0:  # Ensure it's not the same site
                keep_site = False
                break
        if keep_site:
            output_vcf.write_record(query_site)

    output_vcf.close()

# Paths to your VCF files
query_vcf_path = sys.argv[1]#'path/to/query_vcf.vcf'
database_vcf_path = sys.argv[2]#'path/to/database_vcf.vcf'
output_vcf_path = sys.argv[3]#'path/to/output_vcf.vcf'

# Filter query VCF and write results to a new file
filter_and_write_vcf(query_vcf_path, database_vcf_path, output_vcf_path)

