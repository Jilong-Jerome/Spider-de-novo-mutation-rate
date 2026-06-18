import sys
def count_sites_in_chromosomes(vcf_file_path):
    chromosome_counts = {}
    with open(vcf_file_path, 'r') as vcf:
        for line in vcf:
            if line.startswith('##'):
                continue  # Skip header lines
            elif line.startswith('#'):
                # This is the header line with column names, could be used if needed
                pass
            else:
                # Extract chromosome ID
                chrom_id = line.split('\t')[0]
                if chrom_id in chromosome_counts:
                    chromosome_counts[chrom_id] += 1
                else:
                    chromosome_counts[chrom_id] = 1
    return chromosome_counts

# Example usage
vcf_file_path = sys.argv[1]#'path/to/your/file.vcf'
ind = sys.argv[2]# individual id
output = ind+"_callable_sites.tsv" # output file name for sites per chromosome
chromosome_counts = count_sites_in_chromosomes(vcf_file_path)
out = open(output,'w')
for chrom, count in chromosome_counts.items():
    out.write("{chrom}\t{count}\t{ind}\n".format(chrom=chrom,count=count,ind=ind))

