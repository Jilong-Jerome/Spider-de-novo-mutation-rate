import sys
# Define your input (filtered sites) and output (positions list) file paths
filtered_sites_file = sys.argv[1]#'/path/to/filtered_sites.tsv'  # Update this to the path of your filtered sites file
positions_list_file = sys.argv[2]#'/path/to/positions_list.txt'  # Update this to your desired positions list file path
child = sys.argv[3]
father = sys.argv[4]
mother = sys.argv[5]
# Create a positions list file from the filtered sites
with open(filtered_sites_file, 'r') as infile, open(positions_list_file, 'w') as outfile:
    for line in infile:
        cols = line.strip().split('\t')
        chrom, pos, ref, alt= cols[0], cols[1], cols[2],cols[3]
        # Write "chromosome:position" to the positions list file
        outfile.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{child}\t{father}\t{mother}\n")

