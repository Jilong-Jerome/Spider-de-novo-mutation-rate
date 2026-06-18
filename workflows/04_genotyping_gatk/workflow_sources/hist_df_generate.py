import sys

def calculate_ad_frequency_histogram_large_file_by_chromosome_minDP(file_path, output_file, bins=100):
    # Define bin edges and labels
    bin_edges = [i / bins for i in range(bins + 1)]
    bin_labels = [f"{round(bin_edges[i], 2)}-{round(bin_edges[i + 1], 2)}" for i in range(bins)]
    individual_bin_counts = {}

    # Open the file and process line-by-line
    with open(file_path, 'r') as file:
        for line in file:
            # Parse line (assuming tab-separated format)
            fields = line.strip().split('\t')
            individual = fields[0]
            chromosome = fields[1]
            minDP = int(fields[6])  # Minimum depth threshold
            ref_depth = int(fields[3])
            alt_depth = int(fields[4])
            
            # Calculate AD frequency and avoid division by zero
            if ref_depth + alt_depth > 0:
                ad_frequency = alt_depth / (alt_depth + ref_depth)
            else:
                continue  # Skip if both depths are zero
            
            # Determine bin index
            for i in range(bins):
                if bin_edges[i] <= ad_frequency < bin_edges[i + 1]:
                    bin_label = bin_labels[i]
                    break
            else:
                bin_label = bin_labels[-1]  # If ad_frequency == 1.0, it goes to the last bin

            # Initialize structure for individual, chromosome, and minDP if not present
            if individual not in individual_bin_counts:
                individual_bin_counts[individual] = {}
            if chromosome not in individual_bin_counts[individual]:
                individual_bin_counts[individual][chromosome] = {}
            if minDP not in individual_bin_counts[individual][chromosome]:
                individual_bin_counts[individual][chromosome][minDP] = {label: 0 for label in bin_labels}
            
            # Increment the bin count for the specific chromosome and minDP threshold
            individual_bin_counts[individual][chromosome][minDP][bin_label] += 1

    # Write the results to the output file
    with open(output_file, 'w') as out_file:
        out_file.write("individual\tchromosome\tminDP\tAD_frequency_bin\tcount\n")
        for individual, chrom_data in individual_bin_counts.items():
            for chromosome, dp_data in chrom_data.items():
                for minDP, bins in dp_data.items():
                    for bin_label, count in bins.items():
                        out_file.write(f"{individual}\t{chromosome}\t{minDP}\t{bin_label}\t{count}\n")

# Usage example (replace 'data_file.tsv' with the path to your actual data file)
ad_count = sys.argv[1]
output = sys.argv[2]
calculate_ad_frequency_histogram_large_file_by_chromosome_minDP(ad_count,output)

