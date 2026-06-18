import pandas as pd
import sys

input_file = sys.argv[1]
output = sys.argv[2]

def calculate_statistics(values):
    """Calculate mean, median, and quantiles for a list of values."""
    series = pd.Series(values)
    return {
        'mean': series.mean(),
        'median': series.median(),
        'quantile_0.025': series.quantile(0.025),
        'quantile_0.25': series.quantile(0.25),
        'quantile_0.75': series.quantile(0.75),
        'quantile_0.975': series.quantile(0.975)
    }

# Initialize variables
current_chromosome = None
values = []

# Open the output file for writing
with open(output, 'w') as output_file:
    # Write the header
    output_file.write("chromosome\tmean\tmedian\tquantile0.025\tquantile0.975\tquantile0.25\tquantile0.75\n")

    # Read the input file line by line
    with open(input_file, 'r') as file:
        for line in file:
            chromosome, _, depth = line.strip().split('\t')
            depth = int(depth)

            # Check if the chromosome changes
            if chromosome != current_chromosome and current_chromosome is not None:
                # Calculate statistics for the previous chromosome
                stats = calculate_statistics(values)
                # Write the statistics to the output file
                output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(current_chromosome, stats['mean'], stats['median'], stats['quantile_0.025'], stats['quantile_0.975'], stats['quantile_0.25'], stats['quantile_0.75']))
                values = []

            # Update current chromosome and values
            current_chromosome = chromosome
            values.append(depth)

    # Handle the last chromosome
    if current_chromosome is not None:
        stats = calculate_statistics(values)
        output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(current_chromosome, stats['mean'], stats['median'], stats['quantile_0.025'], stats['quantile_0.975'], stats['quantile_0.25'], stats['quantile_0.75']))

# The output file 'output_statistics.tsv' now contains the required statistics for each chromosome

