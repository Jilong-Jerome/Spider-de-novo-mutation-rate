def custom_chromosome_sort(chromosome_name):
    """
    Custom sort function for chromosome names.
    Sorts based on the identifier after the species code and underscore, handling numerical and 'X' suffixed identifiers.
    """
    # Split by underscore and take the part after the species code
    identifier = chromosome_name.split('_')[1]
    if identifier.startswith('X'):
        # Assign a high number to 'X' suffixed identifiers to ensure they come after numerical ones
        return int(identifier[1:]) + 1000
    else:
        return int(identifier)

def generate_intervals_and_save(fasta_index_path, segment_size, species, output_path):
    """
    Reads a fasta index file, sorts chromosomes by their identifiers, generates intervals for GATK,
    and saves them along with database names.
    
    Parameters:
    - fasta_index_path: Path to the fasta index file.
    - segment_size: Size of each interval in bases (e.g., 5Mb = 5,000,000).
    - species_prefix: Prefix of the species code to be used in database naming.
    - output_path: Path to save the output file with intervals and database names.
    """
    intervals = []
    db_index = 0  # To keep track of database naming

    fasta_index_content = []
    with open(fasta_index_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chr_name, chr_length_str = parts[0], parts[1]
                fasta_index_content.append((chr_name, int(chr_length_str)))

    # Sort the chromosomes based on custom logic, focusing on the identifier after the species code
    sorted_chromosomes = sorted(fasta_index_content, key=lambda x: custom_chromosome_sort(x[0]))
    
    for chr_name, chr_length in sorted_chromosomes:
        start = 1
        while start < chr_length:
            end = min(start + segment_size - 1, chr_length)
            interval_str = "{chr_name}:{start}-{end}".format(chr_name=chr_name,start=start,end=end)
            db_name = "{species}_{index}".format(species=species,index = str(db_index).zfill(3))
            intervals.append((interval_str, db_name))
            
            start += segment_size
            db_index += 1  # Increment database index for each interval
    
    # Save the intervals and database names to the output file
    with open(output_path, 'w') as out_f:
        for interval, db_name in intervals:
            out_f.write("{interval}\t{db_name}\n".format(interval=interval,db_name=db_name))



# Example usage parameters
#fasta_index_path = "path/to/your/fasta_index_file.fai"
#segment_size = 5000000  # 5Mb
#species = "your_species"
#output_path = "path/to/your/output_intervals.txt"

# Uncomment the line below to use the function with the example parameters
# generate_intervals_and_save(fasta_index_path, segment_size, species, output_path)

