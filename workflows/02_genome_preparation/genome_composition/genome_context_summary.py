import sys
import pandas as pd
from Bio import SeqIO

def parse_bed(file_path):
    """Parse a BED file to get regions as tuples (chrom, start, end), adjusting for 0-based to 1-based conversion."""
    bed_df = pd.read_csv(file_path, sep='\t', header=None, names=['chrom', 'start', 'end'])
    return [(row['chrom'], row['start'] + 1, row['end']) for _, row in bed_df.iterrows()]  # Adjust start by +1

def load_genome_into_memory(fasta_file):
    """Load entire genome into memory as a dictionary with chromosome names as keys."""
    genome = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        genome[record.id] = record.seq.upper()  # Convert to uppercase for consistency
    return genome

def get_sequence_with_flanking(genome, chrom, start, end):
    """Retrieve a sequence region with one extra base on each side if available."""
    # Handle start and end boundaries to avoid indexing outside the sequence range
    if chrom not in genome:
        return None, False, False  # Chromosome not found
    
    chrom_sequence = genome[chrom]
    left_flank = max(start - 2, 0)  # 1-base before start (0-based)
    right_flank = min(end, len(chrom_sequence))  # 1-base after end
    region_sequence = chrom_sequence[left_flank:right_flank + 1]  # Include the region and potential flanking bases

    # Flags to check if left and right flanks are available
    left_flank_available = (start - 1) > 0
    right_flank_available = (end + 1) <= len(chrom_sequence)

    return region_sequence, left_flank_available, right_flank_available

def count_nucleotides_and_dinucleotides(sequence, left_flank_available, right_flank_available):
    """Count nucleotides and CG/GC dinucleotides, handling cases where flanking bases may not be available."""
    # If no left flank, start counting dinucleotides from the first base
    # If no right flank, stop dinucleotide counting before the last base
    dinucleotide_start = 1 if left_flank_available else 0
    dinucleotide_end = len(sequence) - 1 if right_flank_available else len(sequence)
    
    # Count individual nucleotides (excluding flanking bases for main nucleotide counts)
    main_sequence = sequence[1:-1] if left_flank_available and right_flank_available else sequence
    a_count = main_sequence.count("A")
    t_count = main_sequence.count("T")
    c_count = main_sequence.count("C")
    g_count = main_sequence.count("G")
    
    # Count dinucleotides
    C_followed_by_G = sum(1 for i in range(dinucleotide_start, dinucleotide_end - 1) if sequence[i] == "C" and sequence[i + 1] == "G")
    G_preceded_by_C = sum(1 for i in range(dinucleotide_start, dinucleotide_end - 1) if sequence[i] == "G" and sequence[i - 1] == "C")
    CG_dinucleotide_total = C_followed_by_G + G_preceded_by_C  # Total count of CG dinucleotides in both orientations
    
    return {
        "A": a_count, 
        "T": t_count, 
        "C": c_count, 
        "G": g_count, 
        "C_followed_by_G": C_followed_by_G,
        "G_preceded_by_C": G_preceded_by_C,
        "CG_dinucleotide_total": CG_dinucleotide_total,
        "total_bases": len(main_sequence)  # Count of main sequence bases only
    }

def main(bed_file, fasta_file, output_file, region_name):
    # Load the genome into memory
    genome = load_genome_into_memory(fasta_file)
    
    # Parse regions from BED file
    regions = parse_bed(bed_file)
    
    # Initialize counts
    total_counts = {"A": 0, "T": 0, "C": 0, "G": 0, "C_followed_by_G": 0, "G_preceded_by_C": 0, "CG_dinucleotide_total": 0, "total_bases": 0}
    
    # Process each region independently
    for chrom, start, end in regions:
        seq_with_flanking, left_flank_available, right_flank_available = get_sequence_with_flanking(genome, chrom, start, end)
        if seq_with_flanking:
            # Count nucleotides and dinucleotides, handling flanking cases
            counts = count_nucleotides_and_dinucleotides(seq_with_flanking, left_flank_available, right_flank_available)
            
            # Sum counts across all regions
            for base in total_counts:
                total_counts[base] += counts[base]
    
    # Calculate fractions
    fractions = {base: total_counts[base] / total_counts["total_bases"] for base in ["A", "T", "C", "G"]}
    fractions["C_followed_by_G"] = total_counts["C_followed_by_G"] / (total_counts["total_bases"] - 1) if total_counts["total_bases"] > 1 else 0
    fractions["G_preceded_by_C"] = total_counts["G_preceded_by_C"] / (total_counts["total_bases"] - 1) if total_counts["total_bases"] > 1 else 0
    fractions["CG_dinucleotide_total"] = total_counts["CG_dinucleotide_total"] / (total_counts["total_bases"] - 1) if total_counts["total_bases"] > 1 else 0
    
    # Write results to a tab-separated file
    with open(output_file, "w") as file:
        file.write("Region\tTotal_Sites\tNucleotide\tCount\tFraction\n")
        for base in ["A", "T", "C", "G"]:
            file.write(f"{region_name}\t{total_counts['total_bases']}\t{base}\t{total_counts[base]}\t{fractions[base]:.4f}\n")
        file.write(f"{region_name}\t{total_counts['total_bases']}\tC_followed_by_G\t{total_counts['C_followed_by_G']}\t{fractions['C_followed_by_G']:.4f}\n")
        file.write(f"{region_name}\t{total_counts['total_bases']}\tG_preceded_by_C\t{total_counts['G_preceded_by_C']}\t{fractions['G_preceded_by_C']:.4f}\n")
        file.write(f"{region_name}\t{total_counts['total_bases']}\tCG_dinucleotide_total\t{total_counts['CG_dinucleotide_total']}\t{fractions['CG_dinucleotide_total']:.4f}\n")
    
    print(f"Results written to {output_file}")


# Replace with your file paths and region name
bed_file_path = sys.argv[1]#"path_to_your_bed_file.bed"
fasta_file_path = sys.argv[2]#"path_to_your_genome_file.fasta"
output_file_path = sys.argv[3]#"output_fractions.tsv"
region_name = sys.argv[4]#"User_Specified_Region_Name"  # Replace with the region name you want to use
main(bed_file_path, fasta_file_path, output_file_path, region_name)

