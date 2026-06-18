import sys
from Bio import SeqIO

def calculate_N50(contig_lengths):
    sorted_lengths = sorted(contig_lengths, reverse=True)
    total_length = sum(sorted_lengths)
    cumulative_length = 0
    for length in sorted_lengths:
        cumulative_length += length
        if cumulative_length >= total_length / 2:
            return length
    return 0

def process_fasta(input_file, output_file):
    genome_length = 0
    genome_contig_lengths = []
    total_contigs = 0
    chromosome_data = []

    with open(input_file, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            chromosome = record.id
            sequence = str(record.seq).upper()
            length = len(sequence.replace('N', ''))
            contig_lengths = [len(contig) for contig in sequence.split('N') if contig]
            num_contigs = len(contig_lengths)
            n50 = calculate_N50(contig_lengths)

            # Store chromosome data
            chromosome_data.append((chromosome, length, num_contigs, n50))

            # Accumulate genome statistics
            genome_length += length
            total_contigs += num_contigs
            genome_contig_lengths.extend(contig_lengths)

        # Sort chromosome data by name
        chromosome_data.sort(key=lambda x: x[0])

        # Write sorted data to output
        with open(output_file, 'w') as output:
            output.write("Chromosome\tLength\tContigs\tN50\n")
            for data in chromosome_data:
                output.write(f"{data[0]}\t{data[1]}\t{data[2]}\t{data[3]}\n")

            # Genome N50 and contigs count
            genome_n50 = calculate_N50(genome_contig_lengths)

            # Append genome statistics under "Whole Genome"
            output.write(f"Whole Genome\t{genome_length}\t{total_contigs}\t{genome_n50}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_fasta> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    process_fasta(input_file, output_file)

