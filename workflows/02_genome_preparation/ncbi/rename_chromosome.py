from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

def rename_and_reverse_complement(input_file, genome_fasta, output_fasta, all_output_fasta):
    # Parse the input file to get the renaming and reverse complement information
    rename_dict = {}
    with open(input_file, 'r') as infile:
        for line in infile:
            original_name, new_name, reverse_flag = line.strip().split('\t')
            rename_dict[original_name] = (new_name, reverse_flag)

    # Process the FASTA file
    new_records = []
    all_records = []
    for record in SeqIO.parse(genome_fasta, "fasta"):
        # Rename and reverse complement sequences specified in the input file
        if record.id in rename_dict:
            new_name, reverse_flag = rename_dict[record.id]
            if reverse_flag == '1':
                # Reverse complement the sequence
                new_seq = record.seq.reverse_complement()
            else:
                new_seq = record.seq
            # Rename the sequence
            new_record = SeqRecord(new_seq, id=new_name, description='')
            new_records.append(new_record)
            all_records.append(new_record)
        else:
            # Keep the sequence unchanged for the all_output_fasta
            all_records.append(record)

    # Write the processed sequences to the output FASTA files
    SeqIO.write(new_records, output_fasta, "fasta")
    SeqIO.write(all_records, all_output_fasta, "fasta")

if __name__ == "__main__":

    input_file = sys.argv[1]
    genome_fasta = sys.argv[2]
    output_fasta = sys.argv[3]
    all_output_fasta = sys.argv[4]
    rename_and_reverse_complement(input_file, genome_fasta, output_fasta, all_output_fasta)

