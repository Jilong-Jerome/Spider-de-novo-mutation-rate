from Bio import SeqIO
import sys

trim_info = sys.argv[1]
trim_fa = sys.argv[2]
output_fa = sys.argv[3]

# Step 1: Read trimming information from the tab-separated file
trimming_info = {}  # {chromosome: (start, end, annotation)}
with open(trim_info, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == 5:
            chromosome, _, start, end, annotation = parts
            start, end = int(start), int(end)
            trimming_info[chromosome] = (start, end, annotation)

# Step 2: Parse and trim the genome sequences
trimmed_sequences = []
for record in SeqIO.parse(trim_fa, 'fasta'):
    if record.id in trimming_info:
        start, end, _ = trimming_info[record.id]
        # Be sure to adjust the trimming positions as needed (e.g., 0-based or 1-based index)
        trimmed_seq = record.seq[:start-1] + record.seq[end:]
        record.seq = trimmed_seq
    trimmed_sequences.append(record)

# Step 3: Write the trimmed sequences to a new FASTA file
with open(output_fa, 'w') as output_handle:
    SeqIO.write(trimmed_sequences, output_handle, 'fasta')

