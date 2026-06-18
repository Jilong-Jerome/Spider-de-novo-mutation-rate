#!/usr/bin/env python3
import argparse
from Bio import SeqIO

def is_cg_site(seq, pos):
    """
    Determine if the site at index pos (0-indexed) in the sequence is a CG dinucleotide site.
    A site qualifies if:
      - It is a C and either the preceding or following base is a G, or
      - It is a G and either the preceding or following base is a C.
    """
    base = seq[pos].upper()
    seq_len = len(seq)
    if base == 'C':
        if (pos > 0 and seq[pos-1].upper() == 'G') or (pos < seq_len - 1 and seq[pos+1].upper() == 'G'):
            return True
    elif base == 'G':
        if (pos > 0 and seq[pos-1].upper() == 'C') or (pos < seq_len - 1 and seq[pos+1].upper() == 'C'):
            return True
    return False

def is_cpd_site(seq, pos):
    """
    Determine if the site at index pos (0-indexed) qualifies as a CPD site.
    
    For a C site: count as CPD if the left or right neighbor is in (C, T).
    For a G site (reverse complement of C): count as CPD if the left or right neighbor is in (G, A).
    """
    base = seq[pos].upper()
    seq_len = len(seq)
    
    if base == 'C':
        neighbor_set = {'C', 'T'}
    elif base == 'G':
        neighbor_set = {'G', 'A'}
    else:
        return False

    left_neighbor = seq[pos-1].upper() if pos > 0 else None
    right_neighbor = seq[pos+1].upper() if pos < seq_len - 1 else None
    if (left_neighbor in neighbor_set) or (right_neighbor in neighbor_set):
        return True
    return False

def process_vcf(vcf_file, chrom_seq):
    """
    Process a VCF file and count callable sites based on the provided chromosome sequence.
    Returns:
      total_callable, count_A, count_T, count_C, count_G, count_CG, count_CPD
    """
    total_callable = 0
    count_A = 0
    count_T = 0
    count_C = 0
    count_G = 0
    count_CG = 0
    count_CPD = 0

    with open(vcf_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            # VCF positions are 1-indexed; convert to 0-indexed
            pos = int(parts[1])
            idx = pos - 1
            if idx < 0 or idx >= len(chrom_seq):
                continue  # Skip positions outside the chromosome
            base = chrom_seq[idx].upper()
            if base not in ('A', 'T', 'C', 'G'):
                continue  # Skip ambiguous bases
            total_callable += 1
            if base == 'A':
                count_A += 1
            elif base == 'T':
                count_T += 1
            elif base == 'C':
                count_C += 1
            elif base == 'G':
                count_G += 1

            if is_cg_site(chrom_seq, idx):
                count_CG += 1
            if is_cpd_site(chrom_seq, idx):
                count_CPD += 1

    return total_callable, count_A, count_T, count_C, count_G, count_CG, count_CPD

def load_genome(genome_fasta):
    """
    Load the genome FASTA file into a dictionary.
    Keys are the record IDs (chromosome names) and values are sequences.
    """
    genome_dict = {}
    for record in SeqIO.parse(genome_fasta, "fasta"):
        genome_dict[record.id] = record.seq
    return genome_dict

def main():
    parser = argparse.ArgumentParser(
        description="Summarise nucleotide composition from a VCF file using genome context."
    )
    parser.add_argument("-g", "--genome", required=True, help="Genome FASTA file")
    parser.add_argument("-s", "--species", required=True, help="Species name")
    parser.add_argument("-i", "--ind", required=True, help="Individual ID")
    parser.add_argument("-c", "--chromosome", required=True, help="Chromosome name to be used for output")
    parser.add_argument("-o", "--output", required=True, help="Output summary file")
    parser.add_argument("vcf_file", help="VCF file for the specified chromosome")
    args = parser.parse_args()

    # Load the genome and get the sequence for the specified chromosome
    genome = load_genome(args.genome)
    if args.chromosome not in genome:
        raise ValueError(f"Chromosome '{args.chromosome}' not found in the genome FASTA file.")
    chrom_seq = genome[args.chromosome]

    # Process the VCF file to get nucleotide counts
    total_callable, count_A, count_T, count_C, count_G, count_CG, count_CPD = process_vcf(args.vcf_file, chrom_seq)

    # Write the summary to the specified output file
    header = (
        "species\tind\tchromosome\ttotal_callable_site\ttotal_A_site\ttotal_T_site\t"
        "total_C_site\ttotal_G_site\ttotal_CG_dinucleotide\ttotal_CPD_site\n"
    )
    summary_line = (
        f"{args.species}\t{args.ind}\t{args.chromosome}\t{total_callable}\t{count_A}\t{count_T}"
        f"\t{count_C}\t{count_G}\t{count_CG}\t{count_CPD}\n"
    )

    with open(args.output, "w") as out_f:
        out_f.write(header)
        out_f.write(summary_line)

if __name__ == "__main__":
    main()

