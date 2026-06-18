#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def load_reference_genome(fasta_path):
    """
    Load the reference genome from a FASTA file.
    Returns a dictionary mapping chromosome names to sequences.
    """
    genome = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        genome[record.id] = record.seq
    return genome

def classify_mutation(n_before, ref, n_after, alt):
    """
    Classify the mutation into one of the following types and flags:
      - mut_type: a string that can contain one or more labels.
      - if_CpG: True if the mutation occurs in a CpG context.
      - if_CPD: True if the mutation occurs in a dipyrimidine (CPD) context.
    
    The seven possible labels include:
      "C>A", "C>T", "C>T,CPD", "C>G", "T>A", "T>C", "T>G", "CpG>TpG"
    
    For C>T changes (after converting to a pyrimidine context), the function:
      1. Checks if either neighbor is a G. If so, it flags if_CpG = True.
      2. Checks if either neighbor is a pyrimidine (C or T). If so, it flags if_CPD = True.
      3. Returns a mut_type string that always includes the base "C>T"
         and, if applicable, appends ",CpG>TpG" and/or ",C>T,CPD".
    
    For all other mutations, the base classification is returned and both flags are False.
    """
    # Build the trinucleotide context.
    trinuc = (n_before + ref + n_after).upper()
    
    # If the reference base is not a pyrimidine, convert the trinuc and alt allele
    # to the reverse complement to represent the mutation in a pyrimidine context.
    if ref.upper() not in ["C", "T"]:
        trinuc = str(Seq(trinuc).reverse_complement())
        alt = str(Seq(alt.upper()).complement())
    
    norm_ref = trinuc[1]  # normalized reference (should be C or T)
    norm_alt = alt.upper()
    left_base = trinuc[0]
    right_base = trinuc[2]
    
    # Default flags.
    if_CpG = False
    if_CPD = False
    
    # Classify based on normalized context.
    if norm_ref == "C":
        if norm_alt == "T":
            # For C>T, check both CpG and CPD conditions.
            if_CpG = (left_base == "G" or right_base == "G")
            if_CPD = (left_base in ("C", "T") or right_base in ("C", "T"))
            # Start with base classification.
            mut_type = "C>T"
            if if_CpG:
                mut_type += ",CpG>TpG"
            if if_CPD:
                mut_type += ",C>T_CPD"
            return mut_type, if_CpG, if_CPD
        elif norm_alt == "A":
            return "C>A", False, False
        elif norm_alt == "G":
            return "C>G", False, False
    elif norm_ref == "T":
        if norm_alt == "A":
            return "T>A", False, False
        elif norm_alt == "C":
            return "T>C", False, False
        elif norm_alt == "G":
            return "T>G", False, False
    return "Other", False, False

def process_dnms(genome, dnm_tsv, output_tsv):
    """
    Process the de novo mutations TSV file:
      - For each mutation, extract the flanking nucleotides and the reference nucleotide.
      - Classify the mutation to get mut_type, if_CpG, and if_CPD.
      - Write an output file with the original seven columns plus:
          n_before, nuc_ref, nuc_after, mut_type, if_CpG, if_CPD.
    
    Assumes the input TSV file does not have a header.
    """
    header = ("chrom\tpos\tref\talt\tchild\tfather\tmother\t"
              "n_before\tnuc_ref\tnuc_after\tmut_type\tif_CpG\tif_CPD")
    
    with open(dnm_tsv, "r") as infile, open(output_tsv, "w") as outfile:
        outfile.write(header + "\n")
        
        for line in infile:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            
            chrom, pos, ref, alt, child, father, mother = parts[:7]
            pos = int(pos)  # assume 1-based coordinate
            
            # Retrieve sequence for the given chromosome.
            seq = genome.get(chrom)
            if seq is None:
                print(f"Warning: Chromosome {chrom} not found in the reference genome!")
                continue

            # Extract flanking nucleotides (adjusting for 0-indexing):
            # n_before: nucleotide immediately before mutation (index pos-2)
            # nuc_ref: nucleotide at mutation (index pos-1)
            # nuc_after: nucleotide immediately after mutation (index pos)
            n_before = str(seq[pos - 2]) if pos - 2 >= 0 else "N"
            nuc_ref = str(seq[pos - 1]) if pos - 1 < len(seq) else "N"
            nuc_after = str(seq[pos]) if pos < len(seq) else "N"
            
            # Sanity check.
            if nuc_ref.upper() != ref.upper():
                print(f"Warning: Mismatch at {chrom}:{pos}. Genome nucleotide is {nuc_ref} but TSV ref is {ref}.")

            # Classify the mutation.
            mut_type, flag_CpG, flag_CPD = classify_mutation(n_before, ref, nuc_after, alt)
            
            # Write output line with additional columns.
            outfile.write(
                f"{chrom}\t{pos}\t{ref}\t{alt}\t{child}\t{father}\t{mother}\t"
                f"{n_before}\t{nuc_ref}\t{nuc_after}\t{mut_type}\t"
                f"{str(flag_CpG).lower()}\t{str(flag_CPD).lower()}\n"
            )

def main():
    parser = argparse.ArgumentParser(
        description=("Classify de novo mutations by extracting surrounding nucleotides "
                     "and categorizing mutation types. For C>T changes, both CpG and CPD flags "
                     "are computed, and the mut_type column includes both labels if applicable.")
    )
    parser.add_argument("reference", help="Reference genome FASTA file")
    parser.add_argument("dnm_tsv", help=("Input TSV file with de novo mutation data "
                                          "(chrom, pos, ref, alt, child, father, mother) without header"))
    parser.add_argument("output_tsv", help=("Output TSV file with additional columns: "
                                             "n_before, nuc_ref, nuc_after, mut_type, if_CpG, if_CPD"))
    args = parser.parse_args()
    
    genome = load_reference_genome(args.reference)
    process_dnms(genome, args.dnm_tsv, args.output_tsv)

if __name__ == "__main__":
    main()

