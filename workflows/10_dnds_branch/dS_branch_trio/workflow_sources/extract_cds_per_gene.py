#!/usr/bin/env python3
"""
extract_cds_per_gene.py - Extract per-gene CDS alignments from trio genomes

For each gene in the GFF3 annotation:
  1. Parse CDS regions (coordinates on BIC reference)
  2. Extract CDS sequence from BIC reference, sp1 consensus, sp2 consensus
  3. If strand == '-': reverse complement all 3
  4. QC checks:
     a. CDS length divisible by 3
     b. No premature stop codons in BIC reference
     c. Each species has <= max_n_frac fraction of Ns (coverage check)
  5. Write per-gene FASTA with 3 sequences (BIC, sp1, sp2)
  6. Write passing_genes.txt and gene_stats.tsv

Usage:
    python3 extract_cds_per_gene.py \
        --ref BIC.fasta \
        --sp1_name SAR --sp1_fa SAR_consensus.fa \
        --sp2_name PAC --sp2_fa PAC_consensus.fa \
        --gff3 BIC.gff3 \
        --output_dir steps/06_cds_per_gene/ \
        --min_cov_frac 0.80
"""

import argparse
import os
import sys
from collections import defaultdict

import pysam


# Standard codon table
STOP_CODONS = {"TAA", "TAG", "TGA"}

COMPLEMENT = str.maketrans("ATCGatcgRYSWKMBDHVNryswkmbdhvn",
                           "TAGCtagcYRSWMKVHDBNyrswmkvhdbn")


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def parse_gff3(gff3_path):
    """
    Parse GFF3 and build per-gene CDS region dict.

    Strategy: collect all mRNA->CDS relationships.
    For each gene, pick the mRNA with the longest total CDS.

    Returns:
        dict: {gene_id: {'chrom': str, 'strand': str,
                          'cds_regions': [(start, end), ...], 'mrna_id': str}}
    """
    # First pass: collect gene->mRNA and mRNA->CDS mappings
    gene_to_mrnas = defaultdict(list)
    mrna_to_cds = defaultdict(list)
    mrna_info = {}  # mrna_id -> {chrom, strand}
    gene_info = {}  # gene_id -> {chrom, strand}

    with open(gff3_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue

            chrom = parts[0]
            feature_type = parts[2]
            start = int(parts[3])  # 1-based inclusive in GFF3
            end = int(parts[4])    # 1-based inclusive in GFF3
            strand = parts[6]
            attributes = parts[8]

            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, val = attr.split('=', 1)
                    attr_dict[key] = val

            if feature_type == 'gene':
                gid = attr_dict.get('ID', '')
                if gid:
                    gene_info[gid] = {'chrom': chrom, 'strand': strand}

            elif feature_type == 'mRNA':
                mrna_id = attr_dict.get('ID', '')
                parent = attr_dict.get('Parent', '')
                if mrna_id and parent:
                    gene_to_mrnas[parent].append(mrna_id)
                    mrna_info[mrna_id] = {'chrom': chrom, 'strand': strand}

            elif feature_type == 'CDS':
                parent = attr_dict.get('Parent', '')
                if parent:
                    # Parent can be comma-separated
                    for p in parent.split(','):
                        mrna_to_cds[p].append((start, end))

    # For each gene, pick the mRNA with the longest total CDS
    genes = {}
    for gene_id, mrna_list in gene_to_mrnas.items():
        best_mrna = None
        best_length = 0
        best_cds = []

        for mrna_id in mrna_list:
            cds_regions = mrna_to_cds.get(mrna_id, [])
            if not cds_regions:
                continue
            total_len = sum(e - s + 1 for s, e in cds_regions)
            if total_len > best_length:
                best_length = total_len
                best_mrna = mrna_id
                best_cds = cds_regions

        if best_cds and best_mrna:
            info = mrna_info.get(best_mrna, gene_info.get(gene_id, {}))
            genes[gene_id] = {
                'chrom': info.get('chrom', ''),
                'strand': info.get('strand', '+'),
                'cds_regions': sorted(best_cds, key=lambda x: x[0]),
                'mrna_id': best_mrna
            }

    return genes


def extract_cds_sequence(fasta_handle, chrom, cds_regions, strand):
    """
    Extract concatenated CDS sequence from a pysam FastaFile.

    Args:
        fasta_handle: pysam.FastaFile
        chrom: chromosome name
        cds_regions: list of (start, end) tuples (1-based inclusive)
        strand: '+' or '-'

    Returns:
        str: CDS sequence (reverse complemented if strand == '-')
    """
    # Check if chrom exists in this fasta
    if chrom not in fasta_handle.references:
        return None

    seq_parts = []
    for start, end in cds_regions:
        # pysam fetch is 0-based half-open
        region_seq = fasta_handle.fetch(chrom, start - 1, end)
        seq_parts.append(region_seq.upper())

    full_seq = ''.join(seq_parts)

    if strand == '-':
        full_seq = reverse_complement(full_seq)

    return full_seq


def clean_sequence(seq):
    """Replace any non-ATCG character with N."""
    cleaned = []
    for c in seq:
        if c in 'ATCG':
            cleaned.append(c)
        else:
            cleaned.append('N')
    return ''.join(cleaned)


def n_fraction(seq):
    """Calculate fraction of N characters in sequence."""
    if len(seq) == 0:
        return 1.0
    return seq.count('N') / len(seq)


def has_premature_stop(seq):
    """
    Check if sequence has premature stop codons (excluding the last codon).
    Assumes sequence length is divisible by 3.
    """
    if len(seq) < 6:
        return False

    # Check all codons except the last one
    for i in range(0, len(seq) - 3, 3):
        codon = seq[i:i+3]
        if codon in STOP_CODONS:
            return True
    return False


def main():
    parser = argparse.ArgumentParser(
        description='Extract per-gene CDS alignments from trio genomes')
    parser.add_argument('--ref', required=True,
                        help='Reference genome FASTA')
    parser.add_argument('--ref_name', default='BIC',
                        help='Reference species name to write in output FASTA headers')
    parser.add_argument('--sp1_name', required=True,
                        help='Name of first mapping species (e.g. SAR)')
    parser.add_argument('--sp1_fa', required=True,
                        help='Consensus FASTA for first mapping species')
    parser.add_argument('--sp2_name', required=True,
                        help='Name of second mapping species (e.g. PAC)')
    parser.add_argument('--sp2_fa', required=True,
                        help='Consensus FASTA for second mapping species')
    parser.add_argument('--gff3', required=True,
                        help='GFF3 annotation file')
    parser.add_argument('--output_dir', required=True,
                        help='Output directory for per-gene FASTA files')
    parser.add_argument('--min_cov_frac', type=float, default=0.80,
                        help='Min fraction of non-N bases per species (default: 0.80)')
    args = parser.parse_args()

    sp1_name = args.sp1_name
    sp2_name = args.sp2_name
    ref_name = args.ref_name

    max_n_frac = 1.0 - args.min_cov_frac

    os.makedirs(args.output_dir, exist_ok=True)
    fasta_dir = os.path.join(args.output_dir, 'per_gene_fasta')
    os.makedirs(fasta_dir, exist_ok=True)

    # Parse GFF3
    print("Parsing GFF3 annotation...")
    genes = parse_gff3(args.gff3)
    print(f"  Found {len(genes)} genes with CDS regions")
    print(f"  Min coverage fraction: {args.min_cov_frac} (max N fraction: {max_n_frac})")

    # Open FASTA files
    ref_fa = pysam.FastaFile(args.ref)
    sp1_fa = pysam.FastaFile(args.sp1_fa)
    sp2_fa = pysam.FastaFile(args.sp2_fa)

    # Track stats
    stats = []
    passing_genes = []

    print("Extracting CDS per gene...")
    for gene_id, gene_data in sorted(genes.items()):
        chrom = gene_data['chrom']
        strand = gene_data['strand']
        cds_regions = gene_data['cds_regions']
        n_exons = len(cds_regions)

        # Extract CDS from all 3 genomes
        ref_seq = extract_cds_sequence(ref_fa, chrom, cds_regions, strand)
        sp1_seq = extract_cds_sequence(sp1_fa, chrom, cds_regions, strand)
        sp2_seq = extract_cds_sequence(sp2_fa, chrom, cds_regions, strand)

        # Check extraction success
        if ref_seq is None or sp1_seq is None or sp2_seq is None:
            stats.append((gene_id, 0, n_exons, 'FAIL', 'chrom_not_found',
                          0.0, 0.0, 0.0))
            continue

        cds_length = len(ref_seq)

        # Clean sequences (replace non-ATCG with N)
        ref_clean = clean_sequence(ref_seq)
        sp1_clean = clean_sequence(sp1_seq)
        sp2_clean = clean_sequence(sp2_seq)

        # Calculate N fractions per species
        ref_n_frac = n_fraction(ref_clean)
        sp1_n_frac = n_fraction(sp1_clean)
        sp2_n_frac = n_fraction(sp2_clean)

        # QC checks
        reason = 'PASS'
        passed = True

        if cds_length % 3 != 0:
            reason = 'length_not_divisible_by_3'
            passed = False
        elif has_premature_stop(ref_clean):
            reason = f'premature_stop_in_{ref_name}'
            passed = False
        elif ref_n_frac > max_n_frac:
            reason = f'{ref_name}_low_coverage({ref_n_frac:.2f})'
            passed = False
        elif sp1_n_frac > max_n_frac:
            reason = f'{sp1_name}_low_coverage({sp1_n_frac:.2f})'
            passed = False
        elif sp2_n_frac > max_n_frac:
            reason = f'{sp2_name}_low_coverage({sp2_n_frac:.2f})'
            passed = False

        # Write per-gene FASTA regardless of pass/fail
        gene_fasta = os.path.join(fasta_dir, f"{gene_id}.fa")
        with open(gene_fasta, 'w') as fout:
            fout.write(f">{ref_name}\n{ref_clean}\n")
            fout.write(f">{sp1_name}\n{sp1_clean}\n")
            fout.write(f">{sp2_name}\n{sp2_clean}\n")

        stats.append((gene_id, cds_length, n_exons,
                       'PASS' if passed else 'FAIL', reason,
                       ref_n_frac, sp1_n_frac, sp2_n_frac))

        if passed:
            passing_genes.append(gene_id)

    ref_fa.close()
    sp1_fa.close()
    sp2_fa.close()

    # Write passing genes list
    passing_file = os.path.join(args.output_dir, 'passing_genes.txt')
    with open(passing_file, 'w') as f:
        for gid in passing_genes:
            f.write(f"{gid}\n")

    # Write stats TSV
    stats_file = os.path.join(args.output_dir, 'gene_stats.tsv')
    with open(stats_file, 'w') as f:
        f.write(f"gene_id\tcds_length\tn_exons\tstatus\treason\t"
                f"{ref_name}_n_frac\t{sp1_name}_n_frac\t{sp2_name}_n_frac\n")
        for gene_id, cds_len, n_exons, status, reason, ref_nf, sp1_nf, sp2_nf in stats:
            f.write(f"{gene_id}\t{cds_len}\t{n_exons}\t{status}\t{reason}\t"
                    f"{ref_nf:.4f}\t{sp1_nf:.4f}\t{sp2_nf:.4f}\n")

    # Summary
    n_pass = len(passing_genes)
    n_total = len(stats)
    n_fail = n_total - n_pass
    print(f"\nSummary:")
    print(f"  Total genes processed: {n_total}")
    print(f"  Passing: {n_pass}")
    print(f"  Failing: {n_fail}")
    print(f"  Min coverage fraction: {args.min_cov_frac}")
    print(f"  Output: {args.output_dir}")
    print(f"  Passing genes list: {passing_file}")
    print(f"  Stats: {stats_file}")


if __name__ == '__main__':
    main()
