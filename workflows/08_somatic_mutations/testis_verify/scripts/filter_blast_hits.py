#!/usr/bin/env python3
"""
Filter BLAST results and create homolog mapping.

Parses tblastn output, filters by e-value, identity, and coverage,
and identifies best hits per human gene.
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path


def parse_gene_name_from_query(query: str) -> str:
    """Extract gene name from UniProt query ID.

    Input: sp|A0A087X1C5|CP2D7_HUMAN
    Output: CYP2D7 (from GN= field) or CP2D7 (from entry name)
    """
    # The query in BLAST output is typically just the ID part
    # We need to extract gene name from entry name (e.g., CP2D7_HUMAN -> CP2D7)
    parts = query.split('|')
    if len(parts) >= 3:
        entry_name = parts[2]
        # Remove _HUMAN suffix
        if '_HUMAN' in entry_name:
            return entry_name.replace('_HUMAN', '')
    return query


def extract_mim_gene_id(subject: str) -> str:
    """Extract MIM gene ID from subject hit.

    For CDS extracted by gffread, subject format is typically:
    gene_id or transcript_id
    """
    # Remove transcript suffix if present (e.g., MIM.gene.t1 -> MIM.gene)
    if '.t' in subject:
        # Keep just the gene part
        gene_id = re.sub(r'\.t\d+$', '', subject)
        return gene_id
    return subject


def parse_blast_results(blast_file: Path, min_evalue: float,
                        min_identity: float, min_coverage: float) -> dict:
    """Parse BLAST tabular output and filter hits.

    Expected format: outfmt 6 with columns:
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen

    Returns dict mapping human_gene -> list of (mim_gene, identity, evalue, bitscore, coverage)
    """
    hits_by_gene = defaultdict(list)

    with open(blast_file) as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 14:
                continue

            qseqid = fields[0]
            sseqid = fields[1]
            pident = float(fields[2])
            length = int(fields[3])
            qstart = int(fields[6])
            qend = int(fields[7])
            evalue = float(fields[10])
            bitscore = float(fields[11])
            qlen = int(fields[12])

            # Calculate query coverage
            alignment_len = qend - qstart + 1
            coverage = (alignment_len / qlen) * 100

            # Apply filters
            if evalue > min_evalue:
                continue
            if pident < min_identity:
                continue
            if coverage < min_coverage:
                continue

            # Extract gene names
            human_gene = parse_gene_name_from_query(qseqid)
            mim_gene = extract_mim_gene_id(sseqid)

            hits_by_gene[human_gene].append({
                'mim_gene': mim_gene,
                'identity': pident,
                'evalue': evalue,
                'bitscore': bitscore,
                'coverage': coverage,
                'qseqid': qseqid
            })

    return hits_by_gene


def get_best_hits(hits_by_gene: dict) -> dict:
    """Select best hit (highest bitscore) per human gene."""
    best_hits = {}

    for human_gene, hits in hits_by_gene.items():
        # Sort by bitscore descending
        sorted_hits = sorted(hits, key=lambda x: x['bitscore'], reverse=True)
        best_hit = sorted_hits[0]
        best_hit['hit_count'] = len(hits)
        best_hits[human_gene] = best_hit

    return best_hits


def write_outputs(best_hits: dict, output_dir: Path):
    """Write output files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Homolog mapping (detailed)
    mapping_file = output_dir / 'homolog_mapping.tsv'
    with open(mapping_file, 'w') as f:
        f.write("human_gene\tmim_gene\tidentity\tevalue\tbitscore\tcoverage\thit_count\n")
        for human_gene in sorted(best_hits.keys()):
            hit = best_hits[human_gene]
            f.write(f"{human_gene}\t{hit['mim_gene']}\t{hit['identity']:.1f}\t")
            f.write(f"{hit['evalue']:.2e}\t{hit['bitscore']:.1f}\t")
            f.write(f"{hit['coverage']:.1f}\t{hit['hit_count']}\n")

    # MIM gene list (unique best hits)
    mim_genes = set(hit['mim_gene'] for hit in best_hits.values())
    gene_list_file = output_dir / 'mim_testis_homologs.txt'
    with open(gene_list_file, 'w') as f:
        for gene in sorted(mim_genes):
            f.write(f"{gene}\n")

    # Summary statistics
    summary_file = output_dir / 'blast_summary.txt'

    hit_counts = [hit['hit_count'] for hit in best_hits.values()]
    identities = [hit['identity'] for hit in best_hits.values()]

    with open(summary_file, 'w') as f:
        f.write("BLAST Homolog Search Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Human genes with homologs: {len(best_hits)}\n")
        f.write(f"Unique MIM genes identified: {len(mim_genes)}\n\n")

        f.write("Hit count distribution:\n")
        single_hit = sum(1 for c in hit_counts if c == 1)
        multi_hit = sum(1 for c in hit_counts if c > 1)
        f.write(f"  Genes with single hit: {single_hit}\n")
        f.write(f"  Genes with multiple hits: {multi_hit}\n")
        if hit_counts:
            f.write(f"  Max hits per gene: {max(hit_counts)}\n")
            f.write(f"  Mean hits per gene: {sum(hit_counts)/len(hit_counts):.1f}\n\n")

        f.write("Identity distribution:\n")
        if identities:
            f.write(f"  Min identity: {min(identities):.1f}%\n")
            f.write(f"  Max identity: {max(identities):.1f}%\n")
            f.write(f"  Mean identity: {sum(identities)/len(identities):.1f}%\n")

    return mapping_file, gene_list_file, summary_file


def main():
    parser = argparse.ArgumentParser(
        description='Filter BLAST results and create homolog mapping'
    )
    parser.add_argument('--blast', required=True, type=Path,
                        help='BLAST tabular output file')
    parser.add_argument('--output-dir', required=True, type=Path,
                        help='Output directory')
    parser.add_argument('--min-evalue', type=float, default=1e-5,
                        help='Maximum e-value threshold (default: 1e-5)')
    parser.add_argument('--min-identity', type=float, default=30.0,
                        help='Minimum identity %% (default: 30)')
    parser.add_argument('--min-coverage', type=float, default=50.0,
                        help='Minimum query coverage %% (default: 50)')

    args = parser.parse_args()

    print(f"Parsing BLAST results from {args.blast}...")
    print(f"  E-value threshold: {args.min_evalue}")
    print(f"  Min identity: {args.min_identity}%")
    print(f"  Min coverage: {args.min_coverage}%")

    # Parse and filter
    hits_by_gene = parse_blast_results(
        args.blast, args.min_evalue, args.min_identity, args.min_coverage
    )

    print(f"\nFound hits for {len(hits_by_gene)} human genes")

    # Get best hits
    best_hits = get_best_hits(hits_by_gene)

    # Write outputs
    mapping_file, gene_list_file, summary_file = write_outputs(
        best_hits, args.output_dir
    )

    print(f"\nOutputs written:")
    print(f"  Mapping: {mapping_file}")
    print(f"  Gene list: {gene_list_file}")
    print(f"  Summary: {summary_file}")

    # Print summary
    mim_genes = set(hit['mim_gene'] for hit in best_hits.values())
    print(f"\nSummary:")
    print(f"  Human genes with homologs: {len(best_hits)}")
    print(f"  Unique MIM homologs: {len(mim_genes)}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
