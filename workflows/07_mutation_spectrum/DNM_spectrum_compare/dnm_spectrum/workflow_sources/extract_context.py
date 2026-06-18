#!/usr/bin/env python3
"""
extract_context.py
For each DNM in a species, extract the trinucleotide context from the
reference genome and assign the SBS96 category with strand collapsing
(pyrimidine reference convention).
"""
import argparse
import subprocess
import sys
import os

COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def reverse_complement(seq):
    return ''.join(COMPLEMENT[b] for b in reversed(seq.upper()))


def get_trinucleotides(genome_fa, regions, samtools_path='samtools'):
    """
    Batch-extract trinucleotide sequences from the genome using samtools faidx.
    regions: list of (chrom, pos) tuples with 1-based pos.
    Returns dict mapping (chrom, pos) -> 3-letter uppercase string.
    """
    if not regions:
        return {}
    # Write regions file
    regions_file = genome_fa + '.tmp_regions.txt'
    with open(regions_file, 'w') as f:
        for chrom, pos in regions:
            start = pos - 1
            end = pos + 1
            f.write(f'{chrom}:{start}-{end}\n')
    # Run samtools faidx
    cmd = [samtools_path, 'faidx', genome_fa, '-r', regions_file]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"samtools faidx error: {result.stderr}", file=sys.stderr)
        sys.exit(1)
    os.remove(regions_file)
    # Parse FASTA output
    trinucs = {}
    current_key = None
    current_seq = []
    for line in result.stdout.strip().split('\n'):
        if line.startswith('>'):
            if current_key is not None:
                trinucs[current_key] = ''.join(current_seq).upper()
            header = line[1:].split()[0]  # e.g. afr_1:230325934-230325936
            parts = header.split(':')
            chrom = parts[0]
            coords = parts[1].split('-')
            pos = int(coords[0]) + 1  # convert back to original pos
            current_key = (chrom, pos)
            current_seq = []
        else:
            current_seq.append(line.strip())
    if current_key is not None:
        trinucs[current_key] = ''.join(current_seq).upper()
    return trinucs


def assign_sbs96(ref, alt, trinuc):
    """
    Given ref allele, alt allele, and 3-base trinucleotide context,
    return (collapsed_trinuc, mutation_type, sbs96_category).
    Collapses to pyrimidine reference convention.
    """
    ref = ref.upper()
    alt = alt.upper()
    trinuc = trinuc.upper()
    if ref in ('C', 'T'):
        mut_type = f'{ref}>{alt}'
        category = f'{trinuc[0]}[{mut_type}]{trinuc[2]}'
        return trinuc, mut_type, category
    else:
        # Purine reference: reverse complement
        rc_trinuc = reverse_complement(trinuc)
        rc_ref = COMPLEMENT[ref]
        rc_alt = COMPLEMENT[alt]
        mut_type = f'{rc_ref}>{rc_alt}'
        category = f'{rc_trinuc[0]}[{mut_type}]{rc_trinuc[2]}'
        return rc_trinuc, mut_type, category


def main():
    parser = argparse.ArgumentParser(description='Extract trinucleotide context and assign SBS96 categories')
    parser.add_argument('--dnm_file', required=True, help='Input DNM TSV file')
    parser.add_argument('--genome', required=True, help='Reference genome FASTA (indexed)')
    parser.add_argument('--output', required=True, help='Output annotated TSV')
    parser.add_argument('--x_chromosomes', nargs='*', default=[], help='X chromosome names')
    parser.add_argument('--exclude_trios', nargs='*', default=[], help='Offspring IDs to exclude')
    parser.add_argument('--samtools', default='samtools', help='Path to samtools')
    args = parser.parse_args()

    # Read DNM file
    rows = []
    with open(args.dnm_file) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            row = dict(zip(header, fields))
            rows.append(row)

    # Filter excluded trios
    if args.exclude_trios:
        exclude_set = set(args.exclude_trios)
        rows = [r for r in rows if r['child'] not in exclude_set]

    x_set = set(args.x_chromosomes)

    # Collect regions for batch extraction
    regions = [(r['chrom'], int(r['pos'])) for r in rows]

    # Extract trinucleotides
    trinucs = get_trinucleotides(args.genome, regions, samtools_path=args.samtools)

    # Write output
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, 'w') as out:
        out.write('\t'.join([
            'chrom', 'pos', 'ref', 'alt', 'child', 'father', 'mother',
            'trinuc_context', 'mutation_type', 'sbs96_category', 'chrom_group'
        ]) + '\n')
        for row in rows:
            chrom = row['chrom']
            pos = int(row['pos'])
            ref = row['ref']
            alt = row['alt']
            key = (chrom, pos)
            trinuc = trinucs.get(key, 'NNN')
            if 'N' in trinuc or len(trinuc) != 3:
                print(f"Warning: skipping {chrom}:{pos} — trinucleotide={trinuc}", file=sys.stderr)
                continue
            # Verify the middle base matches ref
            if trinuc[1].upper() != ref.upper():
                print(f"Warning: ref mismatch at {chrom}:{pos} — "
                      f"genome={trinuc[1]}, dnm_ref={ref}", file=sys.stderr)
                continue
            collapsed_trinuc, mut_type, sbs96 = assign_sbs96(ref, alt, trinuc)
            chrom_group = 'x_chromosome' if chrom in x_set else 'autosome'
            out.write('\t'.join([
                chrom, str(pos), ref, alt,
                row['child'], row['father'], row['mother'],
                collapsed_trinuc, mut_type, sbs96, chrom_group
            ]) + '\n')

    print(f"Annotated {len(rows)} DNMs -> {args.output}")


if __name__ == '__main__':
    main()
