#!/usr/bin/env python3
"""
callable_vcf_trinuc.py
For one (species, offspring, chromosome) callable VCF, count the trinucleotide
context of every callable site in pyrimidine-canonical form.

Optimized streaming reader: keeps the chromosome FASTA in memory once and
slices a 3bp window per row. Skips chromosome edges and any window containing
N. Reverse-complements when REF is A/G so the middle base is always C/T.

Input VCF rows look like:
    afr_10  455  .  C  .  .  PASS  .  GT:AD:DP:RGQ  ...
"""
import argparse
import os
import sys

COMPLEMENT = str.maketrans('ACGT', 'TGCA')

PYRIMIDINE_CANONICAL_TRINUCS = []
for first in 'ACGT':
    for mid in 'CT':
        for last in 'ACGT':
            PYRIMIDINE_CANONICAL_TRINUCS.append(first + mid + last)
# 32 trinucs, sorted ACA, ACC, ..., TTT


def reverse_complement(s):
    return s.translate(COMPLEMENT)[::-1]


def load_chromosome(genome_fa, target_chrom):
    """Read a single chromosome from a FASTA into one uppercase string.
    Uses the .fai (samtools faidx) for an O(1) seek+read."""
    fai = genome_fa + '.fai'
    if not os.path.exists(fai):
        sys.exit(f'FASTA index missing: {fai}')

    target_offset = None
    target_length = None
    target_linebases = None
    target_linewidth = None
    with open(fai) as f:
        for line in f:
            name, length, offset, linebases, linewidth = line.rstrip('\n').split('\t')
            if name == target_chrom:
                target_offset = int(offset)
                target_length = int(length)
                target_linebases = int(linebases)
                target_linewidth = int(linewidth)
                break
    if target_offset is None:
        sys.exit(f'Chromosome {target_chrom} not found in {fai}')

    n_full_lines = target_length // target_linebases
    leftover = target_length - n_full_lines * target_linebases
    bytes_to_read = n_full_lines * target_linewidth + leftover
    with open(genome_fa, 'rb') as f:
        f.seek(target_offset)
        raw = f.read(bytes_to_read)
    seq = raw.decode('ascii').replace('\n', '').replace('\r', '').upper()
    if len(seq) != target_length:
        sys.exit(
            f'Chromosome {target_chrom}: read {len(seq)} bp but expected {target_length}'
        )
    return seq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--genome', required=True)
    ap.add_argument('--chrom', required=True,
                    help='Chromosome name (must match the VCF and FASTA).')
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    seq = load_chromosome(args.genome, args.chrom)
    seq_len = len(seq)

    counts = {t: 0 for t in PYRIMIDINE_CANONICAL_TRINUCS}
    n_skipped_n = 0
    n_skipped_edge = 0
    n_skipped_mismatch = 0
    n_processed = 0

    with open(args.vcf) as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            tab1 = line.find('\t')
            if tab1 == -1:
                continue
            tab2 = line.find('\t', tab1 + 1)
            if tab2 == -1:
                continue
            tab3 = line.find('\t', tab2 + 1)
            if tab3 == -1:
                continue
            tab4 = line.find('\t', tab3 + 1)
            if tab4 == -1:
                continue
            try:
                pos = int(line[tab1 + 1:tab2])
            except ValueError:
                continue
            ref = line[tab3 + 1:tab4]
            if len(ref) != 1:
                continue
            ref = ref.upper()
            if ref not in 'ACGT':
                continue
            if pos < 2 or pos > seq_len - 1:
                n_skipped_edge += 1
                continue
            tri = seq[pos - 2:pos + 1]
            if 'N' in tri:
                n_skipped_n += 1
                continue
            mid = tri[1]
            if mid != ref:
                n_skipped_mismatch += 1
                continue
            if mid in 'CT':
                canonical = tri
            else:
                canonical = reverse_complement(tri)
            counts[canonical] = counts.get(canonical, 0) + 1
            n_processed += 1

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, 'w') as out:
        out.write(
            f'# chrom={args.chrom} processed={n_processed} '
            f'skipped_edge={n_skipped_edge} skipped_n={n_skipped_n} '
            f'skipped_ref_mismatch={n_skipped_mismatch}\n'
        )
        out.write('trinuc_canonical\tcount\n')
        for t in PYRIMIDINE_CANONICAL_TRINUCS:
            out.write(f'{t}\t{counts[t]}\n')

    total = sum(counts.values())
    print(
        f'{args.chrom}: {total} callable trinucs '
        f'(processed={n_processed}, skipped_edge={n_skipped_edge}, '
        f'skipped_n={n_skipped_n}, skipped_mismatch={n_skipped_mismatch}) '
        f'-> {args.output}'
    )


if __name__ == '__main__':
    main()
