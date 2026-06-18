#!/usr/bin/env python3
"""
compute_callable_bed.py

For one offspring, parse all autosomal *_callable_sites.vcf files in
chrom_split/ and merge consecutive callable positions into BED segments.
Small intra-segment gaps (≤ merge_gap bp) are bridged so that visualisation
is not broken by single uncallable sites.

Output columns (TSV with header):
    offspring  chrom  start  end
    (start is 0-based, end is exclusive — standard BED half-open coordinates)
"""
import argparse
import os
import glob


def load_fai(fai_path):
    """Return {chrom: length} dict."""
    chrom_lengths = {}
    with open(fai_path) as fh:
        for line in fh:
            parts = line.strip().split('\t')
            chrom_lengths[parts[0]] = int(parts[1])
    return chrom_lengths


def extract_chrom_from_filename(vcf_filename, offspring):
    """
    VCF filename pattern: {offspring}_{chrom}_callable_sites.vcf
    Strip the offspring prefix and '_callable_sites' suffix to get chrom.
    """
    basename = os.path.basename(vcf_filename)
    after_offspring = basename[len(offspring) + 1:]
    chrom = after_offspring.replace('_callable_sites.vcf', '')
    return chrom


def positions_to_segments(positions, merge_gap):
    """
    Convert a sorted list of 1-based positions into 0-based half-open BED
    segments, bridging gaps of ≤ merge_gap bp.
    """
    if not positions:
        return []
    segments = []
    seg_start = positions[0] - 1   # 0-based
    seg_end   = positions[0]        # exclusive end
    for pos in positions[1:]:
        if pos - seg_end <= merge_gap:
            seg_end = pos           # extend current segment
        else:
            segments.append((seg_start, seg_end))
            seg_start = pos - 1
            seg_end   = pos
    segments.append((seg_start, seg_end))
    return segments


def main():
    parser = argparse.ArgumentParser(
        description='Convert per-offspring callable VCFs to merged BED segments.'
    )
    parser.add_argument('--vcf_dir',    required=True,
                        help='Directory containing *_callable_sites.vcf files')
    parser.add_argument('--offspring',  required=True,
                        help='Offspring identifier (e.g. AFR_family1_S1_offspring)')
    parser.add_argument('--genome_fai', required=True,
                        help='Path to genome .fai file')
    parser.add_argument('--x_chroms',  nargs='+', default=[],
                        help='Sex chromosome names to exclude')
    parser.add_argument('--merge_gap', type=int, default=100,
                        help='Bridge intra-segment gaps ≤ this many bp (default: 100)')
    parser.add_argument('--output',    required=True,
                        help='Output TSV path')
    args = parser.parse_args()

    chrom_lengths = load_fai(args.genome_fai)
    x_chroms_set  = set(args.x_chroms)

    vcf_files = glob.glob(os.path.join(args.vcf_dir, '*_callable_sites.vcf'))

    with open(args.output, 'w') as out:
        out.write('offspring\tchrom\tstart\tend\n')

        for vcf_path in sorted(vcf_files):
            chrom = extract_chrom_from_filename(vcf_path, args.offspring)
            if chrom in x_chroms_set:
                continue
            if chrom not in chrom_lengths:
                continue

            # Collect all callable positions for this chromosome
            positions = []
            with open(vcf_path) as fh:
                for line in fh:
                    if line.startswith('#'):
                        continue
                    fields = line.split('\t', 2)
                    if len(fields) < 2:
                        continue
                    positions.append(int(fields[1]))

            positions.sort()
            segments = positions_to_segments(positions, args.merge_gap)

            for start, end in segments:
                out.write(f"{args.offspring}\t{chrom}\t{start}\t{end}\n")

    print(f"Done: {args.output}")


if __name__ == '__main__':
    main()
