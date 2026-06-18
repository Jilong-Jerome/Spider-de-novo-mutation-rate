#!/usr/bin/env python3
"""
Single-chromosome coverage statistics for one individual from a `samtools depth` stream.

Reads `chrom \t pos \t depth` lines on stdin (output of `samtools depth -d 0 -r <chrom> <bam>`,
which by default already excludes unmapped/secondary/qcfail/duplicate reads). Only lines
matching --chrom are counted. Builds a depth histogram (memory-light) and writes a single
row with:

  mean_all_depth = sum_depth / chrom_length              (requirement 1)
  mean_covered   = sum_depth / covered_bases             (requirement 2)
  median_covered = median depth over covered sites only  (requirement 2)

The chromosome length is read from a `chrom \t length` TSV (the contigs step output).
If the chromosome has no covered sites the row is emitted with covered_bases=0,
sum_depth=0, mean_all_depth=0, mean_covered=NA, median_covered=NA.
"""
import argparse
import sys
from collections import defaultdict


def median_from_hist(hist, n):
    """Median of a value-count histogram covering n observations (n > 0)."""
    # lower & upper middle indices (0-based) for an n-length sorted sequence
    lo_idx = (n - 1) // 2
    hi_idx = n // 2
    lo_val = hi_val = None
    cum = 0
    for depth in sorted(hist):
        cum += hist[depth]
        if lo_val is None and cum > lo_idx:
            lo_val = depth
        if cum > hi_idx:
            hi_val = depth
            break
    return (lo_val + hi_val) / 2.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--chrom", required=True, help="the single chromosome this job covers")
    ap.add_argument("--lengths", required=True, help="TSV: chrom<TAB>length")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    lengths = {}
    with open(args.lengths) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            c, ln = line.split("\t")[:2]
            lengths[c] = int(ln)
    if args.chrom not in lengths:
        sys.exit(f"ERROR: chrom {args.chrom} not found in lengths file {args.lengths}")
    length = lengths[args.chrom]

    # depth histogram over covered sites of this one chromosome
    hist = defaultdict(int)
    sum_depth = 0
    covered = 0

    for line in sys.stdin:
        # chrom \t pos \t depth
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        if parts[0] != args.chrom:
            continue
        d = int(parts[2])
        if d <= 0:
            continue
        hist[d] += 1
        sum_depth += d
        covered += 1

    mean_all = sum_depth / length if length > 0 else 0.0
    if covered > 0:
        mean_cov_s = f"{sum_depth / covered:.6f}"
        med_cov_s = f"{median_from_hist(hist, covered):.1f}"
    else:
        mean_cov_s = "NA"
        med_cov_s = "NA"

    with open(args.out, "w") as out:
        out.write("sample\tchrom\tlength\tcovered_bases\tsum_depth\t"
                  "mean_all_depth\tmean_covered\tmedian_covered\n")
        out.write(f"{args.sample}\t{args.chrom}\t{length}\t{covered}\t{sum_depth}\t"
                  f"{mean_all:.6f}\t{mean_cov_s}\t{med_cov_s}\n")


if __name__ == "__main__":
    main()
