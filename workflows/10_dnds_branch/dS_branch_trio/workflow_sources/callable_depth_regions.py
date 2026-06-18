#!/usr/bin/env python3
"""Create callable-depth stats and an uncallable BED from genomecov bedGraph."""

import argparse
import math
from collections import defaultdict


def read_depth_lengths(path):
    covered_depth_lengths = defaultdict(int)
    total_bases = 0
    covered_bases = 0
    zero_depth_bases = 0

    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            chrom, start, end, depth = line.rstrip("\n").split("\t")[:4]
            start = int(start)
            end = int(end)
            depth = int(float(depth))
            length = end - start
            if length <= 0:
                continue
            total_bases += length
            if depth > 0:
                covered_depth_lengths[depth] += length
                covered_bases += length
            else:
                zero_depth_bases += length

    return covered_depth_lengths, total_bases, covered_bases, zero_depth_bases


def weighted_quantile(depth_lengths, n_bases, quantile):
    """Return the smallest depth whose cumulative base count reaches `quantile`."""
    if n_bases == 0:
        raise ValueError("no covered bases found in genomecov bedGraph")

    target = quantile * n_bases
    seen = 0
    last_depth = None
    for depth in sorted(depth_lengths):
        seen += depth_lengths[depth]
        last_depth = depth
        if seen >= target:
            return depth

    return last_depth


def weighted_median(depth_lengths, n_bases):
    return weighted_quantile(depth_lengths, n_bases, 0.5)


def weighted_mean(depth_lengths, n_bases):
    if n_bases == 0:
        raise ValueError("no covered bases found in genomecov bedGraph")
    total = sum(depth * length for depth, length in depth_lengths.items())
    return total / n_bases


def write_depth_distribution(dist_path, depth_lengths, n_bases):
    """Write a per-depth histogram of covered sites (depth > 0)."""
    cumulative = 0
    with open(dist_path, "w") as out:
        out.write("depth\tn_bases\tfrac_of_covered\tcumulative_frac\n")
        for depth in sorted(depth_lengths):
            length = depth_lengths[depth]
            cumulative += length
            frac = length / n_bases if n_bases else 0.0
            cum_frac = cumulative / n_bases if n_bases else 0.0
            out.write(f"{depth}\t{length}\t{frac:.6f}\t{cum_frac:.6f}\n")


def write_stats(stats_path, total_bases, covered_bases, zero_depth_bases,
                covered_median_depth, covered_mean_depth, percentiles,
                min_depth, min_factor, max_factor):
    raw_min = covered_median_depth * min_factor
    raw_max = covered_median_depth * max_factor
    callable_min = max(min_depth, int(math.ceil(raw_min)))
    callable_max = int(math.floor(raw_max))
    covered_fraction = covered_bases / total_bases if total_bases else 0.0

    with open(stats_path, "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"total_genome_bases\t{total_bases}\n")
        out.write(f"covered_bases\t{covered_bases}\n")
        out.write(f"zero_depth_bases\t{zero_depth_bases}\n")
        out.write(f"covered_fraction\t{covered_fraction:.6f}\n")
        out.write(f"covered_median_depth\t{covered_median_depth}\n")
        out.write(f"covered_mean_depth\t{covered_mean_depth:.4f}\n")
        for label, value in percentiles:
            out.write(f"covered_depth_{label}\t{value}\n")
        out.write(f"min_depth_floor\t{min_depth}\n")
        out.write(f"coverage_median_min_factor\t{min_factor}\n")
        out.write(f"coverage_median_max_factor\t{max_factor}\n")
        out.write(f"callable_min_depth\t{callable_min}\n")
        out.write(f"callable_max_depth\t{callable_max}\n")

    return callable_min, callable_max


def flush_interval(out, current):
    if current is not None:
        chrom, start, end = current
        out.write(f"{chrom}\t{start}\t{end}\n")


def write_uncallable_bed(genomecov_path, bed_path, callable_min, callable_max):
    current = None
    with open(genomecov_path) as genomecov, open(bed_path, "w") as out:
        for line in genomecov:
            if not line.strip():
                continue
            chrom, start, end, depth = line.rstrip("\n").split("\t")[:4]
            start = int(start)
            end = int(end)
            depth = int(float(depth))
            if depth < callable_min or depth > callable_max:
                if current is None:
                    current = [chrom, start, end]
                elif chrom == current[0] and start <= current[2]:
                    current[2] = max(current[2], end)
                else:
                    flush_interval(out, current)
                    current = [chrom, start, end]
        flush_interval(out, current)


def main():
    parser = argparse.ArgumentParser(
        description="Compute median-depth callable bounds from bedtools genomecov -bga output")
    parser.add_argument("--genomecov", required=True,
                        help="bedGraph from bedtools genomecov -bga")
    parser.add_argument("--stats", required=True,
                        help="Output TSV of depth summary and callable bounds")
    parser.add_argument("--uncallable-bed", required=True,
                        help="Output BED of intervals outside callable bounds")
    parser.add_argument("--depth-distribution",
                        help="Optional output TSV with the per-depth histogram "
                             "of covered sites")
    parser.add_argument("--min-depth", type=int, required=True,
                        help="Absolute lower depth floor")
    parser.add_argument("--min-factor", type=float, default=0.5,
                        help="Minimum callable median-depth multiplier")
    parser.add_argument("--max-factor", type=float, default=2.0,
                        help="Maximum callable median-depth multiplier")
    args = parser.parse_args()

    covered_depth_lengths, total_bases, covered_bases, zero_depth_bases = (
        read_depth_lengths(args.genomecov)
    )
    covered_median_depth = weighted_median(covered_depth_lengths, covered_bases)
    covered_mean_depth = weighted_mean(covered_depth_lengths, covered_bases)
    percentiles = [
        (label, weighted_quantile(covered_depth_lengths, covered_bases, q))
        for label, q in (("p05", 0.05), ("p25", 0.25), ("p50", 0.50),
                          ("p75", 0.75), ("p95", 0.95))
    ]

    if args.depth_distribution:
        write_depth_distribution(
            args.depth_distribution, covered_depth_lengths, covered_bases)

    callable_min, callable_max = write_stats(
        stats_path=args.stats,
        total_bases=total_bases,
        covered_bases=covered_bases,
        zero_depth_bases=zero_depth_bases,
        covered_median_depth=covered_median_depth,
        covered_mean_depth=covered_mean_depth,
        percentiles=percentiles,
        min_depth=args.min_depth,
        min_factor=args.min_factor,
        max_factor=args.max_factor,
    )
    write_uncallable_bed(
        genomecov_path=args.genomecov,
        bed_path=args.uncallable_bed,
        callable_min=callable_min,
        callable_max=callable_max,
    )

    print(f"Genome bases: {total_bases}")
    print(f"Covered bases: {covered_bases}")
    print(f"Covered median depth: {covered_median_depth}")
    print(f"Covered mean depth: {covered_mean_depth:.4f}")
    print(f"Callable depth: {callable_min}-{callable_max}")


if __name__ == "__main__":
    main()
