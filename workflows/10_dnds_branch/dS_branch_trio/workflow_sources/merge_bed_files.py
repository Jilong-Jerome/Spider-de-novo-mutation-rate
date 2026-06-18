#!/usr/bin/env python3
"""Merge a large sorted BED with a smaller BED without shell sort."""

import argparse
from collections import defaultdict


def parse_bed_line(line):
    chrom, start, end = line.rstrip("\n").split("\t")[:3]
    return chrom, int(start), int(end)


def load_secondary(path):
    intervals = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            chrom, start, end = parse_bed_line(line)
            intervals[chrom].append((start, end))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def merge_intervals(primary_iter, secondary):
    current = None
    seen_chroms = set()
    secondary_idx = defaultdict(int)

    def emit(interval):
        if interval is not None:
            chrom, start, end = interval
            print(f"{chrom}\t{start}\t{end}")

    def add_interval(chrom, start, end):
        nonlocal current
        seen_chroms.add(chrom)
        if current is None:
            current = [chrom, start, end]
        elif chrom == current[0] and start <= current[2]:
            current[2] = max(current[2], end)
        else:
            emit(current)
            current = [chrom, start, end]

    def add_secondary_before(chrom, before_start):
        intervals = secondary.get(chrom, [])
        idx = secondary_idx[chrom]
        while idx < len(intervals) and intervals[idx][0] <= before_start:
            sec_start, sec_end = intervals[idx]
            add_interval(chrom, sec_start, sec_end)
            idx += 1
        secondary_idx[chrom] = idx

    def add_remaining_secondary(chrom):
        intervals = secondary.get(chrom, [])
        idx = secondary_idx[chrom]
        while idx < len(intervals):
            sec_start, sec_end = intervals[idx]
            add_interval(chrom, sec_start, sec_end)
            idx += 1
        secondary_idx[chrom] = idx

    active_chrom = None
    for line in primary_iter:
        if not line.strip():
            continue
        chrom, start, end = parse_bed_line(line)
        if active_chrom is not None and chrom != active_chrom:
            add_remaining_secondary(active_chrom)
        active_chrom = chrom
        add_secondary_before(chrom, start)
        add_interval(chrom, start, end)

    if active_chrom is not None:
        add_remaining_secondary(active_chrom)

    for chrom in sorted(set(secondary) - seen_chroms):
        add_remaining_secondary(chrom)

    emit(current)


def main():
    parser = argparse.ArgumentParser(
        description="Merge a large sorted primary BED with a smaller BED")
    parser.add_argument("--primary", required=True,
                        help="Large sorted BED, usually callable-depth mask")
    parser.add_argument("--secondary", required=True,
                        help="Smaller BED to add, usually rejected SNP sites")
    args = parser.parse_args()

    secondary = load_secondary(args.secondary)
    with open(args.primary) as primary:
        merge_intervals(primary, secondary)


if __name__ == "__main__":
    main()
