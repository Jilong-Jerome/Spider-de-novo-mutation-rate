#!/usr/bin/env python3
"""
Summarize uncompressed read data per (round_id, lane) for one individual.

Output TSV columns (no header): species  family  individual  round_id  lane  size_gb
One row per unique (round_id, lane) plus one TOTAL row.
"""

import argparse
import os
import re
import subprocess
from collections import defaultdict


def parse_round_and_lane(filename):
    """
    Parse round_id and lane from a .fq.gz filename.

    Rules:
    1. Strip read-pair suffix: _1.fq.gz, _2.fq.gz, .1.fq.gz, .2.fq.gz
    2. If starts with PE{digits}_: strip that prefix
    3. Match ^([^_]+)_(L\d+)_ → round_id=group1, lane=group2
    4. Otherwise: round_id=NA, lane=NA
    """
    name = filename

    # Step 1: strip read-pair suffix
    for suffix in ('_1.fq.gz', '_2.fq.gz', '.1.fq.gz', '.2.fq.gz'):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break

    # Step 2: strip PE{digits}_ prefix
    name = re.sub(r'^PE\d+_', '', name)

    # Step 3: match round_id and lane
    m = re.match(r'^([^_]+)_(L\d+)_', name)
    if m:
        return m.group(1), m.group(2)

    return 'NA', 'NA'


def get_uncompressed_bytes(filepath):
    """Return uncompressed byte count of a gzip file via zcat | wc -c."""
    result = subprocess.run(
        f'zcat {filepath} | wc -c',
        shell=True,
        capture_output=True,
        text=True,
        check=True,
    )
    return int(result.stdout.strip())


def main():
    parser = argparse.ArgumentParser(description='Summarize sequencing data for one individual.')
    parser.add_argument('--individual_dir', required=True, help='Path to individual directory')
    parser.add_argument('--species', required=True)
    parser.add_argument('--family', required=True)
    parser.add_argument('--individual', required=True)
    parser.add_argument('--output', required=True, help='Output TSV path')
    args = parser.parse_args()

    # Collect .fq.gz files
    fq_files = [
        f for f in os.listdir(args.individual_dir)
        if f.endswith('.fq.gz')
    ]

    # Group files by (round_id, lane)
    groups = defaultdict(list)
    for fname in fq_files:
        round_id, lane = parse_round_and_lane(fname)
        filepath = os.path.join(args.individual_dir, fname)
        groups[(round_id, lane)].append(filepath)

    # Measure uncompressed bytes per group
    group_bytes = {}
    for key, filepaths in sorted(groups.items()):
        total = sum(get_uncompressed_bytes(fp) for fp in filepaths)
        group_bytes[key] = total

    total_bytes = sum(group_bytes.values())

    # Write output TSV (no header)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    with open(args.output, 'w') as fh:
        for (round_id, lane), nbytes in sorted(group_bytes.items()):
            size_gb = nbytes / 1024**3
            fh.write(f'{args.species}\t{args.family}\t{args.individual}\t'
                     f'{round_id}\t{lane}\t{size_gb:.4f}\n')
        # TOTAL row
        total_gb = total_bytes / 1024**3
        fh.write(f'{args.species}\t{args.family}\t{args.individual}\t'
                 f'TOTAL\tTOTAL\t{total_gb:.4f}\n')


if __name__ == '__main__':
    main()
