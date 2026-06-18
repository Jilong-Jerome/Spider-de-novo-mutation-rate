#!/usr/bin/env python3
"""
select_consensus_snps.py - Resolve SNPs for consensus modes.

Reads a VCF from stdin and writes a VCF to stdout. Records whose selected
allele is REF are omitted, so bcftools consensus keeps the reference base.
Records whose selected allele is ALT are kept.
"""

import argparse
import hashlib
import sys
from collections import Counter

try:
    import pysam
except ImportError:
    pysam = None


VALID_MODES = {"majority", "random", "using_ref", "using_alt"}
VALID_BASES = {"A", "C", "G", "T"}


def parse_info(info):
    values = {}
    for field in info.split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            values[key] = value
        else:
            values[field] = True
    return values


def get_sample_value(format_field, sample_field, key):
    keys = format_field.split(":")
    vals = sample_field.split(":")
    try:
        idx = keys.index(key)
    except ValueError:
        return None
    if idx >= len(vals):
        return None
    return vals[idx]


def normalize_gt(gt):
    if gt is None:
        return None
    return gt.replace("|", "/")


def deterministic_alt_choice(seed, species, chrom, pos):
    token = f"{seed}:{species}:{chrom}:{pos}".encode()
    digest = hashlib.sha256(token).digest()
    return digest[0] % 2 == 1


def dp4_support(info):
    dp4 = info.get("DP4")
    if dp4 is None:
        return None
    try:
        ref_fwd, ref_rev, alt_fwd, alt_rev = [int(x) for x in dp4.split(",")[:4]]
    except ValueError:
        return None
    return ref_fwd + ref_rev, alt_fwd + alt_rev


def count_bases_from_pileup_column(pileup_col, counts=None):
    if counts is None:
        counts = Counter()
    for read in pileup_col.pileups:
        if read.is_del or read.is_refskip:
            continue
        query_pos = read.query_position
        if query_pos is None:
            continue
        base = read.alignment.query_sequence[query_pos].upper()
        if base in VALID_BASES:
            counts[base] += 1
    return counts


def observed_bases_from_pileup(bam_handle, chrom, pos, min_base_quality):
    counts = Counter()
    start = int(pos) - 1
    end = int(pos)

    try:
        pileup_iter = bam_handle.pileup(
            chrom, start, end, truncate=True,
            min_base_quality=min_base_quality)
        for pileup_col in pileup_iter:
            if pileup_col.reference_pos != start:
                continue
            count_bases_from_pileup_column(pileup_col, counts)
    except ValueError:
        return counts

    return counts


def clean_heterozygous_support(base_counts, ref, alt):
    observed = {base for base, count in base_counts.items() if count > 0}
    return (
        observed
        and observed.issubset({ref, alt})
        and ref in observed
        and alt in observed
    )


def keep_alt_for_record(fields, mode, seed, species, bam_handle=None,
                        min_base_quality=20):
    return select_record_action(
        fields, mode, seed, species, bam_handle, min_base_quality) == "alt"


def select_record_action(fields, mode, seed, species, bam_handle=None,
                         min_base_quality=20):
    chrom = fields[0]
    pos = fields[1]
    ref = fields[3].upper()
    alt = fields[4].upper()

    # Keep this helper deliberately restricted to simple biallelic SNPs.
    # Other records are not useful for the CDS/PAML path and are omitted.
    if (
            len(ref) != 1 or len(alt) != 1
            or ref not in VALID_BASES or alt not in VALID_BASES
            or "," in alt):
        return "uncallable"

    info = parse_info(fields[7])
    gt = None
    if len(fields) >= 10:
        gt = normalize_gt(get_sample_value(fields[8], fields[9], "GT"))

    if gt in {"1/1", "1"}:
        return "alt"
    if gt in {"0/0", "0"}:
        return "ref"

    is_het = gt in {"0/1", "1/0"}
    if not is_het:
        return "uncallable"

    if bam_handle is None:
        return "uncallable"
    base_counts = observed_bases_from_pileup(
        bam_handle, chrom, pos, min_base_quality)
    if not clean_heterozygous_support(base_counts, ref, alt):
        return "uncallable"

    if mode == "using_alt":
        return "alt"
    if mode == "using_ref":
        return "ref"
    if mode == "random":
        if deterministic_alt_choice(seed, species, chrom, pos):
            return "alt"
        return "ref"
    if mode == "majority":
        support = dp4_support(info)
        if support is None:
            return "uncallable"
        ref_support, alt_support = support
        if alt_support > ref_support:
            return "alt"
        return "ref"

    raise ValueError(f"unsupported mode: {mode}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=sorted(VALID_MODES))
    parser.add_argument("--seed", type=int, default=20260514)
    parser.add_argument("--species", required=True)
    parser.add_argument("--bam", required=True,
                        help="Final BAM used to validate read support")
    parser.add_argument("--min-base-quality", type=int, default=20,
                        help="Minimum pileup base quality for support validation")
    parser.add_argument("--rejected-bed",
                        help="Optional BED of rejected non-reference SNP sites")
    args = parser.parse_args()

    if pysam is None:
        raise RuntimeError("pysam is required for clean read-support validation")

    rejected_handle = open(args.rejected_bed, "w") if args.rejected_bed else None
    try:
        with pysam.AlignmentFile(args.bam, "rb") as bam_handle:
            for line in sys.stdin:
                if line.startswith("#"):
                    sys.stdout.write(line)
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    continue
                action = select_record_action(
                        fields, args.mode, args.seed, args.species,
                        bam_handle=bam_handle,
                        min_base_quality=args.min_base_quality)
                if action == "alt":
                    sys.stdout.write(line)
                elif action == "uncallable" and rejected_handle is not None:
                    start = int(fields[1]) - 1
                    rejected_handle.write(f"{fields[0]}\t{start}\t{fields[1]}\n")
    finally:
        if rejected_handle is not None:
            rejected_handle.close()


if __name__ == "__main__":
    main()
