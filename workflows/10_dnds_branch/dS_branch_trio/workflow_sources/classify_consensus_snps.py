#!/usr/bin/env python3
"""Classify filtered SNPs once using BAM pileup support.

Uses one streamed pileup per chromosome so each BGZF block is decoded at
most once even when many het sites share overlapping reads.
"""

import argparse
from collections import Counter, defaultdict

try:
    import pysam
except ImportError:
    pysam = None

from select_consensus_snps import (
    VALID_BASES,
    clean_heterozygous_support,
    count_bases_from_pileup_column,
    dp4_support,
    get_sample_value,
    normalize_gt,
    parse_info,
)


def format_counts(counts):
    return ",".join(f"{base}:{counts.get(base, 0)}" for base in sorted(VALID_BASES))


def heterozygous_support_reason(base_counts, ref, alt):
    observed = {base for base, count in base_counts.items() if count > 0}
    if not observed:
        return "no_valid_bases"
    if observed - {ref, alt}:
        return "third_base_support"
    if ref not in observed:
        return "missing_ref_support"
    if alt not in observed:
        return "missing_alt_support"
    return "pass_het_clean"


def classify_genotype(fields):
    """Resolve a VCF record without BAM lookup.

    Returns (row, needs_pileup, het_context).
    - row: list of classification columns when fully resolved, else None.
    - needs_pileup: True for simple het SNPs that still need pileup support.
    - het_context: tuple used to finalize a het row once base counts are
      available: (chrom, pos0, ref, alt, gt_out, ref_support, alt_support).
    """
    chrom = fields[0]
    pos = fields[1]
    ref = fields[3].upper()
    alt = fields[4].upper()
    info = parse_info(fields[7])
    support = dp4_support(info)
    ref_support, alt_support = support if support is not None else (".", ".")

    gt = None
    if len(fields) >= 10:
        gt = normalize_gt(get_sample_value(fields[8], fields[9], "GT"))
    gt_out = gt if gt is not None else "."

    if (
            len(ref) != 1 or len(alt) != 1
            or ref not in VALID_BASES or alt not in VALID_BASES
            or "," in alt):
        return (["uncallable", gt_out, ref_support, alt_support, ".",
                 "non_simple_snp"], False, None)

    if gt in {"0/0", "0"}:
        return (["hom_ref", gt_out, ref_support, alt_support, ".",
                 "pass_hom_ref"], False, None)

    if gt in {"1/1", "1"}:
        return (["hom_alt", gt_out, ref_support, alt_support, ".",
                 "pass_hom_alt"], False, None)

    if gt not in {"0/1", "1/0"}:
        return (["uncallable", gt_out, ref_support, alt_support, ".",
                 "invalid_gt"], False, None)

    pos0 = int(pos) - 1
    return (None, True,
            (chrom, pos0, ref, alt, gt_out, ref_support, alt_support))


def finalize_het_row(het_context, base_counts):
    chrom, pos0, ref, alt, gt_out, ref_support, alt_support = het_context
    observed = format_counts(base_counts)
    reason = heterozygous_support_reason(base_counts, ref, alt)
    label = ("het_clean"
             if clean_heterozygous_support(base_counts, ref, alt)
             else "het_unclean")
    return [label, gt_out, ref_support, alt_support, observed, reason]


def write_stats(path, classification_counts, reason_counts):
    total = sum(classification_counts.values())
    het_clean = classification_counts.get("het_clean", 0)
    het_unclean = classification_counts.get("het_unclean", 0)
    het_total = het_clean + het_unclean
    strict_filtered_fraction = (
        het_unclean / het_total if het_total else 0.0
    )

    metrics = [
        ("total_filtered_snps", total),
        ("hom_ref", classification_counts.get("hom_ref", 0)),
        ("hom_alt", classification_counts.get("hom_alt", 0)),
        ("het_total", het_total),
        ("het_clean", het_clean),
        ("het_unclean", het_unclean),
        ("het_strict_filtered_fraction", f"{strict_filtered_fraction:.6f}"),
        ("uncallable", classification_counts.get("uncallable", 0)),
    ]
    for reason in sorted(reason_counts):
        metrics.append((f"reason_{reason}", reason_counts[reason]))

    with open(path, "w") as out:
        out.write("metric\tvalue\n")
        for metric, value in metrics:
            out.write(f"{metric}\t{value}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Classify filtered SNPs by clean BAM pileup support")
    parser.add_argument("--vcf", required=True,
                        help="Plain VCF produced after depth and quality filtering")
    parser.add_argument("--bam", required=True,
                        help="Final BAM used to validate read support")
    parser.add_argument("--min-base-quality", type=int, default=20,
                        help="Minimum pileup base quality for support validation")
    parser.add_argument("--output", required=True,
                        help="Output TSV with one row per non-header VCF record")
    parser.add_argument("--stats", required=True,
                        help="Output TSV with per-species classification counts")
    args = parser.parse_args()

    if pysam is None:
        raise RuntimeError("pysam is required for clean read-support validation")

    # First pass: stream the VCF, resolve genotype-only rows immediately,
    # buffer hets keyed by chromosome for a single streaming pileup.
    rows = []  # list of [chrom, pos, ref, alt] + classification columns OR None placeholder
    het_contexts = []  # parallel list; populated only for hets
    pending = defaultdict(dict)  # chrom -> {pos0: row_index}

    with open(args.vcf) as vcf:
        for line in vcf:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
            row_prefix = [chrom, pos, ref, alt]

            resolved, needs_pileup, het_context = classify_genotype(fields)
            if not needs_pileup:
                rows.append(row_prefix + resolved)
                het_contexts.append(None)
                continue

            idx = len(rows)
            rows.append(None)
            het_contexts.append((row_prefix, het_context))
            pending[het_context[0]][het_context[1]] = idx

    # Second pass: one streamed pileup per chromosome covering the span
    # of pending het positions.
    with pysam.AlignmentFile(args.bam, "rb") as bam_handle:
        for chrom, position_map in pending.items():
            if not position_map:
                continue
            min_pos0 = min(position_map)
            max_pos0 = max(position_map)
            try:
                pileup_iter = bam_handle.pileup(
                    chrom, min_pos0, max_pos0 + 1, truncate=True,
                    min_base_quality=args.min_base_quality)
                for pileup_col in pileup_iter:
                    ref_pos = pileup_col.reference_pos
                    idx = position_map.get(ref_pos)
                    if idx is None:
                        continue
                    counts = count_bases_from_pileup_column(pileup_col)
                    row_prefix, het_context = het_contexts[idx]
                    rows[idx] = row_prefix + finalize_het_row(het_context, counts)
            except ValueError:
                # Contig missing from BAM or empty region — fall through;
                # remaining unresolved rows are filled below with empty counts.
                pass

    # Any het row still unresolved (position absent from pileup) gets the
    # empty-Counter classification → het_unclean / no_valid_bases.
    empty_counter = Counter()
    classification_counts = Counter()
    reason_counts = Counter()

    with open(args.output, "w") as out:
        out.write(
            "chrom\tpos\tref\talt\tclassification\tgt\t"
            "dp4_ref_support\tdp4_alt_support\tobserved_bases\treason\n")
        for idx, row in enumerate(rows):
            if row is None:
                row_prefix, het_context = het_contexts[idx]
                row = row_prefix + finalize_het_row(het_context, empty_counter)
            classification_counts[row[4]] += 1
            reason_counts[row[-1]] += 1
            out.write("\t".join(str(value) for value in row) + "\n")

    write_stats(args.stats, classification_counts, reason_counts)


if __name__ == "__main__":
    main()
