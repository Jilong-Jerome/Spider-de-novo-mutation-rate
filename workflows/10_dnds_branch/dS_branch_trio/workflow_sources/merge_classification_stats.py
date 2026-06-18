#!/usr/bin/env python3
"""Merge per-chromosome classification_stats TSVs into one whole-genome TSV.

All metrics are summed across input shards except
het_strict_filtered_fraction, which is recomputed from the merged
het_clean / het_unclean totals (per-chrom fractions are not additive).
"""

import argparse
from collections import OrderedDict


FRACTION_METRIC = "het_strict_filtered_fraction"


def read_stats(path):
    metrics = OrderedDict()
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header[:2] != ["metric", "value"]:
            raise ValueError(f"Unexpected header in {path}: {header}")
        for line in handle:
            if not line.strip():
                continue
            metric, value = line.rstrip("\n").split("\t", 1)
            metrics[metric] = value
    return metrics


def main():
    parser = argparse.ArgumentParser(
        description="Merge per-chromosome classification stats shards")
    parser.add_argument("--inputs", required=True, nargs="+",
                        help="Per-chromosome classification_stats TSVs")
    parser.add_argument("--output", required=True,
                        help="Merged whole-genome stats TSV")
    args = parser.parse_args()

    merged = OrderedDict()
    for path in args.inputs:
        for metric, value in read_stats(path).items():
            if metric == FRACTION_METRIC:
                continue
            try:
                summed = merged.get(metric, 0) + int(value)
            except ValueError:
                # Non-integer values (shouldn't occur for the additive
                # metrics) — keep the most recent.
                summed = value
            merged[metric] = summed

    het_clean = int(merged.get("het_clean", 0))
    het_unclean = int(merged.get("het_unclean", 0))
    het_total = het_clean + het_unclean
    fraction = (het_unclean / het_total) if het_total else 0.0

    ordered = OrderedDict()
    for metric in (
            "total_filtered_snps", "hom_ref", "hom_alt",
            "het_total", "het_clean", "het_unclean"):
        if metric in merged:
            ordered[metric] = merged.pop(metric)
    ordered[FRACTION_METRIC] = f"{fraction:.6f}"
    if "uncallable" in merged:
        ordered["uncallable"] = merged.pop("uncallable")
    for metric in sorted(merged):
        ordered[metric] = merged[metric]

    with open(args.output, "w") as out:
        out.write("metric\tvalue\n")
        for metric, value in ordered.items():
            out.write(f"{metric}\t{value}\n")


if __name__ == "__main__":
    main()
