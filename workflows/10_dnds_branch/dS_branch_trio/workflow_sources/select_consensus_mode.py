#!/usr/bin/env python3
"""Select mode-specific consensus SNPs from shared SNP classifications."""

import argparse
import csv
import sys

from select_consensus_snps import VALID_MODES, deterministic_alt_choice

VALID_SUPPORT_POLICIES = {"strict", "relaxed"}


def classification_rows(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            yield row


def select_action(fields, row, mode, seed, species, support_policy):
    classification = row["classification"]
    if classification == "hom_alt":
        return "alt"
    if classification == "hom_ref":
        return "ref"
    if classification == "het_clean":
        is_selectable_het = True
    elif classification == "het_unclean":
        is_selectable_het = support_policy == "relaxed"
    else:
        is_selectable_het = False

    if not is_selectable_het:
        return "uncallable"

    if mode == "using_alt":
        return "alt"
    if mode == "using_ref":
        return "ref"
    if mode == "random":
        if deterministic_alt_choice(seed, species, fields[0], fields[1]):
            return "alt"
        return "ref"
    if mode == "majority":
        try:
            ref_support = int(row["dp4_ref_support"])
            alt_support = int(row["dp4_alt_support"])
        except ValueError:
            return "uncallable"
        if alt_support > ref_support:
            return "alt"
        return "ref"

    raise ValueError(f"unsupported mode: {mode}")


def validate_pair(fields, row):
    expected = (row["chrom"], row["pos"], row["ref"], row["alt"])
    observed = (fields[0], fields[1], fields[3], fields[4])
    if observed != expected:
        raise ValueError(
            "VCF/classification mismatch: "
            f"VCF={observed}, classification={expected}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=sorted(VALID_MODES))
    parser.add_argument("--support-policy", required=True,
                        choices=sorted(VALID_SUPPORT_POLICIES),
                        help="strict masks het_unclean; relaxed resolves it by mode")
    parser.add_argument("--seed", type=int, default=20260514)
    parser.add_argument("--species", required=True)
    parser.add_argument("--classification-tsv", required=True,
                        help="TSV produced by classify_consensus_snps.py")
    parser.add_argument("--rejected-bed", required=True,
                        help="BED of uncallable SNP sites to mask")
    args = parser.parse_args()

    rows = classification_rows(args.classification_tsv)
    seen_records = 0

    with open(args.rejected_bed, "w") as rejected:
        for line in sys.stdin:
            if line.startswith("#"):
                sys.stdout.write(line)
                continue
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            try:
                row = next(rows)
            except StopIteration as exc:
                raise ValueError(
                    "classification TSV ended before VCF records") from exc
            seen_records += 1
            validate_pair(fields, row)
            action = select_action(
                fields, row, args.mode, args.seed, args.species,
                args.support_policy)
            if action == "alt":
                sys.stdout.write(line)
            elif action == "uncallable":
                start = int(fields[1]) - 1
                rejected.write(f"{fields[0]}\t{start}\t{fields[1]}\n")

    try:
        extra = next(rows)
    except StopIteration:
        return
    raise ValueError(
        "classification TSV has extra rows after VCF ended: "
        f"{extra['chrom']}:{extra['pos']}")


if __name__ == "__main__":
    main()
