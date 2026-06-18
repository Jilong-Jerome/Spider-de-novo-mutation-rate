#!/usr/bin/env python3
"""
remap_bim_chrom.py

KING (and many downstream tools) only treat *integer* chromosome codes as
autosomes; non-integer scaffold names (e.g. ``mim_1``) are silently ignored,
which would leave KING with zero usable markers.

This script rewrites the chromosome column of a PLINK ``.bim`` file in place,
mapping each distinct contig name to a sequential integer (1..N, ordered by
first appearance) and writing the mapping to a companion ``.chrommap.tsv`` file.

Because kinship/IBS0 must be computed on AUTOSOMES ONLY (sex chromosomes
distort relatedness as males are hemizygous), the script also asserts that no
sex-chromosome contig leaked through the upstream autosome filter. A contig is
flagged as a sex chromosome if its trailing token (after the last ``_``) is not
purely numeric and begins with ``X``/``Y``/``Z``/``W`` (case-insensitive), e.g.
``mim_X1`` / ``mim_X2``.

Usage:
    python remap_bim_chrom.py <input.bim> <output.chrommap.tsv>
"""
import sys


def _is_sex_contig(name):
    token = name.rsplit('_', 1)[-1]
    return bool(token) and token[0].upper() in ('X', 'Y', 'Z', 'W') and not token.isdigit()


def remap_bim(bim_path, map_path):
    # First pass: read all rows, collect contig order.
    with open(bim_path) as fh:
        rows = [line.rstrip('\n').split('\t') for line in fh if line.strip()]

    order = []
    seen = set()
    for fields in rows:
        chrom = fields[0]
        if chrom not in seen:
            seen.add(chrom)
            order.append(chrom)

    sex = [c for c in order if _is_sex_contig(c)]
    if sex:
        sys.exit(
            "ERROR: sex-chromosome contigs reached the kinship marker set "
            "(autosomes-only expected): %s" % ", ".join(sex)
        )

    chrom_to_int = {c: str(i + 1) for i, c in enumerate(order)}

    # Write the mapping for traceability.
    with open(map_path, 'w') as out:
        out.write("original_contig\tinteger_code\n")
        for c in order:
            out.write("%s\t%s\n" % (c, chrom_to_int[c]))

    # Second pass: rewrite the .bim in place with integer chromosome codes.
    with open(bim_path, 'w') as out:
        for fields in rows:
            fields[0] = chrom_to_int[fields[0]]
            out.write("\t".join(fields) + "\n")

    print("Remapped %d contigs to integer codes (autosomes only)." % len(order))


if __name__ == "__main__":
    bim_path = sys.argv[1]
    map_path = sys.argv[2]
    remap_bim(bim_path, map_path)
