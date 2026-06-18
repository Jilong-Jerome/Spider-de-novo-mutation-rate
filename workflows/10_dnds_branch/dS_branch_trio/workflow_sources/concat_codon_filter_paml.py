#!/usr/bin/env python3
"""
concat_codon_filter_paml.py - Filter a concatenated alignment to fully-callable
codons and write a PAML-compatible phylip file.

Replaces the old region-finding/polymorphism-filtering pipeline (align_filter.py,
find_region.R, filter_region.R, align_filter_region_fa.py, phylip_concat_paml.py).

The trio gene sequences are called per species against a third reference after
alignment QC, so the 3-species alignment is positionally exact. The only
filtering needed is to drop codons that are not fully callable: a codon is kept
only if every species has an unambiguous base (A/C/G/T) at all 3 positions.
Codons containing a gap ('-') or ambiguity (e.g. 'N') in any species are dropped.

Usage:
    python3 concat_codon_filter_paml.py <input.fa.cat> <output.paml.phy>
"""

import sys

VALID_BASES = set("ACGT")


def read_fasta(path):
    """Parse a simple multi-line FASTA into an ordered list of (name, seq)."""
    names = []
    seqs = {}
    current = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                current = line[1:]
                names.append(current)
                seqs[current] = []
            else:
                seqs[current].append(line)
    return [(name, ''.join(seqs[name]).upper()) for name in names]


def main():
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    records = read_fasta(input_file)
    if not records:
        print(f"ERROR: no sequences found in {input_file}")
        sys.exit(1)

    names = [name for name, _ in records]
    sequences = [seq for _, seq in records]

    seq_len = len(sequences[0])
    if any(len(s) != seq_len for s in sequences):
        print("ERROR: input sequences are not all the same length")
        sys.exit(1)
    if seq_len % 3 != 0:
        print(f"ERROR: alignment length {seq_len} is not a multiple of 3")
        sys.exit(1)

    n_codons = seq_len // 3
    kept = {name: [] for name in names}
    n_kept = 0

    for c in range(n_codons):
        i = c * 3
        codons = [s[i:i + 3] for s in sequences]
        if all(set(codon) <= VALID_BASES for codon in codons):
            for name, codon in zip(names, codons):
                kept[name].append(codon)
            n_kept += 1

    kept_bp = n_kept * 3
    print(f"Total codons: {n_codons}")
    print(f"Fully-callable codons kept: {n_kept} ({kept_bp} bp)")

    if n_kept == 0:
        print("ERROR: no fully-callable codons remain after filtering")
        sys.exit(1)

    with open(output_file, 'w') as out:
        out.write(f" {len(names)} {kept_bp}\n")
        for name in names:
            # PAML phylip: name followed by two spaces then sequence
            out.write(f"{name}  {''.join(kept[name])}\n")

    print(f"Written: {output_file}")


if __name__ == '__main__':
    main()
