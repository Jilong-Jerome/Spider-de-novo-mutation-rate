#!/usr/bin/env python3
"""
Add a sequencing coverage column to the sequencing summary TSV.

size_gb is already the uncompressed FASTQ byte count (from zcat | wc -c),
so coverage is computed directly as:

    coverage = (size_gb * 1024^3) / genome_size_bases

Coverage is computed for every row so that per-lane contributions are visible;
the TOTAL row therefore gives each individual's whole-genome coverage estimate.
"""

import argparse
import os

SPECIES_TO_PREFIX = {
    'S_africanus':   'AFR',
    'S_bicolor':     'BIC',
    'S_dumicola':    'DUM',
    'S_lineatus':    'LIN',
    'S_mimosarum':   'MIM',
    'S_sarasinorum': 'SAR',
    'S_tentoriicola':'TEN',
}


def load_genome_sizes(data_dir):
    """Return {species: genome_size_in_bases} from .fai files."""
    sizes = {}
    for species, prefix in SPECIES_TO_PREFIX.items():
        fai_path = os.path.join(data_dir, f'{prefix}_ncbi_chromosome.fa.fai')
        total_bases = 0
        with open(fai_path) as fh:
            for line in fh:
                total_bases += int(line.split('\t')[1])
        sizes[species] = total_bases
    return sizes


def main():
    parser = argparse.ArgumentParser(
        description='Append coverage column to sequencing summary TSV.')
    parser.add_argument('--input',      required=True,
                        help='Path to sequencing_summary.tsv')
    parser.add_argument('--data_dir',   required=True,
                        help='Directory containing {PREFIX}_ncbi_chromosome.fa.fai files')
    parser.add_argument('--output',     required=True,
                        help='Output TSV path')
    args = parser.parse_args()

    genome_sizes = load_genome_sizes(args.data_dir)

    with open(args.input) as fin, open(args.output, 'w') as fout:
        header = fin.readline().rstrip('\n')
        fout.write(header + '\tcoverage\n')

        for line in fin:
            line = line.rstrip('\n')
            parts = line.split('\t')
            species, family, individual, round_id, lane, size_gb = parts

            genome_bases = genome_sizes.get(species)
            if genome_bases is None:
                coverage = 'NA'
            else:
                bytes_total = float(size_gb) * 1024 ** 3
                coverage = f'{bytes_total / genome_bases:.2f}'

            fout.write(line + f'\t{coverage}\n')


if __name__ == '__main__':
    main()
