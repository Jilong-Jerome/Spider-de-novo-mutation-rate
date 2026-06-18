#!/usr/bin/env python3
"""
calculate_dxy.py - Calculate Dxy (absolute divergence) between female and male

Dxy is the average number of pairwise nucleotide differences between two individuals.
For diploid data at each site where both individuals are callable:
- Compare all 4 haplotype pairs (2 female alleles x 2 male alleles)
- Each pair contributes 0 (same) or 0.25 (different) to the site's dxy
- Sum across all callable sites and divide by total callable sites
"""

import argparse
import os
import sys


def parse_genotype_alleles(ref, alt, genotype):
    """
    Parse genotype and return the two alleles.

    Args:
        ref: Reference allele
        alt: Alternate allele (can be '.' for non-variant)
        genotype: Genotype string like '0/0', '0/1', '1/1'

    Returns:
        Tuple of (allele1, allele2) or None if invalid
    """
    if genotype in ['.', './.', '.|.', 'NA']:
        return None

    # Build allele list
    alleles = [ref]
    if alt and alt != '.' and alt != '<NON_REF>':
        alleles.extend(alt.split(','))

    try:
        gt = genotype.replace('|', '/')
        parts = gt.split('/')
        if len(parts) != 2:
            return None

        idx1 = int(parts[0])
        idx2 = int(parts[1])

        # Handle indices that may reference non-existent alt alleles
        if idx1 >= len(alleles):
            idx1 = 0  # Default to ref
        if idx2 >= len(alleles):
            idx2 = 0  # Default to ref

        return (alleles[idx1], alleles[idx2])

    except (ValueError, IndexError):
        return None


def calculate_site_dxy(female_alleles, male_alleles):
    """
    Calculate Dxy for a single site.

    Compare all 4 haplotype pairs:
    - (female_a1, male_a1), (female_a1, male_a2)
    - (female_a2, male_a1), (female_a2, male_a2)

    Returns fraction of differing pairs (0, 0.25, 0.5, 0.75, or 1)
    """
    f1, f2 = female_alleles
    m1, m2 = male_alleles

    pairs = [
        (f1, m1), (f1, m2),
        (f2, m1), (f2, m2)
    ]

    diff_count = sum(1 for f, m in pairs if f != m)
    return diff_count / 4.0


def load_genotypes(genotypes_file):
    """
    Load genotypes file into dictionary indexed by position.

    Returns:
        dict: {pos: {'chrom', 'ref', 'alt', 'genotype', 'is_callable', 'region_type'}}
    """
    data = {}

    with open(genotypes_file, 'r') as f:
        header = f.readline().strip().split('\t')
        col_idx = {name: i for i, name in enumerate(header)}

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < len(header):
                continue

            pos = int(parts[col_idx['pos']])
            data[pos] = {
                'chrom': parts[col_idx['chrom']],
                'ref': parts[col_idx['ref']],
                'alt': parts[col_idx['alt']],
                'genotype': parts[col_idx['genotype']],
                'is_callable': parts[col_idx['is_callable']] == 'True',
                'region_type': parts[col_idx['region_type']]
            }

    return data


def calculate_dxy(female_genotypes_file, male_genotypes_file, gene, family, output_file):
    """
    Calculate Dxy between female and male for a gene.

    Args:
        female_genotypes_file: Path to female genotypes TSV
        male_genotypes_file: Path to male genotypes TSV
        gene: Gene name
        family: Family identifier
        output_file: Output TSV file
    """
    # Load genotype data
    female_data = load_genotypes(female_genotypes_file)
    male_data = load_genotypes(male_genotypes_file)

    # Find common positions
    all_positions = set(female_data.keys()) | set(male_data.keys())

    # Initialize counters for each region
    counts = {
        'whole_gene': {'callable_both': 0, 'dxy_sum': 0.0},
        'exon_only': {'callable_both': 0, 'dxy_sum': 0.0},
        'intron_only': {'callable_both': 0, 'dxy_sum': 0.0}
    }

    for pos in sorted(all_positions):
        # Skip if position not in both datasets
        if pos not in female_data or pos not in male_data:
            continue

        f_data = female_data[pos]
        m_data = male_data[pos]

        # Both must be callable
        if not (f_data['is_callable'] and m_data['is_callable']):
            continue

        # Parse alleles
        female_alleles = parse_genotype_alleles(
            f_data['ref'], f_data['alt'], f_data['genotype']
        )
        male_alleles = parse_genotype_alleles(
            m_data['ref'], m_data['alt'], m_data['genotype']
        )

        if female_alleles is None or male_alleles is None:
            continue

        # Calculate site Dxy
        site_dxy = calculate_site_dxy(female_alleles, male_alleles)

        # Update whole gene counts
        counts['whole_gene']['callable_both'] += 1
        counts['whole_gene']['dxy_sum'] += site_dxy

        # Update region-specific counts (use female's region type, should be same)
        region = f_data['region_type']
        if region == 'exon':
            counts['exon_only']['callable_both'] += 1
            counts['exon_only']['dxy_sum'] += site_dxy
        elif region == 'intron':
            counts['intron_only']['callable_both'] += 1
            counts['intron_only']['dxy_sum'] += site_dxy

    # Calculate final Dxy values
    results = []
    for region in ['whole_gene', 'exon_only', 'intron_only']:
        callable_sites = counts[region]['callable_both']
        dxy_sum = counts[region]['dxy_sum']

        if callable_sites > 0:
            dxy = dxy_sum / callable_sites
        else:
            dxy = 'NA'

        results.append({
            'gene': gene,
            'family': family,
            'region': region,
            'callable_sites_both': callable_sites,
            'dxy_sum': dxy_sum,
            'dxy': dxy
        })

    # Write output
    with open(output_file, 'w') as f:
        f.write("gene\tfamily\tregion\tcallable_sites_both\tdxy_sum\tdxy\n")
        for r in results:
            dxy_str = f"{r['dxy']:.6f}" if isinstance(r['dxy'], float) else r['dxy']
            dxy_sum_str = f"{r['dxy_sum']:.4f}"
            f.write(f"{r['gene']}\t{r['family']}\t{r['region']}\t"
                    f"{r['callable_sites_both']}\t{dxy_sum_str}\t{dxy_str}\n")

    # Print summary
    print(f"Dxy Summary for {family} - {gene}:")
    print("-" * 60)
    for r in results:
        dxy_str = f"{r['dxy']:.6f}" if isinstance(r['dxy'], float) else r['dxy']
        print(f"  {r['region']}: {r['callable_sites_both']} sites, Dxy = {dxy_str}")
    print(f"\nOutput: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Calculate Dxy between female and male'
    )
    parser.add_argument('--female', required=True, help='Female genotypes TSV file')
    parser.add_argument('--male', required=True, help='Male genotypes TSV file')
    parser.add_argument('--gene', required=True, help='Gene name')
    parser.add_argument('--family', required=True, help='Family identifier')
    parser.add_argument('--output', required=True, help='Output TSV file')

    args = parser.parse_args()

    # Create output directory if needed
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    calculate_dxy(
        female_genotypes_file=args.female,
        male_genotypes_file=args.male,
        gene=args.gene,
        family=args.family,
        output_file=args.output
    )


if __name__ == "__main__":
    main()
