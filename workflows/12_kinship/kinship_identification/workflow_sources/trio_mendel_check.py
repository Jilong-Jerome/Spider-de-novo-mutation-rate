#!/usr/bin/env python3
"""
trio_mendel_check.py

Allele-frequency-free, per-trio validation of parent-offspring relationships.

Motivation
----------
Frequency-based kinship estimators (KING robust kinship, PLINK PI_HAT) collapse
under the extreme inbreeding / family structure of social spiders: true
parent-offspring pairs come out with nonsensical (even strongly negative) kinship
because the population allele-frequency baseline they assume is meaningless here.

This test uses pure Mendelian logic instead, one trio at a time. Within a putative
trio (the two parents F + M of a family and one offspring S_k), consider the sites
where the **two parents are fixed for opposite alleles** -- one parent ``0/0`` and
the other ``1/1`` (an "opposite-homozygote" / fixed-difference site). At such a
site a genuine offspring of *both* parents is **obligately heterozygous (0/1)**:
it must inherit the reference allele from the hom-ref parent and the alternate
allele from the hom-alt parent. So the fraction of opposite-homozygote sites at
which the offspring is actually heterozygous ("consistency") is ~1.0 for a real
trio (minus a little genotyping error) and drops sharply if either assigned parent
is wrong. This is independent of any allele frequency and therefore robust to
inbreeding.

Input
-----
A biallelic, autosomal, GQ-masked SNP VCF (the ``{SP}_autosome_snps_nomaf.vcf``
produced by the trio_filter step -- no MAF cutoff, since the informative-site
filter here is the opposite-homozygote condition, not population frequency).
Every genotype used is additionally required to meet a minimum depth (``--min_dp``,
default 26, the de novo mutation callability depth): at low depth a true het can
suffer allelic dropout and be miscalled homozygous, which would fabricate a false
opposite-homozygote parent site or a false offspring Mendelian error. The VCF
carries per-genotype DP (FORMAT ``GT:AD:DP:GQ:PL``), so this gate is applied here
rather than re-running the expensive trio_filter pass.

Trios are read from the VCF sample names ``{SP}_family{N}_{F|M|S#}_{role}``:
per family the two parents are the ``F`` (female) and ``M`` (male) samples and the
offspring are the ``S#`` samples. Each (F, M, S_k) is one trio / one output row.

Output
------
``{SP}_trio_mendel.tsv`` -- one row per offspring:
    species, family, parent_F, parent_M, offspring,
    n_informative, n_het, n_homref, n_homalt, consistency, error_rate, verdict

verdict:
    consistent          consistency >= consistency_min and n_informative >= min_sites
    INCONSISTENT        enough sites but consistency below threshold (trio not solid)
    insufficient-sites  fewer than min_sites opposite-homozygote sites to judge

Usage
-----
    python trio_mendel_check.py \
        --vcf {SP}_autosome_snps_nomaf.vcf \
        --species MIM \
        --output_tsv {SP}_trio_mendel.tsv \
        --min_gq 20 --min_dp 26 --consistency_min 0.9 --min_sites 100
"""
import argparse
import csv
import sys

import pysam


def parse_individual(name):
    """Return (family, role) from ``{SP}_family{N}_{F|M|S#}_{...}``.

    role is one of {'F', 'M', 'offspring', 'unknown'}. Mirrors the naming idiom
    used by classify_kinship.py:parse_individual, but keeps the parents' F/M
    distinction (needed to label the two founders).
    """
    parts = name.split('_')
    family = parts[1] if len(parts) > 1 else name
    token = parts[2] if len(parts) > 2 else ''
    t = token[:1].upper()
    if t == 'S':
        role = 'offspring'
    elif t == 'F':
        role = 'F'
    elif t == 'M':
        role = 'M'
    else:
        role = 'unknown'
    return family, role


def build_trios(samples):
    """Group VCF samples into per-family trios.

    Returns a list of dicts {family, parent_F, parent_M, offspring:[...]} for
    every family that has BOTH parents present. Families missing a parent are
    reported to stderr and skipped (cannot run the opposite-homozygote test).
    """
    families = {}
    for s in samples:
        fam, role = parse_individual(s)
        families.setdefault(fam, {'F': None, 'M': None, 'offspring': []})
        if role == 'F':
            families[fam]['F'] = s
        elif role == 'M':
            families[fam]['M'] = s
        elif role == 'offspring':
            families[fam]['offspring'].append(s)

    trios = []
    for fam in sorted(families):
        info = families[fam]
        if not info['F'] or not info['M']:
            sys.stderr.write(
                "WARN: family %s missing a parent (F=%s, M=%s); skipped\n"
                % (fam, info['F'], info['M']))
            continue
        if not info['offspring']:
            sys.stderr.write("WARN: family %s has no offspring; skipped\n" % fam)
            continue
        trios.append({
            'family': fam,
            'parent_F': info['F'],
            'parent_M': info['M'],
            'offspring': sorted(info['offspring']),
        })
    return trios


def gt_state(gt):
    """Classify a diploid GT tuple -> 'homref' | 'homalt' | 'het' | None.

    None means uncalled or non-informative (missing allele). Only the
    0/0, 1/1 and 0/1 states are recognised (input is biallelic SNPs).
    """
    if gt is None:
        return None
    a, b = gt[0], gt[1]
    if a is None or b is None:
        return None
    if a == b:
        return 'homref' if a == 0 else 'homalt'
    return 'het'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--species', default='')
    ap.add_argument('--output_tsv', required=True)
    ap.add_argument('--min_gq', type=float, default=20.0,
                    help='per-genotype GQ floor; genotypes below this are treated '
                         'as missing (redundant with vcftools --minGQ but defensive)')
    ap.add_argument('--min_dp', type=float, default=26.0,
                    help='per-genotype depth (DP) floor; genotypes below this are '
                         'treated as missing. Guards against low-depth allelic dropout '
                         'miscalling a true het as homozygous (which would fabricate a '
                         'false opposite-homozygote parent site or a false offspring '
                         'Mendelian error). Default 26 matches the de novo mutation '
                         'callability depth.')
    ap.add_argument('--consistency_min', type=float, default=0.9,
                    help='offspring het fraction at opposite-homozygote sites at or '
                         'above which the trio is called consistent')
    ap.add_argument('--min_sites', type=int, default=100,
                    help='minimum opposite-homozygote sites required to issue a verdict')
    args = ap.parse_args()

    vcf = pysam.VariantFile(args.vcf, 'r')
    samples = list(vcf.header.samples)
    trios = build_trios(samples)
    if not trios:
        sys.stderr.write("ERROR: no complete trios found in %s\n" % args.vcf)
        sys.exit(1)

    # Counters keyed by (family, offspring).
    counts = {}
    for t in trios:
        for off in t['offspring']:
            counts[(t['family'], off)] = {
                'n_informative': 0, 'n_het': 0, 'n_homref': 0, 'n_homalt': 0}

    def geno(rec, name):
        """Return GT state for a sample at a record, honouring the GQ and DP floors."""
        smp = rec.samples[name]
        gq = smp.get('GQ')
        if gq is not None and gq < args.min_gq:
            return None
        dp = smp.get('DP')
        if dp is not None and dp < args.min_dp:
            return None
        return gt_state(smp.get('GT'))

    n_records = 0
    for rec in vcf:
        n_records += 1
        for t in trios:
            pf = geno(rec, t['parent_F'])
            pm = geno(rec, t['parent_M'])
            # Opposite homozygotes: one parent homref, the other homalt.
            if not ((pf == 'homref' and pm == 'homalt') or
                    (pf == 'homalt' and pm == 'homref')):
                continue
            for off in t['offspring']:
                state = geno(rec, off)
                if state is None:
                    continue
                c = counts[(t['family'], off)]
                c['n_informative'] += 1
                if state == 'het':
                    c['n_het'] += 1
                elif state == 'homref':
                    c['n_homref'] += 1
                else:
                    c['n_homalt'] += 1
    vcf.close()

    fieldnames = ['species', 'family', 'parent_F', 'parent_M', 'offspring',
                  'n_informative', 'n_het', 'n_homref', 'n_homalt',
                  'consistency', 'error_rate', 'verdict']
    rows = []
    for t in trios:
        for off in t['offspring']:
            c = counts[(t['family'], off)]
            n = c['n_informative']
            if n < args.min_sites:
                consistency = float('nan')
                error_rate = float('nan')
                verdict = 'insufficient-sites'
            else:
                consistency = c['n_het'] / n
                error_rate = (c['n_homref'] + c['n_homalt']) / n
                verdict = ('consistent' if consistency >= args.consistency_min
                           else 'INCONSISTENT')
            rows.append({
                'species': args.species,
                'family': t['family'],
                'parent_F': t['parent_F'],
                'parent_M': t['parent_M'],
                'offspring': off,
                'n_informative': n,
                'n_het': c['n_het'],
                'n_homref': c['n_homref'],
                'n_homalt': c['n_homalt'],
                'consistency': ('%.4f' % consistency) if consistency == consistency else 'NA',
                'error_rate': ('%.4f' % error_rate) if error_rate == error_rate else 'NA',
                'verdict': verdict,
            })

    with open(args.output_tsv, 'w', newline='') as out:
        w = csv.DictWriter(out, fieldnames=fieldnames, delimiter='\t')
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # Console summary.
    n_ok = sum(1 for r in rows if r['verdict'] == 'consistent')
    print("Genotype filters: min_gq=%g, min_dp=%g." % (args.min_gq, args.min_dp))
    print("Scanned %d records; %d trios across %d families."
          % (n_records, len(rows), len(trios)))
    print("Wrote %s" % args.output_tsv)
    print("Consistent trios: %d/%d" % (n_ok, len(rows)))
    flagged = [r for r in rows if r['verdict'] != 'consistent']
    for r in flagged:
        print("  %-18s %-26s n_informative=%s consistency=%s -> %s"
              % (r['family'], r['offspring'], r['n_informative'],
                 r['consistency'], r['verdict']))


if __name__ == '__main__':
    main()
