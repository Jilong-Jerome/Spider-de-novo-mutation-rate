# Methods — Validation of parent–offspring trios by parental opposite-homozygote segregation

## Rationale

Frequency-based relatedness estimators (KING robust kinship; PLINK IBD/PI_HAT) assume
that the founders of the sampled families are drawn from a large, outbred population, so
that population allele frequencies provide a meaningful baseline for expected allele
sharing. This assumption is violated in the social spider system, where colonies are
highly inbred and the founding parents are themselves closely related. Applied to our
data these estimators returned biologically impossible values for known parent–offspring
pairs (robust kinship coefficients ranging from −0.5 to below −13), and could not be used
to confirm the labelled pedigrees. We therefore validated each labelled trio with a
direct, allele-frequency-free Mendelian test that depends only on within-trio segregation
and is consequently unaffected by inbreeding or by the choice of an allele-frequency
reference.

## Principle

For a putative trio consisting of the two parents of a family (the female F and the male
M) and one of their offspring (S), we restrict attention to *parental opposite-homozygote
sites* — biallelic SNPs at which the two parents are fixed for opposite alleles (one
parent homozygous reference, 0/0; the other homozygous alternate, 1/1). At any such site
Mendelian transmission is fully determined: a genuine offspring of *both* parents must
inherit the reference allele from the homozygous-reference parent and the alternate allele
from the homozygous-alternate parent, and is therefore obligately heterozygous (0/1). We
define the trio *consistency* as the fraction of parental opposite-homozygote sites at
which the offspring is observed to be heterozygous. A true trio yields consistency ≈ 1
(departing from unity only through residual genotyping error), whereas an incorrectly
assigned parent produces a sharp drop in consistency. Because the informative sites and
the expected offspring genotype are defined entirely by the two parents' genotypes, no
population allele frequency enters the calculation. Note that this test confirms that
*both* assigned parents are the true parents of the offspring; it is not intended to
distinguish full-sib from parent–offspring relationships, which is not required here.

## Variant input

Genotypes were taken from the per-species GATK joint call set. Sites were restricted to
biallelic SNPs that passed the GATK variant filters (PASS) and lay on autosomes;
sex-chromosome contigs were excluded throughout because hemizygosity in males breaks the
diploid Mendelian expectation. No minor-allele-frequency filter was applied, since the
informative-site criterion is the per-trio opposite-homozygote condition rather than a
population-frequency threshold; consequently low-frequency variants that segregate within
a single family — precisely the most informative sites for this test — are retained.
Genotypes with genotype quality below GQ = 20 were set to missing.

## Per-genotype depth filter

A central source of false signal in this test is allelic dropout at low sequencing depth:
a truly heterozygous genotype that is undersampled can be miscalled as homozygous. Such a
miscall is doubly damaging — at a parent it can fabricate a spurious opposite-homozygote
(informative) site, and at an offspring it can convert a genuinely heterozygous,
Mendelian-consistent genotype into an apparent Mendelian error — in both cases inflating
apparent inconsistency. To guard against this, every genotype used in the test (both
parents and the offspring) was required to have a read depth of at least DP = 26, the same
per-individual callability depth used for the de novo mutation rate estimation in this
study; genotypes below this depth were treated as missing and the site was not counted for
the affected trio. Read depth was evaluated per genotype from the VCF FORMAT field.

## Scoring and verdicts

Trios were enumerated from the structured sample identifiers
(`{species}_family{N}_{F|M|S#}`); each family contributing both parents and at least one
offspring yielded one trio per offspring. The call set was traversed once, and for every
trio at every site at which both parents were called as opposite homozygotes (passing the
GQ and DP filters), the offspring genotype (likewise passing the filters) was classified
as heterozygous (consistent) or homozygous (a Mendelian error). For each offspring we
recorded the number of informative (parental opposite-homozygote) sites, the number of
heterozygous calls, the consistency (heterozygous fraction), and the complementary error
rate. A trio was classified as **consistent** when it was supported by at least 100
informative sites and had consistency ≥ 0.90; as **inconsistent** when it had at least
100 informative sites but consistency below 0.90; and as **insufficient-sites** when fewer
than 100 parental opposite-homozygote sites survived filtering, in which case no verdict on
the pedigree could be issued. The latter category typically reflects parent pairs that are
too similar (too few fixed differences) to certify the offspring rather than evidence
against the assigned relationship.

---

Parameter values (GQ ≥ 20, DP ≥ 26, consistency ≥ 0.90, ≥ 100 informative sites) are the
configured defaults used for all seven species and are exposed as `trio_min_gq`,
`trio_min_dp`, `trio_consistency_min`, and `trio_min_sites` in the workflow
configuration. The test is implemented in `workflow_sources/trio_mendel_check.py`.
