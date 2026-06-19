# Supplementary Methods: Sibling-shared mutations and parental bias in relation to sociality

We used the autosomal germline *de novo* mutation (DNM) call set generated for the
paired social and subsocial *Stegodyphus* trios to address two questions about how
mutations are distributed among full siblings. First, we asked whether the
**fraction of mutation loci that are shared between siblings differs with
sociality** (social vs. subsocial). Second, within the subsocial species that
yielded phaseable shared mutations, we asked whether **sibling-shared and unique
mutations differ in their parental (paternal vs. maternal) origin**. Throughout,
species were grouped by social organisation following the study's species tree:
the **social** species are *S. dumicola* (DUM), *S. mimosarum* (MIM) and
*S. sarasinorum* (SAR), and the **subsocial** species are *S. africanus* (AFR),
*S. bicolor* (BIC), *S. lineatus* (LIN) and *S. tentoriicola* (TEN). Each social
species is the sister of a subsocial species in the working topology
`((((DUM,TEN),(SAR,BIC)),(MIM,AFR)),LIN)`.

## Data and locus definition

Both analyses started from the per-offspring autosomal DNM table
(`autosomal_DNM_spread_sheet.tsv`), in which each row is a candidate DNM observed
in one offspring of a parent–offspring trio, together with its parental phase
assignment (`Phased` ∈ {P = paternal, M = maternal, U = unphased, F =
unphaseable/failed}).

Trios flagged as unreliable were removed before any counting: BIC_family2 and
BIC_family3 in both analyses, and SAR_family3 additionally in the cross-species
sharing comparison, matching the exclusion lists in the analysis configuration.

A **mutation locus** was defined by the combination of family and genomic change,
i.e. `(family, chromosome, position, reference allele, alternate allele)`.
Multiple offspring rows carrying the same locus within a family were collapsed to
a single locus before classification. A locus observed in **two or more offspring**
of the same family was classified as **sibling-shared**; a locus observed in a
single offspring was classified as **unique**. Because shared loci can be observed
in several offspring, a single parental phase was resolved per locus from the set
of offspring-level phase calls using the following priority rule: if all calls
agreed, that phase was kept; a single parental label mixed with U was resolved to
that parent (P+U → P; M+U → M); loci containing both P and M ("conflicts") were
flagged and excluded from phase-based analyses; loci with only U were treated as
unphased (U); and F (unphaseable) calls were dropped before re-evaluating the
remaining labels. No P/M conflicts were observed in the analysed data.

## Test 1 — Fraction of sibling-shared mutations vs. sociality

To test whether sibling sharing of DNMs is associated with sociality, we
classified every retained mutation locus as unique or sibling-shared as defined
above, after removing excluded trios and dropping rows with an unphaseable
(`F`) call. Loci were then pooled within each sociality group across the
constituent species, giving the following 2 × 2 contingency table:

| Sociality | Sibling-shared loci | Unique loci | Total |
|-----------|--------------------:|------------:|------:|
| Subsocial (AFR, BIC, LIN, TEN) | 27 | 452 | 479 |
| Social (DUM, MIM, SAR)         |  0 | 198 | 198 |

We assessed the association between sociality and the shared/unique
classification with a two-sided **Fisher's exact test**
(`scipy.stats.fisher_exact`). Sibling-shared loci were observed only in subsocial
species (27/479 = 5.6% of subsocial loci) and were entirely absent from social
species (0/198); the excess of sharing in subsocial species was significant
(p = 1.1 × 10⁻⁴). Because the social shared cell is zero the uncorrected odds
ratio is undefined (∞), so we report a pseudocount-corrected odds ratio
(1 added to each cell), OR = 12.3; the pseudocount affects only the point
estimate, not the Fisher p-value.

We pooled loci within each sociality group rather than fitting a per-species or
mixed-effects model because all three social species contained zero sibling-shared
loci, which precludes per-species rate estimation and model-based contrasts. The
paired social/subsocial structure of the species tree (each social species sister
to a subsocial species) means the contrast is not driven by a single divergent
lineage.

## Test 2 — Parental bias of sibling-shared vs. unique mutations

This analysis was restricted to the three subsocial species that yielded phased
sibling-shared mutations — BIC, LIN and TEN (BIC_family2 and BIC_family3
excluded); AFR was not included because it contributed too few phased shared loci.
Mutation loci were collapsed and classified as unique or sibling-shared, and each
locus was assigned a single parental phase using the resolution rule above. For
the parental-bias contrast we retained only loci with a resolved parental origin
(P or M); unphased (U), unphaseable (F) and conflicting loci were excluded.

Pooling the three species, parental origin was cross-tabulated against
sharing class:

| Mutation class | Paternal (P) | Maternal (M) | P/M ratio |
|----------------|-------------:|-------------:|----------:|
| Unique         | 106 | 75 | 1.41 |
| Sibling-shared |   5 | 13 | 0.38 |

The difference in parental bias between unique and sibling-shared loci was tested
with a two-sided **Fisher's exact test** on this 2 × 2 table. Unique mutations
were paternally biased (P/M = 1.41) whereas sibling-shared mutations were
maternally biased (P/M = 0.38); the contrast was significant (odds
ratio = 3.67, p = 0.023). The same direction of effect was seen in each species
individually (per-species P/M ratios: unique 1.67, 1.34, 1.37 vs. shared 0.38,
0.50, 0.00 for BIC, LIN, TEN respectively).

The maternal enrichment of sibling-shared mutations is consistent with these loci
arising early — in the maternal germline or during early embryogenesis — such that
the variant is transmitted to, or already present across, multiple offspring of a
clutch. Such shared, maternally biased mutations are detectable in the subsocial
species, whose sampled offspring are full-sib clutches, but not in the social
species (Test 1).

## Software and reproducibility

Both tests were implemented in Python (pandas; `scipy.stats.fisher_exact`,
two-sided). The shared-fraction test (Test 1) is computed by
`scripts/shared_fraction_sociality.py`, which reads the DNM table directly, drops
excluded trios and `F` calls, collapses loci, groups species by sociality, and
writes the per-species counts, the pooled 2×2 table, the Fisher result and the
pseudocount-corrected (+1) odds ratio to
`results/shared_fraction_sociality_summary.txt`; its per-species unique/shared
counts match the DNM_alpha summariser `summarize_shared_mutations.py` (run under
`DNM_alpha/configurations/dnm_alpha.config.yaml`). The parental-bias test
(Test 2) — locus collapsing, phase resolution and the Fisher's exact test — is
implemented in `scripts/paternal_bias_analysis.py`, with the reported counts,
ratios, odds ratio and p-value written to
`results/paternal_bias_unique_shared_summary.txt`. Statistical significance was
assessed at α = 0.05. Note one deliberate convention difference between the tests:
the cross-species sharing summary (Test 1) drops `F` calls before counting loci,
whereas the parental-bias script (Test 2) retains `F` in the locus inventory but
excludes it, together with U and conflicts, from the parental-origin contrast.
