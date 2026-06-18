# Supplemental Methods: Somatic Mutation Rate Estimation

We estimated per-site somatic single-nucleotide mutation rates for each species from the
set of called somatic single-nucleotide variants (SNVs) and a matched count of callable
trinucleotide opportunities. For every species, rates were computed for the overall SNV
class and within seven mutation classes defined by the substituted base and its immediate
sequence context: C>A, C>G, C>T at non-CpG sites, C>T at CpG sites, T>A, T>C, and T>G.

## Somatic mutation counts (numerator)

Somatic SNVs were taken from the per-species filtered mutation table
(`{species}_final_mutation_table_no_clusters.tsv`), from which clustered calls had already
been removed. Each row of this table corresponds to one somatic SNV observed on a single
sequencing molecule (read). Mutations are reported strand-collapsed to a pyrimidine
reference base through the `oriented_ref`, `oriented_alt`, and `oriented_context` fields:
calls referenced to a purine (A or G) and their flanking trinucleotide context were
reverse-complemented upstream so that every call is expressed with a C or T central base.

A call was retained only if it satisfied all of the following conditions: the oriented
reference base was C or T; the oriented trinucleotide context was three bases long with
its central base equal to the oriented reference; and the oriented alternate allele was a
canonical base (A, C, G, or T) different from the reference. Calls failing any condition
were discarded.

Each retained SNV was assigned to one of six base-substitution types from its oriented
reference and alternate alleles (C>A, C>G, C>T, T>A, T>C, T>G). C>T mutations were further
partitioned by the 3′ flanking base of the oriented context: a 3′ G defines a CpG site and
yields the class C>T_CpG, while any other 3′ base yields C>T_nonCpG. Because each row
represents an independent single-molecule observation, no deduplication was applied.

We denote by `n_c` the number of retained somatic SNVs in class `c`, and by `n_overall`
the total number of retained SNVs across all classes.

## Callable trinucleotide opportunities (denominator)

The denominator is the number of callable trinucleotide contexts surveyed for somatic
mutation. Because somatic mutations are detected on individual sequencing molecules, the
opportunity counts were tallied **directly from the raw PacBio HiFi reads** of each
species rather than from the reference genome, so that the denominator enumerates exactly
the trinucleotide contexts that were sequenced and screened for mutation.

For each species, the HiFi reads were split into 64 patches for parallel processing
(`seqkit split2`). Within every read of length `L`, we trimmed `T = 1500` bp from each end
and counted each overlapping trinucleotide whose central base fell in the retained
interval `[T, L − T)`; a counted trinucleotide window was allowed to draw a single flanking
base from the trimmed boundary, but a trimmed base was never counted as a center.
Trinucleotide windows containing any non-ACGT base were skipped, and reads with `L ≤ 2T`
contributed no windows. Each read therefore contributed

```
w = max(0, L − 2T)
```

trinucleotide opportunities. The 1500 bp end-trim mirrors the read-edge exclusion applied
during somatic variant calling (calls within the trimmed margins are not made), so that
the numerator and denominator accumulate over the same per-read region.

Per-patch counts (64 contexts) were summed across all patches and all reads of a species to
give the raw per-species table of 64 trinucleotide counts (`{species}_trinuc_counts_raw.tsv`).
These were then folded to the 32 pyrimidine-canonical contexts (central base C or T) by
summing each purine-central trinucleotide with its reverse complement — for example, the
canonical context ACA combines the raw counts of ACA and TGT — yielding the final
opportunity table (`{species}_trinuc_counts.tsv`); folding preserves the total count.

Crucially, these opportunities are counted **per molecule** (i.e., per read): the
denominator reflects the total number of single-molecule base observations at risk of
mutation, not a number of diploid genotypes. For each mutation class `c` we formed a
class-conditioned denominator `D_c` from the 32 canonical trinucleotide counts `T(xyz)`,
where `x` and `z` range over the four bases and `y` is the central pyrimidine:

| Class            | Denominator `D_c`                                  |
| ---------------- | -------------------------------------------------- |
| overall          | Σ of all 32 trinucleotide counts                   |
| C>A, C>G         | Σ T(N C N)  — all central-C contexts               |
| C>T_CpG          | Σ T(N C G)  — central-C contexts with a 3′ G       |
| C>T_nonCpG       | Σ T(N C N) − Σ T(N C G)                            |
| T>A, T>C, T>G    | Σ T(N T N)  — all central-T contexts               |

## Rate and confidence interval

The somatic mutation rate for class `c` is the per-base-pair, per-read (per-molecule)
ratio

```
r_c = n_c / D_c
```

No factor of two is applied. Each surveyed trinucleotide context contributes a single
molecular observation, so the denominator already enumerates opportunities at the
single-strand / single-molecule level; there is therefore no diploid (two-haplotype)
correction.

Treating `n_c` as Poisson-distributed with mean `r_c · D_c`, we report exact (Garwood)
Poisson 95% confidence limits obtained from chi-square quantiles and divided by the
denominator:

```
r_low  = χ²(0.025, 2 n_c)       / (2 · D_c)        (r_low = 0 when n_c = 0)
r_high = χ²(0.975, 2 (n_c + 1)) / (2 · D_c)
```

where `χ²(p, k)` is the `p`-quantile of the chi-square distribution with `k` degrees of
freedom.

As a worked example, for one subsocial-baseline species (DUM) the overall denominator is
the sum of all 32 canonical trinucleotide counts, 54,147,318,556 molecular opportunities;
with `n_overall = 249` retained somatic SNVs this gives an overall somatic rate of
`249 / 54,147,318,556 = 4.60 × 10⁻⁹` per base pair per read (exact Poisson 95% CI
`4.05 × 10⁻⁹` to `5.21 × 10⁻⁹`).

## Aggregation across species

Group-level somatic rates (e.g., social, subsocial, and all species combined) were
computed by pooling rather than by averaging per-species rates. Numerators and
denominators were summed across the constituent species `s`,

```
r_c^group = (Σ_s n_{c,s}) / (Σ_s D_{c,s})
```

and the same exact Poisson 95% confidence interval was applied to the pooled count and
pooled denominator.

## Trinucleotide-resolved (SBS96) rates

For mutation rates at the full 96-category (SBS96) resolution, each retained oriented SNV
was assigned to its specific `5′[ref>alt]3′` trinucleotide category, and the rate for each
category was the category count divided by the matched canonical trinucleotide opportunity
count for that context. The same exact Poisson 95% confidence interval was applied to each
category. Group-level SBS96 rates were obtained by the same pooling of category counts and
matched opportunity counts described above.
