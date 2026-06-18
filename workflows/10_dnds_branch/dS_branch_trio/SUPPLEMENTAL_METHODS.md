# Supplemental Methods: Branch-specific dN/dS estimation in spider species trios

## Overview

To estimate lineage-specific synonymous (dS) and non-synonymous (dN)
substitution rates we analysed two species trios, each consisting of two query
species and the *B. bicolor* (BIC) genome as the third taxon and shared
mapping reference: ((SAR, PAC), BIC) and ((DUM, TEN), BIC). For each trio,
short reads from the two query species were mapped to the BIC reference
assembly, genotyped, and used to build masked branch-specific consensus
genomes. Per-gene coding-sequence (CDS) alignments were extracted using the BIC
annotation and concatenated, and branch-specific and pairwise dN/dS were
estimated with PAML `codeml`. Confidence intervals were obtained by gene-level
bootstrap resampling. The complete pipeline was implemented as a reproducible
Grid Workflow (`gwf`) DAG run on a SLURM cluster; all steps, parameters, and
helper scripts are available in the project repository.

## Read mapping and alignment processing

Paired-end reads from each query species were aligned to the BIC reference
genome with `bwa-mem2 mem` (v2.2.1) using default parameters. The remaining
processing was performed with samtools v1.17 and Picard v2.25.4. Read groups
were added with Picard `AddOrReplaceReadGroups`. Mate information was corrected
and secondary/unmapped reads removed with `samtools fixmate -rm`; alignments
were then coordinate-sorted (`samtools sort`) and PCR/optical duplicates were
identified and removed with `samtools markdup -r`. Finally, alignments were
filtered with `samtools view` to retain only properly paired, mapped reads
(SAM flags `-f 0x2 -F 0x4`) with a mapping quality of at least 60
(`min_mapping_quality = 60`), and the resulting BAM files were indexed. All
mapping and BAM-processing steps were performed once per query species and
shared across all downstream consensus branches.

## Variant calling

Variants were called per species from the final filtered BAM with `bcftools`
v1.16. Genotype likelihoods were computed with `bcftools mpileup` using the BIC
reference and a minimum base quality of 20 (`-Q 20`), and genotypes were called
with the `bcftools call` consensus caller (`-c`). The minimum base quality of
20 (`min_base_quality`) was applied consistently here and in all downstream
pileup-based steps.

## Callable-depth determination

To exclude regions with anomalous read depth, per-base coverage was computed
from the final BAM with `bedtools genomecov -bga` (bedtools v2.31.1). For each
species we computed the coverage-weighted median depth over all covered sites
(depth > 0) and derived a callable-depth window as a multiple of that median:
the lower bound was `max(min_depth, ceil(0.5 x median))` and the upper bound
was `floor(2.0 x median)`, with an absolute depth floor of `min_depth = 5`.
Genomic intervals with depth outside this callable window were merged into a
per-species "uncallable" BED mask. Coverage statistics (median, mean,
percentiles, covered fraction) and a per-depth histogram were recorded for
quality-control diagnostics.

## SNP filtering and heterozygous-site classification

Single-nucleotide polymorphisms were extracted from each per-species VCF with
`bcftools view`, retaining only biallelic SNP records with a variant quality of
at least 20 (`QUAL >= 20`) and a site depth (`DP`) within the species-specific
callable-depth window defined above.

Each retained SNP was then classified once per species using the read pileup at
that site (pysam v0.23.3, minimum base quality 20). Homozygous-reference and
homozygous-alternate genotypes were accepted directly from the genotype call.
Heterozygous genotypes were sub-classified by their pileup support: a
heterozygous site was labelled `het_clean` only if the observed high-quality
bases consisted exclusively of the reference and alternate alleles, with both
alleles supported by at least one read; otherwise (a third allele observed, or
either expected allele unsupported) it was labelled `het_unclean`. This
classification was computed once and reused across all consensus branches to
avoid repeated pileup traversal.

## Branch-specific consensus construction

Because branch-specific substitution counts depend on how heterozygous sites
are resolved, consensus genomes were built under explicit, factorial
"support-policy x resolution-mode" branches:

- **Support policy** controls which heterozygous sites are eligible for
  resolution. Under the `strict` policy only `het_clean` sites are resolved and
  `het_unclean` sites are masked; under the `relaxed` policy both `het_clean`
  and `het_unclean` sites are resolved.
- **Resolution mode** controls which allele an eligible heterozygous site
  contributes. The mode used for the reported results was `random`, in which
  the allele is chosen deterministically by hashing a fixed seed
  (`random_consensus_seed = 20260514`) with the species name and site
  coordinates, so the choice is reproducible and independent across sites. The
  pipeline also supports `majority` (allele with greater DP4 read support),
  `using_ref`, and `using_alt` modes.

For each branch, homozygous-alternate sites and resolved heterozygous sites
assigned the alternate allele were written to a selected-SNP VCF;
homozygous-reference sites, sites resolved to the reference allele, and
non-resolvable sites were not. Sites that could not be confidently resolved
under the active policy/mode were collected into a rejected-sites BED. The
final per-species mask was the union of the shared uncallable-depth BED and the
branch-specific rejected-sites BED. The branch-specific consensus genome was
generated with `bcftools consensus`, applying the selected-SNP VCF to the BIC
reference and replacing all masked positions with `N`. The BIC reference itself
served as the third (unmodified) genome of each trio.

## Per-gene CDS extraction and quality control

Coding sequences were extracted per gene from the trio of genomes (the BIC
reference and the two query-species consensus genomes) using the BIC GFF3
annotation. For each gene the transcript with the longest total CDS was
selected. CDS exon coordinates were used to extract and concatenate the
corresponding sequence from each of the three genomes; genes on the minus
strand were reverse-complemented. Because all three sequences are defined
against the same BIC coordinate system, the resulting three-sequence alignment
is positional and gap-free by construction. Any non-ACGT character was
converted to `N`. Each gene was retained ("passing") only if (i) its CDS length
was a multiple of 3, (ii) the BIC reference CDS contained no premature stop
codon, and (iii) every one of the three sequences had a callable
(non-`N`) fraction of at least 0.80 (`min_coverage_fraction`). Passing genes
were partitioned into autosomal genes (chromosomes `bic_1`-`bic_14`) and
X-linked genes (chromosomes prefixed `bic_X`); the analyses reported here use
the autosomal gene set.

## Codon alignment filtering and concatenation

For each analysis the per-gene CDS alignments of the passing autosomal genes
were concatenated into a single supermatrix. The supermatrix was then filtered
to fully callable codons: a codon was retained only if all three positions in
all three species were unambiguous bases (A/C/G/T); any codon containing a gap
or ambiguous base in any species was discarded. The filtered supermatrix was
written in PAML sequential PHYLIP format.

## dN/dS estimation with PAML

Branch-specific and pairwise substitution rates were estimated from the
filtered codon supermatrix with PAML `codeml` (v4.9). Two models were run for
each alignment:

- **Free-ratio branch model:** `runmode = 0`, `model = 1` (one independent
  dN/dS ratio per branch), with the trio topology supplied as a fixed
  user tree (`((SAR,PAC),BIC)` or `((DUM,TEN),BIC)`). This yields
  branch-specific dN, dS, and dN/dS for each terminal branch and the internal
  branch.
- **Pairwise model:** `runmode = -2`, `model = 0`, yielding maximum-likelihood
  pairwise dN, dS, and dN/dS for all species pairs.

Both models used codon sequences (`seqtype = 1`), the F3x4 codon-frequency
model (`CodonFreq = 2`), the universal genetic code (`icode = 0`), and
estimated the transition/transversion ratio kappa (`fix_kappa = 0`, initial
`kappa = 2`) and omega (`fix_omega = 0`, initial `omega = 1`). Sites with any
ambiguity were removed by PAML's `cleandata = 1` option, complementing the
upstream codon filter. dN/dS was computed as dN/dS where dS > 0. Branch-level
results were parsed from the `codeml` output with custom scripts.

## Bootstrap confidence intervals

Confidence intervals on dN, dS, and dN/dS were obtained by gene-level
non-parametric bootstrap. The point estimate was computed from the full set of
passing autosomal genes. In addition, 500 bootstrap replicates were generated,
each by resampling genes with replacement (`shuf -r`) to the same number of
genes as the full set. For every replicate the resampled genes were
concatenated, codon-filtered, and analysed through the identical `codeml`
branch and pairwise models. For each branch and each species pair, the
bootstrap distribution across the 500 replicates was summarised by its mean,
median, and the 2.5th and 97.5th percentiles, the latter two defining the 95%
confidence interval.

## Reproducibility and visualisation

The entire analysis was implemented as a `gwf` workflow on a SLURM cluster, with
each step pinned to a dedicated conda environment (bwa-mem2 2.2.1, samtools
1.17, bcftools 1.16, bedtools 2.31.1, PAML 4.9, and Python 3.10 with pysam
0.23.3 and Biopython 1.85). Each trio was processed by an independent
configuration file, and species names and tree topologies were taken from those
configurations rather than hard-coded, so the pipeline is trio-agnostic.
Per-trio results were visualised as pairwise-dS heatmaps and unrooted dS trees;
a combined multi-panel figure was produced across trios using a shared
heatmap colour scale and per-trio tree normalisation.
