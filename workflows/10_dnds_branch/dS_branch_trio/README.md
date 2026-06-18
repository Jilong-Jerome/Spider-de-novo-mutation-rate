# dS Branch Trio Workflow

GWF workflow for branch-specific dN/dS analysis across spider species trios.
For each trio, two species with short reads are mapped to an outgroup/reference
genome, consensus genomes are built, CDS alignments are extracted, and bootstrap
PAML analyses estimate branch and pairwise substitution rates.

## Pipeline Overview

The workflow is defined by `workflow.py`, `workflow_sources/workflow_sources.py`,
and `workflow_sources/workflow_templates.py`.

Shared steps per trio:

```text
01_index
  -> 02_align_{species}
  -> 03a..03f_bam_processing
  -> 04_variant_calling
  -> 05_callable_depth
  -> 05b_filter_snps
  -> 05c_classify_snps
```

Branch-specific steps per trio:

```text
05d_select_snps_{support_policy}_{mode}
  -> 05e_build_mask_{support_policy}_{mode}
  -> 05_consensus_{support_policy}_{mode}
  -> 06_extract_cds
  -> 07_prepare_gene_lists
  -> 08a_sample + 08b_paml for auto_all and auto_bs_1..auto_bs_N
  -> 09_summarize_bootstrap
  -> 10_visualize_pairwise_dS
```

The supported consensus modes are:

- `majority`
- `random`
- `using_ref`
- `using_alt`

Cross-trio combined visualization is also branch-specific.

Each active mode is run under two support policies by default:

- `strict`: heterozygous sites must have clean REF/ALT pileup support.
- `relaxed`: heterozygous sites passing the shared SNP/depth filters are resolved by mode even when pileup support is not clean.

## Consensus Logic

Variant calling is shared by all consensus modes:

```bash
bcftools mpileup -f <reference> -Q <min_base_quality> <bam> |
    bcftools call -c -Oz -o <species>.vcf.gz
```

The base-quality threshold from the config, normally `min_base_quality: 20`,
is therefore the minimum base-quality requirement for all modes.

Before consensus, genome-wide depth is summarized per mapped BAM using
`bedtools genomecov -bga`. The median is calculated only over covered bases
(`DP > 0`), so unmapped reference regions do not force the median to zero.
Callable depth is:

```text
max(min_depth, ceil(0.5 * covered_median_depth)) <= DP <= floor(2.0 * covered_median_depth)
```

SNP records are then filtered once per mapped species with:

```text
TYPE="snp" && QUAL>=20 && DP>=callable_min_depth && DP<=callable_max_depth
```

Read-support cleanliness is then classified once per species. Homozygous calls
that pass the depth and quality filters are accepted from the genotype call.
Heterozygous calls are labelled `het_clean` or `het_unclean` from BAM pileups.
The classification step also writes `{species}_snp_classification_stats.tsv`,
including the number and fraction of heterozygous sites that the strict support
policy would mask. All consensus branches reuse that classification, so failed
branch-specific jobs do not repeat the BAM pileup validation.

Regions from the BAM outside that callable depth range are merged with
branch-specific rejected SNP sites and passed to `bcftools consensus --mask`, so
they become `N` in every mode.

Heterozygous SNPs (`GT=0/1` or `GT=1/0`) are resolved by
`workflow_sources/select_consensus_mode.py`:

- `majority`: choose ALT only when `DP4_alt_forward + DP4_alt_reverse` is greater than `DP4_ref_forward + DP4_ref_reverse`; ties keep REF.
- `random`: choose REF or ALT deterministically from `random_consensus_seed`, species, chromosome, and position.
- `using_ref`: keep the reference/outgroup nucleotide.
- `using_alt`: use the ALT nucleotide.

In the `strict` support policy, a heterozygous site can be resolved by mode only
when pileup support contains only REF and ALT bases, with both represented. In
the `relaxed` support policy, both `het_clean` and `het_unclean` sites are
resolved by mode. Homozygous REF and homozygous ALT calls are not rejected only
because a few reads support another nucleotide: those reads can come from
sequencing error, residual low-quality observations, local alignment noise, or
genotype-likelihood uncertainty.

Homozygous ALT calls remain ALT in all branches. Sites selected as REF are omitted
from the branch-specific VCF passed to `bcftools consensus`, so the reference
base is retained. Sites rejected as unclean are written to a rejected site BED
and included in the `bcftools consensus --mask` BED.

## Configuration

Each trio has one YAML file in `configurations/`. All `*config.y*ml` files are
loaded automatically by `workflow.py`.

Important fields:

```yaml
trio_name: SAR_BIC_TEN
reference_species: TEN
mapping_species:
  SAR:
    read1: data/fastq/SAR/SAR_family1_F_female_R1.fq.gz
    read2: data/fastq/SAR/SAR_family1_F_female_R2.fq.gz
  BIC:
    read1: data/fastq/BIC/BIC_R1.fq.gz
    read2: data/fastq/BIC/BIC_R2.fq.gz
trio_tree: "((SAR,BIC),TEN)"

min_mapping_quality: 60
min_base_quality: 20
min_depth: 5
coverage_median_min_factor: 0.5
coverage_median_max_factor: 2.0
min_coverage_fraction: 0.80

# Supported consensus modes: majority, random, using_ref, using_alt
# To enable all modes, use:
# consensus_modes: [majority, random, using_ref, using_alt]
consensus_modes: [random]
consensus_support_policies: [strict, relaxed]
random_consensus_seed: 20260514
n_bootstrap: 500
```

The default active consensus mode is `random`. The other modes remain supported
and can be re-enabled by editing `consensus_modes`. To run only a subset of
support policies, edit `consensus_support_policies`.

## Output Layout

Shared outputs stay directly under each trio directory:

```text
steps/{TRIO}/
  01_index/
  02_mapping/
  03_bam_processing/
  04_variant_calling/
  05_callable_depth/   # depth mask, filtered SNP VCF, shared SNP classification TSV
```

Branch-specific outputs are nested under `{SUPPORT_POLICY}_{MODE}`:

```text
steps/{TRIO}/{SUPPORT_POLICY}_{MODE}/
  05_consensus/        # selected SNP VCF, rejected-site BED, mask BED, consensus FASTA
  06_cds_per_gene/
  07_gene_lists/
  08_bootstrap/
  09_summary/
  10_visualization/
```

Branch-specific logs and DONE sentinels follow the same pattern:

```text
logs/{TRIO}/{SUPPORT_POLICY}_{MODE}/
```

Combined cross-trio figures are written per branch:

```text
steps/_combined/{SUPPORT_POLICY}_{MODE}/combined_pairwise_dS.{png,pdf}
logs/_combined/{SUPPORT_POLICY}_{MODE}/11_combine_trio_visualizations.DONE
```

## Running

```bash
conda activate gwf_new
gwf status
gwf run
```

Run or inspect a specific branch target, for example:

```bash
gwf status SAR_BIC_TEN_strict_random_bcftools_consensus_SAR
gwf run SAR_BIC_TEN_relaxed_random_bcftools_consensus_SAR
gwf run SAR_BIC_TEN_strict_random_paml_auto_bs_323
gwf run combine_trio_visualizations_strict_random
```

If Slurm accounting is unavailable, a local backend dry run can still be useful
for checking target names and DAG construction:

```bash
gwf -b local run --dry-run SAR_BIC_TEN_strict_random_bcftools_consensus_SAR
```

## Light Checks

```bash
python3 -m py_compile workflow.py workflow_sources/*.py
```

The mode selector can be tested directly with a small VCF and matching
classification TSV:

```bash
python3 workflow_sources/select_consensus_mode.py \
    --mode random \
    --support-policy strict \
    --species SAR \
    --seed 20260514 \
    --classification-tsv SAR_snp_classification.tsv \
    --rejected-bed rejected_sites.bed \
    < input.vcf
```
