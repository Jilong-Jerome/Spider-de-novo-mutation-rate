# coverage_report

GWF workflow computing sequencing coverage **per individual × chromosome × species**,
plus a **per-individual autosomal DP** report and a combined cross-species plot.

Species enabled (one config each): **AFR, BIC, DUM, LIN, MIM, SAR, TEN**.

## Inputs (via `../data/`)
- `data/{SP}/{sample}_final.bam` (+ `.bai`) — dup-removed, coordinate-sorted BAMs.
- `data/{SP}_GATK.vcf.gz` (+ `.tbi`) — **all-sites** GATK VCF, FORMAT has `DP`.

The sample set is read from the VCF (`bcftools query -l`) at graph-build time, so
the report covers only individuals present in the VCF (any BAM-only individuals,
e.g. LIN family5, are intentionally dropped).

## Metrics (per individual × chromosome)
1. `mean_all_depth` = Σ sequenced nucleotides on the chrom ÷ chrom length (from BAM).
2. `mean_covered`, `median_covered` = mean & median depth over sites with depth ≥ 1 (BAM).
3. From the VCF DP: `mean_DP_all`/`median_DP_all` (incl. DP=0 sites) and
   `mean_DP_ge1`/`median_DP_ge1` (only DP ≥ 1 sites).

## Per-individual autosomal DP report
Autosomes = chromosomes not matching `{sp}_X{N}`. From the VCF DP over **all
autosomal sites (DP=0 included)**, `{SP}_autosome_dp` writes
`steps/{SP}/{SP}_autosome_dp_per_individual.tsv` (`sample, n_sites, mean_DP_all,
median_DP_all`) by pooling the per-chromosome DP histograms (`vcf_dp_stats.py`
emits `{chrom}.vcf_dp_hist.tsv` as a side-output; the pooled histogram gives an
exact median). The single top-level `coverage_combined_dp_plot` then reads all
species' TSVs → `steps/combined/all_species_autosome_dp_per_individual.tsv` and a
barplot `..._barplot.png` (one bar per individual = mean DP, coloured by species,
median dot, horizontal grand-mean line across all individuals of all species).

## Structure
- `workflow.py` — globs `configurations/*.config.yaml`, builds per-species graphs,
  then registers the one cross-species `coverage_combined_dp_plot` target.
- `configurations/{SP}.config.yaml` — paths, conda envs, chromosome list.
- `workflow_sources/workflow_sources.py` — builder (reads VCF samples, derives
  autosomes, registers targets; returns metadata for the combined step).
- `workflow_sources/workflow_templates.py` — templates: contigs, bam_depth,
  vcf_dp, merge, autosome_dp, combined_plot.
- `workflow_sources/{bam_depth_stats,vcf_dp_stats,merge_report,autosome_dp_report,combined_dp_plot}.py`
  — analysis scripts.
- `steps/{SP}/` — per-species outputs; `steps/combined/` — cross-species outputs.
- `logs/` — `.DONE` sentinels.

## Targets per species
`{SP}_contigs` → `{SP}_bamcov_{sample}_{chrom}` (one per individual × chromosome;
each scans a single chromosome via the `.bai` index for parallelism) +
`{SP}_vcfdp_{chrom}` (one per chromosome) → `{SP}_merge` (per-chrom table) and
`{SP}_autosome_dp` (per-individual autosomal TSV). Across species: one
`coverage_combined_dp_plot`.

## Run
```bash
cd coverage_report
conda run -n gwf_new gwf status
conda run -n gwf_new gwf run
```

## Add a species
Symlink its data under `../data/`, copy `configurations/LIN.config.yaml` to
`{SP}.config.yaml`, and update `species`, paths, and the `chromosomes` list
(from the BAM `@SQ` / VCF `##contig` header).

## Conda envs
`samtools117` (samtools), `bcftools` (bcftools), `python_phylo` (pandas), `gwf_new` (gwf).
