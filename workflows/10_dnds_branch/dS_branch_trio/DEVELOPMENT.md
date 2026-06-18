# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Does

Branch-specific dN/dS analysis for trios of spider species using PAML codeml free-ratio model. Maps two species' reads to a configured reference/outgroup genome, calls variants, builds consensus sequence sets for each support-policy/mode branch, extracts per-gene CDS alignments, then runs bootstrap PAML to estimate branch-specific substitution rates with confidence intervals. Supports multiple trios via separate config files.

## Running the Workflow

```bash
conda activate gwf_new
gwf status          # check job states
gwf run             # submit pending jobs to SLURM
gwf logs <jobname>  # view SLURM logs for a job
```

The workflow runs on a SLURM cluster (account: `spider2`). All jobs are defined via GWF (Grid WorkFlow) and submitted automatically. There are no unit tests.

## Architecture

**Entry point**: `workflow.py` — loads config YAML files from `configurations/`, calls `trio_dnds_workflow()` for each.

**Three-file workflow pattern**:
- `workflow.py` — GWF workflow instantiation, discovers configs
- `workflow_sources/workflow_sources.py` — DAG orchestration, wires templates together using config values
- `workflow_sources/workflow_templates.py` — individual SLURM job templates (each returns a `gwf.AnonymousTarget`)

**Pipeline DAG**:
```
Shared per trio: Mapping + variant calling
  01_index → 02_align_{sp1,sp2} → 03a-f_bam_processing → 04_variant_call
  05a_callable_depth_{sp} → 05a_coverage_distribution (per-trio diagnostic leaf)

Per consensus support-policy/mode branch: Consensus + bootstrap PAML + visualisation
  05_consensus → 06_extract_cds → 07_prepare_gene_lists
      → [08a_sample + 08b_paml] × 501 replicates
      → 09_summarize → 10_visualize_pairwise_dS

Cross-trio, per consensus branch:
  11_combine_trio_visualizations_{support_policy}_{mode}
```

**Job dependency tracking**: GWF uses `.DONE` sentinel files in `logs/` as inputs/outputs to chain job dependencies.

**Helper scripts** in `workflow_sources/`:
| Script | Purpose |
|--------|---------|
| `extract_cds_per_gene.py` | Extract per-gene CDS from trio consensus FASTAs |
| `prepare_gene_lists.py` | Split genes into configured autosomal vs X-linked chromosome lists |
| `concat_genes.py` | Concatenate per-gene FASTAs into single alignment |
| `concat_codon_filter_paml.py` | Keep fully-callable codons (A/C/G/T in all species) and write PAML phylip |
| `paml2tab.py` | Parse codeml output into TSV |
| `parse_pairwise.py` | Extract pairwise dN/dS from codeml results |
| `branch_dnds.py` | Extract branch-specific dN/dS |
| `summarize_bootstrap.py` | Aggregate bootstrap replicates, compute CIs |
| `visualize_pairwise_dS.py` | Heatmap, unrooted-tree plot, and recalculated branch-length TSV from pairwise dS |
| `combine_trio_visualizations.py` | Multi-trio 2 x N panel figure (heatmap row + unrooted-tree row); shared heatmap colour scale, per-trio tree normalisation |
| `callable_depth_regions.py` | Compute callable-depth bounds from genomecov; write coverage stats (median/mean/percentiles), uncallable BED, and per-depth distribution TSV |
| `plot_coverage_distribution.py` | Per-trio diagnostic: covered-site depth distribution plots + mean-vs-median summary TSV |
| `classify_consensus_snps.py` | Classify SNP read support once per species using BAM pileups and write heterozygous strict-filter stats |
| `select_consensus_mode.py` | Resolve shared SNP classifications for strict/relaxed support policies and the configured consensus modes |

## Configuration

Config files live in `configurations/`. Each `*config.yaml` defines one trio. All configs found by `workflow.py` are processed automatically. Current configs:
- `PAC_SAR_trio_dnds.config.yaml` — SAR_PAC_BIC trio
- `TEN_DUM_trio_dnds.config.yaml` — DUM_TEN_BIC trio

Key parameters:
- `trio_name`: namespaces all output dirs and job names (e.g., `SAR_PAC_BIC`, `DUM_TEN_BIC`)
- `mapping_species`: dict of species → read paths to map against the configured reference genome (species names are used dynamically throughout the pipeline)
- `trio_tree`: Newick tree string for PAML (e.g., `"((SAR,PAC),BIC)"`, `"((DUM,TEN),BIC)"`)
- `consensus_modes`: consensus mode branches to run; default active mode is `random`; supported modes are `majority`, `random`, `using_ref`, and `using_alt`
- `consensus_support_policies`: support-policy branches to run; defaults are `strict`, `relaxed`
- `random_consensus_seed`: deterministic seed for the `random` heterozygous-site mode
- `n_bootstrap`: number of bootstrap replicates (default 500)
- Quality filters: `min_mapping_quality: 60`, `min_base_quality: 20`, `min_depth: 5`, `min_coverage_fraction: 0.80`

## Conda Environments

Jobs use multiple conda environments activated within SLURM scripts:
- `bwa2` — bwa-mem2 for indexing/mapping
- `samtools117` — samtools for BAM processing
- `bcftools` — variant calling and consensus
- `python_phylo` — Python scripts (BioPython, etc.)
- `bedtools` — coverage masking
- `paml` — codeml for dN/dS estimation
- `gwf_new` — for running `gwf` commands locally

## Key Design Decisions

- **Multi-trio support**: Species names are derived from config `mapping_species` keys; tree topology from `trio_tree`. No species names are hardcoded in the pipeline code.
- **Bootstrap at runtime**: The `concat_filter_paml.sh` wrapper reads gene ID files at runtime (not at GWF definition time), which is critical because bootstrap ID lists don't exist when the workflow DAG is defined. The tree is also passed as a runtime argument.
- **Consensus strategy**: Steps 01-04, Step 05 callable-depth masking, SNP filtering, and SNP classification are shared per species. Homozygous calls passing the depth/quality filters are accepted from the genotype call; heterozygous calls are classified as `het_clean` or `het_unclean` by BAM pileup support before branch-specific selection. `strict` branches mask `het_unclean`; `relaxed` branches resolve both heterozygous classes by mode. Branch-specific Step 05 jobs only select SNPs, build the final mask, and run `bcftools consensus`; this avoids repeating BAM pileups across branches and preserves checkpoints after timeouts. Low-coverage regions are masked with `N` in all branches.
- **Base quality**: `min_base_quality` is passed to `bcftools mpileup -Q` and is therefore the shared minimum base-quality threshold for all consensus modes.
- **PAML model**: Free-ratio (`model=1`) with `cleandata=1, CodonFreq=2`. Tree topology defined per trio in config.
- **Bootstrap sampling**: `shuf -r` (with replacement), gene count determined at runtime from the passing genes list.

## Output Structure

Shared outputs under `steps/{TRIO_NAME}/`:
```
01_index/          — bwa-mem2 index files
02_mapping/        — (intermediate, cleaned up)
03_bam_processing/ — final filtered BAMs
04_variant_calling/ — VCF files per species
05_callable_depth/ — per-species covered-depth stats, depth-distribution TSV, shared depth mask, filtered SNP VCF, SNP classification TSV, and classification stats TSV
  coverage_distribution/ — per-trio coverage_distribution.{png,pdf} + coverage_distribution_summary.tsv (diagnostic)
```

Branch-specific outputs under `steps/{TRIO_NAME}/{SUPPORT_POLICY}_{MODE}/`:
```
05_consensus/      — selected SNP VCF, rejected-site BED, final mask BED, and consensus FASTA per species
06_cds_per_gene/   — per_gene_fasta/{gene_id}.fa + passing_genes.txt
07_gene_lists/     — auto_passing_genes.txt, x_passing_genes.txt
08_bootstrap/      — auto_all/ + auto_bs_{1..500}/ replicate dirs
09_summary/        — branch/ and pairwise/ summary TSVs
10_visualization/  — pairwise dS heatmap, unrooted-tree plot, recalculated branch-length TSV
```

Branch-specific cross-trio outputs:
```
_combined/{SUPPORT_POLICY}_{MODE}/
  combined_pairwise_dS.{png,pdf}  — 2 x N multipanel: heatmap row + unrooted-tree row, one column per trio
```
