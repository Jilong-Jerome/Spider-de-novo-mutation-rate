# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository is part of a study on de novo mutation rates in social spiders ("Convergent Evolution of Sociality Causes Reduction of Mutation Rates in Spiders"). It contains three analysis components within the parent `science_resub/` directory:

- **`somatic/`** (this directory) - Filtered somatic mutation data for 7 spider species (AFR, BIC, DUM, LIN, MIM, SAR, TEN). Contains only symlinked TSV files pointing to variant calling results in `/home/jilong/mutationalscanning/`.
- **`ind_div_repair_gene/`** - GWF pipeline analyzing DNA repair gene heterozygosity and divergence (Dxy) across spider families for species AFR, MIM, TEN.
- **`testis_verify/`** - GWF pipeline verifying RNA-seq tissue identity by comparing testis-specific gene expression between ovary (control) and candidate testis samples in MIM.

## Running Workflows

Both pipelines use [GWF (Genome Workflow)](https://gwf.app/) on a SLURM HPC cluster (account: `spider2`).

```bash
# Activate GWF environment
conda activate gwf_new

# Check workflow status (run from the workflow directory, e.g., ind_div_repair_gene/ or testis_verify/)
gwf status

# Run all pending jobs
gwf run

# Check logs for a specific target
gwf logs <target_name>
```

## Architecture

### Workflow Pattern (shared by both pipelines)

Each pipeline follows a consistent structure:

```
workflow.py                  # Entry point: creates GWF Workflow, reads config, wires steps
configurations/*.yaml        # YAML config: paths, species, parameters, conda envs
workflow_sources/
  workflow_sources.py        # Orchestration: reads config, discovers data, creates GWF targets
  workflow_templates.py      # Job templates: defines AnonymousTarget with inputs/outputs/resources/spec
```

GWF tracks job completion via `.DONE` marker files in `logs/`. Each template returns an `AnonymousTarget` with declared `inputs` and `outputs` that GWF uses for dependency resolution.

### ind_div_repair_gene Pipeline

DAG structure (per species, per gene, per family):

```
parse_gff3 (per gene) -> extract_genotypes (per individual)
                             |
                             +-> calculate_heterozygosity (per individual)
                             +-> calculate_dxy (per F-M pair)
                                      |
                                      v
                                  generate_summary (per species)
```

- Config: `configurations/repair_gene_analysis.config.yaml`
- Data layout: `data/{SPECIES}/ind/{SPECIES}_family{N}_{F|M}_{female|male}/` contains per-individual gVCF files
- Gene lists: `data/{SPECIES}/{SPECIES}_repair_gene.tsv` (gene_name, gene_id columns)
- Quality thresholds: `min_depth: 10`, `min_genotype_quality: 20`
- Output: `steps/{SPECIES}/` with per-gene/family heterozygosity, Dxy tables, and species-level summaries

### testis_verify Pipeline

Linear pipeline with per-sample parallelization at alignment stage:

```
extract_human_testis_proteins -> extract_mim_cds -> make_blast_db -> tblastn -> filter_blast_hits
                                                                                       |
star_genome_generate -> [per sample: star_align -> samtools_filter -> featurecounts] ---+
                                                                                       |
                                                                        calculate_tpm -> compare_expression
```

- Config is embedded in `workflow.py` (also mirrored in `configurations/testis_verify.config.yaml`)
- Output directories numbered sequentially: `steps/01_human_proteins/` through `steps/06_analysis/`
- Samples: MIM_ovary_5 (control), MIM_testis_candidate_red, MIM_testis_candidate_white

### Conda Environments

| Environment    | Used for                          |
|---------------|-----------------------------------|
| `gwf_new`     | Running GWF workflow manager      |
| `python_phylo`| Python analysis scripts           |
| `agat`        | AGAT toolkit (CDS extraction)     |
| `ncbi_blast`  | BLAST sequence searches           |
| `STAR`        | RNA-seq alignment                 |
| `samtools117` | BAM processing                    |
| `subread`     | featureCounts                     |

## Species Codes

| Code | Species                |
|------|------------------------|
| AFR  | Solitary spider species  |
| BIC  | Solitary spider species|
| DUM  | Social spider species|
| LIN  | Solitary spider species|
| MIM  | Social spider species  |
| SAR  | Social spider species|
| TEN  | Solitary spider species  |
