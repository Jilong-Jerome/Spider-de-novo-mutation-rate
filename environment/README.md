# Software environments

All workflows were run on a SLURM cluster and orchestrated with
[gwf](https://gwf.app). Each gwf target activates a named **conda** environment
before running its commands (look for `conda activate <env>` inside the
`workflow_templates.py` / `workflow_sources/` of each stage).

The environments below were used across the pipeline. Exact package versions are
whatever was installed in the original compute environment; to capture them as
reproducible `*.yml` files, run [`export_envs.sh`](export_envs.sh) on a machine
that has these conda environments installed (see bottom of this file).

| conda env       | Purpose                                            | Key tools |
|-----------------|----------------------------------------------------|-----------|
| `gwf_new`       | Workflow management / job submission               | gwf |
| `fastqc`        | Read quality control                               | FastQC |
| `bwa2`          | Short-read mapping                                 | bwa-mem2, samtools |
| `pbmm2`         | PacBio HiFi mapping / minimap2 alignment           | pbmm2, minimap2 |
| `samtools117`   | BAM processing (samtools 1.17)                     | samtools |
| `bam_check`     | BAM validation / QC                                | samtools, picard |
| `gatk4`         | Joint genotyping / variant calling                 | GATK4 |
| `bcftools`      | VCF manipulation & filtering                       | bcftools |
| `vcftools`      | VCF statistics & filtering                         | vcftools |
| `bedtools`      | Interval / genome arithmetic                       | bedtools |
| `pysam`         | Custom BAM/VCF parsing in Python                   | pysam |
| `pyvcf`         | VCF parsing in Python                              | PyVCF |
| `biopython`    | Sequence handling / consensus building             | Biopython |
| `pandas`        | Tabular data wrangling & plotting                  | pandas, numpy, matplotlib |
| `python_phylo`  | Phylogenetic / spectrum analysis & plotting        | pandas, matplotlib, ete3/dendropy |
| `paml`          | Branch dN/dS estimation                            | PAML (codeml) |
| `plink`         | Relatedness / kinship checks                       | PLINK |
| `liftoff`       | Annotation lift-over between assemblies            | Liftoff |
| `spider`        | General project environment                        | mixed Python/bioinformatics |
| `base`          | Default conda base for simple steps                | — |

> The table is a best-effort description inferred from the workflow code. Treat the
> tool lists as indicative; the authoritative, version-pinned record is produced by
> `export_envs.sh`.

## Exporting pinned environment files

On a machine where these environments exist:

```bash
cd environment
bash export_envs.sh          # writes <env>.yml for each env that is installed
```

This produces `environment/<env>.yml` files (via `conda env export`) that can be
recreated elsewhere with:

```bash
conda env create -f environment/bwa2.yml
```
