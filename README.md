# Spider de novo mutation rate

Analysis code accompanying the manuscript:

> **Convergent Evolution of Sociality Causes Reduction of Mutation Rates in Spiders**
> Preprint: https://www.biorxiv.org/content/10.1101/2025.07.09.663850v1
> (doi:10.1101/2025.07.09.663850)

This repository deposits the complete analysis pipeline used to estimate germline
and somatic *de novo* mutation (DNM) rates from parent–offspring trios across seven
*Stegodyphus* spider species with independent origins of sociality, and the
downstream analyses of mutation spectrum, callable genome, parental bias, dN/dS,
DNA-repair genes and kinship.

> **Scope.** Only source code and small metadata/result tables are tracked here.
> Raw reads, reference genomes, alignments, variant calls and large intermediate
> files are **not** in git — see [`docs/DATA_AVAILABILITY.md`](docs/DATA_AVAILABILITY.md).

## Study system

Seven *Stegodyphus* species, three social and four subsocial (sociality evolved
independently multiple times):

| Code | Species              | Sociality |
|------|----------------------|-----------|
| DUM  | *S. dumicola*        | Social    |
| MIM  | *S. mimosarum*       | Social    |
| SAR  | *S. sarasinorum*     | Social    |
| LIN  | *S. lineatus*        | Subsocial |
| TEN  | *S. tentoriicola*    | Subsocial |
| AFR  | *S. africanus*       | Subsocial |
| BIC  | *S. bicolor*         | Subsocial |

## Repository layout

```
workflows/                     analysis pipeline, in reproduction order
  01_quality_control/          FastQC; sequencing-batch confound checks
  02_genome_preparation/       NCBI genome cleanup, indexing, genome composition / stats
  03_read_mapping_dnm_scan/    read merging, mapping, germline DNM scanning & filtering
  04_genotyping_gatk/          GATK joint genotyping (autosomes, minDP, X chromosome)
  05_dnm_location/             DNM genomic location & inter-individual variance
  06_mutation_rate/            callable-genome estimation, rate tables, supplement tables
  07_mutation_spectrum/        germline mutation spectrum; germline-vs-somatic comparison
  08_somatic_mutations/        somatic DNM calling, spectrum tables, testis validation
  09_male_mutation_bias/       male mutation bias (alpha); parental-origin / shared DNMs
  10_dnds_branch/              branch dN/dS (PAML) and recombination-associated analyses
  11_dna_repair_genes/         DNA-repair gene divergence & differential expression
  12_kinship/                  trio/kinship verification (PLINK, relatedness)
  13_chromosome_X/             X-chromosome depth & per-chromosome statistics
  14_mating_simulation/        mating / generation-overlap simulation
environment/                   conda environment documentation + export helper
docs/                          collected supplemental methods + data availability
```

Each `workflows/NN_*/` (or its named sub-analyses) is a self-contained
[gwf](https://gwf.app) project. Many sub-directories carry their own `README.md`
and `SUPPLEMENTAL_METHODS*.md` describing that analysis in detail.

## Software environments

Analyses were run on a SLURM cluster and orchestrated with **gwf**. Each step
activates a named conda environment (e.g. `conda activate bwa2`). The full list of
environments and their purpose, plus a helper to export pinned `*.yml` files, is in
[`environment/`](environment/README.md).

## Running a workflow

```bash
conda activate gwf_new          # gwf workflow manager
cd workflows/03_read_mapping_dnm_scan
gwf status                      # inspect the DAG / job states
gwf run                         # submit pending jobs to SLURM
```

Before running, edit the paths and configuration files (`configurations/*.yaml`,
`.gwfconf.json`, and the per-workflow `workflow_dicts.py` / config blocks) to point
at your local copies of the data described in `docs/DATA_AVAILABILITY.md`. The
hard-coded `/faststorage/...` paths reflect the original compute environment and
must be adapted.

## Citation

If you use this code, please cite the manuscript above. Repository citation
metadata is in [`CITATION.cff`](CITATION.cff).

## License

Released under the MIT License — see [`LICENSE`](LICENSE).
