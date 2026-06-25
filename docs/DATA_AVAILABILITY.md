# Data availability

This repository contains **code and small metadata/result tables only**. The data
needed to run the pipeline are large and are deposited in public archives. Fill in
the accession numbers below before/at publication.

## Raw sequencing reads (parent–offspring trios)

Whole-genome short reads (DNBSEQ-G400 / Illumina) for the seven *Stegodyphus*
species and their families.

- Repository: NCBI SRA
- BioProject: PRJNA994315

## Reference genomes & annotations

Chromosome-level assemblies (HiFi + Hi-C) used as mapping references for AFR, BIC, DUM, LIN, MIM, SAR, TEN.

- Repository: NCBI Assembly
- Accessions: 
  - DUM: `GCA_044657685.1`
  - TEN: `GCA_044658465.1`
  - SAR: `GCA_044658445.1`
  - BIC: `GCA_044657725.1`
  - MIM: `GCA_044660205.1`
  - AFR: ``
  - LIN: `GCA_044657705.1`


## RNA-seq (DNA-repair gene expression)

- Repository: NCBI SRA
- BioProject: PRJNA994315

## Processed / intermediate data and results

Alignments (BAM), variant calls (VCF), callable-site BEDs, DNM call sets, mutation
spectra and other large intermediates are too large for git.

## Mapping of file-name codes

Species codes used throughout the code: `AFR, BIC, DUM, LIN, MIM, SAR, TEN`
(see the species table in the top-level [`README`](../README.md)). Note some genome
files use slightly different prefixes (`BI`→BIC, `SARA`→SAR, `TENT`→TEN); the
mapping is recorded in
`workflows/02_genome_preparation/genome_composition/genome_meta.tsv`.

## Reproducing paths

All workflow configs use absolute `/faststorage/...` paths from the original cluster.
After downloading the data above, update the `configurations/*.yaml`,
`.gwfconf.json`, and `workflow_dicts.py` entries in each `workflows/NN_*` stage to
point at your local copies.
