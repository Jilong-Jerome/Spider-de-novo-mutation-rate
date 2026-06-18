# Data availability

This repository contains **code and small metadata/result tables only**. The data
needed to run the pipeline are large and are deposited in public archives. Fill in
the accession numbers below before/at publication.

## Raw sequencing reads (parent–offspring trios)

Whole-genome short reads (DNBSEQ-G400 / Illumina) for the seven *Stegodyphus*
species and their families.

- Repository: NCBI SRA / EBI ENA
- BioProject: `TODO-PRJNAxxxxxx`
- Per-sample accessions: see `TODO` (supplementary table)

## Reference genomes & annotations

Chromosome-level assemblies (HiFi + Hi-C) and gene annotations used as mapping
references (`*_ncbi_chromosome.fa`, `*.gff3`) for AFR, BIC, DUM, LIN, MIM, SAR, TEN.

- Repository: NCBI Assembly / GenBank
- Accessions: `TODO-GCA_xxxxxxxxx` (one per species)

## RNA-seq (DNA-repair gene expression, testis validation)

- Repository: NCBI SRA / EBI ENA
- BioProject: `TODO-PRJNAxxxxxx`

## Processed / intermediate data and results

Alignments (BAM), variant calls (VCF), callable-site BEDs, DNM call sets, mutation
spectra and other large intermediates are too large for git.

- Repository: Zenodo / Figshare / Dryad
- DOI: `TODO-10.5281/zenodo.xxxxxxx`

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
