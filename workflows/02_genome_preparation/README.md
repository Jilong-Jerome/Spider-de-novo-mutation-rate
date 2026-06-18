# 02 · Genome preparation

- `ncbi/` — prepare reference assemblies for NCBI: chromosome renaming
  (`rename_chromosome.py`), contamination trimming (`trim_contamination.py`),
  indexing and genome statistics (gwf workflow `workflow.py`).
- `genome_composition/` — callable-site context and genome composition summaries
  (`genome_context_summary.py`, `callable_site_split_count.py`); `genome_meta.tsv`
  maps species codes to genome files.
- `genome_stats/` — comparative genome statistics across species
  (`species_genomes.txt`, gwf workflow).
