# 10 · Branch dN/dS & DNA-repair gene presence/absence

- `dS_branch_trio/` — branch-specific dN/dS for species trios using PAML codeml
  (gwf workflow; see `README.md`, `SUPPLEMENTAL_METHODS.md`, `DEVELOPMENT.md`,
  per-trio `configurations/*.config.yaml`).
- `hyper/` — screen for unique presence/absence (gene/exon loss) of DNA-repair
  genes (MLH1, PMS2, MSH2, MSH6, OGG1, TDG, MUTYH, MBD4, ERCC2), with a ~600-gene
  random background, across individuals. Uses liftoff-mapped gene annotations and
  per-feature sequencing coverage from the genotyping VCFs (gwf workflow;
  `workflow_sources/check_gene_features` step, `check_unique_absence.py`,
  `feature_summary.py`, `summary_random.py`).
