# 04 · GATK joint genotyping

gwf workflow (entry point `workflow.py` with `workflow_templates.py`,
`workflow_targets.py`, `workflow_dicts.py`, shared code in `workflow_sources/`).

GATK joint genotyping of trios, with per-individual/per-chromosome depth tables
(`*_per_ind_per_chrom_DP.tsv`) and offspring sex metadata
(`offspring_sex_expanded.tsv`). Variant utilities include `igv_batch_create.py`
for manual inspection.

Variants:
- `workflow_minDP/` — minimum-depth genotyping configuration.
- `workflow_X/` — X-chromosome genotyping.
