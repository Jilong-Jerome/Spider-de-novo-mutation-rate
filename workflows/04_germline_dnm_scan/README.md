# 04 · Germline DNM scanning

gwf workflow (entry point `workflow.py` with `workflow_templates.py`,
`workflow_targets.py`, `workflow_dicts.py`, shared code in `workflow_sources/`).

Calls germline *de novo* mutations from the genotyped trio VCFs produced in
[`03_read_mapping_genotyping/`](../03_read_mapping_genotyping/):

- callable-site definition per trio (`run_callable_sites_per_trio`);
- *de novo* mutation scan across parent–offspring trios (`run_dnm_scan`,
  `germline_scan.py`, `germline_scan_AB.py`);
- candidate filtering by genotype quality / depth and clustering distance
  (`callable_GQ_DP_filter*.py`, `distance_filter.py`,
  `bcftools_distance_filter.py`, `vcf_site_count.py`);
- BAM checks and IGV verification of candidate sites (`run_bam_check`,
  `igv_batch_create.py`).

Inputs/metadata: per-individual/per-chromosome depth tables
(`*_per_ind_per_chrom_DP.tsv`) and offspring sex (`offspring_sex_expanded.tsv`).

Variants:
- `workflow_minDP/` — minimum-depth scanning configuration.
- `workflow_X/` — X-chromosome scanning.
