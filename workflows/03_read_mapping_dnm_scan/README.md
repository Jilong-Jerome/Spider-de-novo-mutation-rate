# 03 · Read mapping & germline DNM scanning

gwf workflow (entry point `workflow.py`, with `workflow_templates.py`,
`workflow_targets.py`, `workflow_dicts.py`).

Steps:
- `merge_reads/` — merge per-lane FASTQs into per-sample reads.
- Mapping reads to the species reference (bwa-mem2) and BAM processing.
- Germline *de novo* mutation scanning across parent–offspring trios
  (`germline_scan.py`, `germline_scan_AB.py`) and candidate filtering
  (`callable_GQ_DP_filter*.py`, `distance_filter.py`,
  `bcftools_distance_filter.py`, `vcf_site_count.py`).

Supporting inputs: `intervals/` (per-species scan intervals), `lists/` (read
lists), `sample_maps/` (sample-to-file maps).
