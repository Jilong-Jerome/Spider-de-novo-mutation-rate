# 03 · Read mapping & GATK genotyping

gwf workflow (entry point `workflow.py`, with `workflow_templates.py`,
`workflow_targets.py`, `workflow_dicts.py`).

Steps:
- Genome indexing of each species reference.
- `merge_reads/` — merge per-lane FASTQs into per-sample reads.
- Mapping reads to the species reference (bwa-mem2) and BAM processing.
- GATK joint genotyping: per-individual gVCFs → per-chromosome GenomicsDB →
  joint genotyping → combine/finalise, with per-individual/per-chromosome mean
  depth and per-chromosome VCF splitting. ['A standalone portable pipeline for the GATK step'](https://github.com/Jilong-Jerome/GATK_genotype_calling/tree/main)

Supporting inputs: `intervals/` (per-species intervals), `lists/` (read lists),
`sample_maps/` (sample-to-file maps).

> Germline *de novo* mutation calling from the genotyped VCFs is the next stage,
> [`04_germline_dnm_scan/`](../04_germline_dnm_scan/).
