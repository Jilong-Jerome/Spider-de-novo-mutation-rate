# 01 · Quality control

- `fastqc/` — FastQC quality control of raw reads per individual (gwf workflow;
  entry point `workflow.py`, per-species `*_ind.txt` sample lists).
- `sequencing_batch/` — checks that sequencing flowcell/batch is not confounded
  with the social/subsocial contrast (gwf workflow + summary tables and
  `batch_summary_table.md`).
