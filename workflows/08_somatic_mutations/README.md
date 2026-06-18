# 08 · Somatic mutations

- `somatic_pipeline/` — somatic DNM calling pipeline (per-species gwf configs
  `BIC/`, `LIN/`, `TEN/`; shared code in `workflow_sources/`).
- `somatic_resub/` — resubmission somatic analysis (gwf workflow; see `README.md`).
- `somatic_tables/` — 96-SBS somatic spectrum tables and plotting
  (`somatic_spectrum_table.py`, `plot_somatic_spectrum.py`).
- `testis_verify/` — testis-based validation of somatic/germline calls (gwf
  workflow).
