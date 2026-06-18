#!/bin/env python3
import sys, os, glob

sys.path.insert(0, os.path.realpath(
    os.path.join(os.path.dirname(__file__), 'workflow_sources')
))

from gwf import Workflow
from workflow_sources import coverage_workflow
from workflow_templates import combined_plot_template

gwf = Workflow()

configs = glob.glob(os.path.join(os.path.dirname(__file__), 'configurations', '*.config.yaml'))
metas = []
for config in sorted(configs):
    gwf, meta = coverage_workflow(config_file=config, gwf=gwf)
    metas.append(meta)

# One top-level cross-species barplot of per-individual autosomal DP.
if metas:
    combined_dir = os.path.join(metas[0]["steps_root"], "combined")
    gwf.target_from_template(
        name="coverage_combined_dp_plot",
        template=combined_plot_template(
            species_tsvs=[(m["species"], m["autosome_dp_tsv"]) for m in metas],
            species_dones=[m["autosome_dp_done"] for m in metas],
            combined_dir=combined_dir,
            log_path=metas[0]["log_path"],
            scripts=metas[0]["scripts"],
            conda_env=metas[0]["conda_python"],
            account=metas[0]["account"],
        ),
    )
