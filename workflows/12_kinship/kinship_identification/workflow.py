#!/usr/bin/env python3
"""
workflow.py
GWF entry point for the pairwise kinship identification workflow.

Discovers every configurations/*.config.yaml file and builds the per-species
kinship-inference targets for each one.
"""
import glob
import os
import sys

sys.path.insert(0, os.path.realpath(
    os.path.join(os.path.dirname(__file__), 'workflow_sources')
))

from gwf import Workflow
from workflow_sources import kinship_workflow, combined_trio_plot

gwf = Workflow()

configs = glob.glob(os.path.join(os.path.dirname(__file__), 'configurations', '*.config.yaml'))
for config in configs:
    gwf = kinship_workflow(config_file=config, gwf=gwf)

# One cross-species combined supplemental figure (downstream of every trio_check).
gwf = combined_trio_plot(config_files=configs, gwf=gwf)
