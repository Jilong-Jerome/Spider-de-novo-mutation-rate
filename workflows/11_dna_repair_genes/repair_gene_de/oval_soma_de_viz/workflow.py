#!/usr/bin/env python3
import sys, os, glob

sys.path.insert(0, os.path.realpath(
    os.path.join(os.path.dirname(__file__), 'workflow_sources')
))

from gwf import Workflow
from workflow_sources import (
    pathway_bias_bar_workflow,
    combined_pathway_bias_bar_workflow,
)

gwf = Workflow()

configs = glob.glob(os.path.join(os.path.dirname(__file__), 'configurations', '*.config.yaml'))
for config in configs:
    gwf = pathway_bias_bar_workflow(config_file=config, gwf=gwf)

gwf = combined_pathway_bias_bar_workflow(config_files=configs, gwf=gwf)
