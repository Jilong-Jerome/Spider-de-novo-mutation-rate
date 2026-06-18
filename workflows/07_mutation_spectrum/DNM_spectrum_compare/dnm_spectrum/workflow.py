#!/usr/bin/env python3
import sys
import os
import glob

sys.path.insert(0, os.path.realpath(
    os.path.join(os.path.dirname(__file__), 'workflow_sources')
))

from gwf import Workflow
from workflow_sources import dnm_spectrum_workflow

gwf = Workflow()

configs = glob.glob(os.path.join(os.path.dirname(__file__), 'configurations', '*.config.yaml'))
for config in configs:
    gwf = dnm_spectrum_workflow(config_file=config, gwf=gwf)
