#!/usr/bin/env python3
"""
workflow.py - Entry point for dS branch trio analysis GWF workflow

Maps SAR and PAC reads to BIC reference, builds consensus genomes,
and extracts per-gene CDS alignments for dN/dS analysis.

Usage:
    conda activate gwf_new
    gwf status
    gwf run
"""

import gwf
import sys
import os
import glob

# Add workflow_sources to path
sys.path.insert(0, os.path.realpath('./workflow_sources/'))

from workflow_sources import (
    trio_dnds_workflow,
    combine_trio_visualization_workflow,
)

# Find all configuration files
configs = glob.glob('./configurations/*config.y*ml')

# Create workflow
gwf = gwf.Workflow()

# Process each configuration
for config in configs:
    gwf = trio_dnds_workflow(config_file=config, gwf=gwf)

# Cross-trio combined visualisation (depends on every per-trio Step 09)
gwf = combine_trio_visualization_workflow(config_files=configs, gwf=gwf)
