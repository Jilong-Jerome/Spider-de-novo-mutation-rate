#!/usr/bin/env python3
"""
workflow.py - Entry point for repair gene genotype analysis GWF workflow

This workflow analyzes repair gene genotypes from gVCF files, calculating:
- Heterozygosity per individual (whole gene, exon, intron)
- List of heterozygous sites per individual
- Dxy (absolute divergence) between female-male pairs

Usage:
    conda activate gwf_new
    gwf status
    gwf run
"""

import gwf
import sys
import os
import yaml
import glob

# Add workflow_sources to path
sys.path.insert(0, os.path.realpath('./workflow_sources/'))

from workflow_sources import repair_gene_analysis_workflow

# Find all configuration files
configs = glob.glob('./configurations/*config.y*ml')

# Create workflow
gwf = gwf.Workflow()

# Process each configuration
for config in configs:
    gwf = repair_gene_analysis_workflow(config_file=config, gwf=gwf)
