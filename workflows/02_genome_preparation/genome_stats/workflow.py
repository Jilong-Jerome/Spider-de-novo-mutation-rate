#!/bin/env python3
import sys, os, yaml, glob
from gwf import *
sys.path.insert(0, os.path.realpath('/faststorage/project/spider2/social_spiders_2020/people/jilong/chapter3_dnm/claude/genome_stat_scripts_claude/workflow_sources/')) 
from workflow_sources import *

gwf = Workflow()

# Define species list - modify this based on your available genomes
species_list = ["BIC", "SAR", "TEN", "DUM", "AFR", "MIM", "LIN"]

# Define paths
base_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm"
work_path = f"{base_path}/steps/genome_stats_claude"
data_path = f"{base_path}/data/genome_claude"  # Assuming genomes are stored here
log_path = f"{base_path}/steps/genome_stats_claude/logs"

# Run genome statistics for each species
for species in species_list:
    run_genome_stats(gwf, work_path, data_path, log_path, species)

# Run comparative analysis after all individual analyses
run_genome_comparison(gwf, work_path, log_path, species_list)
