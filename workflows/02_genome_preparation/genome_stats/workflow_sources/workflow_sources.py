from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def run_genome_stats(gwf, work_path, data_path, log_path, species):
    gwf.target_from_template(
        name = f"{species}_genome_stats",
        template = genome_stats_template(work_path, data_path, log_path, species)
    )
    return gwf

def run_genome_comparison(gwf, work_path, log_path, species_list):
    gwf.target_from_template(
        name = "genome_comparison",
        template = genome_comparison_template(work_path, log_path, species_list)
    )
    return gwf
