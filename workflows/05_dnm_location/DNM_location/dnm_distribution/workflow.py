#!/usr/bin/env python3
import sys
import os
import glob
import yaml

sys.path.insert(0, os.path.realpath(
    '/faststorage/project/spider2/social_spiders_2020/people/jilong/chapter3_dnm/science_resub/DNM_location/dnm_distribution/workflow_sources/'
))

from gwf import Workflow
from workflow_sources import dnm_distribution_workflow
from workflow_templates import plot_nearest_neighbor_all_species_template

gwf = Workflow()

configs = sorted(glob.glob('./configurations/*.config.yaml'))
for config in configs:
    gwf = dnm_distribution_workflow(config_file=config, gwf=gwf)

species_results = []
common_output_path = None
common_log_path = None
common_scripts = None

for config in configs:
    with open(config) as fh:
        cfg = yaml.safe_load(fh)

    SP = cfg['project_id']
    sp = cfg['species_prefix']
    output_path = cfg['output_directory_path']
    common_output_path = common_output_path or output_path
    common_log_path = common_log_path or cfg['log_directory_path']
    common_scripts = common_scripts or cfg['scripts_path']

    species_results.append({
        'SP': SP,
        'simulations_tsv': (
            f"{output_path}/nearest_neighbor_test/{SP}/"
            f"{sp}_nearest_neighbor_simulations.tsv"
        ),
        'summary': (
            f"{output_path}/nearest_neighbor_test/{SP}/"
            f"{sp}_nearest_neighbor_summary.txt"
        ),
    })

gwf.target_from_template(
    name="all_species_nearest_neighbor_plot",
    template=plot_nearest_neighbor_all_species_template(
        species_results=species_results,
        output_path=common_output_path,
        log_path=common_log_path,
        scripts=common_scripts,
    ),
)
