#!/usr/bin/env python3
"""
Build GWF targets for final visual presentations and companion datasets.
"""
import yaml
from workflow_templates import combined_rate_visual_template


def dnm_final_visuals_workflow(config_file, gwf):
    with open(config_file) as f:
        config = yaml.safe_load(f)

    gwf.target_from_template(
        name='final_combined_rate_by_mutation_class',
        template=combined_rate_visual_template(
            per_trio_germline_tsv=config['per_trio_germline_tsv'],
            per_species_somatic_tsv=config['per_species_somatic_tsv'],
            spectrum_test_tsv=config['all_species_spectrum_test_tsv'],
            output_path=config['output_directory_path'],
            log_path=config['log_directory_path'],
            scripts=config['scripts_path'],
            python_conda_env=config['python_conda_env'],
            account=config['account'],
            bootstrap_replicates=int(config.get('bootstrap_replicates', 10000)),
            bootstrap_seed=int(config.get('bootstrap_seed', 20260430)),
        ),
    )

    return gwf
