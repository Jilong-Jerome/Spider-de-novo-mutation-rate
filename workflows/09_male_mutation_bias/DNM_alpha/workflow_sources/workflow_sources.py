#!/usr/bin/env python3
"""
workflow_sources.py
Build GWF targets for the DNM alpha mutation-sharing workflow.
"""
import yaml
from workflow_templates import (
    mutation_sharing_summary_template,
    mutation_phasing_summary_template,
    parental_rate_summary_template,
    mutation_sharing_plot_template,
)


def dnm_alpha_workflow(config_file, gwf):
    with open(config_file) as f:
        config = yaml.safe_load(f)

    account = config['account']
    output_path = config['output_directory_path']
    log_path = config['log_directory_path']
    scripts = config['scripts_path']
    python_env = config['python_conda_env']
    dnm_tsv = config['dnm_tsv']
    callable_dir = config['callable_directory_path']

    summary_tsv = f'{output_path}/mutation_sharing/species_mutation_sharing_counts.tsv'
    phasing_tsv = f'{output_path}/mutation_sharing/species_mutation_phasing_counts.tsv'
    parental_rates_tsv = f'{output_path}/mutation_sharing/subsocial_parental_mutation_rates.tsv'

    gwf.target_from_template(
        name='mutation_sharing_summary',
        template=mutation_sharing_summary_template(
            config_file=config_file,
            dnm_tsv=dnm_tsv,
            output_path=output_path,
            log_path=log_path,
            scripts=scripts,
            python_conda_env=python_env,
            account=account,
        ),
    )

    gwf.target_from_template(
        name='mutation_phasing_summary',
        template=mutation_phasing_summary_template(
            config_file=config_file,
            dnm_tsv=dnm_tsv,
            output_path=output_path,
            log_path=log_path,
            scripts=scripts,
            python_conda_env=python_env,
            account=account,
        ),
    )

    gwf.target_from_template(
        name='parental_rate_summary',
        template=parental_rate_summary_template(
            config_file=config_file,
            dnm_tsv=dnm_tsv,
            callable_dir=callable_dir,
            output_path=output_path,
            log_path=log_path,
            scripts=scripts,
            python_conda_env=python_env,
            account=account,
        ),
    )

    gwf.target_from_template(
        name='mutation_sharing_plot',
        template=mutation_sharing_plot_template(
            config_file=config_file,
            summary_tsv=summary_tsv,
            phasing_tsv=phasing_tsv,
            parental_rates_tsv=parental_rates_tsv,
            output_path=output_path,
            log_path=log_path,
            scripts=scripts,
            python_conda_env=python_env,
            account=account,
        ),
    )

    return gwf
