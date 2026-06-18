"""
Config reader and job orchestrator for the sequencing batch summary workflow.
"""

import os
import re

import yaml

from workflow_templates import (
    add_coverage_template,
    concatenate_results_template,
    plot_coverage_scatter_template,
    plot_heatmap_template,
    summarize_individual_template,
)


def load_config(config_file):
    with open(config_file) as fh:
        return yaml.safe_load(fh)


def sequencing_batch_workflow(config_file, gwf):
    config = load_config(config_file)

    data_dir = config['data_directory']
    output_path = config['output_directory_path']
    log_path = config['log_directory_path']
    scripts_path = config['scripts_path']
    final_output = config['final_output']
    fai_directory = config['fai_directory']
    coverage_output = config['coverage_output']
    heatmap_output = config['heatmap_output']
    scatter_output = config['scatter_output']
    account = config['account']

    partial_tsvs = []
    done_logs = []

    for species in sorted(os.listdir(data_dir)):
        species_dir = os.path.join(data_dir, species)
        if not os.path.isdir(species_dir):
            continue

        for family in sorted(os.listdir(species_dir)):
            if family == 'pseudodiploid':
                continue
            family_dir = os.path.join(species_dir, family)
            if not os.path.isdir(family_dir):
                continue

            for individual in sorted(os.listdir(family_dir)):
                individual_dir = os.path.join(family_dir, individual)
                if not os.path.isdir(individual_dir):
                    continue

                # Skip individuals with no .fq.gz files
                fq_files = [
                    f for f in os.listdir(individual_dir)
                    if f.endswith('.fq.gz')
                ]
                if not fq_files:
                    continue

                safe_name = re.sub(r'[^a-zA-Z0-9]', '_',
                                   f'{species}_{family}_{individual}')

                template = summarize_individual_template(
                    individual_dir=individual_dir,
                    species=species,
                    family=family,
                    individual=individual,
                    safe_name=safe_name,
                    output_path=output_path,
                    log_path=log_path,
                    scripts_path=scripts_path,
                    account=account,
                )
                gwf.target_from_template(name=safe_name, template=template)

                partial_tsvs.append(f'{output_path}/{safe_name}_summary.tsv')
                done_logs.append(f'{log_path}/{safe_name}_summarize.DONE')

    # Final concatenation target
    concat_template = concatenate_results_template(
        partial_tsvs=partial_tsvs,
        done_logs=done_logs,
        final_output=final_output,
        log_path=log_path,
        account=account,
    )
    gwf.target_from_template(name='concatenate_results', template=concat_template)

    # Add coverage column (depends on concatenate_results via final_output input)
    cov_template = add_coverage_template(
        summary_tsv=final_output,
        fai_directory=fai_directory,
        coverage_output=coverage_output,
        log_path=log_path,
        scripts_path=scripts_path,
        account=account,
    )
    gwf.target_from_template(name='add_coverage', template=cov_template)

    # Plot sequencing-scheme heatmap (depends on add_coverage via coverage_output input)
    heatmap_template = plot_heatmap_template(
        coverage_tsv=coverage_output,
        heatmap_output=heatmap_output,
        log_path=log_path,
        scripts_path=scripts_path,
        account=account,
    )
    gwf.target_from_template(name='plot_heatmap', template=heatmap_template)

    # Scatter plot of total coverage per individual (depends on add_coverage)
    scatter_template = plot_coverage_scatter_template(
        coverage_tsv=coverage_output,
        scatter_output=scatter_output,
        log_path=log_path,
        scripts_path=scripts_path,
        account=account,
    )
    gwf.target_from_template(name='plot_coverage_scatter', template=scatter_template)
