#!/usr/bin/env python3
"""
workflow_sources.py
Builds gwf targets for the DNM location distribution analysis.
"""
import os
import glob
import yaml

from workflow_templates import (
    compute_window_callable_template,
    merge_window_callable_template,
    dnm_distribution_test_template,
    plot_dnm_genome_template,
    compute_callable_bed_template,
    plot_callable_segments_template,
    nearest_neighbor_test_template,
)


def dnm_distribution_workflow(config_file, gwf):
    """
    Register all gwf targets for one species given a config YAML.
    Returns the gwf Workflow object (modified in place, but returned
    for chaining convenience).
    """
    with open(config_file) as fh:
        cfg = yaml.safe_load(fh)

    SP  = cfg['project_id']          # e.g. AFR
    sp  = cfg['species_prefix']      # e.g. afr
    CALLABLE_GENOME_DIR = cfg['callable_genome_dir']
    CALLABLE_SUMMARY    = cfg['callable_summary_file']
    DNM_FILE            = cfg['dnm_file']
    GENOME_FAI          = cfg['genome_fai']
    OUTPUT_PATH         = cfg['output_directory_path']
    LOG_PATH            = cfg['log_directory_path']
    SCRIPTS             = cfg['scripts_path']
    WINDOW_SIZE         = cfg['window_size']
    X_CHROMS            = cfg['x_chromosomes']

    # ------------------------------------------------------------------ #
    # 1. Per-offspring: compute window callable sites
    # ------------------------------------------------------------------ #
    offspring_dirs = sorted(glob.glob(os.path.join(CALLABLE_GENOME_DIR, '*_offspring')))
    done_files = []

    for offspring_dir in offspring_dirs:
        offspring = os.path.basename(offspring_dir)   # e.g. AFR_family1_S1_offspring
        target_name = f"{SP}_{offspring}_window_callable"

        template = compute_window_callable_template(
            offspring_dir=offspring_dir,
            offspring=offspring,
            sp=sp,
            SP=SP,
            output_path=OUTPUT_PATH,
            log_path=LOG_PATH,
            scripts=SCRIPTS,
            genome_fai=GENOME_FAI,
            window_size=WINDOW_SIZE,
            x_chroms=X_CHROMS,
        )
        gwf.target_from_template(name=target_name, template=template)
        done_files.append(template.outputs['log'])

    # ------------------------------------------------------------------ #
    # 2. Merge per-offspring counts into species-level table
    # ------------------------------------------------------------------ #
    merge_template = merge_window_callable_template(
        offspring_done_files=done_files,
        sp=sp,
        SP=SP,
        output_path=OUTPUT_PATH,
        log_path=LOG_PATH,
        scripts=SCRIPTS,
    )
    gwf.target_from_template(
        name=f"{SP}_merge_window_callable",
        template=merge_template,
    )

    merged_tsv  = merge_template.outputs['merged']
    merge_done  = merge_template.outputs['log']

    # ------------------------------------------------------------------ #
    # 3. Chi-squared distribution test
    # ------------------------------------------------------------------ #
    gwf.target_from_template(
        name=f"{SP}_dnm_distribution_test",
        template=dnm_distribution_test_template(
            merged_tsv=merged_tsv,
            merge_done=merge_done,
            dnm_file=DNM_FILE,
            sp=sp,
            SP=SP,
            output_path=OUTPUT_PATH,
            log_path=LOG_PATH,
            scripts=SCRIPTS,
            window_size=WINDOW_SIZE,
            x_chroms=X_CHROMS,
        ),
    )

    # ------------------------------------------------------------------ #
    # 4. Plot: genome-wide DNM tick-mark visualisation
    # ------------------------------------------------------------------ #
    gwf.target_from_template(
        name=f"{SP}_plot_dnm_genome",
        template=plot_dnm_genome_template(
            dnm_file=DNM_FILE,
            genome_fai=GENOME_FAI,
            sp=sp,
            SP=SP,
            output_path=OUTPUT_PATH,
            log_path=LOG_PATH,
            scripts=SCRIPTS,
            x_chroms=X_CHROMS,
        ),
    )

    # ------------------------------------------------------------------ #
    # 5. Per-offspring: compute callable BED segments from VCFs
    # ------------------------------------------------------------------ #
    bed_done_files = []

    for offspring_dir in offspring_dirs:
        offspring   = os.path.basename(offspring_dir)
        target_name = f"{SP}_{offspring}_callable_bed"

        bed_template = compute_callable_bed_template(
            offspring_dir=offspring_dir,
            offspring=offspring,
            sp=sp,
            SP=SP,
            output_path=OUTPUT_PATH,
            log_path=LOG_PATH,
            scripts=SCRIPTS,
            genome_fai=GENOME_FAI,
            x_chroms=X_CHROMS,
        )
        gwf.target_from_template(name=target_name, template=bed_template)
        bed_done_files.append(bed_template.outputs['log'])

    # ------------------------------------------------------------------ #
    # 6. Plot: callable genome segments per individual per chromosome
    # ------------------------------------------------------------------ #
    gwf.target_from_template(
        name=f"{SP}_plot_callable_segments",
        template=plot_callable_segments_template(
            bed_done_files=bed_done_files,
            sp=sp,
            SP=SP,
            output_path=OUTPUT_PATH,
            log_path=LOG_PATH,
            scripts=SCRIPTS,
            genome_fai=GENOME_FAI,
            dnm_file=DNM_FILE,
            x_chroms=X_CHROMS,
        ),
    )

    # ------------------------------------------------------------------ #
    # 7. Nearest-neighbour permutation test
    # ------------------------------------------------------------------ #
    gwf.target_from_template(
        name=f"{SP}_nearest_neighbor_test",
        template=nearest_neighbor_test_template(
            bed_done_files=bed_done_files,
            sp=sp,
            SP=SP,
            output_path=OUTPUT_PATH,
            log_path=LOG_PATH,
            scripts=SCRIPTS,
            dnm_file=DNM_FILE,
            x_chroms=X_CHROMS,
        ),
    )

    return gwf
