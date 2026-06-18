"""Workflow sources and templates for testis verification pipeline."""

from .workflow_sources import CONDA_ENVS, RESOURCES, get_conda_activate, get_resources
from .workflow_templates import (
    extract_human_testis_proteins,
    extract_mim_cds,
    make_blast_db,
    run_tblastn,
    filter_blast_hits,
    star_genome_generate,
    star_align,
    samtools_filter,
    feature_counts,
    calculate_tpm,
    compare_expression,
)

__all__ = [
    'CONDA_ENVS',
    'RESOURCES',
    'get_conda_activate',
    'get_resources',
    'extract_human_testis_proteins',
    'extract_mim_cds',
    'make_blast_db',
    'run_tblastn',
    'filter_blast_hits',
    'star_genome_generate',
    'star_align',
    'samtools_filter',
    'feature_counts',
    'calculate_tpm',
    'compare_expression',
]
