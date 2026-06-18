"""Workflow sources and templates for somatic mutational spectrum pipeline."""

from .workflow_sources import CONDA_ENVS, RESOURCES, get_conda_activate
from .workflow_templates import mutational_spectrum_template

__all__ = [
    'CONDA_ENVS',
    'RESOURCES',
    'get_conda_activate',
    'mutational_spectrum_template',
]
