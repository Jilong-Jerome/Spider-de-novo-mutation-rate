"""
Somatic Mutational Spectrum GWF Workflow

Summarize somatic mutations per species by counting mutation types (C>A, C>G,
C>T, T>A, T>C, T>G) stratified by trinucleotide context.

Usage:
    gwf status
    gwf run
"""

from pathlib import Path

from gwf import Workflow

from workflow_sources import mutational_spectrum_template

# Initialize workflow
gwf = Workflow(defaults={'account': 'spider2'})

# ============================================================================
# Configuration
# ============================================================================

SPECIES_LIST = ['AFR', 'BIC', 'DUM', 'LIN', 'MIM', 'SAR', 'TEN']
DATA_DIR = Path('data')
STEPS_DIR = Path('steps')
LOGS_DIR = Path('logs')

# ============================================================================
# Create one target per species
# ============================================================================

for species in SPECIES_LIST:
    input_tsv = DATA_DIR / f'{species}_after_first_filter.tsv'
    output_dir = STEPS_DIR / species

    gwf.target_from_template(
        name=f'mutational_spectrum_{species}',
        template=mutational_spectrum_template(
            species=species,
            input_tsv=input_tsv,
            output_dir=output_dir,
            log_dir=LOGS_DIR,
        ),
    )
