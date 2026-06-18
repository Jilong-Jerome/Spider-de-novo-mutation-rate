"""
Workflow source definitions for testis verification pipeline.

Defines conda environments and resource requirements for each tool.
"""

# Conda environment definitions
CONDA_ENVS = {
    'python': 'python_phylo',
    'agat': 'agat',
    'blast': 'ncbi_blast',
    'star': 'STAR',
    'samtools': 'samtools117',
    'featurecounts': 'subread',
}


def get_conda_activate(env_name: str) -> str:
    """Get conda activation command for an environment."""
    return f"CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh && conda activate {env_name}"


# Resource specifications
RESOURCES = {
    'extract_proteins': {
        'cores': 1,
        'memory': '4g',
        'walltime': '00:30:00',
    },
    'agat': {
        'cores': 1,
        'memory': '8g',
        'walltime': '01:00:00',
    },
    'makeblastdb': {
        'cores': 1,
        'memory': '4g',
        'walltime': '00:30:00',
    },
    'tblastn': {
        'cores': 8,
        'memory': '16g',
        'walltime': '04:00:00',
    },
    'filter_blast': {
        'cores': 1,
        'memory': '4g',
        'walltime': '00:30:00',
    },
    'star_build': {
        'cores': 8,
        'memory': '64g',
        'walltime': '04:00:00',
    },
    'star_align': {
        'cores': 8,
        'memory': '64g',
        'walltime': '12:00:00',
    },
    'samtools_filter': {
        'cores': 4,
        'memory': '16g',
        'walltime': '02:00:00',
    },
    'featurecounts': {
        'cores': 4,
        'memory': '8g',
        'walltime': '01:00:00',
    },
    'calculate_tpm': {
        'cores': 1,
        'memory': '4g',
        'walltime': '00:30:00',
    },
    'compare_expression': {
        'cores': 1,
        'memory': '8g',
        'walltime': '00:30:00',
    },
}


def get_resources(task_name: str) -> dict:
    """Get resource specifications for a task."""
    return RESOURCES.get(task_name, {
        'cores': 1,
        'memory': '4g',
        'walltime': '01:00:00',
    })
