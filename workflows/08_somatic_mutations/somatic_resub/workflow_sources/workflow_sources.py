"""Workflow source definitions for somatic mutational spectrum pipeline."""

CONDA_ENVS = {
    'python': 'python_phylo',
}


def get_conda_activate(env_name: str) -> str:
    """Get conda activation command for an environment."""
    return f"CONDA_BASE=$(conda info --base) && source $CONDA_BASE/etc/profile.d/conda.sh && conda activate {env_name}"


RESOURCES = {
    'mutational_spectrum': {
        'cores': 1,
        'memory': '4g',
        'walltime': '01:00:00',
    },
}
