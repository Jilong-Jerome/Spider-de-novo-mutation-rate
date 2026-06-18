"""Command templates for somatic mutational spectrum workflow.

Each function returns an AnonymousTarget with inputs, outputs, options, and spec.
"""

from pathlib import Path
from gwf import AnonymousTarget
from .workflow_sources import get_conda_activate, CONDA_ENVS, RESOURCES


def mutational_spectrum_template(species: str, input_tsv: Path,
                                  output_dir: Path, log_dir: Path):
    """Count mutational spectrum for a single species."""
    conda_cmd = get_conda_activate(CONDA_ENVS['python'])

    inputs = {
        'tsv': str(input_tsv),
    }
    outputs = {
        'spectrum': f"{output_dir}/{species}_mutational_spectrum.tsv",
        'done': f"{log_dir}/mutational_spectrum_{species}.DONE",
    }
    options = {**RESOURCES['mutational_spectrum'], 'account': 'spider2'}

    spec = f"""
{conda_cmd}

mkdir -p {output_dir} {log_dir}

python scripts/count_mutational_spectrum.py \\
    --input {input_tsv} \\
    --species {species} \\
    --output {outputs['spectrum']}

echo done > {outputs['done']}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
