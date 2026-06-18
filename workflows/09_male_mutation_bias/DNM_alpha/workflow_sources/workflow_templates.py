#!/usr/bin/env python3
"""
workflow_templates.py
GWF target templates for the DNM alpha mutation-sharing workflow.
"""
from gwf import AnonymousTarget


def mutation_sharing_summary_template(
    config_file, dnm_tsv, output_path, log_path,
    scripts, python_conda_env, account,
):
    out_dir = f'{output_path}/mutation_sharing'
    out_tsv = f'{out_dir}/species_mutation_sharing_counts.tsv'
    done_file = f'{log_path}/mutation_sharing_summary.DONE'

    inputs = {
        'config': config_file,
        'dnm': dnm_tsv,
    }
    outputs = {
        'summary': out_tsv,
        'log': done_file,
    }
    options = {
        'cores': 1,
        'memory': '2g',
        'walltime': '00:15:00',
        'account': account,
    }

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {python_conda_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/summarize_shared_mutations.py \\
    --config {config_file} \\
    --dnm_tsv {dnm_tsv} \\
    --output_tsv {out_tsv}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def mutation_sharing_plot_template(
    config_file, summary_tsv, phasing_tsv, parental_rates_tsv,
    output_path, log_path,
    scripts, python_conda_env, account,
):
    out_dir = f'{output_path}/mutation_sharing'
    out_pdf = f'{out_dir}/species_mutation_sharing_counts.pdf'
    done_file = f'{log_path}/mutation_sharing_plot.DONE'

    inputs = {
        'config': config_file,
        'summary': summary_tsv,
        'phasing': phasing_tsv,
        'parental_rates': parental_rates_tsv,
    }
    outputs = {
        'plot': out_pdf,
        'log': done_file,
    }
    options = {
        'cores': 1,
        'memory': '2g',
        'walltime': '00:15:00',
        'account': account,
    }

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {python_conda_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/plot_phylo_aligned_mutation_counts.py \\
    --config {config_file} \\
    --summary_tsv {summary_tsv} \\
    --phasing_tsv {phasing_tsv} \\
    --parental_rates_tsv {parental_rates_tsv} \\
    --output_pdf {out_pdf}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def parental_rate_summary_template(
    config_file, dnm_tsv, callable_dir, output_path, log_path,
    scripts, python_conda_env, account,
):
    out_dir = f'{output_path}/mutation_sharing'
    out_tsv = f'{out_dir}/subsocial_parental_mutation_rates.tsv'
    done_file = f'{log_path}/parental_rate_summary.DONE'

    inputs = {
        'config': config_file,
        'dnm': dnm_tsv,
        'callable_dir': callable_dir,
    }
    outputs = {
        'summary': out_tsv,
        'log': done_file,
    }
    options = {
        'cores': 1,
        'memory': '2g',
        'walltime': '00:30:00',
        'account': account,
    }

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {python_conda_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/summarize_parental_rates.py \\
    --config {config_file} \\
    --dnm_tsv {dnm_tsv} \\
    --callable_dir {callable_dir} \\
    --output_tsv {out_tsv}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def mutation_phasing_summary_template(
    config_file, dnm_tsv, output_path, log_path,
    scripts, python_conda_env, account,
):
    out_dir = f'{output_path}/mutation_sharing'
    out_tsv = f'{out_dir}/species_mutation_phasing_counts.tsv'
    done_file = f'{log_path}/mutation_phasing_summary.DONE'

    inputs = {
        'config': config_file,
        'dnm': dnm_tsv,
    }
    outputs = {
        'summary': out_tsv,
        'log': done_file,
    }
    options = {
        'cores': 1,
        'memory': '2g',
        'walltime': '00:15:00',
        'account': account,
    }

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {python_conda_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/summarize_mutation_phasing.py \\
    --config {config_file} \\
    --dnm_tsv {dnm_tsv} \\
    --output_tsv {out_tsv}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
