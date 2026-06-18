"""
GWF AnonymousTarget definitions for the sequencing batch summary workflow.
"""

from gwf import AnonymousTarget


def summarize_individual_template(
    individual_dir,
    species,
    family,
    individual,
    safe_name,
    output_path,
    log_path,
    scripts_path,
    account,
):
    """Measure uncompressed read data per (round_id, lane) for one individual."""
    tsv_out = f'{output_path}/{safe_name}_summary.tsv'
    done_log = f'{log_path}/{safe_name}_summarize.DONE'
    script = f'{scripts_path}/summarize_individual.py'

    inputs = {'individual_dir': individual_dir}
    outputs = {
        'tsv': tsv_out,
        'done': done_log,
    }
    options = {
        'cores': 1,
        'memory': '4g',
        'walltime': '06:00:00',
        'account': account,
    }
    spec = f"""
source $(conda info --base)/etc/profile.d/conda.sh
conda activate python_phylo

python {script} \\
    --individual_dir {individual_dir} \\
    --species {species} \\
    --family {family} \\
    --individual {individual} \\
    --output {tsv_out}

touch {done_log}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def concatenate_results_template(
    partial_tsvs,
    done_logs,
    final_output,
    log_path,
    account,
):
    """Concatenate all per-individual TSVs into a single summary with a header."""
    concat_done = f'{log_path}/concatenate_results.DONE'

    inputs = {
        'tsvs': partial_tsvs,
        'done_logs': done_logs,
    }
    outputs = {
        'final': final_output,
        'done': concat_done,
    }
    options = {
        'cores': 1,
        'memory': '4g',
        'walltime': '00:30:00',
        'account': account,
    }

    # Build the shell loop over partial TSVs
    tsv_list = ' '.join(partial_tsvs)
    spec = f"""
# Write header once
echo -e "species\\tfamily\\tindividual\\tround_id\\tlane\\tsize_gb" > {final_output}

# Append each partial TSV (no header lines to skip)
for tsv in {tsv_list}; do
    cat "$tsv" >> {final_output}
done

touch {concat_done}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def add_coverage_template(
    summary_tsv,
    fai_directory,
    coverage_output,
    log_path,
    scripts_path,
    account,
):
    """Append an approximate sequencing-coverage column to the summary TSV."""
    done_log = f'{log_path}/add_coverage.DONE'
    script = f'{scripts_path}/add_coverage.py'

    inputs = {'summary_tsv': summary_tsv}
    outputs = {
        'coverage_tsv': coverage_output,
        'done': done_log,
    }
    options = {
        'cores': 1,
        'memory': '4g',
        'walltime': '00:15:00',
        'account': account,
    }
    spec = f"""
source $(conda info --base)/etc/profile.d/conda.sh
conda activate python_phylo

python {script} \\
    --input {summary_tsv} \\
    --data_dir {fai_directory} \\
    --output {coverage_output}

touch {done_log}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def plot_heatmap_template(
    coverage_tsv,
    heatmap_output,
    log_path,
    scripts_path,
    account,
):
    """Generate a per-lane coverage heatmap with sociality and species annotations."""
    done_log = f'{log_path}/plot_heatmap.DONE'
    script = f'{scripts_path}/plot_heatmap.py'

    inputs = {'coverage_tsv': coverage_tsv}
    outputs = {
        'heatmap': heatmap_output,
        'done': done_log,
    }
    options = {
        'cores': 1,
        'memory': '8g',
        'walltime': '00:30:00',
        'account': account,
    }
    spec = f"""
source $(conda info --base)/etc/profile.d/conda.sh
conda activate python_phylo

python {script} \\
    --input {coverage_tsv} \\
    --output {heatmap_output}

touch {done_log}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def plot_coverage_scatter_template(
    coverage_tsv,
    scatter_output,
    log_path,
    scripts_path,
    account,
):
    """Scatter plot of total sequencing coverage per individual, grouped by species."""
    done_log = f'{log_path}/plot_coverage_scatter.DONE'
    script = f'{scripts_path}/plot_coverage_scatter.py'

    inputs = {'coverage_tsv': coverage_tsv}
    outputs = {
        'scatter': scatter_output,
        'done': done_log,
    }
    options = {
        'cores': 1,
        'memory': '4g',
        'walltime': '00:15:00',
        'account': account,
    }
    spec = f"""
source $(conda info --base)/etc/profile.d/conda.sh
conda activate python_phylo

python {script} \\
    --input {coverage_tsv} \\
    --output {scatter_output}

touch {done_log}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
