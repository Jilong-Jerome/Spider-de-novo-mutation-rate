#!/usr/bin/env python3
"""
GWF target templates for final visual presentation outputs.
"""
from gwf import AnonymousTarget


def combined_rate_visual_template(
    per_trio_germline_tsv,
    per_species_somatic_tsv,
    spectrum_test_tsv,
    output_path,
    log_path,
    scripts,
    python_conda_env,
    account,
    bootstrap_replicates,
    bootstrap_seed,
):
    out_dir = f'{output_path}/final_visuals/combined_rate_by_mutation_class'
    input_data_tsv = f'{out_dir}/combined_rate_input_data.tsv'
    plot_data_tsv = f'{out_dir}/combined_rate_plot_data.tsv'
    output_pdf = f'{out_dir}/combined_rate_by_mutation_class.pdf'
    done_file = f'{log_path}/final_combined_rate_by_mutation_class.DONE'

    inputs = {
        'per_trio_germline': per_trio_germline_tsv,
        'per_species_somatic': per_species_somatic_tsv,
        'spectrum_test': spectrum_test_tsv,
    }
    outputs = {
        'input_data': input_data_tsv,
        'plot_data': plot_data_tsv,
        'pdf': output_pdf,
        'log': done_file,
    }
    options = {
        'cores': 1,
        'memory': '4g',
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

python3 {scripts}/combined_rate_visual.py \\
    --per_trio_germline {per_trio_germline_tsv} \\
    --per_species_somatic {per_species_somatic_tsv} \\
    --spectrum_test {spectrum_test_tsv} \\
    --output_input_data {input_data_tsv} \\
    --output_plot_data {plot_data_tsv} \\
    --output_pdf {output_pdf} \\
    --bootstrap_replicates {bootstrap_replicates} \\
    --bootstrap_seed {bootstrap_seed}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
