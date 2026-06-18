def final_individual_variation_template(
    config_files,
    cv_config_file,
    outputs,
    scripts_path,
    output_dir,
    account,
    all_dnm_files,
    all_callable_files,
):
    inputs = config_files + [cv_config_file] + all_dnm_files + all_callable_files
    target_outputs = [
        outputs['variance_tsv'],
        outputs['cv_tsv'],
        outputs['cv_simulations_tsv'],
        outputs['pairwise_tsv'],
        outputs['summary_txt'],
        outputs['plot_pdf'],
        outputs['pairwise_plot_pdf'],
        outputs['cv_simulations_plot_pdf'],
        outputs['done_file'],
    ]

    options = {
        'cores': 2,
        'memory': '6g',
        'walltime': '00:30:00',
        'account': account,
    }

    config_files_str = ' '.join(config_files)

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate python_phylo

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {output_dir} $(dirname {outputs['done_file']})

python3 {scripts_path}/make_final_individual_variation_outputs.py \\
    --config_files {config_files_str} \\
    --cv_config_file {cv_config_file} \\
    --output_variance_tsv {outputs['variance_tsv']} \\
    --output_cv_tsv {outputs['cv_tsv']} \\
    --output_cv_simulations_tsv {outputs['cv_simulations_tsv']} \\
    --output_pairwise_tsv {outputs['pairwise_tsv']} \\
    --output_summary {outputs['summary_txt']} \\
    --output_plot_pdf {outputs['plot_pdf']} \\
    --output_pairwise_plot_pdf {outputs['pairwise_plot_pdf']} \\
    --output_cv_simulations_plot_pdf {outputs['cv_simulations_plot_pdf']}

echo "FINISH: $(date)"
echo done > {outputs['done_file']}
"""

    return inputs, target_outputs, options, spec
