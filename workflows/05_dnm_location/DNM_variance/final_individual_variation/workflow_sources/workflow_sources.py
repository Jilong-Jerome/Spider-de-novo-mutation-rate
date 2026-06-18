import os

import yaml

from workflow_templates import final_individual_variation_template


def final_individual_variation_workflow(config_files, cv_config_file, gwf):
    with open(cv_config_file) as f:
        cv_config = yaml.safe_load(f)

    all_dnm_files = []
    all_callable_files = []
    scripts_path = None
    output_directory_path = None
    log_directory_path = None

    for config_file in sorted(config_files):
        with open(config_file) as f:
            config = yaml.safe_load(f)
        all_dnm_files.append(config['dnm_file'])
        all_callable_files.append(config['callable_file'])
        output_directory_path = config['output_directory_path']
        log_directory_path = config['log_directory_path']

    scripts_path = os.path.abspath('./workflow_sources')
    output_dir = os.path.join(output_directory_path, 'final_individual_variation')
    done_file = os.path.join(log_directory_path, 'final_individual_variation_outputs.DONE')

    outputs = {
        'variance_tsv': os.path.join(output_dir, 'individual_variance_per_species.tsv'),
        'cv_tsv': os.path.join(output_dir, 'individual_cv_per_species.tsv'),
        'cv_simulations_tsv': os.path.join(output_dir, 'individual_cv_simulations_per_species.tsv'),
        'pairwise_tsv': os.path.join(output_dir, 'individual_cv_pairwise_tests.tsv'),
        'summary_txt': os.path.join(output_dir, 'individual_cv_pairwise_tests_summary.txt'),
        'plot_pdf': os.path.join(output_dir, 'individual_variation_tree_aligned.pdf'),
        'pairwise_plot_pdf': os.path.join(output_dir, 'individual_cv_pairwise_difference_distributions.pdf'),
        'cv_simulations_plot_pdf': os.path.join(output_dir, 'individual_cv_simulations_per_species.pdf'),
        'done_file': done_file,
    }

    inputs, target_outputs, options, spec = final_individual_variation_template(
        config_files=[os.path.abspath(cf) for cf in sorted(config_files)],
        cv_config_file=os.path.abspath(cv_config_file),
        outputs=outputs,
        scripts_path=scripts_path,
        output_dir=output_dir,
        account=cv_config['account'],
        all_dnm_files=all_dnm_files,
        all_callable_files=all_callable_files,
    )

    gwf.target(
        name='final_individual_variation_outputs',
        inputs=inputs,
        outputs=target_outputs,
        **options,
    ) << spec

    return gwf
