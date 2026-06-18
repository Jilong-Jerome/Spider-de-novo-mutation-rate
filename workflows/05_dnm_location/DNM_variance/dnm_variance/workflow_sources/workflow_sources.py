import yaml
import os
from workflow_templates import (
    dnm_variance_simulation_template,
    dnm_individual_variance_simulation_template,
    dnm_cv_comparison_template,
    dnm_individual_cv_comparison_template,
)


def dnm_variance_workflow(config_file, gwf):
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)

    sp = config['species_prefix']
    SP = config['project_id']

    output_dir = os.path.join(config['output_directory_path'], 'dnm_variance', SP)
    log_dir = config['log_directory_path']

    # --- Family-level variance ---
    simulations_tsv = os.path.join(output_dir, f'{sp}_variance_simulations.tsv')
    summary_txt = os.path.join(output_dir, f'{sp}_variance_summary.txt')
    done_file = os.path.join(log_dir, f'{sp}_dnm_variance_simulation.DONE')

    inputs, outputs, options, spec = dnm_variance_simulation_template(
        dnm_file=config['dnm_file'],
        callable_file=config['callable_file'],
        species_prefix=sp,
        x_chroms=config['x_chromosomes'],
        n_simulations=config.get('n_simulations', 1000),
        random_seed=config.get('random_seed', 42),
        excluded_trios=config.get('excluded_trios', []),
        output_simulations=simulations_tsv,
        output_summary=summary_txt,
        done_file=done_file,
        scripts_path=config['scripts_path'],
        output_dir=output_dir,
        account=config['account'],
    )

    gwf.target(
        name=f'{SP}_dnm_variance_simulation',
        inputs=inputs,
        outputs=outputs,
        **options,
    ) << spec

    # --- Individual-level variance ---
    ind_simulations_tsv = os.path.join(output_dir, f'{sp}_individual_variance_simulations.tsv')
    ind_summary_txt = os.path.join(output_dir, f'{sp}_individual_variance_summary.txt')
    ind_done_file = os.path.join(log_dir, f'{sp}_dnm_individual_variance_simulation.DONE')

    inputs, outputs, options, spec = dnm_individual_variance_simulation_template(
        dnm_file=config['dnm_file'],
        callable_file=config['callable_file'],
        species_prefix=sp,
        x_chroms=config['x_chromosomes'],
        n_simulations=config.get('n_simulations', 1000),
        random_seed=config.get('random_seed', 42),
        excluded_trios=config.get('excluded_trios', []),
        output_simulations=ind_simulations_tsv,
        output_summary=ind_summary_txt,
        done_file=ind_done_file,
        scripts_path=config['scripts_path'],
        output_dir=output_dir,
        account=config['account'],
    )

    gwf.target(
        name=f'{SP}_dnm_individual_variance_simulation',
        inputs=inputs,
        outputs=outputs,
        **options,
    ) << spec

    return gwf


def dnm_cv_comparison_workflow(config_files, cv_config_file, gwf):
    # Load CV comparison config (groupings, simulation params)
    with open(cv_config_file) as f:
        cv_config = yaml.safe_load(f)

    subsocial = cv_config['subsocial']
    social = cv_config['social']
    sister_pairs = cv_config['sister_pairs']
    n_simulations = cv_config.get('n_simulations', 10000)
    random_seed = cv_config.get('random_seed', 42)
    account = cv_config['account']

    # Load per-species configs to gather input data file paths and shared paths
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
        scripts_path = config['scripts_path']
        output_directory_path = config['output_directory_path']
        log_directory_path = config['log_directory_path']

    output_dir = os.path.join(output_directory_path, 'dnm_cv_comparison')
    output_tsv = os.path.join(output_dir, 'dnm_cv_comparison.tsv')
    output_summary = os.path.join(output_dir, 'dnm_cv_comparison_summary.txt')
    done_file = os.path.join(log_directory_path, 'dnm_cv_comparison.DONE')

    # Pass absolute config file paths (config_files from workflow.py are relative)
    abs_config_files = [os.path.abspath(cf) for cf in config_files]

    inputs, outputs, options, spec = dnm_cv_comparison_template(
        config_files=abs_config_files,
        subsocial=subsocial,
        social=social,
        sister_pairs=sister_pairs,
        n_simulations=n_simulations,
        random_seed=random_seed,
        output_tsv=output_tsv,
        output_summary=output_summary,
        done_file=done_file,
        scripts_path=scripts_path,
        output_dir=output_dir,
        account=account,
        all_dnm_files=all_dnm_files,
        all_callable_files=all_callable_files,
    )

    gwf.target(
        name='dnm_cv_comparison',
        inputs=inputs,
        outputs=outputs,
        **options,
    ) << spec

    return gwf


def dnm_individual_cv_comparison_workflow(config_files, cv_config_file, gwf):
    # Load CV comparison config (groupings, simulation params)
    with open(cv_config_file) as f:
        cv_config = yaml.safe_load(f)

    subsocial = cv_config['subsocial']
    social = cv_config['social']
    sister_pairs = cv_config['sister_pairs']
    n_simulations = cv_config.get('n_simulations', 10000)
    random_seed = cv_config.get('random_seed', 42)
    account = cv_config['account']

    # Load per-species configs to gather input data file paths and shared paths
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
        scripts_path = config['scripts_path']
        output_directory_path = config['output_directory_path']
        log_directory_path = config['log_directory_path']

    output_dir = os.path.join(output_directory_path, 'dnm_cv_comparison')
    output_tsv = os.path.join(output_dir, 'dnm_individual_cv_comparison.tsv')
    output_summary = os.path.join(output_dir, 'dnm_individual_cv_comparison_summary.txt')
    done_file = os.path.join(log_directory_path, 'dnm_individual_cv_comparison.DONE')

    abs_config_files = [os.path.abspath(cf) for cf in config_files]

    inputs, outputs, options, spec = dnm_individual_cv_comparison_template(
        config_files=abs_config_files,
        subsocial=subsocial,
        social=social,
        sister_pairs=sister_pairs,
        n_simulations=n_simulations,
        random_seed=random_seed,
        output_tsv=output_tsv,
        output_summary=output_summary,
        done_file=done_file,
        scripts_path=scripts_path,
        output_dir=output_dir,
        account=account,
        all_dnm_files=all_dnm_files,
        all_callable_files=all_callable_files,
    )

    gwf.target(
        name='dnm_individual_cv_comparison',
        inputs=inputs,
        outputs=outputs,
        **options,
    ) << spec

    return gwf
