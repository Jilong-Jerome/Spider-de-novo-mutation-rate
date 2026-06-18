def dnm_individual_cv_comparison_template(
    config_files,
    subsocial,
    social,
    sister_pairs,
    n_simulations,
    random_seed,
    output_tsv,
    output_summary,
    done_file,
    scripts_path,
    output_dir,
    account,
    all_dnm_files,
    all_callable_files,
):
    inputs = all_dnm_files + all_callable_files
    outputs = [output_tsv, output_summary, done_file]

    options = {
        'cores': 2,
        'memory': '4g',
        'walltime': '00:30:00',
        'account': account,
    }

    config_files_str = ' '.join(config_files)
    subsocial_str = ' '.join(subsocial)
    social_str = ' '.join(social)
    sister_pairs_str = ' '.join(sister_pairs)

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate python_phylo

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {output_dir} $(dirname {done_file})

python3 {scripts_path}/dnm_individual_cv_comparison.py \\
    --config_files {config_files_str} \\
    --subsocial {subsocial_str} \\
    --social {social_str} \\
    --sister_pairs {sister_pairs_str} \\
    --n_simulations {n_simulations} \\
    --random_seed {random_seed} \\
    --output_tsv {output_tsv} \\
    --output_summary {output_summary}

echo "FINISH: $(date)"
echo done > {done_file}
"""

    return inputs, outputs, options, spec


def dnm_cv_comparison_template(
    config_files,
    subsocial,
    social,
    sister_pairs,
    n_simulations,
    random_seed,
    output_tsv,
    output_summary,
    done_file,
    scripts_path,
    output_dir,
    account,
    all_dnm_files,
    all_callable_files,
):
    inputs = all_dnm_files + all_callable_files
    outputs = [output_tsv, output_summary, done_file]

    options = {
        'cores': 2,
        'memory': '4g',
        'walltime': '00:30:00',
        'account': account,
    }

    config_files_str = ' '.join(config_files)
    subsocial_str = ' '.join(subsocial)
    social_str = ' '.join(social)
    sister_pairs_str = ' '.join(sister_pairs)

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate python_phylo

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {output_dir} $(dirname {done_file})

python3 {scripts_path}/dnm_cv_comparison.py \\
    --config_files {config_files_str} \\
    --subsocial {subsocial_str} \\
    --social {social_str} \\
    --sister_pairs {sister_pairs_str} \\
    --n_simulations {n_simulations} \\
    --random_seed {random_seed} \\
    --output_tsv {output_tsv} \\
    --output_summary {output_summary}

echo "FINISH: $(date)"
echo done > {done_file}
"""

    return inputs, outputs, options, spec


def dnm_individual_variance_simulation_template(
    dnm_file,
    callable_file,
    species_prefix,
    x_chroms,
    n_simulations,
    random_seed,
    excluded_trios,
    output_simulations,
    output_summary,
    done_file,
    scripts_path,
    output_dir,
    account,
):
    inputs = [dnm_file, callable_file]
    outputs = [output_simulations, output_summary, done_file]

    options = {
        'cores': 2,
        'memory': '8g',
        'walltime': '01:00:00',
        'account': account,
    }

    x_chroms_str = ' '.join(x_chroms)

    if excluded_trios:
        excluded_trios_arg = f' \\\n    --excluded_trios {" ".join(excluded_trios)}'
    else:
        excluded_trios_arg = ' \\\n    --excluded_trios'

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate python_phylo

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {output_dir} $(dirname {done_file})

python3 {scripts_path}/dnm_individual_variance_simulation.py \\
    --dnm_file {dnm_file} \\
    --callable_file {callable_file} \\
    --species_prefix {species_prefix} \\
    --x_chroms {x_chroms_str} \\
    --n_simulations {n_simulations} \\
    --random_seed {random_seed} \\
    --output_simulations {output_simulations} \\
    --output_summary {output_summary}{excluded_trios_arg}

echo "FINISH: $(date)"
echo done > {done_file}
"""

    return inputs, outputs, options, spec


def dnm_variance_simulation_template(
    dnm_file,
    callable_file,
    species_prefix,
    x_chroms,
    n_simulations,
    random_seed,
    excluded_trios,
    output_simulations,
    output_summary,
    done_file,
    scripts_path,
    output_dir,
    account,
):
    inputs = [dnm_file, callable_file]
    outputs = [output_simulations, output_summary, done_file]

    options = {
        'cores': 2,
        'memory': '8g',
        'walltime': '01:00:00',
        'account': account,
    }

    x_chroms_str = ' '.join(x_chroms)

    if excluded_trios:
        excluded_trios_arg = f' \\\n    --excluded_trios {" ".join(excluded_trios)}'
    else:
        excluded_trios_arg = ' \\\n    --excluded_trios'

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate python_phylo

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {output_dir} $(dirname {done_file})

python3 {scripts_path}/dnm_variance_simulation.py \\
    --dnm_file {dnm_file} \\
    --callable_file {callable_file} \\
    --species_prefix {species_prefix} \\
    --x_chroms {x_chroms_str} \\
    --n_simulations {n_simulations} \\
    --random_seed {random_seed} \\
    --output_simulations {output_simulations} \\
    --output_summary {output_summary}{excluded_trios_arg}

echo "FINISH: $(date)"
echo done > {done_file}
"""

    return inputs, outputs, options, spec
