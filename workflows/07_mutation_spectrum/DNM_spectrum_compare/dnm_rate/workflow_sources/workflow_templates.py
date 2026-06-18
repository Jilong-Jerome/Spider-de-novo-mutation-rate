#!/usr/bin/env python3
"""
workflow_templates.py
GWF target templates for the DNM mutation-rate analysis.
"""
from gwf import AnonymousTarget


def callable_vcf_trinuc_template(
    vcf_file, genome_fa, SP, offspring_id, chrom,
    output_path, log_path, scripts, samtools_conda_env, account,
):
    """Per (species, offspring, chrom) callable VCF -> 32-trinuc count TSV."""
    out_dir = f'{output_path}/rate/germline_callable_trinuc/{SP}/{offspring_id}'
    out_tsv = f'{out_dir}/{offspring_id}_{chrom}_trinuc.tsv'
    done_file = f'{log_path}/{SP}_{offspring_id}_{chrom}_callable_trinuc.DONE'
    index_done = f'{log_path}/{SP}_index_genome.DONE'

    inputs = {'vcf': vcf_file, 'index_done': index_done}
    outputs = {'trinuc': out_tsv, 'log': done_file}
    options = {
        'cores': 1, 'memory': '4g', 'walltime': '02:00:00',
        'account': account,
    }

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {samtools_conda_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/callable_vcf_trinuc.py \\
    --vcf {vcf_file} \\
    --genome {genome_fa} \\
    --chrom {chrom} \\
    --output {out_tsv}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def aggregate_callable_trinuc_template(
    SP, offspring_chrom_pairs, output_path, log_path,
    scripts, python_conda_env, account,
):
    """Sum per-(offspring, chrom) trinuc TSVs into a single species-level table."""
    out_dir = f'{output_path}/rate/germline_callable_trinuc/{SP}'
    out_tsv = f'{out_dir}/{SP}_autosome_callable_trinuc.tsv'
    done_file = f'{log_path}/{SP}_germline_callable_trinuc.DONE'

    inputs = {}
    input_paths = []
    for offspring_id, chrom in offspring_chrom_pairs:
        p = f'{out_dir}/{offspring_id}/{offspring_id}_{chrom}_trinuc.tsv'
        inputs[f'{offspring_id}_{chrom}'] = p
        input_paths.append(p)
    outputs = {'aggregated': out_tsv, 'log': done_file}
    options = {
        'cores': 1, 'memory': '2g', 'walltime': '00:30:00',
        'account': account,
    }

    files_arg = ' '.join(input_paths)

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {python_conda_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/aggregate_callable_trinuc.py \\
    --species {SP} \\
    --output {out_tsv} \\
    --input_files {files_arg}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def germline_rate_summary_template(
    SP, dnm_file, callable_file, callable_trinuc_tsv,
    x_chromosomes, exclude_trios,
    output_path, log_path, scripts, python_conda_env, account,
):
    """Per-species germline rate summary (overall + 7 mutation classes)."""
    out_dir = f'{output_path}/rate/germline/{SP}'
    out_tsv = f'{out_dir}/{SP}_germline_rate_summary.tsv'
    done_file = f'{log_path}/{SP}_germline_rate_summary.DONE'

    inputs = {
        'dnm': dnm_file,
        'callable_file': callable_file,
        'callable_trinuc': callable_trinuc_tsv,
    }
    outputs = {'summary': out_tsv, 'log': done_file}
    options = {
        'cores': 1, 'memory': '2g', 'walltime': '00:15:00',
        'account': account,
    }

    x_arg = ' '.join(x_chromosomes) if x_chromosomes else ''
    x_flag = f'--x_chromosomes {x_arg}' if x_arg else ''
    excl_arg = ' '.join(exclude_trios) if exclude_trios else ''
    excl_flag = f'--exclude_trios {excl_arg}' if excl_arg else ''

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {python_conda_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/germline_rate_summary.py \\
    --species {SP} \\
    --dnm_file {dnm_file} \\
    --callable_file {callable_file} \\
    --callable_trinuc {callable_trinuc_tsv} \\
    {x_flag} \\
    {excl_flag} \\
    --output {out_tsv}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def somatic_rate_summary_template(
    SP, somatic_dnm_file, somatic_trinuc_counts_file,
    output_path, log_path, scripts, python_conda_env, account,
):
    """Per-species somatic rate summary (overall + 7 mutation classes)."""
    out_dir = f'{output_path}/rate/somatic/{SP}'
    out_tsv = f'{out_dir}/{SP}_somatic_rate_summary.tsv'
    done_file = f'{log_path}/{SP}_somatic_rate_summary.DONE'

    inputs = {
        'somatic_dnm': somatic_dnm_file,
        'trinuc_counts': somatic_trinuc_counts_file,
    }
    outputs = {'summary': out_tsv, 'log': done_file}
    options = {
        'cores': 1, 'memory': '2g', 'walltime': '00:15:00',
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

python3 {scripts}/somatic_rate_summary.py \\
    --species {SP} \\
    --somatic_dnm {somatic_dnm_file} \\
    --trinuc_counts {somatic_trinuc_counts_file} \\
    --output {out_tsv}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def aggregate_rates_template(
    species_list, spectrum_config_path, output_path, log_path,
    scripts, python_conda_env, account,
):
    """Combine all per-species summaries + per-(offspring, chrom) trinuc tables
    into the master rate TSV, per-trio germline TSV, per-species somatic TSV,
    and somatic-vs-germline fold-change TSV. Reads species/groups from the
    spectrum config to stay a single source of truth."""
    out_dir = f'{output_path}/rate/master'
    master_tsv = f'{out_dir}/mutation_rate_master.tsv'
    per_trio_tsv = f'{out_dir}/per_trio_germline_mutations_callable.tsv'
    per_species_somatic_tsv = f'{out_dir}/per_species_somatic_mutations_callable.tsv'
    fold_change_tsv = f'{out_dir}/somatic_vs_germline_fold_change.tsv'
    fold_change_test_tsv = f'{out_dir}/social_vs_subsocial_fold_change_test.tsv'
    master_96_tsv = f'{out_dir}/mutation_rate_96_context_master.tsv'
    fold_change_96_tsv = f'{out_dir}/somatic_vs_germline_96_context_fold_change.tsv'
    done_file = f'{log_path}/mutation_rate_aggregate.DONE'

    inputs = {'spectrum_config': spectrum_config_path}
    for SP in species_list:
        inputs[f'{SP}_germline'] = f'{output_path}/rate/germline/{SP}/{SP}_germline_rate_summary.tsv'
        inputs[f'{SP}_somatic'] = f'{output_path}/rate/somatic/{SP}/{SP}_somatic_rate_summary.tsv'
        inputs[f'{SP}_germline_trinuc'] = f'{output_path}/rate/germline_callable_trinuc/{SP}/{SP}_autosome_callable_trinuc.tsv'

    outputs = {
        'master': master_tsv,
        'per_trio': per_trio_tsv,
        'per_species_somatic': per_species_somatic_tsv,
        'fold_change': fold_change_tsv,
        'fold_change_test': fold_change_test_tsv,
        'master_96': master_96_tsv,
        'fold_change_96': fold_change_96_tsv,
        'log': done_file,
    }
    options = {
        'cores': 1, 'memory': '4g', 'walltime': '00:30:00',
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

python3 {scripts}/aggregate_rates.py \\
    --spectrum_config {spectrum_config_path} \\
    --output_path {output_path} \\
    --output_master {master_tsv} \\
    --output_per_trio {per_trio_tsv} \\
    --output_per_species_somatic {per_species_somatic_tsv} \\
    --output_fold_change {fold_change_tsv} \\
    --output_fold_change_test {fold_change_test_tsv} \\
    --output_master_96 {master_96_tsv} \\
    --output_fold_change_96 {fold_change_96_tsv}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def plot_rates_template(
    output_path, log_path, scripts, python_conda_env, account,
):
    """Render rate plots from the master TSV + fold-change TSV."""
    in_dir = f'{output_path}/rate/master'
    out_dir = f'{output_path}/rate/plots'
    master_tsv = f'{in_dir}/mutation_rate_master.tsv'
    fold_change_tsv = f'{in_dir}/somatic_vs_germline_fold_change.tsv'
    fold_change_test_tsv = f'{in_dir}/social_vs_subsocial_fold_change_test.tsv'
    master_96_tsv = f'{in_dir}/mutation_rate_96_context_master.tsv'
    fold_change_96_tsv = f'{in_dir}/somatic_vs_germline_96_context_fold_change.tsv'
    pdf_per_species = f'{out_dir}/per_species_rates.pdf'
    pdf_group = f'{out_dir}/group_rates.pdf'
    pdf_fc = f'{out_dir}/fold_change.pdf'
    pdf_per_species_96 = f'{out_dir}/per_species_96_context_rates.pdf'
    pdf_group_96 = f'{out_dir}/group_96_context_rates.pdf'
    pdf_fc_96 = f'{out_dir}/fold_change_96_context.pdf'
    done_file = f'{log_path}/mutation_rate_plot.DONE'
    aggregate_done = f'{log_path}/mutation_rate_aggregate.DONE'

    inputs = {
        'master': master_tsv,
        'fold_change': fold_change_tsv,
        'fold_change_test': fold_change_test_tsv,
        'master_96': master_96_tsv,
        'fold_change_96': fold_change_96_tsv,
        'aggregate_done': aggregate_done,
    }
    outputs = {
        'per_species': pdf_per_species,
        'group': pdf_group,
        'fold_change': pdf_fc,
        'per_species_96': pdf_per_species_96,
        'group_96': pdf_group_96,
        'fold_change_96': pdf_fc_96,
        'log': done_file,
    }
    options = {
        'cores': 1, 'memory': '2g', 'walltime': '00:15:00',
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

python3 {scripts}/plot_rates.py \\
    --master {master_tsv} \\
    --fold_change {fold_change_tsv} \\
    --fold_change_test {fold_change_test_tsv} \\
    --output_per_species {pdf_per_species} \\
    --output_group {pdf_group} \\
    --output_fold_change {pdf_fc} \\
    --master_96 {master_96_tsv} \\
    --fold_change_96 {fold_change_96_tsv} \\
    --output_per_species_96 {pdf_per_species_96} \\
    --output_group_96 {pdf_group_96} \\
    --output_fold_change_96 {pdf_fc_96}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def phylo_aligned_rate_plot_template(
    spectrum_config_path, output_path, log_path,
    scripts, python_conda_env, account,
    n_bootstrap=10000, seed=42,
):
    """Per-species overall germline (point + bootstrap 95% CI) and somatic
    (point) rates with somatic/germline fold change, arranged as a 4-column
    figure aligned to the species phylogeny."""
    in_dir = f'{output_path}/rate/master'
    plot_dir = f'{output_path}/rate/plots'
    per_trio_tsv = f'{in_dir}/per_trio_germline_mutations_callable.tsv'
    per_species_somatic_tsv = f'{in_dir}/per_species_somatic_mutations_callable.tsv'
    out_tsv = f'{in_dir}/phylo_aligned_overall_rates.tsv'
    out_pdf = f'{plot_dir}/phylo_aligned_overall_rates.pdf'
    done_file = f'{log_path}/phylo_aligned_rate_plot.DONE'
    aggregate_done = f'{log_path}/mutation_rate_aggregate.DONE'

    inputs = {
        'per_trio': per_trio_tsv,
        'per_species_somatic': per_species_somatic_tsv,
        'spectrum_config': spectrum_config_path,
        'aggregate_done': aggregate_done,
    }
    outputs = {'tsv': out_tsv, 'pdf': out_pdf, 'log': done_file}
    options = {
        'cores': 1, 'memory': '4g', 'walltime': '00:30:00',
        'account': account,
    }

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {python_conda_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {in_dir} {plot_dir} {log_path}

python3 {scripts}/plot_phylo_aligned_rates.py \\
    --per_trio_tsv {per_trio_tsv} \\
    --per_species_somatic_tsv {per_species_somatic_tsv} \\
    --spectrum_config {spectrum_config_path} \\
    --output_tsv {out_tsv} \\
    --output_pdf {out_pdf} \\
    --n_bootstrap {n_bootstrap} \\
    --seed {seed}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
