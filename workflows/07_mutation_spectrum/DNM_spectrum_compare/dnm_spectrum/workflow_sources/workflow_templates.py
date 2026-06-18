#!/usr/bin/env python3
"""
workflow_templates.py
GWF target templates for the DNM spectrum analysis.
"""
import os
from gwf import AnonymousTarget


def index_genome_template(genome_fa, SP, log_path, samtools_conda_env, account):
    """Create .fai index for a genome FASTA."""
    fai_file = genome_fa + '.fai'
    done_file = f'{log_path}/{SP}_index_genome.DONE'

    inputs = {'genome': genome_fa}
    outputs = {'fai': fai_file, 'log': done_file}
    options = {
        'cores': 1, 'memory': '4g', 'walltime': '00:30:00',
        'account': account,
    }

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {samtools_conda_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {log_path}

samtools faidx {genome_fa}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def extract_context_template(
    dnm_file, genome_fa, SP, output_path, log_path, scripts,
    x_chromosomes, exclude_trios, samtools_conda_env, python_conda_env, account
):
    """Extract trinucleotide context and assign SBS96 categories for one species."""
    out_dir = f'{output_path}/annotated_dnm/{SP}'
    out_tsv = f'{out_dir}/{SP}_dnm_annotated.tsv'
    done_file = f'{log_path}/{SP}_extract_context.DONE'
    index_done = f'{log_path}/{SP}_index_genome.DONE'

    inputs = {'dnm_file': dnm_file, 'index_done': index_done}
    outputs = {'annotated': out_tsv, 'log': done_file}
    options = {
        'cores': 1, 'memory': '4g', 'walltime': '00:30:00',
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
conda activate {samtools_conda_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/extract_context.py \\
    --dnm_file {dnm_file} \\
    --genome {genome_fa} \\
    --output {out_tsv} \\
    {x_flag} \\
    {excl_flag}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def compute_somatic_spectrum_template(
    somatic_dnm_file, SP, output_path, log_path,
    scripts, python_conda_env, account
):
    """Build a 96-category somatic SBS spectrum locally from the raw
    *_final_mutation_table_no_clusters.tsv (already strand-collapsed)."""
    out_dir = f'{output_path}/somatic_spectrum/{SP}'
    out_tsv = f'{out_dir}/{SP}_somatic_spectrum.tsv'
    done_file = f'{log_path}/{SP}_somatic_spectrum.DONE'

    inputs = {'somatic_dnm': somatic_dnm_file}
    outputs = {'spectrum': out_tsv, 'log': done_file}
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

python3 {scripts}/compute_somatic_spectrum.py \\
    --input {somatic_dnm_file} \\
    --species {SP} \\
    --output {out_tsv}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def compute_spectrum_template(
    input_files, chrom_group, group_name, output_path, log_path,
    scripts, dependency_done_files, python_conda_env, account
):
    """Compute SBS96 spectrum for a group of species and chromosome group."""
    out_dir = f'{output_path}/spectrum/{group_name}'
    out_tsv = f'{out_dir}/{group_name}_{chrom_group}_spectrum.tsv'
    done_file = f'{log_path}/{group_name}_{chrom_group}_spectrum.DONE'

    inputs = {f'done_{i}': d for i, d in enumerate(dependency_done_files)}
    outputs = {'spectrum': out_tsv, 'log': done_file}
    options = {
        'cores': 1, 'memory': '4g', 'walltime': '00:30:00',
        'account': account,
    }

    files_arg = ' '.join(input_files)

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {python_conda_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/compute_spectrum.py \\
    --input_files {files_arg} \\
    --chrom_group {chrom_group} \\
    --output {out_tsv}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def plot_spectrum_template(
    spectrum_tsv, group_name, chrom_group, output_path, log_path,
    scripts, python_conda_env, account
):
    """Render the SBS96 bar plot for a spectrum TSV."""
    out_dir = f'{output_path}/spectrum/{group_name}'
    out_pdf = f'{out_dir}/{group_name}_{chrom_group}_spectrum.pdf'
    done_file = f'{log_path}/{group_name}_{chrom_group}_plot.DONE'
    spectrum_done = f'{log_path}/{group_name}_{chrom_group}_spectrum.DONE'

    inputs = {'spectrum': spectrum_tsv, 'spectrum_done': spectrum_done}
    outputs = {'plot': out_pdf, 'log': done_file}
    options = {
        'cores': 1, 'memory': '4g', 'walltime': '00:30:00',
        'account': account,
    }

    title = f'{group_name} {chrom_group}'

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {python_conda_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/plot_spectrum.py \\
    --input {spectrum_tsv} \\
    --output {out_pdf} \\
    --title "{title}"

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def plot_comparison_template(
    germline_spectrum_tsv, somatic_tsvs, group_name,
    output_path, log_path, scripts, python_conda_env, account
):
    """Render germline (autosome) vs merged-somatic SBS96 comparison for a group."""
    out_dir = f'{output_path}/comparison/{group_name}'
    out_pdf = f'{out_dir}/{group_name}_germline_vs_somatic.pdf'
    out_tsv = f'{out_dir}/{group_name}_germline_vs_somatic_test.tsv'
    done_file = f'{log_path}/{group_name}_comparison_plot.DONE'
    germline_done = f'{log_path}/{group_name}_autosome_spectrum.DONE'

    inputs = {
        'germline': germline_spectrum_tsv,
        'germline_done': germline_done,
    }
    for i, s in enumerate(somatic_tsvs):
        inputs[f'somatic_{i}'] = s
    outputs = {'plot': out_pdf, 'tsv': out_tsv, 'log': done_file}
    options = {
        'cores': 1, 'memory': '1g', 'walltime': '00:05:00',
        'account': account,
    }

    somatic_arg = ' '.join(somatic_tsvs)

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {python_conda_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/plot_comparison.py \\
    --germline {germline_spectrum_tsv} \\
    --somatic {somatic_arg} \\
    --group_name {group_name} \\
    --output {out_pdf} \\
    --output_tsv {out_tsv}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def compare_groups_test_template(
    germline_social_tsv, germline_subsocial_tsv,
    somatic_social_tsvs, somatic_subsocial_tsvs,
    output_path, log_path, scripts, python_conda_env, account
):
    """Chi-square + Fisher's exact test across 7 collapsed categories,
    comparing social vs subsocial species groups, for germline and somatic."""
    out_dir = f'{output_path}/tests/social_vs_subsocial'
    out_tsv = f'{out_dir}/social_vs_subsocial_7category_test.tsv'
    out_pdf = f'{out_dir}/social_vs_subsocial_7category_test.pdf'
    done_file = f'{log_path}/social_vs_subsocial_7category_test.DONE'
    soc_germ_done = f'{log_path}/social_autosome_spectrum.DONE'
    sub_germ_done = f'{log_path}/subsocial_autosome_spectrum.DONE'

    inputs = {
        'germline_social': germline_social_tsv,
        'germline_subsocial': germline_subsocial_tsv,
        'soc_germ_done': soc_germ_done,
        'sub_germ_done': sub_germ_done,
    }
    for i, s in enumerate(somatic_social_tsvs):
        inputs[f'somatic_social_{i}'] = s
    for i, s in enumerate(somatic_subsocial_tsvs):
        inputs[f'somatic_subsocial_{i}'] = s
    outputs = {'tsv': out_tsv, 'pdf': out_pdf, 'log': done_file}
    options = {
        'cores': 1, 'memory': '4g', 'walltime': '00:30:00',
        'account': account,
    }

    soc_somatic_arg = ' '.join(somatic_social_tsvs)
    sub_somatic_arg = ' '.join(somatic_subsocial_tsvs)

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {python_conda_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/compare_groups_test.py \\
    --germline_social {germline_social_tsv} \\
    --germline_subsocial {germline_subsocial_tsv} \\
    --somatic_social {soc_somatic_arg} \\
    --somatic_subsocial {sub_somatic_arg} \\
    --output_tsv {out_tsv} \\
    --output_pdf {out_pdf}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def pairwise_somatic_test_template(
    pairs, output_path, log_path, scripts, python_conda_env, account
):
    """Per-pair social vs subsocial somatic spectrum test + CMH across pairs.

    pairs: list of dicts {name, social_path, subsocial_path}
    """
    out_dir = f'{output_path}/tests/pairwise_somatic'
    out_tsv = f'{out_dir}/pairwise_somatic_test.tsv'
    out_pdf = f'{out_dir}/pairwise_somatic_test.pdf'
    done_file = f'{log_path}/pairwise_somatic_test.DONE'

    inputs = {}
    for p in pairs:
        inputs[f"{p['name']}_social"] = p['social_path']
        inputs[f"{p['name']}_subsocial"] = p['subsocial_path']
    outputs = {'tsv': out_tsv, 'pdf': out_pdf, 'log': done_file}
    options = {
        'cores': 1, 'memory': '4g', 'walltime': '00:30:00',
        'account': account,
    }

    pair_args = ' '.join(
        f"--pair {p['name']} {p['social_path']} {p['subsocial_path']}"
        for p in pairs
    )

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {python_conda_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/pairwise_somatic_test.py \\
    {pair_args} \\
    --output_tsv {out_tsv} \\
    --output_pdf {out_pdf}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def spectrum_phylo_mantel_template(
    species_paths, tree_newick, output_path, log_path,
    scripts, python_conda_env, account,
):
    """Mantel test: per-species somatic spectrum cosine distance vs
    phylogenetic distance (tree topology)."""
    out_dir = f'{output_path}/tests/spectrum_phylo_mantel'
    out_tsv = f'{out_dir}/spectrum_phylo_mantel.tsv'
    out_pdf = f'{out_dir}/spectrum_phylo_mantel.pdf'
    done_file = f'{log_path}/spectrum_phylo_mantel.DONE'

    inputs = {sp: path for sp, path in species_paths}
    outputs = {'tsv': out_tsv, 'pdf': out_pdf, 'log': done_file}
    options = {
        'cores': 1, 'memory': '4g', 'walltime': '00:30:00',
        'account': account,
    }

    sp_args = ' '.join(f'--species {sp} {path}' for sp, path in species_paths)

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {python_conda_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/spectrum_phylo_mantel.py \\
    {sp_args} \\
    --tree_newick '{tree_newick}' \\
    --output_tsv {out_tsv} \\
    --output_pdf {out_pdf}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def baseline_somatic_test_template(
    social_species, subsocial_species, output_path, log_path,
    scripts, python_conda_env, account,
):
    """Test each species' somatic spectrum vs merged subsocial baseline
    (leave-one-out within subsocial)."""
    out_dir = f'{output_path}/tests/baseline_somatic'
    out_tsv = f'{out_dir}/baseline_somatic_test.tsv'
    out_pdf = f'{out_dir}/baseline_somatic_test.pdf'
    done_file = f'{log_path}/baseline_somatic_test.DONE'

    inputs = {}
    for sp, path in social_species:
        inputs[f'social_{sp}'] = path
    for sp, path in subsocial_species:
        inputs[f'subsocial_{sp}'] = path
    outputs = {'tsv': out_tsv, 'pdf': out_pdf, 'log': done_file}
    options = {
        'cores': 1, 'memory': '4g', 'walltime': '00:30:00',
        'account': account,
    }

    soc_args = ' '.join(f'--social {sp} {path}' for sp, path in social_species)
    sub_args = ' '.join(f'--subsocial {sp} {path}' for sp, path in subsocial_species)

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {python_conda_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/baseline_somatic_test.py \\
    {soc_args} \\
    {sub_args} \\
    --output_tsv {out_tsv} \\
    --output_pdf {out_pdf}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def deviation_interaction_test_template(
    germline_social_tsv, germline_subsocial_tsv,
    somatic_social_tsvs, somatic_subsocial_tsvs,
    output_path, log_path, scripts, python_conda_env, account
):
    """Per-category Breslow-Day + global 3-way log-linear interaction test of
    whether the germline-vs-somatic deviation differs between social and
    subsocial species groups."""
    out_dir = f'{output_path}/tests/social_vs_subsocial_deviation'
    out_tsv = f'{out_dir}/social_vs_subsocial_deviation_interaction_test.tsv'
    out_pdf = f'{out_dir}/social_vs_subsocial_deviation_interaction_test.pdf'
    done_file = f'{log_path}/social_vs_subsocial_deviation_interaction_test.DONE'
    soc_germ_done = f'{log_path}/social_autosome_spectrum.DONE'
    sub_germ_done = f'{log_path}/subsocial_autosome_spectrum.DONE'

    inputs = {
        'germline_social': germline_social_tsv,
        'germline_subsocial': germline_subsocial_tsv,
        'soc_germ_done': soc_germ_done,
        'sub_germ_done': sub_germ_done,
    }
    for i, s in enumerate(somatic_social_tsvs):
        inputs[f'somatic_social_{i}'] = s
    for i, s in enumerate(somatic_subsocial_tsvs):
        inputs[f'somatic_subsocial_{i}'] = s
    outputs = {'tsv': out_tsv, 'pdf': out_pdf, 'log': done_file}
    options = {
        'cores': 1, 'memory': '4g', 'walltime': '00:30:00',
        'account': account,
    }

    soc_somatic_arg = ' '.join(somatic_social_tsvs)
    sub_somatic_arg = ' '.join(somatic_subsocial_tsvs)

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {python_conda_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/deviation_interaction_test.py \\
    --germline_social {germline_social_tsv} \\
    --germline_subsocial {germline_subsocial_tsv} \\
    --somatic_social {soc_somatic_arg} \\
    --somatic_subsocial {sub_somatic_arg} \\
    --output_tsv {out_tsv} \\
    --output_pdf {out_pdf}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
