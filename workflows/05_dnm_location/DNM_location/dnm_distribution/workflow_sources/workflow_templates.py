#!/usr/bin/env python3
"""
workflow_templates.py
gwf template functions for the DNM location distribution analysis.

Output hierarchy under steps/:
  window_callable/{SP}/        — per-offspring TSVs and merged species TSV
  callable_bed/{SP}/           — per-offspring callable BED segment TSVs
  dnm_distribution_test/{SP}/  — chi-squared results, per-ind counts, summary
  plots/{SP}/                  — PDF figures
"""
import os
from gwf import AnonymousTarget


def compute_window_callable_template(
    offspring_dir,
    offspring,
    sp,
    SP,
    output_path,
    log_path,
    scripts,
    genome_fai,
    window_size,
    x_chroms,
):
    """
    Count callable sites per 100 Mb window for a single offspring.
    Parses all *_callable_sites.vcf files in chrom_split/.
    """
    vcf_dir    = os.path.join(offspring_dir, 'chrom_split')
    out_dir    = f"{output_path}/window_callable/{SP}"
    out_tsv    = f"{out_dir}/{sp}_{offspring}_window_callable.tsv"
    done_file  = f"{log_path}/{sp}_{offspring}_window_callable.DONE"
    x_arg      = ' '.join(x_chroms)

    inputs  = {"callable_genome_dir": vcf_dir}
    outputs = {"window_counts": out_tsv, "log": done_file}
    options = {
        "cores": 1, "memory": "16g", "walltime": "08:00:00",
        "account": "spider2",
    }

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate python_phylo

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/compute_window_callable.py \\
    --vcf_dir {vcf_dir} \\
    --offspring {offspring} \\
    --genome_fai {genome_fai} \\
    --window_size {window_size} \\
    --x_chroms {x_arg} \\
    --output {out_tsv}

echo "FINISH: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def merge_window_callable_template(
    offspring_done_files,
    sp,
    SP,
    output_path,
    log_path,
    scripts,
):
    """
    Merge per-offspring callable window counts into a species-level table.
    Depends on ALL per-offspring DONE sentinels.
    """
    out_dir    = f"{output_path}/window_callable/{SP}"
    merged_tsv = f"{out_dir}/{sp}_window_callable.tsv"
    done_file  = f"{log_path}/{sp}_merge_window_callable.DONE"

    inputs  = {f"offspring_{i}_done": f for i, f in enumerate(offspring_done_files)}
    outputs = {"merged": merged_tsv, "log": done_file}
    options = {
        "cores": 1, "memory": "4g", "walltime": "00:30:00",
        "account": "spider2",
    }

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate python_phylo

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/merge_window_callable.py \\
    --input_dir {out_dir} \\
    --species {sp} \\
    --output {merged_tsv}

echo "FINISH: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def dnm_distribution_test_template(
    merged_tsv,
    merge_done,
    dnm_file,
    sp,
    SP,
    output_path,
    log_path,
    scripts,
    window_size,
    x_chroms,
):
    """
    Chi-squared GOF test: are DNMs proportional to callable coverage?
    """
    out_dir     = f"{output_path}/dnm_distribution_test/{SP}"
    result_tsv  = f"{out_dir}/{sp}_dnm_distribution_test.tsv"
    per_ind_tsv = f"{out_dir}/{sp}_dnm_per_individual_per_window.tsv"
    summary_txt = f"{out_dir}/{sp}_dnm_distribution_test_summary.txt"
    done_file   = f"{log_path}/{sp}_dnm_distribution_test.DONE"
    x_arg       = ' '.join(x_chroms)

    inputs  = {"window_callable": merged_tsv, "merge_done": merge_done}
    outputs = {
        "result_tsv":         result_tsv,
        "per_ind_counts_tsv": per_ind_tsv,
        "result_summary":     summary_txt,
        "log":                done_file,
    }
    options = {
        "cores": 1, "memory": "4g", "walltime": "00:30:00",
        "account": "spider2",
    }

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate python_phylo

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/dnm_distribution_test.py \\
    --window_callable {merged_tsv} \\
    --dnm_file {dnm_file} \\
    --species {sp} \\
    --window_size {window_size} \\
    --x_chroms {x_arg} \\
    --output_tsv {result_tsv} \\
    --output_per_ind {per_ind_tsv} \\
    --output_summary {summary_txt}

echo "FINISH: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def plot_dnm_genome_template(
    dnm_file,
    genome_fai,
    sp,
    SP,
    output_path,
    log_path,
    scripts,
    x_chroms,
):
    """
    Genome-wide DNM tick-mark visualisation (one subplot per autosome).
    """
    out_dir   = f"{output_path}/plots/{SP}"
    out_pdf   = f"{out_dir}/{sp}_dnm_genome.pdf"
    done_file = f"{log_path}/{sp}_plot_dnm_genome.DONE"
    x_arg     = ' '.join(x_chroms)

    inputs  = {"dnm_file": dnm_file}
    outputs = {"plot": out_pdf, "log": done_file}
    options = {
        "cores": 1, "memory": "4g", "walltime": "00:30:00",
        "account": "spider2",
    }

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate python_phylo

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/plot_dnm_genome.py \\
    --dnm_file {dnm_file} \\
    --genome_fai {genome_fai} \\
    --species {sp} \\
    --x_chroms {x_arg} \\
    --output {out_pdf}

echo "FINISH: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def compute_callable_bed_template(
    offspring_dir,
    offspring,
    sp,
    SP,
    output_path,
    log_path,
    scripts,
    genome_fai,
    x_chroms,
):
    """
    Convert per-offspring callable VCF files into merged BED segments.
    Output: one TSV per offspring with columns offspring/chrom/start/end.
    """
    vcf_dir   = os.path.join(offspring_dir, 'chrom_split')
    out_dir   = f"{output_path}/callable_bed/{SP}"
    out_bed   = f"{out_dir}/{sp}_{offspring}_callable.bed"
    done_file = f"{log_path}/{sp}_{offspring}_callable_bed.DONE"
    x_arg     = ' '.join(x_chroms)

    inputs  = {"callable_genome_dir": vcf_dir}
    outputs = {"bed": out_bed, "log": done_file}
    options = {
        "cores": 1, "memory": "16g", "walltime": "08:00:00",
        "account": "spider2",
    }

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate python_phylo

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/compute_callable_bed.py \\
    --vcf_dir {vcf_dir} \\
    --offspring {offspring} \\
    --genome_fai {genome_fai} \\
    --x_chroms {x_arg} \\
    --output {out_bed}

echo "FINISH: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def plot_callable_segments_template(
    bed_done_files,
    sp,
    SP,
    output_path,
    log_path,
    scripts,
    genome_fai,
    dnm_file,
    x_chroms,
    bridge_gap=20000,
):
    """
    Combined callable segment + de novo mutation visualisation.
    One panel per autosome; per-individual callable tracks with DNM dots overlaid.
    Depends on ALL per-offspring callable BED DONE sentinels.
    """
    bed_dir   = f"{output_path}/callable_bed/{SP}"
    out_dir   = f"{output_path}/plots/{SP}"
    out_pdf   = f"{out_dir}/{sp}_callable_segments.pdf"
    done_file = f"{log_path}/{sp}_plot_callable_segments.DONE"
    x_arg     = ' '.join(x_chroms)

    inputs  = {f"offspring_{i}_bed_done": f for i, f in enumerate(bed_done_files)}
    inputs["dnm_file"] = dnm_file
    outputs = {"plot": out_pdf, "log": done_file}
    options = {
        "cores": 1, "memory": "32g", "walltime": "01:00:00",
        "account": "spider2",
    }

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate python_phylo

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/plot_callable_segments.py \\
    --bed_dir {bed_dir} \\
    --species_prefix {sp} \\
    --genome_fai {genome_fai} \\
    --dnm_file {dnm_file} \\
    --x_chroms {x_arg} \\
    --bridge_gap {bridge_gap} \\
    --output {out_pdf}

echo "FINISH: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def nearest_neighbor_test_template(
    bed_done_files,
    sp,
    SP,
    output_path,
    log_path,
    scripts,
    dnm_file,
    x_chroms,
    n_simulations=1000,
    seed=42,
):
    """
    Permutation test: average nearest-neighbour distance of DNMs vs. random
    draws from the union of all individuals' callable genome.
    Depends on ALL per-offspring callable BED DONE sentinels.
    """
    bed_dir     = f"{output_path}/callable_bed/{SP}"
    out_dir     = f"{output_path}/nearest_neighbor_test/{SP}"
    out_tsv     = f"{out_dir}/{sp}_nearest_neighbor_simulations.tsv"
    out_summary = f"{out_dir}/{sp}_nearest_neighbor_summary.txt"
    out_plot    = f"{out_dir}/{sp}_nearest_neighbor_plot.pdf"
    done_file   = f"{log_path}/{sp}_nearest_neighbor_test.DONE"
    x_arg       = ' '.join(x_chroms)

    inputs  = {f"offspring_{i}_bed_done": f for i, f in enumerate(bed_done_files)}
    inputs["dnm_file"] = dnm_file
    outputs = {
        "simulations_tsv": out_tsv,
        "summary":         out_summary,
        "plot":            out_plot,
        "log":             done_file,
    }
    options = {
        "cores": 1, "memory": "32g", "walltime": "02:00:00",
        "account": "spider2",
    }

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate python_phylo

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/nearest_neighbor_test.py \\
    --bed_dir {bed_dir} \\
    --species_prefix {sp} \\
    --dnm_file {dnm_file} \\
    --x_chroms {x_arg} \\
    --n_simulations {n_simulations} \\
    --seed {seed} \\
    --output_tsv {out_tsv} \\
    --output_summary {out_summary} \\
    --output_plot {out_plot}

echo "FINISH: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def plot_nearest_neighbor_all_species_template(
    species_results,
    output_path,
    log_path,
    scripts,
):
    """
    Multi-panel nearest-neighbour result plot across all species.
    Depends on all per-species nearest-neighbour summaries and simulation TSVs.
    """
    out_dir   = f"{output_path}/nearest_neighbor_test"
    out_pdf   = f"{out_dir}/all_species_nearest_neighbor_plot.pdf"
    done_file = f"{log_path}/all_species_nearest_neighbor_plot.DONE"

    inputs = {}
    for species in species_results:
        SP = species['SP']
        inputs[f"{SP}_simulations"] = species['simulations_tsv']
        inputs[f"{SP}_summary"] = species['summary']

    outputs = {"plot": out_pdf, "log": done_file}
    options = {
        "cores": 1, "memory": "8g", "walltime": "00:30:00",
        "account": "spider2",
    }

    spec = f"""
# Setting conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate python_phylo

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

mkdir -p {out_dir} {log_path}

python3 {scripts}/plot_nearest_neighbor_all_species.py \\
    --input_dir {out_dir} \\
    --output {out_pdf}

echo "FINISH: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
