#!/usr/bin/env python3
"""
workflow_templates.py - GWF job templates for repair gene genotype analysis

This module defines job templates for:
1. parse_gff3 - Extract gene/exon/intron coordinates
2. extract_genotypes - Extract genotypes from gVCF
3. calculate_heterozygosity - Calculate het metrics per individual
4. calculate_dxy - Calculate pairwise Dxy between pairs
5. generate_summary - Compile summary reports
"""

from gwf import AnonymousTarget
import os


def parse_gff3_template(output_path: str, script_path: str, log_path: str,
                        gff3_file: str, gene_id: str, gene_name: str,
                        species: str):
    """
    Template for extracting gene/exon/intron coordinates from GFF3.

    Args:
        output_path: Base output directory
        script_path: Path to helper scripts
        log_path: Path for log files
        gff3_file: Path to GFF3 annotation file
        gene_id: Gene/transcript ID to extract
        gene_name: Common gene name for output naming
        species: Species identifier
    """
    gene_dir = f"{output_path}/{species}/{gene_name}"

    inputs = {"gff3": gff3_file}
    outputs = {
        "gene_bed": f"{gene_dir}/gene_coords.bed",
        "exon_bed": f"{gene_dir}/exon_coords.bed",
        "intron_bed": f"{gene_dir}/intron_coords.bed",
        "log": f"{log_path}/{species}_{gene_name}_parse_gff3.DONE"
    }
    options = {
        'cores': 1,
        'memory': '2g',
        'walltime': '00:30:00',
        'account': 'spider2'
    }

    spec = f"""
    # Setting conda environment
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate python_phylo

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    mkdir -p {gene_dir}
    mkdir -p {log_path}

    python3 {script_path}/parse_gff3.py \\
        --gff3 {gff3_file} \\
        --gene_id {gene_id} \\
        --output_dir {gene_dir}

    echo "FINISH: $(date)"
    echo done > {outputs["log"]}
    """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def extract_genotypes_template(output_path: str, script_path: str, log_path: str,
                               gvcf_file: str, gene_bed: str, exon_bed: str,
                               intron_bed: str, individual: str, gene_name: str,
                               species: str, family: str, min_dp: int, min_gq: int):
    """
    Template for extracting genotypes from gVCF for a gene region.

    Args:
        output_path: Base output directory
        script_path: Path to helper scripts
        log_path: Path for log files
        gvcf_file: Path to individual's gVCF file
        gene_bed: Gene coordinates BED file
        exon_bed: Exon coordinates BED file
        intron_bed: Intron coordinates BED file
        individual: Individual identifier
        gene_name: Gene name
        species: Species identifier
        family: Family identifier
        min_dp: Minimum depth threshold
        min_gq: Minimum genotype quality threshold
    """
    out_dir = f"{output_path}/{species}/{family}/{gene_name}"

    inputs = {
        "gvcf": gvcf_file,
        "gene_bed": gene_bed,
        "exon_bed": exon_bed,
        "intron_bed": intron_bed
    }
    outputs = {
        "genotypes": f"{out_dir}/{individual}_genotypes.tsv",
        "log": f"{log_path}/{species}_{family}_{gene_name}_{individual}_genotypes.DONE"
    }
    options = {
        'cores': 1,
        'memory': '8g',
        'walltime': '08:00:00',
        'account': 'spider2'
    }

    spec = f"""
    # Setting conda environment
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate python_phylo

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    mkdir -p {out_dir}
    mkdir -p {log_path}

    python3 {script_path}/extract_genotypes.py \\
        --gvcf {gvcf_file} \\
        --gene_bed {gene_bed} \\
        --exon_bed {exon_bed} \\
        --intron_bed {intron_bed} \\
        --output {outputs["genotypes"]} \\
        --min_dp {min_dp} \\
        --min_gq {min_gq}

    echo "FINISH: $(date)"
    echo done > {outputs["log"]}
    """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def calculate_heterozygosity_template(output_path: str, script_path: str,
                                      log_path: str, genotypes_file: str,
                                      individual: str, gene_name: str,
                                      species: str, family: str):
    """
    Template for calculating heterozygosity metrics.

    Args:
        output_path: Base output directory
        script_path: Path to helper scripts
        log_path: Path for log files
        genotypes_file: Input genotypes TSV file
        individual: Individual identifier
        gene_name: Gene name
        species: Species identifier
        family: Family identifier
    """
    out_dir = f"{output_path}/{species}/{family}/{gene_name}"

    inputs = {"genotypes": genotypes_file}
    outputs = {
        "heterozygosity": f"{out_dir}/{individual}_heterozygosity.tsv",
        "het_sites": f"{out_dir}/{individual}_het_sites.tsv",
        "log": f"{log_path}/{species}_{family}_{gene_name}_{individual}_het.DONE"
    }
    options = {
        'cores': 1,
        'memory': '4g',
        'walltime': '01:00:00',
        'account': 'spider2'
    }

    spec = f"""
    # Setting conda environment
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate python_phylo

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    mkdir -p {log_path}

    python3 {script_path}/calculate_heterozygosity.py \\
        --genotypes {genotypes_file} \\
        --individual {individual} \\
        --gene {gene_name} \\
        --output_het {outputs["heterozygosity"]} \\
        --output_sites {outputs["het_sites"]}

    echo "FINISH: $(date)"
    echo done > {outputs["log"]}
    """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def calculate_dxy_template(output_path: str, script_path: str, log_path: str,
                           female_genotypes: str, male_genotypes: str,
                           gene_name: str, species: str, family: str):
    """
    Template for calculating Dxy between female and male.

    Args:
        output_path: Base output directory
        script_path: Path to helper scripts
        log_path: Path for log files
        female_genotypes: Female genotypes TSV file
        male_genotypes: Male genotypes TSV file
        gene_name: Gene name
        species: Species identifier
        family: Family identifier
    """
    out_dir = f"{output_path}/{species}/{family}/{gene_name}"

    inputs = {
        "female_genotypes": female_genotypes,
        "male_genotypes": male_genotypes
    }
    outputs = {
        "dxy": f"{out_dir}/dxy_{family}.tsv",
        "log": f"{log_path}/{species}_{family}_{gene_name}_dxy.DONE"
    }
    options = {
        'cores': 1,
        'memory': '4g',
        'walltime': '01:00:00',
        'account': 'spider2'
    }

    spec = f"""
    # Setting conda environment
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate python_phylo

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    mkdir -p {log_path}

    python3 {script_path}/calculate_dxy.py \\
        --female {female_genotypes} \\
        --male {male_genotypes} \\
        --gene {gene_name} \\
        --family {family} \\
        --output {outputs["dxy"]}

    echo "FINISH: $(date)"
    echo done > {outputs["log"]}
    """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def generate_summary_template(output_path: str, script_path: str, log_path: str,
                              species: str, het_files: list, dxy_files: list,
                              het_site_files: list):
    """
    Template for generating summary reports.

    Args:
        output_path: Base output directory
        script_path: Path to helper scripts
        log_path: Path for log files
        species: Species identifier
        het_files: List of heterozygosity TSV files
        dxy_files: List of Dxy TSV files
        het_site_files: List of het sites TSV files
    """
    out_dir = f"{output_path}/{species}"

    # Convert lists to space-separated strings for command line
    het_files_str = " ".join(het_files)
    dxy_files_str = " ".join(dxy_files)
    sites_files_str = " ".join(het_site_files)

    # Build inputs dict
    inputs = {}
    for i, f in enumerate(het_files):
        inputs[f"het_{i}"] = f
    for i, f in enumerate(dxy_files):
        inputs[f"dxy_{i}"] = f
    for i, f in enumerate(het_site_files):
        inputs[f"sites_{i}"] = f

    outputs = {
        "het_summary": f"{out_dir}/{species}_heterozygosity_summary.tsv",
        "dxy_summary": f"{out_dir}/{species}_dxy_summary.tsv",
        "sites_catalog": f"{out_dir}/{species}_het_sites_catalog.tsv",
        "total_het": f"{out_dir}/{species}_individual_total_heterozygosity.tsv",
        "total_dxy": f"{out_dir}/{species}_family_total_dxy.tsv",
        "log": f"{log_path}/{species}_summary.DONE"
    }
    options = {
        'cores': 1,
        'memory': '8g',
        'walltime': '01:00:00',
        'account': 'spider2'
    }

    spec = f"""
    # Setting conda environment
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate python_phylo

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    mkdir -p {out_dir}
    mkdir -p {log_path}

    python3 {script_path}/generate_summary.py \\
        --species {species} \\
        --het_files {het_files_str} \\
        --dxy_files {dxy_files_str} \\
        --sites_files {sites_files_str} \\
        --output_het {outputs["het_summary"]} \\
        --output_dxy {outputs["dxy_summary"]} \\
        --output_sites {outputs["sites_catalog"]} \\
        --output_total_het {outputs["total_het"]} \\
        --output_total_dxy {outputs["total_dxy"]}

    echo "FINISH: $(date)"
    echo done > {outputs["log"]}
    """

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
