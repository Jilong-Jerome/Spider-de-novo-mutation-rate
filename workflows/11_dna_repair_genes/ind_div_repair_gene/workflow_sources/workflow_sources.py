#!/usr/bin/env python3
"""
workflow_sources.py - Workflow orchestration for repair gene genotype analysis

This module reads configuration and orchestrates the execution of:
1. parse_gff3 - Extract gene/exon/intron coordinates (per gene)
2. extract_genotypes - Extract genotypes from gVCF (per individual per gene)
3. calculate_heterozygosity - Calculate het metrics (per individual per gene)
4. calculate_dxy - Calculate pairwise Dxy (per family per gene)
5. generate_summary - Compile summary reports (per species)
"""

from gwf import Workflow
import os
import yaml
import glob
import re
from workflow_templates import (
    parse_gff3_template,
    extract_genotypes_template,
    calculate_heterozygosity_template,
    calculate_dxy_template,
    generate_summary_template
)


def sanitize_name(name):
    """
    Sanitize a name to be valid for gwf target names.
    GWF only allows alphanumeric characters and underscores.

    Args:
        name: Original name string

    Returns:
        Sanitized name string
    """
    # Replace hyphens and other special characters with underscores
    sanitized = re.sub(r'[^a-zA-Z0-9_]', '_', name)
    # Remove leading numbers (not allowed as first character)
    sanitized = re.sub(r'^[0-9]+', '', sanitized)
    # Remove consecutive underscores
    sanitized = re.sub(r'_+', '_', sanitized)
    # Remove leading/trailing underscores
    sanitized = sanitized.strip('_')
    return sanitized


def parse_repair_genes(repair_gene_file):
    """
    Parse repair gene TSV file.

    Args:
        repair_gene_file: Path to TSV file with gene_name and gene_id columns

    Returns:
        dict: {gene_name: {'gene_id': gene_id, 'safe_name': sanitized_name}}
    """
    repair_genes = {}

    with open(repair_gene_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) >= 2:
                gene_name = parts[0]
                gene_id = parts[1]
                safe_name = sanitize_name(gene_name)
                repair_genes[gene_name] = {
                    'gene_id': gene_id,
                    'safe_name': safe_name
                }

    return repair_genes


def discover_families(ind_dir, species):
    """
    Discover families and individuals from directory structure.

    Args:
        ind_dir: Path to individuals directory
        species: Species identifier

    Returns:
        dict: {family_key: {'F': female_name, 'M': male_name}}
    """
    families = {}

    if not os.path.exists(ind_dir):
        return families

    for entry in os.listdir(ind_dir):
        # Pattern: {species}_family{N}_{sex_code}_{sex_label}
        # Example: AFR_family1_F_female
        match = re.match(rf'{species}_family(\d+)_(F|M)_(female|male)', entry)

        if match:
            family_num = match.group(1)
            sex_code = match.group(2)
            family_key = f"family{family_num}"

            if family_key not in families:
                families[family_key] = {}

            families[family_key][sex_code] = entry

    return families


def repair_gene_analysis_workflow(config_file: str, gwf):
    """
    Main workflow function for repair gene genotype analysis.

    Pipeline DAG:
    parse_gff3 (per gene)
        |
        v
    extract_genotypes (per individual per gene)
        |
        +---> calculate_heterozygosity (per individual per gene)
        |
        +---> calculate_dxy (per family per gene)
                  |
                  v
              generate_summary (per species)

    Args:
        config_file: Path to configuration YAML file
        gwf: GWF Workflow object

    Returns:
        Updated gwf Workflow object
    """
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    CONFIG = yaml.safe_load(open(config_file))

    ACCOUNT = CONFIG['account']
    PROJECT_ID = CONFIG['project_id']
    BASE_DIR = CONFIG['base_directory']
    OUTPUT_DIR = CONFIG['output_directory_path']
    LOG_DIR = CONFIG['log_directory_path']
    SCRIPTS_PATH = CONFIG['scripts_path']
    DATA_DIR = CONFIG['data_directory_path']
    SPECIES_LIST = CONFIG['species']
    MIN_DP = CONFIG.get('min_depth', 10)
    MIN_GQ = CONFIG.get('min_genotype_quality', 20)

    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------

    for species in SPECIES_LIST:
        # Species-specific paths
        species_data_dir = f"{DATA_DIR}/{species}"
        gff3_file = f"{species_data_dir}/{species}_gene.gff3"
        repair_gene_file = f"{species_data_dir}/{species}_repair_gene.tsv"
        ind_dir = f"{species_data_dir}/ind"

        # Validate required files exist
        if not os.path.exists(gff3_file):
            print(f"Warning: GFF3 file not found for {species}: {gff3_file}")
            continue

        if not os.path.exists(repair_gene_file):
            print(f"Warning: Repair gene file not found for {species}: {repair_gene_file}")
            continue

        # Parse repair genes
        repair_genes = parse_repair_genes(repair_gene_file)
        if not repair_genes:
            print(f"Warning: No repair genes found for {species}")
            continue

        # Discover families
        families = discover_families(ind_dir, species)
        if not families:
            print(f"Warning: No families found for {species}")
            continue

        print(f"Processing species: {species}")
        print(f"  Repair genes: {len(repair_genes)}")
        print(f"  Families: {len(families)}")

        # Track output files for summary
        all_het_files = []
        all_dxy_files = []
        all_het_site_files = []

        # --------------------------------------------------
        # Step 1: Parse GFF3 for each gene (once per species)
        # --------------------------------------------------
        gene_outputs = {}

        for gene_name, gene_info in repair_genes.items():
            gene_id = gene_info['gene_id']
            safe_name = gene_info['safe_name']

            task_name = f"{species}_{safe_name}_parse_gff3"
            gene_dir = f"{OUTPUT_DIR}/{species}/{gene_name}"

            gwf.target_from_template(
                name=task_name,
                template=parse_gff3_template(
                    output_path=OUTPUT_DIR,
                    script_path=SCRIPTS_PATH,
                    log_path=LOG_DIR,
                    gff3_file=gff3_file,
                    gene_id=gene_id,
                    gene_name=gene_name,
                    species=species
                )
            )

            gene_outputs[gene_name] = {
                "gene_bed": f"{gene_dir}/gene_coords.bed",
                "exon_bed": f"{gene_dir}/exon_coords.bed",
                "intron_bed": f"{gene_dir}/intron_coords.bed",
                "safe_name": safe_name
            }

        # --------------------------------------------------
        # Steps 2-4: For each family, for each gene
        # --------------------------------------------------
        for family, members in families.items():
            # Skip incomplete families
            if 'F' not in members or 'M' not in members:
                print(f"  Warning: Incomplete family {family}, skipping")
                continue

            female = members['F']
            male = members['M']

            for gene_name in repair_genes.keys():
                gene_coords = gene_outputs[gene_name]
                safe_name = gene_coords['safe_name']
                family_gene_dir = f"{OUTPUT_DIR}/{species}/{family}/{gene_name}"

                # Sanitize individual names for task naming
                def safe_ind_name(ind):
                    return sanitize_name(ind)

                # Track genotype files for Dxy calculation
                genotype_files = {}

                for ind, sex_label in [(female, 'F'), (male, 'M')]:
                    # Individual gVCF path
                    ind_gvcf = f"{ind_dir}/{ind}/{ind}.g.vcf"
                    ind_safe = safe_ind_name(ind)

                    # --------------------------------------------------
                    # Step 2: Extract genotypes
                    # --------------------------------------------------
                    extract_task = f"{species}_{family}_{safe_name}_{ind_safe}_genotypes"

                    gwf.target_from_template(
                        name=extract_task,
                        template=extract_genotypes_template(
                            output_path=OUTPUT_DIR,
                            script_path=SCRIPTS_PATH,
                            log_path=LOG_DIR,
                            gvcf_file=ind_gvcf,
                            gene_bed=gene_coords["gene_bed"],
                            exon_bed=gene_coords["exon_bed"],
                            intron_bed=gene_coords["intron_bed"],
                            individual=ind,
                            gene_name=gene_name,
                            species=species,
                            family=family,
                            min_dp=MIN_DP,
                            min_gq=MIN_GQ
                        )
                    )

                    genotypes_file = f"{family_gene_dir}/{ind}_genotypes.tsv"
                    genotype_files[sex_label] = genotypes_file

                    # --------------------------------------------------
                    # Step 3: Calculate heterozygosity
                    # --------------------------------------------------
                    het_task = f"{species}_{family}_{safe_name}_{ind_safe}_het"

                    gwf.target_from_template(
                        name=het_task,
                        template=calculate_heterozygosity_template(
                            output_path=OUTPUT_DIR,
                            script_path=SCRIPTS_PATH,
                            log_path=LOG_DIR,
                            genotypes_file=genotypes_file,
                            individual=ind,
                            gene_name=gene_name,
                            species=species,
                            family=family
                        )
                    )

                    # Track output files for summary
                    all_het_files.append(
                        f"{family_gene_dir}/{ind}_heterozygosity.tsv"
                    )
                    all_het_site_files.append(
                        f"{family_gene_dir}/{ind}_het_sites.tsv"
                    )

                # --------------------------------------------------
                # Step 4: Calculate Dxy between pair
                # --------------------------------------------------
                dxy_task = f"{species}_{family}_{safe_name}_dxy"

                gwf.target_from_template(
                    name=dxy_task,
                    template=calculate_dxy_template(
                        output_path=OUTPUT_DIR,
                        script_path=SCRIPTS_PATH,
                        log_path=LOG_DIR,
                        female_genotypes=genotype_files['F'],
                        male_genotypes=genotype_files['M'],
                        gene_name=gene_name,
                        species=species,
                        family=family
                    )
                )

                all_dxy_files.append(
                    f"{family_gene_dir}/dxy_{family}.tsv"
                )

        # --------------------------------------------------
        # Step 5: Generate summary for species
        # --------------------------------------------------
        summary_task = f"{species}_summary"

        gwf.target_from_template(
            name=summary_task,
            template=generate_summary_template(
                output_path=OUTPUT_DIR,
                script_path=SCRIPTS_PATH,
                log_path=LOG_DIR,
                species=species,
                het_files=all_het_files,
                dxy_files=all_dxy_files,
                het_site_files=all_het_site_files
            )
        )

    return gwf
