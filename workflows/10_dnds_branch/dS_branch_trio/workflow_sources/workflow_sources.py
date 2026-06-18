#!/usr/bin/env python3
"""
workflow_sources.py - Workflow orchestration for dS branch trio analysis

Pipeline DAG:
    Part 1:
    01_bwa2_index
        ├──> 02_align_sp1 ──> 03a..03f_sp1 ──> 04_call_sp1 ──> 05a..05f_sp1 ──┐
        │                                                                     ├──> 06_extract_cds
        └──> 02_align_sp2 ──> 03a..03f_sp2 ──> 04_call_sp2 ──> 05a..05f_sp2 ──┘

    Part 2:
    06_extract_cds ──> 07_prepare_gene_lists
        └──> 08a_sample_auto_all  ──> 08b_paml_auto_all
             08a_sample_auto_bs_1 ──> 08b_paml_auto_bs_1
             ...
             08a_sample_auto_bs_500 ──> 08b_paml_auto_bs_500

Total: 20 (Part 1) + 1 + 501*2 (Part 2) = 1023 SLURM jobs
"""

import os
import yaml
from workflow_templates import (
    bwa2_index_template,
    bwa2_align_template,
    picard_addRG_template,
    samtools_fixmate_template,
    samtools_sort_template,
    samtools_markdup_template,
    samtools_filter_template,
    samtools_index_template,
    bcftools_call_template,
    callable_depth_template,
    coverage_distribution_plot_template,
    consensus_filter_snps_template,
    consensus_classify_snps_chrom_template,
    consensus_classify_snps_merge_template,
    consensus_select_snps_template,
    consensus_mask_template,
    bcftools_consensus_template,
    extract_cds_template,
    prepare_gene_lists_template,
    bootstrap_sample_template,
    bootstrap_paml_template,
    summarize_bootstrap_template,
    visualize_pairwise_dS_template,
    combine_trio_visualizations_template,
)


def trio_dnds_workflow(config_file, gwf):
    """
    Main workflow function for dS branch trio analysis.

    Args:
        config_file: Path to configuration YAML file
        gwf: GWF Workflow object

    Returns:
        Updated gwf Workflow object
    """
    # --------------------------------------------------
    #  Configuration
    # --------------------------------------------------
    CONFIG = yaml.safe_load(open(config_file))

    BASE_DIR = CONFIG['base_directory']
    SCRIPTS_PATH = CONFIG['scripts_path']

    # Trio identifier for namespacing output dirs and job names
    TRIO = CONFIG['trio_name']
    OUTPUT_DIR = os.path.join(CONFIG['output_directory_path'], TRIO)
    LOG_DIR = os.path.join(CONFIG['log_directory_path'], TRIO)

    REF_GENOME = os.path.join(BASE_DIR, CONFIG['reference_genome'])
    REF_ANNOTATION = os.path.join(BASE_DIR, CONFIG['reference_annotation'])
    REF_SPECIES = CONFIG.get('reference_species', 'BIC')

    MAPPING_SPECIES = CONFIG['mapping_species']
    species_list = list(MAPPING_SPECIES.keys())
    sp1_name, sp2_name = species_list[0], species_list[1]
    TRIO_TREE = CONFIG['trio_tree']

    MIN_MQ = CONFIG['min_mapping_quality']
    MIN_BQ = CONFIG['min_base_quality']
    MIN_DP = CONFIG['min_depth']
    MIN_COV_FRAC = CONFIG.get('min_coverage_fraction', 0.80)
    COVERAGE_MEDIAN_MIN_FACTOR = CONFIG.get('coverage_median_min_factor', 0.5)
    COVERAGE_MEDIAN_MAX_FACTOR = CONFIG.get('coverage_median_max_factor', 2.0)

    # Step output directories (all under steps/{TRIO}/)
    INDEX_DIR = os.path.join(OUTPUT_DIR, '01_index')
    MAP_DIR = os.path.join(OUTPUT_DIR, '02_mapping')
    BAM_DIR = os.path.join(OUTPUT_DIR, '03_bam_processing')
    VCF_DIR = os.path.join(OUTPUT_DIR, '04_variant_calling')
    CALLABLE_DEPTH_DIR = os.path.join(OUTPUT_DIR, '05_callable_depth')
    CONS_DIR = os.path.join(OUTPUT_DIR, '05_consensus')
    CDS_DIR = os.path.join(OUTPUT_DIR, '06_cds_per_gene')

    # Contig list used to fan out step 05c (classify SNPs) by chromosome.
    # Autosomes come straight from config; X-linked contigs are resolved by
    # prefix-matching the reference .fai so configs only need to name the
    # prefix (e.g., bic_X → bic_X1, bic_X2).
    autosomal_chroms = list(CONFIG.get(
        'autosomal_chromosomes',
        [f'bic_{i}' for i in range(1, 15)]
    ))
    x_prefixes = list(CONFIG.get('x_chromosome_prefixes', ['bic_X']))
    fai_path = REF_GENOME + '.fai'
    x_chroms = []
    if os.path.exists(fai_path):
        with open(fai_path) as fai:
            for line in fai:
                name = line.split('\t', 1)[0]
                if any(name.startswith(p) for p in x_prefixes):
                    x_chroms.append(name)
    classify_chroms = autosomal_chroms + sorted(x_chroms)

    # --------------------------------------------------
    #  Step 01: Reference Indexing (1 job)
    # --------------------------------------------------
    gwf.target_from_template(
        name=f'{TRIO}_bwa2_index_{REF_SPECIES}',
        template=bwa2_index_template(
            ref_genome=REF_GENOME,
            ref_species=REF_SPECIES,
            output_dir=INDEX_DIR,
            log_dir=LOG_DIR
        )
    )

    # --------------------------------------------------
    #  Per-species steps (SAR, PAC run in parallel)
    # --------------------------------------------------
    for species, reads in MAPPING_SPECIES.items():
        read1 = os.path.join(BASE_DIR, reads['read1'])
        read2 = os.path.join(BASE_DIR, reads['read2'])
        ref_index_prefix = os.path.join(INDEX_DIR, REF_SPECIES)

        # Step 02: Read Mapping (1 job per species)
        gwf.target_from_template(
            name=f'{TRIO}_bwa2_align_{species}',
            template=bwa2_align_template(
                ref_index_prefix=ref_index_prefix,
                read1=read1,
                read2=read2,
                species=species,
                output_dir=BAM_DIR,
                log_dir=LOG_DIR
            )
        )

        # Step 03: BAM Processing chain (6 jobs per species)
        gwf.target_from_template(
            name=f'{TRIO}_picard_addRG_{species}',
            template=picard_addRG_template(
                species=species,
                bam_dir=BAM_DIR,
                log_dir=LOG_DIR
            )
        )

        gwf.target_from_template(
            name=f'{TRIO}_samtools_fixmate_{species}',
            template=samtools_fixmate_template(
                species=species,
                bam_dir=BAM_DIR,
                log_dir=LOG_DIR
            )
        )

        gwf.target_from_template(
            name=f'{TRIO}_samtools_sort_{species}',
            template=samtools_sort_template(
                species=species,
                bam_dir=BAM_DIR,
                log_dir=LOG_DIR
            )
        )

        gwf.target_from_template(
            name=f'{TRIO}_samtools_markdup_{species}',
            template=samtools_markdup_template(
                species=species,
                bam_dir=BAM_DIR,
                log_dir=LOG_DIR
            )
        )

        gwf.target_from_template(
            name=f'{TRIO}_samtools_filter_{species}',
            template=samtools_filter_template(
                species=species,
                bam_dir=BAM_DIR,
                log_dir=LOG_DIR,
                min_mq=MIN_MQ
            )
        )

        gwf.target_from_template(
            name=f'{TRIO}_samtools_index_{species}',
            template=samtools_index_template(
                species=species,
                bam_dir=BAM_DIR,
                log_dir=LOG_DIR
            )
        )

        # Step 04: Variant Calling (1 job per species)
        gwf.target_from_template(
            name=f'{TRIO}_bcftools_call_{species}',
            template=bcftools_call_template(
                ref_genome=REF_GENOME,
                species=species,
                bam_dir=BAM_DIR,
                output_dir=VCF_DIR,
                log_dir=LOG_DIR,
                min_bq=MIN_BQ
            )
        )

        # Step 05a: Shared callable-depth mask (1 job per species)
        gwf.target_from_template(
            name=f'{TRIO}_callable_depth_{species}',
            template=callable_depth_template(
                species=species,
                bam_dir=BAM_DIR,
                output_dir=CALLABLE_DEPTH_DIR,
                log_dir=LOG_DIR,
                script_path=SCRIPTS_PATH,
                min_dp=MIN_DP,
                coverage_median_min_factor=COVERAGE_MEDIAN_MIN_FACTOR,
                coverage_median_max_factor=COVERAGE_MEDIAN_MAX_FACTOR
            )
        )

        # Step 05b: Shared SNP filter (1 job per species)
        gwf.target_from_template(
            name=f'{TRIO}_filter_consensus_snps_{species}',
            template=consensus_filter_snps_template(
                species=species,
                vcf_dir=VCF_DIR,
                output_dir=CALLABLE_DEPTH_DIR,
                callable_depth_dir=CALLABLE_DEPTH_DIR,
                shared_log_dir=LOG_DIR
            )
        )

        # Step 05c: Shared clean-support classification.
        # Fan out per chromosome (one job per contig) and merge into the
        # whole-genome classification TSV that downstream branches consume.
        for chrom in classify_chroms:
            gwf.target_from_template(
                name=f'{TRIO}_classify_consensus_snps_{species}_{chrom}',
                template=consensus_classify_snps_chrom_template(
                    species=species,
                    chrom=chrom,
                    bam_dir=BAM_DIR,
                    output_dir=CALLABLE_DEPTH_DIR,
                    shared_log_dir=LOG_DIR,
                    script_path=SCRIPTS_PATH,
                    min_bq=MIN_BQ
                )
            )
        gwf.target_from_template(
            name=f'{TRIO}_classify_consensus_snps_{species}',
            template=consensus_classify_snps_merge_template(
                species=species,
                chroms=classify_chroms,
                output_dir=CALLABLE_DEPTH_DIR,
                shared_log_dir=LOG_DIR,
                script_path=SCRIPTS_PATH
            )
        )

    # Step 05a (diagnostic): per-trio coverage-distribution TSV + plots.
    # Reference-only leaf step; not consumed by any downstream target.
    gwf.target_from_template(
        name=f'{TRIO}_coverage_distribution',
        template=coverage_distribution_plot_template(
            species_list=species_list,
            callable_depth_dir=CALLABLE_DEPTH_DIR,
            log_dir=LOG_DIR,
            script_path=SCRIPTS_PATH,
            trio_name=TRIO
        )
    )

    # ==========================================================
    #  Branch-specific consensus and downstream analysis
    # ==========================================================

    valid_modes = ['majority', 'random', 'using_ref', 'using_alt']
    default_modes = ['random']
    consensus_modes = CONFIG.get('consensus_modes', default_modes)
    unknown_modes = [m for m in consensus_modes if m not in valid_modes]
    if unknown_modes:
        raise ValueError(f"{config_file}: unknown consensus modes: {unknown_modes}")
    valid_support_policies = ['strict', 'relaxed']
    support_policies = CONFIG.get(
        'consensus_support_policies', valid_support_policies)
    unknown_policies = [
        p for p in support_policies if p not in valid_support_policies]
    if unknown_policies:
        raise ValueError(
            f"{config_file}: unknown consensus support policies: "
            f"{unknown_policies}")
    consensus_branches = [
        (policy, mode, f'{policy}_{mode}')
        for policy in support_policies
        for mode in consensus_modes
    ]
    random_seed = CONFIG.get('random_consensus_seed', 20260514)

    # Bootstrap config
    N_BOOTSTRAP = CONFIG.get('n_bootstrap', 500)
    AUTOSOMAL_CHROMS = ','.join(CONFIG.get(
        'autosomal_chromosomes',
        [f'bic_{i}' for i in range(1, 15)]
    ))
    X_CHROM_PREFIXES = ','.join(CONFIG.get('x_chromosome_prefixes', ['bic_X']))

    # Build list of replicate IDs: auto_all + auto_bs_1..auto_bs_500
    replicate_ids = ['auto_all'] + [f'auto_bs_{i}' for i in range(1, N_BOOTSTRAP + 1)]

    for support_policy, mode, branch_label in consensus_branches:
        MODE_OUTPUT_DIR = os.path.join(OUTPUT_DIR, branch_label)
        MODE_LOG_DIR = os.path.join(LOG_DIR, branch_label)
        CONS_DIR = os.path.join(MODE_OUTPUT_DIR, '05_consensus')
        CDS_DIR = os.path.join(MODE_OUTPUT_DIR, '06_cds_per_gene')
        GENE_LIST_DIR = os.path.join(MODE_OUTPUT_DIR, '07_gene_lists')
        BOOTSTRAP_DIR = os.path.join(MODE_OUTPUT_DIR, '08_bootstrap')
        SUMMARY_DIR = os.path.join(MODE_OUTPUT_DIR, '09_summary')
        VIS_DIR = os.path.join(MODE_OUTPUT_DIR, '10_visualization')
        FASTA_DIR = os.path.join(CDS_DIR, 'per_gene_fasta')
        PASSING_GENES = os.path.join(CDS_DIR, 'passing_genes.txt')

        # --------------------------------------------------
        #  Step 05: Branch-specific consensus chain
        # --------------------------------------------------
        for species in MAPPING_SPECIES:
            gwf.target_from_template(
                name=f'{TRIO}_{branch_label}_select_consensus_snps_{species}',
                template=consensus_select_snps_template(
                    species=species,
                    output_dir=CONS_DIR,
                    shared_consensus_dir=CALLABLE_DEPTH_DIR,
                    shared_log_dir=LOG_DIR,
                    log_dir=MODE_LOG_DIR,
                    script_path=SCRIPTS_PATH,
                    consensus_mode=mode,
                    support_policy=support_policy,
                    branch_label=branch_label,
                    random_seed=random_seed
                )
            )

            gwf.target_from_template(
                name=f'{TRIO}_{branch_label}_build_consensus_mask_{species}',
                template=consensus_mask_template(
                    species=species,
                    output_dir=CONS_DIR,
                    callable_depth_dir=CALLABLE_DEPTH_DIR,
                    shared_log_dir=LOG_DIR,
                    log_dir=MODE_LOG_DIR,
                    script_path=SCRIPTS_PATH,
                    branch_label=branch_label
                )
            )

            gwf.target_from_template(
                name=f'{TRIO}_{branch_label}_bcftools_consensus_{species}',
                template=bcftools_consensus_template(
                    ref_genome=REF_GENOME,
                    species=species,
                    output_dir=CONS_DIR,
                    log_dir=MODE_LOG_DIR,
                    branch_label=branch_label
                )
            )

        # --------------------------------------------------
        #  Step 06: CDS Extraction (1 job per branch)
        # --------------------------------------------------
        sp1_consensus = os.path.join(CONS_DIR, f'{sp1_name}_consensus.fa')
        sp2_consensus = os.path.join(CONS_DIR, f'{sp2_name}_consensus.fa')

        gwf.target_from_template(
            name=f'{TRIO}_{branch_label}_extract_cds_per_gene',
            template=extract_cds_template(
                ref_genome=REF_GENOME,
                ref_name=REF_SPECIES,
                sp1_name=sp1_name,
                sp1_consensus=sp1_consensus,
                sp2_name=sp2_name,
                sp2_consensus=sp2_consensus,
                gff3=REF_ANNOTATION,
                output_dir=CDS_DIR,
                log_dir=MODE_LOG_DIR,
                script_path=SCRIPTS_PATH,
                min_cov_frac=MIN_COV_FRAC
            )
        )

        # --------------------------------------------------
        #  Step 07: Prepare Gene Lists (1 job per branch)
        # --------------------------------------------------
        gwf.target_from_template(
            name=f'{TRIO}_{branch_label}_prepare_gene_lists',
            template=prepare_gene_lists_template(
                gff3=REF_ANNOTATION,
                passing_genes=PASSING_GENES,
                output_dir=GENE_LIST_DIR,
                log_dir=MODE_LOG_DIR,
                script_path=SCRIPTS_PATH,
                autosomal_chromosomes=AUTOSOMAL_CHROMS,
                x_chromosome_prefixes=X_CHROM_PREFIXES
            )
        )

        AUTO_GENE_LIST = os.path.join(GENE_LIST_DIR, 'auto_passing_genes.txt')

        # --------------------------------------------------
        #  Step 08: Bootstrap sampling + PAML for each replicate
        # --------------------------------------------------
        paml_done_logs = []
        for rep_id in replicate_ids:
            rep_work_dir = os.path.join(BOOTSTRAP_DIR, rep_id)
            id_list_file = os.path.join(rep_work_dir, f'{rep_id}_id.txt')

            gwf.target_from_template(
                name=f'{TRIO}_{branch_label}_sample_{rep_id}',
                template=bootstrap_sample_template(
                    gene_list=AUTO_GENE_LIST,
                    replicate_id=rep_id,
                    work_dir=rep_work_dir,
                    log_dir=MODE_LOG_DIR
                )
            )

            paml_done_log = f'{MODE_LOG_DIR}/08b_paml_{rep_id}.DONE'
            paml_done_logs.append(paml_done_log)
            gwf.target_from_template(
                name=f'{TRIO}_{branch_label}_paml_{rep_id}',
                template=bootstrap_paml_template(
                    replicate_id=rep_id,
                    id_list_file=id_list_file,
                    fasta_dir=FASTA_DIR,
                    work_dir=rep_work_dir,
                    log_dir=MODE_LOG_DIR,
                    script_path=SCRIPTS_PATH,
                    trio_tree=TRIO_TREE
                )
            )

        # --------------------------------------------------
        #  Step 09: Summarize Bootstrap Results (1 job per branch)
        # --------------------------------------------------
        gwf.target_from_template(
            name=f'{TRIO}_{branch_label}_summarize_bootstrap',
            template=summarize_bootstrap_template(
                bootstrap_dir=BOOTSTRAP_DIR,
                n_bootstrap=N_BOOTSTRAP,
                output_dir=SUMMARY_DIR,
                log_dir=MODE_LOG_DIR,
                script_path=SCRIPTS_PATH,
                paml_done_logs=paml_done_logs
            )
        )

        # --------------------------------------------------
        #  Step 10: Visualize Pairwise dS (1 job per branch)
        # --------------------------------------------------
        gwf.target_from_template(
            name=f'{TRIO}_{branch_label}_visualize_pairwise_dS',
            template=visualize_pairwise_dS_template(
                summary_dir=SUMMARY_DIR,
                output_dir=VIS_DIR,
                log_dir=MODE_LOG_DIR,
                script_path=SCRIPTS_PATH,
                trio_name=f'{TRIO}_{branch_label}',
                trio_tree=TRIO_TREE,
                summary_done_log=f'{MODE_LOG_DIR}/09_summarize_bootstrap.DONE',
            )
        )

    return gwf


def combine_trio_visualization_workflow(config_files, gwf):
    """Cross-trio Step 11: build a single 2 x N combined figure across trios.

    Depends on every per-trio Step 09 having produced its pairwise summary +
    bootstrap TSVs (the per-trio Step 10 sentinel is not strictly required, so
    we wire the dependency to the Step 09 sentinels — that lets this run as
    soon as all summaries exist).
    """
    if not config_files:
        return gwf

    output_root = None
    log_root = None
    scripts_path = None
    valid_modes = ['majority', 'random', 'using_ref', 'using_alt']
    default_modes = ['random']
    valid_support_policies = ['strict', 'relaxed']
    branches = []
    config_data = []
    for cf in config_files:
        with open(cf) as f:
            cfg = yaml.safe_load(f)
        config_data.append(cfg)
        consensus_modes = cfg.get('consensus_modes', default_modes)
        support_policies = cfg.get(
            'consensus_support_policies', valid_support_policies)
        for policy in support_policies:
            if policy not in valid_support_policies:
                raise ValueError(
                    f"{cf}: unknown consensus support policy: {policy}")
            for mode in consensus_modes:
                if mode not in valid_modes:
                    raise ValueError(f"{cf}: unknown consensus mode: {mode}")
                branch_label = f'{policy}_{mode}'
                if branch_label not in branches:
                    branches.append(branch_label)
        output_root = output_root or cfg['output_directory_path']
        log_root = log_root or cfg['log_directory_path']
        scripts_path = scripts_path or cfg['scripts_path']

    for branch_label in branches:
        trios = []
        for cfg in config_data:
            consensus_modes = cfg.get('consensus_modes', default_modes)
            support_policies = cfg.get(
                'consensus_support_policies', valid_support_policies)
            cfg_branches = {
                f'{policy}_{mode}'
                for policy in support_policies
                for mode in consensus_modes
            }
            if branch_label not in cfg_branches:
                continue
            trio_name = cfg['trio_name']
            trio_tree = cfg['trio_tree']
            summary_dir = os.path.join(cfg['output_directory_path'],
                                       trio_name, branch_label, '09_summary')
            log_dir = os.path.join(
                cfg['log_directory_path'], trio_name, branch_label)
            trios.append({
                'trio_name': f'{trio_name}_{branch_label}',
                'trio_tree': trio_tree,
                'summary_dir': summary_dir,
                'summary_done_log': f'{log_dir}/09_summarize_bootstrap.DONE',
            })

        combined_out = os.path.join(output_root, '_combined', branch_label)
        combined_log = os.path.join(log_root, '_combined', branch_label)

        gwf.target_from_template(
            name=f'combine_trio_visualizations_{branch_label}',
            template=combine_trio_visualizations_template(
                trios=trios,
                output_dir=combined_out,
                log_dir=combined_log,
                script_path=scripts_path,
            )
        )

    return gwf
