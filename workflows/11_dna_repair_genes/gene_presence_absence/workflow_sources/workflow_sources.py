from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *
sp_dict = {
        "AFR":2,
        "BIC":3,
        "DUM":4,
        "LIN":5,
        "MIM":6,
        "SAR":7,
        "TEN":8
        }
def gene_variant_check_workflow(config_file: str, gwf):
    """
    """
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    CONFIG = yaml.safe_load(open(config_file))
    ACCOUNT: str = CONFIG['account']
    SPECIES_ID: str = CONFIG['species_id']
    GENES : list = CONFIG['gene']
    RANDOM_GENES: list = CONFIG['random_gene']
    VCF: str = CONFIG['species_vcf']
    OLD_REF: str = CONFIG['old_ref']
    NEW_REF: str = CONFIG['new_ref']
    WORK_DIR: str = CONFIG['working_directory_path']
    LOG_DIR: str = CONFIG['log_directory_path']
    HELP_SCRIPTS_PATH: str = CONFIG['help_scripts_path']
    GENE_BLAST: str = CONFIG['homologs_blast']
    GENE_BLAST_SP: str = CONFIG['homologs_blast_species']
    GENE_BLAST_SP_GFF: str = CONFIG['homologs_blast_species_gff']
    TARGET_SP_GFF: str = CONFIG['target_species_gff']
    ORTHOGROUPS: str = CONFIG['orthogroups']
    DEPTH_PATH: str = CONFIG['depth_path']
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    gwf.target_from_template(
        name = f"{SPECIES_ID}_depth_summary_per_chrom_per_ind",
        template =  prepare_depth_stat(
            work_path = WORK_DIR,
            log_path = LOG_DIR,
            script_path = HELP_SCRIPTS_PATH,
            vcf = VCF,
            sp = SPECIES_ID)
        )
    for gene in GENES:
        gwf.target_from_template(
            name = f"{SPECIES_ID}_{gene}_prepare",
            template = prepare_gene(
                work_path = WORK_DIR,
                gene = gene, 
                target_sp = SPECIES_ID, 
                target_sp_id = sp_dict[SPECIES_ID],
                gene_blast = GENE_BLAST,
                gene_blast_sp = GENE_BLAST_SP,
                blast_sp_gff = GENE_BLAST_SP_GFF, 
                target_sp_gff = TARGET_SP_GFF,
                orthogroups = ORTHOGROUPS,
                log_path = LOG_DIR)
            )
        gwf.target_from_template(
            name = f"{SPECIES_ID}_{gene}_liftoff",
            template = liftoff_gff(
                work_path = WORK_DIR,
                old_fasta = OLD_REF,
                new_fasta = NEW_REF,
                sp = SPECIES_ID,
                gene = gene,
                log_path = LOG_DIR) 
        )
        gwf.target_from_template(
            name = f"{SPECIES_ID}_{gene}_vcf_retrival",
            template = retrieve_vcf_gene(
                work_path = WORK_DIR,
                vcf = VCF,
                sp = SPECIES_ID,
                gene = gene, 
                log_path = LOG_DIR)
            )
        gwf.target_from_template(
            name = f"{SPECIES_ID}_{gene}_feature_summary",
            template = check_gene_features(
                work_path = WORK_DIR,
                script_path = HELP_SCRIPTS_PATH,
                mean_depth_path = DEPTH_PATH,
                sp = SPECIES_ID,
                gene = gene,
                log_path = LOG_DIR)
            )
    for gene in RANDOM_GENES:
        gwf.target_from_template(
            name = f"{SPECIES_ID}_{gene}_prepare_random",
            template = prepare_gene(
                work_path = f"{WORK_DIR}/random",
                gene = gene, 
                target_sp = SPECIES_ID, 
                target_sp_id = sp_dict[SPECIES_ID],
                gene_blast = GENE_BLAST,
                gene_blast_sp = GENE_BLAST_SP,
                blast_sp_gff = GENE_BLAST_SP_GFF, 
                target_sp_gff = TARGET_SP_GFF,
                orthogroups = ORTHOGROUPS,
                log_path = LOG_DIR)
            )
        gwf.target_from_template(
            name = f"{SPECIES_ID}_{gene}_liftoff_random",
            template = liftoff_gff(
                work_path = f"{WORK_DIR}/random",
                old_fasta = OLD_REF,
                new_fasta = NEW_REF,
                sp = SPECIES_ID,
                gene = gene,
                log_path = LOG_DIR) 
        )
        gwf.target_from_template(
            name = f"{SPECIES_ID}_{gene}_vcf_retrival_random",
            template = retrieve_vcf_gene(
                work_path = f"{WORK_DIR}/random",
                vcf = VCF,
                sp = SPECIES_ID,
                gene = gene, 
                log_path = LOG_DIR)
            )
        gwf.target_from_template(
            name = f"{SPECIES_ID}_{gene}_feature_summary_random",
            template = check_gene_features(
                work_path = f"{WORK_DIR}/random",
                script_path = HELP_SCRIPTS_PATH,
                mean_depth_path = DEPTH_PATH,
                sp = SPECIES_ID,
                gene = gene,
                log_path = LOG_DIR)
            )

    return gwf
