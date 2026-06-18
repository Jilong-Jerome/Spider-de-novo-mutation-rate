from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *
from workflow_dicts import *

def call_dnm_workflow(config_file: str, gwf):

    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    CONFIG = yaml.safe_load(open(config_file))
    ACCOUNT: str = CONFIG['account']
    SPECIES_ID: str = CONFIG['species_id']
    WORK_DIR: str = CONFIG['working_directory_path']
    OUTPUT_DIR: str = CONFIG['output_directory_path']
    LOG_DIR: str = CONFIG['log_directory_path']
    DATA_DIR: str = CONFIG['data_directory_path']
    BAM_PATH: str = CONFIG['bam_directory_path']
    REF_FASTA: str = CONFIG['ref_fasta']
    HELP_SCRIPTS_PATH: str = CONFIG['help_scripts_path']
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------


    # --------------------------------------------------
    #             Callable sites per trio
    # --------------------------------------------------
    for fam in sp_fam_dict[SPECIES_ID]:
        female = f"{SPECIES_ID}_{fam}_F_female"
        male = f"{SPECIES_ID}_{fam}_M_male"
        for child in child_list:
            offspring = f"{SPECIES_ID}_{fam}_{child}_offspring"
            for chrom in sp_chrom_dict_autosome[SPECIES_ID]:
                for minDP in [20,22,24,25,26,28,30]:
                    gwf.target_from_template(
                        name = f"{offspring}_{chrom}_minDP_{minDP}_callable_sites",
                        template = callable_sites(
                        work_path = WORK_DIR,
                        script_path = HELP_SCRIPTS_PATH,
                        data_path = DATA_DIR,
                        log_path = LOG_DIR,
                        sp = SPECIES_ID,
                        female = female,
                        male = male,
                        offspring = offspring,
                        chrom = chrom,
                        minDP = minDP))
                    gwf.target_from_template(
                        name = f"{offspring}_{chrom}_minDP_{minDP}_ind_het",
                        template= ind_het_est(
                        work_path = WORK_DIR,
                        script_path =  HELP_SCRIPTS_PATH,
                        log_path = LOG_DIR,
                        sp = SPECIES_ID,
                        offspring = offspring,
                        chrom = chrom,
                        minDP = minDP)
                        )
                    gwf.target_from_template(
                        name = f"{offspring}_{chrom}_minDP_{minDP}_DNM_scan",
                        template = dnm_scan(
                        work_path = WORK_DIR,
                        script_path = HELP_SCRIPTS_PATH,
                        data_path = DATA_DIR,
                        log_path = LOG_DIR,
                        sp = SPECIES_ID,
                        female = female,
                        male = male,
                        offspring = offspring,
                        chrom = chrom,
                        minDP = minDP))
                    gwf.target_from_template(
                        name = f"{offspring}_{chrom}_minDP_{minDP}_het_AD",
                        template = het_AD_retrive(
                        work_path = WORK_DIR,
                        script_path = HELP_SCRIPTS_PATH,
                        data_path = DATA_DIR,
                        log_path = LOG_DIR,
                        female = female,
                        male = male,
                        offspring = offspring,
                        chrom = chrom,
                        sp = SPECIES_ID,
                        minDP = minDP)
                        )
                    gwf.target_from_template(
                        name = f"{offspring}_{chrom}_minDP_{minDP}_DNM_check",
                        template = dnm_check(
                        work_path = WORK_DIR,
                        script_path = HELP_SCRIPTS_PATH,
                        data_path = DATA_DIR,
                        log_path = LOG_DIR,
                        sp = SPECIES_ID,
                        female = female,
                        male = male,
                        offspring = offspring,
                        chrom = chrom,
                        minDP = minDP,
                        bam_path = BAM_PATH,
                        ref = REF_FASTA))
    for  minDP in [20,22,24,25,26,28,30]:
        gwf.target_from_template(
            name = f"{SPECIES_ID}_minDP_{minDP}_IGV_script",
            template = dnm_igv_script(
            work_path = WORK_DIR.replace("/home/jilong","/Users/au688344/GenomeDK"),
            script_path = HELP_SCRIPTS_PATH.replace("/home/jilong","/Users/au688344/GenomeDK"),
            bam_path = BAM_PATH.replace("/home/jilong","/Users/au688344/GenomeDK"),
            log_path = LOG_DIR.replace("/home/jilong","/Users/au688344/GenomeDK"),
            sp = SPECIES_ID,
            minDP = minDP
            )
        )
    return gwf
