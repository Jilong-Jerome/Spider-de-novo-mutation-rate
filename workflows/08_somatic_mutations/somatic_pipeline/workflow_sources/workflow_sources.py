from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def somatic_mutation_workflow(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    """
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    CONFIG = yaml.safe_load(open(config_file))
    ACCOUNT: str = CONFIG['account']
    DIPLOID_ASSEMBLY: str = CONFIG['assembly_file_path']
    PACB_READS: str = CONFIG['pacbio_seq_file_path']
    TASK_ID: str = CONFIG['task_id']
    SPECIES_ID: str = CONFIG['species_id']
    WORK_DIR: str = CONFIG['working_directory_path']
    OUTPUT_DIR: str = CONFIG['output_directory_path']
    LOG_DIR: str = CONFIG['log_directory_path']


    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------

    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    # Align pacb reads to diploid assemblies #
    gwf.target_from_template(
        name = f'pbmm2_align_{TASK_ID}',
        template = diploid_pbmm(
            work_path = WORK_DIR,
            asm_fa = DIPLOID_ASSEMBLY,
            pacb_file = PACB_READS,
            name = TASK_ID,
            spid = SPECIES_ID,
            log_path = LOG_DIR
        )
    )
    return gwf
