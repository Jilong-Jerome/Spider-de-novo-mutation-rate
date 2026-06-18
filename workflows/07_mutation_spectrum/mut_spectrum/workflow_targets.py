from workflow_templates import *
import os
from gwf import *

def run_class_dnm(gwf, work_path, script_path, log_path, sp, minDP, genome, tsv):
    gwf.target_from_template(
            name = f"classify_DNM_{sp}_minDP_{minDP}",
            template = classify_DNM(work_path,script_path,log_path,sp,minDP,genome,tsv)
            )


