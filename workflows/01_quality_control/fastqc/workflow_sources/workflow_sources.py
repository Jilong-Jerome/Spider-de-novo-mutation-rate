from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *
def run_fastqc(gwf,work_path,data_path,log_path,ind,sp):
    gwf.target_from_template(
        name = f"{ind}_fastqc",
        template = fastqc_fastq(work_path,data_path,log_path,ind,sp)
        )
    return gwf
