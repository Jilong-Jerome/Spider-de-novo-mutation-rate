from workflow_templates import *
import os
from gwf import *

def run_asm_stat(gwf, work_path, script_path, log_path, sp, genome, bed):
    gwf.target_from_template(
            name = f"genome_stat_{sp}",
            template = genome_composition(work_path,script_path,log_path,sp,genome,bed) 
            )

def run_callable_stat(gwf, work_path, script_path, log_path, sp, ind, chrom, vcf):
    gwf.target_from_template(
            name = f"callable_stat_{ind}_{chrom}",
            template = split_callable_genome_ind(work_path,script_path,log_path,sp,ind,chrom,vcf) 
            )

def run_callable_split_count(gwf,work_path,script_path,log_path,sp,ind,chrom,vcf,genome):
    gwf.target_from_template(
            name = f"callable_split_count_{ind}_{chrom}",
            template = split_callable_genome_count(work_path,script_path,log_path,sp,ind,chrom,vcf,genome)
            )
