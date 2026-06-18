#!/bin/env python3
from gwf import AnonymousTarget
import os, glob

def callable_sites(work_path,script_path,data_path,log_path,sp,female,male,offspring,chrom,minDP):
    inputs = {"vcf":f"{data_path}/germline/{sp}/chrom_splits/{sp}_{chrom}.vcf.gz",# The vcf called per species split by chromosome, gzipped
            "depth":f"{data_path}/{sp}_per_ind_per_chrom_DP.tsv",
            "sex":f"{data_path}/offspring_sex_expanded.tsv"}
    outputs = {"callable_vcf":f"{work_path}/{sp}/callable/minDP_{minDP}/{offspring}/chrom_split/{offspring}_{chrom}_callable_sites.vcf",
            "callable_tsv":f"{work_path}/{sp}/callable/minDP_{minDP}/{offspring}/chrom_split/{offspring}_{chrom}_callable_sites.tsv",
            "log": f"{log_path}/germline_step2/minDP_{minDP}/{offspring}_{chrom}_minDP_{minDP}_callable_sites.DONE"}
    options = {
              'cores': 1,
              'memory': '4g',
              'walltime':'24:00:00',
              'account':"spider2"
              }
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate vcftools
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {log_path}/germline_step2/minDP_{minDP}
    mkdir -p {work_path}/{sp}/callable/minDP_{minDP}/{offspring}/chrom_split
    cd {work_path}/{sp}/callable/minDP_{minDP}/{offspring}/chrom_split

    cp {inputs["vcf"]} /scratch/$SLURM_JOBID/{sp}_{chrom}.vcf.gz
    if [ ! -f \"{offspring}_{chrom}_raw.recode.vcf\" ]; then
        echo "Retriving VCF file for the selected trio"
        vcftools --gzvcf /scratch/$SLURM_JOBID/{sp}_{chrom}.vcf.gz --indv {female} --indv {male} --indv {offspring} --chr {chrom} --min-alleles 1 --max-alleles 2  --max-missing 0 --recode --out {offspring}_{chrom}_raw
    else
        echo "VCF file already retrived, go directly for callable site scanning"
    fi
    conda activate pandas
    python {script_path}/callable_GQ_DP_filter_minDP_X.py {offspring}_{chrom}_raw.recode.vcf {data_path}/{sp}_per_ind_per_chrom_DP.tsv {offspring}_{chrom}_callable_sites.vcf {inputs["sex"]} {minDP}
    rm {offspring}_{chrom}_raw.recode.vcf
    python {script_path}/vcf_site_count.py {offspring}_{chrom}_callable_sites.vcf {offspring}_{chrom}
    echo done > {outputs["log"]}
    date
    jobinfo $SLURM_JOBID
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def dnm_scan(work_path,script_path,data_path,log_path,sp,female,male,offspring,chrom,minDP):
    inputs = {"callable_sites":f"{work_path}/{sp}/callable/minDP_{minDP}/{offspring}/chrom_split/{offspring}_{chrom}_callable_sites.vcf",
            "callable_log":f"{log_path}/germline_step2/minDP_{minDP}/{offspring}_{chrom}_minDP_{minDP}_callable_sites.DONE",
            "sex":f"{data_path}/offspring_sex_expanded.tsv"}
    outputs = {
            "log":f"{log_path}/germline_step2/minDP_{minDP}/{offspring}_{chrom}_minDP_{minDP}_DNM_scan.DONE",
            "vcf_out":f"{work_path}/{sp}/DNMs/minDP_{minDP}/{offspring}/chrom_splits/{offspring}_{chrom}_putative_DNMs_distance_filtered_AB.vcf"
            }
    options = {
              'cores': 1,
              'memory': '4g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate vcftools
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {work_path}/{sp}/DNMs/minDP_{minDP}/{offspring}/chrom_splits
    cd {work_path}/{sp}/DNMs/minDP_{minDP}/{offspring}/chrom_splits
    vcftools --vcf {inputs["callable_sites"]} --remove-indels --recode --recode-INFO-all --min-alleles 2 --max-alleles 2 --out {offspring}_{chrom}_SNPs
    conda activate pyvcf
    python {script_path}/germline_scan.py {offspring}_{chrom}_SNPs.recode.vcf {male} {female} {offspring} {inputs["sex"]} {offspring}_{chrom}_putative_DNMs.vcf
    python {script_path}/distance_filter.py {offspring}_{chrom}_putative_DNMs.vcf {offspring}_{chrom}_putative_DNMs_distance_filtered.vcf
    python {script_path}/germline_scan_AB.py {offspring}_{chrom}_putative_DNMs_distance_filtered.vcf {male} {female} {offspring} {inputs["sex"]} {offspring}_{chrom}_putative_DNMs_distance_filtered_AB.vcf
    echo done > {outputs["log"]}
    date
    jobinfo $SLURM_JOBID
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def dnm_check(work_path,script_path,data_path,bam_path,log_path,sp,ref,female,male,offspring,chrom,minDP):
    inputs = {"vcf":f"{work_path}/{sp}/DNMs/minDP_{minDP}/{offspring}/chrom_splits/{offspring}_{chrom}_putative_DNMs_distance_filtered_AB.vcf",
            "log":f"{log_path}/germline_step2/minDP_{minDP}/{offspring}_{chrom}_minDP_{minDP}_DNM_scan.DONE",
            "sex":f"{data_path}/offspring_sex_expanded.tsv"}
    outputs = {
            "log":f"{log_path}/germline_step2/minDP_{minDP}/{offspring}_{chrom}_minDP_{minDP}_DNM_check.DONE",
            "tsv":f"{work_path}/{sp}/DNMs/minDP_{minDP}/{offspring}/chrom_splits/{offspring}_{chrom}_automate.tsv"
            }
    options = {
              'cores': 1,
              'memory': '4g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bcftools
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {work_path}/{sp}/DNMs/minDP_{minDP}/{offspring}/chrom_splits
    cd {work_path}/{sp}/DNMs/minDP_{minDP}/{offspring}/chrom_splits
    python {script_path}/bcftools_check_region_generate.py {inputs["vcf"]} {offspring}_{chrom}_DNM_regions.tsv
    bcftools mpileup -f {ref} -a AD,DP -R {offspring}_{chrom}_DNM_regions.tsv -Ou {bam_path}/{sp}/{offspring}_final.bam {bam_path}/{sp}/{male}_final.bam {bam_path}/{sp}/{female}_final.bam | bcftools call -c -Ov -o {offspring}_{chrom}_DNM_region_bcftools.vcf
    python {script_path}/germline_scan.py {offspring}_{chrom}_DNM_region_bcftools.vcf {male} {female} {offspring} {inputs["sex"]} {offspring}_{chrom}_DNM_surrounding_region.vcf
    python {script_path}/bcftools_distance_filter.py {inputs["vcf"]} {offspring}_{chrom}_DNM_surrounding_region.vcf {offspring}_{chrom}_DNM_bcftools_dist_checked.vcf
     bcftools mpileup -f {ref} -a AD,DP -R {offspring}_{chrom}_DNM_bcftools_dist_checked.vcf {bam_path}/{sp}/{offspring}_final.bam {bam_path}/{sp}/{male}_final.bam {bam_path}/{sp}/{female}_final.bam| bcftools call -c - | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP\t%AD]\n' > {offspring}_{chrom}_DNM_bcftools_AD_check.tsv
     python {script_path}/bam_check.py {offspring}_{chrom}_DNM_bcftools_AD_check.tsv {inputs["sex"]} {offspring} {offspring}_{chrom}_bam_check_pass.tsv
     python {script_path}/automate_filter_res_parse.py {offspring}_{chrom}_bam_check_pass.tsv {offspring}_{chrom}_automate.tsv {offspring} {male} {female}
     echo done > {outputs["log"]}
    date
    jobinfo $SLURM_JOBID
     """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def dnm_igv_script(work_path,script_path,bam_path,log_path,sp,minDP):
    cluster_work_path = f"{work_path}".replace("/Users/au688344/GenomeDK","/home/jilong")
    cluster_script_path = f"{script_path}".replace("/Users/au688344/GenomeDK","/home/jilong")
    possible_file = f"{work_path}/{sp}/DNMs/minDP_{minDP}/{sp}_minDP_{minDP}_automate.tsv".replace("/Users/au688344/GenomeDK","/home/jilong")
    inputs = [possible_file]
    outputs = {"log":f"{log_path}/germline_step2/minDP_{minDP}/{sp}_minDP_{minDP}_IGV.DONE".replace("/Users/au688344/GenomeDK","/home/jilong")}
    options = {
              'cores': 1,
              'memory': '1g',
              'walltime':'2:00:00',
              'account':"spider2"
              }
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate pyvcf
    echo jobinfo $SLURM_JOBID
    date
    export DISPLAY=localhost:10.0
    mkdir -p {cluster_work_path}/{sp}/pngs/minDP_{minDP}
    cd {cluster_work_path}/{sp}/pngs/minDP_{minDP}
    python {cluster_script_path}/igv_batch_create.py {cluster_work_path}/{sp}/DNMs/minDP_{minDP}/{sp}_minDP_{minDP}_automate.tsv {work_path}/{sp}/pngs/minDP_{minDP} {sp} {bam_path}
    echo done > {outputs["log"]}
    date
    jobinfo $SLURM_JOBID
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
