from workflow_dicts import *
from gwf import *
import os
LOG_PATH = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/logs"

def callable_sites(path,sp,input_vcf,female,male,offspring,dp_path,log_list,chrom):
    inputs = [LOG_PATH + "/align_step1/{sp}_{chrom}_vcf_gzipped.DONE".format(sp = sp,chrom=chrom)]#log_list
    outputs = [LOG_PATH + '/germline_step2/{child}_{chrom}_callable_sites.DONE'.format(child = offspring,chrom=chrom)]
    script_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/align"
    options = {
              'cores': 1,
              'memory': '4g',
              'walltime':'24:00:00',
              'account':"spider2"
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate vcftools
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}

    cp {vcf} /scratch/$SLURM_JOBID/{sp}_{chrom}.vcf.gz

    vcftools --gzvcf /scratch/$SLURM_JOBID/{sp}_{chrom}.vcf.gz --indv {female} --indv {male} --indv {child} --chr {chrom} --min-alleles 1 --max-alleles 2  --max-missing 0 --recode --out {child}_{chrom}_raw

    /home/jilong/anaconda3/bin/python {script_path}/callable_GQ_DP_filter.py {child}_{chrom}_raw.recode.vcf /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/germline_call/{sp}_per_ind_per_chrom_DP.tsv {child}_{chrom}_callable_sites.vcf
    rm {child}_{chrom}_raw.recode.vcf
    python {script_path}/vcf_site_count.py {child}_{chrom}_callable_sites.vcf {child}_{chrom}
    echo done > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(sp=sp,path=path,script_path=script_path,vcf=input_vcf,female=female,male=male,child= offspring,chrom=chrom,log =outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def dnm_scan(path,sp,input_vcf,female,male,offspring,chrom):
    inputs = [LOG_PATH + '/germline_step2/{child}_{chrom}_callable_sites.DONE'.format(child = offspring,chrom=chrom)]
    outputs = [LOG_PATH + '/germline_step2/{child}_{chrom}_DNM.DONE'.format(child = offspring,chrom=chrom)]
    options = {
              'cores': 1,
              'memory': '4g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate vcftools
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    vcftools --vcf {vcf} --remove-indels --recode --recode-INFO-all --min-alleles 2 --max-alleles 2 --out {child}_{chrom}_SNPs
    conda activate pyvcf
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/align/germline_scan.py {child}_{chrom}_SNPs.recode.vcf {father} {mother} {child} {child}_{chrom}_putative_DNMs.vcf 
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/align/distance_filter.py {child}_{chrom}_putative_DNMs.vcf {child}_{chrom}_putative_DNMs_distance_filtered.vcf 
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/align/germline_scan_AB.py {child}_{chrom}_putative_DNMs_distance_filtered.vcf {father} {mother} {child} {child}_{chrom}_putative_DNMs_distance_filtered_AB.vcf 
    echo done > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,vcf=input_vcf,child=offspring,father = male, mother = female,log = outputs[0],chrom=chrom)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def dnm_bcftools_check(path,sp,ref,chrom,vcf,female,male,offspring):
    inputs = [LOG_PATH + '/germline_step2/{child}_{chrom}_DNM.DONE'.format(child = offspring,chrom=chrom)]
    outputs = [LOG_PATH + '/germline_step2/{child}_{chrom}_DNM_check.DONE'.format(child = offspring,chrom=chrom)]
    options = {
              'cores': 1,
              'memory': '4g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bcftools
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/align/bcftools_check_region_generate.py {vcf} {child}_{chrom}_DNM_regions.tsv
    bcftools mpileup -f {ref} -a AD,DP -R {child}_{chrom}_DNM_regions.tsv -Ou /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/align/bams/{sp}/{child}_final.bam /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/align/bams/{sp}/{father}_final.bam /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/align/bams/{sp}/{mother}_final.bam | bcftools call -c -Ov -o {child}_{chrom}_DNM_region_bcftools.vcf
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/align/germline_scan.py {child}_{chrom}_DNM_region_bcftools.vcf {father} {mother} {child} {child}_{chrom}_DNM_surrounding_region.vcf
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/align/bcftools_distance_filter.py {vcf} {child}_{chrom}_DNM_surrounding_region.vcf {child}_{chrom}_DNM_bcftools_dist_checked.vcf
    bcftools mpileup -f {ref} -a AD,DP -R {child}_{chrom}_DNM_bcftools_dist_checked.vcf /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/align/bams/{sp}/{child}_final.bam /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/align/bams/{sp}/{father}_final.bam /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/align/bams/{sp}/{mother}_final.bam| bcftools call -c - | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP\t%AD]\n' > {child}_{chrom}_DNM_bcftools_AD_check.tsv
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/align/bam_check.py {child}_{chrom}_DNM_bcftools_AD_check.tsv {child}_{chrom}_bam_check_pass.tsv
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/align/automate_filter_res_parse.py {child}_{chrom}_bam_check_pass.tsv {child}_{chrom}_automate.tsv {child} {father} {mother}
    echo done > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,sp=sp,ref=ref,vcf=vcf,child=offspring,father = male, mother = female,log = outputs[0],chrom=chrom)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def create_IGV_check(sp,path):
    inputs = ["/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/germline/{sp}/DNMs/{sp}_DNMs_automate.tsv".format(sp=sp)]
    outputs = [LOG_PATH + '/germline_step2/{sp}_IGV_check.DONE'.format(sp=sp)]
    options = {
              'cores': 1,
              'memory': '32g',
              'walltime':'8:00:00',
              'account':"spider2"
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate pyvcf
    echo jobinfo $SLURM_JOBID
    date
    export DISPLAY=localhost:10.0
    mkdir -p {path}/pngs
    cd {path}
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/germline_call/igv_batch_create.py /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/germline/{sp}/DNMs/{sp}_DNMs_automate.tsv {path}/pngs {sp} /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/align/bams 
    /home/jilong/software/IGV_Linux_2.9.4/igv.sh -b {path}/igv_{sp}_script.txt
    echo done > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,sp=sp,log = outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

   
