from workflow_dicts import *
from gwf import *
import os
LOG_PATH = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/logs/align_step1"
def bwa2_index(ref_genome,path,basename):
    """ Template for indexing a genome with bwa2-mem"""
    inputs = [ref_genome]
    outputs = [LOG_PATH+"/{basename}_index.DONE".format(basename=basename)]
    options = {
               'cores': 12,
               'memory': '128g',
               'walltime':"12:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bwa2
    echo jobinfo $SLURM_JOBID
    echo "start indexing bwa2"
    date
    mkdir -p {path}
    cd {path}
    #/home/jilong/software/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index -p {basename} {ref_genome}
    bwa-mem2 index -p {basename} {ref_genome}
    echo "finish indexing bwa2"
    echo done > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(ref_genome=ref_genome,path=path,basename=basename,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def bwa2_align(path,sp,indname,read1,read2):
    """ Template for aligning paired read to a genome with bwa2-mem"""
    inputs = [LOG_PATH+"/{sp}_index.DONE".format(sp=sp)]
    outputs = [LOG_PATH+'/{indname}_align_s0.DONE'.format(indname=indname)]
    options = {
               'cores': 20,
               'memory': '64g',
               'walltime':"24:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    #conda activate spider
    conda activate bwa2
    echo jobinfo $SLURM_JOBID
    echo "start aligning bwa-mem2"
    date
    mkdir -p {path}
    cd {path}
    bwa-mem2 mem -t 20 /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/align/index/{sp} {read1} {read2} | samtools view -Sb > {indname}_s0.bam
    cp {indname}_s0.bam {indname}_s0_backup.bam
    echo "finish aligning with bwa-mem2"
    echo done > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(sp=sp,indname=indname,path=path,read1=read1,read2=read2,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def picard_addRG(path,indname):
    """GATK4.0 preprocessing, read in bam file need RGs"""
    inputs = [LOG_PATH+'/{indname}_align_s0.DONE'.format(indname=indname)]
    outputs = [LOG_PATH+'/{indname}_addRG_s1.DONE'.format(indname=indname)]
    options = {
              'cores': 8,
              'memory': '16g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = """    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate spider
    echo "start gatk markduplicatesSpark"
    echo jobinfo $SLURM_JOBID
    date
    cd {path}
    java -jar /home/jilong/software/picard.jar AddOrReplaceReadGroups --I {indname}_s0.bam --O {indname}_s1.bam -RGID {indname} -RGPU unknown -RGSM {indname} -RGPL illumina -RGLB lib0
    rm {indname}_s0.bam
    echo "finished gatk markduplicatesSpark"
    date
    jobinfo $SLURM_JOBID
    echo done > {log}
    """.format(path=path,log=outputs[0],indname=indname)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def samtool_fixmate(path,indname):
    """Samtools fixtmate for remove duplicates"""
    inputs = [LOG_PATH+'/{indname}_addRG_s1.DONE'.format(indname=indname)]
    outputs = [LOG_PATH+'/{indname}_fixmate_s2.DONE'.format(indname=indname)]
    options = {
              'cores': 16,
              'memory': '32g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate spider
    echo "start samtools fixmate, remove secondary and unmapped reads"
    echo jobinfo $SLURM_JOBID
    date
    cd {path}
    samtools fixmate -rm {indname}_s1.bam {indname}_s2.bam -@ 16
    rm {indname}_s1.bam
    echo "finished samtools fixmate"
    date
    jobinfo $SLURM_JOBID
    echo done > {log}
    """.format(path=path,indname=indname,log = outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def samtool_sort(path,indname):
    """Samtools sort bam file"""
    inputs = [LOG_PATH+'/{indname}_fixmate_s2.DONE'.format(indname=indname)]
    outputs = [LOG_PATH+'/{indname}_sort_s3.DONE'.format(indname=indname)]
    options = {
              'cores': 16,
              'memory': '32g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate spider
    echo "start samtools sort bam"
    echo jobinfo $SLURM_JOBID
    date
    cd {path}
    samtools sort {indname}_s2.bam -o {indname}_s3.bam -@ 16
    echo "finished samtools sort"
    date
    rm {indname}_s2.bam
    echo done > {log}
    jobinfo $SLURM_JOBID
    """.format(path=path,indname=indname,log = outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def samtool_markdup(path,indname):
    """Samtools markdup for remove duplicates"""
    inputs = [LOG_PATH+'/{indname}_sort_s3.DONE'.format(indname=indname)]
    outputs = [LOG_PATH+'/{indname}_markdup_s4.DONE'.format(indname=indname)]
    options = {
              'cores': 16,
              'memory': '32g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate spider
    echo "start samtools markdup to remove duplicate reads"
    echo jobinfo $SLURM_JOBID
    date
    cd {path}
    samtools markdup -r -f {indname}.stat -s  {indname}_s3.bam {indname}_s4.bam -@ 16
    echo "finished samtools markdup"
    rm {indname}_s3.bam
    date
    jobinfo $SLURM_JOBID
    echo done > {log}
    """.format(path=path,indname=indname,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def samtool_mapped(path,indname):
    """Samtools take only the mapped reads"""
    inputs = [LOG_PATH+'/{indname}_markdup_s4.DONE'.format(indname=indname)]
    outputs = [LOG_PATH+'/{indname}_MQ_s5.DONE'.format(indname=indname)]
    options = {
              'cores': 16,
              'memory': '32g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = """    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate spider
    echo "start samtools to remove unmapped reads"
    echo jobinfo $SLURM_JOBID
    date
    cd {path}
    samtools view -@ 16 -bq 60  -f 0x2 -F 0x4 {indname}_s4.bam > {indname}_final.bam
    rm {indname}_s4.bam
    echo "finished samtools get mapped reads"
    date
    echo done > {log}
    jobinfo $SLURM_JOBID
    """.format(path=path,indname=indname,log = outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def picard_bai(path,indname):
    """GATK4.0 preprocessing, index bam files"""
    inputs = [LOG_PATH+'/{indname}_MQ_s5.DONE'.format(indname=indname)]
    outputs = [LOG_PATH+'/{indname}_bai_s6.DONE'.format(indname=indname)]
    options = {
              'cores': 12,
              'memory': '32g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate spider
    echo "start picar build bam index"
    echo jobinfo $SLURM_JOBID
    date
    cd {path}
    java -jar /home/jilong/software/picard.jar BuildBamIndex --INPUT {indname}_final.bam --OUTPUT {indname}_final.bam.bai
    echo "finished picard index bam"
    date
    echo done > {log}
    jobinfo $SLURM_JOBID
    """.format(path=path,indname=indname,log = outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def gatk_haplotype_call(path,bam,assembly,indname,chrom):
    """GATK4 haplotype caller by each sample each chrom"""
    inputs = [LOG_PATH+'/{indname}_bai_s6.DONE'.format(indname=indname)]
    outputs = [LOG_PATH+'/{indname}_{chrom}_gvcf.DONE'.format(indname=indname,chrom=chrom)]
    options = {
              'cores': 16,
              'memory': '36g',
              'walltime':'36:00:00',
              'account':"spider2"
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate spider
    echo "start gatk haplotype caller by sample with bam output"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    /home/jilong/software/gatk-4.2.0.0/gatk --java-options "-Xmx36g" HaplotypeCaller -R {assembly} -I {ibam} -O {indname}_{chrom}.g.vcf -L {chrom} -ERC BP_RESOLUTION -bamout {indname}_{chrom}_remap.bam --tmp-dir /scratch/$SLURM_JOBID/  --native-pair-hmm-threads 16 --disable-tool-default-read-filters true --max-reads-per-alignment-start 0 --do-not-run-physical-phasing
    echo "finished gatk call vcf for sample" > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,assembly=assembly,ibam=bam,indname=indname,chrom=chrom,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def merge_gvcf_by_chrom(log_list,merge_list,path,indname):
    inputs = log_list
    outputs = [LOG_PATH+'/{indname}_gvcf.DONE'.format(indname=indname)]
    options = {
              'cores': 6,
              'memory': '24g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    concat_string = ""
    for chrom_file in merge_list:
        concat_string = concat_string+chrom_file+" "
    spec = """    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bcftools
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    bcftools concat {concat_string} -o {indname}.g.vcf -O v --threads 6
    conda activate gatk4
    gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx24g -Xms24g" IndexFeatureFile --input {indname}.g.vcf
    echo done > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,concat_string=concat_string,indname = indname, log = outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def merge_bam_by_chrom(log_list,merge_list,path,indname):
    inputs = log_list
    outputs = [LOG_PATH+'/{indname}_remap_bam.DONE'.format(indname=indname)]
    options = {
              'cores': 6,
              'memory': '24g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    concat_string = ""
    for chrom_file in merge_list:
        concat_string = concat_string+chrom_file+" "
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate samtools117
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    samtools merge {indname}_remap.bam {concat_string} -@ 6
    samtools index -b {indname}_remap.bam {indname}_remap.bam.bai -@ 6
    date
    echo done > {log}
    jobinfo $SLURM_JOBID
    """.format(path=path,concat_string=concat_string,indname=indname,log = outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def gatk_consolidate(log_list,ref,vcf_path,sp,interval_range,db_id,batch_size,sample_map):
    inputs = log_list
    outputs = [LOG_PATH + "/gatk_geno/GATK4_genotype_segment_{db_id}.DONE".format(db_id=db_id)]
    options = {
              'cores': 6,
              'memory': '64g',
              'walltime':'2:00:00',
              'account':"spider2"
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate gatk4
    mkdir -p {vcf_path}
    echo jobinfo $SLURM_JOBID
    gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx60g -Xms60g" GenomicsDBImport \
            --genomicsdb-workspace-path /scratch/$SLURM_JOBID/{db_id} \
            --batch-size {batch_size} \
            --sample-name-map {sample_map} \
            --tmp-dir /scratch/$SLURM_JOBID/ \
            --reader-threads 6 \
            -L {interval} \
            --genomicsdb-shared-posixfs-optimizations true \
            --genomicsdb-vcf-buffer-size 4194304
    date
    gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx60g -Xms60g" GenotypeGVCFs \
            -R {ref} \
            -V gendb:///scratch/$SLURM_JOBID/{db_id} \
            --tmp-dir /scratch/$SLURM_JOBID/ \
            -O /scratch/$SLURM_JOBID/{db_id}.vcf.gz \
            -L {interval} \
            -all-sites true
    gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx60g -Xms60g" VariantFiltration \
            -V /scratch/$SLURM_JOBID/{db_id}.vcf.gz \
            -O /scratch/$SLURM_JOBID/{db_id}_GATK_filtered.vcf.gz \
            --filter-name "gatk_germline" \
            --tmp-dir /scratch/$SLURM_JOBID/ \
            --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
            --verbosity ERROR
    mv /scratch/$SLURM_JOBID/{db_id}_GATK_filtered.vcf.gz {vcf_path}
    jobinfo $SLURM_JOBID
    echo done > {log}
    """.format(db_id=db_id,ref=ref,interval=interval_range,vcf_path=vcf_path,batch_size=batch_size,sample_map=sample_map,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def merge_vcfs(log_list,merge_list,path,outname):
    inputs = log_list
    outputs = [LOG_PATH+'/{outname}_vcf.DONE'.format(outname=outname)]
    options = {
              'cores': 6,
              'memory': '24g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    concat_string = ""
    for vcf_file in merge_list:
        concat_string = concat_string + vcf_file + " "
    spec = """    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bcftools
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    bcftools concat {concat_string} -o {outname}.vcf.gz -O v --threads 6
    conda activate gatk4
    gatk --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx24g -Xms24g" IndexFeatureFile --input {outname}.vcf.gz
    echo done > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,concat_string=concat_string,outname = outname, log = outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)


def vcf_depths(sp,vcf,chrom,path):
    inputs = [LOG_PATH+'/{sp}_GATK_vcf.DONE'.format(sp=sp)]
    outputs = [LOG_PATH + '{sp}_{chrom}_depth.DONE'.format(sp=sp,chrom=chrom)]
    options = {
              'cores': 1,
              'memory': '1g',
              'walltime':'6:00:00',
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
    ln -s -f {vcf} .
    vcftools --gzvcf {vcf} --depth --chr {chrom} --out {sp}_{chrom}
    echo done > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(sp=sp,log = outputs[0],path=path,vcf = vcf,chrom=chrom)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def sp_chrom_vcf_split(vcf,path,sp):
    inputs = [LOG_PATH+'/{sp}_GATK_vcf.DONE'.format(sp=sp)]
    outputs = [LOG_PATH+'/{sp}_GATK_vcf_split.DONE'.format(sp=sp)]
    options = {
              'cores': 1,
              'memory': '1g',
              'walltime':'48:00:00',
              'account':"spider2"
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate pysam
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/align/chrom_split.py {vcf} {path} {sp}  
    echo done > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(sp=sp,log = outputs[0],path=path,vcf = vcf)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def gzip_file(sp,chrom,path,vcf):
    inputs = [LOG_PATH+'/{sp}_GATK_vcf_split.DONE'.format(sp=sp)]
    outputs = [LOG_PATH+'/{sp}_{chrom}_vcf_gzipped.DONE'.format(sp=sp,chrom=chrom)]
    options = {
              'cores': 1,
              'memory': '4g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    gzip {vcf}
    echo done > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(sp=sp,log = outputs[0],path=path,vcf = vcf)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)



def callable_sites(path,sp,input_vcf,female,male,offspring,dp_path,log_list,chrom):
    inputs = [LOG_PATH + "/{sp}_{chrom}_vcf_gzipped.DONE".format(sp = sp,chrom=chrom)]#log_list
    outputs = [LOG_PATH + '/{child}_{chrom}_callable_sites.DONE'.format(child = offspring,chrom=chrom)]
    script_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/align"
    options = {
              'cores': 1,
              'memory': '8g',
              'walltime':'8:00:00',
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
    
    python {script_path}/callable_GQ_DP_filter.py {child}_{chrom}_raw.recode.vcf {script_path}/{sp}_per_ind_per_chrom_DP.tsv {child}_{chrom}_callable_sites.vcf
    rm {child}_{chrom}_raw.recode.vcf
    python {script_path}/vcf_site_count.py {child}_{chrom}_callable_sites.vcf {child}_{chrom}
    echo done > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(sp=sp,path=path,script_path=script_path,vcf=input_vcf,female=female,male=male,child= offspring,chrom=chrom,log =outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)


def dnm_scan(path,sp,input_vcf,female,male,offspring,chrom):
    inputs = [LOG_PATH + '/{child}_{chrom}_callable_sites.DONE'.format(child = offspring,chrom=chrom)]
    outputs = [LOG_PATH + '/{child}_{chrom}_DNM.DONE'.format(child = offspring,chrom=chrom)]
    options = {
              'cores': 1,
              'memory': '4g',
              'walltime':'2:00:00',
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
    inputs = [LOG_PATH + '/{child}_{chrom}_DNM.DONE'.format(child = offspring,chrom=chrom)]
    outputs = [LOG_PATH + '/{child}_{chrom}_DNM_check.DONE'.format(child = offspring,chrom=chrom)]
    options = {
              'cores': 1,
              'memory': '4g',
              'walltime':'4:00:00',
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


