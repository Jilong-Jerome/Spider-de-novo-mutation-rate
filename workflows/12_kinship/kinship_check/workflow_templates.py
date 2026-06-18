from gwf import *
LOG_PATH = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/logs/kinship"
def kinship_vcf_filter(sp,path,sp_vcf):
    """ Template for checking kinship among individulas in a vcf file"""
    inputs = [sp_vcf]
    outputs = [LOG_PATH+'/{sp}_kinship_filter.DONE'.format(sp=sp)]
    options = {
               'cores': 1,
               'memory': '12g',
               'walltime':"12:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate vcftools
    echo jobinfo $SLURM_JOBID
    echo "start finding X kinship, filter vcf files for biallelic sites, minor allele frequency and LD prunening"
    date
    mkdir -p {path}
    cd {path}
    vcftools --gzvcf {sp_vcf} --min-alleles 2 --max-alleles 2 --remove-filtered-all --recode --recode-INFO-all --out {sp}_biallelic
    conda activate pysam
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/kinship_check/vcf_filter_maf_snps.py {sp}_biallelic.recode.vcf {sp}_Autosome_SNPs_maf_filtered.vcf
    rm {sp}_biallelic.recode.vcf
    #conda activate plink
    #plink --vcf {sp}_Autosome_SNPs_maf_filtered.vcf --double-id --allow-extra-chr --genome -indep-pairwise 50 5 0.2 --out {sp}_Autosome_SNPs_maf_filtered
    echo done > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(sp=sp,path = path,sp_vcf=sp_vcf,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def kinship_prune_and_test(sp,path):
    """ Template for checking kinship among individulas in a vcf file"""
    inputs = [LOG_PATH+'/{sp}_kinship_filter.DONE'.format(sp=sp)]
    outputs = [LOG_PATH+'/{sp}_kinship_test.DONE'.format(sp=sp)]
    options = {
               'cores': 1,
               'memory': '4g',
               'walltime':"4:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    #python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/kinship_check/vcf_filter_maf_snps.py {sp}_biallelic.recode.vcf {sp}_Autosome_SNPs_maf_filtered.vcf
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    conda activate plink
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/kinship_check/vcf_fill_id.py {sp}_Autosome_SNPs_maf_filtered.vcf {sp}_Autosome_SNPs_maf_filtered_named.vcf
    plink --vcf {sp}_Autosome_SNPs_maf_filtered_named.vcf --double-id --allow-extra-chr -indep-pairwise 50 5 0.2 --out {sp}_Autosome_SNPs_maf_filtered_named
    plink --vcf {sp}_Autosome_SNPs_maf_filtered_named.vcf --double-id --allow-extra-chr --extract {sp}_Autosome_SNPs_maf_filtered_named.prune.in --genome --out {sp}_kinship
    echo done > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(sp=sp,path = path,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

