from gwf import *
from workflow_templates import *
from workflow_targets import *
from workflow_dicts import *

gwf = Workflow()

### Genome indexing
for sp in ["LIN","SARA","BI","MIM","DUM","TENT","AFR"]:
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/align/index"
    if sp == "AFR":
        genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/data/genome/{sp}_ncbi_chromosome.fa".format(sp=sp)
    else:
        genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/data/genome/{sp}_ncbi_chromosome.fa".format(sp=sp)
    outname = sp
    index_genome(gwf,path,genome,outname)

### Aligning
for sp in ["LIN","SAR","BIC","DUM","TEN","MIM"]:
    fam_list = ["family1","family2","family3","family4","family5"]
    if sp == "BIC":
        fam_list = ["family1","family2","family3","family4","family5","family6"]
    for fam in fam_list:
        for child in ["S1","S2","S3","S4","S5","S6"]:
            indname = sp+"_"+fam+"_"+child+"_"+"offspring"
            align_ind(gwf,indname)
        male = sp+"_"+fam+"_"+"M"+"_"+"male"
        align_ind(gwf,male)
        female = sp+"_"+fam+"_"+"F"+"_"+"female"
        align_ind(gwf,female)
for sp in ["AFR"]:
    for fam in ["family1","family2","family3","family4"]:
        for child in ["S1","S2","S3","S4","S5","S6"]:
            indname = sp+"_"+fam+"_"+child+"_"+"offspring"
            align_ind(gwf,indname)
        male = sp+"_"+fam+"_"+"M"+"_"+"male"
        align_ind(gwf,male)
        female = sp+"_"+fam+"_"+"F"+"_"+"female"
        align_ind(gwf,female)
#for sp in ["MIM","SAR","BIC","DUM","TEN"]:
#    for fam in ["family1","family2","family3","family4","family5"]:
#        male = sp+"_"+fam+"_"+"M"+"_"+"male"
#        align_ind(gwf,male)

### GATK
for sp in ["LIN","SAR","BIC","TEN","DUM","MIM","AFR"]:
    if sp in ["LIN","MIM","AFR"]:
        fam_list = ["family1","family2","family3","family4"]
    elif sp in ["BIC"]:
        fam_list = ["family1","family2","family3","family4","family5","family6"]
    else:
        fam_list = ["family1","family2","family3","family4","family5"] 
    for fam in fam_list:
        for child in ["S1","S2","S3","S4","S5","S6"]:
            indname = sp+"_"+fam+"_"+child+"_"+"offspring"
            gatk_gvcf_ind(gwf,indname)
        male = sp+"_"+fam+"_"+"M"+"_"+"male"
        gatk_gvcf_ind(gwf,male)
        female = sp+"_"+fam+"_"+"F"+"_"+"female"
        gatk_gvcf_ind(gwf,female)
for sp in ["LIN","SAR","BIC","DUM","TEN","MIM","AFR"]:
    gatk_DB_build(gwf,sp) #Build species GATK databse per chromosome for joint genotyping
for sp in ["LIN","SAR","BIC","DUM","TEN","MIM","AFR"]:
    gatk_finalise(gwf,sp) #Combine the vcf files called from GATK pipeline per specei
    mean_DP_per_chrom_per_ind(gwf,sp)
    run_chrom_split(gwf,sp)
    run_chrom_vcf_gzip(gwf,sp)
