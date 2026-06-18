from workflow_templates import *
from workflow_targets import *
from gwf import *

gwf = Workflow()

for sp in ["BIC","SAR","DUM","TEN","MIM","AFR"]:
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/kinship/{sp}".format(sp=sp) 
    vcf = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/gatk/genotype/{sp}/{sp}_GATK.vcf.gz".format(sp=sp)
    cal_IBD_snps(gwf,path,sp,vcf)
