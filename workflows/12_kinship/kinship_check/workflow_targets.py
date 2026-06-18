from workflow_templates import *
from gwf import *

def cal_IBD_snps(gwf,path,sp,vcf):
    gwf.target_from_template(
        name = "IBD_snps_{sp}".format(sp=sp),
        template = kinship_vcf_filter(sp,path,vcf)
        )
    gwf.target_from_template(
        name = "IBD_cal_{sp}".format(sp=sp),
        template = kinship_prune_and_test(sp,path) 
        )


