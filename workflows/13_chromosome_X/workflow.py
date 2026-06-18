from workflow_templates import *
from workflow_targets import *
from gwf import *

gwf = Workflow()
# get coverage information per window in bam files
for sp in ["LIN","DUM","TEN","SAR","MIM","BIC"]:
    if sp == "BIC":
        ref_sp = "BI"
    elif sp == "SAR":
        ref_sp = "SARA"
    elif sp == "TEN":
        ref_sp = "TENT"
    else:
        ref_sp = sp
    indname = "{sp}_family1_M_male".format(sp=sp)
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/align/chromX"
    cal_coverage(gwf,path,indname,sp,ref_sp)
    viz_coverage(gwf,path,indname)

inds = open("offsprings_not_scanned.tsv")
for ind_info in inds:
    indname = ind_info.strip("\n")
    sp = indname.split("_")[0]
    ref_sp = sp
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/align/chromX"
    cal_coverage(gwf,path,indname,sp,ref_sp)
    viz_coverage(gwf,path,indname)




