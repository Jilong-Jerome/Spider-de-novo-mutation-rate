from gwf import *
from workflow_templates import *
from workflow_targets import *
from workflow_dicts import *

gwf = Workflow()

# Germline call
for sp in ["LIN","SAR","BIC","DUM","TEN","MIM","AFR"]:
    if sp in ["LIN","MIM","AFR"]:
       fam_list = ["family1","family2","family3","family4"]
    elif sp == "BIC":
       fam_list = ["family1","family2","family3","family4","family5","family6"]
    else:
       fam_list = ["family1","family2","family3","family4","family5"]
    for fam in fam_list :#,"family2","family3","family4","family5"]:
        if sp == "MIM" and fam == "family4":
                child_list = ["S1","S2","S3","S4"]
        elif sp == "DUM" and fam == "family3":
                child_list = ["S1","S2","S4","S5","S6"]
        else:
                child_list = ["S1","S2","S3","S4","S5","S6"]
        for child in child_list:
            female = "{sp}_{fam}_F_female".format(sp=sp,fam=fam)
            male = "{sp}_{fam}_M_male".format(sp=sp,fam=fam)
            offspring = "{sp}_{fam}_{child}_offspring".format(sp=sp,fam=fam,child=child)
            run_callable_sites_per_trio(gwf,sp,female,male,offspring)
            run_dnm_scan(gwf,sp,female,male,offspring)
            run_bam_check(gwf,sp,female,male,offspring)

# IGV check
for sp in ["LIN","BIC","AFR","MIM","TEN","DUM","SAR"]:
    run_igv_check(gwf,sp)
