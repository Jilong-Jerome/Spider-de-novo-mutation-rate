from gwf import *
from workflow_templates import *
from workflow_targets import *

gwf = Workflow()

for sp in ["AFR","MIM","TEN","DUM","BIC","SAR","LIN"]:
    for minDP in [26]:
        if sp == "BIC":
            ref_sp = "BI"
        elif sp == "TEN":
            ref_sp = "TENT"
        elif sp == "SAR":
            ref_sp = "SARA"
        else:
            ref_sp = sp
        script_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/mut_spectrum"
        work_path = f"/faststorage/project/spider2/social_spiders_2020/people/jilong/chapter3_dnm/steps/germline/minDP/{sp}/DNMs/minDP_{minDP}"
        log_path  = "/faststorage/project/spider2/social_spiders_2020/people/jilong/chapter3_dnm/script/mut_spectrum/logs"
        genome = f"/faststorage/project/spider2/social_spiders_2020/people/jilong/chapter3_dnm/data/genome/{ref_sp}_ncbi_chromosome.fa"
        tsv = f"/faststorage/project/spider2/social_spiders_2020/people/jilong/chapter3_dnm/steps/germline/minDP/{sp}/DNMs/minDP_{minDP}/{sp}_minDP_{minDP}_automate.tsv"
        run_class_dnm(gwf, work_path, script_path, log_path, sp, minDP, genome, tsv)
