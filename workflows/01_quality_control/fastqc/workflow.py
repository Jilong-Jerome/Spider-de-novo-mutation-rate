#!/bin/env python3
import sys, os, yaml, glob
from gwf import *
sys.path.insert(0, os.path.realpath('/faststorage/project/spider2/social_spiders_2020/people/jilong/chapter3_dnm/script/fastqc/workflow_sources/'))
from workflow_sources import *

gwf = Workflow()
for filename in ["BIC_ind.txt","SAR_ind.txt","TEN_ind.txt","DUM_ind.txt",
        "AFR_ind.txt","MIM_ind.txt","LIN_ind.txt"]:
    inds = open(filename)
    for line in inds:
        ind = line.strip("\n")
        infos = ind.split("_")
        sp = infos[0]
        work_path = f"/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/fastqc/{sp}"
        data_path = f"/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/data/raw_reads"
        log_path = f"/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/fastqc/logs/{sp}"
        run_fastqc(gwf,work_path,data_path,log_path,ind,sp)
