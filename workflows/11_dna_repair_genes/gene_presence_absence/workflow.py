#!/bin/env python3
import sys, os, yaml, glob
from gwf import *
sys.path.insert(0, os.path.realpath('/faststorage/project/spider2/social_spiders_2020/people/jilong/chapter3_dnm/script/hyper/workflow_sources/'))
from workflow_sources import *
configs = glob.glob('./configs/*config.y*ml')
gwf = Workflow()
for config in configs:
    gwf = gene_variant_check_workflow(config, gwf)
