#!/bin/env python3
import gwf
import sys, os, yaml, glob
sys.path.insert(0, os.path.realpath('/faststorage/project/spider2/social_spiders_2020/people/jilong/chapter3_dnm/script/germline_call/workflow_X/workflow_sources/'))
from workflow_sources import *
configs = glob.glob('./configurations/*config.y*ml')

gwf = Workflow()
for config in configs:
    gwf = call_dnm_X_workflow(config_file = config,gwf = gwf)
