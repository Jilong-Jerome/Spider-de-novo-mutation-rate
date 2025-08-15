#!/bin/env python3
from gwf import AnonymousTarget
import os, glob

def fastqc_fastq(work_path,data_path,log_path,ind,sp):
    inputs = {"R1":f"{data_path}/{sp}/{ind}_R1.fq.gz",
            "R2":f"{data_path}/{sp}/{ind}_R2.fq.gz"}
    outputs = {
            "log":f"{log_path}/{ind}_fastqc.DONE",
            }
    options = {
              'cores':4,
              'memory':'36g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate fastqc
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {work_path}/{sp}
    cd {work_path}/{sp}
    fastqc -t 4 -o {work_path}/{sp} {inputs["R1"]} {inputs["R2"]}
    mkdir -p {log_path}
    echo done > {log_path}/{ind}_fastqc.DONE
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

