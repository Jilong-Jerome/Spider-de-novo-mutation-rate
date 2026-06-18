from gwf import AnonymousTarget
import os, glob

def diploid_pbmm(work_path: str,asm_fa: str,pacb_file: str, name: str, spid:str,log_path: str):
    """
    """
    inputs = {"pacb":pacb_file,
              "asm":asm_fa}
    outputs = {"log":log_path+"/"+name+"_pbmm_diploid.log"}
    options = {
              'cores':32,
              'memory':'128g',
              'walltime':'12:00:00',
              'account':"spider2"}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate pbmm2
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {work_path}/align
    mkdir -p {log_path}
    cd {work_path}/align
    pbmm2 align --preset HIFI  -j 32 {inputs["asm"]} {inputs["pacb"]} {spid}_hifi_pbmm2.bam
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
