from gwf import *

def classify_DNM(work_path,script_path,log_path,sp,minDP,genome,tsv):
    inputs = [tsv,genome]
    outputs = [f"{work_path}/{sp}_minDP_{minDP}_automate_classified.tsv",
            f"{log_path}/{sp}_minDP_{minDP}_DNM_classify.DONE"]
    options = {
               'cores': 1,
               'memory': '4g',
               'walltime':"1:00:00",
               'account':"spider2"
    }
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate biopython
    echo jobinfo $SLURM_JOBID
    mkdir -p {work_path}
    mkdir -p {log_path}
    cd {work_path}
    python {script_path}/mut_class_check.py {genome} {tsv} {outputs[0]}
    echo done > {outputs[1]}
    echo jobinfo $SLURM_JOBID
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
