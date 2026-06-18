from gwf import *
import os

def genome_composition(work_path,script_path,log_path,sp,genome,bed):
    inputs = [genome,
            bed]
    outputs = [log_path + f"/{sp}_genome_composition_stat.DONE"]
    options = {
               'cores': 1,
               'memory': '8g',
               'walltime':"4:00:00",
               'account':"spider2"
    }
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate biopython
    mkdir -p {work_path}/{sp}
    mkdir -p {log_path}
    cd {work_path}/{sp}
    python {script_path}/genome_context_summary.py {bed} {genome} {sp}_genome_fraction.tsv {sp}_wgs
    echo done > {outputs[0]}
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def split_callable_genome_ind(work_path,script_path,log_path,sp,ind,chrom,vcf):
    inputs = [vcf]
    outputs = [f"{work_path}/{sp}/{ind}_{chrom}_callable_spread.tsv",
            f"{log_path}/{sp}/{ind}_{chrom}_callable_spread.DONE"]
    options = {
               'cores': 1,
               'memory': '1g',
               'walltime':"1:00:00",
               'account':"spider2"
    }
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate biopython
    mkdir -p {work_path}/{sp}/{ind}
    mkdir -p {log_path}/{sp}/{ind}
    cd {work_path}/{sp}
    bash {script_path}/ind_callable_split.sh {vcf} {sp} {ind} {chrom} {ind}_{chrom}_callable_spread.tsv
    echo done > {outputs[1]}
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)


def split_callable_genome_count(work_path,script_path,log_path,sp,ind,chrom,vcf,genome):
    inputs = [vcf,genome]
    outputs = [f"{work_path}/{sp}/{ind}/{ind}_{chrom}_callable_split_count.tsv",
            f"{log_path}/{sp}/{ind}/{ind}_{chrom}_callable_split_count.DONE"]
    options = {
               'cores': 1,
               'memory': '8g',
               'walltime':"1:00:00",
               'account':"spider2"
    }
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate biopython
    echo jobinfo $SLURM_JOBID
    mkdir -p {work_path}/{sp}/{ind}
    mkdir -p {log_path}/{sp}/{ind}
    cd {work_path}/{sp}/{ind}
    python {script_path}/callable_site_split_count.py -g {genome} -s {sp} -i {ind} -c {chrom} -o {ind}_{chrom}_callable_split_count.tsv {vcf}
    echo done > {outputs[1]}
    echo jobinfo $SLURM_JOBID
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
