from gwf import*
from workflow_dicts import *
LOG_PATH = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/logs/ncbi"

def trim_ncbi(path,trim,genome,outname):
    inputs = [trim]
    outputs = [LOG_PATH+"/trim_{outname}.DONE".format(outname=outname)]
    options = {
               'cores': 4,
               'memory': '16g',
               'walltime':"2:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate biopython
    echo jobinfo $SLURM_JOBID
    echo "start trimming contamination found from NCBI"
    mkdir -p {path}
    cd {path}
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/ncbi/trim_contamination.py {trim} {genome} {outname}.fa
    echo done > {log}
    """.format(path=path,trim=trim,genome=genome,outname=outname,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def rename_ncbi(path,rename,genome,outname):
    inputs = [LOG_PATH+"/trim_{outname}.DONE".format(outname=outname)]
    outputs = [LOG_PATH+"/rename_{outname}.DONE".format(outname=outname)]
    options = {
               'cores': 4,
               'memory': '16g',
               'walltime':"2:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate biopython
    echo jobinfo $SLURM_JOBID
    echo "start trimming contamination found from NCBI"
    mkdir -p {path}
    cd {path}
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/ncbi/rename_chromosome.py {rename} {genome} {outname}_chromosome.fa {outname}_full.fa
    rm {genome}
    echo done > {log}
    """.format(path=path,rename=rename,genome=genome,outname=outname,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)



gwf = Workflow()
# Trim contamination sequences
for sp in ["DUM","TENT","SARA","BI","MIM","LIN"]:
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/data/genome"
    trim = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/data/genome/{sp}_ncbi.tsv".format(sp=sp)
    genome = sp_raw_genome[sp] 
    outname = "{sp}_ncbi".format(sp=sp)
    gwf.target_from_template(
        name = "trim_{name}".format(name=outname),
        template = trim_ncbi(path,trim,genome,outname)
    )
# Rename chromosomes
for sp in ["DUM","TENT","SARA","BI","MIM","LIN"]:
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/data/genome"
    rename = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/data/genome/{sp}_rename.tsv".format(sp=sp)
    outname = "{sp}_ncbi".format(sp=sp)
    genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/data/genome/{outname}.fa".format(outname=outname)
    gwf.target_from_template(
        name = "rename_{name}".format(name=outname),
        template = rename_ncbi(path,rename,genome,outname)
    )
