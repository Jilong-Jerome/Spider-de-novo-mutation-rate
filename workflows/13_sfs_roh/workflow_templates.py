from gwf import *
LOG_PATH = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/logs/SFS"

def process_vcfs(path,vcf_list,ind_list,outname):
    inputs = [vcf_list,ind_list]
    outputs = [LOG_PATH+"/{outname}_SNPs.DONE".format(outname=outname)]
    options = {
               'cores': 1,
               'memory': '2g',
               'walltime':"2:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bcftools
    echo jobinfo $SLURM_JOBID
    echo "start merging and filtering for SNPs vcfs"
    date
    mkdir -p {path}
    cd {path}
    /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/SFS/process_vcfs.sh {vcfs} {inds} {outname}.vcf
    echo "finish merging vcfs"
    echo done > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,vcfs=vcf_list,inds=ind_list,log=outputs[0],outname=outname)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
