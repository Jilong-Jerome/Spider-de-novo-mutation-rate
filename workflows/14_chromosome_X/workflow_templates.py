from gwf import *
LOG_PATH = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/logs/chromX"
def samtools_depth(path,indname,sp,ref_sp):
    """ Template for checking sex chromosomes of a single individual bam"""
    inputs = ["/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/align/bams/{sp}/{indname}_final.bam".format(sp=sp,indname=indname)]
    outputs = [LOG_PATH+'/{indname}_chromX.DONE'.format(indname=indname)]
    options = {
               'cores': 1,
               'memory': '96g',
               'walltime':"12:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bam_check 
    echo jobinfo $SLURM_JOBID
    echo "start finding X chromosomes"
    date
    mkdir -p {path}
    cd {path}
    samtools depth -X /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/align/bams/{sp}/{indname}_final.bam /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/align/bams/{sp}/{indname}_final.bam.bai -o {indname}.depth
    echo done > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(sp=sp,indname=indname,path=path,ref_sp = ref_sp,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def viz_depth(path,indname):
    inputs = [LOG_PATH+'/{indname}_chromX.DONE'.format(indname=indname)]
    outputs = [LOG_PATH+'/{indname}_chromX_viz.DONE'.format(indname=indname)]
    options = {
               'cores': 1,
               'memory': '96g',
               'walltime':"12:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate pandas
    echo jobinfo $SLURM_JOBID
    echo "start finding X chromosomes"
    date
    mkdir -p {path}
    cd {path}
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/chromX/stat_per_chrom.py {indname}.depth {indname}.dp.stat
    echo done > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(indname=indname,path=path,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
