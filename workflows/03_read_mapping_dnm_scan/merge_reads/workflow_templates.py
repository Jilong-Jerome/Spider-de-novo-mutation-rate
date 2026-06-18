from gwf import *
LOG_PATH = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/logs/data_prepare"
def combine_fq(combine_list,path,outname):
    """ Combine reads from several Lane to a single file"""
    inputs = combine_list
    catfiles=""
    for i in inputs:
        catfiles=catfiles+str(i)+" "
    outputs = [LOG_PATH+'/merge_{outname}.DONE'.format(outname=outname)]
    options = {
              'cores': 1,
              'memory': "2g",
              'walltime':"08:00:00",
              'account':"spider2"
              }
    spec = """
    echo jobinfo $SLURM_JOBID
    echo "start combining fq"
    date
    mkdir -p {path}
    cd {path}
    cat {catfiles} > {combined}.fq.gz
    echo done > {log}
    echo "finished combination"
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,catfiles=catfiles,combined=outname,log = outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
