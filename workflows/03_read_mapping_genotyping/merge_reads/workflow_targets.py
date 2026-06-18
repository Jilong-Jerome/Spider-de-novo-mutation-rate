from workflow_templates import *
from gwf import *
def merge_fq(gwf,outname,fq_list,path):
    gwf.target_from_template(
            name = "merge_fq_{task}".format(task=outname),
        template = combine_fq(
            combine_list = sorted(fq_list),
            path = path,
            outname = outname
        )
    )
