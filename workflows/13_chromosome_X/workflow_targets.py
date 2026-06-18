from workflow_templates import *
from gwf import *

def cal_coverage(gwf,path,indname,sp,ref_sp):
    gwf.target_from_template(
        name = "cov_per_window_{indname}".format(indname = indname),
        template = samtools_depth(path,indname,sp,ref_sp)

        )

def viz_coverage(gwf,path,indname):
    gwf.target_from_template(
        name = "viz_depth_{indname}".format(indname = indname),
        template = viz_depth(path,indname)
        )

