from workflow_dicts import *
from workflow_templates import *
from gwf import *
import os

def run_callable_sites_per_trio(gwf,sp,female,male,offspring):
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/germline/{sp}/callable/{offspring}/chrom_splits".format(sp=sp,offspring=offspring)
    dp_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/germline/{sp}/filters".format(sp=sp)
    log_list = []
    for chrom in sp_chrom_dict[sp]:
        log_list.append("/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/logs/align_step1{sp}_{chrom}_depth.DONE".format(sp=sp,chrom=chrom))
    outfile_name = "{sp}_per_ind_per_chrom_DP.tsv".format(sp = sp)
    if not os.path.exists(outfile_name):
        with open(outfile_name,'w') as outfile:
            outfile.write("sample_id\tN_sites\tmean_depth\tchromosome\n")
            for chrom_id in sp_chrom_dict[sp]:
                chrom_dp_file = dp_path+"/{sp}_{chrom}.idepth".format(sp=sp,chrom=chrom_id)
                with open(chrom_dp_file) as infile:
                    next(infile)
                    for line in infile:
                        outfile.write("{line_info}\t{chrom}\n".format(line_info=line.strip("\n"),chrom = chrom_id))
    for chrom in sp_chrom_dict[sp]:
        input_vcf = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/germline/{sp}/chrom_splits/{sp}_{chrom}.vcf.gz".format(sp=sp,chrom=chrom)
        gwf.target_from_template(
                name = "{child}_{chrom}_callable_sites".format(chrom=chrom,child = offspring),
                template = callable_sites(path,sp,input_vcf,female,male,offspring,dp_path,log_list,chrom)
        )

def run_dnm_scan(gwf,sp,female,male,offspring):
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/germline/{sp}/DNMs/{offspring}/chrom_splits".format(sp=sp,offspring=offspring)
    for chrom in sp_chrom_dict[sp]:
        vcf = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/germline/{sp}/callable/{offspring}/chrom_splits/{offspring}_{chrom}_callable_sites.vcf".format(sp=sp,offspring=offspring,chrom=chrom)
        gwf.target_from_template(
            name = "{child}_{chrom}_DNMs".format(child=offspring,chrom=chrom),
            template = dnm_scan(path,sp,vcf,female,male,offspring,chrom)
        )

def run_bam_check(gwf,sp,female,male,offspring):
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/germline/{sp}/DNMs/{offspring}/chrom_splits".format(sp=sp,offspring=offspring)
    if sp == "SAR":
        ref_sp = "SARA"
    elif sp == "BIC":
        ref_sp = "BI"
    elif sp == "TEN":
        ref_sp = "TENT"
    else:
        ref_sp = sp
    ref = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/data/genome/{sp}_ncbi_chromosome.fa".format(sp=ref_sp)
    for chrom in sp_chrom_dict[sp]:
        vcf = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/germline/{sp}/DNMs/{offspring}/chrom_splits/{offspring}_{chrom}_putative_DNMs_distance_filtered_AB.vcf".format(sp=sp,offspring=offspring,chrom=chrom)
        gwf.target_from_template(
            name = "{child}_{chrom}_DNMs_bam_check".format(child=offspring,chrom=chrom),
            template = dnm_bcftools_check(path,sp,ref,chrom,vcf,female,male,offspring)
         )


def run_igv_check(gwf,sp):
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/germline/{sp}/IGV".format(sp=sp)
    gwf.target_from_template(
        name = "{sp}_IGV_check".format(sp=sp),
        template = create_IGV_check(sp,path)
    )
