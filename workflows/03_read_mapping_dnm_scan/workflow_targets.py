from workflow_dicts import *
from workflow_templates import *
from gwf import *
from split_segments import *
import os
def index_genome(gwf,path,genome,outname):
    gwf.target_from_template(
            name = "index_{task}".format(task=outname),
            template = bwa2_index(
            ref_genome = genome,
            path = path,
            basename = outname
            )
    )

def align_ind(gwf,indname):
    sp = indname.split("_")[0]
    if sp == "TEN":
        ref_sp = "TENT"
    elif sp == "SAR":
        ref_sp =  "SARA"
    elif sp  == "BIC":
        ref_sp = "BI"
    else:
        ref_sp = sp
    read1 = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/data/raw_reads/{sp}/{indname}_R1.fq.gz".format(sp=sp,indname=indname)
    read2 = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/data/raw_reads/{sp}/{indname}_R2.fq.gz".format(sp=sp,indname=indname)
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/align/bams/{sp}".format(sp=sp)
    gwf.target_from_template(
        name = "align_bwa2_{}".format(indname),
        template = bwa2_align(
            path = path,
            sp=ref_sp,
            indname = indname,
            read1 = read1,
            read2 = read2
            )
        )
    gwf.target_from_template(
        name = "picard_addRG_{}".format(indname),
        template = picard_addRG(path,indname)
        )
    gwf.target_from_template(
        name = "samtools_fixmate_{}".format(indname),
        template = samtool_fixmate(path,indname)
        )
    gwf.target_from_template(
        name = "samtools_sort_{}".format(indname),
        template = samtool_sort(path,indname)
        )
    gwf.target_from_template(
        name = "samtools_markdup_{}".format(indname),
        template = samtool_markdup(path,indname)
        )
    gwf.target_from_template(
        name = "samtools_mapped_MQ_{}".format(indname),
        template = samtool_mapped(path,indname)
    )
    gwf.target_from_template(
        name = "picard_bai_{}".format(indname),
        template = picard_bai(path,indname)
    )

def gatk_gvcf_ind(gwf,indname):
    sp = indname.split("_")[0]
    chrom_list = sp_chrom_dict[sp]
    if sp == "SAR":
        ref_sp = "SARA"
    elif sp == "TEN":
        ref_sp = "TENT"
    elif sp == "BIC":
        ref_sp = "BI"
    else:
        ref_sp = sp
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/gatk/{sp}/{indname}/chroms".format(sp=sp,indname=indname)
    assembly = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/data/genome/{sp}_ncbi_chromosome.fa".format(sp=ref_sp)
    bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/align/bams/{sp}/{indname}_final.bam".format(sp=sp,indname=indname)
    log_list = []
    chrom_gvcf_list = []
    chrom_remap_list = []
    for chrom in chrom_list:
        gwf.target_from_template(
            name = "gatk_haplocall_{indname}_{chrom}".format(indname=indname,chrom=chrom),
            template = gatk_haplotype_call(path,bam,assembly,indname,chrom)
            )
        log_list.append(LOG_PATH+'/{indname}_{chrom}_gvcf.DONE'.format(indname=indname,chrom=chrom))
        chrom_gvcf_list.append(path+"/{indname}_{chrom}.g.vcf".format(indname=indname,chrom=chrom))
        chrom_remap_list.append(path+"/{indname}_{chrom}_remap.bam".format(indname=indname,chrom=chrom))
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/gatk/{sp}/{indname}".format(sp=sp,indname=indname)
    merge_list = chrom_gvcf_list
    gwf.target_from_template(
        name = "gatk_gvcf_{indname}".format(indname=indname),
        template = merge_gvcf_by_chrom(log_list,merge_list,path,indname)
        )
    merge_list = chrom_remap_list
    gwf.target_from_template(
        name = "gatk_remap_{indname}".format(indname=indname),
        template = merge_bam_by_chrom(log_list,merge_list,path,indname)
        )

def write_map(sp):
    log_list = []
    sample_map = open("./sample_maps/{sp}_sample_map.tsv".format(sp=sp),"w")
    if sp in ["LIN","MIM","AFR"]:
        fam_list = ["family1","family2","family3","family4"]
    elif sp in ["BIC"]:
        fam_list = ["family1","family2","family3","family4","family5","family6"]
    else:
        fam_list = ["family1","family2","family3","family4","family5"]
    for fam in fam_list:
            ind = sp+"_"+fam+"_"+"M"+"_"+"male"
            ind_vcf = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/gatk/{sp}/{ind}/{ind}.g.vcf".format(sp=sp,ind=ind)
            ind_log = LOG_PATH + '/{indname}_gvcf.DONE'.format(indname=ind)
            log_list.append(ind_log)
            sample_map.write("{ind}\t{vcf}\n".format(ind=ind,vcf=ind_vcf))
            ind = sp+"_"+fam+"_"+"F"+"_"+"female"
            ind_vcf = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/gatk/{sp}/{ind}/{ind}.g.vcf".format(sp=sp,ind=ind)
            ind_log = LOG_PATH + '/{indname}_gvcf.DONE'.format(indname=ind)
            log_list.append(ind_log)
            sample_map.write("{ind}\t{vcf}\n".format(ind=ind,vcf=ind_vcf))
            if sp == "MIM" and fam == "family4":
                child_list = ["S1","S2","S3","S4"]
            elif sp == "DUM" and fam == "family3":
                child_list = ["S1","S2","S4","S5","S6"]
            else:
                child_list = ["S1","S2","S3","S4","S5","S6"]
            for child in child_list:
                ind = sp+"_"+fam+"_"+child+"_"+"offspring"
                ind_vcf = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/gatk/{sp}/{ind}/{ind}.g.vcf".format(sp=sp,ind=ind)
                ind_log = LOG_PATH + '/{indname}_gvcf.DONE'.format(indname=ind)
                log_list.append(ind_log)
                sample_map.write("{ind}\t{vcf}\n".format(ind=ind,vcf=ind_vcf))
    return log_list

def write_gatk_logs_merges(sp):
    log_list = []
    merge_list = []
    intervals = open("./intervals/{sp}_intervals.tsv".format(sp=sp))
    for line in intervals:
        db_id = line.strip("\n").split("\t")[1]
        log = LOG_PATH + "/gatk_geno/GATK4_genotype_segment_{db_id}.DONE".format(db_id = db_id)
        merge = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/gatk/genotype/{sp}/splits/{db_id}_GATK_filtered.vcf.gz".format(sp=sp,db_id=db_id)
        log_list.append(log)
        merge_list.append(merge)
    return log_list,merge_list

def gatk_DB_build(gwf,sp):
    # Example usage parameters
    if sp == "TEN":
        genome_sp = "TENT"
    elif sp == "SAR":
        genome_sp = "SARA"
    elif sp == "BIC":
        genome_sp = "BI"
    else:
        genome_sp = sp
    fasta_index_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/data/genome/{sp}_ncbi_chromosome.fa.fai".format(sp=genome_sp)
    ref = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/data/genome/{sp}_ncbi_chromosome.fa".format(sp=genome_sp)
    segment_size = 5000000  # 5Mb
    output_path = "./intervals/{sp}_intervals.tsv".format(sp=sp)
    generate_intervals_and_save(fasta_index_path, segment_size, sp, output_path)
    log_list = write_map(sp)
    vcf_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/gatk/genotype/{sp}/splits".format(sp=sp)
    batch_size = 50
    sample_map = "./sample_maps/{sp}_sample_map.tsv".format(sp=sp)
    interval_file = open( "./intervals/{sp}_intervals.tsv".format(sp=sp))
    for interval in interval_file:
        infos = interval.strip("\n").split("\t")
        interval_range = infos[0]
        db_id = infos[1]
        gwf.target_from_template(
        name = "gatk_DB_genotype_segment_{db_id}".format(db_id=db_id),
        template = gatk_consolidate(log_list,ref,vcf_path,sp,interval_range,db_id,batch_size,sample_map)
        )

def gatk_finalise(gwf,sp):
    log_list,merge_list = write_gatk_logs_merges(sp)
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/gatk/genotype/{sp}".format(sp=sp)
    outname = "{sp}_GATK".format(sp=sp)
    gwf.target_from_template(
        name = "{sp}_GATK_finalise".format(sp=sp),
        template = merge_vcfs(log_list,merge_list,path,outname)
    )

def mean_DP_per_chrom_per_ind(gwf,sp):
    vcf = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/gatk/genotype/{sp}/{sp}_GATK.vcf.gz".format(sp=sp)
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/germline/{sp}/filters".format(sp=sp)
    for chrom in sp_chrom_dict[sp]:
        gwf.target_from_template(
            name = "{sp}_{chrom}_DP_per_ind".format(sp=sp,chrom=chrom),
            template = vcf_depths(sp,vcf,chrom,path)
        )
def run_chrom_split(gwf,sp):
    out_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/germline/{sp}/chrom_splits".format(sp=sp)
    vcf = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/gatk/genotype/{sp}/{sp}_GATK.vcf.gz".format(sp=sp)
    gwf.target_from_template(
        name = "{sp}_GATK_chrom_split".format(sp=sp),
        template =sp_chrom_vcf_split(vcf,out_path,sp)
    )

def run_chrom_vcf_gzip(gwf,sp):
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/steps/germline/{sp}/chrom_splits".format(sp=sp)
    for chrom in sp_chrom_dict[sp]:
        vcf = "{sp}_{chrom}.vcf".format(sp=sp,chrom=chrom)
        gwf.target_from_template(
            name = "{sp}_{chrom}_vcf_gz".format(sp=sp,chrom=chrom),
            template = gzip_file(sp,chrom,path,vcf)
        )

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
