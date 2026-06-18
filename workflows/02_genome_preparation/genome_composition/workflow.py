from gwf import *
from workflow_templates import *
from workflow_targets import *

gwf = Workflow()

genomes = open("genome_meta.tsv")
data_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/data/genome"
work_path = "/faststorage/project/spider2/social_spiders_2020/people/jilong/chapter3_dnm/steps/mut_spectrum/genome_stat"
log_path = "/faststorage/project/spider2/social_spiders_2020/people/jilong/chapter3_dnm/steps/mut_spectrum/logs"
script_path = "/faststorage/project/spider2/social_spiders_2020/people/jilong/chapter3_dnm/script/genome_composition"
for genome in genomes:
    infos = genome.strip("\n").split("\t")
    sp = infos[0]
    bed = f"{data_path}/{infos[1]}"
    genome = f"{data_path}/{infos[2]}"
    run_asm_stat(gwf, work_path, script_path, log_path, sp, genome, bed)


vcfs = open("callable_vcf.txt")
work_path = "/faststorage/project/spider2/social_spiders_2020/people/jilong/chapter3_dnm/steps/mut_spectrum/callable_stat"
log_path = "/faststorage/project/spider2/social_spiders_2020/people/jilong/chapter3_dnm/steps/mut_spectrum/logs"
script_path = "/faststorage/project/spider2/social_spiders_2020/people/jilong/chapter3_dnm/script/genome_composition"
for vcf_file in vcfs:
    vcf = vcf_file.strip("\n")
    sp = vcf.split("/")[-1].split("_")[0]
    infos = vcf.split("/")[-1].split("_")
    ind = f"{infos[0]}_{infos[1]}_{infos[2]}_{infos[3]}"
    chrom = f"{infos[4]}_{infos[5]}"
    run_callable_stat(gwf, work_path, script_path, log_path, sp, ind, chrom, vcf)

vcfs = open("callable_vcf.txt")
work_path = "/faststorage/project/spider2/social_spiders_2020/people/jilong/chapter3_dnm/steps/mut_spectrum/callable_stat"
log_path = "/faststorage/project/spider2/social_spiders_2020/people/jilong/chapter3_dnm/steps/mut_spectrum/logs"
script_path = "/faststorage/project/spider2/social_spiders_2020/people/jilong/chapter3_dnm/script/genome_composition"
for vcf_file in vcfs:
    vcf = vcf_file.strip("\n")
    sp = vcf.split("/")[-1].split("_")[0]
    infos = vcf.split("/")[-1].split("_")
    ind = f"{infos[0]}_{infos[1]}_{infos[2]}_{infos[3]}"
    chrom = f"{infos[4]}_{infos[5]}"
    if sp == "BIC":
        sp_ref = "BI"
    elif sp == "TEN":
        sp_ref = "TENT"
    elif sp == "SAR":
        sp_ref = "SARA"
    else:
        sp_ref = sp
    genome = f"/faststorage/project/spider2/social_spiders_2020/people/jilong/chapter3_dnm/data/genome/{sp_ref}_ncbi_chromosome.fa"
    run_callable_split_count(gwf,work_path,script_path,log_path,sp,ind,chrom,vcf,genome)
