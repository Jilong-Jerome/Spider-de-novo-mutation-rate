from merge_reads import *
from workflow_targets import *
from gwf import * 

gwf = Workflow()
for species in ["LIN","MIM","DUM","TEN","SAR","BIC","AFR"]:
    species_file = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/script/align/lists/{sp}_reads.txt".format(sp=species)
    sp_dict = generate_merge_dict(species_file,species)
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/data/raw_reads/{species}".format(species=species) 
    for ind in sp_dict:
        R1_combine = sp_dict[ind][0]
        R2_combine = sp_dict[ind][1]
        merge_fq(gwf,"{ind}_R1".format(ind=ind),R1_combine,path)
        merge_fq(gwf,"{ind}_R2".format(ind=ind),R2_combine,path)
