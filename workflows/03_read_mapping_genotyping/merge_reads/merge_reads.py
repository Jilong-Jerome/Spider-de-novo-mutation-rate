import re

def generate_merge_dict(species_file,species):
    reads = open(species_file,"r")
    reads_dict={}
    for read in reads:
        read_pair = read.strip("\n").split("/")[-1].split(".")[0][-1]
        family = read.strip("\n").split("/")[-3]
        ind_full = read.strip("\n").split("/")[-2]
        ind = re.split("_|-|\.",ind_full)[-1]
        if (ind == "F" or ind == "female" or ind == "Fema"):
            suffix = "female"
            ind = "F"
        elif (ind == "M" or ind == "male" or ind == "Ma" or ind == "Male"):
            suffix = "male"
            ind = "M"
        else:
            suffix = "offspring"
            if "S" not in ind:
                ind = "S"+ind
        newname = species+"_"+family+"_"+ind+"_"+suffix
        if newname in reads_dict:
            if read_pair == "1":
                reads_dict[newname][0].append(read.strip("\n"))
            if read_pair == "2":
                reads_dict[newname][1].append(read.strip("\n"))
        else:
            reads_dict[newname]=[[],[]]
            if read_pair == "1":
                reads_dict[newname][0].append(read.strip("\n"))
            if read_pair == "2":
                reads_dict[newname][1].append(read.strip("\n"))
    return reads_dict

