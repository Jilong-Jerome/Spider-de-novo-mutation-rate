out = open("Table_S1.tsv","w")
for sp in ["BIC","TEN","LIN","MIM","DUM","SAR"]:
    for family in ["family1","family2","family3","family4","family5"]:
        for off in [1,2,3,4,5,6]:
            out.write(f"{sp}\t{sp}_{family}\t{sp}_{family}_F_female\t{sp}_{family}_M_male\t{sp}_{family}_S{off}_offspring\n")

