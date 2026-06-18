import sys
from Bio import AlignIO
fasta_aln = sys.argv[1]
outname = sys.argv[2]
out = open(outname,"w")
out.write("site\tgap\tpoly\tremove\n")
alignment = AlignIO.read(open(fasta_aln),"fasta")
codon_dict = {}
for i in range(alignment.get_alignment_length()):
    codon_pos = i % 3
    nuc_list = []
    for record in alignment:
        nuc_list.append(record.seq[i])
    #print(set(nuc_list.remove("-")))
    if "-" in nuc_list:
        site_gap = True
    else:
        site_gap = False
    if ("-" in nuc_list):
        nuc_list.remove('-')
    if len(set(nuc_list)) > 1:
        site_poly = True
    else:
        site_poly = False
    codon_dict[codon_pos]=[i,site_gap,site_poly]
    if codon_pos == 2:
        gap_stage = []
        for key in codon_dict:
            gap_stage.append(codon_dict[key][site_gap])
        if gap_stage.count(True) > 0:
            codon_remove = True
        else:
            codon_remove = False
        for n in range(3):
            out.write("{site}\t{gap}\t{poly}\t{remove}\n".format(site = codon_dict[n][0], gap = codon_dict[n][1], poly = codon_dict[n][2], remove = codon_remove))
