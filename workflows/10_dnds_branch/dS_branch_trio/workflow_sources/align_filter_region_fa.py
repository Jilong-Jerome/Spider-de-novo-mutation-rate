import sys
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

fasta_aln = sys.argv[1]
regions = sys.argv[2]
outname = sys.argv[3]
alignment = AlignIO.read(open(fasta_aln),"fasta")
region_list = []
regions_file = open(regions)
for region in regions_file:
    region_info = region.split("\t")
    if region_info[0] != "region_n":
        region_start = int(region_info[1])
        region_end = int(region_info[2])
        region_list.append((region_start,region_end))
alignment_list = []
for region in region_list:
    records = []
    for record in alignment:
        seq = record.seq[region[0]:region[1]+1]
        seq_id = record.id
        records.append(SeqRecord(Seq(seq),id = seq_id))
    multi_align = MultipleSeqAlignment(records)
    alignment_list.append(multi_align)

AlignIO.write(alignment_list,outname,"phylip")
