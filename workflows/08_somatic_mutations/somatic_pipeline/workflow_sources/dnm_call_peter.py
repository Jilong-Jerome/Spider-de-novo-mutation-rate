import pysam
import sys
from var_call_pacbio import candidate_mutations_first_run, candidate_mutations_second_run, build_fastq_single_read, fastq_dicts

bam_in = sys.argv[1]
fa_in =  sys.argv[2]
contig_name = sys.argv[3]
contig_length = sys.argv[4]
quality_threshold = sys.argv[5]
len_of_interval_around_mut = sys.argv[6]
min_base_qual = sys.argv[7]
min_map_qual = sys.argv[8]
reference = pysam.FastaFile(fa_in)
bamfile = pysam.AlignmentFile(bam_in, "rb")

mutation_reads_fastq_path = sys.argv[9]

reference = pysam.FastaFile(fa_in)
bam_alignment = pysam.AlignmentFile(bam_in, "rb")

mutation_reads_fastq = open(mutation_reads_fastq_path, 'a')
candidate_mutations_first_run(bamfile, reference, contig_name, int(quality_threshold), int(min_base_qual), int(min_map_qual))

print('contig', 
      'start', 
      'end', 
      'ref', 
      'alt', 
      'context', 
      'oriented_ref',
      'oriented_alt', 
      'oriented_context', 
      'mutation_qual', 
      'strand', 
      'reads_with_alt', 
      'coverage', 
      'sm', 
      'sx', 
      'read_name', 
      'avg_read_lengths',
      'position_in_read',
      sep='\t')

mutation_reads = candidate_mutations_second_run(bamfile, reference, contig_name, int(quality_threshold), int(len_of_interval_around_mut), int(min_base_qual), int(min_map_qual))

output_path = '/home/pepo/mutationalscanning/Workspaces/pepo/de_novo_assembly/results/variant_calling/human/ob007/pacbio/diploid/'

for read in mutation_reads:
    position = read[0]
    read_name = read[1]
    for pileupcolumn in bam_alignment.pileup(contig_name, start = position-1,  end = position):
        build_fastq_single_read(pileupcolumn, read_name)

single_read_dict = fastq_dicts()
for read in single_read_dict:
    mutation_reads_fastq.write('\n'.join(single_read_dict[read])+'\n')

