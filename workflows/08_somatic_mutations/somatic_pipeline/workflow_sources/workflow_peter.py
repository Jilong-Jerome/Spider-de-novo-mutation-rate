from gwf import *
gwf = Workflow()

#assembly = 'haploid'
assembly = 'diploid'

#human_id = 'ob2v' 
#human_id = 'ob003h'
human_id = 'ob006'
#human_id = 'la4'
#human_id = 'da4'
#human_id = 'tsaed'
#human_id = 'da1'
#human_id = 'ob007'
#human_id = 'aa_sperm'
inpath = "/home/pepo/mutationalscanning/DerivedData/"
# Use chimp varaiant caller when there is no position
# ############################################################################################################################################################################################
# #                                                                          VARIANT CALLING                                                                                                 #
# ############################################################################################################################################################################################

def var_call(bam_in, fa_in, fai_in, bed_out):
    inputs = [fai_in]
    outputs = [bed_out]
    options = {
              'cores':1,
              'memory':'2g',
              'walltime':'04:00:00', 
              'account':"mutationalscanning"}
    spec =""" 
        date
        echo jobinfo $SLURM_JOBID
        python /home/pepo/mutationalscanning/Workspaces/pepo/de_novo_assembly/scripts/variant_calling/human/dnm_call.py {bam_in} {fa_in} {contig_name} {contig_length} {quality_threshold} {len_of_interval_around_mut} {min_base_qual} {min_map_qual} {mutation_reads_fastq} > {bed_out}.tmp
        mv {bed_out}.tmp {bed_out}
        date
    """.format(bam_in = bam_in, 
               fa_in = fa_in, 
               fai_in = fai_in, 
               contig_name = contig_name, 
               contig_length = contig_dict[contig_name], 
               quality_threshold = quality_threshold, 
               len_of_interval_around_mut = len_of_interval_around_mut, 
               min_base_qual = min_base_qual, 
               min_map_qual = min_map_qual, 
               bed_out = bed_out, 
               mutation_reads_fastq = mutation_reads_fastq)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

bam_in = '/home/pepo/mutationalscanning/DerivedData/bam/HiFi/human/{human_id}/kinetics/{human_id}_kinetics_{assembly}.bam'.format(human_id = human_id, assembly = assembly)
fa_in = '/home/pepo/mutationalscanning/DerivedData/fasta/HiFi/human/{human_id}/{assembly}_assembly/{human_id}_{assembly}.fa'.format(human_id = human_id, assembly = assembly)
fai_in = '/home/pepo/mutationalscanning/DerivedData/fasta/HiFi/human/{human_id}/{assembly}_assembly/{human_id}_{assembly}.fa.fai'.format(human_id = human_id, assembly = assembly)
quality_threshold = 93  #40
len_of_interval_around_mut = 10
min_base_qual = 55 # 40
min_map_qual = 1
contig_dict = {}
with open('{fai_in}'.format(fai_in = fai_in)) as contigs:
    for line in contigs:
        contig_dict[line.split()[0]] = line.split()[1]

file_names = []
for contig_name, length in contig_dict.items():
    if int(length) > 10**6:
        bed_out = '../../../../results/variant_calling/human/{human_id}/pacbio/{assembly}/{contig_name}_pacbio_variants.bed'.format(human_id = human_id, assembly = assembly, contig_name = contig_name)
        mutation_reads_fastq = '../../../../results/variant_calling/human/{human_id}/pacbio/{assembly}/mutation_reads_{contig_name}.fastq'.format(human_id = human_id, assembly = assembly, contig_name = contig_name)
        file_names.append(bed_out)
        gwf.target_from_template(
            name = "{human_id}_{contig_name}_variant_{assembly}".format(human_id = human_id, assembly = assembly, contig_name = contig_name),
            template = var_call(bam_in, fa_in, fai_in, bed_out)
)

# ############################################################################################################################################################################################
# #                                                                        CONCATENATING                                                                                                     #
# ############################################################################################################################################################################################
def concat():
    inputs = file_names
    outputs = [bed_out+'.bed']
    options = {
              'cores':12, 
              'memory':'128g', 
              'walltime':'02:00:00',
              'account':"mutationalscanning"}
    spec =""" 
        echo jobinfo $SLURM_JOBID
        date
        cat ../../../../results/variant_calling/human/{human_id}/pacbio/{assembly}/h*.bed > {bed_out}.bed
        rm ../../../../results/variant_calling/human/{human_id}/pacbio/{assembly}/h*
        cat ../../../../results/variant_calling/human/{human_id}/pacbio/{assembly}/mutation_reads*.fastq > {fastq_out}.fastq
        rm ../../../../results/variant_calling/human/{human_id}/pacbio/{assembly}/mutation_reads*.fastq*
        pbmm2 align /home/pepo/mutationalscanning/DerivedData/fasta/HiFi/human/{human_id}/haploid_assembly/{human_id}_haploid.fa {fastq_out}.fastq /home/pepo/mutationalscanning/Workspaces/pepo/de_novo_assembly/results/variant_calling/human/{human_id}/pacbio/{assembly}/{human_id}_single_read_ref_pos.bam --preset HIFI --sort -j 12 -J 4 -m 8G 
        samtools view /home/pepo/mutationalscanning/Workspaces/pepo/de_novo_assembly/results/variant_calling/human/{human_id}/pacbio/{assembly}/{human_id}_single_read_ref_pos.bam | awk -v OFS='\t' '{{print $1,$2,$3,$4,$5,$6}}' > /home/pepo/mutationalscanning/Workspaces/pepo/de_novo_assembly/results/variant_calling/human/{human_id}/pacbio/{assembly}/{human_id}_read_locations.bed
        date
    """.format(bed_out = bed_out, fastq_out = fastq_out, assembly = assembly, human_id = human_id)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

bed_out = '../../../../results/variant_calling/human/{human_id}/pacbio/{assembly}/{human_id}_true_mutations'.format(human_id = human_id, assembly = assembly)
fastq_out = '../../../../results/variant_calling/human/{human_id}/pacbio/{assembly}/{human_id}_mutation_reads'.format(human_id = human_id, assembly = assembly)

gwf.target_from_template(
    name = "{human_id}_concat".format(human_id = human_id),
    template = concat()
)


