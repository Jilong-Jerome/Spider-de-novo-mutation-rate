#!/bin/env python3
from gwf import AnonymousTarget
import os, glob

def prepare_depth_stat(work_path,log_path,script_path, vcf, sp):
    inputs = {"vcf":vcf}
    outputs = {
            "log":f"{log_path}/{sp}_depth_summary_per_ind_per_chrom.DONE",
            "depth":f"{work_path}/depth/{sp}_depth_summary.tsv"
            }
    options = {
              'cores':4,
              'memory':'200g',
              'walltime':'72:00:00',
              'account':"spider2"
              }
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate pysam
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {work_path}/depth/{sp}
    cd {work_path}/depth/{sp}
    python {script_path}/depth_summary_per_ind_per_chrom.py {vcf} {sp}_depth_summary.tsv 
    mkdir -p {log_path}
    echo done > {log_path}/{sp}_depth_summary_per_ind_per_chrom.DONE
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def prepare_gene(work_path, gene, target_sp, target_sp_id, gene_blast,gene_blast_sp,blast_sp_gff, target_sp_gff,orthogroups,log_path):
    inputs = {
            "blast":gene_blast,
            "blast_gff":blast_sp_gff,
            "target_gff":target_sp_gff,
            "orthogroups":orthogroups
            }
    outputs = {
            "log":f"{log_path}/{target_sp}_{gene}_prepare.DONE"
            }
    options = {
              'cores':1,
              'memory':'1g',
              'walltime':'1:00:00',
              'account':"spider2"
              }
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate biopython
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {work_path}/prepare/{target_sp}
    cd {work_path}/prepare/{target_sp}
    less {gene_blast}| grep {gene}_ |cut -f1 >  {gene_blast_sp}_{gene}_homologus.txt
    read line < {gene_blast_sp}_{gene}_homologus.txt
    CHROM=$(echo "$line" | awk -F'_' '{{print $1"_"$2"_"$3}}')
    HIC_ID=$(echo "$line" | awk -F'_' '{{print $3}}')
    GENE_ID=$(echo "$line" | awk -F'[_.]' '{{print $4}}')
    echo "Find gene {gene} in {gene_blast_sp}"
    echo "CHROM = $CHROM"
    echo "GENE_ID = $GENE_ID"
    less {blast_sp_gff}|grep $CHROM| grep $GENE_ID > {gene_blast_sp}_{gene}_blast.gff
    echo "Finding orthologs"
    echo "{gene_blast_sp}.HiC_${{HIC_ID}}_file_1_file_1_${{GENE_ID}}"
    less {orthogroups} | grep {gene_blast_sp}.HiC_${{HIC_ID}}_file_1_file_1_${{GENE_ID}} > {target_sp}_{gene}_orthogroups.tsv
    cut -f {target_sp_id} {target_sp}_{gene}_orthogroups.tsv |sed 's/\\n//g'|sed  's/, /\\n/g' > {target_sp}_{gene}.list
    index=1
    while read full_name; do
        echo $full_name
        full_name=$(echo $full_name | tr -d '\\r\\n')
        echo {target_sp_gff}
        cp {target_sp_gff}  {target_sp}_{gene}_temp.gff
        less {target_sp}_{gene}_temp.gff | grep ${{full_name}}> {target_sp}_{gene}_$index.gff
        rm  {target_sp}_{gene}_temp.gff
        index=$((index + 1))
        echo "$index: full_name: $cleaned"
    done < {target_sp}_{gene}.list
    mkdir -p {log_path}
    echo done > {log_path}/{target_sp}_{gene}_prepare.DONE
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def liftoff_gff(work_path, old_fasta, new_fasta, sp, gene ,log_path):
    inputs = {
            "prepare_log":f"{log_path}/{sp}_{gene}_prepare.DONE",
            "old_genome":old_fasta,
            "new_genome":new_fasta
            }
    outputs = {
            "log":f"{log_path}/{sp}_{gene}_lift.DONE",
            "gff":f"{work_path}/lift/{sp}/{gene}/{sp}_{gene}.gff"
            }
    options = {
              'cores':8,
              'memory':'16g',
              'walltime':'12:00:00',
              'account':"spider2"}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate liftoff
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {work_path}/lift/{sp}/{gene}
    cd {work_path}/lift/{sp}/{gene}
    mkdir -p {log_path}
    liftoff {new_fasta} {old_fasta} -g {work_path}/prepare/{sp}/{sp}_{gene}_1.gff -o {sp}_{gene}.gff -p 8 -u unmapped.gff -polish
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def retrieve_vcf_gene(work_path, vcf, sp, gene, log_path):
    inputs = {
            "vcf": vcf,
            "gff":f"{work_path}/lift/{sp}/{gene}/{sp}_{gene}.gff"
            }
    outputs = {
            "log":f"{log_path}/{sp}_{gene}_vcf_retrival.DONE",
            "gene_vcf":f"{work_path}/vcfs/{sp}/{gene}/{sp}_{gene}.vcf"}
    options = {
              'cores':4,
              'memory':'200g',
              'walltime':'2:00:00',
              'account':"spider2"}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bcftools
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {work_path}/vcfs/{sp}/{gene}
    cd {work_path}/vcfs/{sp}/{gene}
    mkdir -p {log_path}
    awk 'BEGIN{{OFS="\t"}} $3=="exon" {{print $1, $4-1, $5}}' {inputs["gff"]} > {sp}_{gene}.bed
    bcftools view -R {sp}_{gene}.bed -Ov -o {sp}_{gene}.vcf {inputs["vcf"]}
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def check_gene_features(work_path,script_path,mean_depth_path,sp,gene,log_path):
    inputs = {
            "retreive_log":f"{log_path}/{sp}_{gene}_vcf_retrival.DONE",
            }
    outputs = {
            "log":f"{log_path}/{sp}_{gene}_feature_summary.DONE"
            }
    options = {
              'cores':1,
              'memory':'4g',
              'walltime':'1:00:00',
              'account':"spider2"}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate pysam
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {work_path}/vcfs/{sp}/{gene}
    cd {work_path}/vcfs/{sp}/{gene}
    CHROM=$(head -n 1 {sp}_{gene}.bed| cut -f 1)
    cp {mean_depth_path}/{sp}_${{CHROM}}.idepth {sp}_{gene}_mean_DP.idepth
    cp {work_path}/lift/{sp}/{gene}/{sp}_{gene}.gff {sp}_{gene}.gff
    python {script_path}/feature_summary.py {work_path}/vcfs/{sp}/{gene}/{sp}_{gene}.vcf {sp}_{gene}_mean_DP.idepth {sp}_{gene}.gff {sp}_{gene}.pdf {sp}_{gene}_summary.tsv
    mkdir -p {log_path}
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)


