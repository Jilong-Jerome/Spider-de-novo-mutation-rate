#!/usr/bin/env python3
"""
workflow_templates.py - GWF job templates for dS branch trio analysis

Pipeline steps:
  01 - Reference indexing (bwa-mem2 + samtools faidx)
  02 - Read mapping (bwa-mem2 mem)
  03 - BAM processing (addRG, fixmate, sort, markdup, filter, index)
  04 - Variant calling (bcftools mpileup + call)
  05 - Consensus building with configurable heterozygous-site modes
  06 - CDS extraction per gene (custom Python script)
  07 - Prepare gene lists (autosomal vs X-linked)
  08 - Bootstrap PAML analysis (sample + concat + filter + codeml + parse)
"""

from gwf import AnonymousTarget


# ===========================================================
#  Step 01 - Reference Indexing
# ===========================================================

def bwa2_index_template(ref_genome, ref_species, output_dir, log_dir):
    """Index reference genome with bwa-mem2 and samtools faidx."""
    inputs = [ref_genome]
    outputs = [f"{log_dir}/01_bwa2_index.DONE"]
    options = {
        'cores': 12,
        'memory': '128g',
        'walltime': '12:00:00',
        'account': 'spider2'
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bwa2
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {output_dir} {log_dir}
    cd {output_dir}
    bwa-mem2 index -p {ref_species} {ref_genome}
    conda activate samtools117
    samtools faidx {ref_genome}
    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        ref_genome=ref_genome,
        ref_species=ref_species,
        output_dir=output_dir,
        log_dir=log_dir,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ===========================================================
#  Step 02 - Read Mapping
# ===========================================================

def bwa2_align_template(ref_index_prefix, read1, read2, species, output_dir, log_dir):
    """Align paired-end reads to reference with bwa-mem2."""
    inputs = [
        f"{log_dir}/01_bwa2_index.DONE",
        read1,
        read2
    ]
    outputs = [f"{log_dir}/02_bwa2_align_{species}.DONE"]
    options = {
        'cores': 20,
        'memory': '64g',
        'walltime': '24:00:00',
        'account': 'spider2'
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bwa2
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {output_dir} {log_dir}
    cd {output_dir}
    bwa-mem2 mem -t 20 {ref_index_prefix} {read1} {read2} | samtools view -Sb > {species}_s0.bam
    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        ref_index_prefix=ref_index_prefix,
        read1=read1,
        read2=read2,
        species=species,
        output_dir=output_dir,
        log_dir=log_dir,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ===========================================================
#  Step 03 - BAM Processing (sequential chain per species)
# ===========================================================

def picard_addRG_template(species, bam_dir, log_dir):
    """Add read groups with picard AddOrReplaceReadGroups."""
    inputs = [f"{log_dir}/02_bwa2_align_{species}.DONE"]
    outputs = [f"{log_dir}/03a_addRG_{species}.DONE"]
    options = {
        'cores': 8,
        'memory': '16g',
        'walltime': '12:00:00',
        'account': 'spider2'
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate python_phylo
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {log_dir}
    cd {bam_dir}
    java -jar /home/jilong/software/picard.jar AddOrReplaceReadGroups \
        --I {species}_s0.bam --O {species}_s1.bam \
        -RGID {species} -RGPU unknown -RGSM {species} -RGPL illumina -RGLB lib0
    rm {species}_s0.bam
    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        species=species,
        bam_dir=bam_dir,
        log_dir=log_dir,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def samtools_fixmate_template(species, bam_dir, log_dir):
    """Fix mate information and remove secondary/unmapped reads."""
    inputs = [f"{log_dir}/03a_addRG_{species}.DONE"]
    outputs = [f"{log_dir}/03b_fixmate_{species}.DONE"]
    options = {
        'cores': 16,
        'memory': '32g',
        'walltime': '12:00:00',
        'account': 'spider2'
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate samtools117
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {log_dir}
    cd {bam_dir}
    samtools fixmate -rm {species}_s1.bam {species}_s2.bam -@ 16
    rm {species}_s1.bam
    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        species=species,
        bam_dir=bam_dir,
        log_dir=log_dir,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def samtools_sort_template(species, bam_dir, log_dir):
    """Sort BAM file by coordinate."""
    inputs = [f"{log_dir}/03b_fixmate_{species}.DONE"]
    outputs = [f"{log_dir}/03c_sort_{species}.DONE"]
    options = {
        'cores': 16,
        'memory': '32g',
        'walltime': '12:00:00',
        'account': 'spider2'
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate samtools117
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {log_dir}
    cd {bam_dir}
    samtools sort {species}_s2.bam -o {species}_s3.bam -@ 16
    rm {species}_s2.bam
    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        species=species,
        bam_dir=bam_dir,
        log_dir=log_dir,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def samtools_markdup_template(species, bam_dir, log_dir):
    """Mark and remove duplicate reads."""
    inputs = [f"{log_dir}/03c_sort_{species}.DONE"]
    outputs = [f"{log_dir}/03d_markdup_{species}.DONE"]
    options = {
        'cores': 16,
        'memory': '64g',
        'walltime': '12:00:00',
        'account': 'spider2'
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate samtools117
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {log_dir}
    cd {bam_dir}
    samtools markdup -r -f {species}.stat -s {species}_s3.bam {species}_s4.bam -@ 16
    rm {species}_s3.bam
    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        species=species,
        bam_dir=bam_dir,
        log_dir=log_dir,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def samtools_filter_template(species, bam_dir, log_dir, min_mq):
    """Filter BAM: mapped, properly paired, min mapping quality."""
    inputs = [f"{log_dir}/03d_markdup_{species}.DONE"]
    outputs = [f"{log_dir}/03e_filter_{species}.DONE"]
    options = {
        'cores': 16,
        'memory': '64g',
        'walltime': '12:00:00',
        'account': 'spider2'
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate samtools117
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {log_dir}
    cd {bam_dir}
    samtools view -@ 16 -bq {min_mq} -f 0x2 -F 0x4 {species}_s4.bam > {species}_final.bam
    rm {species}_s4.bam
    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        species=species,
        bam_dir=bam_dir,
        min_mq=min_mq,
        log_dir=log_dir,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def samtools_index_template(species, bam_dir, log_dir):
    """Index final BAM file."""
    inputs = [f"{log_dir}/03e_filter_{species}.DONE"]
    outputs = [f"{log_dir}/03f_index_{species}.DONE"]
    options = {
        'cores': 8,
        'memory': '16g',
        'walltime': '04:00:00',
        'account': 'spider2'
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate samtools117
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {log_dir}
    cd {bam_dir}
    samtools index {species}_final.bam {species}_final.bam.bai
    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        species=species,
        bam_dir=bam_dir,
        log_dir=log_dir,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ===========================================================
#  Step 04 - Variant Calling (bcftools mpileup + call)
# ===========================================================

def bcftools_call_template(ref_genome, species, bam_dir, output_dir, log_dir, min_bq):
    """Call variants with bcftools mpileup + call (consensus calling)."""
    inputs = [f"{log_dir}/03f_index_{species}.DONE"]
    outputs = [f"{log_dir}/04_bcftools_call_{species}.DONE"]
    options = {
        'cores': 8,
        'memory': '32g',
        'walltime': '24:00:00',
        'account': 'spider2'
    }
    bam = f"{bam_dir}/{species}_final.bam"
    vcf_out = f"{output_dir}/{species}.vcf.gz"
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bcftools
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {output_dir} {log_dir}
    bcftools mpileup -f {ref_genome} -Q {min_bq} --threads 8 {bam} | \
        bcftools call -c --threads 8 -Oz -o {vcf_out}
    bcftools index {vcf_out}
    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        ref_genome=ref_genome,
        min_bq=min_bq,
        bam=bam,
        output_dir=output_dir,
        vcf_out=vcf_out,
        log_dir=log_dir,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ===========================================================
#  Step 05a - Shared Callable-Depth Mask
# ===========================================================

def callable_depth_template(species, bam_dir, output_dir, log_dir, script_path,
                            min_dp, coverage_median_min_factor,
                            coverage_median_max_factor):
    """Build per-species callable-depth stats and merged uncallable BED."""
    bam = f"{bam_dir}/{species}_final.bam"
    inputs = [f"{log_dir}/03f_index_{species}.DONE"]
    genomecov_bedgraph = f"{output_dir}/{species}_genomecov.bedgraph"
    coverage_stats = f"{output_dir}/{species}_coverage_stats.tsv"
    uncallable_bed = f"{output_dir}/{species}_uncallable_depth.bed"
    depth_distribution = f"{output_dir}/{species}_depth_distribution.tsv"
    outputs = [
        f"{log_dir}/05a_callable_depth_{species}.DONE",
        coverage_stats,
        uncallable_bed,
        depth_distribution
    ]
    options = {
        'cores': 1,
        'memory': '16g',
        'walltime': '12:00:00',
        'account': 'spider2'
    }
    spec = """
    set -eo pipefail
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {output_dir} {log_dir}

    conda activate bedtools
    echo "STEP genomecov START: $(date)"
    bedtools genomecov -ibam {bam} -bga > {genomecov_bedgraph}
    echo "STEP genomecov FINISH: $(date)"
    conda activate python_phylo
    echo "STEP callable_depth_regions START: $(date)"
    python3 {script_path}/callable_depth_regions.py \
        --genomecov {genomecov_bedgraph} \
        --stats {coverage_stats} \
        --uncallable-bed {uncallable_bed} \
        --depth-distribution {depth_distribution} \
        --min-depth {min_dp} \
        --min-factor {coverage_median_min_factor} \
        --max-factor {coverage_median_max_factor}
    echo "STEP callable_depth_regions FINISH: $(date)"

    CALLABLE_MIN_DP=$(awk -F '\\t' '$1=="callable_min_depth" {{print $2}}' {coverage_stats})
    CALLABLE_MAX_DP=$(awk -F '\\t' '$1=="callable_max_depth" {{print $2}}' {coverage_stats})
    echo "Callable depth range: $CALLABLE_MIN_DP-$CALLABLE_MAX_DP"
    echo "Uncallable-depth regions: $(wc -l < {uncallable_bed})"
    rm {genomecov_bedgraph}
    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        bam=bam,
        output_dir=output_dir,
        log_dir=log_dir,
        genomecov_bedgraph=genomecov_bedgraph,
        coverage_stats=coverage_stats,
        uncallable_bed=uncallable_bed,
        depth_distribution=depth_distribution,
        script_path=script_path,
        min_dp=min_dp,
        coverage_median_min_factor=coverage_median_min_factor,
        coverage_median_max_factor=coverage_median_max_factor,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ===========================================================
#  Step 05a (diagnostic) - Coverage distribution plots per trio
# ===========================================================

def coverage_distribution_plot_template(species_list, callable_depth_dir,
                                        log_dir, script_path, trio_name):
    """Reference-only per-trio coverage-distribution TSV summary and plots."""
    inputs = [
        f"{log_dir}/05a_callable_depth_{species}.DONE"
        for species in species_list
    ]
    out_dir = f"{callable_depth_dir}/coverage_distribution"
    summary_tsv = f"{out_dir}/coverage_distribution_summary.tsv"
    outputs = [
        f"{log_dir}/05a_coverage_distribution.DONE",
        summary_tsv,
        f"{out_dir}/coverage_distribution.png",
        f"{out_dir}/coverage_distribution.pdf"
    ]
    options = {
        'cores': 1,
        'memory': '4g',
        'walltime': '00:30:00',
        'account': 'spider2'
    }
    distributions = ','.join(
        f"{sp}={callable_depth_dir}/{sp}_depth_distribution.tsv"
        for sp in species_list
    )
    stats = ','.join(
        f"{sp}={callable_depth_dir}/{sp}_coverage_stats.tsv"
        for sp in species_list
    )
    spec = """
    set -eo pipefail
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate python_phylo
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {out_dir} {log_dir}

    python3 {script_path}/plot_coverage_distribution.py \
        --species-distributions {distributions} \
        --species-stats {stats} \
        --output-dir {out_dir} \
        --trio-name {trio_name}

    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        out_dir=out_dir,
        log_dir=log_dir,
        script_path=script_path,
        distributions=distributions,
        stats=stats,
        trio_name=trio_name,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ===========================================================
#  Step 05b - Shared SNP Filtering
# ===========================================================

def consensus_filter_snps_template(species, vcf_dir, output_dir,
                                   callable_depth_dir, shared_log_dir):
    """Filter SNPs by shared callable-depth bounds once per species."""
    inputs = [
        f"{shared_log_dir}/04_bcftools_call_{species}.DONE",
        f"{shared_log_dir}/05a_callable_depth_{species}.DONE",
    ]
    filtered_vcf = f"{output_dir}/{species}_snps_filtered.vcf.gz"
    filtered_vcf_plain = f"{output_dir}/{species}_snps_filtered.vcf"
    outputs = [
        f"{shared_log_dir}/05b_filter_snps_{species}.DONE",
        filtered_vcf,
        filtered_vcf_plain,
    ]
    options = {
        'cores': 1,
        'memory': '16g',
        'walltime': '08:00:00',
        'account': 'spider2'
    }
    vcf_in = f"{vcf_dir}/{species}.vcf.gz"
    coverage_stats = f"{callable_depth_dir}/{species}_coverage_stats.tsv"
    spec = """
    set -eo pipefail
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bcftools
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {output_dir} {shared_log_dir}

    CALLABLE_MIN_DP=$(awk -F '\\t' '$1=="callable_min_depth" {{print $2}}' {coverage_stats})
    CALLABLE_MAX_DP=$(awk -F '\\t' '$1=="callable_max_depth" {{print $2}}' {coverage_stats})
    echo "Callable depth range: $CALLABLE_MIN_DP-$CALLABLE_MAX_DP"

    echo "STEP filter_vcf START: $(date)"
    bcftools view -i 'TYPE="snp" && QUAL>=20 && DP>='$CALLABLE_MIN_DP' && DP<='$CALLABLE_MAX_DP'' {vcf_in} -Oz -o {filtered_vcf}
    echo "STEP filter_vcf FINISH: $(date)"
    echo "STEP index_filtered_vcf START: $(date)"
    bcftools index {filtered_vcf}
    echo "STEP index_filtered_vcf FINISH: $(date)"
    echo "STEP write_plain_filtered_vcf START: $(date)"
    bcftools view -Ov -o {filtered_vcf_plain} {filtered_vcf}
    echo "STEP write_plain_filtered_vcf FINISH: $(date)"
    grep -q '^#CHROM' {filtered_vcf_plain}
    echo "Filtered SNPs: $(bcftools view -H {filtered_vcf} | wc -l)"

    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        vcf_in=vcf_in,
        output_dir=output_dir,
        shared_log_dir=shared_log_dir,
        coverage_stats=coverage_stats,
        filtered_vcf=filtered_vcf,
        filtered_vcf_plain=filtered_vcf_plain,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ===========================================================
#  Step 05c - Shared SNP Support Classification
# ===========================================================

def consensus_classify_snps_chrom_template(species, chrom, bam_dir, output_dir,
                                           shared_log_dir, script_path, min_bq):
    """Classify clean SNP support for one chromosome of one species."""
    bam = f"{bam_dir}/{species}_final.bam"
    filtered_vcf_plain = f"{output_dir}/{species}_snps_filtered.vcf"
    shard_dir = f"{output_dir}/classify_shards/{species}"
    chrom_vcf = f"{shard_dir}/{species}_{chrom}_snps_filtered.vcf"
    chrom_tsv = f"{shard_dir}/{species}_{chrom}_classification.tsv"
    chrom_stats = f"{shard_dir}/{species}_{chrom}_classification_stats.tsv"
    inputs = [
        f"{shared_log_dir}/05b_filter_snps_{species}.DONE",
    ]
    outputs = [
        f"{shared_log_dir}/05c_classify_snps_{species}_{chrom}.DONE",
        chrom_tsv,
        chrom_stats,
    ]
    options = {
        'cores': 1,
        'memory': '8g',
        'walltime': '04:00:00',
        'account': 'spider2'
    }
    spec = """
    set -eo pipefail
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate python_phylo
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {shard_dir} {shared_log_dir}

    echo "STEP slice_vcf START: $(date)"
    awk -v c={chrom} '/^#/ || $1==c' {filtered_vcf_plain} > {chrom_vcf}
    echo "STEP slice_vcf FINISH: $(date)"

    echo "STEP classify_snps START: $(date)"
    python3 {script_path}/classify_consensus_snps.py \
        --vcf {chrom_vcf} \
        --bam {bam} \
        --min-base-quality {min_bq} \
        --output {chrom_tsv} \
        --stats {chrom_stats}
    echo "STEP classify_snps FINISH: $(date)"
    test -s {chrom_tsv}
    test -s {chrom_stats}
    echo "Classified SNPs ({chrom}): $(awk 'NR>1 {{n++}} END {{print n+0}}' {chrom_tsv})"

    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        bam=bam,
        shard_dir=shard_dir,
        shared_log_dir=shared_log_dir,
        filtered_vcf_plain=filtered_vcf_plain,
        chrom=chrom,
        chrom_vcf=chrom_vcf,
        chrom_tsv=chrom_tsv,
        chrom_stats=chrom_stats,
        script_path=script_path,
        min_bq=min_bq,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def consensus_classify_snps_merge_template(species, chroms, output_dir,
                                           shared_log_dir, script_path):
    """Merge per-chromosome classification shards into whole-genome TSVs.

    The shards are concatenated in the order chromosomes first appear in
    the input filtered VCF so the merged TSV pairs row-for-row with the
    VCF that select_consensus_mode.py consumes downstream.
    """
    shard_dir = f"{output_dir}/classify_shards/{species}"
    filtered_vcf_plain = f"{output_dir}/{species}_snps_filtered.vcf"
    chrom_stats = [f"{shard_dir}/{species}_{c}_classification_stats.tsv" for c in chroms]
    classification_tsv = f"{output_dir}/{species}_snp_classification.tsv"
    classification_stats = f"{output_dir}/{species}_snp_classification_stats.tsv"
    inputs = [
        f"{shared_log_dir}/05c_classify_snps_{species}_{c}.DONE" for c in chroms
    ]
    outputs = [
        f"{shared_log_dir}/05c_classify_snps_{species}.DONE",
        classification_tsv,
        classification_stats,
    ]
    options = {
        'cores': 1,
        'memory': '2g',
        'walltime': '08:00:00',
        'account': 'spider2'
    }
    chrom_stats_args = " ".join(chrom_stats)
    spec = """
    set -eo pipefail
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate python_phylo
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {output_dir} {shared_log_dir}

    echo "STEP merge_classification_tsv START: $(date)"
    VCF_CHROM_ORDER=$(awk '!/^#/ && !seen[$1]++ {{print $1}}' {filtered_vcf_plain})
    FIRST_CHROM=$(echo "$VCF_CHROM_ORDER" | head -n 1)
    if [ -z "$FIRST_CHROM" ]; then
        echo "ERROR: no records in {filtered_vcf_plain}" >&2
        exit 1
    fi
    head -n 1 {shard_dir}/{species}_${{FIRST_CHROM}}_classification.tsv > {classification_tsv}
    for chrom in $VCF_CHROM_ORDER; do
        shard={shard_dir}/{species}_${{chrom}}_classification.tsv
        if [ ! -f "$shard" ]; then
            echo "ERROR: missing shard $shard for chromosome $chrom" >&2
            exit 1
        fi
        tail -n +2 "$shard" >> {classification_tsv}
    done
    echo "STEP merge_classification_tsv FINISH: $(date)"

    echo "STEP merge_classification_stats START: $(date)"
    python3 {script_path}/merge_classification_stats.py \
        --inputs {chrom_stats_args} \
        --output {classification_stats}
    echo "STEP merge_classification_stats FINISH: $(date)"

    test -s {classification_tsv}
    test -s {classification_stats}
    VCF_RECORDS=$(awk '!/^#/ {{n++}} END {{print n+0}}' {filtered_vcf_plain})
    TSV_RECORDS=$(awk 'NR>1 {{n++}} END {{print n+0}}' {classification_tsv})
    if [ "$VCF_RECORDS" != "$TSV_RECORDS" ]; then
        echo "ERROR: VCF has $VCF_RECORDS records but merged TSV has $TSV_RECORDS" >&2
        exit 1
    fi
    echo "Classified SNPs (merged): $TSV_RECORDS"
    awk -F '\\t' '$1=="het_strict_filtered_fraction" {{print "Strict-filtered het fraction: "$2}}' {classification_stats}

    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        output_dir=output_dir,
        shared_log_dir=shared_log_dir,
        shard_dir=shard_dir,
        species=species,
        filtered_vcf_plain=filtered_vcf_plain,
        chrom_stats_args=chrom_stats_args,
        classification_tsv=classification_tsv,
        classification_stats=classification_stats,
        script_path=script_path,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ===========================================================
#  Step 05d - Mode-Specific SNP Selection
# ===========================================================

def consensus_select_snps_template(species, output_dir, shared_consensus_dir,
                                   shared_log_dir, log_dir, script_path,
                                   consensus_mode, support_policy,
                                   branch_label, random_seed):
    """Select branch-specific ALT SNPs from shared SNP classifications."""
    inputs = [
        f"{shared_log_dir}/05b_filter_snps_{species}.DONE",
        f"{shared_log_dir}/05c_classify_snps_{species}.DONE",
    ]
    filtered_vcf_plain = f"{shared_consensus_dir}/{species}_snps_filtered.vcf"
    selected_vcf_plain = f"{output_dir}/{species}_{branch_label}_selected_snps.vcf"
    mode_vcf = f"{output_dir}/{species}_{branch_label}_selected_snps.vcf.gz"
    classification_tsv = f"{shared_consensus_dir}/{species}_snp_classification.tsv"
    rejected_sites_bed = f"{output_dir}/{species}_{branch_label}_rejected_sites.bed"
    outputs = [
        f"{log_dir}/05d_select_snps_{species}.DONE",
        mode_vcf,
        rejected_sites_bed,
    ]
    options = {
        'cores': 1,
        'memory': '8g',
        'walltime': '08:00:00',
        'account': 'spider2'
    }
    spec = """
    set -eo pipefail
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {output_dir} {log_dir}

    conda activate python_phylo
    echo "STEP select_consensus_mode START: $(date)"
    python3 {script_path}/select_consensus_mode.py \
        --mode {consensus_mode} \
        --support-policy {support_policy} \
        --seed {random_seed} \
        --species {species} \
        --classification-tsv {classification_tsv} \
        --rejected-bed {rejected_sites_bed} \
        < {filtered_vcf_plain} \
        > {selected_vcf_plain}
    echo "STEP select_consensus_mode FINISH: $(date)"

    test -s {selected_vcf_plain}
    grep -q '^#CHROM' {selected_vcf_plain}
    SELECTED_COUNT=$(awk 'BEGIN {{n=0}} !/^#/ {{n++}} END {{print n}}' {selected_vcf_plain})
    echo "Selected SNPs ({branch_label}): $SELECTED_COUNT"

    conda activate bcftools
    echo "STEP compress_selected_vcf START: $(date)"
    bcftools view -Oz -o {mode_vcf} {selected_vcf_plain}
    echo "STEP compress_selected_vcf FINISH: $(date)"
    echo "STEP index_selected_vcf START: $(date)"
    bcftools index {mode_vcf}
    echo "STEP index_selected_vcf FINISH: $(date)"
    echo "Rejected SNP sites ({branch_label}): $(wc -l < {rejected_sites_bed})"

    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        output_dir=output_dir,
        log_dir=log_dir,
        filtered_vcf_plain=filtered_vcf_plain,
        selected_vcf_plain=selected_vcf_plain,
        classification_tsv=classification_tsv,
        mode_vcf=mode_vcf,
        rejected_sites_bed=rejected_sites_bed,
        script_path=script_path,
        consensus_mode=consensus_mode,
        support_policy=support_policy,
        branch_label=branch_label,
        random_seed=random_seed,
        species=species,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ===========================================================
#  Step 05e - Mode-Specific Mask Building
# ===========================================================

def consensus_mask_template(species, output_dir, callable_depth_dir,
                            shared_log_dir, log_dir, script_path,
                            branch_label):
    """Merge shared depth mask and branch-specific rejected SNP sites."""
    inputs = [
        f"{shared_log_dir}/05a_callable_depth_{species}.DONE",
        f"{log_dir}/05d_select_snps_{species}.DONE",
    ]
    uncallable_bed = f"{callable_depth_dir}/{species}_uncallable_depth.bed"
    rejected_sites_bed = f"{output_dir}/{species}_{branch_label}_rejected_sites.bed"
    mask_bed = f"{output_dir}/{species}_{branch_label}_mask.bed"
    outputs = [
        f"{log_dir}/05e_build_mask_{species}.DONE",
        mask_bed,
    ]
    options = {
        'cores': 1,
        'memory': '32g',
        'walltime': '12:00:00',
        'account': 'spider2'
    }
    spec = """
    set -eo pipefail
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {output_dir} {log_dir}

    conda activate python_phylo
    echo "STEP merge_mask_bed START: $(date)"
    python3 {script_path}/merge_bed_files.py \
        --primary {uncallable_bed} \
        --secondary {rejected_sites_bed} \
        > {mask_bed}
    echo "STEP merge_mask_bed FINISH: $(date)"
    echo "Final mask regions ({branch_label}): $(wc -l < {mask_bed})"

    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        output_dir=output_dir,
        log_dir=log_dir,
        uncallable_bed=uncallable_bed,
        rejected_sites_bed=rejected_sites_bed,
        mask_bed=mask_bed,
        script_path=script_path,
        branch_label=branch_label,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ===========================================================
#  Step 05f - Mode-Specific Consensus FASTA
# ===========================================================

def bcftools_consensus_template(ref_genome, species, output_dir, log_dir,
                                branch_label):
    """Build a branch-specific masked consensus FASTA."""
    inputs = [
        f"{log_dir}/05d_select_snps_{species}.DONE",
        f"{log_dir}/05e_build_mask_{species}.DONE",
    ]
    consensus_fa = f"{output_dir}/{species}_consensus.fa"
    outputs = [f"{log_dir}/05_consensus_{species}.DONE", consensus_fa]
    options = {
        'cores': 1,
        'memory': '8g',
        'walltime': '12:00:00',
        'account': 'spider2'
    }
    mode_vcf = f"{output_dir}/{species}_{branch_label}_selected_snps.vcf.gz"
    mask_bed = f"{output_dir}/{species}_{branch_label}_mask.bed"
    spec = """
    set -eo pipefail
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {output_dir} {log_dir}

    conda activate bcftools
    echo "STEP bcftools_consensus START: $(date)"
    bcftools consensus --mask {mask_bed} -f {ref_genome} {mode_vcf} > {consensus_fa}
    echo "STEP bcftools_consensus FINISH: $(date)"

    conda activate samtools117
    echo "STEP samtools_faidx START: $(date)"
    samtools faidx {consensus_fa}
    echo "STEP samtools_faidx FINISH: $(date)"

    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        ref_genome=ref_genome,
        mode_vcf=mode_vcf,
        output_dir=output_dir,
        log_dir=log_dir,
        mask_bed=mask_bed,
        consensus_fa=consensus_fa,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ===========================================================
#  Step 06 - CDS Extraction Per Gene
# ===========================================================

def extract_cds_template(ref_genome, ref_name, sp1_name, sp1_consensus, sp2_name, sp2_consensus,
                         gff3, output_dir, log_dir, script_path, min_cov_frac):
    """Extract per-gene CDS alignments from trio genomes."""
    inputs = [
        f"{log_dir}/05_consensus_{sp1_name}.DONE",
        f"{log_dir}/05_consensus_{sp2_name}.DONE",
        ref_genome,
        sp1_consensus,
        sp2_consensus,
        gff3
    ]
    outputs = [f"{log_dir}/06_extract_cds.DONE"]
    options = {
        'cores': 1,
        'memory': '16g',
        'walltime': '08:00:00',
        'account': 'spider2'
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate python_phylo
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {output_dir} {log_dir}
    python3 {script_path}/extract_cds_per_gene.py \
        --ref {ref_genome} \
        --ref_name {ref_name} \
        --sp1_name {sp1_name} --sp1_fa {sp1_consensus} \
        --sp2_name {sp2_name} --sp2_fa {sp2_consensus} \
        --gff3 {gff3} \
        --output_dir {output_dir} \
        --min_cov_frac {min_cov_frac}
    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        ref_genome=ref_genome,
        ref_name=ref_name,
        sp1_name=sp1_name,
        sp1_consensus=sp1_consensus,
        sp2_name=sp2_name,
        sp2_consensus=sp2_consensus,
        gff3=gff3,
        output_dir=output_dir,
        log_dir=log_dir,
        script_path=script_path,
        min_cov_frac=min_cov_frac,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ===========================================================
#  Step 07 - Prepare Gene Lists (autosomal vs X-linked)
# ===========================================================

def prepare_gene_lists_template(gff3, passing_genes, output_dir, log_dir, script_path,
                                autosomal_chromosomes, x_chromosome_prefixes):
    """Split passing genes into autosomal vs X-linked lists."""
    inputs = [
        f"{log_dir}/06_extract_cds.DONE"
    ]
    outputs = [
        f"{log_dir}/07_prepare_gene_lists.DONE",
        f"{output_dir}/auto_passing_genes.txt"
    ]
    options = {
        'cores': 1,
        'memory': '4g',
        'walltime': '00:30:00',
        'account': 'spider2'
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate python_phylo
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {output_dir} {log_dir}
    python3 {script_path}/prepare_gene_lists.py \
        --gff3 {gff3} \
        --passing {passing_genes} \
        --output_dir {output_dir} \
        --autosomal_chromosomes {autosomal_chromosomes} \
        --x_chromosome_prefixes {x_chromosome_prefixes}
    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        gff3=gff3,
        passing_genes=passing_genes,
        output_dir=output_dir,
        log_dir=log_dir,
        script_path=script_path,
        autosomal_chromosomes=autosomal_chromosomes,
        x_chromosome_prefixes=x_chromosome_prefixes,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ===========================================================
#  Step 08a - Bootstrap Sampling (per replicate)
# ===========================================================

def bootstrap_sample_template(gene_list, replicate_id, work_dir, log_dir):
    """Sample genes for a bootstrap replicate (with replacement) or copy for auto_all."""
    inputs = [
        f"{log_dir}/07_prepare_gene_lists.DONE",
        gene_list
    ]
    id_file = f"{work_dir}/{replicate_id}_id.txt"
    outputs = [
        f"{log_dir}/08a_sample_{replicate_id}.DONE",
        id_file
    ]
    options = {
        'cores': 1,
        'memory': '100m',
        'walltime': '00:05:00',
        'account': 'spider2'
    }
    # For auto_all: copy the full gene list (no sampling)
    # For bootstrap replicates: sample with replacement using shuf -r
    # n_genes is determined at runtime from the gene list (not hardcoded)
    if replicate_id == 'auto_all':
        sample_cmd = "cp {gene_list} {id_file}".format(
            gene_list=gene_list, id_file=id_file
        )
    else:
        sample_cmd = "N_GENES=$(wc -l < {gene_list}) && shuf -r -n $N_GENES {gene_list} > {id_file}".format(
            gene_list=gene_list, id_file=id_file
        )
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate paml
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {work_dir} {log_dir}
    cd {work_dir}
    {sample_cmd}
    echo "Sampled $(wc -l < {id_file}) genes for {replicate_id}"
    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        work_dir=work_dir,
        log_dir=log_dir,
        sample_cmd=sample_cmd,
        id_file=id_file,
        replicate_id=replicate_id,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ===========================================================
#  Step 08b - Concat + Filter + PAML (per replicate)
# ===========================================================

def bootstrap_paml_template(replicate_id, id_list_file, fasta_dir, work_dir, log_dir, script_path, trio_tree):
    """Run full concat+filter+phylip+codeml+parse pipeline via wrapper script."""
    inputs = [
        f"{log_dir}/08a_sample_{replicate_id}.DONE",
        id_list_file
    ]
    outputs = [
        f"{log_dir}/08b_paml_{replicate_id}.DONE"
    ]
    options = {
        'cores': 1,
        'memory': '8g',
        'walltime': '08:00:00',
        'account': 'spider2'
    }
    spec = """
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {work_dir} {log_dir}
    bash {script_path}/concat_filter_paml.sh \
        {replicate_id} \
        {id_list_file} \
        {fasta_dir} \
        {work_dir} \
        {script_path} \
        '{trio_tree}'
    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        replicate_id=replicate_id,
        id_list_file=id_list_file,
        fasta_dir=fasta_dir,
        work_dir=work_dir,
        log_dir=log_dir,
        script_path=script_path,
        trio_tree=trio_tree,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ===========================================================
#  Step 09 - Summarize Bootstrap Results
# ===========================================================

def summarize_bootstrap_template(bootstrap_dir, n_bootstrap, output_dir, log_dir, script_path,
                                  paml_done_logs):
    """Summarize bootstrap dN/dS results with confidence intervals."""
    inputs = paml_done_logs
    outputs = [
        f"{log_dir}/09_summarize_bootstrap.DONE",
        f"{output_dir}/branch/auto_branch_summary.tsv",
        f"{output_dir}/branch/auto_branch_bootstrap_values.tsv",
        f"{output_dir}/pairwise/auto_pairwise_summary.tsv",
        f"{output_dir}/pairwise/auto_pairwise_bootstrap_values.tsv"
    ]
    options = {
        'cores': 1,
        'memory': '4g',
        'walltime': '01:00:00',
        'account': 'spider2'
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate python_phylo
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {output_dir} {log_dir}
    python3 {script_path}/summarize_bootstrap.py \
        --bootstrap_dir {bootstrap_dir} \
        --n_bootstrap {n_bootstrap} \
        --output_dir {output_dir}
    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        bootstrap_dir=bootstrap_dir,
        n_bootstrap=n_bootstrap,
        output_dir=output_dir,
        log_dir=log_dir,
        script_path=script_path,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ===========================================================
#  Step 10 - Visualize Pairwise dS
# ===========================================================

def visualize_pairwise_dS_template(summary_dir, output_dir, log_dir, script_path,
                                    trio_name, trio_tree, summary_done_log):
    """Heatmap, unrooted-tree plot, and recalculated branch-length TSV from pairwise dS."""
    inputs = [
        summary_done_log,
        f"{summary_dir}/pairwise/auto_pairwise_summary.tsv",
        f"{summary_dir}/pairwise/auto_pairwise_bootstrap_values.tsv",
    ]
    outputs = [
        f"{log_dir}/10_visualize_pairwise_dS.DONE",
        f"{output_dir}/pairwise_dS_heatmap.png",
        f"{output_dir}/pairwise_dS_heatmap.pdf",
        f"{output_dir}/branch_dS_tree.png",
        f"{output_dir}/branch_dS_tree.pdf",
        f"{output_dir}/branch_dS_lengths.tsv",
        f"{output_dir}/social_lineage_shortening.tsv",
    ]
    options = {
        'cores': 1,
        'memory': '2g',
        'walltime': '00:30:00',
        'account': 'spider2'
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate python_phylo
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {output_dir} {log_dir}
    python3 {script_path}/visualize_pairwise_dS.py \
        --pairwise_summary {summary_dir}/pairwise/auto_pairwise_summary.tsv \
        --pairwise_bootstrap {summary_dir}/pairwise/auto_pairwise_bootstrap_values.tsv \
        --trio_tree '{trio_tree}' \
        --trio_name {trio_name} \
        --output_dir {output_dir}
    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        summary_dir=summary_dir,
        output_dir=output_dir,
        log_dir=log_dir,
        script_path=script_path,
        trio_name=trio_name,
        trio_tree=trio_tree,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# ===========================================================
#  Step 11 - Combined Cross-Trio Visualisation
# ===========================================================

def combine_trio_visualizations_template(trios, output_dir, log_dir, script_path):
    """Build a 2 x N panel figure across multiple trios.

    `trios` is a list of dicts, each with keys:
        trio_name, trio_tree, summary_dir, summary_done_log
    summary_dir is the per-trio Step 09 dir (containing pairwise/...).
    """
    inputs = []
    trio_args = []
    for t in trios:
        inputs.append(t['summary_done_log'])
        inputs.append(f"{t['summary_dir']}/pairwise/auto_pairwise_summary.tsv")
        inputs.append(f"{t['summary_dir']}/pairwise/auto_pairwise_bootstrap_values.tsv")
        trio_args.append(
            f"--trio_name {t['trio_name']} "
            f"--trio_tree '{t['trio_tree']}' "
            f"--pairwise_summary {t['summary_dir']}/pairwise/auto_pairwise_summary.tsv "
            f"--pairwise_bootstrap {t['summary_dir']}/pairwise/auto_pairwise_bootstrap_values.tsv"
        )
    trio_args_str = " \\\n        ".join(trio_args)

    outputs = [
        f"{log_dir}/11_combine_trio_visualizations.DONE",
        f"{output_dir}/combined_pairwise_dS.png",
        f"{output_dir}/combined_pairwise_dS.pdf",
    ]
    options = {
        'cores': 1,
        'memory': '2g',
        'walltime': '00:30:00',
        'account': 'spider2'
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate python_phylo
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {output_dir} {log_dir}
    python3 {script_path}/combine_trio_visualizations.py \\
        {trio_args} \\
        --output_dir {output_dir}
    echo "FINISH: $(date)"
    echo done > {log}
    """.format(
        output_dir=output_dir,
        log_dir=log_dir,
        script_path=script_path,
        trio_args=trio_args_str,
        log=outputs[0]
    )
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
