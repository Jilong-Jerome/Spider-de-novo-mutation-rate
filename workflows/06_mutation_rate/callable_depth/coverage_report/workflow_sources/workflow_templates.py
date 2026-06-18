#!/bin/env python3
from gwf import AnonymousTarget
import os


def contigs_template(species, bam, out_dir, log_path, scripts, conda_env, account):
    """Extract chrom<TAB>length from the BAM @SQ header for one representative BAM."""
    lengths_tsv = f"{out_dir}/{species}.chrom_lengths.tsv"
    done = f"{log_path}/{species}_contigs.DONE"
    inputs = {"bam": bam}
    outputs = {"lengths": lengths_tsv, "done": done}
    options = {"cores": 1, "memory": "2g", "walltime": "00:30:00", "account": account}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate {conda_env}
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    mkdir -p {out_dir} {log_path}

    samtools view -H {bam} \\
        | awk -F'\\t' '/^@SQ/ {{ sn=""; ln=""; for (i=1;i<=NF;i++) {{ if ($i ~ /^SN:/) sn=substr($i,4); if ($i ~ /^LN:/) ln=substr($i,4) }} if (sn!="" && ln!="") print sn"\\t"ln }}' \\
        > {lengths_tsv}

    echo "DONE: $(date)"
    echo done > {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def bam_depth_template(species, sample, chrom, bam, lengths_tsv, lengths_done,
                       out_dir, log_path, scripts, conda_env, account):
    """One `samtools depth -d 0 -r {chrom}` pass per (individual x chromosome).

    Uses the .bai index for random access so each job only scans a single
    chromosome (32 x 23 = 736 concurrent jobs for LIN). `-@` adds BGZF
    decompression threads.
    """
    bam_dir = f"{out_dir}/bam"
    out_tsv = f"{bam_dir}/{sample}.{chrom}.bam_cov.tsv"
    done = f"{log_path}/{species}_bamcov_{sample}_{chrom}.DONE"
    inputs = {"bam": bam, "lengths": lengths_tsv, "lengths_done": lengths_done}
    outputs = {"tsv": out_tsv, "done": done}
    options = {"cores": 4, "memory": "4g", "walltime": "06:00:00", "account": account}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate {conda_env}
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    mkdir -p {bam_dir} {log_path}

    samtools depth -d 0 -@ 3 -r {chrom} {bam} \\
        | python3 {scripts}/bam_depth_stats.py \\
            --sample {sample} \\
            --chrom {chrom} \\
            --lengths {lengths_tsv} \\
            --out {out_tsv}

    echo "DONE: $(date)"
    echo done > {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def vcf_dp_template(species, chrom, vcf, samples, out_dir, log_path,
                    scripts, conda_env, account):
    """One job per chromosome: extract per-sample DP via tabix region -> DP stats."""
    vcf_dir = f"{out_dir}/vcf"
    out_tsv = f"{vcf_dir}/{chrom}.vcf_dp.tsv"
    hist_tsv = f"{vcf_dir}/{chrom}.vcf_dp_hist.tsv"
    done = f"{log_path}/{species}_vcfdp_{chrom}.DONE"
    samples_csv = ",".join(samples)
    inputs = {"vcf": vcf}
    outputs = {"tsv": out_tsv, "hist": hist_tsv, "done": done}
    options = {"cores": 1, "memory": "4g", "walltime": "12:00:00", "account": account}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate {conda_env}
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    mkdir -p {vcf_dir} {log_path}

    # Resolve symlink so the adjacent .tbi index is discoverable.
    VCF=$(readlink -f {vcf})

    bcftools query -r {chrom} -f '[%DP\\t]\\n' "$VCF" \\
        | python3 {scripts}/vcf_dp_stats.py \\
            --chrom {chrom} \\
            --samples '{samples_csv}' \\
            --out {out_tsv} \\
            --hist-out {hist_tsv}

    echo "DONE: $(date)"
    echo done > {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def merge_template(species, samples, chromosomes, out_dir, log_path,
                   scripts, conda_env, account):
    """Inner-join BAM and VCF stats into the per-species report."""
    bam_dir = f"{out_dir}/bam"
    vcf_dir = f"{out_dir}/vcf"
    report = f"{out_dir}/{species}_coverage_report.tsv"
    done = f"{log_path}/{species}_merge.DONE"
    inputs = {}
    inputs.update({f"bam_{s}_{c}": f"{log_path}/{species}_bamcov_{s}_{c}.DONE"
                   for s in samples for c in chromosomes})
    inputs.update({f"vcf_{c}": f"{log_path}/{species}_vcfdp_{c}.DONE" for c in chromosomes})
    outputs = {"report": report, "done": done}
    options = {"cores": 1, "memory": "4g", "walltime": "01:00:00", "account": account}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate {conda_env}
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    mkdir -p {out_dir} {log_path}

    python3 {scripts}/merge_report.py \\
        --bam-dir {bam_dir} \\
        --vcf-dir {vcf_dir} \\
        --out {report}

    echo "DONE: $(date)"
    echo done > {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def autosome_dp_template(species, autosomes, out_dir, log_path,
                         scripts, conda_env, account):
    """Per-individual autosomal DP report (all sites, DP=0 included).

    Pools the per-sample all-sites DP histograms across the autosomal chromosomes
    only, then writes one row per individual with n_sites, mean_DP_all and the
    exact pooled median_DP_all. Depends on the autosomal vcfdp jobs.
    """
    vcf_dir = f"{out_dir}/vcf"
    report = f"{out_dir}/{species}_autosome_dp_per_individual.tsv"
    done = f"{log_path}/{species}_autosome_dp.DONE"
    autosomes_csv = ",".join(autosomes)
    inputs = {f"vcf_{c}": f"{log_path}/{species}_vcfdp_{c}.DONE" for c in autosomes}
    outputs = {"report": report, "done": done}
    options = {"cores": 1, "memory": "4g", "walltime": "01:00:00", "account": account}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate {conda_env}
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    mkdir -p {out_dir} {log_path}

    python3 {scripts}/autosome_dp_report.py \\
        --vcf-dir {vcf_dir} \\
        --autosomes '{autosomes_csv}' \\
        --out {report}

    echo "DONE: $(date)"
    echo done > {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def combined_plot_template(species_tsvs, species_dones, combined_dir, log_path,
                           scripts, conda_env, account):
    """One cross-species barplot of per-individual autosomal DP.

    `species_tsvs` is a list of (species, tsv_path); `species_dones` is the list of
    per-species `{SP}_autosome_dp.DONE` files this step waits on. Produces a merged
    table and a barplot (one bar per individual, coloured by species, median dot,
    and a horizontal grand-mean line across all individuals of all species).
    """
    merged_tsv = f"{combined_dir}/all_species_autosome_dp_per_individual.tsv"
    plot_png = f"{combined_dir}/all_species_autosome_dp_per_individual_barplot.png"
    done = f"{log_path}/coverage_combined_dp_plot.DONE"
    spec_args = " ".join(f"{sp}:{path}" for sp, path in species_tsvs)
    inputs = {f"done_{i}": d for i, d in enumerate(species_dones)}
    outputs = {"merged": merged_tsv, "plot": plot_png, "done": done}
    options = {"cores": 1, "memory": "4g", "walltime": "01:00:00", "account": account}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate {conda_env}
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    mkdir -p {combined_dir} {log_path}

    python3 {scripts}/combined_dp_plot.py \\
        --inputs {spec_args} \\
        --merged-out {merged_tsv} \\
        --plot-out {plot_png}

    echo "DONE: $(date)"
    echo done > {done}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
