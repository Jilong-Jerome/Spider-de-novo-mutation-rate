from gwf import Workflow
import os
import re
import subprocess
import yaml
from workflow_templates import (
    contigs_template,
    bam_depth_template,
    vcf_dp_template,
    merge_template,
    autosome_dp_template,
)


def get_vcf_samples(vcf, conda_env):
    """Read the sample list from the VCF at graph-build time via bcftools query -l.

    Returns samples in VCF order. This is the canonical sample set for the report:
    BAM-only individuals (not in the VCF) are never created as targets.
    """
    out = subprocess.run(
        ["conda", "run", "-n", conda_env, "bcftools", "query", "-l", vcf],
        check=True, capture_output=True, text=True,
    ).stdout
    return [s for s in out.splitlines() if s.strip()]


def coverage_workflow(config_file, gwf):
    with open(config_file) as f:
        cfg = yaml.safe_load(f)

    sp = cfg["species"]
    account = cfg["account"]
    bam_dir = cfg["bam_dir"]
    vcf = cfg["vcf"]
    scripts = cfg["scripts_path"]
    chromosomes = cfg["chromosomes"]
    out_dir = os.path.join(cfg["output_directory_path"], sp)
    log_path = cfg["log_directory_path"]
    env_sam = cfg["conda_samtools"]
    env_bcf = cfg["conda_bcftools"]
    env_py = cfg["conda_python"]

    samples = get_vcf_samples(vcf, env_bcf)

    # Autosomes = chromosomes that are not sex chromosomes ({sp}_X{N}). The
    # per-individual DP report pools only these.
    autosomes = [c for c in chromosomes if not re.search(r"_X\d+$", c)]

    # representative BAM for the contigs/lengths step (first VCF sample)
    rep_bam = os.path.join(bam_dir, f"{samples[0]}_final.bam")
    lengths_tsv = f"{out_dir}/{sp}.chrom_lengths.tsv"
    lengths_done = f"{log_path}/{sp}_contigs.DONE"

    gwf.target_from_template(
        name=f"{sp}_contigs",
        template=contigs_template(sp, rep_bam, out_dir, log_path, scripts, env_sam, account),
    )

    for s in samples:
        for c in chromosomes:
            gwf.target_from_template(
                name=f"{sp}_bamcov_{s}_{c}",
                template=bam_depth_template(
                    sp, s, c, os.path.join(bam_dir, f"{s}_final.bam"),
                    lengths_tsv, lengths_done, out_dir, log_path, scripts, env_sam, account),
            )

    for c in chromosomes:
        gwf.target_from_template(
            name=f"{sp}_vcfdp_{c}",
            template=vcf_dp_template(sp, c, vcf, samples, out_dir, log_path,
                                     scripts, env_bcf, account),
        )

    gwf.target_from_template(
        name=f"{sp}_merge",
        template=merge_template(sp, samples, chromosomes, out_dir, log_path,
                                scripts, env_py, account),
    )

    gwf.target_from_template(
        name=f"{sp}_autosome_dp",
        template=autosome_dp_template(sp, autosomes, out_dir, log_path,
                                      scripts, env_py, account),
    )

    # Metadata the top-level combined plot step needs (wired in workflow.py).
    meta = {
        "species": sp,
        "autosome_dp_tsv": f"{out_dir}/{sp}_autosome_dp_per_individual.tsv",
        "autosome_dp_done": f"{log_path}/{sp}_autosome_dp.DONE",
        "log_path": log_path,
        "scripts": scripts,
        "conda_python": env_py,
        "account": account,
        "steps_root": cfg["output_directory_path"],
    }
    return gwf, meta
