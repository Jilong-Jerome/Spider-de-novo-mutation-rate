#!/usr/bin/env python3
"""
workflow_templates.py
GWF target templates for the pairwise kinship identification workflow.

Per species the pipeline is:
    filter_snps -> prepare_plink -> { plink_ibd, king_kinship } -> classify -> plot

Each template writes a DONE sentinel as its last action; downstream targets list
the upstream sentinel(s) in their inputs so gwf wires the dependency DAG.
"""
from gwf import AnonymousTarget


def filter_snps_template(sp, vcf, out_dir, log_path, scripts,
                         vcftools_env, python_env, maf_threshold, min_gq,
                         max_missing, account):
    """Keep biallelic autosomal SNPs that are well-genotyped and common.

    vcftools masks low-confidence genotypes (--minGQ) and drops sites by call
    rate (--max-missing) across all individuals; vcf_filter_maf_snps.py then
    keeps autosomal SNPs with MAF > threshold. Autosomes only: that script keeps
    contigs named <name>_<digits>, so sex chromosomes (mim_X1/mim_X2) are dropped.
    """
    snp_vcf = f'{out_dir}/{sp}_autosome_snps_maf.vcf'
    done_file = f'{log_path}/{sp}_filter_snps.DONE'
    inputs = {'vcf': vcf}
    outputs = {'snp_vcf': snp_vcf, 'log': done_file}
    options = {'cores': 1, 'memory': '12g', 'walltime': '12:00:00', 'account': account}
    spec = f"""
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"
mkdir -p {out_dir} {log_path}
cd {out_dir}

conda activate {vcftools_env}
vcftools --gzvcf {vcf} \\
    --min-alleles 2 --max-alleles 2 --remove-filtered-all \\
    --minGQ {min_gq} --max-missing {max_missing} \\
    --recode --recode-INFO-all \\
    --out {sp}_biallelic

conda activate {python_env}
python {scripts}/vcf_filter_maf_snps.py {sp}_biallelic.recode.vcf {snp_vcf}
rm -f {sp}_biallelic.recode.vcf

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def prepare_plink_template(sp, snp_vcf, filter_done, out_dir, log_path, scripts,
                           python_env, plink_env, ld_prune, account):
    """Fill variant IDs, LD-prune, make PLINK binary, remap chrom to integers for KING."""
    named_vcf = f'{out_dir}/{sp}_named.vcf'
    bed = f'{out_dir}/{sp}_pruned.bed'
    bim = f'{out_dir}/{sp}_pruned.bim'
    fam = f'{out_dir}/{sp}_pruned.fam'
    chrommap = f'{out_dir}/{sp}_pruned.chrommap.tsv'
    done_file = f'{log_path}/{sp}_prepare_plink.DONE'
    inputs = {'snp_vcf': snp_vcf, 'filter_done': filter_done}
    outputs = {'bed': bed, 'bim': bim, 'fam': fam, 'chrommap': chrommap, 'log': done_file}
    options = {'cores': 1, 'memory': '8g', 'walltime': '04:00:00', 'account': account}
    spec = f"""
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"
mkdir -p {out_dir} {log_path}
cd {out_dir}

conda activate {python_env}
python {scripts}/vcf_fill_id.py {snp_vcf} {named_vcf}

conda activate {plink_env}
plink --vcf {named_vcf} --double-id --allow-extra-chr \\
    --indep-pairwise {ld_prune} --out {sp}_named
plink --vcf {named_vcf} --double-id --allow-extra-chr \\
    --extract {sp}_named.prune.in --make-bed --out {sp}_pruned

# Remap scaffold names to integer chromosome codes so KING treats them as
# autosomes (asserts no sex chromosome leaked in).
python {scripts}/remap_bim_chrom.py {bim} {chrommap}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def plink_ibd_template(sp, bed, bim, fam, prepare_done, out_dir, log_path,
                       plink_env, account):
    """PLINK --genome IBD estimation (Z0/Z1/Z2/PI_HAT) as an independent cross-check."""
    genome = f'{out_dir}/{sp}_plink.genome'
    done_file = f'{log_path}/{sp}_plink_ibd.DONE'
    inputs = {'bed': bed, 'bim': bim, 'fam': fam, 'prepare_done': prepare_done}
    outputs = {'genome': genome, 'log': done_file}
    options = {'cores': 1, 'memory': '8g', 'walltime': '04:00:00', 'account': account}
    spec = f"""
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {plink_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"
mkdir -p {out_dir} {log_path}
cd {out_dir}

plink --bfile {out_dir}/{sp}_pruned --allow-extra-chr --genome --out {sp}_plink

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def king_kinship_template(sp, bed, bim, fam, prepare_done, out_dir, log_path,
                          king_binary, account):
    """KING kinship coefficient + IBS0, and --related relationship inference."""
    kin0 = f'{out_dir}/{sp}_king.kin0'
    rel_kin0 = f'{out_dir}/{sp}_king_rel.kin0'
    done_file = f'{log_path}/{sp}_king_kinship.DONE'
    inputs = {'bed': bed, 'bim': bim, 'fam': fam, 'prepare_done': prepare_done}
    outputs = {'kin0': kin0, 'log': done_file}
    options = {'cores': 4, 'memory': '8g', 'walltime': '04:00:00', 'account': account}
    spec = f"""
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"
mkdir -p {out_dir} {log_path}
cd {out_dir}

# Kinship coefficient + IBS0 for every pair (.kin0; all pairs are cross-family
# because PLINK --double-id makes every FID unique).
{king_binary} -b {bed} --kinship --ibs --cpus 4 --prefix {sp}_king

# Built-in relationship inference (adds InfType column: PO/FS/2nd/3rd/UN/...).
{king_binary} -b {bed} --related --degree 2 --cpus 4 --prefix {sp}_king_rel

# Guarantee the .kin0 outputs exist even if KING emitted only .kin (defensive).
touch {kin0} {rel_kin0}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def classify_template(sp, kin0, genome, king_done, plink_done, out_dir, log_path,
                      scripts, python_env, thresholds, account):
    """Merge KING + PLINK estimates, classify pairs, compare to known truth."""
    rel_kin0 = f'{out_dir}/{sp}_king_rel.kin0'
    classified = f'{out_dir}/{sp}_kinship_classified.tsv'
    done_file = f'{log_path}/{sp}_classify.DONE'
    inputs = {'kin0': kin0, 'genome': genome,
              'king_done': king_done, 'plink_done': plink_done}
    outputs = {'classified': classified, 'log': done_file}
    options = {'cores': 1, 'memory': '2g', 'walltime': '00:30:00', 'account': account}
    spec = f"""
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {python_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"
mkdir -p {out_dir} {log_path}

python {scripts}/classify_kinship.py \\
    --king_kin0 {kin0} \\
    --king_related {rel_kin0} \\
    --plink_genome {genome} \\
    --output_tsv {classified} \\
    --kinship_dup {thresholds['kinship_dup']} \\
    --kinship_first {thresholds['kinship_first']} \\
    --kinship_second {thresholds['kinship_second']} \\
    --kinship_third {thresholds['kinship_third']} \\
    --ibs0_po_max {thresholds['ibs0_po_max']}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def trio_filter_template(sp, vcf, out_dir, log_path, scripts,
                         vcftools_env, python_env, min_gq, trio_max_missing, account):
    """No-MAF biallelic autosomal SNP VCF for the per-trio Mendelian check.

    Independent of filter_snps (whose MAF>0.1 output cannot recover the
    low-frequency, single-family-segregating sites this test needs). vcftools
    masks low-confidence genotypes (--minGQ) but keeps all sites regardless of
    missingness (--max-missing {trio_max_missing}, 0.0 = keep all); per-trio
    completeness is enforced downstream. vcf_filter_maf_snps.py is run in no-MAF
    mode (3rd arg -1) so every biallelic autosomal SNP is kept (sex chromosomes
    dropped by is_autosome).
    """
    snp_vcf = f'{out_dir}/{sp}_autosome_snps_nomaf.vcf'
    done_file = f'{log_path}/{sp}_trio_filter.DONE'
    inputs = {'vcf': vcf}
    outputs = {'snp_vcf': snp_vcf, 'log': done_file}
    options = {'cores': 1, 'memory': '12g', 'walltime': '12:00:00', 'account': account}
    spec = f"""
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"
mkdir -p {out_dir} {log_path}
cd {out_dir}

conda activate {vcftools_env}
vcftools --gzvcf {vcf} \\
    --min-alleles 2 --max-alleles 2 --remove-filtered-all \\
    --minGQ {min_gq} --max-missing {trio_max_missing} \\
    --recode --recode-INFO-all \\
    --out {sp}_trio_biallelic

conda activate {python_env}
python {scripts}/vcf_filter_maf_snps.py {sp}_trio_biallelic.recode.vcf {snp_vcf} -1
rm -f {sp}_trio_biallelic.recode.vcf

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def trio_check_template(sp, snp_vcf, trio_filter_done, out_dir, log_path, scripts,
                        python_env, min_gq, min_dp, consistency_min, min_sites, account):
    """Per-trio opposite-homozygote Mendelian test (allele-frequency-free)."""
    trio_tsv = f'{out_dir}/{sp}_trio_mendel.tsv'
    done_file = f'{log_path}/{sp}_trio_check.DONE'
    inputs = {'snp_vcf': snp_vcf, 'trio_filter_done': trio_filter_done}
    outputs = {'trio_tsv': trio_tsv, 'log': done_file}
    options = {'cores': 1, 'memory': '4g', 'walltime': '04:00:00', 'account': account}
    spec = f"""
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {python_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"
mkdir -p {out_dir} {log_path}

python {scripts}/trio_mendel_check.py \\
    --vcf {snp_vcf} \\
    --species {sp} \\
    --output_tsv {trio_tsv} \\
    --min_gq {min_gq} \\
    --min_dp {min_dp} \\
    --consistency_min {consistency_min} \\
    --min_sites {min_sites}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def trio_plot_template(sp, trio_tsv, trio_check_done, out_dir, log_path, scripts,
                       python_env, consistency_min, account):
    """Figure for the per-trio Mendelian validation."""
    pdf = f'{out_dir}/{sp}_trio_mendel.pdf'
    done_file = f'{log_path}/{sp}_trio_plot.DONE'
    inputs = {'trio_tsv': trio_tsv, 'trio_check_done': trio_check_done}
    outputs = {'pdf': pdf, 'log': done_file}
    options = {'cores': 1, 'memory': '2g', 'walltime': '00:30:00', 'account': account}
    spec = f"""
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {python_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"
mkdir -p {out_dir} {log_path}

python {scripts}/plot_trio_mendel.py \\
    --input_tsv {trio_tsv} \\
    --output_pdf {pdf} \\
    --species {sp} \\
    --consistency_min {consistency_min}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def combined_trio_plot_template(species_tsvs, check_dones, combined_dir, log_path,
                                scripts, python_env, consistency_min, exclude, account):
    """Combined all-species supplemental figure for the per-trio Mendelian test.

    One target across every species: reads all {SP}_trio_mendel.tsv files and
    writes a single multi-row PDF (one row per species, two columns; the
    per-species histogram panel is dropped). `exclude` is a list of 'SP:family'
    tokens to omit (the supplement drops BIC family6).
    """
    pdf = f'{combined_dir}/all_species_trio_mendel.pdf'
    done_file = f'{log_path}/combined_trio_plot.DONE'
    inputs = dict(check_dones)
    inputs.update(species_tsvs)
    outputs = {'pdf': pdf, 'log': done_file}
    options = {'cores': 1, 'memory': '2g', 'walltime': '00:30:00', 'account': account}
    tsv_args = ' '.join(species_tsvs[sp] for sp in sorted(species_tsvs))
    exclude_args = ' '.join(f'--exclude {e}' for e in exclude)
    spec = f"""
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {python_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"
mkdir -p {combined_dir} {log_path}

python {scripts}/plot_combined_trio_mendel.py \\
    --input_tsvs {tsv_args} \\
    --output_pdf {pdf} \\
    --consistency_min {consistency_min} \\
    {exclude_args}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def plot_template(sp, classified, classify_done, out_dir, log_path,
                  scripts, python_env, thresholds, account):
    """Validation figure: kinship-vs-IBS0, Z0-vs-Z1, and confusion matrix."""
    pdf = f'{out_dir}/{sp}_kinship_validation.pdf'
    done_file = f'{log_path}/{sp}_plot.DONE'
    inputs = {'classified': classified, 'classify_done': classify_done}
    outputs = {'pdf': pdf, 'log': done_file}
    options = {'cores': 1, 'memory': '2g', 'walltime': '00:30:00', 'account': account}
    spec = f"""
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate {python_env}

echo "START: $(date)"
echo "JobID: $SLURM_JOBID"
mkdir -p {out_dir} {log_path}

python {scripts}/plot_kinship.py \\
    --input_tsv {classified} \\
    --output_pdf {pdf} \\
    --species {sp} \\
    --kinship_first {thresholds['kinship_first']} \\
    --kinship_second {thresholds['kinship_second']} \\
    --kinship_third {thresholds['kinship_third']} \\
    --ibs0_po_max {thresholds['ibs0_po_max']}

echo "DONE: $(date)"
echo done > {done_file}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
