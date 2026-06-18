#!/usr/bin/env python3
"""
workflow_sources.py
Build GWF targets for the pairwise kinship identification workflow.

For each species in the config, registers the 6-step pipeline:
    {SP}_filter_snps -> {SP}_prepare_plink
        -> {SP}_plink_ibd / {SP}_king_kinship
        -> {SP}_classify -> {SP}_plot
"""
import yaml

from workflow_templates import (
    filter_snps_template,
    prepare_plink_template,
    plink_ibd_template,
    king_kinship_template,
    classify_template,
    plot_template,
    trio_filter_template,
    trio_check_template,
    trio_plot_template,
    combined_trio_plot_template,
)


def kinship_workflow(config_file, gwf):
    with open(config_file) as f:
        config = yaml.safe_load(f)

    account = config['account']
    output_path = config['output_directory_path']
    log_path = config['log_directory_path']
    scripts = config['scripts_path']

    python_env = config['python_conda_env']
    plink_env = config['plink_conda_env']
    vcftools_env = config['vcftools_conda_env']
    king_binary = config['king_binary']

    # Toggle the frequency-based KING/PLINK all-pairs branch. Default True keeps
    # backward compatibility (the completed MIM/BIC runs). When False, only the
    # allele-frequency-free per-trio Mendelian branch is registered -- used for
    # the species where the inbreeding-sensitive estimators are not wanted.
    run_king_plink = config.get('run_king_plink', True)

    maf_threshold = config['maf_threshold']
    min_gq = config['min_gq']
    max_missing = config['max_missing']
    ld_prune = config['ld_prune']

    # Per-trio Mendelian validation branch.
    trio_max_missing = config['trio_max_missing']
    trio_min_gq = config['trio_min_gq']
    trio_min_dp = config['trio_min_dp']
    trio_consistency_min = config['trio_consistency_min']
    trio_min_sites = config['trio_min_sites']
    thresholds = {
        'kinship_dup': config['kinship_dup'],
        'kinship_first': config['kinship_first'],
        'kinship_second': config['kinship_second'],
        'kinship_third': config['kinship_third'],
        'ibs0_po_max': config['ibs0_po_max'],
    }

    for sp, sp_cfg in config['species'].items():
        vcf = sp_cfg['vcf']
        out_dir = f'{output_path}/{sp}'

        bed = f'{out_dir}/{sp}_pruned.bed'
        bim = f'{out_dir}/{sp}_pruned.bim'
        fam = f'{out_dir}/{sp}_pruned.fam'
        kin0 = f'{out_dir}/{sp}_king.kin0'
        genome = f'{out_dir}/{sp}_plink.genome'
        classified = f'{out_dir}/{sp}_kinship_classified.tsv'

        filter_done = f'{log_path}/{sp}_filter_snps.DONE'
        prepare_done = f'{log_path}/{sp}_prepare_plink.DONE'
        plink_done = f'{log_path}/{sp}_plink_ibd.DONE'
        king_done = f'{log_path}/{sp}_king_kinship.DONE'
        classify_done = f'{log_path}/{sp}_classify.DONE'

        snp_vcf = f'{out_dir}/{sp}_autosome_snps_maf.vcf'

        # Trio branch paths.
        trio_snp_vcf = f'{out_dir}/{sp}_autosome_snps_nomaf.vcf'
        trio_tsv = f'{out_dir}/{sp}_trio_mendel.tsv'
        trio_filter_done = f'{log_path}/{sp}_trio_filter.DONE'
        trio_check_done = f'{log_path}/{sp}_trio_check.DONE'

        # --- Frequency-based KING/PLINK all-pairs branch (optional) ---
        if run_king_plink:
            gwf.target_from_template(
                name=f'{sp}_filter_snps',
                template=filter_snps_template(
                    sp=sp, vcf=vcf, out_dir=out_dir, log_path=log_path, scripts=scripts,
                    vcftools_env=vcftools_env, python_env=python_env,
                    maf_threshold=maf_threshold, min_gq=min_gq, max_missing=max_missing,
                    account=account,
                ),
            )

            gwf.target_from_template(
                name=f'{sp}_prepare_plink',
                template=prepare_plink_template(
                    sp=sp, snp_vcf=snp_vcf, filter_done=filter_done,
                    out_dir=out_dir, log_path=log_path, scripts=scripts,
                    python_env=python_env, plink_env=plink_env,
                    ld_prune=ld_prune, account=account,
                ),
            )

            gwf.target_from_template(
                name=f'{sp}_plink_ibd',
                template=plink_ibd_template(
                    sp=sp, bed=bed, bim=bim, fam=fam, prepare_done=prepare_done,
                    out_dir=out_dir, log_path=log_path,
                    plink_env=plink_env, account=account,
                ),
            )

            gwf.target_from_template(
                name=f'{sp}_king_kinship',
                template=king_kinship_template(
                    sp=sp, bed=bed, bim=bim, fam=fam, prepare_done=prepare_done,
                    out_dir=out_dir, log_path=log_path,
                    king_binary=king_binary, account=account,
                ),
            )

            gwf.target_from_template(
                name=f'{sp}_classify',
                template=classify_template(
                    sp=sp, kin0=kin0, genome=genome,
                    king_done=king_done, plink_done=plink_done,
                    out_dir=out_dir, log_path=log_path, scripts=scripts,
                    python_env=python_env, thresholds=thresholds, account=account,
                ),
            )

            gwf.target_from_template(
                name=f'{sp}_plot',
                template=plot_template(
                    sp=sp, classified=classified, classify_done=classify_done,
                    out_dir=out_dir, log_path=log_path, scripts=scripts,
                    python_env=python_env, thresholds=thresholds, account=account,
                ),
            )

        # --- Per-trio Mendelian validation branch (always registered) ---
        gwf.target_from_template(
            name=f'{sp}_trio_filter',
            template=trio_filter_template(
                sp=sp, vcf=vcf, out_dir=out_dir, log_path=log_path, scripts=scripts,
                vcftools_env=vcftools_env, python_env=python_env,
                min_gq=min_gq, trio_max_missing=trio_max_missing, account=account,
            ),
        )

        gwf.target_from_template(
            name=f'{sp}_trio_check',
            template=trio_check_template(
                sp=sp, snp_vcf=trio_snp_vcf, trio_filter_done=trio_filter_done,
                out_dir=out_dir, log_path=log_path, scripts=scripts,
                python_env=python_env, min_gq=trio_min_gq, min_dp=trio_min_dp,
                consistency_min=trio_consistency_min, min_sites=trio_min_sites,
                account=account,
            ),
        )

        gwf.target_from_template(
            name=f'{sp}_trio_plot',
            template=trio_plot_template(
                sp=sp, trio_tsv=trio_tsv, trio_check_done=trio_check_done,
                out_dir=out_dir, log_path=log_path, scripts=scripts,
                python_env=python_env, consistency_min=trio_consistency_min,
                account=account,
            ),
        )

    return gwf


def combined_trio_plot(config_files, gwf):
    """Register a single cross-species combined trio-validation figure.

    Loops over every config, collects each species' trio_mendel TSV and its
    trio_check sentinel, and registers ONE `combined_trio_plot` target (not in
    the per-species loop). BIC family6 is excluded from the supplement.
    """
    species_tsvs = {}
    check_dones = {}
    common = None
    for config_file in config_files:
        with open(config_file) as f:
            config = yaml.safe_load(f)
        output_path = config['output_directory_path']
        log_path = config['log_directory_path']
        common = {
            'combined_dir': f'{output_path}/combined',
            'log_path': log_path,
            'scripts': config['scripts_path'],
            'python_env': config['python_conda_env'],
            'consistency_min': config['trio_consistency_min'],
            'account': config['account'],
        }
        for sp in config['species']:
            species_tsvs[sp] = f'{output_path}/{sp}/{sp}_trio_mendel.tsv'
            check_dones[sp] = f'{log_path}/{sp}_trio_check.DONE'

    if not species_tsvs:
        return gwf

    gwf.target_from_template(
        name='combined_trio_plot',
        template=combined_trio_plot_template(
            species_tsvs=species_tsvs, check_dones=check_dones,
            combined_dir=common['combined_dir'], log_path=common['log_path'],
            scripts=common['scripts'], python_env=common['python_env'],
            consistency_min=common['consistency_min'], exclude=['BIC:family6'],
            account=common['account'],
        ),
    )
    return gwf
