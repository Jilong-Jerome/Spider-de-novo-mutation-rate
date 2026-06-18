#!/usr/bin/env python3
"""
workflow_sources.py
Build GWF targets for the DNM mutation-rate analysis.
Reads the species/group blocks from the spectrum.config.yaml referenced via
`species_config:` in rate.config.yaml, so configuration stays in one place.
"""
import os
import yaml
from workflow_templates import (
    callable_vcf_trinuc_template,
    aggregate_callable_trinuc_template,
    germline_rate_summary_template,
    somatic_rate_summary_template,
    aggregate_rates_template,
    plot_rates_template,
    phylo_aligned_rate_plot_template,
)


def _load_offspring_chrom_pairs(callable_file, x_chromosomes, exclude_trios):
    """Return ordered list of (offspring_id, chrom) pairs from the species
    callable TSV after filtering out X chromosomes and excluded trios."""
    x_set = set(x_chromosomes or [])
    excl_set = set(exclude_trios or [])
    seen = set()
    pairs = []
    with open(callable_file) as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 3:
                continue
            chrom = fields[0]
            tag = fields[2]
            # tag is "{offspring_id}_{chrom}" with offspring_id ending in "_offspring"
            suffix = f'_{chrom}'
            if not tag.endswith(suffix):
                # fallback: split off the trailing chromosome token
                continue
            offspring_id = tag[:-len(suffix)]
            if chrom in x_set:
                continue
            if offspring_id in excl_set:
                continue
            key = (offspring_id, chrom)
            if key in seen:
                continue
            seen.add(key)
            pairs.append(key)
    return pairs


def _vcf_path(vcf_root, SP, offspring_id, chrom):
    return (
        f'{vcf_root}/{SP}/callable/minDP_26/{offspring_id}/chrom_split/'
        f'{offspring_id}_{chrom}_callable_sites.vcf'
    )


def dnm_rate_workflow(config_file, gwf):
    with open(config_file) as f:
        rate_cfg = yaml.safe_load(f)

    spectrum_cfg_path = rate_cfg['species_config']
    with open(spectrum_cfg_path) as f:
        spectrum_cfg = yaml.safe_load(f)

    account = rate_cfg['account']
    output_path = rate_cfg['output_directory_path']
    log_path = rate_cfg['log_directory_path']
    scripts = rate_cfg['scripts_path']
    samtools_env = rate_cfg['samtools_conda_env']
    python_env = rate_cfg['python_conda_env']

    vcf_root = spectrum_cfg['germline_callable_vcf_root']
    species_configs = spectrum_cfg['species']
    spectrum_groups = spectrum_cfg['spectrum_groups']

    species_list = list(species_configs.keys())

    for SP in species_list:
        sp_conf = species_configs[SP]
        genome_fa = sp_conf['genome']
        dnm_file = sp_conf['dnm_file']
        somatic_dnm_file = sp_conf.get('somatic_dnm_file')
        callable_file = sp_conf.get('germline_callable_file')
        somatic_trinuc_counts_file = sp_conf.get('somatic_trinuc_counts_file')
        x_chroms = sp_conf.get('x_chromosomes', [])
        exclude = sp_conf.get('exclude_trios', [])

        # Step 1: per-(offspring, chrom) callable trinuc extraction.
        # Skip species without callable inputs.
        if not callable_file or not os.path.exists(callable_file):
            continue

        offspring_chrom_pairs = _load_offspring_chrom_pairs(
            callable_file, x_chroms, exclude
        )

        for offspring_id, chrom in offspring_chrom_pairs:
            vcf_path = _vcf_path(vcf_root, SP, offspring_id, chrom)
            target_name = f'{SP}_{offspring_id}_{chrom}_callable_trinuc'
            gwf.target_from_template(
                name=target_name,
                template=callable_vcf_trinuc_template(
                    vcf_file=vcf_path,
                    genome_fa=genome_fa,
                    SP=SP,
                    offspring_id=offspring_id,
                    chrom=chrom,
                    output_path=output_path,
                    log_path=log_path,
                    scripts=scripts,
                    samtools_conda_env=samtools_env,
                    account=account,
                ),
            )

        # Step 2: per-species aggregator
        if offspring_chrom_pairs:
            gwf.target_from_template(
                name=f'{SP}_germline_callable_trinuc',
                template=aggregate_callable_trinuc_template(
                    SP=SP,
                    offspring_chrom_pairs=offspring_chrom_pairs,
                    output_path=output_path,
                    log_path=log_path,
                    scripts=scripts,
                    python_conda_env=python_env,
                    account=account,
                ),
            )

            # Step 3: per-species germline rate summary
            callable_trinuc_tsv = (
                f'{output_path}/rate/germline_callable_trinuc/{SP}/'
                f'{SP}_autosome_callable_trinuc.tsv'
            )
            gwf.target_from_template(
                name=f'{SP}_germline_rate_summary',
                template=germline_rate_summary_template(
                    SP=SP,
                    dnm_file=dnm_file,
                    callable_file=callable_file,
                    callable_trinuc_tsv=callable_trinuc_tsv,
                    x_chromosomes=x_chroms,
                    exclude_trios=exclude,
                    output_path=output_path,
                    log_path=log_path,
                    scripts=scripts,
                    python_conda_env=python_env,
                    account=account,
                ),
            )

        # Step 4: per-species somatic rate summary
        if somatic_dnm_file and somatic_trinuc_counts_file:
            gwf.target_from_template(
                name=f'{SP}_somatic_rate_summary',
                template=somatic_rate_summary_template(
                    SP=SP,
                    somatic_dnm_file=somatic_dnm_file,
                    somatic_trinuc_counts_file=somatic_trinuc_counts_file,
                    output_path=output_path,
                    log_path=log_path,
                    scripts=scripts,
                    python_conda_env=python_env,
                    account=account,
                ),
            )

    # Step 5: master aggregator
    gwf.target_from_template(
        name='mutation_rate_aggregate',
        template=aggregate_rates_template(
            species_list=species_list,
            spectrum_config_path=spectrum_cfg_path,
            output_path=output_path,
            log_path=log_path,
            scripts=scripts,
            python_conda_env=python_env,
            account=account,
        ),
    )

    # Step 6: visualization
    gwf.target_from_template(
        name='mutation_rate_plot',
        template=plot_rates_template(
            output_path=output_path,
            log_path=log_path,
            scripts=scripts,
            python_conda_env=python_env,
            account=account,
        ),
    )

    # Step 7: phylogeny-aligned per-species overall rate plot
    gwf.target_from_template(
        name='phylo_aligned_rate_plot',
        template=phylo_aligned_rate_plot_template(
            spectrum_config_path=spectrum_cfg_path,
            output_path=output_path,
            log_path=log_path,
            scripts=scripts,
            python_conda_env=python_env,
            account=account,
        ),
    )

    return gwf
