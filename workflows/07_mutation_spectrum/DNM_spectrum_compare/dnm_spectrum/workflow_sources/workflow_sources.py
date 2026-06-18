#!/usr/bin/env python3
"""
workflow_sources.py
Build GWF targets for the DNM spectrum analysis from a config file.
"""
import yaml
from workflow_templates import (
    index_genome_template,
    extract_context_template,
    compute_somatic_spectrum_template,
    compute_spectrum_template,
    plot_spectrum_template,
    plot_comparison_template,
    compare_groups_test_template,
    deviation_interaction_test_template,
    pairwise_somatic_test_template,
    baseline_somatic_test_template,
    spectrum_phylo_mantel_template,
)


def dnm_spectrum_workflow(config_file, gwf):
    """Register all targets for the DNM spectrum workflow."""
    with open(config_file) as f:
        config = yaml.safe_load(f)

    account = config['account']
    output_path = config['output_directory_path']
    log_path = config['log_directory_path']
    scripts = config['scripts_path']
    samtools_env = config['samtools_conda_env']
    python_env = config['python_conda_env']

    species_configs = config['species']
    spectrum_groups = config['spectrum_groups']
    somatic_spectra = config.get('somatic_spectra', {})
    social_subsocial_pairs = config.get('social_subsocial_pairs', [])
    species_tree_newick = config.get('species_tree_newick')

    # Step 1 & 2: Per-species index genome and extract context
    for SP, sp_conf in species_configs.items():
        genome_fa = sp_conf['genome']
        dnm_file = sp_conf['dnm_file']
        x_chroms = sp_conf.get('x_chromosomes', [])
        exclude = sp_conf.get('exclude_trios', [])

        # Step 1: Index genome
        gwf.target_from_template(
            name=f'{SP}_index_genome',
            template=index_genome_template(
                genome_fa=genome_fa,
                SP=SP,
                log_path=log_path,
                samtools_conda_env=samtools_env,
                account=account,
            ),
        )

        # Step 2: Extract trinucleotide context
        gwf.target_from_template(
            name=f'{SP}_extract_context',
            template=extract_context_template(
                dnm_file=dnm_file,
                genome_fa=genome_fa,
                SP=SP,
                output_path=output_path,
                log_path=log_path,
                scripts=scripts,
                x_chromosomes=x_chroms,
                exclude_trios=exclude,
                samtools_conda_env=samtools_env,
                python_conda_env=python_env,
                account=account,
            ),
        )

        # Step 2b: Build the per-species 96-cat somatic spectrum locally
        # from the raw oriented_ref/oriented_alt/oriented_context table.
        somatic_dnm_file = sp_conf.get('somatic_dnm_file')
        if somatic_dnm_file:
            gwf.target_from_template(
                name=f'{SP}_somatic_spectrum',
                template=compute_somatic_spectrum_template(
                    somatic_dnm_file=somatic_dnm_file,
                    SP=SP,
                    output_path=output_path,
                    log_path=log_path,
                    scripts=scripts,
                    python_conda_env=python_env,
                    account=account,
                ),
            )

    # Step 3: Compute spectrum per group per chromosome group
    for group_name, group_conf in spectrum_groups.items():
        group_species = group_conf['species']

        # Collect annotated DNM files and their DONE sentinels
        input_files = []
        done_files = []
        for SP in group_species:
            annotated = f'{output_path}/annotated_dnm/{SP}/{SP}_dnm_annotated.tsv'
            done = f'{log_path}/{SP}_extract_context.DONE'
            input_files.append(annotated)
            done_files.append(done)

        for chrom_group in ['autosome', 'x_chromosome']:
            spectrum_tsv = f'{output_path}/spectrum/{group_name}/{group_name}_{chrom_group}_spectrum.tsv'

            gwf.target_from_template(
                name=f'{group_name}_{chrom_group}_spectrum',
                template=compute_spectrum_template(
                    input_files=input_files,
                    chrom_group=chrom_group,
                    group_name=group_name,
                    output_path=output_path,
                    log_path=log_path,
                    scripts=scripts,
                    dependency_done_files=done_files,
                    python_conda_env=python_env,
                    account=account,
                ),
            )

            gwf.target_from_template(
                name=f'{group_name}_{chrom_group}_plot',
                template=plot_spectrum_template(
                    spectrum_tsv=spectrum_tsv,
                    group_name=group_name,
                    chrom_group=chrom_group,
                    output_path=output_path,
                    log_path=log_path,
                    scripts=scripts,
                    python_conda_env=python_env,
                    account=account,
                ),
            )

        # Germline (autosome) vs merged-somatic comparison plot
        somatic_tsvs = [somatic_spectra[SP] for SP in group_species if SP in somatic_spectra]

        if somatic_tsvs:
            germline_autosome_tsv = f'{output_path}/spectrum/{group_name}/{group_name}_autosome_spectrum.tsv'
            gwf.target_from_template(
                name=f'{group_name}_comparison_plot',
                template=plot_comparison_template(
                    germline_spectrum_tsv=germline_autosome_tsv,
                    somatic_tsvs=somatic_tsvs,
                    group_name=group_name,
                    output_path=output_path,
                    log_path=log_path,
                    scripts=scripts,
                    python_conda_env=python_env,
                    account=account,
                ),
            )

    # Social vs subsocial 7-category statistical test (germline + somatic)
    if 'social' in spectrum_groups and 'subsocial' in spectrum_groups:
        soc_species = spectrum_groups['social']['species']
        sub_species = spectrum_groups['subsocial']['species']
        soc_somatic = [somatic_spectra[SP] for SP in soc_species if SP in somatic_spectra]
        sub_somatic = [somatic_spectra[SP] for SP in sub_species if SP in somatic_spectra]

        germline_social_tsv = f'{output_path}/spectrum/social/social_autosome_spectrum.tsv'
        germline_subsocial_tsv = f'{output_path}/spectrum/subsocial/subsocial_autosome_spectrum.tsv'

        if soc_somatic and sub_somatic:
            gwf.target_from_template(
                name='social_vs_subsocial_7category_test',
                template=compare_groups_test_template(
                    germline_social_tsv=germline_social_tsv,
                    germline_subsocial_tsv=germline_subsocial_tsv,
                    somatic_social_tsvs=soc_somatic,
                    somatic_subsocial_tsvs=sub_somatic,
                    output_path=output_path,
                    log_path=log_path,
                    scripts=scripts,
                    python_conda_env=python_env,
                    account=account,
                ),
            )

            # Pairwise somatic test across closest social-subsocial pairs
            pairs = []
            for pair in social_subsocial_pairs:
                soc_sp, sub_sp = pair[0], pair[1]
                if soc_sp in somatic_spectra and sub_sp in somatic_spectra:
                    pairs.append({
                        'name': f'{soc_sp}_vs_{sub_sp}',
                        'social_path': somatic_spectra[soc_sp],
                        'subsocial_path': somatic_spectra[sub_sp],
                    })
            if pairs:
                gwf.target_from_template(
                    name='pairwise_somatic_test',
                    template=pairwise_somatic_test_template(
                        pairs=pairs,
                        output_path=output_path,
                        log_path=log_path,
                        scripts=scripts,
                        python_conda_env=python_env,
                        account=account,
                    ),
                )

            # Per-species somatic shift vs merged subsocial baseline
            social_species_paths = [(sp, somatic_spectra[sp])
                                    for sp in soc_species if sp in somatic_spectra]
            subsocial_species_paths = [(sp, somatic_spectra[sp])
                                       for sp in sub_species if sp in somatic_spectra]
            if social_species_paths and len(subsocial_species_paths) >= 2:
                gwf.target_from_template(
                    name='baseline_somatic_test',
                    template=baseline_somatic_test_template(
                        social_species=social_species_paths,
                        subsocial_species=subsocial_species_paths,
                        output_path=output_path,
                        log_path=log_path,
                        scripts=scripts,
                        python_conda_env=python_env,
                        account=account,
                    ),
                )

            gwf.target_from_template(
                name='social_vs_subsocial_deviation_interaction_test',
                template=deviation_interaction_test_template(
                    germline_social_tsv=germline_social_tsv,
                    germline_subsocial_tsv=germline_subsocial_tsv,
                    somatic_social_tsvs=soc_somatic,
                    somatic_subsocial_tsvs=sub_somatic,
                    output_path=output_path,
                    log_path=log_path,
                    scripts=scripts,
                    python_conda_env=python_env,
                    account=account,
                ),
            )

    # Spectrum-phylogeny Mantel test (global, uses all species with somatic spectra)
    if species_tree_newick and somatic_spectra:
        mantel_species = [(sp, somatic_spectra[sp]) for sp in sorted(somatic_spectra)]
        if len(mantel_species) >= 3:
            gwf.target_from_template(
                name='spectrum_phylo_mantel',
                template=spectrum_phylo_mantel_template(
                    species_paths=mantel_species,
                    tree_newick=species_tree_newick,
                    output_path=output_path,
                    log_path=log_path,
                    scripts=scripts,
                    python_conda_env=python_env,
                    account=account,
                ),
            )

    return gwf
