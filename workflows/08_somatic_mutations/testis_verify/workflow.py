"""
Testis Verification GWF Workflow

Verify whether RNA-seq samples are testes by comparing expression
of testis-specific gene homologs between ovary (control) and candidate samples.

Usage:
    gwf status
    gwf run
"""

from pathlib import Path

from gwf import Workflow

from workflow_sources import (
    extract_human_testis_proteins,
    extract_mim_cds,
    make_blast_db,
    run_tblastn,
    filter_blast_hits,
    star_genome_generate,
    star_align,
    samtools_filter,
    feature_counts,
    calculate_tpm,
    compare_expression,
)

# Initialize workflow
gwf = Workflow(defaults={'account': 'spider2'})

# ============================================================================
# Configuration
# ============================================================================

# Data paths
DATA_DIR = Path('data')
GENE_LIST = DATA_DIR / 'gene' / 'human_testis_specific_genes.txt'
HUMAN_PROTEOME = DATA_DIR / 'gene' / 'uniprotkb_proteome_UP000005640_AND_revi_2026_01_19.fasta'
GENOME = DATA_DIR / 'genome' / 'MIM' / 'MIM_ncbi_chromosome.fa'
ANNOTATION = DATA_DIR / 'genome' / 'MIM' / 'MIM_gene.gff3'

# Sample definitions
SAMPLES = {
    'MIM_ovary_5': {
        'r1': DATA_DIR / 'sample' / 'MIM_ovary_5_R1.fastq.gz',
        'r2': DATA_DIR / 'sample' / 'MIM_ovary_5_R2.fastq.gz',
        'type': 'control',
    },
    'MIM_testis_candidate_red': {
        'r1': DATA_DIR / 'sample' / 'MIM_testis_candidate_red_R1.fastq.gz',
        'r2': DATA_DIR / 'sample' / 'MIM_testis_candidate_red_R2.fastq.gz',
        'type': 'candidate',
    },
    'MIM_testis_candidate_white': {
        'r1': DATA_DIR / 'sample' / 'MIM_testis_candidate_white_R1.fastq.gz',
        'r2': DATA_DIR / 'sample' / 'MIM_testis_candidate_white_R2.fastq.gz',
        'type': 'candidate',
    },
}

# Output directories
STEPS_DIR = Path('steps')
HUMAN_PROTEINS_DIR = STEPS_DIR / '01_human_proteins'
MIM_CDS_DIR = STEPS_DIR / '02_mim_cds'
BLAST_DIR = STEPS_DIR / '03_blast'
ALIGNMENT_DIR = STEPS_DIR / '04_alignment'
STAR_INDEX_DIR = ALIGNMENT_DIR / 'STAR_index'
COUNTS_DIR = STEPS_DIR / '05_counts'
ANALYSIS_DIR = STEPS_DIR / '06_analysis'
LOGS_DIR = Path('logs')

# Parameters
BLAST_EVALUE = 1e-5
BLAST_MAX_TARGETS = 5
BLAST_MIN_IDENTITY = 30.0
BLAST_MIN_COVERAGE = 50.0
STAR_THREADS = 8
FEATURECOUNTS_THREADS = 4

# ============================================================================
# Step 1: Extract Human Testis Gene Protein Sequences
# ============================================================================

step1 = gwf.target_from_template(
    name='extract_human_proteins',
    template=extract_human_testis_proteins(
        proteome=HUMAN_PROTEOME,
        gene_list=GENE_LIST,
        output_dir=HUMAN_PROTEINS_DIR,
        log_dir=LOGS_DIR,
    ),
)

# ============================================================================
# Step 2: Extract MIM CDS from Genome
# ============================================================================

step2 = gwf.target_from_template(
    name='extract_mim_cds',
    template=extract_mim_cds(
        genome=GENOME,
        annotation=ANNOTATION,
        output_dir=MIM_CDS_DIR,
        log_dir=LOGS_DIR,
    ),
)

# ============================================================================
# Step 3: Create BLAST Database
# ============================================================================

step3 = gwf.target_from_template(
    name='make_blast_db',
    template=make_blast_db(
        input_fasta=step2.outputs['cds'],
        output_dir=BLAST_DIR,
        log_dir=LOGS_DIR,
        dbtype='nucl',
        upstream_done=step2.outputs['done'],
    ),
)

# ============================================================================
# Step 4: BLAST Human Proteins vs MIM CDS (tblastn)
# ============================================================================

step4 = gwf.target_from_template(
    name='run_tblastn',
    template=run_tblastn(
        query=step1.outputs['fasta'],
        db_prefix=BLAST_DIR / 'MIM_cds_db',
        output_dir=BLAST_DIR,
        log_dir=LOGS_DIR,
        evalue=BLAST_EVALUE,
        max_targets=BLAST_MAX_TARGETS,
        threads=8,
        upstream_done=step3.outputs['done'],
    ),
)

# ============================================================================
# Step 5: Filter BLAST and Create Homolog Mapping
# ============================================================================

step5 = gwf.target_from_template(
    name='filter_blast_hits',
    template=filter_blast_hits(
        blast_file=step4.outputs['results'],
        output_dir=BLAST_DIR,
        log_dir=LOGS_DIR,
        min_evalue=BLAST_EVALUE,
        min_identity=BLAST_MIN_IDENTITY,
        min_coverage=BLAST_MIN_COVERAGE,
        upstream_done=step4.outputs['done'],
    ),
)

# ============================================================================
# Step 6: Build STAR Index
# ============================================================================

step6 = gwf.target_from_template(
    name='star_build_index',
    template=star_genome_generate(
        genome=GENOME,
        annotation=ANNOTATION,
        index_dir=STAR_INDEX_DIR,
        log_dir=LOGS_DIR,
        threads=STAR_THREADS,
    ),
)

# ============================================================================
# Steps 7-9: Per-sample alignment, filtering, and counting
# ============================================================================

sample_targets = {}
featurecounts_dones = []

for sample_name, sample_info in SAMPLES.items():
    # Step 7: STAR alignment
    align_target = gwf.target_from_template(
        name=f'star_align_{sample_name}',
        template=star_align(
            index_dir=STAR_INDEX_DIR,
            r1=sample_info['r1'],
            r2=sample_info['r2'],
            output_dir=ALIGNMENT_DIR,
            log_dir=LOGS_DIR,
            sample_name=sample_name,
            threads=STAR_THREADS,
            upstream_done=step6.outputs['done'],
        ),
    )

    # Step 8: Filter and index BAM
    filter_target = gwf.target_from_template(
        name=f'samtools_filter_{sample_name}',
        template=samtools_filter(
            star_bam=align_target.outputs['bam'],
            output_dir=ALIGNMENT_DIR,
            log_dir=LOGS_DIR,
            sample_name=sample_name,
            min_mapq=10,
            threads=4,
            upstream_done=align_target.outputs['done'],
        ),
    )

    # Step 9: Feature counting
    counts_target = gwf.target_from_template(
        name=f'featurecounts_{sample_name}',
        template=feature_counts(
            bam_file=filter_target.outputs['filtered_bam'],
            annotation=ANNOTATION,
            output_dir=COUNTS_DIR,
            log_dir=LOGS_DIR,
            sample_name=sample_name,
            feature_type='gene',
            attribute='ID',
            threads=FEATURECOUNTS_THREADS,
            upstream_done=filter_target.outputs['done'],
        ),
    )

    sample_targets[sample_name] = {
        'align': align_target,
        'filter': filter_target,
        'counts': counts_target,
    }
    featurecounts_dones.append(counts_target.outputs['done'])

# ============================================================================
# Step 10: Calculate TPM and Merge
# ============================================================================

sample_names = list(SAMPLES.keys())
count_files = [sample_targets[s]['counts'].outputs['counts'] for s in sample_names]

step10 = gwf.target_from_template(
    name='calculate_tpm',
    template=calculate_tpm(
        count_files=count_files,
        sample_names=sample_names,
        output_dir=ANALYSIS_DIR,
        log_dir=LOGS_DIR,
        upstream_dones=featurecounts_dones,
    ),
)

# ============================================================================
# Step 11: Compare Expression & Generate Report
# ============================================================================

control_samples = [name for name, info in SAMPLES.items() if info['type'] == 'control']
candidate_samples = [name for name, info in SAMPLES.items() if info['type'] == 'candidate']

step11 = gwf.target_from_template(
    name='compare_expression',
    template=compare_expression(
        tpm_file=step10.outputs['tpm'],
        testis_genes=step5.outputs['genes'],
        control_samples=control_samples,
        candidate_samples=candidate_samples,
        output_dir=ANALYSIS_DIR,
        log_dir=LOGS_DIR,
        upstream_done=step10.outputs['done'],
    ),
)
