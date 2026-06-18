"""
Command templates for testis verification workflow.

Each function returns an AnonymousTarget with inputs, outputs, options, and spec.
All templates create .DONE marker files for proper GWF dependency tracking.
"""

from pathlib import Path
from gwf import AnonymousTarget
from .workflow_sources import get_conda_activate, CONDA_ENVS, RESOURCES


def extract_human_testis_proteins(proteome: Path, gene_list: Path,
                                   output_dir: Path, log_dir: Path):
    """Extract testis-specific proteins from human proteome."""
    conda_cmd = get_conda_activate(CONDA_ENVS['python'])

    inputs = {
        'proteome': str(proteome),
        'gene_list': str(gene_list),
    }
    outputs = {
        'fasta': f"{output_dir}/human_testis_proteins.fa",
        'stats': f"{output_dir}/extraction_stats.txt",
        'done': f"{log_dir}/extract_human_proteins.DONE",
    }
    options = {**RESOURCES['extract_proteins'], 'account': 'spider2'}

    spec = f"""
{conda_cmd}

mkdir -p {output_dir} {log_dir}

python scripts/extract_human_testis_proteins.py \\
    --proteome {proteome} \\
    --genes {gene_list} \\
    --output {outputs['fasta']} \\
    --stats {outputs['stats']}

echo done > {outputs['done']}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def extract_mim_cds(genome: Path, annotation: Path, output_dir: Path,
                    log_dir: Path, upstream_done: str = None):
    """Extract CDS sequences from MIM genome using AGAT."""
    conda_cmd = get_conda_activate(CONDA_ENVS['agat'])

    inputs = {
        'genome': str(genome),
        'annotation': str(annotation),
    }
    if upstream_done:
        inputs['upstream_done'] = upstream_done

    outputs = {
        'cds': f"{output_dir}/MIM_cds.fa",
        'done': f"{log_dir}/extract_mim_cds.DONE",
    }
    options = {**RESOURCES['agat'], 'account': 'spider2'}

    spec = f"""
{conda_cmd}

mkdir -p {output_dir} {log_dir}

agat_sp_extract_sequences.pl \\
    -g {annotation} \\
    -f {genome} \\
    -t cds \\
    -o {outputs['cds']}

echo done > {outputs['done']}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def make_blast_db(input_fasta: str, output_dir: Path, log_dir: Path,
                  dbtype: str = 'nucl', upstream_done: str = None):
    """Create BLAST database."""
    conda_cmd = get_conda_activate(CONDA_ENVS['blast'])
    db_prefix = f"{output_dir}/MIM_cds_db"

    inputs = {'fasta': input_fasta}
    if upstream_done:
        inputs['upstream_done'] = upstream_done

    outputs = {
        'nhr': f"{db_prefix}.nhr",
        'nin': f"{db_prefix}.nin",
        'nsq': f"{db_prefix}.nsq",
        'done': f"{log_dir}/make_blast_db.DONE",
    }
    options = {**RESOURCES['makeblastdb'], 'account': 'spider2'}

    spec = f"""
{conda_cmd}

mkdir -p {output_dir} {log_dir}

makeblastdb -in {input_fasta} -dbtype {dbtype} -out {db_prefix} -parse_seqids

echo done > {outputs['done']}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def run_tblastn(query: str, db_prefix: Path, output_dir: Path, log_dir: Path,
                evalue: float = 1e-5, max_targets: int = 5,
                threads: int = 8, upstream_done: str = None):
    """Run tblastn search."""
    conda_cmd = get_conda_activate(CONDA_ENVS['blast'])
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

    inputs = {
        'query': query,
        'db_nhr': f"{db_prefix}.nhr",
        'db_nin': f"{db_prefix}.nin",
        'db_nsq': f"{db_prefix}.nsq",
    }
    if upstream_done:
        inputs['upstream_done'] = upstream_done

    outputs = {
        'results': f"{output_dir}/blast_results.tsv",
        'done': f"{log_dir}/run_tblastn.DONE",
    }
    options = {**RESOURCES['tblastn'], 'account': 'spider2'}

    spec = f"""
{conda_cmd}

mkdir -p {output_dir} {log_dir}

tblastn -query {query} -db {db_prefix} -out {outputs['results']} \\
    -evalue {evalue} -max_target_seqs {max_targets} \\
    -num_threads {threads} \\
    -outfmt "{outfmt}"

echo done > {outputs['done']}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def filter_blast_hits(blast_file: str, output_dir: Path, log_dir: Path,
                      min_evalue: float = 1e-5,
                      min_identity: float = 30.0,
                      min_coverage: float = 50.0,
                      upstream_done: str = None):
    """Filter BLAST results and create homolog mapping."""
    conda_cmd = get_conda_activate(CONDA_ENVS['python'])

    inputs = {'blast': blast_file}
    if upstream_done:
        inputs['upstream_done'] = upstream_done

    outputs = {
        'mapping': f"{output_dir}/homolog_mapping.tsv",
        'genes': f"{output_dir}/mim_testis_homologs.txt",
        'summary': f"{output_dir}/blast_summary.txt",
        'done': f"{log_dir}/filter_blast_hits.DONE",
    }
    options = {**RESOURCES['filter_blast'], 'account': 'spider2'}

    spec = f"""
{conda_cmd}

mkdir -p {output_dir} {log_dir}

python scripts/filter_blast_hits.py \\
    --blast {blast_file} \\
    --output-dir {output_dir} \\
    --min-evalue {min_evalue} \\
    --min-identity {min_identity} \\
    --min-coverage {min_coverage}

echo done > {outputs['done']}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def star_genome_generate(genome: Path, annotation: Path, index_dir: Path,
                         log_dir: Path, threads: int = 8, upstream_done: str = None):
    """Build STAR genome index.

    Note: Annotation is GFF3 format with 'Parent' attribute (not 'transcript_id').
    """
    conda_cmd = get_conda_activate(CONDA_ENVS['star'])

    inputs = {
        'genome': str(genome),
        'annotation': str(annotation),
    }
    if upstream_done:
        inputs['upstream_done'] = upstream_done

    outputs = {
        'sa': f"{index_dir}/SA",
        'sa_index': f"{index_dir}/SAindex",
        'genome': f"{index_dir}/Genome",
        'genome_params': f"{index_dir}/genomeParameters.txt",
        'done': f"{log_dir}/star_build_index.DONE",
    }
    options = {**RESOURCES['star_build'], 'account': 'spider2'}

    spec = f"""
{conda_cmd}

mkdir -p {index_dir} {log_dir}

STAR --runMode genomeGenerate \\
    --genomeDir {index_dir} \\
    --genomeFastaFiles {genome} \\
    --sjdbGTFfile {annotation} \\
    --sjdbGTFfeatureExon exon \\
    --sjdbGTFtagExonParentTranscript Parent \\
    --runThreadN {threads}

echo done > {outputs['done']}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def star_align(index_dir: Path, r1: Path, r2: Path, output_dir: Path,
               log_dir: Path, sample_name: str, threads: int = 8,
               upstream_done: str = None):
    """Align reads with STAR.

    Outputs coordinate-sorted BAM directly.
    """
    conda_cmd = get_conda_activate(CONDA_ENVS['star'])
    output_prefix = f"{output_dir}/{sample_name}."

    inputs = {
        'r1': str(r1),
        'r2': str(r2),
        'index_sa': f"{index_dir}/SA",
        'index_genome': f"{index_dir}/Genome",
    }
    if upstream_done:
        inputs['upstream_done'] = upstream_done

    outputs = {
        'bam': f"{output_prefix}Aligned.sortedByCoord.out.bam",
        'log_final': f"{output_prefix}Log.final.out",
        'done': f"{log_dir}/star_align_{sample_name}.DONE",
    }
    options = {**RESOURCES['star_align'], 'account': 'spider2'}

    spec = f"""
{conda_cmd}

mkdir -p {output_dir} {log_dir}

STAR --runMode alignReads \\
    --genomeDir {index_dir} \\
    --readFilesIn {r1} {r2} \\
    --readFilesCommand zcat \\
    --outSAMtype BAM SortedByCoordinate \\
    --outFileNamePrefix {output_prefix} \\
    --runThreadN {threads} \\
    --outFilterMultimapNmax 1 \\
    --outFilterMismatchNmax 2

echo done > {outputs['done']}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def samtools_filter(star_bam: str, output_dir: Path, log_dir: Path,
                    sample_name: str, min_mapq: int = 10,
                    threads: int = 4, upstream_done: str = None):
    """Filter BAM for proper pairs and primary alignments, then index.

    STAR produces coordinate-sorted BAM, so no sorting needed.

    Filters:
    - Require proper pairs (-f 2)
    - Remove unmapped reads (-F 4)
    - Remove secondary alignments (-F 256)
    - Remove supplementary alignments (-F 2048)
    - Combined: -F 2308
    - Require minimum mapping quality (-q)
    """
    conda_cmd = get_conda_activate(CONDA_ENVS['samtools'])

    inputs = {'bam': star_bam}
    if upstream_done:
        inputs['upstream_done'] = upstream_done

    outputs = {
        'filtered_bam': f"{output_dir}/{sample_name}.filtered.bam",
        'bai': f"{output_dir}/{sample_name}.filtered.bam.bai",
        'flagstat': f"{output_dir}/{sample_name}.flagstat.txt",
        'done': f"{log_dir}/samtools_filter_{sample_name}.DONE",
    }
    options = {**RESOURCES['samtools_filter'], 'account': 'spider2'}

    spec = f"""
{conda_cmd}

mkdir -p {output_dir} {log_dir}

# Filter: proper pairs, primary alignments, MAPQ >= {min_mapq}
samtools view -@ {threads} -b -f 2 -F 2308 -q {min_mapq} {star_bam} -o {outputs['filtered_bam']}

# Index the filtered BAM
samtools index {outputs['filtered_bam']}

# Generate alignment statistics
samtools flagstat {outputs['filtered_bam']} > {outputs['flagstat']}

echo done > {outputs['done']}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def feature_counts(bam_file: str, annotation: Path, output_dir: Path,
                   log_dir: Path, sample_name: str,
                   feature_type: str = 'gene', attribute: str = 'ID',
                   threads: int = 4, upstream_done: str = None):
    """Count features with featureCounts."""
    conda_cmd = get_conda_activate(CONDA_ENVS['featurecounts'])

    inputs = {
        'bam': bam_file,
        'annotation': str(annotation),
    }
    if upstream_done:
        inputs['upstream_done'] = upstream_done

    outputs = {
        'counts': f"{output_dir}/{sample_name}_counts.txt",
        'summary': f"{output_dir}/{sample_name}_counts.txt.summary",
        'done': f"{log_dir}/featurecounts_{sample_name}.DONE",
    }
    options = {**RESOURCES['featurecounts'], 'account': 'spider2'}

    spec = f"""
{conda_cmd}

mkdir -p {output_dir} {log_dir}

featureCounts -T {threads} -p -B -C \\
    -t {feature_type} -g {attribute} \\
    -a {annotation} -o {outputs['counts']} {bam_file}

echo done > {outputs['done']}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def calculate_tpm(count_files: list, sample_names: list, output_dir: Path,
                  log_dir: Path, upstream_dones: list = None):
    """Calculate TPM and merge samples."""
    conda_cmd = get_conda_activate(CONDA_ENVS['python'])

    inputs = {f'counts_{i}': str(f) for i, f in enumerate(count_files)}
    if upstream_dones:
        for i, done in enumerate(upstream_dones):
            inputs[f'upstream_done_{i}'] = done

    outputs = {
        'tpm': f"{output_dir}/tpm_matrix.tsv",
        'counts': f"{output_dir}/count_matrix.tsv",
        'done': f"{log_dir}/calculate_tpm.DONE",
    }
    options = {**RESOURCES['calculate_tpm'], 'account': 'spider2'}

    counts_str = ' '.join(str(f) for f in count_files)
    names_str = ' '.join(sample_names)

    spec = f"""
{conda_cmd}

mkdir -p {output_dir} {log_dir}

python scripts/calculate_tpm.py \\
    --counts {counts_str} \\
    --sample-names {names_str} \\
    --output {outputs['tpm']}

echo done > {outputs['done']}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def compare_expression(tpm_file: str, testis_genes: str,
                       control_samples: list, candidate_samples: list,
                       output_dir: Path, log_dir: Path,
                       upstream_done: str = None):
    """Compare expression and generate report."""
    conda_cmd = get_conda_activate(CONDA_ENVS['python'])

    inputs = {
        'tpm': tpm_file,
        'genes': testis_genes,
    }
    if upstream_done:
        inputs['upstream_done'] = upstream_done

    outputs = {
        'report': f"{output_dir}/testis_verification_report.txt",
        'heatmap': f"{output_dir}/expression_heatmap.png",
        'barplot': f"{output_dir}/expression_barplot.png",
        'done': f"{log_dir}/compare_expression.DONE",
    }
    options = {**RESOURCES['compare_expression'], 'account': 'spider2'}

    control_str = ' '.join(control_samples)
    candidate_str = ' '.join(candidate_samples)

    spec = f"""
{conda_cmd}

mkdir -p {output_dir} {log_dir}

python scripts/compare_expression.py \\
    --tpm {tpm_file} \\
    --testis-genes {testis_genes} \\
    --control-samples {control_str} \\
    --candidate-samples {candidate_str} \\
    --output-dir {output_dir}

echo done > {outputs['done']}
"""
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
