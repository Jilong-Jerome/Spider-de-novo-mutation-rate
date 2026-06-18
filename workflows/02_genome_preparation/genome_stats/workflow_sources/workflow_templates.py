#!/bin/env python3
from gwf import AnonymousTarget
import os, glob

def genome_stats_template(work_path, data_path, log_path, species):
    inputs = {"genome": f"{data_path}/{species}.fasta"}
    outputs = {
        "stats": f"{work_path}/{species}/{species}_genome_stats.txt",
        "length_dist": f"{work_path}/{species}/{species}_length_distribution.txt",
        "summary": f"{work_path}/{species}/{species}_genome_summary.txt",
        "log": f"{log_path}/{species}_genome_stats.DONE"
    }
    options = {
        'cores': 2,
        'memory': '16g',
        'walltime': '02:00:00',
        'account': "spider2"
    }
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate seqkit
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    mkdir -p {work_path}/{species}
    mkdir -p {log_path}
    
    # Basic genome statistics using seqkit
    echo "=== BASIC GENOME STATISTICS ===" > {outputs["stats"]}
    echo "Species: {species}" >> {outputs["stats"]}
    echo "Genome file: {inputs["genome"]}" >> {outputs["stats"]}
    echo "Analysis date: $(date)" >> {outputs["stats"]}
    echo "" >> {outputs["stats"]}
    
    # Get basic stats
    seqkit stats {inputs["genome"]} >> {outputs["stats"]}
    echo "" >> {outputs["stats"]}
    
    # Get number of sequences/contigs/chromosomes
    echo "=== SEQUENCE COUNT ===" >> {outputs["stats"]}
    echo "Total sequences: $(grep -c '^>' {inputs["genome"]})" >> {outputs["stats"]}
    echo "" >> {outputs["stats"]}
    
    # Length distribution
    echo "=== LENGTH DISTRIBUTION ===" > {outputs["length_dist"]}
    seqkit fx2tab -nl {inputs["genome"]} | awk '{{print $2}}' | sort -n > {outputs["length_dist"]}
    
    # Summary statistics
    echo "=== GENOME SUMMARY ===" > {outputs["summary"]}
    echo "Species: {species}" >> {outputs["summary"]}
    
    # Calculate N50, N90, etc.
    python3 -c "
import sys
lengths = []
with open('{outputs["length_dist"]}', 'r') as f:
    lengths = [int(line.strip()) for line in f if line.strip()]

lengths.sort(reverse=True)
total_length = sum(lengths)
num_contigs = len(lengths)

# Calculate Nx statistics
def calculate_nx(lengths, x):
    target = total_length * x / 100
    cumsum = 0
    for length in lengths:
        cumsum += length
        if cumsum >= target:
            return length
    return 0

n50 = calculate_nx(lengths, 50)
n90 = calculate_nx(lengths, 90)
n95 = calculate_nx(lengths, 95)

print(f'Total genome length: {{total_length:,}} bp', file=open('{outputs["summary"]}', 'a'))
print(f'Number of contigs/scaffolds: {{num_contigs:,}}', file=open('{outputs["summary"]}', 'a'))
print(f'Longest contig: {{max(lengths):,}} bp', file=open('{outputs["summary"]}', 'a'))
print(f'Shortest contig: {{min(lengths):,}} bp', file=open('{outputs["summary"]}', 'a'))
print(f'Mean contig length: {{int(total_length/num_contigs):,}} bp', file=open('{outputs["summary"]}', 'a'))
print(f'N50: {{n50:,}} bp', file=open('{outputs["summary"]}', 'a'))
print(f'N90: {{n90:,}} bp', file=open('{outputs["summary"]}', 'a'))
print(f'N95: {{n95:,}} bp', file=open('{outputs["summary"]}', 'a'))

# Length distribution bins
bins = [0, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, float('inf')]
bin_labels = ['<1kb', '1-5kb', '5-10kb', '10-50kb', '50-100kb', '100-500kb', '500kb-1Mb', '>1Mb']
bin_counts = [0] * len(bin_labels)

for length in lengths:
    for i, bin_max in enumerate(bins[1:]):
        if length <= bin_max:
            bin_counts[i] += 1
            break

print('\\nLength distribution:', file=open('{outputs["summary"]}', 'a'))
for label, count in zip(bin_labels, bin_counts):
    print(f'  {{label}}: {{count}} contigs', file=open('{outputs["summary"]}', 'a'))
"
    
    # GC content analysis
    echo "" >> {outputs["summary"]}
    echo "=== GC CONTENT ===" >> {outputs["summary"]}
    seqkit fx2tab -g {inputs["genome"]} | awk '{{sum+=$3; count++}} END {{if(count>0) print "Average GC content: " sum/count "%"}}' >> {outputs["summary"]}
    
    echo "Analysis completed: $(date)" >> {outputs["summary"]}
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def genome_comparison_template(work_path, log_path, species_list):
    inputs = {sp: f"{work_path}/{sp}/{sp}_genome_summary.txt" for sp in species_list}
    outputs = {
        "comparison": f"{work_path}/comparative_genome_stats.txt",
        "log": f"{log_path}/genome_comparison.DONE"
    }
    options = {
        'cores': 1,
        'memory': '4g',
        'walltime': '01:00:00',
        'account': "spider2"
    }
    
    input_files = " ".join([f"{work_path}/{sp}/{sp}_genome_summary.txt" for sp in species_list])
    
    spec = f"""
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    
    mkdir -p {log_path}
    
    echo "=== COMPARATIVE GENOME STATISTICS ===" > {outputs["comparison"]}
    echo "Analysis date: $(date)" >> {outputs["comparison"]}
    echo "Species compared: {', '.join(species_list)}" >> {outputs["comparison"]}
    echo "" >> {outputs["comparison"]}
    
    # Create comparison table
    echo "Species\tGenome_Size(bp)\tNum_Contigs\tN50(bp)\tLongest_Contig(bp)\tGC_Content(%)" >> {outputs["comparison"]}
    
    for summary_file in {input_files}; do
        if [ -f "$summary_file" ]; then
            species=$(basename $(dirname $summary_file))
            genome_size=$(grep "Total genome length:" $summary_file | awk '{{print $4}}' | tr -d ',')
            num_contigs=$(grep "Number of contigs/scaffolds:" $summary_file | awk '{{print $4}}' | tr -d ',')
            n50=$(grep "N50:" $summary_file | awk '{{print $2}}' | tr -d ',')
            longest=$(grep "Longest contig:" $summary_file | awk '{{print $3}}' | tr -d ',')
            gc_content=$(grep "Average GC content:" $summary_file | awk '{{print $4}}' | tr -d '%')
            
            echo "$species\t$genome_size\t$num_contigs\t$n50\t$longest\t$gc_content" >> {outputs["comparison"]}
        fi
    done
    
    echo "" >> {outputs["comparison"]}
    echo "Summary completed: $(date)" >> {outputs["comparison"]}
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
