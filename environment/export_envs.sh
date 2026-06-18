#!/usr/bin/env bash
# Export pinned conda environment files for the environments used by this project.
# Run on a machine that has these conda environments installed:
#     bash export_envs.sh
# Produces <env>.yml next to this script for each environment that exists.
set -uo pipefail

ENVS=(
  gwf_new fastqc bwa2 pbmm2 samtools117 bam_check gatk4 bcftools vcftools
  bedtools pysam pyvcf biopython pandas python_phylo paml plink liftoff spider
)

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# List installed envs once.
installed="$(conda env list | awk '{print $1}')"

for env in "${ENVS[@]}"; do
  if grep -qxF "$env" <<< "$installed"; then
    echo "Exporting $env -> $HERE/$env.yml"
    conda env export -n "$env" > "$HERE/$env.yml"
  else
    echo "Skipping $env (not installed)"
  fi
done

echo "Done. Recreate with: conda env create -f <env>.yml"
