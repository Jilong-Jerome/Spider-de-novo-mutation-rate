#!/bin/bash
# concat_filter_paml.sh - Runtime wrapper for concat+filter+phylip+codeml+parse pipeline
#
# This script reads the gene ID list at runtime (not workflow definition time),
# which is necessary for bootstrap replicates whose ID files are generated
# by earlier jobs in the DAG.
#
# Runs both branch model (model=1) and pairwise model (runmode=-2) in separate subfolders.
#
# Usage:
#   bash concat_filter_paml.sh <replicate_id> <id_list_file> <fasta_dir> <work_dir> <script_path> <trio_tree>
#
# Args:
#   replicate_id  - e.g. auto_all, auto_bs_1, auto_bs_2, ...
#   id_list_file  - path to file with gene IDs (one per line)
#   fasta_dir     - directory containing per-gene FASTA files ({gene_id}.fa)
#   work_dir      - working directory for this replicate
#   script_path   - path to workflow_sources/ with helper scripts
#   trio_tree     - Newick tree string, e.g. "((SAR,PAC),BIC)"

set -euo pipefail

REP_ID=$1
ID_LIST=$2
FASTA_DIR=$3
WORK_DIR=$4
SCRIPT_PATH=$5
TRIO_TREE=$6

echo "START: $(date)"
echo "Replicate: ${REP_ID}"
echo "ID list: ${ID_LIST}"
echo "FASTA dir: ${FASTA_DIR}"
echo "Work dir: ${WORK_DIR}"
echo "Script path: ${SCRIPT_PATH}"
echo "Trio tree: ${TRIO_TREE}"

mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"

CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# ==============================================================
#  Shared steps: concat + codon filter + phylip conversion
# ==============================================================

# ---- Step 1: Concatenate all gene alignments into supermatrix ----
echo "Concatenating gene alignments..."
conda activate python_phylo
python3 "${SCRIPT_PATH}/concat_genes.py" "${ID_LIST}" "${FASTA_DIR}" "${REP_ID}.fa.cat"
echo "Concatenation done: $(wc -c < "${REP_ID}.fa.cat") bytes"

# ---- Step 2: Keep fully-callable codons + convert to PAML phylip ----
# The 3-species alignment is positionally exact (per-species CDS called against
# a shared reference after alignment QC), so the only filtering needed is to
# drop codons not callable in all species.
echo "Filtering to fully-callable codons..."
python3 "${SCRIPT_PATH}/concat_codon_filter_paml.py" \
    "${REP_ID}.fa.cat" "${REP_ID}.fa.filtered.concat.paml.phy"
echo "Codon filtering done"

PAML_PHY="$(pwd)/${REP_ID}.fa.filtered.concat.paml.phy"

# ==============================================================
#  Branch model (model=1, free-ratio) in branch/ subfolder
# ==============================================================
echo "===== Branch model ====="
mkdir -p branch
cd branch
cp "${SCRIPT_PATH}/codeml_ofree.ctl" codeml.ctl
echo "${TRIO_TREE};" > tree.txt
ln -f -s "${PAML_PHY}" seq.txt

echo "Running codeml (branch model)..."
conda activate paml
codeml
echo "codeml branch done"

echo "Parsing branch results..."
conda activate ete3
python3 "${SCRIPT_PATH}/paml2tab.py" results.txt "${REP_ID}.tab" "${REP_ID}"
python3 "${SCRIPT_PATH}/branch_dnds.py" "${REP_ID}.tab" "${REP_ID}_branch.tab"
echo "Branch parsing done"

cd "${WORK_DIR}"

# ==============================================================
#  Pairwise model (runmode=-2) in pairwise/ subfolder
# ==============================================================
echo "===== Pairwise model ====="
mkdir -p pairwise
cd pairwise
cp "${SCRIPT_PATH}/codeml_pairwise.ctl" codeml.ctl
echo "${TRIO_TREE};" > tree.txt
ln -f -s "${PAML_PHY}" seq.txt

echo "Running codeml (pairwise model)..."
conda activate paml
codeml
echo "codeml pairwise done"

echo "Parsing pairwise results..."
conda activate python_phylo
python3 "${SCRIPT_PATH}/parse_pairwise.py" "$(pwd)" "${REP_ID}_pairwise.tab" "${REP_ID}"
echo "Pairwise parsing done"

cd "${WORK_DIR}"

echo "FINISH: $(date)"
