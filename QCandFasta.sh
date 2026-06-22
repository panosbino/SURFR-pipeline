#!/bin/bash
#SBATCH -A naiss2024-22-713
#SBATCH -p memory
#SBATCH -c 30
#SBATCH -t 20:00:00
#SBATCH -J QCandFasta
#SBATCH --begin=now
#SBATCH --mail-type=ALL

# =============================================================================
# SURFR Pipeline - Step 2: QCandFasta.sh
#
# Runs miRTrace quality control on cancer and adjacent-normal FASTQ files,
# then compresses the resulting uncollapsed FASTA outputs for use in Step 3.
#
# Usage (direct):
#   bash QCandFasta.sh <project_id> <dataset> <projPath> <mirtrace_path>
#
# Usage (SLURM):
#   sbatch QCandFasta.sh <project_id> <dataset> <projPath> <mirtrace_path>
#
# Arguments:
#   project_id    - Project identifier, e.g. LUAD
#   dataset       - Dataset name, e.g. TCGA
#   projPath      - Absolute path to the project root directory
#   mirtrace_path - Full path to the miRTrace executable
#
# Inputs (from Step 1):
#   <outdir>/all_<project_id>_cancer.fastq.gz
#   <outdir>/all_<project_id>_adjacent.fastq.gz
#
# Outputs:
#   <outdir>/all_<project_id>_cancer_miRTrace/qc_passed_reads.all.uncollapsed/all_<project_id>_cancer.fasta.gz
#   <outdir>/all_<project_id>_adjacent_miRTrace/qc_passed_reads.all.uncollapsed/all_<project_id>_adjacent.fasta.gz
# =============================================================================

set -euo pipefail

# ----------------------------
# INPUT ARGUMENTS
# ----------------------------
if [ "$#" -ne 4 ]; then
    echo "ERROR: Expected 4 arguments, got $#"
    echo "Usage: $0 <project_id> <dataset> <projPath> <mirtrace_path>"
    exit 1
fi

project_id=$1       # e.g. LUAD
dataset=$2          # e.g. TCGA
projPath=$3         # e.g. /proj/myproject
mirtrace_path=$4    # e.g. /path/to/mirtrace

# ----------------------------
# PATHS
# ----------------------------
outdir=${projPath}/Data/${project_id}/${dataset}

# Validate miRTrace executable
if [ ! -x "${mirtrace_path}" ]; then
    echo "ERROR: miRTrace executable not found or not executable: ${mirtrace_path}"
    exit 1
fi

# Validate input FASTQ files from Step 1
for condition in cancer adjacent; do
    fastq="${outdir}/all_${project_id}_${condition}.fastq.gz"
    if [ ! -f "${fastq}" ]; then
        echo "ERROR: Input FASTQ not found (did Step 1 complete?): ${fastq}"
        exit 1
    fi
done

# ----------------------------
# LOAD MODULES
# ----------------------------
ml PDC
ml pigz

# ----------------------------
# miRTrace QC — CANCER
# ----------------------------
echo "[$(date)] Running miRTrace QC on cancer samples..."

# FIX: input file must come AFTER all flags in miRTrace syntax
"${mirtrace_path}" qc \
    --species hsa \
    --output-dir "${outdir}/all_${project_id}_cancer_miRTrace" \
    --write-fasta \
    --uncollapse-fasta \
    --num-threads 30 \
    --force \
    "${outdir}/all_${project_id}_cancer.fastq.gz"

# ----------------------------
# miRTrace QC — ADJACENT-NORMAL
# ----------------------------
echo "[$(date)] Running miRTrace QC on adjacent-normal samples..."

"${mirtrace_path}" qc \
    --species hsa \
    --output-dir "${outdir}/all_${project_id}_adjacent_miRTrace" \
    --write-fasta \
    --uncollapse-fasta \
    --num-threads 30 \
    --force \
    "${outdir}/all_${project_id}_adjacent.fastq.gz"

# ----------------------------
# COMPRESS OUTPUT FASTA FILES
# ----------------------------
echo "[$(date)] Compressing output FASTA files..."

cancer_fasta="${outdir}/all_${project_id}_cancer_miRTrace/qc_passed_reads.all.uncollapsed/all_${project_id}_cancer.fasta"
adjacent_fasta="${outdir}/all_${project_id}_adjacent_miRTrace/qc_passed_reads.all.uncollapsed/all_${project_id}_adjacent.fasta"

if [ ! -f "${cancer_fasta}" ]; then
    echo "ERROR: Expected cancer FASTA not produced by miRTrace: ${cancer_fasta}"
    exit 1
fi
if [ ! -f "${adjacent_fasta}" ]; then
    echo "ERROR: Expected adjacent-normal FASTA not produced by miRTrace: ${adjacent_fasta}"
    exit 1
fi

pigz --processes 30 "${cancer_fasta}"
pigz --processes 30 "${adjacent_fasta}"

echo "[$(date)] Step 2 (QCandFasta) complete."
