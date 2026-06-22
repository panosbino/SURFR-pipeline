#!/bin/bash
#SBATCH -A naiss2024-22-713
#SBATCH -p shared
#SBATCH -c 30
#SBATCH -t 20:00:00
#SBATCH -J MergeAndConvertBAMs
#SBATCH --begin=now
#SBATCH --mail-type=ALL

# =============================================================================
# SURFR Pipeline - Step 1: MergeAndConvertBAMs.sh
#
# Merges all BAM files per condition (cancer / adjacent-normal), converts each
# merged BAM to FASTQ, and compresses the output with pigz.
#
# NOTE: If your input files are already in FASTQ format, skip this script and
#       pass the FASTQ files directly to QCandFasta.sh (Step 2).
#
# Usage (direct):
#   bash MergeAndConvertBAMs.sh <project_id> <dataset> <projPath>
#
# Usage (SLURM):
#   sbatch MergeAndConvertBAMs.sh <project_id> <dataset> <projPath>
#
# Arguments:
#   project_id  - Project identifier, e.g. LUAD
#   dataset     - Dataset name, e.g. TCGA
#   projPath    - Absolute path to the project root directory
#
# Expected input directory structure:
#   <projPath>/Data/<project_id>/<dataset>/bams/cancer_bams/   (*.bam files)
#   <projPath>/Data/<project_id>/<dataset>/bams/adjacent_bams/ (*.bam files)
#
# Outputs:
#   <outdir>/all_<project_id>_cancer.fastq.gz
#   <outdir>/all_<project_id>_adjacent.fastq.gz
# =============================================================================

set -euo pipefail   # Exit on error, unset variable, or pipe failure

# ----------------------------
# INPUT ARGUMENTS
# ----------------------------
if [ "$#" -ne 3 ]; then
    echo "ERROR: Expected 3 arguments, got $#"
    echo "Usage: $0 <project_id> <dataset> <projPath>"
    exit 1
fi

project_id=$1   # e.g. LUAD
dataset=$2      # e.g. TCGA
projPath=$3     # e.g. /proj/myproject

# ----------------------------
# PATHS
# ----------------------------
outdir=${projPath}/Data/${project_id}/${dataset}

# Verify the output directory exists
if [ ! -d "${outdir}" ]; then
    echo "ERROR: Output directory does not exist: ${outdir}"
    exit 1
fi

# ----------------------------
# LOAD MODULES
# (PDC module is system-specific; adjust for your cluster environment)
# ----------------------------
ml PDC
ml samtools
ml pigz

# ----------------------------
# PROCESSING CANCER SAMPLES
# ----------------------------
echo "[$(date)] Merging cancer BAMs..."

cancer_bam_dir="${outdir}/bams/cancer_bams"
if [ ! -d "${cancer_bam_dir}" ]; then
    echo "ERROR: Cancer BAM directory not found: ${cancer_bam_dir}"
    exit 1
fi

# Merge all cancer BAM files into one; --write-index creates .bai automatically
samtools merge \
    --threads 30 \
    --write-index \
    -o "${outdir}/all_${project_id}_cancer.bam" \
    $(find "${cancer_bam_dir}" -name "*.bam" | sort)

echo "[$(date)] Converting cancer BAM to FASTQ..."
samtools fastq \
    --threads 30 \
    "${outdir}/all_${project_id}_cancer.bam" \
    > "${outdir}/all_${project_id}_cancer.fastq"

echo "[$(date)] Compressing cancer FASTQ..."
pigz --processes 30 "${outdir}/all_${project_id}_cancer.fastq"

# ----------------------------
# PROCESSING ADJACENT-NORMAL SAMPLES
# ----------------------------
echo "[$(date)] Merging adjacent-normal BAMs..."

adjacent_bam_dir="${outdir}/bams/adjacent_bams"
if [ ! -d "${adjacent_bam_dir}" ]; then
    echo "ERROR: Adjacent-normal BAM directory not found: ${adjacent_bam_dir}"
    exit 1
fi

samtools merge \
    --threads 30 \
    --write-index \
    -o "${outdir}/all_${project_id}_adjacent.bam" \
    $(find "${adjacent_bam_dir}" -name "*.bam" | sort)

echo "[$(date)] Converting adjacent-normal BAM to FASTQ..."
samtools fastq \
    --threads 30 \
    "${outdir}/all_${project_id}_adjacent.bam" \
    > "${outdir}/all_${project_id}_adjacent.fastq"

echo "[$(date)] Compressing adjacent-normal FASTQ..."
pigz --processes 30 "${outdir}/all_${project_id}_adjacent.fastq"

echo "[$(date)] Step 1 (MergeAndConvertBAMs) complete."
