#!/bin/bash
#SBATCH -A naiss2024-22-713
#SBATCH -p memory
#SBATCH -c 30
#SBATCH -t 5:00:00
#SBATCH -J countKmers
#SBATCH --begin=now
#SBATCH --mail-type=ALL

# =============================================================================
# SURFR Pipeline - Step 3: countKmers.sh
#
# Uses KMC to count k-mers in cancer and adjacent-normal FASTA files.
# Filters cancer k-mers to those with ≥30 counts, dumps both tables to text,
# sorts them by k-mer sequence, then left-joins on cancer to produce a merged
# count table for use in Step 4.
#
# Usage (direct):
#   bash countKmers.sh <project_id> <dataset> <projPath> <kmc_bin_dir>
#
# Usage (SLURM):
#   sbatch countKmers.sh <project_id> <dataset> <projPath> <kmc_bin_dir>
#
# Arguments:
#   project_id  - Project identifier, e.g. LUAD
#   dataset     - Dataset name, e.g. TCGA
#   projPath    - Absolute path to the project root directory
#   kmc_bin_dir - Directory containing both 'kmc' and 'kmc_tools' executables
#
# Inputs (from Step 2):
#   <outdir>/all_<project_id>_cancer_miRTrace/qc_passed_reads.all.uncollapsed/all_<project_id>_cancer.fasta.gz
#   <outdir>/all_<project_id>_adjacent_miRTrace/qc_passed_reads.all.uncollapsed/all_<project_id>_adjacent.fasta.gz
#
# Outputs:
#   <outdir>/all_<project_id>_<dataset>_17mers_merged.txt
#     Columns: kmer <TAB> cancer_count <TAB> adjacent_count
# =============================================================================

set -euo pipefail

# ----------------------------
# INPUT ARGUMENTS
# ----------------------------
if [ "$#" -ne 4 ]; then
    echo "ERROR: Expected 4 arguments, got $#"
    echo "Usage: $0 <project_id> <dataset> <projPath> <kmc_bin_dir>"
    exit 1
fi

project_id=$1       # e.g. LUAD
dataset=$2          # e.g. TCGA
projPath=$3         # e.g. /proj/myproject
kmc_bin_dir=$4      # e.g. /path/to/KMC/bin

k_length=17         # K-mer length (fixed at 17 for the SURFR pipeline)

# ----------------------------
# DEFINE PROGRAM PATHS
# FIX: both kmc and kmc_tools live in the same bin directory
# ----------------------------
kmc="${kmc_bin_dir}/kmc"
kmc_tools="${kmc_bin_dir}/kmc_tools"

# Validate executables
for exe in "${kmc}" "${kmc_tools}"; do
    if [ ! -x "${exe}" ]; then
        echo "ERROR: Executable not found or not executable: ${exe}"
        exit 1
    fi
done

# ----------------------------
# PATHS
# ----------------------------
outdir=${projPath}/Data/${project_id}/${dataset}

cancer_fasta="${outdir}/all_${project_id}_cancer_miRTrace/qc_passed_reads.all.uncollapsed/all_${project_id}_cancer.fasta.gz"
adjacent_fasta="${outdir}/all_${project_id}_adjacent_miRTrace/qc_passed_reads.all.uncollapsed/all_${project_id}_adjacent.fasta.gz"

# Validate inputs from Step 2
for f in "${cancer_fasta}" "${adjacent_fasta}"; do
    if [ ! -f "${f}" ]; then
        echo "ERROR: Input FASTA not found (did Step 2 complete?): ${f}"
        exit 1
    fi
done

# ----------------------------
# K-MER COUNTING — CANCER
# ----------------------------
echo "[$(date)] Counting cancer k-mers..."

# -fa        : input is FASTA format
# -t30       : use 30 threads
# -b         : RAM-only mode (no temporary disk storage)
# -ci30      : minimum k-mer occurrence cutoff (≥30 reads required)
# -cs4294967296 : maximum k-mer count (large upper bound)
# -k17       : k-mer length
# Final '.'  : working directory for KMC temporary files
"${kmc}" -fa -t30 -b -ci30 -cs4294967296 -k${k_length} \
    "${cancer_fasta}" \
    "${outdir}/all_${project_id}_cancer_${k_length}mers_FILTER30" \
    .

# Dump cancer k-mer counts to a tab-delimited text file
# -s : simple output format (kmer <TAB> count)
"${kmc_tools}" -t30 transform \
    "${outdir}/all_${project_id}_cancer_${k_length}mers_FILTER30" \
    dump -s \
    "${outdir}/all_${project_id}_cancer_${k_length}mers_FILTER30.txt"

# ----------------------------
# K-MER COUNTING — ADJACENT-NORMAL
# ----------------------------
echo "[$(date)] Counting adjacent-normal k-mers..."

# -ci1 : include all k-mers with ≥1 occurrence (no minimum filter for adjacent)
"${kmc}" -fa -t30 -b -ci1 -cs4294967296 -k${k_length} \
    "${adjacent_fasta}" \
    "${outdir}/all_${project_id}_adjacent_${k_length}mers" \
    .

"${kmc_tools}" -t30 transform \
    "${outdir}/all_${project_id}_adjacent_${k_length}mers" \
    dump -s \
    "${outdir}/all_${project_id}_adjacent_${k_length}mers.txt"

# ----------------------------
# SORT K-MER FILES
# ----------------------------
echo "[$(date)] Sorting k-mer count files..."

sort -t$'\t' -k1,1 \
    "${outdir}/all_${project_id}_cancer_${k_length}mers_FILTER30.txt" \
    > "${outdir}/all_${project_id}_cancer_${k_length}mers_FILTER30_sorted.txt"

sort -t$'\t' -k1,1 \
    "${outdir}/all_${project_id}_adjacent_${k_length}mers.txt" \
    > "${outdir}/all_${project_id}_adjacent_${k_length}mers_sorted.txt"

# ----------------------------
# MERGE K-MER COUNTS (left join on cancer)
# ----------------------------
echo "[$(date)] Merging cancer and adjacent-normal k-mer tables..."

# -t$'\t'    : tab delimiter
# -a1        : keep all cancer k-mers even if absent in adjacent-normal
# -e 0       : fill missing values with "0"
# -o 0,1.2,2.2 : output: kmer, cancer_count, adjacent_count
# -1 1 -2 1  : join on column 1 of both files
join -t$'\t' -a1 -e 0 -o 0,1.2,2.2 -1 1 -2 1 \
    "${outdir}/all_${project_id}_cancer_${k_length}mers_FILTER30_sorted.txt" \
    "${outdir}/all_${project_id}_adjacent_${k_length}mers_sorted.txt" \
    > "${outdir}/all_${project_id}_${dataset}_${k_length}mers_merged.txt"

echo "[$(date)] Step 3 (countKmers) complete."
echo "[$(date)] Output: ${outdir}/all_${project_id}_${dataset}_${k_length}mers_merged.txt"
