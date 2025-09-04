#!/bin/sh

#SBATCH -A naiss2024-22-713
#SBATCH -p memory
#SBATCH -c 30
#SBATCH -t 30:00
#SBATCH -J make_kmc_kmers
#SBATCH --begin=now
#SBATCH --mail-type=ALL

# ----------------------------
# INPUT ARGUMENTS
# ----------------------------
project_id=$1   # Project identifier (first script argument)
dataset=$2      # Dataset name (second script argument)
projPath=$3     # Project base path (third script argument)
k_length=17     # K-mer length (here fixed at 17)

# ----------------------------
# DEFINE PROGRAM PATHS
# ----------------------------
# Paths to KMC and KMC tools executables
kmc=/path/to/KMC/bin/kmc
kmc_tools=/path/to/KMC/KMC/bin/kmc_tools

# Output directory for results
outdir=${projPath}/Data/${project_id}/${dataset}

# ----------------------------
# K-MER COUNTING (CANCER)
# ----------------------------
# Count k-mers in the cancer FASTA file
#   -fa            : input is FASTA format
#   -t30           : use 30 threads
#   -b             : use RAM-only mode (no temporary disk storage)
#   -ci30          : minimum k-mer occurrence cutoff (≥30 reads required)
#   -cs4294967296  : maximum k-mer count value (very large upper bound)
#   -k${k_length}  : length of k-mers to count
# Output: binary KMC database named all_${project_id}_cancer_${k_length}mers_FILTER30
${kmc} -fa -t30 -b -ci30 -cs4294967296 -k${k_length} \
${outdir}/all_${project_id}_cancer_miRTrace/qc_passed_reads.all.uncollapsed/all_${project_id}_cancer.fasta.gz \
${outdir}/all_${project_id}_cancer_${k_length}mers_FILTER30 .

# Dump k-mer counts to a human-readable text file (sorted by KMC’s internal order)
#   -s : output in a simple tab-delimited format (kmer \t count)
${kmc_tools} -t30 transform ${outdir}/all_${project_id}_cancer_${k_length}mers_FILTER30 dump -s \
${outdir}/all_${project_id}_cancer_${k_length}mers_FILTER30.txt

# ----------------------------
# K-MER COUNTING (HEALTHY)
# ----------------------------
# Same as above, but with a looser cutoff (-ci1 means include all k-mers with ≥1 occurrence)
${kmc} -fa -t30 -b -ci1 -cs4294967296 -k${k_length} \
${outdir}/all_${project_id}_healthy_miRTrace/qc_passed_reads.all.uncollapsed/all_${project_id}_healthy.fasta.gz \
${outdir}/all_${project_id}_healthy_${k_length}mers .

# Dump healthy k-mer counts to text file
${kmc_tools} -t30 transform ${outdir}/all_${project_id}_healthy_${k_length}mers dump -s \
${outdir}/all_${project_id}_healthy_${k_length}mers.txt

# ----------------------------
# SORT K-MER FILES
# ----------------------------
# Sort both cancer and healthy k-mer count files by k-mer sequence (field 1)
#   -t$'\t' : tab delimiter
#   -k1,1   : sort by first column only
sort -t$'\t' -k1,1 ${outdir}/all_${project_id}_cancer_${k_length}mers_FILTER30.txt \
    > ${outdir}/all_${project_id}_cancer_${k_length}mers_FILTER30_sorted.txt

sort -t$'\t' -k1,1 ${outdir}/all_${project_id}_healthy_${k_length}mers.txt \
    > ${outdir}/all_${project_id}_healthy_${k_length}mers_sorted.txt

# ----------------------------
# MERGE K-MER COUNTS
# ----------------------------
# Use 'join' to merge cancer and healthy k-mer counts by k-mer sequence
#   -t$'\t' : tab delimiter
#   -a1     : include all cancer k-mers, even if not present in healthy
#   -e 0    : fill missing values with "0"
#   -o 0,1.2,2.2 : output columns: (kmer, cancer_count, healthy_count)
#   -1 1 -2 1    : join on first column of both files
join -t$'\t' -a1 -e 0 -o 0,1.2,2.2 -1 1 -2 1 \
${outdir}/all_${project_id}_cancer_${k_length}mers_FILTER30_sorted.txt \
${outdir}/all_${project_id}_healthy_${k_length}mers_sorted.txt \
> ${outdir}/all_${project_id}_${dataset}_${k_length}mers_merged.txt
