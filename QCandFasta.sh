#!/bin/sh

#SBATCH -A naiss2024-22-713
#SBATCH -p memory
#SBATCH -c 30
#SBATCH -t 20:00:00
#SBATCH -J mirtrace
#SBATCH --begin=now
#SBATCH --mail-type=ALL

# This script runs miRTrace quality control on cancer and healthy datasets,
# then compresses the resulting FASTA files.

# ----------------------------
# INPUT ARGUMENTS
# ----------------------------
project_id=$1   # Project identifier (passed as first argument) e.g LUAD
dataset=$2      # Dataset name (passed as second argument) e.g TCGA
projPath=$3     # Path to the project directory (passed as third argument)

# Define the base output directory for all results
outdir=${projPath}/Data/${project_id}/${dataset}

# ----------------------------
# RUN miRTrace QC (CANCER)
# ----------------------------
# Run miRTrace quality control on the cancer FASTQ file
#   --species hsa : assumes human (hsa) reads
#   --output-dir : directory to store QC results
#   --write-fasta : outputs passed reads in FASTA format
#   --uncollapse-fasta : keeps duplicate reads uncollapsed
#   --num-threads 30 : use 30 threads for faster execution
#   --force : overwrite existing results if present
/path/to/mirtrace qc \
    --species hsa \
    ${outdir}/all_${project_id}_cancer.fastq.gz \
    --output-dir ${outdir}/all_${project_id}_cancer_miRTrace \
    --write-fasta \
    --uncollapse-fasta \
    --num-threads 30 \
    --force

# ----------------------------
# RUN miRTrace QC (HEALTHY)
# ----------------------------
/path/to/mirtrace qc \
    --species hsa \
    ${outdir}/all_${project_id}_healthy.fastq.gz \
    --output-dir ${outdir}/all_${project_id}_healthy_miRTrace \
    --write-fasta \
    --uncollapse-fasta \
    --num-threads 30 \
    --force

# ----------------------------
# COMPRESS OUTPUT FASTA FILES
# ----------------------------
# Compress the uncollapsed FASTA files from miRTrace QC using pigz
#   --processes 30 : parallel compression with 30 threads
pigz --processes 30 ${outdir}/all_${project_id}_cancer_miRTrace/qc_passed_reads.all.uncollapsed/all_${project_id}_cancer.fasta
pigz --processes 30 ${outdir}/all_${project_id}_healthy_miRTrace/qc_passed_reads.all.uncollapsed/all_${project_id}_healthy.fasta
