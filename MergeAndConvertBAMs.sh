#!/bin/sh

#SBATCH -A naiss2024-22-713
#SBATCH -p shared
#SBATCH -c 30
#SBATCH -t 20:00:00
#SBATCH -J Download
#SBATCH --begin=now
#SBATCH --mail-type=ALL

# This script merges the different samples in each condition, converts the merged bam file to a fastq and compresses it. 

# Assign input arguments to variables
project_id=$1   # First argument: project identifier e.g LUAD
dataset=$2      # Second argument: dataset name e.g TCGA
projPath=$3     # Third argument: base path to project directory

# Define the output directory where results will be stored
outdir=${projPath}/Data/${project_id}/${dataset}

# Load required modules, the PDC module is system specific
ml PDC
ml samtools
ml pigz

# --------------------------
# PROCESSING CANCER SAMPLES
# --------------------------

# Merge all BAM files from the cancer_bams directory into one BAM file
#   --threads 30 : use 30 threads for parallel processing
#   --write-index : automatically generate a BAM index file (.bai)
#   -o : specify output BAM file path
#   $(find ...) : find all BAM files in the cancer_bams folder
samtools merge --threads 30 --write-index -o ${outdir}/all_${project_id}_cancer.bam \
    $(find ${outdir}/bams/cancer_bams/ -name "*.bam")

# Convert the merged cancer BAM file into a FASTQ file
#   --threads 30 : use 30 threads for faster conversion
samtools fastq --threads 30 ${outdir}/all_${project_id}_cancer.bam \
    > ${outdir}/all_${project_id}_cancer.fastq

# Compress the FASTQ file using pigz (parallel gzip)
#   --processes 30 : use 30 parallel threads for compression
pigz --processes 30 ${outdir}/all_${project_id}_cancer.fastq

# --------------------------
# PROCESSING HEALTHY SAMPLES
# --------------------------

# Merge all BAM files from the healthy_bams directory into one BAM file
samtools merge --threads 30 --write-index -o ${outdir}/all_${project_id}_healthy.bam \
    $(find ${outdir}/bams/healthy_bams/ -name "*.bam")

# Convert the merged healthy BAM file into a FASTQ file
samtools fastq --threads 30 ${outdir}/all_${project_id}_healthy.bam \
    > ${outdir}/all_${project_id}_healthy.fastq

# Compress the FASTQ file with pigz using 30 threads
pigz --processes 30 ${outdir}/all_${project_id}_healthy.fastq
