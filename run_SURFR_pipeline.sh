#!/bin/bash

# =============================================================================
# SURFR Pipeline — Master Script
#
# Runs the full SURFR pipeline sequentially across one or more cohorts:
#
#   Step 1: MergeAndConvertBAMs.sh  — Merge BAMs, convert to FASTQ
#   Step 2: QCandFasta.sh           — miRTrace QC, produce FASTA
#   Step 3: countKmers.sh           — KMC k-mer counting and merging
#   Step 4: FindCancerSpecificRNAs.R — Identify cancer-specific sequences
#
# SLURM NOTE:
#   Each step submits its own sbatch job and chains to the next via dependencies.
#   The master script itself should be run on the login node or in a screen/tmux 
#   session, NOT submitted as an sbatch job.
#
# Usage:
#   bash run_SURFR_pipeline.sh [--skip-step1] [--dry-run]
#
# =============================================================================

set -euo pipefail

# =============================================================================
# >>>>>>>>>>>>  USER CONFIGURATION — edit everything in this block  <<<<<<<<<<<
# =============================================================================

# --- Project settings --------------------------------------------------------
PROJECT_ID="LUAD"           # Project identifier (e.g. LUAD, BRCA, COAD)
PROJ_PATH="/proj/myproject" # Absolute path to your project root

# Datasets to run through Steps 1–3 (space-separated; Step 4 uses both)
DATASETS="TCGA CPTAC"

# --- Container settings -------------------------------------------------------
# Path to the Singularity sandbox on Dardel
SANDBOX="/cfs/klemming/projects/snic/naiss2024-6-235/programs/surfr_pipeline"

# Every pipeline command runs through this wrapper so it executes inside the
# container. -B /cfs/klemming is required on Dardel
SING_EXEC="singularity exec -B /cfs/klemming ${SANDBOX}"

# --- Tool paths (INSIDE the container — do not change unless the image changes)
MIRTRACE_PATH="/opt/mirtrace/mirtrace"             # miRTrace executable
KMC_BIN_DIR="/opt/kmc/bin"                         # Directory with kmc and kmc_tools
DEKUPL_PATH="/opt/dekupl/bin/mergeTags"            # dekupl-mergeTags executable
RSCRIPT_PATH="Rscript"                             # Rscript inside the container

# --- Step 4 input paths ------------------------------------------------------
MERGED_TABLES_DIR="${PROJ_PATH}/Data/${PROJECT_ID}/merged_tables"
SRA_KMER_TABLE="${PROJ_PATH}/Data/SRA_kmer_counts.txt"
METADATA_DIR="${PROJ_PATH}/Data/${PROJECT_ID}/metadata"
OUTPUT_DIR="${PROJ_PATH}/Results/${PROJECT_ID}"

# --- Script directory --------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# =============================================================================
# END OF USER CONFIGURATION
# =============================================================================

# --- Parse flags -------------------------------------------------------------
SKIP_STEP1=false
DRY_RUN=false

for arg in "$@"; do
  case $arg in
    --skip-step1) SKIP_STEP1=true ;;
    --dry-run)    DRY_RUN=true ;;
    *) echo "WARNING: Unknown argument: $arg" ;;
  esac
done

# --- Helper Function for Slurm Orchestration / Dry Runs ----------------------
# Handles job submissions. If --dry-run is active, it logs the sbatch structure 
# to stderr and passes a dummy job ID to stdout to preserve downstream dependencies.
run_sbatch() {
  if [ "${DRY_RUN}" = true ]; then
    echo -e "\n[DRY-RUN] sbatch command:" >&2
    for arg in "$@"; do
      echo "  $arg" >&2
    done
    # Generate an incremental-like dummy Job ID for tracking tracking dependencies
    echo "mock_job_$((RANDOM % 89999 + 10000))"
  else
    sbatch "$@"
  fi
}

# --- Validate script paths ---------------------------------------------------
for script in MergeAndConvertBAMs.sh QCandFasta.sh countKmers.sh FindCancerSpecificRNAs.R; do
  if [ ! -f "${SCRIPT_DIR}/${script}" ]; then
    echo "ERROR: Pipeline script not found: ${SCRIPT_DIR}/${script}"
    exit 1
  fi
done

echo "============================================================"
echo " SURFR Pipeline Launcher"
echo " Project:  ${PROJECT_ID}"
echo " Datasets: ${DATASETS}"
echo " ProjPath: ${PROJ_PATH}"
echo " Dry Run:  ${DRY_RUN}"
echo " Started:  $(date)"
echo "============================================================"

# =============================================================================
# STEP 1 — Merge BAMs and convert to FASTQ
# =============================================================================
if [ "${SKIP_STEP1}" = false ]; then
  echo ""
  echo "--- Step 1: MergeAndConvertBAMs ---"

  step1_job_ids=()
  for dataset in ${DATASETS}; do
    echo "[$(date)] Submitting Step 1 for dataset: ${dataset}"

    job_id=$(run_sbatch \
      --parsable \
      --account=naiss2026-3-153 \
      --partition=shared \
      --cpus-per-task=30 \
      --time=20:00:00 \
      --job-name=SURFR_Step1_MergeBAMs \
      --wrap="ml PDC singularity && ${SING_EXEC} bash ${SCRIPT_DIR}/MergeAndConvertBAMs.sh ${PROJECT_ID} ${dataset} ${PROJ_PATH}")

    echo "[$(date)] Step 1 job ID for ${dataset}: ${job_id}"
    step1_job_ids+=("${job_id}")
  done

  step1_deps=$(IFS=:; echo "afterok:${step1_job_ids[*]}")
  echo "[$(date)] Step 1 dependencies mapped successfully."
else
  echo ""
  echo "--- Step 1 skipped (--skip-step1 flag set) ---"
  step1_deps=""
fi

# =============================================================================
# STEP 2 — miRTrace QC and FASTA conversion
# =============================================================================
echo ""
echo "--- Step 2: QCandFasta ---"

step2_job_ids=()
for dataset in ${DATASETS}; do
  echo "[$(date)] Submitting Step 2 for dataset: ${dataset}"

  sbatch_args=(
    --parsable
    --account=naiss2026-3-153
    --partition=memory
    --cpus-per-task=30
    --time=20:00:00
    --job-name=SURFR_Step2_QC
  )
  if [ -n "${step1_deps}" ]; then
    sbatch_args+=(--dependency="${step1_deps}")
  fi

  job_id=$(run_sbatch \
    "${sbatch_args[@]}" \
    --wrap="ml PDC singularity && ${SING_EXEC} bash ${SCRIPT_DIR}/QCandFasta.sh ${PROJECT_ID} ${dataset} ${PROJ_PATH} ${MIRTRACE_PATH}")

  echo "[$(date)] Step 2 job ID for ${dataset}: ${job_id}"
  step2_job_ids+=("${job_id}")
done

step2_deps=$(IFS=:; echo "afterok:${step2_job_ids[*]}")
echo "[$(date)] Step 2 dependencies mapped successfully."

# =============================================================================
# STEP 3 — k-mer counting and merging
# =============================================================================
echo ""
echo "--- Step 3: countKmers ---"

step3_job_ids=()
for dataset in ${DATASETS}; do
  echo "[$(date)] Submitting Step 3 for dataset: ${dataset}"

  step3_sbatch_args=(
    --parsable
    --account=naiss2026-3-153
    --partition=memory
    --cpus-per-task=30
    --time=5:00:00
    --job-name=SURFR_Step3_Kmers
  )
  if [ -n "${step2_deps}" ]; then
    step3_sbatch_args+=(--dependency="${step2_deps}")
  fi

  job_id=$(run_sbatch \
    "${step3_sbatch_args[@]}" \
    --wrap="ml PDC singularity && ${SING_EXEC} bash ${SCRIPT_DIR}/countKmers.sh ${PROJECT_ID} ${dataset} ${PROJ_PATH} ${KMC_BIN_DIR}")

  echo "[$(date)] Step 3 job ID for ${dataset}: ${job_id}"
  step3_job_ids+=("${job_id}")
done

step3_deps=$(IFS=:; echo "afterok:${step3_job_ids[*]}")
echo "[$(date)] Step 3 dependencies mapped successfully."

# =============================================================================
# STEP 4 — Identify cancer-specific RNAs (R script)
# =============================================================================
echo ""
echo "--- Step 4: FindCancerSpecificRNAs ---"
echo "[$(date)] Generating delayed workspace linking logic for compute node execution..."

# Dynamically build a execution string string for Step 4.
# This builds environment directories and runs symlinks ON the compute node.
wrap_cmd="mkdir -p '${OUTPUT_DIR}' && mkdir -p '${MERGED_TABLES_DIR}'"

for dataset in ${DATASETS}; do
  src="${PROJ_PATH}/Data/${PROJECT_ID}/${dataset}/all_${PROJECT_ID}_${dataset}_17mers_merged.txt"
  dst="${MERGED_TABLES_DIR}/all_${PROJECT_ID}_${dataset}_17mers_merged.txt"
  if [ "${src}" != "${dst}" ]; then
    wrap_cmd="${wrap_cmd} && ln -sf '${src}' '${dst}'"
  fi
done

# Append the actual R execution script pipeline invocation 
wrap_cmd="${wrap_cmd} && ml PDC singularity && ${SING_EXEC} ${RSCRIPT_PATH} ${SCRIPT_DIR}/FindCancerSpecificRNAs.R \
  ${PROJECT_ID} \
  ${MERGED_TABLES_DIR} \
  ${SRA_KMER_TABLE} \
  ${METADATA_DIR} \
  ${OUTPUT_DIR} \
  ${DEKUPL_PATH}"

step4_sbatch_args=(
  --parsable
  --account=naiss2026-3-153
  --partition=memory
  --cpus-per-task=8
  --time=4:00:00
  --job-name=SURFR_Step4_RScript
)
if [ -n "${step3_deps}" ]; then
  step4_sbatch_args+=(--dependency="${step3_deps}")
fi

echo "[$(date)] Submitting Step 4..."
step4_job_id=$(run_sbatch "${step4_sbatch_args[@]}" --wrap="${wrap_cmd}")
echo "[$(date)] Step 4 job ID: ${step4_job_id}"

# =============================================================================
# SUMMARY
# =============================================================================
echo ""
echo "============================================================"
echo " All tasks processed."
echo ""
echo " Scheduled Job Chain Layout:"
if [ "${SKIP_STEP1}" = false ]; then
  echo "   Step 1 (MergeAndConvertBAMs): ${step1_job_ids[*]}"
fi
echo "   Step 2 (QCandFasta):          ${step2_job_ids[*]}"
echo "   Step 3 (countKmers):          ${step3_job_ids[*]}"
echo "   Step 4 (FindCancerSpecificRNAs): ${step4_job_id}"
echo "============================================================"