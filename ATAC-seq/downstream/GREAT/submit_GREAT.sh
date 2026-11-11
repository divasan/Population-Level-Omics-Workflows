#!/usr/bin/env bash
#SBATCH --job-name=rGREAT_GO_BP
#SBATCH --output=GREAT_%j.out
#SBATCH --error=GREAT_%j.err
#SBATCH --time=72:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=standard

# Example SLURM submission script; customize modules and resource requests.
# Submit from this directory with:
# sbatch --export=ALL,DAR_CSV=/path/to/DAR_results.csv,RESULTS_DIR=/path/to/GREAT_results submit_GREAT.sh

set -euo pipefail

: "${DAR_CSV:?Set DAR_CSV to the PyDESeq2 results CSV}"
: "${RESULTS_DIR:?Set RESULTS_DIR to the desired output directory}"

# Adjust for environment setup
module purge
module load gcc/11.4.0 openmpi/4.1.4 R/4.4.1
module load zlib/1.3.1
# Optional: export R_LIBS_USER=/path/to/personal/R/library

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

printf 'Job ID: %s\nHost: %s\nStarted: %s\n' "${SLURM_JOB_ID:-NA}" "$(hostname)" "$(date)"
Rscript "${SCRIPT_DIR}/run_GREAT_GO_BP.R" "${DAR_CSV}" "${RESULTS_DIR}" "${FDR_CUTOFF:-0.05}" "${LOG2FC_CUTOFF:-0.5}"
printf 'Completed: %s\n' "$(date)"
