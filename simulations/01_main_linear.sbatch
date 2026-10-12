#!/usr/bin/env bash
#SBATCH --job-name=bcf-gibbs-linear
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=11:50:00
#SBATCH --array=1-600%20
#SBATCH --output=%x-%A_%a.out
set -euo pipefail
# The dispatcher submits offsets 0 and 600 to cover all 1200 replicate batches.
# Internal 600-task chunk; BCF_TASK_OFFSET must be explicit for manual use.
# Use sbatch 00_run_all.sbatch for the full study and evaluation.
bcf_entry_root="${BCF_CODE_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
exec bash "$bcf_entry_root/slurm/worker.sh" main-linear 01_main_linear.sbatch
