#!/usr/bin/env bash
#SBATCH --job-name=bcf-gibbs-ihdp
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=11:50:00
#SBATCH --array=1-8%8
#SBATCH --output=%x-%A_%a.out
set -euo pipefail
bcf_entry_root="${BCF_CODE_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
exec bash "$bcf_entry_root/slurm/worker.sh" ihdp 04_ihdp.sbatch
