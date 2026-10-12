#!/usr/bin/env bash
# Analysis only. Default: normal linear/nonlinear and orthogonal simulations.
# Simulation pipeline: afterok on all six arrays; automatic continuation is OFF.
# Manual recovery/IHDP: wait for ALL required fits, including any continuations.
# Use "ihdp" or "all" explicitly if wanted.
#SBATCH --job-name=bcf-gibbs-analysis
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=11:50:00
#SBATCH --output=%x-%j.out
set -euo pipefail
bcf_entry_root="${BCF_CODE_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
exec bash "$bcf_entry_root/slurm/rscript.sh" run.R analyse "${1:-simulations}"
