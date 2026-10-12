#!/usr/bin/env bash
# One submission: normal + orthogonal simulations, IHDP and automatic exports.
# Submit from this extracted folder. No environment/path setup is required.
#SBATCH --job-name=bcf-run-all
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=%x-%j.out
set -euo pipefail
# Slurm runs a spool copy: use the actual submission directory, not a stale
# BCF_CODE_ROOT inherited from a previous workflow.
bcf_easy_root="${SLURM_SUBMIT_DIR:-$PWD}"
[[ -f "$bcf_easy_root/run.R" && -f "$bcf_easy_root/slurm/run_all.sh" ]] || {
  echo 'Run sbatch 00_run_all.sbatch from the simulations folder.' >&2
  exit 2
}
exec bash "$bcf_easy_root/slurm/run_all.sh" "$bcf_easy_root"
