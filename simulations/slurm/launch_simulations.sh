#!/usr/bin/env bash
# Internal fixed-layout submitter. The public 00 script runs all tiny checks
# first. The new all-studies entry explicitly requests the additional IHDP stage.
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/paths.sh"
bcf_launch_scope="${1:-simulations}"
[[ "$bcf_launch_scope" == simulations || "$bcf_launch_scope" == all ]] || { echo 'Invalid launch scope.' >&2; exit 2; }
if [[ "$bcf_launch_scope" == all ]]; then
  [[ "${BCF_EASY_RUN_ID:-}" =~ ^[1-9][0-9]*$ && "${BCF_AUTO_IHDP_ANALYSIS:-}" == 1 ]] || {
    echo 'Use 00_run_all.sbatch for the all-studies workflow.' >&2; exit 2;
  }
fi
command -v sbatch >/dev/null 2>&1 || { echo 'sbatch is unavailable; nothing submitted.' >&2; exit 2; }
# Keep all submissions in one transaction. Advisory flock can be node-local on
# BeeGFS, so also claim a same-filesystem hard link. Never steal another owner.
mkdir -p -- "$BCF_OUTPUT/runtime-benchmark"
exec 8>"$BCF_OUTPUT/runtime-benchmark/launch.lock"
flock -n 8 || { echo 'Another simulation launch transaction is active; nothing submitted.' >&2; exit 2; }
bcf_launch_claim="$BCF_OUTPUT/runtime-benchmark/launch.owner"
bcf_launch_token=$(mktemp "$BCF_OUTPUT/runtime-benchmark/.launch_owner.XXXXXX")
printf '%s\t%s\t%s\n' "${SLURM_JOB_ID:-manual}" "${HOSTNAME:-unknown}" "${bcf_launch_token##*/}" > "$bcf_launch_token"
if ! ln "$bcf_launch_token" "$bcf_launch_claim" 2>/dev/null; then
  rm -f -- "$bcf_launch_token"
  echo "Launch owner already exists: $bcf_launch_claim" >&2
  echo 'Verify its recorded job has ended before removing that exact owner file after node failure/SIGKILL.' >&2
  exit 2
fi
bcf_launch_cleanup() {
  if [[ -f "$bcf_launch_claim" ]] && cmp -s "$bcf_launch_token" "$bcf_launch_claim"; then rm -f -- "$bcf_launch_claim"; fi
  rm -f -- "$bcf_launch_token"
}
trap bcf_launch_cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
echo 'Fixed layout: six 600-task arrays, each limited to 20 simultaneous tasks, 4 CPUs and 16G per task.'
echo 'Maximum running fit allocation: 120 tasks / 480 CPUs / 1920G, subject to shared account and scheduler limits.'
echo 'Inspect per-quantity and grouped diagnostics after fitting; see ../README.md for the rerun protocol.'
echo 'The 11:50 job limit is per allocation, not a guarantee that the entire study finishes in 12 hours.'
bcf_production_ids=()
for bcf_entry in 01_main_linear.sbatch 02_main_nonlinear.sbatch 03_orthogonalisation.sbatch; do
  for bcf_offset in 0 600; do
    if ! bcf_submitted=$(sbatch --parsable --chdir="$BCF_CODE_ROOT" \
      --array=1-600%20 \
      --export="ALL,BCF_CONTINUATION=0,BCF_TASK_OFFSET=$bcf_offset,BCF_AUTOCONTINUE=0" \
      "$BCF_CODE_ROOT/$bcf_entry"); then
      echo "Submission failed for $bcf_entry offset=$bcf_offset. Previously printed jobs remain submitted; inspect them before retrying." >&2
      exit 2
    fi
    if [[ ! "$bcf_submitted" =~ ^[1-9][0-9]*(\;[^[:space:]]+)?$ ]]; then
      echo "Unexpected sbatch reply for $bcf_entry offset=$bcf_offset: $bcf_submitted. Inspect squeue before retrying." >&2
      exit 2
    fi
    printf 'Submitted %-30s offset=%-3s array=1-600%%20 job=%s\n' "$bcf_entry" "$bcf_offset" "$bcf_submitted"
    bcf_production_ids+=("${bcf_submitted%%;*}")
    if [[ "$bcf_launch_scope" == all ]]; then
      printf '%s\t%s\t%s\n' "$bcf_entry" "$bcf_offset" "$bcf_submitted" >> "$BCF_OUTPUT/SUBMITTED_JOBS.tsv"
    fi
  done
done
# Auto-continuation is OFF and budget stops exit75. All six afterok dependencies
# therefore require genuine completion; a failed batch blocks auto-evaluation.
bcf_analysis_dependency=$(IFS=:; printf '%s' "${bcf_production_ids[*]}")
if ! bcf_analysis=$(sbatch --parsable --chdir="$BCF_CODE_ROOT" \
  --dependency="afterok:$bcf_analysis_dependency" --export=ALL,BCF_AUTOCONTINUE=0 \
  "$BCF_CODE_ROOT/05_analyse.sbatch" simulations); then
  echo 'Analysis submission failed; the six printed fitting arrays remain submitted. After every fit completes, submit 05_analyse.sbatch manually.' >&2
  exit 2
fi
if [[ ! "$bcf_analysis" =~ ^[1-9][0-9]*(\;[^[:space:]]+)?$ ]]; then
  echo "Unexpected analysis sbatch reply: $bcf_analysis. Inspect squeue before retrying." >&2
  exit 2
fi
printf 'Submitted automatic evaluation after all six arrays: %s\n' "$bcf_analysis"
if [[ "$bcf_launch_scope" == all ]]; then
  printf 'simulation-analysis\t-\t%s\n' "$bcf_analysis" >> "$BCF_OUTPUT/SUBMITTED_JOBS.tsv"
  # Depend on successful FITS, not analysis: convergence-draft export exits2
  # and must not prevent the independent IHDP comparison from fitting.
  if ! bcf_ihdp=$(sbatch --parsable --chdir="$BCF_CODE_ROOT" \
    --dependency="afterok:$bcf_analysis_dependency" \
    --export=ALL,BCF_CONTINUATION=0,BCF_TASK_OFFSET=0,BCF_AUTOCONTINUE=1,BCF_AUTO_IHDP_ANALYSIS=1 \
    "$BCF_CODE_ROOT/04_ihdp.sbatch"); then
    echo 'IHDP submission failed; the printed simulation/analysis jobs remain submitted. Inspect SUBMITTED_JOBS.tsv before recovery.' >&2
    exit 2
  fi
  [[ "$bcf_ihdp" =~ ^[1-9][0-9]*(\;[^[:space:]]+)?$ ]] || {
    echo "Unexpected IHDP sbatch reply: $bcf_ihdp. Inspect squeue before retrying." >&2; exit 2;
  }
  printf 'ihdp-initial\t-\t%s\n' "$bcf_ihdp" >> "$BCF_OUTPUT/SUBMITTED_JOBS.tsv"
  printf 'Submitted IHDP after all six simulation FIT arrays: %s\n' "$bcf_ihdp"
  echo 'IHDP analysis will be submitted once all eight shards, including continuations, have truly finished.'
fi
echo "Simulation output: $BCF_OUTPUT"
printf 'Status: BCF_OUTPUT=%q bash slurm/rscript.sh run.R status %s\n' "$BCF_OUTPUT" "$bcf_launch_scope"
echo 'If a fit fails or reaches its budget, automatic evaluation remains blocked. Recover missing fits, verify status, then submit 05_analyse.sbatch manually.'
