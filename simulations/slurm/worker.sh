#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/paths.sh"
bcf_suite="${1:?suite required}"
bcf_entry="${2:?entry script required}"
bcf_local_task="${SLURM_ARRAY_TASK_ID:?This is an sbatch array worker}"
case "$bcf_suite" in
  main-linear|main-nonlinear|ortho)
    [[ -n "${BCF_TASK_OFFSET:-}" ]] || {
  echo 'Simulation entries are 600-task chunks, not a full study. Use sbatch 00_run_all.sbatch for a fresh all-studies run, or set BCF_TASK_OFFSET=0 or 600 for an explicit partial/manual resume.' >&2
      exit 2
    }
    bcf_offset="$BCF_TASK_OFFSET"
    [[ "$bcf_offset" == 0 || "$bcf_offset" == 600 ]] || { echo "BCF_TASK_OFFSET must be 0 or 600" >&2; exit 2; }
    [[ "$bcf_local_task" =~ ^[1-9][0-9]{0,2}$ ]] && (( bcf_local_task <= 600 )) || {
      echo "Simulation array index must be 1..600" >&2; exit 2;
    }
    bcf_task=$((bcf_offset + bcf_local_task))
    bcf_tasks=1200
    bcf_lock_key="batch_$bcf_task"
    bcf_r_args=(run.R batch "$bcf_suite" "$bcf_task")
    bcf_autocontinue="${BCF_AUTOCONTINUE:-0}"
    ;;
  ihdp)
    [[ "$bcf_local_task" =~ ^[1-8]$ ]] || { echo "IHDP array index must be 1..8" >&2; exit 2; }
    bcf_offset=0
    bcf_task="$bcf_local_task"
    bcf_tasks=8
    bcf_lock_key="shard_$bcf_task"
    bcf_r_args=(run.R worker "$bcf_suite" "$bcf_task" "$bcf_tasks")
    bcf_autocontinue="${BCF_AUTOCONTINUE:-1}"
    ;;
  runtime-benchmark)
    [[ "$bcf_local_task" =~ ^([1-9]|1[0-8])$ ]] || { echo "Runtime benchmark array index must be 1..18" >&2; exit 2; }
    bcf_offset=0
    bcf_task="$bcf_local_task"
    bcf_tasks=18
    bcf_lock_key="task_$bcf_task"
    bcf_r_args=(run.R benchmark-task "$bcf_task")
    # The optional report needs complete full-length timings. A budget stop
    # remains a nonzero exit, never a hidden continuation or shortened timing.
    bcf_autocontinue=0
    ;;
  *) echo "Unknown study suite: $bcf_suite" >&2; exit 2 ;;
esac
[[ "${BCF_CONTINUATION:-0}" =~ ^(0|[1-9][0-9]{0,2})$ ]] || {
  echo "BCF_CONTINUATION must be a nonnegative integer below 1000" >&2; exit 2;
}
[[ -f "$BCF_CODE_ROOT/$bcf_entry" ]] || exit 2
mkdir -p -- "$BCF_OUTPUT/$bcf_suite/locks"
# The OS releases the local advisory lock on timeout/cancellation. The separate
# cross-node claim below can need manual recovery after node failure/SIGKILL.
exec 9>"$BCF_OUTPUT/$bcf_suite/locks/${bcf_lock_key}.lock"
flock -n 9 || {
  echo "This study task already has an active worker; duplicate did not complete its fits." >&2
  # A skipped duplicate must not satisfy afterok and trigger early evaluation.
  exit 2
}
# BeeGFS flock can be node-local. A server-atomic, no-overwrite hard link adds
# cross-node ownership; no automatic stealing of another allocation's claim.
bcf_claim="$BCF_OUTPUT/$bcf_suite/locks/${bcf_lock_key}.owner"
bcf_token=$(mktemp "$BCF_OUTPUT/$bcf_suite/locks/.owner_${bcf_lock_key}.XXXXXX")
printf '%s\t%s\t%s\n' "${SLURM_JOB_ID}" "${HOSTNAME:-unknown}" "${bcf_token##*/}" > "$bcf_token"
if ! ln "$bcf_token" "$bcf_claim" 2>/dev/null; then
  rm -f -- "$bcf_token"
  echo "Cross-node owner exists for this task: $bcf_claim"
  echo "Do not resubmit an active task. After a node failure/SIGKILL, verify its recorded Slurm job has ended before removing that exact owner file."
  exit 2
fi
bcf_cleanup() {
  if [[ -f "$bcf_claim" ]] && cmp -s "$bcf_token" "$bcf_claim"; then rm -f -- "$bcf_claim"; fi
  rm -f -- "$bcf_token"
}
trap bcf_cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
cd "$BCF_CODE_ROOT"
echo "Study=$bcf_suite task=$bcf_task/$bcf_tasks local_index=$bcf_local_task offset=$bcf_offset output=$BCF_OUTPUT"
set +e
# Leave 35 minutes below the Slurm limit for termination/continuation handling.
# A single interrupted chain restarts from its seed; completed chains survive.
timeout --signal=TERM --kill-after=30s 40500 \
  bash "$BCF_CODE_ROOT/slurm/rscript.sh" "${bcf_r_args[@]}"
bcf_rc=$?
set -e
if [[ "$bcf_rc" == 75 ]]; then
  if [[ "$bcf_autocontinue" != 1 ]]; then
    echo "Time budget reached. Resubmit $bcf_entry with --array=$bcf_local_task and BCF_TASK_OFFSET=$bcf_offset to continue from checkpoints."
    exit 75
  fi
  bcf_next=$(( ${BCF_CONTINUATION:-0} + 1 ))
  if (( bcf_next > 200 )); then echo "Continuation safety limit reached; inspect status before resubmitting." >&2; exit 2; fi
  bcf_dependency="${SLURM_ARRAY_JOB_ID:-$SLURM_JOB_ID}_${bcf_local_task}"
  bcf_new=$(sbatch --parsable --chdir="$BCF_CODE_ROOT" --array="$bcf_local_task" \
    --dependency="afterany:$bcf_dependency" \
    --export="ALL,BCF_CONTINUATION=$bcf_next,BCF_TASK_OFFSET=$bcf_offset" "$BCF_CODE_ROOT/$bcf_entry")
  [[ "$bcf_new" =~ ^[1-9][0-9]*(\;[^[:space:]]+)?$ ]] || {
    echo "Unexpected continuation sbatch reply: $bcf_new. Check squeue before retrying; completed checkpoints remain." >&2
    exit 2
  }
  printf '%s\t%s\t%s\t%s\n' "$bcf_suite" "$bcf_task" "$bcf_next" "$bcf_new" \
    >> "$BCF_OUTPUT/$bcf_suite/locks/continuations_${bcf_task}.tsv"
  echo "Checkpointed continuation submitted: $bcf_new (same task, after this allocation)."
  exit 0
fi
if [[ "$bcf_rc" == 124 || "$bcf_rc" == 137 ]]; then
  echo "Hard timeout: retained completed chains. Inspect chain logs and runtime before resubmitting; no blind retry loop." >&2
fi
if [[ "$bcf_rc" == 0 && "$bcf_suite" == ihdp && "${BCF_AUTO_IHDP_ANALYSIS:-0}" == 1 ]]; then
  bash "$BCF_CODE_ROOT/slurm/ihdp_completed.sh"
fi
exit "$bcf_rc"
