#!/usr/bin/env bash
# Called ONLY after an IHDP R worker exits0. Markers certify completed R grids,
# not Slurm job states. Analysis independently checks every expected result.
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/paths.sh"
[[ "${BCF_AUTO_IHDP_ANALYSIS:-}" == 1 && "${BCF_EASY_RUN_ID:-}" =~ ^[1-9][0-9]*$ ]] || {
  echo 'Automatic IHDP joining requires an all-studies run identity.' >&2; exit 2;
}
bcf_join_shard="${SLURM_ARRAY_TASK_ID:-}"
bcf_join_array="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-}}"
[[ "$bcf_join_shard" =~ ^[1-8]$ && "$bcf_join_array" =~ ^[1-9][0-9]*$ ]] || {
  echo 'Invalid final IHDP array element identity.' >&2; exit 2;
}
[[ "${BCF_OUTPUT##*/}" == "run-$BCF_EASY_RUN_ID" ]] || {
  echo 'IHDP completion output/run identity mismatch.' >&2; exit 2;
}
bcf_join_job="${bcf_join_array}_${bcf_join_shard}"
bcf_join_dir="$BCF_OUTPUT/ihdp/completion/$BCF_EASY_RUN_ID"
mkdir -p -- "$bcf_join_dir"
bcf_join_token=$(mktemp "$bcf_join_dir/.completed.XXXXXX")
trap 'rm -f -- "$bcf_join_token"' EXIT
printf 'interactions-v5\t%s\t%s\t%s\n' "$BCF_EASY_RUN_ID" "$bcf_join_shard" "$bcf_join_job" > "$bcf_join_token"
# Same-filesystem hard links publish a complete immutable record atomically,
# including across BeeGFS nodes. Never overwrite another successful attempt.
bcf_join_marker="$bcf_join_dir/shard_$bcf_join_shard.tsv"
if ! ln "$bcf_join_token" "$bcf_join_marker" 2>/dev/null; then
  [[ -f "$bcf_join_marker" ]] || { echo 'Could not publish IHDP completion marker.' >&2; exit 2; }
fi
# Do not reuse the linked inode: modifying it would corrupt the immutable marker.
rm -f -- "$bcf_join_token"
bcf_join_token=$(mktemp "$bcf_join_dir/.analysis_claim.XXXXXX")
bcf_join_ids=()
bcf_join_missing=0
for bcf_join_i in {1..8}; do
  bcf_join_file="$bcf_join_dir/shard_$bcf_join_i.tsv"
  if [[ ! -f "$bcf_join_file" ]]; then bcf_join_missing=1; continue; fi
  IFS=$'\t' read -r bcf_join_version bcf_join_run bcf_join_index bcf_join_id bcf_join_extra < "$bcf_join_file" || {
    echo "Unreadable IHDP completion marker: $bcf_join_file" >&2; exit 2;
  }
  [[ "$bcf_join_version" == interactions-v5 && "$bcf_join_run" == "$BCF_EASY_RUN_ID" &&
     "$bcf_join_index" == "$bcf_join_i" && "$bcf_join_id" =~ ^[1-9][0-9]*_[1-8]$ &&
     "${bcf_join_id##*_}" == "$bcf_join_i" && -z "$bcf_join_extra" &&
     "$(wc -l < "$bcf_join_file")" -eq 1 ]] || {
    echo "Invalid IHDP completion marker: $bcf_join_file. Inspect before recovery." >&2; exit 2;
  }
  bcf_join_ids+=("$bcf_join_id")
done
if (( bcf_join_missing )); then
  echo "IHDP shard $bcf_join_shard complete; automatic analysis waits for all eight completion markers."
  exit 0
fi
# Persistent, server-atomic exactly-once submission claim. Leave it in place
# after success, sbatch failure, or ambiguous reply; automatic retry could
# duplicate a job already accepted by Slurm. No node-local flock assumption.
printf 'interactions-v5\t%s\twinner=%s\tcompleted=%s\n' "$BCF_EASY_RUN_ID" "$bcf_join_job" "${bcf_join_ids[*]}" > "$bcf_join_token"
if ! ln "$bcf_join_token" "$bcf_join_dir/analysis.claim" 2>/dev/null; then
  [[ -f "$bcf_join_dir/analysis.claim" ]] || { echo 'Could not claim IHDP analysis submission.' >&2; exit 2; }
  echo "IHDP analysis already claimed; inspect $bcf_join_dir/analysis.status. No duplicate submitted."
  exit 0
fi
printf 'SUBMITTING\nWinner: %s\n' "$bcf_join_job" > "$bcf_join_dir/analysis.status"
# Only the current, still-running winner is an afterok dependency. Earlier
# final shard IDs may have aged out of Slurm (MinJobAge) during continuations;
# their immutable markers above prove their R completion instead.
if ! bcf_join_analysis=$(sbatch --parsable --chdir="$BCF_CODE_ROOT" \
  --dependency="afterok:$bcf_join_job" --export=ALL,BCF_AUTOCONTINUE=0 \
  "$BCF_CODE_ROOT/05_analyse.sbatch" ihdp); then
  printf 'SUBMISSION_FAILED_OR_UNCERTAIN\nInspect squeue and logs before manually submitting 05_analyse.sbatch ihdp. Persistent claim retained.\n' > "$bcf_join_dir/analysis.status"
  echo "IHDP analysis submission failed; all fit checkpoints remain. Inspect $bcf_join_dir/analysis.status; do not blindly retry." >&2
  exit 2
fi
if [[ ! "$bcf_join_analysis" =~ ^[1-9][0-9]*(\;[^[:space:]]+)?$ ]]; then
  printf 'SUBMISSION_REPLY_UNCERTAIN\nReply: %s\nInspect squeue before recovery. Persistent claim retained.\n' "$bcf_join_analysis" > "$bcf_join_dir/analysis.status"
  echo 'Unexpected IHDP analysis sbatch reply; inspect squeue before recovery.' >&2
  exit 2
fi
printf 'SUBMITTED_NOT_FINISHED\nAnalysis job: %s\nWinner dependency: %s\n' "$bcf_join_analysis" "$bcf_join_job" > "$bcf_join_dir/analysis.status"
printf 'ihdp-analysis\t-\t%s\n' "$bcf_join_analysis" >> "$BCF_OUTPUT/SUBMITTED_JOBS.tsv"
echo "Submitted IHDP analysis after all eight R shards completed: $bcf_join_analysis (afterok current winner $bcf_join_job)."
