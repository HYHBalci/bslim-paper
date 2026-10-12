#!/usr/bin/env bash
# Internal functions; the public entry always uses /net/beegfs/users. The second
# function argument exists only so isolated tests can supply temporary storage.
bcf_easy_account() {
  local bcf_path="$1" bcf_name
  if [[ "$bcf_path" =~ ^/(home|appdata/users|net/beegfs/users)/([^/]+)(/|$) ]]; then
    bcf_name="${BASH_REMATCH[2]}"
  else
    bcf_name="${SLURM_JOB_USER:-${USER:-$(id -un)}}"
    bcf_name="${bcf_name##*\\}"
  fi
  [[ "$bcf_name" =~ ^[A-Za-z0-9._-]+$ ]] || { echo 'Cannot resolve a safe Helios account name.' >&2; return 2; }
  printf '%s' "$bcf_name"
}
bcf_setup_easy_run() {
  local bcf_root="$1" bcf_storage="${2:-/net/beegfs/users}" bcf_account_root bcf_parent
  [[ "${SLURM_JOB_ID:-}" =~ ^[1-9][0-9]*$ ]] || { echo 'This entry must run through sbatch.' >&2; return 2; }
  export BCF_CODE_ROOT="$(cd -- "$bcf_root" && pwd -P)"
  [[ -f "$BCF_CODE_ROOT/run.R" ]] || { echo 'The submitted folder has no run.R.' >&2; return 2; }
  export BCF_ACCOUNT="$(bcf_easy_account "$BCF_CODE_ROOT")"
  bcf_account_root="$bcf_storage/$BCF_ACCOUNT"
  if [[ ! -d "$bcf_account_root" && -d "$bcf_storage/${BCF_ACCOUNT^^}" ]]; then
    export BCF_ACCOUNT="${BCF_ACCOUNT^^}"; bcf_account_root="$bcf_storage/$BCF_ACCOUNT"
  fi
  [[ -d "$bcf_account_root" && -w "$bcf_account_root" ]] || {
    echo "Required writable BeeGFS account folder missing: $bcf_account_root. Refusing to store this study in home." >&2
    return 2
  }
  bcf_parent="$bcf_account_root/spbart_8.3.5_interactions_v5"
  mkdir -p -- "$bcf_parent"
  export BCF_EASY_RUN_ID="$SLURM_JOB_ID"
  export BCF_OUTPUT="$bcf_parent/run-$BCF_EASY_RUN_ID"
  # No -p here: a rerun of this same dispatcher job must not dispatch twice.
  mkdir -- "$BCF_OUTPUT" || { echo "Run folder already exists or cannot be created: $BCF_OUTPUT" >&2; return 2; }
  export BCF_AUTO_IHDP_ANALYSIS=1
  export BCF_CONTINUATION=0
  unset BCF_TASK_OFFSET BCF_AUTOCONTINUE
}
