#!/usr/bin/env bash
# Internal entry implementation. Argument2 is exclusively an isolated-test
# storage root; 00_run_all.sbatch never passes it or exposes a storage knob.
set -euo pipefail
bcf_easy_root="${1:?code folder required}"
source "$bcf_easy_root/slurm/easy_paths.sh"
bcf_setup_easy_run "$bcf_easy_root" "${2:-/net/beegfs/users}"
cd "$BCF_CODE_ROOT"
printf 'RUN_ID=%s\nCODE=%s\nOUTPUT=%s\n' "$BCF_EASY_RUN_ID" "$BCF_CODE_ROOT" "$BCF_OUTPUT"
printf 'CHECKING_PACKAGE_AND_INPUTS\nCode: %s\nResults: %s\n' "$BCF_CODE_ROOT" "$BCF_OUTPUT" > "$BCF_OUTPUT/RUN_STATUS.txt"
printf 'stage\toffset\tjob_id\n' > "$BCF_OUTPUT/SUBMITTED_JOBS.tsv"
bcf_easy_failure() {
  local bcf_status=$?
  if (( bcf_status != 0 )); then
    printf 'SETUP_OR_SUBMISSION_FAILED\nExit: %s\nInspect dispatcher log and SUBMITTED_JOBS.tsv. Already submitted jobs remain submitted; do not resubmit blindly.\n' "$bcf_status" > "$BCF_OUTPUT/RUN_STATUS.txt"
  fi
}
trap bcf_easy_failure EXIT
[[ -f SHA256SUMS.txt ]] || { echo 'Package checksum file is missing.' >&2; exit 2; }
command -v sha256sum >/dev/null 2>&1 || { echo 'sha256sum is required to verify the uploaded package.' >&2; exit 2; }
sha256sum --check --quiet SHA256SUMS.txt
# Check IHDP dependencies/public inputs BEFORE any simulation arrays. The
# following helper also runs all three tiny simulation integration checks.
bash "$BCF_CODE_ROOT/slurm/rscript.sh" run.R check ihdp
bash "$BCF_CODE_ROOT/slurm/check_simulations.sh"
bash "$BCF_CODE_ROOT/slurm/launch_simulations.sh" all
printf 'SUBMITTED_NOT_FINISHED\nResults: %s\nSee SUBMITTED_JOBS.tsv and each suite analysis/STATUS.txt for completion and convergence.\n' "$BCF_OUTPUT" > "$BCF_OUTPUT/RUN_STATUS.txt"
echo "All study stages scheduled. Results: $BCF_OUTPUT"
