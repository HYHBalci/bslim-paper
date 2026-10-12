#!/usr/bin/env bash
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/paths.sh"
cd "$BCF_CODE_ROOT"
if [[ "${BCF_NATIVE_R:-0}" == 1 ]]; then
  exec "${BCF_RSCRIPT:-Rscript}" --no-init-file --no-site-file "$@"
fi
bcf_container="${BCF_CONTAINER:-/trinity/shared/apps/ood_rstudio/bioconductor_3.20-R-4.4.2.sif}"
[[ -f "$bcf_container" ]] || { echo "Container missing: $bcf_container" >&2; exit 2; }
if ! command -v apptainer >/dev/null 2>&1 && type module >/dev/null 2>&1; then module load apptainer; fi
command -v apptainer >/dev/null 2>&1 || { echo "Load the apptainer module before submission." >&2; exit 2; }
bcf_binds="$BCF_CODE_ROOT,$BCF_OUTPUT"
for bcf_dir in "/appdata/users/$BCF_ACCOUNT" "/net/beegfs/users/$BCF_ACCOUNT" "${BCF_RLIB:-}" "${BCF_IHDP_DATA:-}"; do
  if [[ -n "$bcf_dir" && -d "$bcf_dir" ]]; then bcf_binds+=",$bcf_dir"; fi
done
if [[ -f /trinity/shared/apps/ood_rstudio/Renviron ]]; then
  bcf_binds+=",/trinity/shared/apps/ood_rstudio/Renviron:/usr/local/lib/R/etc/Renviron.site"
fi
exec apptainer exec --bind "${BCF_BINDPATHS:-$bcf_binds}" "$bcf_container" Rscript --no-init-file --no-site-file "$@"
