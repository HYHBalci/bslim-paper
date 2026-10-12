#!/usr/bin/env bash
# Source only. No installation, automatic Git updates, or restricted-data copies.
export BCF_CODE_ROOT="${BCF_CODE_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)}"
[[ -f "$BCF_CODE_ROOT/run.R" ]] || { echo "Submit from the simulations folder." >&2; return 2; }
if [[ -z "${BCF_ACCOUNT:-}" ]]; then
  if [[ "$BCF_CODE_ROOT" =~ ^/(appdata|net/beegfs)/users/([^/]+)/ ]]; then
    BCF_ACCOUNT="${BASH_REMATCH[2]}"
  else
    BCF_ACCOUNT="${USER:-$(id -un)}"; BCF_ACCOUNT="${BCF_ACCOUNT##*\\}"
    if [[ ! -d "/net/beegfs/users/$BCF_ACCOUNT" && -d "/net/beegfs/users/${BCF_ACCOUNT^^}" ]]; then BCF_ACCOUNT="${BCF_ACCOUNT^^}"; fi
  fi
fi
export BCF_ACCOUNT
if [[ -z "${BCF_OUTPUT:-}" ]]; then
  if [[ -d "/net/beegfs/users/$BCF_ACCOUNT" ]]; then
    export BCF_OUTPUT="/net/beegfs/users/$BCF_ACCOUNT/spbart_8.3.5_interactions_v5"
  else
    export BCF_OUTPUT="$BCF_CODE_ROOT/output/interactions-v5"
  fi
fi
[[ "$BCF_OUTPUT" == /* ]] || { echo "BCF_OUTPUT must be absolute" >&2; return 2; }
mkdir -p -- "$BCF_OUTPUT"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
if [[ -n "${BCF_RLIB:-}" ]]; then export R_LIBS_USER="$BCF_RLIB"; fi
