#!/usr/bin/env bash
# Only dependency checks and tiny actual fits, using the production container
# and R library. These tests do not establish convergence or benchmark runtime.
set -euo pipefail
source "$(dirname -- "${BASH_SOURCE[0]}")/paths.sh"
for bcf_suite in main-linear main-nonlinear ortho; do
  echo "Checking simulation dependencies: $bcf_suite"
  bash "$BCF_CODE_ROOT/slurm/rscript.sh" run.R check "$bcf_suite"
done
for bcf_suite in main-linear main-nonlinear ortho; do
  echo "Running tiny simulation preflight: $bcf_suite"
  bash "$BCF_CODE_ROOT/slurm/rscript.sh" run.R preflight "$bcf_suite"
done
echo 'Normal and fixed-external-propensity orthogonal simulation preflights passed.'
