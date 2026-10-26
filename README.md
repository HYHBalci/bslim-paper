# B-SLIM paper — version 8.3.6

Latest manuscript, supplement, and analysis scripts for *A BART Sparse Linear Interactions Model for Heterogeneous Treatment Effects*.

- `manuscript/`: version 8.3.6 main paper and supplement, as PDF and LaTeX, with their bibliography and figures.
- `simulations/`: synthetic simulations, orthogonalisation, and IHDP analyses.
- `canvas/`: CANVAS preparation, external validation, and interpretation scripts.

Compile both documents with `manuscript/build.ps1` on Windows or `bash manuscript/build.sh` on Linux/macOS, using a LaTeX installation with `pdflatex` and `bibtex`.

Install the R package from the source checkout, including its C++ dependencies. Follow the [package installation instructions](https://github.com/HYHBalci/bslim#r-package) to install the required R packages and C++17 toolchain.

```sh
git clone --recurse-submodules https://github.com/HYHBalci/bslim.git
git -C bslim checkout 5ea43eeffcbd8efc387eb485b71f88f854be5e59
git -C bslim submodule update --init --recursive
R CMD INSTALL bslim
```

The analysis scripts use R 4.4.2. Install the additional packages required by the dependency checks in `simulations/run.R`; IHDP preparation also requires Python with NumPy.

To reproduce the fitting-time benchmark, run `Rscript simulations/benchmark_fit.R benchmark-local.json` from the repository root, choosing a new output path if that file already exists. It uses one simulated data set with N=500, linear treatment effects and non-linear prognosis, and the four-chain settings in `simulations/config.R`. The chains run independently and sequentially. Only calls to `bslim::bcf_linear_probit` are timed, including their internal preprocessing and sampling; data generation, external propensity estimation, package loading, validation and file output are excluded. The recorded run in `simulations/benchmark_fit.json`, under the earlier package name `spbart`, took 132.13 seconds in total on an Intel Core i5-1335U (10 cores, 16 GB RAM; Windows 11, R 4.4.2), with single-thread limits for OpenMP and BLAS. `BCF_RLIB` can specify the R library containing `bslim`.

The updated scripts use four independent package fits per Bayesian model, each with 25 grow-from-root initialization iterations, 2000 burn-in iterations, and 4000 retained draws. This applies to synthetic simulations, orthogonalisation, IHDP, and both CANVAS analyses. Observational propensity models use Bayesian logistic regression with main effects and pairwise interactions. The basis is selected using training rows only; redundant columns are recorded. IHDP shares the cross-fitted propensity estimates across methods. CANVAS passes the same constant treatment proportion from the full CANVAS-R training pool to both models in fitting and prediction, and both include it in their prognostic component.

The main synthetic simulation tables have been updated from the full rerun results; the paired orthogonalisation results, IHDP results and CANVAS figures retain their corresponding study outputs. The updated methods and sampling table describe the rerun scripts. Per-quantity and grouped posterior diagnostics allow treatment effects, coefficients, residual variance, and shrinkage scales to be assessed separately; they do not exclude fits from performance summaries. Simulations retain aligned local/global scale and coefficient draws in the compact result objects for joint inspection (`keep_scale_posteriors` in `config.R`).

From `simulations/`, prepare the public IHDP inputs and inspect the study configuration:

```sh
python prepare_ihdp.py
Rscript run.R manifest
Rscript run.R check main-linear
Rscript run.R check ortho
Rscript run.R check ihdp
```

On Helios, submit `sbatch 00_run_all.sbatch` from `simulations/`. The supplied launchers use Slurm, BeeGFS storage, and the configured R 4.4.2 Apptainer image. Install the R dependencies in the library available inside that environment. `BCF_RLIB`, `BCF_CONTAINER`, and `BCF_IHDP_DATA` can specify the library, image, and input directory. The `spbart-8.3.5-interactions-v5` protocol uses fresh output directories and does not reuse the preceding propensity-model checkpoints. The analysis exports include `diagnostic_groups.csv` alongside the per-quantity diagnostics.

CANVAS participant data require authorized access through the [YODA Project](https://yoda.yale.edu/) and are not included. Run the CANVAS scripts in the authorized workspace. `canvas_prepare.R` uses `BCF_WORKSPACE` and `BCF_DATA_ROOT`, the installed `bslim` package, and the supplied `canvas_mcmc.R` and `shapley_aux.R`. No separate `linear.R` is required. The script defaults to existing prepared analysis sets; raw-data rebuilding is controlled in its configuration. From the repository root, run:

```sh
Rscript canvas/canvas_prepare.R
Rscript canvas/new_OOS.R
Rscript canvas/interpretation_all_covariates.R
```

`new_OOS.R` performs external validation; the final script reads the saved four-chain coefficient arrays and encoded design for interpretation summaries. Outputs default to `BCF_WORKSPACE/spbart-8.3.5-rerun-mu-propensity-20261010/{interpretation,out-of-sample}`. Set `BCF_CANVAS_RUN` to another fresh run directory if needed. Existing outputs are protected; `BCF_CANVAS_RESUME=true` permits only matching external-validation caches to be resumed. Set `BCF_KEEP_SCALE_DRAWS=true` before the interpretation fit to retain aligned sampled scales alongside coefficients for joint inspection; this uses additional workspace storage. `BCF_SHAPLEY_SRC` can override the supplied Shapley helper.
