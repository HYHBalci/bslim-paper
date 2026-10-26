#!/usr/bin/env Rscript
# One n=500 data set, four independent sequential chains, fitting calls only.
# Run: Rscript simulations/benchmark_fit.R [output.json]
# BCF_RLIB can identify the library containing the installed spbart build.

args <- commandArgs(trailingOnly = FALSE)
script <- normalizePath(sub("^--file=", "", args[grepl("^--file=", args)]),
                        mustWork = TRUE)
root <- dirname(script)
output_args <- commandArgs(trailingOnly = TRUE)
if (length(output_args) > 1L) stop("Supply at most one output JSON path")
output <- if (length(output_args)) output_args[[1L]] else
  file.path(root, "benchmark_fit.json")
if (file.exists(output)) stop("Output already exists; choose a new path: ", output)
if (!dir.exists(dirname(output))) stop("Output directory does not exist")
rlib <- Sys.getenv("BCF_RLIB")
if (nzchar(rlib)) .libPaths(c(rlib, .libPaths()))
for (file in c("common.R", "data_generator.R", "main.R"))
  source(file.path(root, "R", file))
source(file.path(root, "config.R"))
need_packages(c("bslim", "arm", "jsonlite", "digest"))
check_variance_repair()
stopifnot(config$chains == 4L, config$gfr == 25L,
          config$burnin == 2000L, config$draws == 4000L)

machine_info <- function() {
  if (.Platform$OS.type == "windows") {
    query <- paste(
      "$c=Get-CimInstance Win32_Processor;",
      "$m=Get-CimInstance Win32_PhysicalMemory;",
      "$o=Get-CimInstance Win32_OperatingSystem;",
      "[pscustomobject]@{cpu=($c.Name -join '; ');",
      "physical_cores=($c.NumberOfCores|Measure-Object -Sum).Sum;",
      "logical_processors=($c.NumberOfLogicalProcessors|Measure-Object -Sum).Sum;",
      "installed_ram_bytes=($m.Capacity|Measure-Object -Sum).Sum;",
      "os=$o.Caption;os_version=$o.Version;architecture=$o.OSArchitecture}",
      "|ConvertTo-Json -Compress")
    answer <- system2("powershell.exe",
      c("-NoProfile", "-NonInteractive", "-Command", shQuote(query)),
      stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(answer, "status"))) stop("Could not read machine specifications")
    return(jsonlite::fromJSON(paste(answer, collapse = "\n")))
  }
  cpu <- if (file.exists("/proc/cpuinfo")) {
    line <- grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[1L]
    sub("^[^:]+:\\s*", "", line)
  } else NA_character_
  ram <- if (file.exists("/proc/meminfo")) {
    line <- grep("^MemTotal:", readLines("/proc/meminfo"), value = TRUE)
    as.numeric(sub("^MemTotal:\\s*([0-9]+).*$", "\\1", line)) * 1024
  } else NA_real_
  list(cpu = cpu, physical_cores = parallel::detectCores(logical = FALSE),
       logical_processors = parallel::detectCores(), installed_ram_bytes = ram,
       os = unname(Sys.info()["sysname"]), os_version = unname(Sys.info()["release"]),
       architecture = R.version$arch)
}

# All preparation, package loading, and argument construction precede timing.
hardware <- machine_info()
job <- data.frame(suite = "main-nonlinear", n = 500L, replicate = 1L,
  scenario = "default", heterogeneity = TRUE, mu_form = "nonlinear",
  model = "SemiParam", fold = 0L, stringsAsFactors = FALSE)
data <- prepare_main(job, config)
X <- data$X_train
y <- data$y_train
z <- data$z_binary - .5
pi <- data$pi_train
gfr <- config$gfr
burnin <- config$burnin
draws <- config$draws
fit_function <- bslim::bcf_linear_probit
seeds <- vapply(seq_len(config$chains), function(k) fit_seed(job, k), integer(1))
general_parameters <- lapply(seeds, gp_hchs, cfg = config)
timings <- vector("list", config$chains)
validations <- vector("list", config$chains)
started_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
cat("Benchmark: N=500; linear treatment effect; non-linear prognosis;",
    "four independent sequential chains.\n")

for (chain in seq_len(config$chains)) {
  gp <- general_parameters[[chain]]
  gc(verbose = FALSE)
  cat("Starting chain", chain, "of", config$chains, "\n")
  flush.console()
  # The measured expression contains only the model-fitting function call.
  # Its own preprocessing, sampling, and returned fit object are included.
  elapsed <- system.time(
    fit <- fit_function(X_train = X, y_train = y, Z_train = z,
      propensity_train = pi, num_gfr = gfr, num_burnin = burnin,
      num_mcmc = draws, general_params = gp), gcFirst = FALSE)
  timings[[chain]] <- data.frame(chain = chain, seed = seeds[[chain]],
    elapsed_seconds = unname(elapsed[["elapsed"]]),
    user_cpu_seconds = unname(elapsed[["user.self"]]),
    system_cpu_seconds = unname(elapsed[["sys.self"]]))
  # Validation and disposal are deliberately outside the measured expression.
  stopifnot(identical(dim(fit$Beta)[1:2], c(1L, draws)),
            length(fit$sigma2_global_samples) == draws,
            all(is.finite(fit$Beta)), all(is.finite(fit$sigma2_global_samples)))
  validations[[chain]] <- list(coefficient_draws = dim(fit$Beta),
    interaction_draws = dim(fit$Beta_int),
    residual_variance_draws = length(fit$sigma2_global_samples),
    model_params = fit$model_params)
  cat("Chain", chain, "fitting time:", timings[[chain]]$elapsed_seconds, "seconds\n")
  rm(fit)
}

timings <- do.call(rbind, timings)
stopifnot(nrow(timings) == 4L, all(is.finite(timings$elapsed_seconds)),
          all(timings$elapsed_seconds > 0), !anyDuplicated(timings$seed))
pd <- utils::packageDescription("bslim")
source_files <- c("benchmark_fit.R", "config.R", "R/common.R",
                  "R/data_generator.R", "R/main.R")
record <- list(
  complete = TRUE, started_utc = started_utc,
  finished_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  timing_scope = "Only bslim::bcf_linear_probit calls; excludes data generation, external propensity estimation, package loading, validation, posterior processing, and file output",
  execution = "Four independent chains run sequentially in one R process",
  fit_groups = 1L, job = job, chains = config$chains,
  gfr_per_chain = gfr, burnin_per_chain = burnin, retained_draws_per_chain = draws,
  data_seed = job$replicate, general_parameters = general_parameters,
  chain_timings = timings, total_fitting_seconds = sum(timings$elapsed_seconds),
  hardware = hardware,
  software = list(R = R.version.string, platform = R.version$platform,
    bslim = as.character(utils::packageVersion("bslim")),
    arm = as.character(utils::packageVersion("arm")),
    variance_fix = pd[["Config/VarianceFix"]],
    tree_sampler_fix = pd[["Config/TreeSamplerFix"]],
    BLAS = unname(extSoftVersion()["BLAS"]),
    dll_sha256 = digest::digest(file = getLoadedDLLs()[["bslim"]][["path"]], algo = "sha256"),
    fitting_function_sha256 = hash_object(body(fit_function))),
  source_sha256 = setNames(lapply(source_files, function(path)
    digest::digest(file = file.path(root, path), algo = "sha256")), source_files),
  input_sha256 = hash_object(list(X = X, y = y, z = z, propensity = pi)),
  propensity = data$metadata$propensity_model, validation = validations)
jsonlite::write_json(record, output, pretty = TRUE, auto_unbox = TRUE,
                    digits = NA, null = "null", na = "null")
cat("Total fitting time:", record$total_fitting_seconds, "seconds\n")
cat("Saved:", normalizePath(output), "\n")
