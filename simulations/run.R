#!/usr/bin/env Rscript
# One entry point; launchers determine paths. No package installation or download.
arguments <- commandArgs(trailingOnly=TRUE)
file_arg <- grep("^--file=",commandArgs(),value=TRUE)
ROOT <- dirname(normalizePath(sub("^--file=","",file_arg[1]),winslash="/",mustWork=TRUE))
setwd(ROOT)
source("config.R")
config$root <- ROOT
config$output_root <- Sys.getenv("BCF_OUTPUT",file.path(ROOT,"output","interactions-v5"))
config$ihdp_data_dir <- Sys.getenv("BCF_IHDP_DATA",file.path(ROOT,config$ihdp_data_dir))
rlib <- Sys.getenv("BCF_RLIB")
if(nzchar(rlib)) .libPaths(c(rlib,.libPaths()))
# Child R processes inherit these exact resolved library locations.
Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
for(file in c("common.R","data_generator.R","main.R","metrics.R","orthogonalisation.R","ihdp.R","engine.R","runtime.R"))
  source(file.path("R",file))
if(file.exists("R/ihdp_calibration.R")) source("R/ihdp_calibration.R")
usage <- function() cat(paste0(
 "Usage: Rscript run.R COMMAND [arguments]\n",
 "  check [SUITE]                Dependency, package and input checks; no fits\n",
 "  preflight SUITE              Tiny actual fits; not a convergence study\n",
 "  manifest                     Print exact study sizes\n",
 "  status [SUITE|GROUP]          Completed/missing fits and diagnostic failures\n",
 "  analyse [SUITE|GROUP]         Export new results only; require complete grid\n",
 "  worker SUITE SHARD SHARDS     Internal resumable Slurm worker\n",
 "  batch SUITE TASK              Small simulation batch: one n and replicate\n",
 "  benchmark-task TASK           Full-length Helios runtime batch (18 cases)\n",
 "  runtime-report                Evaluate the fixed layout; never submit jobs\n",
 "  chain REQUEST CHAIN          Internal isolated-chain worker\n",
 "SUITE: main-linear, main-nonlinear, ortho, ihdp.\n",
 "GROUP: normal (both main suites), simulations (normal + ortho), all (also IHDP).\n",
 "From Helios use: bash slurm/rscript.sh run.R COMMAND ...\n"))
if(!length(arguments) || arguments[1] %in% c("help","--help","-h")) { usage(); quit(status=0) }
command <- arguments[1]; arguments <- arguments[-1]
status <- tryCatch({
  switch(command,
    check={ check_environment(if(length(arguments)) arguments[1] else "main-linear",config); 0L },
    preflight={ if(length(arguments)!=1) stop("Specify one suite"); preflight(arguments[1],config); 0L },
    manifest={ print_manifest(config); 0L },
    status={ study_status(if(length(arguments)) arguments[1] else "all",config); 0L },
    analyse={ analyse_studies(if(length(arguments)) arguments[1] else "all",config) },
    worker={ if(length(arguments)!=3) stop("worker requires SUITE SHARD SHARDS")
      worker(arguments[1],as.integer(arguments[2]),as.integer(arguments[3]),config) },
    batch={ if(length(arguments)!=2) stop("batch requires SUITE TASK")
      batch_worker(arguments[1],as.integer(arguments[2]),config) },
    `benchmark-task`={ if(length(arguments)!=1) stop("benchmark-task requires TASK")
      benchmark_task(as.integer(arguments[1]),config); 0L },
    `runtime-report`={ runtime_report(config) },
    chain={ if(length(arguments)!=2) stop("chain requires REQUEST CHAIN")
      run_chain(arguments[1],as.integer(arguments[2]),config); 0L },
    stop("Unknown command: ",command))
},error=function(e) { message("ERROR: ",conditionMessage(e)); 1L })
quit(status=as.integer(status))
