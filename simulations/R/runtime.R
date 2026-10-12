# Small, exhaustive simulation batches and an optional fixed-resource runtime report.
# Timing is hardware-specific. A padded compute forecast is NOT a queue guarantee.
simulation_batch_grid <- function(cfg) {
  expand.grid(n=cfg$n_values,replicate=cfg$replicates,KEEP.OUT.ATTRS=FALSE)
}
simulation_batch_selection <- function(suite,task,cfg) {
  if(!suite %in% simulation_suites) stop("Only main/ortho use simulation batches")
  grid <- simulation_batch_grid(cfg)
  if(length(task)!=1L || !is.finite(task) || task!=as.integer(task) ||
     task<1L || task>nrow(grid)) stop("Invalid simulation batch ID")
  jobs <- get_jobs(suite,cfg); cell <- grid[task,]
  which(jobs$n==cell$n & jobs$replicate==cell$replicate)
}
batch_worker <- function(suite,task,cfg) {
  selected <- simulation_batch_selection(suite,task,cfg)
  cat("Small simulation batch",task,"of",nrow(simulation_batch_grid(cfg)),
      "; fit groups:",length(selected),"\n")
  worker(suite,1L,1L,cfg,selected_override=selected)
}
runtime_case_grid <- function(cfg) {
  x <- expand.grid(n=cfg$n_values,suite=simulation_suites,stringsAsFactors=FALSE,
                   KEEP.OUT.ATTRS=FALSE)
  x$task <- seq_len(nrow(x)); x
}
runtime_benchmark_dir <- function(cfg) file.path(cfg$output_root,"runtime-benchmark")
runtime_record_path <- function(task,cfg)
  file.path(runtime_benchmark_dir(cfg),sprintf("timing_%02d.rds",task))
validate_runtime_record <- function(record,cell,identity,cfg) {
  count <- if(cell$suite=="ortho") 4L else 8L
  valid <- isTRUE(record$complete) && identical(record$identity,identity) &&
    identical(record$case,cell) && identical(as.integer(record$chains),as.integer(cfg$chains)) &&
    identical(as.integer(c(record$gfr,record$burnin,record$draws)),
              as.integer(c(cfg$gfr,cfg$burnin,cfg$draws))) &&
    identical(record$host$system,"Linux") &&
    isTRUE(suppressWarnings(as.integer(record$host$allocated_cpus))>=cfg$chains) &&
    is.character(record$host$slurm_job) && length(record$host$slurm_job)==1L &&
    nzchar(record$host$slurm_job) && identical(as.integer(record$fit_groups),count) &&
    length(record$fit_seconds)==count && all(is.finite(record$fit_seconds)) &&
    all(record$fit_seconds>0) && length(record$elapsed_seconds)==1L &&
    is.finite(record$elapsed_seconds) && record$elapsed_seconds>=sum(record$fit_seconds)*.99 &&
    length(record$diagnostic_passes)==count && !anyNA(record$diagnostic_passes)
  if(!isTRUE(valid)) stop("Incompatible/non-Helios or incomplete timing record ",cell$task)
  invisible(TRUE)
}
runtime_host <- function() {
  cpu <- Sys.getenv("PROCESSOR_IDENTIFIER","")
  if(file.exists("/proc/cpuinfo")) {
    lines <- readLines("/proc/cpuinfo",warn=FALSE)
    model <- grep("^model name",lines,value=TRUE)
    if(length(model)) cpu <- sub("^[^:]*:[ ]*","",model[1])
  }
  list(host=unname(Sys.info()["nodename"]),system=unname(Sys.info()["sysname"]),
       cpu=cpu,allocated_cpus=Sys.getenv("SLURM_CPUS_PER_TASK",""),
       slurm_job=Sys.getenv("SLURM_JOB_ID",""),R=R.version.string)
}
benchmark_task <- function(task,cfg) {
  cases <- runtime_case_grid(cfg)
  if(length(task)!=1L || !is.finite(task) || task!=as.integer(task) ||
     task<1L || task>nrow(cases)) stop("Invalid runtime benchmark task")
  if(.Platform$OS.type!="unix" || !nzchar(Sys.getenv("SLURM_JOB_ID")) ||
     as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","0"))<cfg$chains)
    stop("Helios benchmark needs Unix and one allocated CPU per independent chain")
  if(!identical(as.integer(cfg$chains),4L)) stop("Launchers require four independent chains")
  started <- Sys.time()
  cell <- cases[task,]; check_environment(cell$suite,cfg)
  identity <- protocol_identity(cfg); path <- runtime_record_path(task,cfg)
  if(file.exists(path)) {
    record <- readRDS(path)
    validate_runtime_record(record,cell,identity,cfg)
    cat("Reusing complete full-length benchmark",task,"; seconds:",record$elapsed_seconds,"\n")
    return(invisible(record))
  }
  preflight(cell$suite,cfg)
  dir.create(runtime_benchmark_dir(cfg),recursive=TRUE,showWarnings=FALSE)
  # Every interrupted attempt gets a NEW output directory. Timing a partial
  # restart would understate full-batch cost. Valid completed records may be reused.
  attempt <- tempfile(sprintf("attempt_%02d_",task),tmpdir=runtime_benchmark_dir(cfg))
  dir.create(attempt); bench_cfg <- cfg; bench_cfg$output_root <- attempt
  jobs <- get_jobs(cell$suite,cfg)
  selected <- which(jobs$n==cell$n & jobs$replicate==cfg$replicates[1])
  if(!length(selected)) stop("Empty timing batch")
  seconds <- numeric(length(selected))
  passes <- logical(length(selected)); diagnostic_rows <- vector("list",length(selected))
  for(i in seq_along(selected)) {
    before <- proc.time()[3]
    result <- run_fit_group(jobs[selected[i],,drop=FALSE],bench_cfg,identity)
    gc(verbose=FALSE)
    seconds[i] <- unname(proc.time()[3]-before)
    passes[i] <- isTRUE(result$diagnostic_pass)
    d <- result$diagnostics
    finite_extreme <- function(x,fun) if(any(is.finite(x))) fun(x[is.finite(x)]) else NA_real_
    diagnostic_rows[[i]] <- data.frame(job=job_key(jobs[selected[i],,drop=FALSE]),
      pass=passes[i],failing_parameters=sum(!d$pass),
      max_rhat=finite_extreme(d$rhat,max),min_bulk_ess=finite_extreme(d$ess_bulk,min),
      min_tail_ess=finite_extreme(d$ess_tail,min))
  }
  record <- list(complete=TRUE,case=cell,identity=identity,host=runtime_host(),
    started=as.numeric(started),finished=as.numeric(Sys.time()),
    elapsed_seconds=as.numeric(difftime(Sys.time(),started,units="secs")),
    fit_seconds=seconds,fit_groups=length(selected),diagnostic_passes=passes,
    diagnostic_summary=do.call(rbind,diagnostic_rows),
    chains=cfg$chains,gfr=cfg$gfr,burnin=cfg$burnin,draws=cfg$draws,
    attempt=attempt)
  validate_runtime_record(record,cell,identity,cfg)
  atomic_rds(record,path)
  cat("Measured full-length batch:",cell$suite,"n=",cell$n,
      "groups=",length(selected),"seconds=",round(record$elapsed_seconds,1),"\n")
  invisible(record)
}
queue_free_schedule_seconds <- function(costs,slots=20L) {
  load <- numeric(slots)
  for(cost in costs) { next_slot <- which.min(load); load[next_slot] <- load[next_slot]+cost }
  max(load)
}
fixed_runtime_forecast <- function(timings,cfg) {
  needed <- c("suite","n","elapsed_seconds")
  if(!all(needed %in% names(timings))) stop("Missing benchmark columns")
  expected <- runtime_case_grid(cfg)
  key <- function(x) paste(x$suite,x$n,sep="|")
  if(anyDuplicated(key(timings)) || !setequal(key(timings),key(expected)) ||
     any(!is.finite(timings$elapsed_seconds) | timings$elapsed_seconds<=0))
    stop("Incomplete, duplicate or invalid full-length benchmark grid")
  grid <- simulation_batch_grid(cfg)
  if(nrow(grid)!=1200L || length(cfg$n_values)!=6L || cfg$chains!=4L)
    stop("Fixed layout requires six sample sizes, 200 replicates and four chains")
  rows <- list(); batch_limit <- cfg$worker_seconds-cfg$reserve_seconds
  for(suite in simulation_suites) for(offset in c(0L,600L)) {
    task_n <- grid$n[offset+seq_len(600L)]
    measured <- timings$elapsed_seconds[match(paste(suite,task_n,sep="|"),key(timings))]
    padded <- 2*measured+60
    rows[[length(rows)+1L]] <- data.frame(suite=suite,offset=offset,tasks=600L,
      concurrent_tasks=20L,cpus_per_task=4L,maximum_cpus=80L,
      measured_work_seconds=sum(measured),max_measured_batch_seconds=max(measured),
      measured_queue_free_seconds=queue_free_schedule_seconds(measured+60),
      padded_queue_free_seconds=queue_free_schedule_seconds(padded),
      # List-scheduling upper bound, conditional on twenty continuously
      # available identical allocations. It is not a scheduler guarantee.
      padded_compute_bound_seconds=sum(padded)/20+max(padded),
      max_padded_batch_seconds=max(padded),batch_within_budget=max(padded)<batch_limit)
  }
  report <- do.call(rbind,rows)
  estimated <- max(report$padded_compute_bound_seconds)+3600
  list(report=report,maximum_cpus=480L,analysis_reserve_seconds=3600,
    padded_total_seconds=estimated,
    timing_within_12h=estimated<=43200 && all(report$batch_within_budget))
}
runtime_report <- function(cfg) {
  out <- runtime_benchmark_dir(cfg); dir.create(out,recursive=TRUE,showWarnings=FALSE)
  status_path <- file.path(out,"RUNTIME_STATUS.txt")
  writeLines("INVALID_OR_IN_PROGRESS: any older runtime report is stale",status_path)
  identity <- protocol_identity(cfg); cases <- runtime_case_grid(cfg)
  records <- lapply(cases$task,function(task) {
    path <- runtime_record_path(task,cfg)
    if(!file.exists(path)) stop("Missing full-length Helios timing task ",task)
    r <- readRDS(path)
    validate_runtime_record(r,cases[task,,drop=FALSE],identity,cfg)
    r
  })
  timings <- do.call(rbind,lapply(records,function(r) data.frame(suite=r$case$suite,n=r$case$n,
    elapsed_seconds=r$elapsed_seconds,fit_groups=r$fit_groups,host=r$host$host,cpu=r$host$cpu,
    convergence_failures=sum(!r$diagnostic_passes))))
  write.csv(timings,file.path(out,"measured_timings.csv"),row.names=FALSE)
  diagnostic_rows <- lapply(records,function(r) {
    if(is.null(r$diagnostic_summary)) return(NULL)
    cbind(suite=r$case$suite,n=r$case$n,r$diagnostic_summary)
  })
  if(any(lengths(diagnostic_rows)>0)) write.csv(do.call(rbind,diagnostic_rows),
    file.path(out,"pilot_convergence.csv"),row.names=FALSE)
  forecast <- fixed_runtime_forecast(timings,cfg)
  write.csv(forecast$report,file.path(out,"runtime_forecast.csv"),row.names=FALSE)
  failures <- sum(timings$convergence_failures)
  lines <- c(
    if(forecast$timing_within_12h) "TIMING_ESTIMATE_WITHIN_12H_NOT_A_GUARANTEE" else "TIMING_ESTIMATE_EXCEEDS_12H_OR_BATCH_BUDGET",
    if(failures) "PILOT_CONVERGENCE_FAILED: do not treat these settings as publication-ready" else "PILOT_DIAGNOSTICS_PASSED",
    "Fixed settings: four CPUs and 16 GiB per task; twenty concurrent tasks per array; six arrays; maximum 480 CPUs.",
    paste("Padded estimated production plus analysis hours:",round(forecast$padded_total_seconds/3600,2)),
    paste("Pilot fit groups failing convergence:",failures),
    "Two times measured batch time plus sixty seconds per task; unmeasured one-hour analysis allowance.",
    "One full-length replicate per sample size/suite is not a statistical runtime upper bound.",
    "The optional calibration run and all scheduler queue waits are excluded from this estimate.",
    "This report never submits jobs, changes CPU settings, or changes the MCMC budget.",
    "Inspect the full pilot diagnostics before deciding to launch the production study.")
  writeLines(lines,status_path)
  print(forecast$report,row.names=FALSE); cat(paste(lines,collapse="\n"),"\n")
  if(forecast$timing_within_12h && failures==0) 0L else 2L
}
