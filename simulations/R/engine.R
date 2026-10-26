suite_names <- c("main-linear","main-nonlinear","ortho","ihdp")
simulation_suites <- c("main-linear","main-nonlinear","ortho")
resolve_suites <- function(suite) {
  if(length(suite)!=1L || is.na(suite)) stop("Specify one suite or group")
  if(suite=="all") return(suite_names)
  if(suite=="simulations") return(simulation_suites)
  if(suite=="normal") return(c("main-linear","main-nonlinear"))
  if(!suite %in% suite_names) stop("Unknown suite or group: ",suite)
  suite
}
get_jobs <- function(suite,cfg) {
  if(!suite %in% suite_names) stop("Unknown suite: ",suite)
  jobs <- if(suite=="ortho") orthogonalisation_jobs(cfg) else if(suite=="ihdp") ihdp_jobs(cfg) else main_jobs(cfg,suite)
  # Keep paired model/arm fits close in the manifest, and distribute sizes.
  jobs <- jobs[order(jobs$replicate,jobs$scenario,jobs$heterogeneity,jobs$mu_form,jobs$fold,jobs$n,jobs$model),]
  rownames(jobs) <- NULL
  if(anyDuplicated(vapply(seq_len(nrow(jobs)),function(i)job_key(jobs[i,]),character(1)))) stop("Duplicate jobs")
  jobs
}
prepare_job <- function(job,cfg) {
  if(job$suite=="ortho") prepare_orthogonalisation(job,cfg) else if(job$suite=="ihdp") ihdp_prepare(job,cfg) else prepare_main(job,cfg)
}
fit_job_chain <- function(job,data,chain,cfg) {
  if(job$suite=="ortho") fit_orthogonalisation_chain(job,data,chain,cfg) else
    if(job$suite=="ihdp") ihdp_fit_chain(job,data,chain,cfg) else fit_main_chain(job,data,chain,cfg)
}
is_bayesian <- function(job) job$model %in% c("SemiParam","BCF","Standard","Orthogonal")
result_path <- function(job,cfg) file.path(cfg$output_root,job$suite,"results",paste0(job_key(job),".rds"))
code_identity <- function(cfg) {
  files <- c("config.R","run.R",sort(list.files(file.path(cfg$root,"R"),pattern="[.]R$",full.names=FALSE)))
  paths <- c(file.path(cfg$root,files[1:2]),file.path(cfg$root,"R",files[-(1:2)]))
  setNames(vapply(paths,function(p)digest::digest(file=p,algo="sha256"),character(1)),files)
}
protocol_identity <- function(cfg) {
  input_manifest <- file.path(cfg$ihdp_data_dir,"input_manifest.csv")
  list(settings=cfg[setdiff(names(cfg),c("root","output_root","ihdp_data_dir"))],
       code=code_identity(cfg),software=software_identity(),
       ihdp_input_manifest=if(file.exists(input_manifest))
         digest::digest(file=input_manifest,algo="sha256") else NULL)
}
check_environment <- function(suite,cfg) {
  if(!suite %in% suite_names) stop("Unknown suite")
  need_packages(c("bslim","arm","digest","posterior"))
  check_variance_repair()
  if(suite=="ihdp") {
    need_packages(c("SuperLearner","gam","glmnet","xgboost","ranger","grf"))
    for(r in cfg$ihdp_replicates) ihdp_load(r,cfg$ihdp_data_dir)
  }
  if(cfg$chains<2 || cfg$chains>9 || cfg$draws<4 || cfg$burnin<1) stop("Invalid MCMC configuration")
  if(!identical(as.integer(cfg$ihdp_folds),10L)) stop("This manuscript protocol requires ten IHDP outer folds")
  if(any(!cfg$replicates %in% 1:200) || any(!cfg$ihdp_replicates %in% 1:100)) stop("Replicate outside protocol grid")
  cat("R libraries:\n",paste(.libPaths(),collapse="\n"),"\n",sep="")
  cat("spbart:",as.character(packageVersion("bslim")),"; Gibbs=TRUE, hn_scale=1\n")
  if(suite=="ortho") {
    check_fixed_external_package()
    cat("Fixed external orthogonalisation: no propensity in the forest and no internal propensity draws.\n")
  }
  invisible(TRUE)
}
preflight <- function(suite,cfg) {
  check_environment(suite,cfg)
  if(suite=="ortho") return(orthogonalisation_preflight(cfg))
  tiny <- cfg; tiny$chains <- 2L; tiny$gfr <- 2L; tiny$burnin <- 3L; tiny$draws <- 6L
  # A small Gaussian integration check, independent of expensive IHDP nuisances.
  jobs <- main_jobs(tiny,if(suite %in% c("main-linear","main-nonlinear")) suite else "main-linear")[1:2,]
  jobs$n <- 60L; jobs$replicate <- 1L; jobs$scenario <- "default"; jobs$heterogeneity <- TRUE
  data <- prepare_main(jobs[1,],tiny)
  for(model in c("SemiParam","BCF")) {
    fun <- if(model=="SemiParam") fit_sp_chain else fit_bcf_chain
    z <- fun(data$X_train,data$y_train,data$z_binary,data$X_eval,
      data$pi_train,data$pi_eval,seed=901L,cfg=tiny,z_eval=data$z_eval)
    strict_draws(z$tau,60,6)
    if(any(!is.finite(z$y_hat))) stop("Non-finite preflight predictions")
  }
  cat("Tiny Gaussian Gibbs/BCF integration fits passed. This is NOT convergence or coverage validation.\n")
  invisible(TRUE)
}
print_manifest <- function(cfg) {
  for(suite in suite_names) {
    jobs <- get_jobs(suite,cfg)
    fits <- sum(ifelse(vapply(seq_len(nrow(jobs)),function(i)is_bayesian(jobs[i,]),logical(1)),cfg$chains,1L))
    cat(sprintf("%-16s %6d fit groups, %7d independent chain/point fits\n",suite,nrow(jobs),fits))
  }
}
run_chain <- function(request_file,chain,cfg) {
  need_packages(c("digest","bslim"))
  request <- readRDS(request_file)
  validate_request(request)
  if(!identical(protocol_identity(cfg),request$identity)) stop("Code/package/config changed after request creation")
  if(request$job$suite=="ihdp") {
    current <- ihdp_load(request$job$replicate,cfg$ihdp_data_dir)
    if(!identical(current$source_sha256,request$data$input$source_sha256)) stop("IHDP input changed after preparation")
  }
  maximum <- if(is_bayesian(request$job)) cfg$chains else 1L
  if(chain<1 || chain>maximum) stop("Invalid chain number")
  path <- file.path(dirname(request_file),sprintf("chain_%02d.rds",chain))
  if(file.exists(path)) {
    z <- readRDS(path)
    if(!identical(z$request_hash,request$request_hash) || z$chain!=chain) stop("Incompatible chain checkpoint")
    cat("Reusing completed chain",chain,"\n"); return(invisible(NULL))
  }
  t0 <- proc.time()[3]
  z <- fit_job_chain(request$job,request$data,chain,cfg)
  z$seed <- fit_seed(request$job,chain); z$chain <- chain
  z$request_hash <- request$request_hash; z$elapsed_seconds <- unname(proc.time()[3]-t0)
  strict_draws(z$tau,length(request$data$truth_eval),if(is_bayesian(request$job)) cfg$draws else 1L)
  if(!is.null(z$diagnostic_draws) && any(!is.finite(z$diagnostic_draws))) stop("Non-finite diagnostic draws")
  atomic_rds(z,path)
  cat("Saved chain",chain,"in",round(z$elapsed_seconds),"seconds\n")
}
run_fit_group <- function(job,cfg,identity) {
  target <- result_path(job,cfg); id <- hash_object(identity)
  if(file.exists(target)) {
    result <- readRDS(target)
    if(!identical(result$protocol_hash,id) || !identical(result$job,job)) stop("Existing fit has incompatible settings/code/package")
    return(invisible(result))
  }
  checkpoint <- file.path(cfg$output_root,job$suite,"checkpoints",job_key(job))
  dir.create(checkpoint,recursive=TRUE,showWarnings=FALSE)
  request_file <- file.path(checkpoint,"request.rds")
  if(file.exists(request_file)) {
    request <- readRDS(request_file)
    validate_request(request)
    if(!identical(request$identity,identity) || !identical(request$job,job)) stop("Incompatible request checkpoint")
    # Frozen prepared inputs are reused; no fresh nuisance randomness on resume.
  } else {
    data <- prepare_job(job,cfg)
    request <- list(job=job,data=data,identity=identity,data_hash=hash_object(data))
    request$request_hash <- hash_object(request)
    atomic_rds(request,request_file)
  }
  nc <- if(is_bayesian(job)) cfg$chains else 1L
  cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","1"))
  launch <- function(k) {
    log <- file.path(checkpoint,sprintf("chain_%02d.log",k))
    status <- system2(file.path(R.home("bin"),"Rscript"),
      c("--no-init-file","--no-site-file",shQuote(file.path(cfg$root,"run.R")),
        "chain",shQuote(request_file),as.character(k)),stdout=log,stderr=log)
    if(status!=0) stop("Chain ",k," failed; see ",log)
    k
  }
  if(.Platform$OS.type=="unix" && cpus>1 && nc>1) {
    launched <- parallel::mclapply(seq_len(nc),launch,mc.cores=min(nc,cpus),mc.set.seed=FALSE)
    if(any(vapply(launched,inherits,logical(1),"try-error"))) stop("A chain failed; completed checkpoints are preserved")
  } else for(k in seq_len(nc)) launch(k)
  chains <- lapply(seq_len(nc),function(k) {
    z <- readRDS(file.path(checkpoint,sprintf("chain_%02d.rds",k)))
    if(!identical(z$request_hash,request$request_hash) || z$chain!=k) stop("Mismatched chain")
    z
  })
  result <- pool_fit(job,request$data,chains,cfg)
  result$protocol_hash <- id; result$request_hash <- request$request_hash
  result$data_hash <- request$data_hash; result$identity <- identity
  atomic_rds(result,target)
  # Remove ONLY this worker's verified generated posterior checkpoints, after
  # reopening and hashing the complete compact result. Never touch input data.
  if(!cfg$keep_chain_posteriors && identical(hash_object(readRDS(target)),hash_object(result))) {
    exact <- file.path(checkpoint,sprintf("chain_%02d.rds",seq_len(nc)))
    for(p in exact) {
      z <- readRDS(p)
      if(!identical(z$request_hash,request$request_hash)) stop("Refusing checkpoint compaction: identity mismatch")
      unlink(p,recursive=FALSE)
    }
  }
  cat(job_key(job),": saved; convergence pass =",result$diagnostic_pass,"\n")
  invisible(result)
}
worker <- function(suite,shard,shards,cfg,selected_override=NULL) {
  if(anyNA(c(shard,shards)) || shard<1 || shard>shards) stop("Invalid shard")
  check_environment(suite,cfg)
  identity <- protocol_identity(cfg)
  manifest <- list(identity=identity,jobs=get_jobs(suite,cfg))
  dir.create(file.path(cfg$output_root,suite),recursive=TRUE,showWarnings=FALSE)
  manifest_file <- file.path(cfg$output_root,suite,"manifest.rds")
  atomic_rds(manifest,manifest_file)
  # Check installed integration once per worker process. Orthogonalisation must
  # fail here, not silently become another model when the build lacks support.
  preflight(suite,cfg)
  jobs <- manifest$jobs
  selected <- which((seq_len(nrow(jobs))-1L)%%shards+1L==shard)
  if(!is.null(selected_override)) {
    if(!length(selected_override) || anyNA(selected_override) ||
       any(selected_override<1L | selected_override>nrow(jobs)) ||
       anyDuplicated(selected_override)) stop("Invalid batch selection")
    selected <- selected_override
  }
  elapsed <- numeric(); start <- Sys.time(); remaining <- FALSE
  for(j in selected) {
    job <- jobs[j,,drop=FALSE]; target <- result_path(job,cfg)
    if(file.exists(target)) {
      old <- readRDS(target)
      if(!identical(old$protocol_hash,hash_object(identity)) || !identical(old$job,job)) stop("Incompatible saved fit")
      next
    }
    used <- as.numeric(difftime(Sys.time(),start,units="secs"))
    reserve <- max(cfg$reserve_seconds,if(length(elapsed)) 1.5*max(tail(elapsed,5)) else 0)
    if(used>0 && used+reserve>=cfg$worker_seconds) { remaining <- TRUE; break }
    t0 <- Sys.time()
    run_fit_group(job,cfg,identity)
    elapsed <- c(elapsed,as.numeric(difftime(Sys.time(),t0,units="secs")))
    gc(verbose=FALSE)
  }
  if(remaining) { cat("TIME_BUDGET: checkpoints saved; continue this shard in a new allocation.\n"); return(75L) }
  cat("Shard",shard,"of",shards,"complete\n"); 0L
}
study_status <- function(suite,cfg) {
  suites <- resolve_suites(suite)
  rows <- lapply(suites,function(s) {
    jobs <- get_jobs(s,cfg)
    paths <- vapply(seq_len(nrow(jobs)),function(i)result_path(jobs[i,,drop=FALSE],cfg),character(1))
    present <- file.exists(paths)
    failed <- sum(vapply(paths[present],function(p)!isTRUE(readRDS(p)$diagnostic_pass),logical(1)))
    data.frame(suite=s,expected=nrow(jobs),complete=sum(present),missing=sum(!present),diagnostic_failures=failed)
  })
  ans <- do.call(rbind,rows); print(ans,row.names=FALSE); invisible(ans)
}
analysis_status <- function(out,label,details=character(),stale=FALSE) {
  files <- list.files(out,pattern="[.](csv|rds)$",recursive=TRUE)
  prior_tables <- setdiff(files,"missing_fits.csv")
  caution <- if(stale && length(prior_tables))
    "STALE TABLES: existing CSV/RDS outputs are preserved but are not valid for this analysis invocation." else
      if(stale) "No tables from this invocation are valid until analysis completes successfully." else character()
  writeLines(c(label,details,caution),file.path(out,"STATUS.txt"))
}

export_main_analysis <- function(raw,out) {
  summary <- aggregate_metrics(raw)
  write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE)
  main_columns <- c("n","scenario","heterogeneity","mu_form","model","replicates",
    "rmse_ate_true","cover_ate","len_ate","rmse_cate","rmse_cate_mcse",
    "spearman","cover_cate","cover_cate_mcse","len_cate","convergence_failures")
  modifier_columns <- grep("^modifier_",names(summary),value=TRUE)
  columns <- unique(c(main_columns,modifier_columns))
  write.csv(summary[,intersect(columns,names(summary)),drop=FALSE],
    file.path(out,"manuscript_metrics.csv"),row.names=FALSE)
  homogeneous <- summary[!summary$heterogeneity,,drop=FALSE]
  if(nrow(homogeneous)) {
    # Report the requested false-positive counts alongside ordinary accuracy
    # and coverage. The _n columns are usable replicate counts, not the number
    # of coefficients; the latter has its own explicit denominator column.
    coverage <- grep("^(cover|len)(50|80|90)?_(ate|cate)(_(n|sd|mcse))?$",
                     names(summary),value=TRUE)
    columns <- unique(c(main_columns,coverage,modifier_columns))
    write.csv(homogeneous[,intersect(columns,names(summary)),drop=FALSE],
      file.path(out,"homogeneous_metrics.csv"),row.names=FALSE)
  }
  invisible(summary)
}

analyse_studies <- function(suite,cfg) {
  suites <- resolve_suites(suite)
  outputs <- setNames(lapply(suites,function(s) {
    out <- file.path(cfg$output_root,s,"analysis")
    dir.create(out,recursive=TRUE,showWarnings=FALSE)
    out
  }),suites)
  # Invalidate every selected suite before any completeness/package checks.
  # On errors or interruption a previous COMPLETE label must never survive.
  for(s in suites) analysis_status(outputs[[s]],"INVALID_ANALYSIS_IN_PROGRESS",
    "Current analysis has not yet validated the saved result grid.",stale=TRUE)
  pending <- suites
  any_problem <- FALSE
  tryCatch({
  need_packages("digest")
  for(s in suites) {
    jobs <- get_jobs(s,cfg); paths <- vapply(seq_len(nrow(jobs)),function(i)result_path(jobs[i,,drop=FALSE],cfg),character(1))
    out <- outputs[[s]]
    missing <- jobs[!file.exists(paths),,drop=FALSE]
    write.csv(missing,file.path(out,"missing_fits.csv"),row.names=FALSE)
    if(nrow(missing)) {
      analysis_status(out,"INCOMPLETE_MISSING_FITS",
        paste("Missing fits:",nrow(missing),"; see missing_fits.csv. No current manuscript tables produced."),stale=TRUE)
      message(s,": ",nrow(missing)," fits missing; no incomplete manuscript table produced.")
      any_problem<-TRUE; pending<-setdiff(pending,s); next
    }
    manifest_file <- file.path(cfg$output_root,s,"manifest.rds")
    if(!file.exists(manifest_file)) stop("Missing immutable manifest: ",manifest_file)
    manifest <- readRDS(manifest_file)
    if(!identical(jobs,manifest$jobs)) stop("Analysis grid differs from saved manifest")
    if(!identical(manifest$identity$settings$protocol,cfg$protocol))
      stop("Analysis protocol differs from the saved manifest; do not mix historical and fixed-external runs")
    expected_hash <- hash_object(manifest$identity)
    raw <- list(); diagnostic_rows <- list(); diagnostic_group_rows <- list(); failures <- 0L
    for(i in seq_along(paths)) {
      r <- readRDS(paths[i])
      if(!identical(r$protocol_hash,expected_hash) || !identical(r$job,jobs[i,,drop=FALSE])) stop("Mixed or mismatched analysis inputs")
      if(s!="ihdp" && !isTRUE(as.logical(r$job$heterogeneity))) {
        # This also recovers the newly exposed denominator from older complete
        # fixed-external results, without refitting or modifying saved results.
        coefficients <- if(r$job$model=="BCF") NULL else r$coefficient_summary
        if(r$job$model!="BCF" && is.null(coefficients))
          stop("Missing saved coefficient intervals for homogeneous-effect analysis")
        modifiers <- homogeneous_modifier_metrics(coefficients)
        r$metrics[names(modifiers)] <- modifiers
      }
      raw[[i]] <- cbind(r$job,as.data.frame(as.list(r$metrics)),diagnostic_pass=r$diagnostic_pass)
      diagnostic_rows[[i]] <- data.frame(job=job_key(r$job),pass=r$diagnostic_pass,
        failing_parameters=if(nrow(r$diagnostics)) sum(!r$diagnostics$pass) else 0,
        max_rhat=if(nrow(r$diagnostics) && any(is.finite(r$diagnostics$rhat))) max(r$diagnostics$rhat,na.rm=TRUE) else NA_real_,
        gate_scope="all_monitored_quantities")
      grouped <- posterior_diagnostic_groups(r$diagnostics,cfg)
      if(nrow(grouped)) diagnostic_group_rows[[length(diagnostic_group_rows)+1L]] <-
        cbind(job=job_key(r$job),r$job[rep(1L,nrow(grouped)),,drop=FALSE],grouped)
      failures <- failures+!r$diagnostic_pass
    }
    raw <- bind_rows_base(raw)
    write.csv(raw,file.path(out,"per_fit_metrics.csv"),row.names=FALSE)
    write.csv(bind_rows_base(diagnostic_rows),file.path(out,"convergence.csv"),row.names=FALSE)
    write.csv(bind_rows_base(diagnostic_group_rows),file.path(out,"diagnostic_groups.csv"),row.names=FALSE)
    if(s=="ihdp") {
      # Aggregate one replicate at a time: never hold all saved draws in memory.
      combined <- list()
      for(rep in cfg$ihdp_replicates) {
        selected <- which(jobs$replicate==rep)
        rs <- lapply(paths[selected],readRDS)
        cfg_rep <- cfg; cfg_rep$ihdp_replicates <- rep
        combined[[as.character(rep)]] <- ihdp_aggregate(rs,cfg_rep)
      }
      saveRDS(combined,file.path(out,"ihdp_replicate_summaries.rds"))
      # Module supplies its study-specific table collection.
      export_ihdp_analysis(combined,out,cfg)
    } else {
      export_main_analysis(raw,out)
      if(s=="ortho") export_ortho_analysis(raw,out)
    }
    label <- if(failures) "DRAFT_CONVERGENCE_FAILED" else "COMPLETE_DIAGNOSTICS_PASSED"
    analysis_status(out,label,c(paste("Fits:",nrow(jobs)),paste("Fits with at least one flagged quantity:",failures),
      "The legacy status and pass field apply to all monitored quantities; a flag is not a diagnosis of every estimand.",
      "Inspect diagnostic_groups.csv and the saved per-quantity diagnostics for location and severity of flags.",
      "All requested replicates included. No replicate was discarded for a diagnostic flag.",
      paste("Protocol hash:",expected_hash)))
    pending <- setdiff(pending,s)
    if(failures) any_problem<-TRUE
    cat(s,": exported",nrow(jobs),"fits to",out,";",label,"\n")
  }
  if(any_problem) 2L else 0L
  },error=function(error) {
    for(s in pending) analysis_status(outputs[[s]],"INVALID_ANALYSIS",
      paste("Analysis failed:",conditionMessage(error)),stale=TRUE)
    stop(error)
  })
}
export_ortho_analysis <- function(raw,out) {
  a <- raw[raw$model=="Standard",]; b <- raw[raw$model=="Orthogonal",]
  keys <- c("n","replicate","mu_form")
  pairs <- merge(a,b,by=keys,suffixes=c("_std","_ortho"),all=TRUE)
  if(nrow(pairs)!=nrow(a) || anyNA(pairs$rmse_cate_std) || anyNA(pairs$rmse_cate_ortho)) stop("Incomplete paired orthogonalisation grid")
  pairs$difference <- pairs$rmse_cate_ortho-pairs$rmse_cate_std
  tab <- lapply(split(pairs,paste(pairs$n,pairs$mu_form)),function(d) data.frame(
    n=d$n[1],mu_form=d$mu_form[1],replicates=nrow(d),rmse_std=mean(d$rmse_cate_std),
    rmse_ortho=mean(d$rmse_cate_ortho),difference=mean(d$difference),
    difference_mcse=sd(d$difference)/sqrt(nrow(d)),cover_std=mean(d$cover_cate_std),
    cover_ortho=mean(d$cover_cate_ortho),len_std=mean(d$len_cate_std),len_ortho=mean(d$len_cate_ortho)))
  write.csv(do.call(rbind,tab),file.path(out,"ortho_paired.csv"),row.names=FALSE)
  write.csv(pairs,file.path(out,"ortho_paired_replicates.csv"),row.names=FALSE)
  # External scores are shared by the two paired arms: report each generated
  # dataset only once, without fabricating propensity posterior intervals.
  propensity_columns <- grep("^pi_",names(a),value=TRUE)
  if(length(propensity_columns)) {
    propensity <- a[,c(keys,propensity_columns),drop=FALSE]
    write.csv(propensity,file.path(out,"ortho_propensity_replicates.csv"),row.names=FALSE)
    summary <- aggregate_metrics(a,groups=c("n","mu_form"))
    columns <- c("n","mu_form","replicates",grep("^pi_",names(summary),value=TRUE))
    write.csv(summary[,columns,drop=FALSE],file.path(out,"ortho_propensity_summary.csv"),row.names=FALSE)
  }
}
