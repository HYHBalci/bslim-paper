frac <- function(num,den) if(den>0) num/den else NA_real_
finite_mean <- function(x) if(any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
postprocess_map <- function(items,fun) {
  # Called only while pooling, after all independent-chain workers have joined.
  # Fork numeric, read-only summaries on Unix; do not create Windows PSOCK
  # copies of the potentially large pooled posterior matrix.
  allocated <- suppressWarnings(as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK","1")))
  if(length(allocated)!=1L || !is.finite(allocated) || allocated<1 || allocated!=floor(allocated))
    allocated <- 1
  cores <- min(4L,allocated,length(items))
  if(.Platform$OS.type!="unix" || cores<2L) return(lapply(items,fun))
  result <- parallel::mclapply(items,fun,mc.cores=cores,mc.preschedule=TRUE,
                             mc.set.seed=FALSE,mc.allow.recursive=FALSE)
  failed <- vapply(result,function(x)inherits(x,"try-error") || is.null(x),logical(1L))
  if(any(failed)) stop("Parallel posterior postprocessing failed: ",
    paste(as.character(result[[which(failed)[1L]]]),collapse=" "))
  result
}

empirical_crps_rows <- function(tau,truth) {
  if(!is.matrix(tau) || !nrow(tau) || !ncol(tau) || length(truth)!=nrow(tau) ||
     any(!is.finite(tau)) || any(!is.finite(truth))) stop("Invalid empirical CRPS inputs")
  # All draws and the same per-person expression are retained. Chunking only
  # distributes rows; no approximation, thinning, RNG or changed summation.
  chunks <- split(seq_len(nrow(tau)),ceiling(seq_len(nrow(tau))/64L))
  values <- postprocess_map(chunks,function(rows) vapply(rows,function(i) {
    s <- sort(tau[i,]); m <- length(s)
    mean(abs(s-truth[i]))-sum((2*seq_len(m)-m-1)*s)/m^2
  },numeric(1L)))
  unname(unlist(values,use.names=FALSE))
}

homogeneous_modifier_metrics <- function(coefficient_summary = NULL) {
  fields <- c("modifier_null_coefficient_count", "modifier_false_positive_count",
              "modifier_false_positive_fraction", "modifier_any_false_positive")
  # BCF has no directly corresponding linear modifier coefficients. Missing is
  # not zero false positives, and the nonzero homogeneous intercept is not null.
  if(is.null(coefficient_summary)) return(setNames(rep(NA_real_,length(fields)),fields))
  required <- c("parameter","lower","upper")
  if(!is.data.frame(coefficient_summary) || !all(required %in% names(coefficient_summary)) ||
     anyNA(coefficient_summary$parameter) || any(!nzchar(coefficient_summary$parameter)) ||
     anyDuplicated(coefficient_summary$parameter) ||
     !"alpha" %in% coefficient_summary$parameter ||
     any(!is.finite(coefficient_summary$lower)) ||
     any(!is.finite(coefficient_summary$upper)) ||
     any(coefficient_summary$lower>coefficient_summary$upper))
    stop("Invalid coefficient intervals for homogeneous-effect false positives")
  null <- coefficient_summary[coefficient_summary$parameter!="alpha",,drop=FALSE]
  flags <- null$lower>0 | null$upper<0
  setNames(c(nrow(null),sum(flags),if(length(flags)) mean(flags) else NA_real_,
             as.numeric(any(flags))),fields)
}
safe_cor <- function(x,y,method="pearson") {
  if(length(x)<3 || sd(x)<1e-10 || sd(y)<1e-10) NA_real_ else cor(x,y,method=method)
}
diagnostic_parameter_group <- function(variable) {
  group <- rep("other",length(variable))
  group[variable=="ATE"] <- "ATE"
  group[grepl("^CATE_",variable)] <- "CATE"
  group[variable=="alpha"] <- "treatment_intercept"
  group[grepl("^beta_",variable)] <- "main_coefficients"
  group[grepl("^interaction_",variable)] <- "interaction_coefficients"
  group[variable=="sigma2"] <- "residual_variance"
  group[variable=="lambda_global"] <- "global_shrinkage"
  group[grepl("^lambda_local_",variable)] <- "local_shrinkage"
  group
}
posterior_diagnostic_groups <- function(diagnostics,cfg) {
  if(!nrow(diagnostics)) return(data.frame())
  required <- c("variable","sd","rhat","ess_bulk","ess_tail","mcse_mean","pass")
  if(!all(required %in% names(diagnostics))) stop("Incomplete per-quantity diagnostics")
  # Reconstruct labels when collecting older saved objects. Never change their
  # per-quantity flags or the legacy all-quantities gate.
  diagnostics$parameter_group <- diagnostic_parameter_group(diagnostics$variable)
  qfinite <- function(x,p) {
    x <- x[is.finite(x)]
    if(length(x)) as.numeric(quantile(x,p,names=FALSE)) else NA_real_
  }
  rows <- lapply(unique(diagnostics$parameter_group),function(group) {
    d <- diagnostics[diagnostics$parameter_group==group,,drop=FALSE]
    relative_mcse <- ifelse(is.finite(d$sd)&d$sd>0,d$mcse_mean/d$sd,NA_real_)
    flagged <- is.na(d$pass)|!d$pass
    data.frame(parameter_group=group,quantities=nrow(d),flagged_quantities=sum(flagged),
      flagged_fraction=mean(flagged),any_quantity_flagged=any(flagged),
      rhat_limit=cfg$rhat_limit,ess_min=cfg$ess_min,
      rhat_flags=sum(!is.finite(d$rhat)|d$rhat>cfg$rhat_limit),
      bulk_ess_flags=sum(!is.finite(d$ess_bulk)|d$ess_bulk<cfg$ess_min),
      tail_ess_flags=sum(!is.finite(d$ess_tail)|d$ess_tail<cfg$ess_min),
      nonfinite_rhat=sum(!is.finite(d$rhat)),
      nonfinite_bulk_ess=sum(!is.finite(d$ess_bulk)),nonfinite_tail_ess=sum(!is.finite(d$ess_tail)),
      zero_posterior_sd=sum(is.finite(d$sd)&d$sd==0),
      rhat_median=qfinite(d$rhat,.5),rhat_p90=qfinite(d$rhat,.9),rhat_max=qfinite(d$rhat,1),
      bulk_ess_min=qfinite(d$ess_bulk,0),bulk_ess_median=qfinite(d$ess_bulk,.5),
      tail_ess_min=qfinite(d$ess_tail,0),tail_ess_median=qfinite(d$ess_tail,.5),
      mcse_mean_median=qfinite(d$mcse_mean,.5),mcse_mean_max=qfinite(d$mcse_mean,1),
      mcse_mean_over_sd_median=qfinite(relative_mcse,.5),mcse_mean_over_sd_max=qfinite(relative_mcse,1))
  })
  do.call(rbind,rows)
}
posterior_diagnostics <- function(chains,cfg) {
  if(is.null(chains[[1]]$diagnostic_draws)) return(data.frame())
  need_packages("posterior")
  p <- ncol(chains[[1]]$diagnostic_draws); n <- nrow(chains[[1]]$tau)
  variables <- c(colnames(chains[[1]]$diagnostic_draws),paste0("CATE_",seq_len(n)))
  # Bound working memory when there are thousands of CATE parameters.
  rows <- postprocess_map(seq.int(1L,length(variables),by=100L),function(start) {
    idx <- seq.int(start,min(start+99L,length(variables)))
    a <- array(NA_real_,c(cfg$draws,length(chains),length(idx)),
      dimnames=list(NULL,NULL,variables[idx]))
    for(k in seq_along(chains)) for(j in seq_along(idx)) {
      v <- idx[j]
      a[,k,j] <- if(v<=p) chains[[k]]$diagnostic_draws[,v] else chains[[k]]$tau[v-p,]
    }
    d <- as.data.frame(posterior::summarise_draws(posterior::as_draws_array(a),
      "mean","sd","rhat","ess_bulk","ess_tail","mcse_mean"))
    d$pass <- length(chains)>=2 & is.finite(d$rhat) & d$rhat<=cfg$rhat_limit &
      is.finite(d$ess_bulk) & d$ess_bulk>=cfg$ess_min &
      is.finite(d$ess_tail) & d$ess_tail>=cfg$ess_min
    d$parameter_group <- diagnostic_parameter_group(d$variable)
    d$mcse_mean_over_sd <- ifelse(is.finite(d$sd)&d$sd>0,d$mcse_mean/d$sd,NA_real_)
    d
  })
  do.call(rbind,rows)
}
individual_metrics <- function(individual) {
  truth <- individual$truth; pred <- individual$tau_mean
  lo <- individual$tau_lower; hi <- individual$tau_upper
  Bayesian <- all(is.finite(lo)) && all(is.finite(hi))
  error <- pred-truth; err_ate <- mean(error)
  ans <- c(est_ate=mean(pred),true_ate=mean(truth),err_ate=err_ate,
    sq_err_ate=err_ate^2,mae_ate=abs(err_ate),rmse_cate=sqrt(mean(error^2)),
    bias_cate=mean(error),spearman=safe_cor(pred,truth,"spearman"),
    r2_true=if(sd(truth)>1e-10) 1-sum(error^2)/sum((truth-mean(truth))^2) else NA_real_,
    calib_slope=if(sd(pred)>1e-10 && sd(truth)>1e-10) unname(coef(lm(truth~pred))[2]) else NA_real_,
    sd_ratio=if(sd(truth)>1e-10) sd(pred)/sd(truth) else NA_real_,
    rmse_y=if(all(is.finite(individual$y_hat))) sqrt(mean((individual$y_hat-individual$y)^2)) else NA_real_)
  for(level in c(95,80,90,50)) {
    suffix <- if(level==95) "" else as.character(level)
    l <- individual[[if(level==95) "tau_lower" else paste0("lo",level)]]
    h <- individual[[if(level==95) "tau_upper" else paste0("hi",level)]]
    ok <- !is.null(l) && length(l)==length(truth) && all(is.finite(l)) && all(is.finite(h))
    ans[paste0("cover",suffix,"_cate")] <- if(ok) mean(l<=truth & h>=truth) else NA_real_
    ans[paste0("len",suffix,"_cate")] <- if(ok) mean(h-l) else NA_real_
  }
  if(Bayesian) {
    sig <- lo>0|hi<0; direction <- ifelse(lo>0,1,ifelse(hi<0,-1,0))
    ans <- c(ans,int_score_cate=mean(hi-lo+40*pmax(0,lo-truth)+40*pmax(0,truth-hi)),
      frac_excl0=mean(sig),n_sig=sum(sig),
      sign_correct=frac(sum(direction[sig]==sign(truth[sig])),sum(sig)),
      type_S=frac(sum(direction[sig]!=sign(truth[sig])),sum(sig)),
      sens_benefit=frac(sum(direction==1 & truth>0),sum(truth>0)),
      sens_harm=frac(sum(direction== -1 & truth<0),sum(truth<0)),
      prec_benefit=frac(sum(direction==1 & truth>0),sum(direction==1)),
      prec_harm=frac(sum(direction== -1 & truth<0),sum(direction== -1)))
  } else ans <- c(ans,setNames(rep(NA_real_,9),c("int_score_cate","frac_excl0","n_sig",
    "sign_correct","type_S","sens_benefit","sens_harm","prec_benefit","prec_harm")))
  for(top in c(10,25,50)) {
    cut <- if(sd(truth)>1e-10) quantile(truth,1-top/100,names=FALSE) else NA_real_
    target <- truth>=cut; selected <- lo>=cut
    ans[paste0("sens_top",top)] <- if(Bayesian && is.finite(cut)) frac(sum(selected & target),sum(target)) else NA_real_
    ans[paste0("prec_top",top)] <- if(Bayesian && is.finite(cut)) frac(sum(selected & target),sum(selected)) else NA_real_
  }
  ans
}
pool_fit <- function(job,data,chains,cfg) {
  Bayesian <- job$model %in% c("SemiParam","BCF","Standard","Orthogonal")
  expected <- if(Bayesian) cfg$draws else 1L
  for(z in chains) strict_draws(z$tau,length(data$truth_eval),expected)
  if(Bayesian && length(chains)!=cfg$chains) stop("Incomplete independent chain collection")
  tau <- do.call(cbind,lapply(chains,`[[`,"tau"))
  probs <- c(.025,.975,.1,.9,.05,.95,.25,.75)
  q <- if(Bayesian) t(apply(tau,1,quantile,probs=probs,names=FALSE)) else matrix(NA_real_,nrow(tau),8)
  yhat <- Reduce(`+`,lapply(chains,function(z) if(is.null(z$y_hat)) rep(NA_real_,nrow(tau)) else z$y_hat))/length(chains)
  ind <- data.frame(index=data$eval_index,tau_mean=rowMeans(tau),tau_lower=q[,1],tau_upper=q[,2],
    lo80=q[,3],hi80=q[,4],lo90=q[,5],hi90=q[,6],lo50=q[,7],hi50=q[,8],
    tau_sd=if(Bayesian) apply(tau,1,sd) else NA_real_,
    p_positive=if(Bayesian) rowMeans(tau>0) else NA_real_,
    truth=data$truth_eval,y_hat=yhat,y=data$y_eval,z=data$z_eval)
  metrics <- individual_metrics(ind)
  ate <- colMeans(tau)
  for(level in c(95,80,90,50)) {
    suffix <- if(level==95) "" else as.character(level)
    a <- (1-level/100)/2
    ci <- if(Bayesian) quantile(ate,c(a,1-a),names=FALSE) else c(NA_real_,NA_real_)
    metrics[paste0("cover",suffix,"_ate")] <- if(Bayesian) as.numeric(ci[1]<=mean(data$truth_eval) && ci[2]>=mean(data$truth_eval)) else NA_real_
    metrics[paste0("len",suffix,"_ate")] <- diff(ci)
    if(level==95) metrics <- c(metrics,ate_lower=ci[1],ate_upper=ci[2])
  }
  if(Bayesian) {
    # Exact empirical CRPS using all retained draws, no posterior thinning.
    crps <- empirical_crps_rows(tau,data$truth_eval)
    metrics["crps_cate"] <- mean(crps)
  } else metrics["crps_cate"] <- NA_real_
  coefficients <- if(!is.null(chains[[1]]$coefficients)) do.call(rbind,lapply(chains,`[[`,"coefficients")) else NULL
  coefficient_summary <- NULL
  if(!is.null(coefficients)) {
    cq <- apply(coefficients,2,quantile,probs=c(.025,.975),names=FALSE)
    coefficient_summary <- data.frame(parameter=colnames(coefficients),mean=colMeans(coefficients),
       lower=cq[1,],upper=cq[2,],excludes_zero=cq[1,]>0|cq[2,]<0)
  }
  # Compute the denominator from the retained coefficient design, not a fixed
  # magic number. Here alpha + six main effects + fourteen interactions gives
  # twenty truly null modifiers when tau=3. BCF reports NA, never a false zero.
  if(!job$heterogeneity && job$suite!="ihdp")
    metrics <- c(metrics,homogeneous_modifier_metrics(coefficient_summary))
  extras <- unique(unlist(lapply(chains,function(z) names(z$extra_metrics))))
  for(name in extras) metrics[name] <- finite_mean(vapply(chains,function(z)
    if(!name %in% names(z$extra_metrics)) NA_real_ else as.numeric(z$extra_metrics[[name]]),numeric(1)))
  if(job$suite=="ortho" &&
     exists("orthogonalisation_pooled_metrics",mode="function")) {
    pooled_extra <- orthogonalisation_pooled_metrics(chains,data)
    metrics[names(pooled_extra)] <- pooled_extra
  }
  diagnostics <- if(Bayesian) posterior_diagnostics(chains,cfg) else data.frame()
  pass <- if(Bayesian) nrow(diagnostics)>0 && all(diagnostics$pass) else TRUE
  compact <- lapply(chains,function(z) list(seed=z$seed,
    diagnostic_draws=if(is.null(z$diagnostic_draws)) NULL else z$diagnostic_draws[,
      intersect(c("ATE","sigma2","lambda_global",
        if(isTRUE(cfg$keep_scale_posteriors))
          grep("^lambda_local_",colnames(z$diagnostic_draws),value=TRUE) else character()),
        colnames(z$diagnostic_draws)),drop=FALSE],
    coefficients=z$coefficients,model_params=z$model_params))
  list(job=job,individual=ind,metrics=metrics,diagnostics=diagnostics,
       tau_draws=if(job$suite=="ihdp" && Bayesian && isTRUE(cfg$keep_ihdp_posteriors)) tau else NULL,
       diagnostic_pass=pass,diagnostic_groups=posterior_diagnostic_groups(diagnostics,cfg),
       coefficient_summary=coefficient_summary,
       compact_chains=compact,preprocessing=data$preprocess,
       orthogonalisation=chains[[1]]$orthogonalisation,
       fixed_external=chains[[1]]$fixed_external,
       extra_metadata=data$metadata,elapsed_seconds=sum(vapply(chains,function(z)
         if(is.null(z$elapsed_seconds)) 0 else z$elapsed_seconds,numeric(1))))
}
bind_rows_base <- function(rows) {
  if(!length(rows)) return(data.frame())
  columns <- unique(unlist(lapply(rows,names)))
  do.call(rbind,lapply(rows,function(x) {
    x <- as.data.frame(x,stringsAsFactors=FALSE)
    for(name in setdiff(columns,names(x))) x[[name]] <- NA
    x[,columns,drop=FALSE]
  }))
}
aggregate_metrics <- function(raw,groups=c("suite","n","scenario","heterogeneity","mu_form","model")) {
  if(!nrow(raw)) return(raw)
  split_key <- do.call(paste,c(raw[groups],sep="|"))
  metrics <- setdiff(names(raw)[vapply(raw,is.numeric,logical(1))],c("n","replicate","fold"))
  rows <- lapply(split(raw,split_key),function(d) {
    row <- as.list(d[1,groups,drop=FALSE]); row$replicates <- nrow(d)
    row$convergence_failures <- sum(!d$diagnostic_pass)
    for(m in metrics) {
      v <- d[[m]]; good <- is.finite(v)
      row[[m]] <- finite_mean(v)
      row[[paste0(m,"_n")]] <- sum(good)
      row[[paste0(m,"_sd")]] <- if(sum(good)>1) sd(v[good]) else NA_real_
      row[[paste0(m,"_mcse")]] <- if(sum(good)>1) sd(v[good])/sqrt(sum(good)) else NA_real_
    }
    if("modifier_null_coefficient_count" %in% names(d)) {
      counts <- d$modifier_null_coefficient_count
      finite <- counts[is.finite(counts)]
      row$modifier_null_coefficient_count_min <- if(length(finite)) min(finite) else NA_real_
      row$modifier_null_coefficient_count_max <- if(length(finite)) max(finite) else NA_real_
    }
    if("sq_err_ate" %in% names(d)) row$rmse_ate_true <- sqrt(mean(d$sq_err_ate))
    if("rmse_cate" %in% names(d)) row$pehe_rms <- sqrt(mean(d$rmse_cate^2))
    as.data.frame(row,stringsAsFactors=FALSE)
  })
  bind_rows_base(rows)
}
