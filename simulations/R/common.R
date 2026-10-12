need_packages <- function(packages) {
  missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if (length(missing)) stop("Missing R packages in THIS job's library: ",paste(missing,collapse=", "),
                           ". Install them before submitting; launchers never install packages.")
}

# One external propensity specification for every analysis: main effects and
# all pairwise products, with Bayesian logistic shrinkage. Select an ordered
# full-rank basis using the training rows only, before fitting or prediction.
pairwise_logistic_propensity <- function(X, z, X_new = X, trim = NULL) {
  need_packages("arm")
  X <- as.matrix(X); X_new <- as.matrix(X_new)
  if (!is.numeric(X) || !is.numeric(X_new) || nrow(X) < 2L || !nrow(X_new) ||
      !ncol(X) || ncol(X_new) != ncol(X) ||
      any(!is.finite(X)) || any(!is.finite(X_new)) ||
      length(z) != nrow(X) || any(!is.finite(z)) || any(!z %in% c(0, 1)) ||
      length(unique(z)) != 2L)
    stop("Propensity fitting requires finite aligned covariates and both binary treatment arms")
  if (!is.null(trim) && (length(trim) != 1L || !is.finite(trim) || trim <= 0 || trim >= .5))
    stop("Propensity trimming must be a single value between zero and one half")
  training_names <- colnames(X)
  if (is.null(training_names)) training_names <- paste0("x", seq_len(ncol(X)))
  if (anyNA(training_names) || any(!nzchar(training_names)) || anyDuplicated(training_names))
    stop("Propensity training columns must have unique nonempty names")
  if (!is.null(colnames(X_new))) {
    if (anyDuplicated(colnames(X_new)) || !setequal(colnames(X_new), training_names))
      stop("Propensity prediction columns do not match the training columns")
    X_new <- X_new[, match(training_names, colnames(X_new)), drop = FALSE]
  }
  colnames(X) <- colnames(X_new) <- training_names
  pairs <- if (ncol(X) > 1L) utils::combn(seq_len(ncol(X)), 2L) else matrix(integer(), 2L, 0L)
  design <- function(A) {
    products <- matrix(0, nrow(A), ncol(pairs))
    if (ncol(pairs)) for (j in seq_len(ncol(pairs)))
      products[, j] <- A[, pairs[1L, j]] * A[, pairs[2L, j]]
    colnames(products) <- if (ncol(pairs))
      paste(training_names[pairs[1L, ]], training_names[pairs[2L, ]], sep = ":") else character()
    cbind("(Intercept)" = 1, A, products)
  }
  training_design <- design(X); prediction_design <- design(X_new)
  if (any(!is.finite(training_design)) || any(!is.finite(prediction_design)))
    stop("Non-finite pairwise propensity design")
  # Normalisation makes the rank tolerance independent of measurement units.
  # The ordered QR retains the intercept and identifiable main effects first.
  column_scale <- sqrt(colSums(training_design^2))
  column_scale[column_scale == 0] <- 1
  basis <- qr(sweep(training_design, 2L, column_scale, "/"), tol = 1e-8, LAPACK = FALSE)
  retained <- sort(basis$pivot[seq_len(basis$rank)])
  if (!1L %in% retained) stop("Propensity design lost its intercept")
  slopes <- setdiff(retained, 1L)
  training_frame <- as.data.frame(training_design[, slopes, drop = FALSE])
  prediction_frame <- as.data.frame(prediction_design[, slopes, drop = FALSE])
  names(training_frame) <- names(prediction_frame) <-
    if (length(slopes)) paste0("basis", seq_along(slopes)) else character()
  training_frame$.treatment <- as.numeric(z)
  formula <- if (length(slopes)) .treatment ~ . else .treatment ~ 1
  fit_warnings <- character()
  fit <- withCallingHandlers(
    arm::bayesglm(formula, family = stats::binomial(link = "logit"),
      data = training_frame, prior.scale = 2.5, prior.df = 1,
      prior.scale.for.intercept = 10, prior.df.for.intercept = 1,
      scaled = TRUE, maxit = 1000L),
    warning = function(warning) {
      fit_warnings <<- unique(c(fit_warnings, conditionMessage(warning)))
    })
  training_probability <- as.numeric(stats::fitted(fit))
  probability <- as.numeric(stats::predict(fit, newdata = prediction_frame, type = "response"))
  if (!isTRUE(fit$converged) || any(!is.finite(stats::coef(fit))) ||
      length(training_probability) != nrow(X) || length(probability) != nrow(X_new) ||
      any(!is.finite(c(training_probability, probability))))
    stop("Bayesian logistic propensity with main effects and pairwise interactions failed")
  raw_range <- range(c(training_probability, probability))
  training_trimmed <- prediction_trimmed <- 0L
  if (is.null(trim)) {
    if (any(c(training_probability, probability) <= 0 | c(training_probability, probability) >= 1))
      stop("Bayesian logistic propensity produced a numerical boundary probability")
  } else {
    training_trimmed <- sum(training_probability < trim | training_probability > 1 - trim)
    prediction_trimmed <- sum(probability < trim | probability > 1 - trim)
    training_probability <- pmax(trim, pmin(1 - trim, training_probability))
    probability <- pmax(trim, pmin(1 - trim, probability))
  }
  list(training_probability = training_probability, probability = probability,
       metadata = list(method = "Bayesian logistic regression with main effects and all pairwise interactions",
         training_columns = training_names, proposed_columns = colnames(training_design),
         retained_columns = colnames(training_design)[retained],
         dropped_columns = colnames(training_design)[setdiff(seq_len(ncol(training_design)), retained)],
         design_rank = basis$rank, training_rows = nrow(X),
         basis = "ordered QR from training rows only; intercept and main effects precede interactions",
         rank_tolerance = 1e-8, prior = "scaled Cauchy(0, 2.5); intercept Cauchy(0, 10)",
         trim = trim, training_trimmed = training_trimmed, prediction_trimmed = prediction_trimmed,
         raw_probability_range = raw_range, warnings = fit_warnings,
         converged = fit$converged, iterations = fit$iter))
}

hash_object <- function(x) digest::digest(x,algo="sha256",serialize=TRUE)
safe_name <- function(x) gsub("[^A-Za-z0-9_.-]","_",x)
job_key <- function(job) paste(safe_name(unlist(job[c("suite","n","replicate","scenario",
                                                "heterogeneity","mu_form","model","fold")])),collapse="_")
fit_seed <- function(job,chain) {
  # Unique across the full fixed protocol, not an order-dependent worker RNG.
  suite <- match(job$suite,c("main-linear","main-nonlinear","ortho","ihdp"))-1L
  ni <- if(job$suite=="ihdp") 0L else match(job$n,c(250,500,750,1000,1500,3000))-1L
  sc <- match(job$scenario,c("default","medical_saturation","mild_saturation","ihdp"))-1L
  model <- match(job$model,c("SemiParam","BCF","Standard","Orthogonal","DR","R","S"))-1L
  mu <- match(job$mu_form,c("linear","nonlinear","ihdp"))-1L
  if(anyNA(c(suite,ni,sc,model,mu)) || chain<1 || chain>9 || job$replicate<1 ||
     job$replicate>200 || job$fold<0 || job$fold>10) stop("Invalid seed key")
  as.integer(1000000 + ((((((((suite*6+ni)*200+job$replicate-1)*4+sc)*2+
     as.integer(job$heterogeneity))*3+mu)*7+model)*11+job$fold)*10)+chain)
}
atomic_rds <- function(object,path,allow_identical=TRUE) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  if(file.exists(path)) {
    if(allow_identical && identical(readRDS(path),object)) return(invisible(path))
    stop("Refusing to overwrite incompatible output: ",path)
  }
  tmp <- tempfile(pattern=paste0(basename(path),".",safe_name(Sys.info()[["nodename"]]),".",Sys.getpid(),"."),
                  tmpdir=dirname(path),fileext=".tmp")
  on.exit(unlink(tmp),add=TRUE)
  saveRDS(object,tmp,compress="gzip",version=3)
  if(file.exists(path)) {
    if(allow_identical && identical(readRDS(path),object)) return(invisible(path))
    stop("Cannot publish output over a different existing file: ",path)
  }
  # Hard-link publication is atomic and NEVER replaces an existing path, on
  # POSIX too. The temporary file is on the same filesystem as its destination.
  if(!suppressWarnings(file.link(tmp,path))) {
    if(allow_identical && file.exists(path) && identical(readRDS(path),object)) return(invisible(path))
    stop("Cannot publish output without replacement: ",path,
         ". The output filesystem must support same-directory hard links.")
  }
  invisible(path)
}
validate_request <- function(request) {
  if(!identical(hash_object(request$data),request$data_hash) ||
     !identical(hash_object(request[setdiff(names(request),"request_hash")]),request$request_hash))
    stop("Prepared request checkpoint is corrupt or has been modified")
  invisible(TRUE)
}
software_identity <- function() {
  need_packages(c("spbart","digest","posterior","arm"))
  loadNamespace("spbart")
  pd <- packageDescription("spbart")
  dll <- getLoadedDLLs()[["spbart"]][["path"]]
  pk <- c("spbart","posterior","arm","digest","SuperLearner","glmnet","gam","xgboost","ranger","grf")
  list(R=R.version.string,platform=R.version$platform,
    version=as.character(pd$Version),commit=pd$RemoteSha,
    variance_fix=pd[["Config/VarianceFix"]],
    orthogonalization=pd[["Config/Orthogonalization"]],
    fixed_external_helpers=setNames(lapply(c(".predict_linear_bcf_fixed_external",
      ".fixed_external_linear_draws", ".validate_external_propensity", ".validate_binary_treatment"),
      function(name) if(exists(name,asNamespace("spbart"),inherits=FALSE))
        hash_object(body(get(name,asNamespace("spbart")))) else NULL),
      c("prediction","draws","propensity_validation","treatment_validation")),
    packages=setNames(vapply(pk,function(p)tryCatch(as.character(packageVersion(p)),
      error=function(e) NA_character_),character(1)),pk),
    dll_sha256=digest::digest(file=dll,algo="sha256"),
    sp_function=hash_object(body(spbart::bcf_linear_probit)),
    bcf_function=hash_object(body(spbart::bcf)))
}
check_variance_repair <- function() {
  ns <- asNamespace("spbart")
  if(!exists(".hchs_variance_parameters",ns,inherits=FALSE))
    stop("Installed spbart does not expose the repaired Gaussian hc-hs variance helper. ",
         "Check the library used INSIDE the Helios container; updating the source folder is not installation.")
  ans <- get(".hchs_variance_parameters",ns)(c(1,2),c(.5),c(.8),c(2,4),3,1,.001)
  expected <- .001+(5+(.5/2)^2+(.8/12)^2)/2
  if(abs(ans$shape-3)>1e-12 || abs(ans$scale-expected)>1e-12)
    stop("Installed Gaussian variance helper fails its algebraic check")
  invisible(TRUE)
}
check_fixed_external_package <- function() {
  pd <- utils::packageDescription("spbart")
  if(!identical(pd[["Config/Orthogonalization"]],"fixed-external-propensity-20260921") ||
     utils::packageVersion("spbart") < "0.1.1.9004" ||
     !exists(".predict_linear_bcf_fixed_external",asNamespace("spbart"),inherits=FALSE))
    stop("Orthogonalisation requires the installed fixed-external-propensity spbart build ",
         "0.1.1.9004 (commit 85d2a845) or a compatible newer build. Reinstall it in the ",
         "R library visible INSIDE the Helios container; a Git pull alone is not installation.")
  invisible(TRUE)
}
gp_hchs <- function(seed,cfg) list(cutpoint_grid_size=100,standardize=TRUE,
  sample_sigma2_global=TRUE,sigma2_global_init=1,sigma2_global_shape=1,
  sigma2_global_scale=.001,propensity_covariate="mu",adaptive_coding=FALSE,
  control_coding_init=-.5,treated_coding_init=.5,random_seed=seed,
  num_chains=1L,keep_burnin=FALSE,keep_gfr=FALSE,keep_every=1L,verbose=FALSE,
  sample_global_prior="hc-hs",unlink=TRUE,gibbs=TRUE,hn_scale=1,
  propensity_seperate="none",interaction_rule="continuous_or_binary",
  standardize_cov=FALSE,simple_prior=FALSE,regularize_ATE=FALSE,
  sigma_residual=0,use_ncp=FALSE,orthogonalize=FALSE,robust=FALSE,p0=0)
strict_draws <- function(x,n,draws=NULL) {
  if(!is.matrix(x) || nrow(x)!=n || !ncol(x) || any(!is.finite(x)) ||
     (!is.null(draws) && ncol(x)!=draws)) stop("Invalid or misaligned posterior draws")
  x
}
flatten_coefficients <- function(x) {
  if(length(dim(x))!=3L || dim(x)[1]!=1L) stop("Expected one independently fitted chain")
  matrix(x[1,,,drop=FALSE],nrow=dim(x)[2],ncol=dim(x)[3])
}
semiparam_draws <- function(fit,X) {
  beta <- flatten_coefficients(fit$Beta)
  alpha <- as.vector(t(fit$alpha))
  bint <- if(is.null(fit$Beta_int)) matrix(numeric(),length(alpha),0L) else flatten_coefficients(fit$Beta_int)
  pairs <- fit$interaction_pairs
  if(is.null(pairs) && !ncol(bint)) pairs <- matrix(integer(),2,0)
  scale <- if(isTRUE(fit$model_params$standardize)) fit$model_params$outcome_scale else 1
  if(ncol(beta)!=ncol(X) || nrow(beta)!=length(alpha) || nrow(bint)!=length(alpha) ||
     !is.matrix(pairs) || nrow(pairs)!=2 || ncol(pairs)!=ncol(bint) ||
     any(!pairs %in% seq_len(ncol(X))) || length(scale)!=1 || !is.finite(scale) || scale<=0)
    stop("Coefficient metadata does not match evaluation design")
  W <- matrix(0,nrow(X),ncol(pairs))
  if(ncol(pairs)) for(j in seq_len(ncol(pairs))) W[,j] <- X[,pairs[1,j]]*X[,pairs[2,j]]
  coefficients <- cbind(alpha,beta,bint)*scale
  colnames(coefficients) <- c("alpha",paste0("beta_",seq_len(ncol(beta))),
    if(ncol(pairs)) paste0("interaction_",pairs[1,],"_",pairs[2,]) else character())
  list(tau=strict_draws(cbind(1,X,W)%*%t(coefficients),nrow(X)),
       coefficients=coefficients,interaction_pairs=pairs)
}
sp_prognosis <- function(fit,X,pi) {
  if(!is.null(fit$model_params$propensity_covariate) &&
     fit$model_params$propensity_covariate %in% c("mu","both")) X <- cbind(X,pi)
  ds <- spbart::createForestDataset(X)
  m <- fit$forests_mu$predict(ds)
  s <- if(isTRUE(fit$model_params$standardize)) fit$model_params$outcome_scale else 1
  b <- if(isTRUE(fit$model_params$standardize)) fit$model_params$outcome_mean else 0
  as.numeric(rowMeans(m)*s+b)
}
extract_sp <- function(fit,X_eval,cfg,z_eval=NULL,pi_eval=NULL) {
  X_eval <- get("preprocessPredictionData",asNamespace("spbart"))(X_eval,fit$train_set_metadata)
  out <- semiparam_draws(fit,X_eval)
  strict_draws(out$tau,nrow(X_eval),cfg$draws)
  sigma <- as.vector(fit$sigma2_global_samples)
  glob <- as.vector(t(fit$Tau_glob))
  if(length(sigma)!=cfg$draws || length(glob)!=cfg$draws || any(!is.finite(glob)|glob<=0))
    stop("Missing residual variance or positive global-scale draws")
  out$diagnostic_draws <- cbind(ATE=colMeans(out$tau),sigma2=sigma,
                              out$coefficients,lambda_global=glob)
  if(!is.null(fit$Tau)) {
    local <- flatten_coefficients(fit$Tau)
    colnames(local) <- paste0("lambda_local_",seq_len(ncol(local)))
    out$diagnostic_draws <- cbind(out$diagnostic_draws,local)
  }
  if(!is.null(fit$theta_samples))
    stop("Internal propensity posterior draws are incompatible with this fixed-external protocol")
  orthogonal <- isTRUE(fit$model_params$orthogonalize)
  if(orthogonal && (!identical(fit$model_params$orthogonalization,"fixed-external") ||
     !identical(fit$model_params$propensity_covariate,"none") ||
     isTRUE(fit$model_params$adaptive_coding) || isTRUE(fit$model_params$internal_propensity_model)))
    stop("Fitted orthogonalisation model does not use fixed external propensity exclusively in treatment coding")
  if(orthogonal && (is.null(z_eval) || length(z_eval)!=nrow(X_eval) ||
     any(!is.finite(z_eval)) || any(!z_eval %in% c(0,1)) ||
     length(pi_eval)!=nrow(X_eval) || any(!is.finite(pi_eval)) || any(pi_eval<0|pi_eval>1)))
    stop("Fixed orthogonalisation predictions require binary treatment and external propensity per evaluation row")
  coding <- if(is.null(z_eval)) NULL else if(orthogonal) z_eval-pi_eval else z_eval-.5
  out$y_hat <- if(is.null(z_eval)) rep(NA_real_,nrow(X_eval)) else
    sp_prognosis(fit,X_eval,pi_eval)+rowMeans(out$tau)*coding
  if(orthogonal) out$fixed_external <- list(
    propensity_train=as.numeric(fit$propensity_train),
    treatment_coding_train=as.numeric(fit$treatment_coding_train),
    propensity_eval=as.numeric(pi_eval),treatment_coding_eval=as.numeric(coding))
  out$model_params <- fit$model_params
  out
}
fit_sp_chain <- function(X_train,y_train,z_binary,X_eval,pi_train,pi_eval,seed,cfg,
                         extra_gp=list(),z_eval=NULL,z_train_override=NULL) {
  gp <- modifyList(gp_hchs(seed,cfg),extra_gp)
  if(!isTRUE(gp$gibbs) || gp$hn_scale!=1 || gp$sample_global_prior!="hc-hs" || !gp$unlink)
    stop("This protocol requires unit-scale unlinked hc-hs Gibbs sampling")
  z <- if(is.null(z_train_override)) z_binary-.5 else z_train_override
  if(isTRUE(gp$orthogonalize)) {
    check_fixed_external_package()
    if(!identical(as.numeric(z),as.numeric(z_binary)))
      stop("Supply raw binary treatment to the fixed external orthogonalisation arm")
  }
  fit <- spbart::bcf_linear_probit(X_train=X_train,y_train=y_train,Z_train=z,
    propensity_train=pi_train,num_gfr=cfg$gfr,num_burnin=cfg$burnin,
    num_mcmc=cfg$draws,general_params=gp)
  if(isTRUE(gp$orthogonalize) &&
     (!identical(as.numeric(fit$propensity_train),as.numeric(pi_train)) ||
      !identical(as.numeric(fit$treatment_coding_train),as.numeric(z_binary-pi_train))))
    stop("The installed sampler altered the external propensity or treatment coding")
  extract_sp(fit,X_eval,cfg,z_eval,pi_eval)
}
fit_bcf_chain <- function(X_train,y_train,z_binary,X_eval,pi_train,pi_eval,seed,cfg,
                          extra_gp=list(),z_eval=NULL) {
  gp <- modifyList(list(cutpoint_grid_size=100,standardize=TRUE,propensity_covariate="mu",
    adaptive_coding=TRUE,random_seed=seed,num_chains=1L,keep_burnin=FALSE,
    keep_gfr=FALSE,keep_every=1L,verbose=FALSE,probit_outcome_model=FALSE),extra_gp)
  fit <- spbart::bcf(X_train=X_train,y_train=y_train,Z_train=z_binary,
    propensity_train=pi_train,num_gfr=cfg$gfr,num_burnin=cfg$burnin,
    num_mcmc=cfg$draws,general_params=gp)
  if(is.null(z_eval)) z_eval <- rep(0,nrow(X_eval))
  pr <- predict(fit,X=X_eval,Z=as.matrix(z_eval),propensity=pi_eval)
  tau <- strict_draws(pr$tau_hat,nrow(X_eval),cfg$draws)
  sig <- as.vector(fit$sigma2_global_samples)
  if(length(sig)!=cfg$draws) stop("Missing BCF variance draws")
  list(tau=tau,y_hat=rowMeans(strict_draws(pr$y_hat,nrow(X_eval),cfg$draws)),
       diagnostic_draws=cbind(ATE=colMeans(tau),sigma2=sig),model_params=fit$model_params)
}
