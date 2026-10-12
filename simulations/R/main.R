main_jobs <- function(cfg, suite) {
  stopifnot(suite %in% c("main-linear","main-nonlinear"))
  cells <- data.frame(scenario=c("default","default","medical_saturation","mild_saturation"),
                      heterogeneity=c(TRUE,FALSE,TRUE,TRUE))
  x <- expand.grid(n=cfg$n_values, replicate=cfg$replicates,
                   cell=seq_len(nrow(cells)), model=c("SemiParam","BCF"),
                   stringsAsFactors=FALSE)
  x$scenario <- cells$scenario[x$cell]; x$heterogeneity <- cells$heterogeneity[x$cell]
  x$suite <- suite; x$mu_form <- sub("main-","",suite); x$fold <- 0L
  x$cell <- NULL
  x[,c("suite","n","replicate","scenario","heterogeneity","mu_form","model","fold")]
}
prepare_main <- function(job, cfg) {
  d <- generate_data_medical(n=job$n, is_te_hetero=job$heterogeneity,
    is_mu_nonlinear=job$mu_form=="nonlinear", seed=job$replicate,
    scenario=job$scenario, z_diff=.5, BCF=FALSE)
  X <- as.matrix(d[,1:6]); z <- as.numeric(d$z>0)
  propensity <- pairwise_logistic_propensity(X,z)
  pi <- propensity$training_probability
  list(X_train=X,y_train=d$y,z_binary=z,X_eval=X,
       truth_eval=d$tau,y_eval=d$y,z_eval=z,eval_index=seq_len(nrow(X)),
       pi_train=pi,pi_eval=pi,pi_true=d$pi_x,
       metadata=list(propensity_model=propensity$metadata))
}
fit_main_chain <- function(job, data, chain, cfg) {
  f <- if (job$model=="SemiParam") fit_sp_chain else fit_bcf_chain
  f(data$X_train,data$y_train,data$z_binary,data$X_eval,
    data$pi_train,data$pi_eval,fit_seed(job,chain),cfg,z_eval=data$z_eval)
}
