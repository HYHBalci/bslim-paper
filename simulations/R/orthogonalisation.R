# Paired standard / fixed external-propensity orthogonalisation study.
# Both arms exclude propensity from the prognostic forest. The orthogonal arm
# supplies raw binary treatment and freezes externally estimated pi throughout
# all chains; spbart >= 0.1.1.9004 performs the coding Z - pi internally.

orthogonalisation_jobs <- function(cfg) {
  jobs <- expand.grid(
    n = as.integer(cfg$n_values), replicate = as.integer(cfg$replicates),
    mu_form = c("linear", "nonlinear"), model = c("Standard", "Orthogonal"),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  jobs$suite <- "ortho"
  jobs$scenario <- "default"
  jobs$heterogeneity <- TRUE
  jobs$fold <- 0L
  jobs[, c("suite", "n", "replicate", "scenario", "heterogeneity", "mu_form", "model", "fold")]
}

orthogonalisation_auc <- function(score, label) {
  if (length(score) != length(label) || any(!is.finite(score))) return(NA_real_)
  positive <- score[label == 1]
  negative <- score[label == 0]
  if (!length(positive) || !length(negative)) return(NA_real_)
  ranks <- rank(c(positive, negative), ties.method = "average")
  (sum(ranks[seq_along(positive)]) - length(positive) * (length(positive) + 1) / 2) /
    (length(positive) * length(negative))
}

orthogonalisation_propensity_metrics <- function(probability, truth, z, prefix) {
  valid <- length(probability) == length(truth) && length(z) == length(truth) &&
    all(is.finite(probability)) && all(is.finite(truth))
  values <- if (valid) c(
    rmse = sqrt(mean((probability - truth)^2)),
    cor = suppressWarnings(stats::cor(probability, truth)),
    auc = orthogonalisation_auc(probability, z)) else
      c(rmse = NA_real_, cor = NA_real_, auc = NA_real_)
  stats::setNames(as.list(values), paste0(prefix, "_", names(values)))
}

orthogonalisation_external_propensity <- function(X, z) {
  pairwise_logistic_propensity(X, z)
}

prepare_orthogonalisation <- function(job, cfg) {
  if (!job$model %in% c("Standard", "Orthogonal") ||
      !job$mu_form %in% c("linear", "nonlinear") ||
      job$scenario != "default" || !isTRUE(as.logical(job$heterogeneity)))
    stop("Invalid orthogonalisation job: use the paired heterogeneous default DGP")
  generated <- generate_data_medical(
    n = as.integer(job$n), is_te_hetero = TRUE,
    is_mu_nonlinear = identical(as.character(job$mu_form), "nonlinear"),
    seed = as.integer(job$replicate), RCT = FALSE, scenario = "default",
    z_diff = 0.5, BCF = FALSE, sigma_sq = 1)
  X <- as.matrix(generated[, c("x1", "x2", "x3", "x4", "x5_1", "x5_2")])
  storage.mode(X) <- "double"
  z_binary <- as.numeric(generated$z > 0)
  propensity <- orthogonalisation_external_propensity(X, z_binary)
  metrics <- orthogonalisation_propensity_metrics(
    propensity$training_probability, generated$pi_x, z_binary, "pi_int")
  list(
    X_train = X, y_train = as.numeric(generated$y), z_binary = z_binary,
    X_eval = X, truth_eval = as.numeric(generated$tau),
    y_eval = as.numeric(generated$y), z_eval = z_binary,
    eval_index = seq_len(nrow(X)),
    pi_train = propensity$training_probability, pi_eval = propensity$training_probability,
    propensity_truth = as.numeric(generated$pi_x),
    propensity_metrics = metrics,
    propensity_diagnostic_error = NULL,
    noiseless_eval = as.numeric(generated$y_hat),
    metadata = list(propensity_model = propensity$metadata,
      generator = "historical paired default DGP", outcome_treatment_coding = "Z - 0.5",
      orthogonalisation = "fixed-external", propensity_used_for_coding = "external main-effect and pairwise-interaction bayesglm",
      propensity_in_prognostic_forest = FALSE, propensity_fixed_across_chains = TRUE,
      propensity_diagnostic_models = "external main-effect and pairwise-interaction bayesglm"),
    dgp_metadata = list(scenario = "default", heterogeneity = TRUE,
      mu_form = as.character(job$mu_form), generator_seed = as.integer(job$replicate),
      outcome_treatment_coding = "Z - 0.5", sigma_sq = 1,
      historical_pairing = "same generated data in the two arms"))
}

fit_orthogonalisation_chain <- function(job, data, chain, cfg) {
  orthogonal <- identical(as.character(job$model), "Orthogonal")
  if (!job$model %in% c("Standard", "Orthogonal")) stop("Unknown orthogonalisation arm")
  result <- fit_sp_chain(
    X_train = data$X_train, y_train = data$y_train, z_binary = data$z_binary,
    X_eval = data$X_eval, pi_train = data$pi_train, pi_eval = data$pi_eval,
    seed = fit_seed(job, chain), cfg = cfg,
    extra_gp = list(propensity_covariate = "none", orthogonalize = orthogonal,
                    adaptive_coding = FALSE, save_output = FALSE),
    z_eval = data$z_eval,
    # The fixed-external API requires raw binary treatment and subtracts pi
    # itself. The standard arm retains the helper's centred coding Z - 0.5.
    z_train_override = if (orthogonal) data$z_binary else NULL)
  result$propensity_metrics <- data$propensity_metrics
  result$propensity_diagnostic_error <- data$propensity_diagnostic_error
  result$orthogonalisation <- list(
    arm = as.character(job$model),
    propensity_method = "external main-effect and pairwise-interaction bayesglm; fixed across all chains",
    supplied_treatment_coding = if (orthogonal) "binary 0/1" else "Z - 0.5",
    fitted_treatment_coding = if (orthogonal) "Z - external propensity" else "Z - 0.5",
    propensity_covariate = "none",
    # Store the input identity even for Standard, where pi is not used by the
    # outcome model, to validate the paired study and pooled diagnostics.
    propensity_train = as.numeric(data$pi_train),
    propensity_eval = as.numeric(data$pi_eval))
  if (!is.null(result$propensity_draws) || !is.null(result$theta_samples))
    stop("The fixed-external study must not contain internal propensity draws")
  if (orthogonal) {
    orthogonalisation_validate_fixed_result(result, data)
  }
  result$extra_metrics <- result$propensity_metrics
  result
}

orthogonalisation_validate_fixed_result <- function(result, data) {
  params <- result$model_params
  if (!identical(params$orthogonalization, "fixed-external") ||
      !isTRUE(params$orthogonalize) ||
      !identical(params$propensity_covariate, "none") ||
      !identical(params$adaptive_coding, FALSE) ||
      !identical(params$internal_propensity_model, FALSE))
    stop("Orthogonalisation requires verified fixed-external package metadata")
  expected <- list(propensity_train = as.numeric(data$pi_train),
                   treatment_coding_train = as.numeric(data$z_binary - data$pi_train),
                   propensity_eval = as.numeric(data$pi_eval),
                   treatment_coding_eval = as.numeric(data$z_eval - data$pi_eval))
  if (!identical(result$fixed_external, expected))
    stop("Fixed external propensity or treatment coding differs from the frozen study inputs")
  invisible(TRUE)
}

orthogonalisation_pooled_metrics <- function(chains, data) {
  if (!length(chains)) stop("Cannot pool an empty orthogonalisation chain collection")
  arms <- vapply(chains, function(chain) as.character(chain$orthogonalisation$arm), character(1L))
  if (length(unique(arms)) != 1L || !arms[1L] %in% c("Standard", "Orthogonal"))
    stop("Cannot pool different orthogonalisation arms")
  for (chain in chains) {
    if (!is.null(chain$propensity_draws) || !is.null(chain$theta_samples))
      stop("Internal propensity draws are incompatible with fixed-external pooling")
    identity <- chain$orthogonalisation
    if (!identical(identity$propensity_train, as.numeric(data$pi_train)) ||
        !identical(identity$propensity_eval, as.numeric(data$pi_eval)))
      stop("External propensity changed between orthogonalisation chains")
    if (arms[1L] == "Orthogonal") orthogonalisation_validate_fixed_result(chain, data)
    if (!identical(chain$propensity_metrics, data$propensity_metrics))
      stop("External propensity diagnostics differ between orthogonalisation chains")
  }
  # These are diagnostics for one fixed external point estimate, not posterior
  # draws. In particular, do not manufacture propensity credible intervals.
  data$propensity_metrics
}

orthogonalisation_preflight <- function(cfg, fit_function = NULL) {
  # An actual tiny Gaussian Gibbs fit verifies fixed-input semantics and native
  # capability. This is not convergence or frequentist coverage validation.
  if (is.null(fit_function)) fit_function <- getExportedValue("bslim", "bcf_linear_probit")
  generated <- generate_data_medical(n = 48L, is_te_hetero = TRUE,
    is_mu_nonlinear = TRUE, seed = 81033L, RCT = FALSE, scenario = "default",
    z_diff = 0.5, BCF = FALSE, sigma_sq = 1)
  X <- as.matrix(generated[, c("x1", "x2", "x3", "x4", "x5_1", "x5_2")])
  z <- as.numeric(generated$z > 0)
  # This deterministic toy external score needs no arm dependency. Production
  # preparation fits main-effect and pairwise-interaction bayesglm once,
  # outside every chain.
  pi <- as.numeric(0.1 + 0.8 * stats::plogis(X[, 1L] - 0.3 * X[, 2L]))
  test <- seq_len(12L)
  gp <- list(
    standardize = TRUE, standardize_cov = FALSE, sample_sigma2_global = TRUE,
    sigma2_global_init = 1, sigma2_global_shape = 1, sigma2_global_scale = 0.001,
    adaptive_coding = FALSE, control_coding_init = -0.5, treated_coding_init = 0.5,
    sample_global_prior = "hc-hs", unlink = TRUE, gibbs = TRUE, hn_scale = 1,
    p0 = 0, use_ncp = FALSE, regularize_ATE = FALSE, robust = FALSE,
    simple_prior = FALSE, probit_outcome_model = FALSE, sigma_residual = 0,
    interaction_rule = "continuous_or_binary", propensity_covariate = "none",
    propensity_seperate = "none", orthogonalize = TRUE, save_output = FALSE,
    random_seed = 81034L, num_chains = 1L, keep_gfr = FALSE,
    keep_burnin = FALSE, keep_every = 1L, verbose = FALSE)
  fit <- tryCatch(fit_function(
    X_train = X, y_train = generated$y, Z_train = z, propensity_train = pi,
    X_test = X[test, , drop = FALSE], Z_test = z[test], propensity_test = pi[test],
    num_gfr = 1L, num_burnin = 2L, num_mcmc = 3L,
    general_params = gp, prognostic_forest_params = list(num_trees = 5L)),
    error = function(error) error)
  if (inherits(fit, "error")) stop(
    "Orthogonalisation preflight failed for the installed spbart package: ",
    conditionMessage(fit), "\n",
    "This study requires the fixed-external implementation in spbart 0.1.1.9004 or later ",
    "(commit 85d2a845), orthogonalize=TRUE, gibbs=TRUE, and sampled residual variance. ",
    "Install that implementation before submitting the orthogonalisation array. ",
    "The main and IHDP suites can be preflighted independently.", call. = FALSE)
  params <- fit$model_params
  if (!identical(params$orthogonalization, "fixed-external") ||
      !isTRUE(params$orthogonalize) || !identical(params$propensity_covariate, "none") ||
      !identical(params$propensity_seperate, "none") ||
      !identical(params$adaptive_coding, FALSE) ||
      !identical(params$internal_propensity_model, FALSE) ||
      !isTRUE(params$sample_sigma2_global) || !identical(params$probit_outcome_model, FALSE) ||
      !identical(params$variance_update, "hchs-full-conditional-rng-20260911") ||
      !identical(as.integer(params$num_prognostic_covariates), as.integer(ncol(X))))
    stop("Orthogonalisation preflight: incompatible fixed-external Gaussian variance metadata", call. = FALSE)
  if (!is.null(fit$theta_samples) || !is.null(fit$bart_propensity_model) ||
      !identical(fit$propensity_train, pi) || !identical(fit$propensity_test, pi[test]) ||
      !identical(fit$treatment_coding_train, z - pi) ||
      !identical(fit$treatment_coding_test, z[test] - pi[test]))
    stop("Orthogonalisation preflight: supplied external propensity was not preserved exactly", call. = FALSE)
  coefficients <- fit$orthogonal_coefficients
  if (!is.matrix(coefficients) || nrow(coefficients) != 3L ||
      ncol(coefficients) < 1L + ncol(X) || any(!is.finite(coefficients)))
    stop("Orthogonalisation preflight: missing finite retained coefficient draws", call. = FALSE)
  for (field in c("sigma2_global_samples", "Tau_glob"))
    if (length(fit[[field]]) != 3L || any(!is.finite(fit[[field]])) || any(fit[[field]] <= 0))
      stop("Orthogonalisation preflight: missing positive variance/global-scale draws", call. = FALSE)
  for (part in c("train", "test")) {
    nr <- if (part == "train") nrow(X) else length(test)
    for (field in c("mu_hat_", "tau_hat_", "y_hat_")) {
      draws <- fit[[paste0(field, part)]]
      if (!is.matrix(draws) || !identical(dim(draws), c(nr, 3L)) || any(!is.finite(draws)))
        stop("Orthogonalisation preflight: invalid retained prediction draws", call. = FALSE)
    }
    expected <- fit[[paste0("mu_hat_", part)]] + fit[[paste0("tau_hat_", part)]] *
      fit[[paste0("treatment_coding_", part)]]
    if (!isTRUE(all.equal(fit[[paste0("y_hat_", part)]], expected, tolerance = 1e-12)))
      stop("Orthogonalisation preflight: outcome predictions do not use Z - external propensity", call. = FALSE)
  }
  invisible(list(ok = TRUE, mode = "hc-hs Gibbs with fixed external propensity",
                 retained_draws = 3L, propensity_in_mu = FALSE))
}
