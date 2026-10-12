# Shared CANVAS rerun protocol. Each chain starts a separate package fit.
# Package num_chains is deliberately one: no sampler state passes between calls.
canvas_mcmc_settings <- function() list(chains = 4L, gfr = 25L, burnin = 2000L,
                                      draws = 4000L, rhat_limit = 1.01, ess_min = 400)

canvas_run_root <- function(workspace) {
  Sys.getenv("BCF_CANVAS_RUN", unset = file.path(workspace, "spbart-8.3.5-rerun-mu-propensity-20261010"))
}

canvas_output_directory <- function(workspace, component, resume = FALSE) {
  path <- file.path(canvas_run_root(workspace), component)
  if (dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE)) && !resume)
    stop("Rerun output directory is not empty: ", path,
         ". Set BCF_CANVAS_RUN to a new directory; existing results are not overwritten.")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop("Cannot create CANVAS output directory: ", path)
  path
}

canvas_chain_seeds <- function(seed, model, chains = 4L) {
  model_index <- match(model, c("BCF", "SemiParam"))
  # Ten slots per input seed keep both model/chain sequences disjoint.
  values <- as.double(seed) * 10 + (model_index - 1L) * 4 + seq_len(chains)
  if (length(seed) != 1L || !is.finite(seed) || seed != floor(seed) || seed < 1 ||
      is.na(model_index) || chains != 4L || any(values > .Machine$integer.max))
    stop("Invalid CANVAS seed/model/chain configuration")
  as.integer(values)
}

canvas_gp <- function(model, seed) {
  common <- list(cutpoint_grid_size = 100, standardize = TRUE, propensity_covariate = "mu",
    num_chains = 1L, probit_outcome_model = FALSE, random_seed = seed,
    keep_burnin = FALSE, keep_gfr = FALSE, keep_every = 1L, verbose = FALSE)
  if (model == "BCF") return(modifyList(common, list(adaptive_coding = TRUE)))
  if (model != "SemiParam") stop("Unknown CANVAS model")
  modifyList(common, list(sample_sigma2_global = TRUE, sigma2_global_init = 1,
    sigma2_global_shape = 1, sigma2_global_scale = 0.001,
    unlink = TRUE, gibbs = TRUE, hn_scale = 1, adaptive_coding = FALSE,
    control_coding_init = -0.5, treated_coding_init = 0.5,
    orthogonalize = FALSE, simple_prior = FALSE, regularize_ATE = FALSE,
    use_ncp = FALSE, p0 = 0, robust = FALSE, sigma_residual = 0,
    interaction_rule = "continuous_or_binary", step_out = 0.5, max_steps = 150,
    sample_global_prior = "hc-hs",
    propensity_seperate = "none", standardize_cov = FALSE))
}

canvas_draw_matrix <- function(value, observations, draws, label) {
  if (!is.matrix(value) || !identical(dim(value), c(as.integer(observations), as.integer(draws))) ||
      any(!is.finite(value))) stop("Invalid retained ", label, " draws")
  value
}

canvas_coefficient_matrix <- function(value, draws, label) {
  dims <- dim(value)
  if (length(dims) != 3L || dims[1] != 1L || dims[2] != draws || any(!is.finite(value)))
    stop("Invalid single-chain ", label, " array")
  matrix(value, nrow = draws, ncol = dims[3])
}

canvas_parameter_group <- function(variable) {
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
canvas_diagnostic_groups <- function(diagnostics,cfg) {
  if(!nrow(diagnostics)) return(data.frame())
  required <- c("variable","sd","rhat","ess_bulk","ess_tail","mcse_mean","pass")
  if(!all(required %in% names(diagnostics))) stop("Incomplete per-quantity diagnostics")
  # Reconstruct labels when collecting older saved objects. Never change their
  # per-quantity flags or the legacy all-quantities gate.
  diagnostics$parameter_group <- canvas_parameter_group(diagnostics$variable)
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

canvas_diagnostics <- function(chains, settings) {
  if (!requireNamespace("posterior", quietly = TRUE)) stop("Install the posterior R package")
  variables <- colnames(chains[[1]])
  if (length(chains) != settings$chains || any(vapply(chains, nrow, integer(1)) != settings$draws) ||
      !all(vapply(chains, function(x) identical(colnames(x), variables), logical(1))))
    stop("Diagnostic chains have inconsistent dimensions or variables")
  # Work in blocks rather than making one extra array of every patient CATE.
  rows <- lapply(seq.int(1L, length(variables), by = 100L), function(start) {
    index <- seq.int(start, min(start + 99L, length(variables)))
    values <- array(NA_real_, c(settings$draws, length(chains), length(index)),
                    dimnames = list(NULL, NULL, variables[index]))
    for (chain in seq_along(chains)) values[, chain, ] <- chains[[chain]][, index, drop = FALSE]
    result <- as.data.frame(posterior::summarise_draws(posterior::as_draws_array(values),
      "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean"))
    result$parameter_group <- canvas_parameter_group(result$variable)
    result$mcse_mean_over_sd <- ifelse(is.finite(result$sd) & result$sd > 0, result$mcse_mean / result$sd, NA_real_)
    result$chains <- length(chains)
    result$draws_per_chain <- settings$draws
    result$pass <- is.finite(result$rhat) & result$rhat <= settings$rhat_limit &
      is.finite(result$ess_bulk) & result$ess_bulk >= settings$ess_min &
      is.finite(result$ess_tail) & result$ess_tail >= settings$ess_min
    result
  })
  do.call(rbind, rows)
}

canvas_fit_chains <- function(model, Xtr, Ytr, Ztr, Xte, Zte, propensity, seed,
                              settings = canvas_mcmc_settings(), fit_function = NULL,
                              predict_function = NULL, forest_params = NULL) {
  if (!requireNamespace("posterior", quietly = TRUE)) stop("Install the posterior R package")
  if (length(propensity) != 1L || !is.finite(propensity) || propensity <= 0 || propensity >= 1)
    stop("CANVAS propensity must be one training-pool treatment proportion between zero and one")
  seeds <- canvas_chain_seeds(seed, model, settings$chains)
  if (is.null(fit_function)) fit_function <- if (model == "BCF") spbart::bcf else spbart::bcf_linear_probit
  if (is.null(predict_function)) predict_function <- if (model == "BCF") stats::predict else spbart::predict_linear_bcf
  diagnostic_chains <- vector("list", settings$chains)
  coefficients <- vector("list", settings$chains)
  keep_scale_draws <- identical(tolower(Sys.getenv("BCF_KEEP_SCALE_DRAWS", "false")), "true")
  ate <- vector("list", settings$chains)
  tau_mean <- y_mean <- numeric(nrow(Xte))
  reference <- NULL
  for (chain in seq_len(settings$chains)) {
    set.seed(seeds[chain])
    args <- list(X_train = Xtr, y_train = Ytr,
      Z_train = if (model == "BCF") Ztr else Ztr - 0.5,
      num_gfr = settings$gfr, num_burnin = settings$burnin, num_mcmc = settings$draws,
      general_params = canvas_gp(model, seeds[chain]))
    # The same training-pool fraction enters only the prognostic forest in both models.
    args$propensity_train <- rep(propensity, length(Ztr))
    if (!is.null(forest_params)) args$prognostic_forest_params <- forest_params
    fit <- do.call(fit_function, args)
    prediction_args <- list(object = fit, X = Xte,
      Z = as.matrix(if (model == "BCF") Zte else Zte - 0.5))
    prediction_args$propensity <- rep(propensity, nrow(Xte))
    prediction <- do.call(predict_function, prediction_args)
    tau <- canvas_draw_matrix(prediction$tau_hat, nrow(Xte), settings$draws, "CATE")
    yhat <- canvas_draw_matrix(prediction$y_hat, nrow(Xte), settings$draws, "outcome")
    sigma <- as.numeric(fit$sigma2_global_samples)
    if (length(sigma) != settings$draws || any(!is.finite(sigma) | sigma <= 0))
      stop("Missing positive residual-variance draws")
    ate[[chain]] <- colMeans(tau)
    tau_mean <- tau_mean + rowMeans(tau) / settings$chains
    y_mean <- y_mean + rowMeans(yhat) / settings$chains
    diagnostic <- cbind(ATE = ate[[chain]], sigma2 = sigma)
    if (model == "SemiParam") {
      beta <- canvas_coefficient_matrix(fit$Beta, settings$draws, "main coefficient")
      bint <- if (is.null(fit$Beta_int)) matrix(numeric(), settings$draws, 0L) else
        canvas_coefficient_matrix(fit$Beta_int, settings$draws, "interaction coefficient")
      alpha <- as.numeric(fit$alpha)
      if (length(alpha) != settings$draws || any(!is.finite(alpha))) stop("Invalid intercept draws")
      metadata <- list(model_params = fit$model_params, interaction_pairs = fit$interaction_pairs)
      if (is.null(reference)) reference <- metadata
      # The seed is intentionally different; encoding and outcome scales must agree.
      if (!identical(metadata$interaction_pairs, reference$interaction_pairs) ||
          !identical(metadata$model_params$outcome_scale, reference$model_params$outcome_scale) ||
          !identical(metadata$model_params$outcome_mean, reference$model_params$outcome_mean))
        stop("Independent chains used different outcome scales or interaction columns")
      coefficients[[chain]] <- list(alpha = alpha, beta = beta, beta_int = bint, sigma = sigma)
      colnames(beta) <- paste0("beta_", seq_len(ncol(beta)))
      if (ncol(bint)) colnames(bint) <- paste0("interaction_", seq_len(ncol(bint)))
      diagnostic <- cbind(diagnostic, alpha = alpha, beta, bint)
      # Tau_int is fixed at one in this unlinked hc-hs implementation; it is not a sampled scale.
      for (field in c("Tau", "Tau_glob")) {
        value <- fit[[field]]
        if (!is.null(value)) {
          value <- if (length(dim(value)) == 3L) canvas_coefficient_matrix(value, settings$draws, field) else
            matrix(as.numeric(value), nrow = settings$draws)
          if (any(!is.finite(value))) stop("Non-finite shrinkage draws")
          colnames(value) <- if (field == "Tau_glob") "lambda_global" else paste0("lambda_local_", seq_len(ncol(value)))
          if (keep_scale_draws) coefficients[[chain]]$scales[[field]] <- value
          diagnostic <- cbind(diagnostic, value)
        }
      }
    }
    cate <- t(tau); colnames(cate) <- paste0("CATE_", seq_len(nrow(Xte)))
    diagnostic_chains[[chain]] <- cbind(diagnostic, cate)
    rm(fit, prediction, tau, yhat, diagnostic, cate); gc(verbose = FALSE)
  }
  diagnostics <- canvas_diagnostics(diagnostic_chains, settings)
  diagnostics$model <- model
  diagnostic_groups <- canvas_diagnostic_groups(diagnostics, settings)
  diagnostic_groups$model <- model
  pooled <- NULL
  if (model == "SemiParam") {
    stack <- function(field) {
      p <- ncol(coefficients[[1]][[field]])
      result <- array(NA_real_, c(settings$chains, settings$draws, p))
      for (chain in seq_len(settings$chains)) result[chain, , ] <- coefficients[[chain]][[field]]
      result
    }
    pooled <- list(alpha = do.call(rbind, lapply(coefficients, `[[`, "alpha")),
      Beta = stack("beta"), Beta_int = stack("beta_int"),
      sigma2_global_samples = do.call(rbind, lapply(coefficients, `[[`, "sigma")),
      interaction_pairs = reference$interaction_pairs, model_params = reference$model_params)
    pooled$model_params$num_chains <- settings$chains
    pooled$model_params$num_samples <- settings$chains * settings$draws
    if (keep_scale_draws) pooled$scale_draws <- lapply(names(coefficients[[1]]$scales), function(field) {
      values <- array(NA_real_, c(settings$chains, settings$draws, ncol(coefficients[[1]]$scales[[field]])))
      for (chain in seq_len(settings$chains)) values[chain, , ] <- coefficients[[chain]]$scales[[field]]
      values
    })
    if (keep_scale_draws) names(pooled$scale_draws) <- names(coefficients[[1]]$scales)
  }
  list(tau_mean = tau_mean, y_mean = y_mean, ate_draws = unlist(ate, use.names = FALSE),
    posterior = pooled, diagnostics = diagnostics, diagnostic_groups = diagnostic_groups,
    mcmc = list(protocol = "canvas-independent-4chain-mu-propensity-8.3.5", settings = settings,
                model = model, seeds = seeds, independent_package_calls = TRUE,
                propensity = list(value = propensity, covariate = "mu", treatment_design = "none"),
                scale_draws_retained = model == "SemiParam" && keep_scale_draws))
}
