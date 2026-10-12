# IHDP analysis: the five manuscript methods, ten outer folds and full-sample fits.
# Outer-fold outcomes never enter preprocessing, nuisances or learner training.
# Source this file; it has no side effects and never installs/downloads packages.

ihdp_value <- function(cfg, name, default) {
  value <- cfg[[name]]
  if (is.null(value)) default else value
}

ihdp_jobs <- function(cfg) {
  jobs <- expand.grid(
    replicate = as.integer(ihdp_value(cfg, "ihdp_replicates", 1:100)),
    fold = 0:as.integer(ihdp_value(cfg, "ihdp_folds", 10L)),
    model = c("SemiParam", "BCF", "DR", "R", "S"),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  jobs$suite <- "ihdp"
  jobs$n <- 747L
  jobs$scenario <- "ihdp"
  jobs$heterogeneity <- TRUE
  jobs$mu_form <- "ihdp"
  jobs[, c("suite", "n", "replicate", "scenario", "heterogeneity", "mu_form", "model", "fold")]
}

ihdp_load <- function(replicate, directory) {
  if (length(replicate) != 1L || is.na(replicate) || replicate < 1L || replicate > 100L)
    stop("IHDP replicate must be in 1..100")
  path <- file.path(directory, sprintf("ihdp_npci_%d.csv", replicate))
  if (!file.exists(path)) stop("Missing IHDP input: ", path,
    ". Prepare the 100 public inputs before submitting compute jobs.")
  manifest_path <- file.path(directory, "input_manifest.csv")
  if (!file.exists(manifest_path)) stop("Missing IHDP input manifest: ", manifest_path)
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  if (nrow(manifest) != 100L || anyDuplicated(manifest$replicate) ||
      !setequal(manifest$replicate, 1:100)) stop("Incomplete IHDP input manifest")
  entry <- manifest[manifest$replicate == replicate, , drop = FALSE]
  if (!requireNamespace("digest", quietly = TRUE)) stop("IHDP input verification requires digest")
  source_sha256 <- digest::digest(file = path, algo = "sha256")
  if (entry$file != basename(path) || entry$sha256 != source_sha256)
    stop("IHDP input does not match its SHA-256 manifest: ", path)
  raw <- as.matrix(utils::read.csv(path, header = FALSE, check.names = FALSE))
  storage.mode(raw) <- "double"
  if (!identical(dim(raw), c(747L, 30L)) || any(!is.finite(raw)))
    stop("IHDP input must be a finite 747 x 30 numeric matrix: ", path)
  if (any(!raw[, 1] %in% c(0, 1)) || sum(raw[, 1]) != 139)
    stop("Unexpected IHDP treatment assignment: ", path)
  X <- raw[, 6:30, drop = FALSE]
  colnames(X) <- paste0("x", 1:25)
  list(X = X, y = raw[, 2], z = raw[, 1], truth = raw[, 5] - raw[, 4],
       source_file = normalizePath(path, winslash = "/", mustWork = TRUE),
       source_md5 = unname(tools::md5sum(path)), source_sha256 = source_sha256)
}

ihdp_fold_ids <- function(n, replicate, folds = 10L) {
  if (folds < 2L || folds >= n) stop("Invalid IHDP fold count")
  # Retain the historical outer assignment rule and its dependence on input order.
  set.seed(1000L + as.integer(replicate))
  sample(rep(seq_len(folds), length.out = n))
}

ihdp_design <- function(X_train, X_eval) {
  X_train <- as.matrix(X_train); X_eval <- as.matrix(X_eval)
  if (ncol(X_train) != 25L || ncol(X_eval) != 25L ||
      any(!is.finite(X_train)) || any(!is.finite(X_eval)))
    stop("IHDP design requires 25 finite raw covariates")
  # The benchmark's original covariate order: x4 is four-level categorical;
  # x1,x2,x3,x5,x6 are continuous; x7..x25 are binary, with x14 coded 1/2.
  continuous <- c(1L, 2L, 3L, 5L, 6L)
  centres <- colMeans(X_train[, continuous, drop = FALSE])
  scales <- apply(X_train[, continuous, drop = FALSE], 2, stats::sd)
  if (any(!is.finite(scales) | scales <= 0)) stop("Constant continuous IHDP training covariate")
  category_levels <- sort(unique(X_train[, 4]))
  if (length(category_levels) != 4L || any(!X_eval[, 4] %in% category_levels))
    stop("All four x4 levels must be represented in the IHDP training partition")
  contrasts <- stats::contr.sum(4L)
  transform <- function(X) {
    columns <- list()
    for (j in seq_len(25L)) {
      name <- paste0("x", j)
      if (j %in% continuous) {
        k <- match(j, continuous)
        columns[[name]] <- (X[, j] - centres[k]) / scales[k]
      } else if (j == 4L) {
        encoded <- contrasts[match(X[, j], category_levels), , drop = FALSE]
        for (k in 1:3) columns[[paste0(name, k)]] <- encoded[, k]
      } else {
        expected <- if (j == 14L) c(1, 2) else c(0, 1)
        if (any(!X[, j] %in% expected)) stop("Unexpected values for binary ", name)
        columns[[name]] <- ifelse(X[, j] == expected[1], -1, 1)
      }
    }
    matrix <- do.call(cbind, columns)
    storage.mode(matrix) <- "double"
    matrix
  }
  train <- transform(X_train); evaluation <- transform(X_eval)
  list(train = train, evaluation = evaluation,
       metadata = list(continuous = paste0("x", continuous), centres = centres,
         scales = scales, category_levels = category_levels,
         columns = colnames(train), coding = "sum contrasts; binary -1/+1",
         moments_from = "outer training partition only",
         non_continuous_index = which(!colnames(train) %in% paste0("x", continuous)) - 1L))
}

ihdp_inner_folds <- function(z, folds, seed) {
  folds <- min(as.integer(folds), min(table(z)))
  if (folds < 2L) stop("Too few observations per treatment arm for nuisance cross-fitting")
  set.seed(seed)
  ids <- integer(length(z))
  # Stratification avoids empty treatment arms in an inner training partition.
  for (arm in c(0, 1)) {
    index <- which(z == arm)
    ids[index] <- sample(rep(seq_len(folds), length.out = length(index)))
  }
  ids
}

ihdp_propensity_fit <- function(X, z, X_new, trim) {
  fit <- pairwise_logistic_propensity(X, z, X_new, trim = trim)
  structure(fit$probability, propensity_model = fit$metadata)
}

ihdp_sl_library <- function() {
  packages <- c("SuperLearner", "glmnet", "gam", "xgboost", "ranger")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing IHDP learner packages: ", paste(missing, collapse = ", "))
  environment <- new.env(parent = asNamespace("SuperLearner"))
  gam <- SuperLearner::create.Learner("SL.gam", tune = list(deg.gam = c(2, 3, 4)),
                                    name_prefix = "IHDP_GAM", env = environment)
  list(names = c("SL.glm", "SL.glmnet", gam$names, "SL.xgboost", "SL.ranger"),
       environment = environment)
}

ihdp_sl_fit <- function(X, y, X_new, seed, library = NULL) {
  if (is.null(library)) library <- ihdp_sl_library()
  set.seed(seed)
  fit <- SuperLearner::SuperLearner(Y = y, X = as.data.frame(X),
    newX = as.data.frame(X_new), SL.library = library$names,
    family = stats::gaussian(), env = library$environment)
  prediction <- as.numeric(fit$SL.predict)
  if (length(prediction) != nrow(X_new) || any(!is.finite(prediction)))
    stop("IHDP SuperLearner returned incomplete predictions")
  list(prediction = prediction, weights = fit$coef, cv_risk = fit$cvRisk,
       library = library$names)
}

ihdp_nuisance <- function(X_train, y_train, z_train, X_eval, seed,
                          folds = 5L, trim = 0.01, outcomes = FALSE,
                          propensity_fit = ihdp_propensity_fit, outcome_fit = ihdp_sl_fit,
                          outcome_stages = c("m", "mu0", "mu1")) {
  if (!is.finite(trim) || trim <= 0 || trim >= 0.5) stop("Invalid propensity trimming")
  inner <- ihdp_inner_folds(z_train, folds, seed)
  pi_train <- rep(NA_real_, length(y_train))
  m_train <- mu0_train <- mu1_train <- pi_train
  weights <- list()
  propensity_models <- list()
  if (any(!outcome_stages %in% c("m", "mu0", "mu1"))) stop("Unknown IHDP nuisance outcome stage")
  for (fold in sort(unique(inner))) {
    training <- which(inner != fold); validation <- which(inner == fold)
    propensity <- propensity_fit(X_train[training, , drop = FALSE],
      z_train[training], X_train[validation, , drop = FALSE], trim)
    pi_train[validation] <- as.numeric(propensity)
    propensity_models[[paste0("inner_", fold)]] <- attr(propensity, "propensity_model")
    if (outcomes) {
      sets <- list(m = training, mu0 = training[z_train[training] == 0],
                   mu1 = training[z_train[training] == 1])
      for (stage in outcome_stages) {
        index <- sets[[stage]]
        value <- outcome_fit(X_train[index, , drop = FALSE], y_train[index],
          X_train[validation, , drop = FALSE], seed + 100L * fold + match(stage, names(sets)))
        if (stage == "m") m_train[validation] <- value$prediction
        if (stage == "mu0") mu0_train[validation] <- value$prediction
        if (stage == "mu1") mu1_train[validation] <- value$prediction
        weights[[paste0(stage, "_", fold)]] <- value[c("weights", "cv_risk", "library")]
      }
    }
  }
  propensity <- propensity_fit(X_train, z_train, X_eval, trim)
  propensity_models$outer_training <- attr(propensity, "propensity_model")
  pi_eval <- as.numeric(propensity)
  if (any(!is.finite(c(pi_train, pi_eval)))) stop("Incomplete IHDP propensity predictions")
  psi <- NULL
  if (outcomes) {
    values <- list(m = m_train, mu0 = mu0_train, mu1 = mu1_train)
    if (any(!is.finite(unlist(values[outcome_stages])))) stop("Incomplete IHDP nuisances")
    if (all(c("mu0", "mu1") %in% outcome_stages)) {
      mu_observed <- ifelse(z_train == 1, mu1_train, mu0_train)
      psi <- mu1_train - mu0_train + (z_train - pi_train) *
        (y_train - mu_observed) / (pi_train * (1 - pi_train))
    }
  }
  list(pi_train = pi_train, pi_eval = pi_eval, m_train = m_train,
       mu0_train = mu0_train, mu1_train = mu1_train, psi_train = psi,
       inner_folds = inner, trim = trim, weights = weights, propensity_models = propensity_models,
       description = "Nuisances cross-fitted strictly within the outer training partition")
}

ihdp_prepare <- function(job, cfg) {
  replicate <- as.integer(job$replicate)
  fold <- as.integer(job$fold)
  folds <- as.integer(ihdp_value(cfg, "ihdp_folds", 10L))
  if (length(fold) != 1L || is.na(fold) || !fold %in% 0:folds) stop("Invalid IHDP outer fold")
  raw <- ihdp_load(replicate, cfg$ihdp_data_dir)
  fold_ids <- ihdp_fold_ids(length(raw$y), replicate, folds)
  training <- if (fold == 0L) seq_along(raw$y) else which(fold_ids != fold)
  evaluation <- if (fold == 0L) seq_along(raw$y) else which(fold_ids == fold)
  design <- ihdp_design(raw$X[training, , drop = FALSE], raw$X[evaluation, , drop = FALSE])
  nuisance_seed <- as.integer(ihdp_value(cfg, "seed_base", 20260921L) + 1000L * replicate + fold)
  nuisance <- ihdp_nuisance(design$train, raw$y[training], raw$z[training], design$evaluation,
    seed = nuisance_seed, folds = ihdp_value(cfg, "ihdp_inner_folds", 5L),
    trim = ihdp_value(cfg, "ihdp_trim", 0.01), outcomes = job$model %in% c("DR", "R"),
    outcome_stages = if (job$model == "DR") c("mu0", "mu1") else "m")
  identity <- list(replicate = replicate, fold = fold, source_sha256 = raw$source_sha256,
    indices = training, nuisance_seed = nuisance_seed)
  list(X_train = design$train, y_train = raw$y[training], z_train = raw$z[training],
       X_eval = design$evaluation, pi_train = nuisance$pi_train, pi_eval = nuisance$pi_eval,
       truth_eval = raw$truth[evaluation], y_eval = raw$y[evaluation], z_eval = raw$z[evaluation],
       eval_index = evaluation, train_index = training, fold_ids = fold_ids,
       training_identity = identity,
       input = raw[c("source_file", "source_md5", "source_sha256")], preprocessing = design$metadata,
       preprocess = design$metadata,
       metadata = list(input = raw[c("source_file", "source_md5", "source_sha256")], training_identity = identity,
         fold_ids = fold_ids, inner_folds = nuisance$inner_folds, trim = nuisance$trim,
         nuisance_weights = nuisance$weights, propensity_models = nuisance$propensity_models,
         protocol = nuisance$description),
       nuisance = nuisance, sample = if (fold == 0L) "IS" else "OOS")
}

ihdp_fit_chain <- function(job, data, chain, cfg) {
  model <- as.character(job$model)
  seed <- fit_seed(job, chain)
  budget <- cfg
  if (isTRUE(ihdp_value(cfg, "ihdp_historical_budget", FALSE)))
    stop("This release uses the common MCMC budget for fitting and diagnostics; change that budget only in a fresh run")
  if (model %in% c("SemiParam", "BCF")) {
    fit <- if (model == "SemiParam") fit_sp_chain else fit_bcf_chain
    return(fit(X_train = data$X_train, y_train = data$y_train, z_binary = data$z_train,
      X_eval = data$X_eval, pi_train = data$pi_train, pi_eval = data$pi_eval,
      seed = seed, cfg = budget, z_eval = data$z_eval))
  }
  if (chain != 1L) stop("Point-estimation learners have one fit, not MCMC chains")
  if (model == "DR") {
    value <- ihdp_sl_fit(data$X_train, data$nuisance$psi_train, data$X_eval, seed)
    tau <- value$prediction
    y_hat <- rep(NA_real_, length(tau))
    metadata <- value[c("weights", "cv_risk", "library")]
  } else if (model == "R") {
    if (!requireNamespace("grf", quietly = TRUE)) stop("IHDP R-learner requires grf")
    set.seed(seed)
    fit <- grf::causal_forest(data$X_train, data$y_train, data$z_train,
      Y.hat = data$nuisance$m_train, W.hat = data$pi_train, seed = seed, num.threads = 1L)
    tau <- as.numeric(stats::predict(fit, data$X_eval, num.threads = 1L)$predictions)
    y_hat <- rep(NA_real_, length(tau))
    metadata <- list(method = "grf::causal_forest", num.trees = fit$num.trees)
  } else if (model == "S") {
    train <- cbind(data$X_train, A = data$z_train)
    eval1 <- cbind(data$X_eval, A = 1)
    eval0 <- cbind(data$X_eval, A = 0)
    value <- ihdp_sl_fit(train, data$y_train, rbind(eval1, eval0), seed)
    n <- nrow(data$X_eval)
    mu1 <- value$prediction[seq_len(n)]; mu0 <- value$prediction[n + seq_len(n)]
    tau <- mu1 - mu0
    y_hat <- ifelse(data$z_eval == 1, mu1, mu0)
    metadata <- value[c("weights", "cv_risk", "library")]
  } else stop("Unknown IHDP model: ", model)
  if (length(tau) != nrow(data$X_eval) || any(!is.finite(tau)))
    stop("Incomplete IHDP learner predictions")
  list(tau = matrix(tau, ncol = 1L), y_hat = y_hat,
       diagnostic_draws = NULL, learner = metadata)
}

ihdp_aggregate <- function(results, cfg) {
  if (!length(results)) stop("No IHDP results supplied")
  jobs <- do.call(rbind, lapply(results, function(x) x$job))
  keys <- paste(jobs$replicate, jobs$fold, jobs$model, sep = ":")
  if (anyDuplicated(keys)) stop("Duplicate IHDP fit identities in aggregation")
  expected <- ihdp_jobs(cfg)
  expected_keys <- paste(expected$replicate, expected$fold, expected$model, sep = ":")
  if (!setequal(keys, expected_keys)) stop("IHDP aggregation requires every requested replicate, fold and model")
  rows <- list(); individuals <- list()
  models <- c("SemiParam", "BCF", "DR", "R", "S")
  for (replicate in as.integer(ihdp_value(cfg, "ihdp_replicates", 1:100))) {
    for (model in models) for (sample in c("IS", "OOS")) {
      selected <- which(jobs$replicate == replicate & jobs$model == model &
        if (sample == "IS") jobs$fold == 0L else jobs$fold > 0L)
      selected <- selected[order(jobs$fold[selected])]
      combined <- do.call(rbind, lapply(selected, function(k) {
        value <- results[[k]]$individual
        value$fold <- jobs$fold[k]
        value
      }))
      if (nrow(combined) != 747L || anyDuplicated(combined$index) ||
          !setequal(combined$index, seq_len(747L)))
        stop("Incomplete or overlapping IHDP evaluation rows: ", replicate, "/", model, "/", sample)
      combined <- combined[order(combined$index), , drop = FALSE]
      if (any(!is.finite(c(combined$tau_mean, combined$truth)))) stop("Incomplete IHDP CATE estimates")
      error <- combined$tau_mean - combined$truth
      bayesian <- model %in% c("SemiParam", "BCF")
      if (bayesian && any(!is.finite(c(combined$tau_lower, combined$tau_upper))))
        stop("Missing IHDP Bayesian intervals")
      if (bayesian && (!"p_positive" %in% names(combined) ||
          length(combined$p_positive) != 747L ||
          any(!is.finite(combined$p_positive)) ||
          any(combined$p_positive < 0 | combined$p_positive > 1)))
        stop("Missing or invalid retained IHDP posterior sign probabilities")
      lower <- if (bayesian) combined$tau_lower else rep(NA_real_, 747L)
      upper <- if (bayesian) combined$tau_upper else rep(NA_real_, 747L)
      if (bayesian && any(lower > upper)) stop("Reversed IHDP Bayesian intervals")
      exclude <- lower > 0 | upper < 0
      # All evaluated children enter the AUC, not only intervals excluding zero.
      # Average ranks give tied positive/negative scores half a concordance.
      positive <- combined$truth > 0
      n_positive <- sum(positive); n_nonpositive <- sum(!positive)
      auc_sign <- if (bayesian && n_positive > 0L && n_nonpositive > 0L) {
        ranks <- rank(combined$p_positive, ties.method = "average")
        (sum(ranks[positive]) - n_positive * (n_positive + 1) / 2) /
          (n_positive * n_nonpositive)
      } else NA_real_
      slope <- if (stats::sd(combined$tau_mean) > 0)
        unname(stats::coef(stats::lm(combined$truth ~ combined$tau_mean))[2]) else NA_real_
      row <- data.frame(replicate = replicate, model = model, sample = sample, n = 747L,
        cate_rmse = sqrt(mean(error^2)), cate_mae = mean(abs(error)),
        cate_spearman = if (stats::sd(combined$tau_mean) > 0 && stats::sd(combined$truth) > 0)
          stats::cor(combined$tau_mean, combined$truth, method = "spearman") else NA_real_,
        calibration_slope = slope, ate_error = mean(error),
        sd_tau_true = stats::sd(combined$truth),
        outcome_rmse = if (all(is.finite(combined$y_hat))) sqrt(mean((combined$y_hat - combined$y)^2)) else NA_real_,
        cate_coverage = if (bayesian) mean(lower <= combined$truth & upper >= combined$truth) else NA_real_,
        cate_width = if (bayesian) mean(upper - lower) else NA_real_,
        exclude_zero = if (bayesian) mean(exclude) else NA_real_,
        sign_error_given_exclusion = if (bayesian && any(exclude))
          mean((lower > 0 & combined$truth <= 0 | upper < 0 & combined$truth >= 0)[exclude]) else NA_real_,
        sign_correct = if (bayesian && any(exclude))
          mean((lower > 0 & combined$truth > 0 | upper < 0 & combined$truth < 0)[exclude]) else NA_real_,
        auc_sign = auc_sign,
        stringsAsFactors = FALSE)
      # Point ATE error is meaningful across held-out rows. A joint OOS ATE
      # credible interval is not: different folds correspond to different posteriors.
      row$ate_lower <- row$ate_upper <- NA_real_
      row$convergence_failures <- sum(vapply(selected, function(k)
        identical(results[[k]]$diagnostic_pass, FALSE), logical(1)))
      for (level in c(80L, 90L, 50L)) {
        lo <- combined[[paste0("lo", level)]]; hi <- combined[[paste0("hi", level)]]
        available <- bayesian && length(lo) == 747L && length(hi) == 747L &&
          all(is.finite(c(lo, hi)))
        row[[paste0("cate_coverage", level)]] <- if (available)
          mean(lo <= combined$truth & hi >= combined$truth) else NA_real_
        row[[paste0("cate_width", level)]] <- if (available) mean(hi - lo) else NA_real_
      }
      if (sample == "IS" && bayesian) {
        metrics <- results[[selected]]$metrics
        if (all(c("ate_lower", "ate_upper") %in% names(metrics))) {
          row$ate_lower <- unname(metrics[["ate_lower"]])
          row$ate_upper <- unname(metrics[["ate_upper"]])
        }
      }
      rows[[length(rows) + 1L]] <- row
      combined$replicate <- replicate; combined$model <- model; combined$sample <- sample
      individuals[[length(individuals) + 1L]] <- combined
    }
  }
  list(metrics = do.call(rbind, rows), individual = do.call(rbind, individuals),
       note = "OOS posterior intervals are individual; no joint interval is formed across outer folds.")
}
