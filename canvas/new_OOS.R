# ==============================================================================
# CANVAS-R (train) -> CANVAS (external validation), sample-size sweep.
# CANVAS external-validation rerun for methodology-8.3.5.tex.
# Revision 2026-09-21: explicit Gibbs hc-hs settings, validation-size-normalized
# AUUC/net AUUC, and run-identified caches. Historical hchs results are NOT
# Gibbs results and are never resumed or overwritten by this output tag.
#
# Revision 2026-08-23. Changes relative to the previous version (kept as
# new_OOS.R.bak), each marked FIX n at its call site:
#   FIX 1  Prior. sample_global_prior was "half-cauchy", which in
#          stochtree/src/horseshoe_samplers.cpp applies the global scale to the
#          main effects as well (tau_glob_main = tau_glob). The paper's
#          hierarchy (eq. hshc_main has no lambda_glob) and every simulation /
#          IHDP arm use "hc-hs" (tau_glob_main = 1). Now "hc-hs".
#   FIX 2  n_draws was set to 1 inside the loop and never restored; N = 5000 of
#          a 5018 pool is not a Monte Carlo cell. Dropped; one fit on the full
#          CANVAS-R pool is run separately and labelled "full".
#   FIX 3  Subsampling was unseeded. Every draw is seeded from
#          seed_base + 1000*n + draw (subsample and both MCMC chains).
#   FIX 4  model.matrix() drops NA rows, so train/test indices would have been
#          applied to a shorter matrix. Complete cases are taken first and
#          nrow is asserted after the design is built.
#   FIX 5  Continuous covariates were scaled with train+test pooled moments.
#          Stacking is kept (it aligns factor levels); moments come from the
#          training pool only.
#   FIX 6  Dashboard panels 2-4 were drawn on the draws pooled into one frame
#          (n_draws copies of the same 3004 patients): uplift y-axis inflated
#          by n_draws, global ranking mixed draws with different shrinkage,
#          rolling-window SE divided by 500 rows holding ~25 distinct patients,
#          title correlation was the pooled (upward-biased) value. Every curve
#          is now computed WITHIN a draw and averaged across draws.
#   FIX 7  Both models receive the same constant treatment proportion from
#          the full training pool at fit and prediction, only in the prognostic
#          forest. It is not added to either treatment design.
#   FIX 8  Raw AUUC includes the ATE * n / 2 baseline as well as ranking gain.
#          Net AUUC subtracts this random-targeting baseline; raw is kept.
#   FIX 9  "008080" is not a colour (missing #).
#   FIX 10 Posterior-mean CATE predictions are saved per draw (.rds) and the
#          metric table is flushed after every draw. Slim files do NOT retain
#          ATE posterior draws or all Y predictions; complete matching metrics
#          are therefore required when resuming a saved prediction file.
#   FIX 11 (2026-08-24) The YODA workspace quota is ~25 MB and the first run
#          filled it (2 MB per .rds, 1.6 MB per EPS, 0.5 MB per TIFF), so the
#          N = 4000 and full-pool figures were never written. Export is now PNG
#          only on the VM (TIFF/EPS for submission are rebuilt off-VM from the
#          .rds), the .rds are stored as 1e-4-resolution integers under xz
#          (~0.3 MB each). Resume now requires matching run identities in the
#          manifest, prediction file and complete metric rows for that size.
#          Delete the stray Balci_*.eps / Balci_*.tiff on the VM by hand.
# ==============================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(spbart)

# A source checkout update is not an installed-package update. Stop before
# touching participant inputs if this R session lacks the Gaussian repair.
check_hchs_gibbs_package <- function(helper = NULL) {
  if (is.null(helper)) {
    ns <- asNamespace("spbart")
    if (!exists(".hchs_variance_parameters", ns, inherits = FALSE))
      stop("Install the repaired Gaussian hc-hs spbart build in this R session's library before running CANVAS.")
    helper <- get(".hchs_variance_parameters", ns)
  }
  value <- helper(c(1, 2), c(.5), c(.8), c(2, 4), 3, 1, .001)
  expected <- .001 + (5 + (.5 / 2)^2 + (.8 / 12)^2) / 2
  if (!isTRUE(all.equal(value$shape, 3, tolerance = 1e-12)) ||
      !isTRUE(all.equal(value$scale, expected, tolerance = 1e-12)))
    stop("Installed hc-hs residual-variance helper fails its algebraic check")
  invisible(TRUE)
}
check_hchs_gibbs_package()

# ------------------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------------------
# Shared sampler protocol; keep canvas_mcmc.R beside this script.
canvas_source_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(canvas_source_path)) {
  canvas_cli_path <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(canvas_cli_path)) canvas_source_path <- sub("^--file=", "", canvas_cli_path[1])
}
CANVAS_CODE_DIR <- if (is.null(canvas_source_path)) getwd() else dirname(normalizePath(canvas_source_path, mustWork = TRUE))
source(file.path(CANVAS_CODE_DIR, "canvas_mcmc.R"))
if (!requireNamespace("posterior", quietly = TRUE)) stop("Install the posterior R package before the CANVAS rerun")

WORKSPACE_DIR <- Sys.getenv("BCF_WORKSPACE", unset = "X:/Workspace")
CFG <- list(
  train_file = file.path(WORKSPACE_DIR, "analsis.RData"),    # CANVAS-R, n = 5018
  test_file  = file.path(WORKSPACE_DIR, "analsis_2.RData"),  # CANVAS,   n = 3004
  out_dir    = file.path(canvas_run_root(WORKSPACE_DIR), "out-of-sample"),
  tag        = "spbart_8.3.5_four_chains_mu_propensity_20261010",
  n_draws    = 20,
  n_targets  = c(500, 1000, 1500, 2000, 3000, 4000),
  fit_full   = TRUE,          # one fit on the whole CANVAS-R pool (the r = 0.791 sentence)
  seed_base  = 20260823,
  num_chains = 4L, num_gfr = 25L, num_burnin = 2000L, num_mcmc = 4000L,
  benefit_sign = -1,          # Y = change in HbA1c: lower is better
  resume = identical(tolower(Sys.getenv("BCF_CANVAS_RESUME", "false")), "true"),
  # --- journal artwork (see FIGURES section) ---
  fig_numbers = list(),       # e.g. list(N4000 = "3", compound = "S2"); unnumbered -> N label in the file name
  # FIX 11: the YODA workspace has a ~25 MB quota. On the VM export "png" only
  # and keep the slim .rds; the submission TIFF/EPS are rebuilt OFF the VM from
  # those .rds by re-running the FIGURES + COMPOUND sections (or just
  # make_dashboard / make_compound_figure) with export = c("tiff", "eps").
  export      = c("png"),
  dpi         = 800,          # TIFF: graphs/drawings, 800 dpi preferred, 600 minimum
  png_dpi     = 300,          # PNG draft for the tex; small enough for the quota
  fig_width   = 7.5, fig_height = 6.2,   # inches at print size (two-column width)
  dashboard_height = 3.6,     # one two-panel row; avoids wasting vertical space
  base_pt     = 8,            # lettering at print size; identical in every panel
  dashboard_panels = c("cate", "uplift"),          # two panes per training-size scenario
  compound_sizes  = c(1000, 2000, 3000, 4000),   # rows of the one compound figure
  compound_panels = c("cate", "uplift")          # columns of the compound figure
)
CFG$out_dir <- canvas_output_directory(WORKSPACE_DIR, "out-of-sample", CFG$resume)
out_path <- function(...) file.path(CFG$out_dir, paste0(...))
MODEL_COLS <- c("BCF" = "tomato", "Semi-parametric BART" = "#008080")   # FIX 9

# ------------------------------------------------------------------------------
# DESIGN MATRIX (script-side copy; the package's internal copy is not used
# because standardize_cov = FALSE below)
# ------------------------------------------------------------------------------
standardize_X_by_index <- function(X_initial,
                                   process_data = TRUE,
                                   interaction_rule = c("continuous", "continuous_or_binary", "all"),
                                   cat_coding_method = c("sum", "difference"),
                                   scale_idx = NULL) {   # FIX 5: rows whose moments are used for scaling
  interaction_rule <- match.arg(interaction_rule)
  cat_coding_method <- match.arg(cat_coding_method)
  if (cat_coding_method == "difference" && !requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for difference coding (contr.sdif).", call. = FALSE)
  }
  if (!is.data.frame(X_initial) && !is.matrix(X_initial)) stop("X_initial must be a data.frame or a matrix.")
  if (anyNA(X_initial)) stop("X_initial contains NA; take complete cases first (model.matrix would drop rows).")  # FIX 4

  original_colnames <- colnames(X_initial)
  if (is.null(original_colnames)) original_colnames <- paste0("V", 1:ncol(X_initial))
  X_df <- as.data.frame(X_initial); colnames(X_df) <- original_colnames
  if (is.null(scale_idx)) scale_idx <- seq_len(nrow(X_df))

  if (process_data) {
    all_cols_info <- lapply(1:ncol(X_df), function(j_idx) {
      col_data <- X_df[[j_idx]]
      n_unique <- length(unique(col_data[!is.na(col_data)]))
      is_continuous <- FALSE; is_binary <- FALSE; is_categorical <- FALSE
      processed_col_data <- NULL; should_include <- FALSE
      if (n_unique > 1) {
        if (is.logical(col_data)) col_data <- as.character(col_data)
        if (n_unique == 2) {
          factor_col <- factor(col_data)
          processed_col_data <- ifelse(factor_col == levels(factor_col)[1], -1, 1)
          is_binary <- TRUE; should_include <- TRUE
        } else if (is.character(col_data) || is.factor(col_data) || (is.numeric(col_data) && n_unique > 2 && n_unique < 4)) {
          col_as_factor <- factor(col_data)
          contrasts(col_as_factor) <- if (cat_coding_method == "difference") MASS::contr.sdif(nlevels(col_as_factor)) else contr.sum(nlevels(col_as_factor))
          processed_col_data <- col_as_factor
          is_categorical <- TRUE; should_include <- TRUE
        } else if (is.numeric(col_data)) {
          processed_col_data <- col_data
          is_continuous <- TRUE; should_include <- TRUE
        }
      }
      list(processed_data = processed_col_data, is_continuous = is_continuous,
           is_binary = is_binary, is_categorical = is_categorical,
           original_idx = j_idx, col_name = original_colnames[j_idx],
           include_in_intermediate = should_include)
    })
    valid_cols_info <- all_cols_info[sapply(all_cols_info, `[[`, "include_in_intermediate")]
    if (length(valid_cols_info) == 0) stop("No valid columns to process.")
    X_intermediate_list <- lapply(valid_cols_info, `[[`, "processed_data")
    names(X_intermediate_list) <- sapply(valid_cols_info, `[[`, "col_name")
    X_intermediate_df <- as.data.frame(X_intermediate_list, stringsAsFactors = FALSE)
    intermediate_map <- data.frame(
      col_name_intermediate = sapply(valid_cols_info, `[[`, "col_name"),
      original_idx = sapply(valid_cols_info, `[[`, "original_idx"),
      is_continuous = sapply(valid_cols_info, `[[`, "is_continuous"),
      is_binary = sapply(valid_cols_info, `[[`, "is_binary"),
      is_categorical = sapply(valid_cols_info, `[[`, "is_categorical"),
      stringsAsFactors = FALSE)

    # FIX 5: centre/scale with the moments of scale_idx (the training pool) only
    cols_to_scale <- intermediate_map$col_name_intermediate[intermediate_map$is_continuous]
    for (cn in cols_to_scale) {
      m <- mean(X_intermediate_df[[cn]][scale_idx]); s <- sd(X_intermediate_df[[cn]][scale_idx])
      X_intermediate_df[[cn]] <- (X_intermediate_df[[cn]] - m) / s
    }

    X_final_with_int <- model.matrix(~ ., data = X_intermediate_df)
    X_final <- X_final_with_int[, -1, drop = FALSE]
    attr_assign <- attr(X_final_with_int, "assign")[-1]
    matched_rows <- match(colnames(X_intermediate_df)[attr_assign], intermediate_map$col_name_intermediate)
    X_final_var_info <- data.frame(
      original_idx = intermediate_map$original_idx[matched_rows],
      is_continuous = intermediate_map$is_continuous[matched_rows],
      is_binary = intermediate_map$is_binary[matched_rows],
      is_categorical = intermediate_map$is_categorical[matched_rows],
      col_name_final = colnames(X_final), stringsAsFactors = FALSE)
    non_continous_idx_cpp <- which(!sapply(all_cols_info, `[[`, "is_continuous")) - 1
  } else {
    X_final <- as.matrix(X_df)
    all_cols_info <- lapply(1:ncol(X_df), function(j_idx) {
      col_data <- X_df[[j_idx]]
      n_unique <- length(unique(col_data[!is.na(col_data)]))
      is_continuous <- FALSE; is_binary <- FALSE; is_categorical <- FALSE
      if (n_unique > 1) {
        if (n_unique == 2) is_binary <- TRUE
        else if (is.character(col_data) || is.factor(col_data) || (is.numeric(col_data) && n_unique < min(20, 0.2 * nrow(X_df)))) is_categorical <- TRUE
        else if (is.numeric(col_data)) is_continuous <- TRUE
      }
      list(is_continuous = is_continuous, is_binary = is_binary, is_categorical = is_categorical)
    })
    X_final_var_info <- data.frame(
      original_idx = 1:ncol(X_final),
      is_continuous = sapply(all_cols_info, `[[`, "is_continuous"),
      is_binary = sapply(all_cols_info, `[[`, "is_binary"),
      is_categorical = sapply(all_cols_info, `[[`, "is_categorical"),
      col_name_final = original_colnames, stringsAsFactors = FALSE)
    non_continous_idx_cpp <- which(!X_final_var_info$is_continuous) - 1
  }
  stopifnot(nrow(X_final) == nrow(X_initial))   # FIX 4

  p_int_calculated <- 0
  n_final_cols <- ncol(X_final)
  if (n_final_cols > 1) {
    is_candidate <- switch(interaction_rule,
                           "all" = rep(TRUE, n_final_cols),
                           "continuous_or_binary" = X_final_var_info$is_continuous | X_final_var_info$is_binary,
                           "continuous" = X_final_var_info$is_continuous)
    for (i in 1:(n_final_cols - 1)) for (j in (i + 1):n_final_cols)
      if (is_candidate[i] || is_candidate[j]) p_int_calculated <- p_int_calculated + 1
  }
  list(X_final = X_final, p_int = p_int_calculated,
       non_continous_idx_cpp = non_continous_idx_cpp, X_final_var_info = X_final_var_info)
}

# ------------------------------------------------------------------------------
# METRIC HELPERS. All operate on ONE draw (n = n_test rows). Never pool draws.
# ------------------------------------------------------------------------------
# Rolling-window calibration: exactly as the supplement describes (window 500,
# >= 6 per arm, SE of the within-window difference in means recorded).
get_rolling_calibration <- function(df, model_col, model_name, window_size = 500, step = 20) {
  df_sorted <- df %>% arrange(.data[[model_col]])
  n <- nrow(df_sorted)
  if (n < window_size) { window_size <- max(10, floor(n / 4)); step <- max(1, floor(window_size / 5)) }
  indices <- seq(1, max(1, n - window_size + 1), by = step)
  res <- lapply(indices, function(i) {
    w <- df_sorted[i:(i + window_size - 1), ]
    n1 <- sum(w$Z == 1); n0 <- sum(w$Z == 0)
    if (n1 > 5 && n0 > 5) {
      data.frame(Mean_Pred = mean(w[[model_col]]),
                 Obs_Diff  = mean(w$Y[w$Z == 1]) - mean(w$Y[w$Z == 0]),
                 SE        = sqrt(var(w$Y[w$Z == 1]) / n1 + var(w$Y[w$Z == 0]) / n0))
    } else NULL
  })
  res_df <- bind_rows(res); res_df$Model <- model_name; res_df
}

# Non-overlapping quintile bins: used for the tabulated calibration MAE
# (rolling windows overlap, so an MAE from them has uninterpretable uncertainty).
get_calibration <- function(df, model_col, model_name) {
  df %>% mutate(Bin = ntile(.data[[model_col]], 5)) %>% group_by(Bin) %>%
    summarise(Mean_Pred = mean(.data[[model_col]]),
              Obs_Diff = mean(Y[Z == 1]) - mean(Y[Z == 0]),
              SE = sqrt(var(Y[Z == 1]) / (sum(Z == 1) + 1e-6) + var(Y[Z == 0]) / (sum(Z == 0) + 1e-6)),
              Model = model_name, .groups = "drop")
}

# Uplift curve. Patients are sorted by predicted benefit (benefit_sign * CATE,
# largest first), then cumulative benefit in the two randomised arms is compared.
# Equal scores are evaluated as one block.  Their order is otherwise arbitrary;
# treating them as a sequence would manufacture ranking performance when a model
# has collapsed to a homogeneous CATE.  The explicit (0, 0) origin makes the
# trapezoidal AUUC and the random-targeting chord use the same definition.
# Keep the historical count-scaled `uplift` for explicit legacy comparisons.
# The reported curve is uplift_normalized = uplift / n_validation = fraction
# targeted * within-prefix mean benefit difference, in outcome units.
get_uplift <- function(df, model_col, model_name, benefit_sign = CFG$benefit_sign) {
  if (!all(c("Y", "Z", model_col) %in% names(df))) stop("uplift data are missing Y, Z, or the score column")
  if (!all(df$Z %in% c(0, 1)) || !all(is.finite(df$Y)) || !all(is.finite(df[[model_col]])))
    stop("uplift requires finite outcomes/scores and Z coded 0/1")
  if (length(unique(df$Z)) != 2) stop("uplift requires both treatment arms")
  if (length(benefit_sign) != 1L || !benefit_sign %in% c(-1, 1)) stop("benefit_sign must be -1 or +1")

  blocks <- df %>%
    transmute(.b = benefit_sign * Y, .s = benefit_sign * .data[[model_col]], Z = Z) %>%
    arrange(desc(.s)) %>%
    mutate(.tie = cumsum(c(TRUE, diff(.s) != 0))) %>%
    group_by(.tie) %>%
    summarise(n = n(), n_1 = sum(Z), n_0 = sum(1 - Z),
              b_1 = sum(.b * Z), b_0 = sum(.b * (1 - Z)), .groups = "drop") %>%
    mutate(k = cumsum(n), cum_n_1 = cumsum(n_1), cum_n_0 = cumsum(n_0),
           cum_b_1 = cumsum(b_1), cum_b_0 = cumsum(b_0),
           uplift = ifelse(cum_n_1 > 0 & cum_n_0 > 0,
                           (cum_b_1 / cum_n_1 - cum_b_0 / cum_n_0) * k, 0),
           frac = k / sum(n), uplift_normalized = uplift / sum(n),
           n_validation = sum(n), Model = model_name) %>%
    select(frac, uplift, uplift_normalized, n_validation, Model)

  bind_rows(data.frame(frac = 0, uplift = 0, uplift_normalized = 0,
                       n_validation = nrow(df), Model = model_name), blocks)
}
get_auuc <- function(u, scale = c("normalized", "count")) {
  scale <- match.arg(scale)
  u <- u[order(u$frac), , drop = FALSE]
  y <- u[[if (scale == "normalized") "uplift_normalized" else "uplift"]]
  stopifnot(nrow(u) >= 2, length(y) == nrow(u), all(is.finite(y)),
            u$frac[1] == 0, tail(u$frac, 1) == 1, all(diff(u$frac) > 0))
  sum(diff(u$frac) * (head(y, -1) + tail(y, -1)) / 2)
}
# FIX 8: random-targeting baseline is the chord from (0,0) to (1, uplift(1)),
# whose area is endpoint/2 on the selected scale. A homogeneous score has zero
# net AUUC. Normalization divides by n_validation, NOT by the ATE or max gain.
get_net_auuc <- function(u, scale = c("normalized", "count")) {
  scale <- match.arg(scale)
  y <- u[[if (scale == "normalized") "uplift_normalized" else "uplift"]]
  get_auuc(u, scale) - y[which.max(u$frac)] / 2
}

# Small deterministic metric audit.  Run this on the VM with
# BCF_RUN_METRIC_TESTS=true before a costly fit; it proves that tied scores have
# zero ranking gain and that a correctly ordered two-risk example gains AUUC.
audit_uplift_metrics <- function() {
  d <- data.frame(Y = c(4, 0, 0, 0), Z = c(1, 0, 1, 0), score = c(2, 2, 0, 0))
  tied <- get_uplift(transform(d, score = 0), "score", "audit", benefit_sign = 1)
  ranked <- get_uplift(d, "score", "audit", benefit_sign = 1)
  stopifnot(isTRUE(all.equal(get_net_auuc(tied), 0, tolerance = 1e-12)),
            get_net_auuc(ranked) > 0,
            isTRUE(all.equal(get_auuc(tied), tail(tied$uplift_normalized, 1) / 2, tolerance = 1e-12)))
  invisible(TRUE)
}
if (identical(tolower(Sys.getenv("BCF_RUN_METRIC_TESTS")), "true")) audit_uplift_metrics()

# ------------------------------------------------------------------------------
# PREDICT for the semi-parametric model. Same logic as the exported
# stochtree::predict_linear_bcf (kept inline so the script does not depend on
# the package version installed on the VM). Coefficients are fitted on the
# standardised residual (y - ybar)/sd, hence the * outcome_scale. Interaction
# columns follow object$interaction_pairs, which is built by the same double
# loop as generate_interaction_matrix() at training.
# ------------------------------------------------------------------------------
predict_linear_bcf_patched <- function(object, X, Z, propensity = NULL, ...) {
  if ((!is.data.frame(X)) && (!is.matrix(X))) stop("X must be a matrix or dataframe")
  X_processed <- preprocessPredictionData(X, object$train_set_metadata)
  X_linear <- X_processed
  if (!is.null(Z) && is.null(dim(Z))) Z <- as.matrix(as.numeric(Z))
  if (!is.null(propensity) && is.null(dim(propensity))) propensity <- as.matrix(propensity)

  if ((object$model_params$propensity_covariate != "none") || (object$model_params$propensity_seperate != "none")) {
    if (is.null(propensity)) {
      if (object$model_params$internal_propensity_model) {
        propensity <- as.matrix(rowMeans(predict(object$bart_propensity_model, X_processed)$y_hat))
      } else stop("Propensity scores must be provided for this model.")
    }
  }
  X_forest <- X_processed
  if (object$model_params$propensity_covariate != "none") X_forest <- cbind(X_processed, propensity)
  forest_dataset_pred <- createForestDataset(X_forest, Z)
  prop_sep <- object$model_params$propensity_seperate
  if (!is.null(prop_sep) && prop_sep == "tau") X_linear <- cbind(X_processed, propensity)

  num_chains <- dim(object$Beta)[1]; num_mcmc <- dim(object$Beta)[2]
  total_samples <- num_chains * num_mcmc
  alpha_samples <- as.vector(t(object$alpha))
  beta_samples <- matrix(aperm(object$Beta, c(2, 1, 3)), nrow = total_samples)
  has_interactions <- !is.null(object$Beta_int) && (dim(object$Beta_int)[3] > 0)
  linear_pred <- matrix(rep(alpha_samples, each = nrow(X)), nrow = nrow(X)) + as.matrix(X_linear) %*% t(beta_samples)
  if (has_interactions) {
    beta_int_samples <- matrix(aperm(object$Beta_int, c(2, 1, 3)), nrow = total_samples)
    int_pairs <- object$interaction_pairs
    X_int <- matrix(0, nrow = nrow(X_linear), ncol = ncol(int_pairs))
    for (k in seq_len(ncol(int_pairs))) X_int[, k] <- X_linear[, int_pairs[1, k]] * X_linear[, int_pairs[2, k]]
    linear_pred <- linear_pred + (X_int %*% t(beta_int_samples))
  }
  y_std <- object$model_params$outcome_scale; y_bar <- object$model_params$outcome_mean
  mu_hat <- object$forests_mu$predict(forest_dataset_pred) * y_std + y_bar
  tau_hat <- if (object$model_params$standardize) linear_pred * y_std else linear_pred
  list(mu_hat = mu_hat, tau_hat = tau_hat, y_hat = mu_hat + tau_hat * as.vector(Z))
}

# ------------------------------------------------------------------------------
# MODEL SETTINGS. Explicit intended hc-hs Gibbs hierarchy; no adaptive coding
# in the semi-parametric arm, whose supplied Z is already centred at +/-0.5.
# ------------------------------------------------------------------------------
gp_semiparam <- function(seed) canvas_gp("SemiParam", seed)
gp_bcf <- function(seed) canvas_gp("BCF", seed)

# Fit each chain from a separate initialization and pool posterior predictions.
fit_both <- function(Xtr, Ytr, Ztr, Xte, Zte, pi_tr, seed) {
  settings <- canvas_mcmc_settings()
  settings$chains <- CFG$num_chains; settings$gfr <- CFG$num_gfr
  settings$burnin <- CFG$num_burnin; settings$draws <- CFG$num_mcmc
  standard <- canvas_fit_chains("BCF", Xtr, Ytr, Ztr, Xte, Zte, pi_tr, seed, settings)
  linear <- canvas_fit_chains("SemiParam", Xtr, Ytr, Ztr, Xte, Zte, pi_tr, seed, settings)
  list(Std_CATE = standard$tau_mean, Lin_CATE = linear$tau_mean,
    Std_Y = standard$y_mean, Lin_Y = linear$y_mean,
    Std_ATE_draws = standard$ate_draws, Lin_ATE_draws = linear$ate_draws,
    diagnostics = rbind(standard$diagnostics, linear$diagnostics),
    diagnostic_groups = rbind(standard$diagnostic_groups, linear$diagnostic_groups),
    chain_seeds = list(BCF = standard$mcmc$seeds, SemiParam = linear$mcmc$seeds))
}

save_fit_diagnostics <- function(pr, n, draw, seed) {
  result <- cbind(data.frame(Run_ID = RUN_ID, N_Train = n, Draw = draw, Seed = seed), pr$diagnostics)
  group_result <- cbind(data.frame(Run_ID = RUN_ID, N_Train = n, Draw = draw, Seed = seed), pr$diagnostic_groups)
  write.csv(group_result, out_path("mcmc_groups_", CFG$tag, "_N", n, "_draw", draw, ".csv"), row.names = FALSE)
  saveRDS(list(diagnostics = result, groups = group_result, chain_seeds = pr$chain_seeds),
    out_path("mcmc_", CFG$tag, "_N", n, "_draw", draw, ".rds"), compress = "xz")
  summary <- do.call(rbind, lapply(split(result, result$model), function(d) {
    data.frame(Run_ID = RUN_ID, N_Train = n, Draw = draw, Seed = seed,
      Model = d$model[1], Chains = CFG$num_chains, Draws_Per_Chain = CFG$num_mcmc,
      Parameters = nrow(d), Failing_Parameters = sum(!d$pass),
      Max_Rhat = if (any(is.finite(d$rhat))) max(d$rhat[is.finite(d$rhat)]) else NA_real_,
      Min_Bulk_ESS = if (any(is.finite(d$ess_bulk))) min(d$ess_bulk[is.finite(d$ess_bulk)]) else NA_real_,
      Min_Tail_ESS = if (any(is.finite(d$ess_tail))) min(d$ess_tail[is.finite(d$ess_tail)]) else NA_real_)
  }))
  write.csv(summary, out_path("mcmc_summary_", CFG$tag, "_N", n, "_draw", draw, ".csv"), row.names = FALSE)
  invisible(NULL)
}

# All per-draw metrics from one set of predictions.
draw_metrics <- function(pr, Y_test, Z_test) {
  tp <- data.frame(Y = Y_test, Z = Z_test, Std_CATE_OOS = pr$Std_CATE, Lin_CATE_OOS = pr$Lin_CATE)
  emp_ate <- mean(Y_test[Z_test == 1]) - mean(Y_test[Z_test == 0])
  emp_se  <- sqrt(var(Y_test[Z_test == 1]) / sum(Z_test == 1) + var(Y_test[Z_test == 0]) / sum(Z_test == 0))
  cov_fun <- function(d) {
    ci <- quantile(d, c(0.025, 0.975), na.rm = TRUE)
    c(raw = as.numeric(emp_ate >= ci[1] & emp_ate <= ci[2]),
      # empirical ATE is itself estimated: compare the gap against both uncertainties
      adj = as.numeric(abs(mean(d) - emp_ate) <= 1.96 * sqrt(sd(d)^2 + emp_se^2)))
  }
  cs <- get_calibration(tp, "Std_CATE_OOS", "BCF"); cl <- get_calibration(tp, "Lin_CATE_OOS", "Semi-parametric BART")
  us <- get_uplift(tp, "Std_CATE_OOS", "BCF");      ul <- get_uplift(tp, "Lin_CATE_OOS", "Semi-parametric BART")
  covs <- cov_fun(pr$Std_ATE_draws); covl <- cov_fun(pr$Lin_ATE_draws)
  data.frame(
    N_Validation = length(Y_test),
    RMSE_Y_Std = sqrt(mean((Y_test - pr$Std_Y)^2)), RMSE_Y_Lin = sqrt(mean((Y_test - pr$Lin_Y)^2)),
    CATE_Correlation = suppressWarnings(cor(pr$Std_CATE, pr$Lin_CATE)),
    CATE_Spearman    = suppressWarnings(cor(pr$Std_CATE, pr$Lin_CATE, method = "spearman")),
    SD_CATE_Std = sd(pr$Std_CATE), SD_CATE_Lin = sd(pr$Lin_CATE),
    Calib_MAE_Std = mean(abs(cs$Mean_Pred - cs$Obs_Diff)), Calib_MAE_Lin = mean(abs(cl$Mean_Pred - cl$Obs_Diff)),
    AUUC_Normalized_Std = get_auuc(us), AUUC_Normalized_Lin = get_auuc(ul),
    NetAUUC_Normalized_Std = get_net_auuc(us), NetAUUC_Normalized_Lin = get_net_auuc(ul),
    AUUC_CountScaled_Std = get_auuc(us, "count"), AUUC_CountScaled_Lin = get_auuc(ul, "count"),
    NetAUUC_CountScaled_Std = get_net_auuc(us, "count"), NetAUUC_CountScaled_Lin = get_net_auuc(ul, "count"),
    # Backward-compatible names retain their old count-scaled meaning.
    AUUC_Std = get_auuc(us, "count"), AUUC_Lin = get_auuc(ul, "count"),
    NetAUUC_Std = get_net_auuc(us, "count"), NetAUUC_Lin = get_net_auuc(ul, "count"),
    ATE_Bias_Std = mean(pr$Std_CATE) - emp_ate, ATE_Bias_Lin = mean(pr$Lin_CATE) - emp_ate,
    ATE_Cov_Std = covs[["raw"]], ATE_Cov_Lin = covl[["raw"]],
    ATE_CovAdj_Std = covs[["adj"]], ATE_CovAdj_Lin = covl[["adj"]]
  )
}

# ------------------------------------------------------------------------------
# FIGURES. Every curve is computed within a draw, then averaged (FIX 6), and the
# artwork follows the journal's illustration rules (FIG 1-6 below):
#   FIG 1  No tints: no alpha, no light fills. Bands are drawn as thin
#          dotted/dashed lines in the model's own colour; the CATE support
#          density is a line, not a filled ribbon; scatter points are solid,
#          small, and open/closed shapes distinguish models in greyscale.
#   FIG 2  Drawn at print size (CFG$fig_width x CFG$fig_height inches) with
#          CFG$base_pt lettering, the same in every panel and figure, so nothing
#          has to survive a reduction.
#   FIG 3  One key, inside the artwork (shared legend at the bottom), colour +
#          linetype + shape so the figure also reads in black and white.
#   FIG 4  Panels lettered a, b, c, d; compound figures are one file.
#   FIG 5  Exported as TIFF at CFG$dpi (800, LZW) and EPS (cairo; no
#          transparency, so EPS is exact), plus a small PNG for the tex draft.
#   FIG 6  File names carry the figure number where configured. The manuscript
#          supplies the caption; the artwork carries no editorial masthead.
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SLIM PREDICTION FILES (FIX 11). Full-precision doubles made each .rds ~2 MB,
# which alone filled most of the 25 MB YODA quota (the N = 4000 and full-pool
# figures then failed to write). Predictions are stored as integers at 1e-4
# resolution (ample on the HbA1c scale) with xz compression, and only draw 1 of
# the Y predictions is kept (the outcome panel shows one fit). ~0.3 MB per size.
# Legacy/unidentified files are deliberately rejected, never migrated in place.
# A manifest records the installed package build, sampler settings, input/design
# fingerprints and analysis helper code. The same ID is stored in EVERY CSV row
# and prediction file; resume also checks dimensions, draw numbers and seeds.
# ------------------------------------------------------------------------------
object_md5 <- function(x) {
  f <- tempfile("canvas-identity-", fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  saveRDS(x, f, version = 2, compress = FALSE)
  unname(tools::md5sum(f))
}

package_fingerprint <- function(package) {
  root <- find.package(package)
  files <- sort(list.files(root, recursive = TRUE, full.names = FALSE))
  # Includes installed R databases, DESCRIPTION and compiled DLL/shared library,
  # not just a version number shared by different local builds.
  files <- files[!file.info(file.path(root, files))$isdir]
  hashes <- tools::md5sum(file.path(root, files))
  if (!length(hashes) || anyNA(hashes)) stop("Cannot fingerprint installed package: ", package)
  names(hashes) <- files
  list(package = package, version = as.character(packageVersion(package)), files = hashes)
}

function_signature <- function(fun) list(formals = formals(fun), body = deparse(body(fun), width.cutoff = 500L))

make_run_identity <- function(config, package_builds, input_hashes, design_hash, analysis_code,
                              runtime = list(R = R.version.string, RNG = RNGkind(), platform = R.version$platform)) {
  # Artwork/export paths are deliberately excluded: identical fits can be
  # re-exported elsewhere or at another resolution without a new MCMC run.
  settings <- config[c("n_draws", "n_targets", "fit_full", "seed_base", "num_gfr",
                       "num_burnin", "num_mcmc", "num_chains", "benefit_sign")]
  specification <- list(schema = "canvas-independent-4chain-mu-propensity-8.3.5", settings = settings,
    semiparametric = gp_semiparam(0), bcf = gp_bcf(0),
    package_builds = package_builds, input_hashes = input_hashes,
    design_hash = design_hash, analysis_code = analysis_code, runtime = runtime)
  list(run_id = object_md5(specification), specification = specification)
}

ensure_run_manifest <- function(identity, file, cache_files = character()) {
  if (file.exists(file)) {
    old <- readRDS(file)
    if (!identical(old$run_id, identity$run_id) ||
        !identical(old$specification, identity$specification))
      stop("Run identity differs from the existing manifest. Use a NEW CFG$tag; do not mix or relabel historical results.")
  } else {
    if (any(file.exists(cache_files)))
      stop("Output files exist without a run manifest. Use a NEW CFG$tag; unidentified caches cannot be resumed.")
    saveRDS(identity, file, compress = "xz")
  }
  invisible(identity$run_id)
}

validate_metrics <- function(df, run_id, n_validation, required_columns, n_train = NULL,
                             expected_draws = NULL, expected_seeds = NULL) {
  required <- unique(c("Run_ID", "N_Train", "Draw", "Seed", "N_Validation", required_columns))
  if (!all(required %in% names(df))) stop("Metric cache has missing identity/schema columns; no legacy resume is allowed.")
  ids <- df[, c("N_Train", "Draw", "Seed", "N_Validation"), drop = FALSE]
  if (!all(vapply(ids, function(x) is.numeric(x) && all(is.finite(x)) && all(x == floor(x)), logical(1))) ||
      any(df$N_Train <= 0) || any(df$Draw <= 0) || any(df$N_Validation <= 0))
    stop("Metric cache has invalid numeric draw identifiers or sample sizes.")
  if (anyNA(df$Run_ID) || any(df$Run_ID != run_id) ||
      anyNA(df$N_Validation) || any(df$N_Validation != n_validation))
    stop("Metric cache contains a different/missing run identity or validation size.")
  if (anyNA(df[, c("N_Train", "Draw", "Seed")]) || anyDuplicated(df[, c("N_Train", "Draw")]))
    stop("Metric cache contains missing or duplicate draw identifiers.")
  if (!is.null(n_train)) {
    rows <- df[df$N_Train == n_train, , drop = FALSE]
    rows <- rows[order(rows$Draw), , drop = FALSE]
    if (!identical(as.integer(rows$Draw), as.integer(seq_len(expected_draws))) ||
        !identical(as.numeric(rows$Seed), as.numeric(expected_seeds)))
      stop("Prediction cache lacks complete matching metric draws/seeds; slim predictions cannot reconstruct all metrics. Use a new tag or explicitly refit this run.")
  }
  invisible(TRUE)
}

.enc <- function(m) matrix(as.integer(round(as.matrix(m) * 1e4)), nrow = nrow(as.matrix(m)))
slim_preds <- function(p, run_id, cache_key, seeds) {
  stopifnot(is.character(run_id), length(run_id) == 1L, nzchar(run_id),
            ncol(p$Std_CATE) == length(seeds), identical(dim(p$Std_CATE), dim(p$Lin_CATE)),
            identical(dim(p$Std_Y), dim(p$Std_CATE)), identical(dim(p$Lin_Y), dim(p$Std_CATE)),
            all(vapply(p[c("Std_CATE", "Lin_CATE", "Std_Y", "Lin_Y")],
                       function(m) is.matrix(m) && all(is.finite(m)) && all(abs(m * 1e4) < .Machine$integer.max), logical(1))))
  list(fmt = "i4-run-v1", run_id = run_id, cache_key = as.character(cache_key),
       seeds = as.numeric(seeds), n_validation = nrow(p$Std_CATE), n_draws = ncol(p$Std_CATE),
       Std_CATE = .enc(p$Std_CATE), Lin_CATE = .enc(p$Lin_CATE),
       Std_Y1 = .enc(cbind(p$Std_Y[, 1])), Lin_Y1 = .enc(cbind(p$Lin_Y[, 1])))
}
save_preds <- function(p, f, run_id, cache_key, seeds) {
  saveRDS(slim_preds(p, run_id, cache_key, seeds), f, compress = "xz")
  cat("  -> Saved", f, "\n")
}
load_preds <- function(f, run_id, cache_key, n_validation, seeds) {
  p <- readRDS(f)
  if (!identical(p$fmt, "i4-run-v1") || !identical(p$run_id, run_id) ||
      !identical(p$cache_key, as.character(cache_key)) ||
      !identical(p$seeds, as.numeric(seeds)))
    stop("Prediction cache identity/settings/seeds do not match this run; legacy caches cannot be resumed.")
  dims <- c(as.integer(n_validation), as.integer(length(seeds)))
  if (!identical(as.integer(p$n_validation), dims[1]) ||
      !identical(as.integer(p$n_draws), dims[2]) ||
      !identical(dim(p$Std_CATE), dims) || !identical(dim(p$Lin_CATE), dims) ||
      !identical(dim(p$Std_Y1), c(dims[1], 1L)) || !identical(dim(p$Lin_Y1), c(dims[1], 1L)) ||
      anyNA(p$Std_CATE) || anyNA(p$Lin_CATE) || anyNA(p$Std_Y1) || anyNA(p$Lin_Y1))
    stop("Prediction cache has incomplete or inconsistent dimensions/values.")
  list(Std_CATE = p$Std_CATE / 1e4, Lin_CATE = p$Lin_CATE / 1e4,
       Std_Y = p$Std_Y1 / 1e4, Lin_Y = p$Lin_Y1 / 1e4)
}

MODEL_LTY <- c("BCF" = "solid", "Semi-parametric BART" = "dashed")
MODEL_SHP <- c("BCF" = 16, "Semi-parametric BART" = 2)   # closed circle / open triangle

fig_theme <- function() {
  theme_classic(base_size = CFG$base_pt) +
    theme(plot.title = element_text(face = "bold", size = CFG$base_pt + 1),
          plot.subtitle = element_text(size = CFG$base_pt),
          plot.tag = element_text(face = "bold", size = CFG$base_pt + 2),
          axis.title = element_text(size = CFG$base_pt),
          axis.text = element_text(size = CFG$base_pt),
          legend.position = "none",
          axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
          plot.margin = margin(4, 6, 2, 4))
}
scale_models <- function(shape = FALSE, linetype = TRUE) c(
  list(scale_colour_manual(values = MODEL_COLS, name = NULL)),
  if (linetype) list(scale_linetype_manual(values = MODEL_LTY, name = NULL)),
  if (shape) list(scale_shape_manual(values = MODEL_SHP, name = NULL)))

# Honest comparison plots need the same numerical limits on both axes. Padding
# is computed once from the combined values, so neither model receives a more
# favourable visual scale. The same helper is also used across scenario rows.
shared_plot_limits <- function(x, padding = 0.04) {
  x <- x[is.finite(x)]
  if (!length(x)) stop("Cannot determine plot limits from non-finite values.")
  lim <- range(x)
  span <- diff(lim)
  if (span == 0) span <- max(abs(lim), 1) * 0.1
  as.numeric(lim + c(-1, 1) * span * padding)
}
shared_axis_breaks <- function(limits, n = 5) {
  b <- pretty(limits, n = n)
  b[b >= limits[1] & b <= limits[2]]
}

# Shared key placed once at the bottom (FIG 3): line type + marker + colour per
# model, drawn from a dummy plot so every panel reads from the same key.
make_key <- function() {
  d <- data.frame(Model = rep(names(MODEL_COLS), each = 2), x = rep(1:2, 2), y = 1)
  p <- ggplot(d, aes(x, y, colour = Model, linetype = Model, shape = Model)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.5, stroke = 0.4) +
    scale_models(shape = TRUE) +
    theme_void(base_size = CFG$base_pt) +
    theme(legend.position = "bottom", legend.key.width = unit(2.2, "lines"),
          legend.text = element_text(size = CFG$base_pt), legend.key.spacing.x = unit(1.5, "lines"))
  g <- ggplotGrob(p)
  i <- which(g$layout$name == "guide-box-bottom")
  if (length(i) == 0) i <- grep("guide-box", g$layout$name)[1]
  g$grobs[[i]]
}

# Write one figure object in every requested format.
save_fig <- function(grob, stem, fig_number = NULL,
                     width = CFG$fig_width, height = CFG$fig_height) {
  label <- if (is.null(fig_number)) stem else paste0("Fig", fig_number)
  base  <- out_path("canvas_", CFG$tag, "_", label)
  for (fmt in CFG$export) {
    f <- paste0(base, ".", fmt)
    if (fmt == "tiff") {
      tiff(f, width = width, height = height, units = "in", res = CFG$dpi, compression = "lzw")
    } else if (fmt == "eps") {
      cairo_ps(f, width = width, height = height, fallback_resolution = CFG$dpi)
    } else {
      png(f, width = width, height = height, units = "in", res = CFG$png_dpi)
    }
    grid.newpage(); grid.draw(grob); dev.off()
    cat("  -> Saved", f, "\n")
  }
}

# Assemble lettered panels with one shared key.  Figure number and caption belong
# in the manuscript, not in a line printed over the artwork.
assemble_figure <- function(panels, ncol, title, fig_number = NULL) {
  panels <- lapply(seq_along(panels), function(i) panels[[i]] + labs(tag = letters[i]))
  key <- make_key()
  body <- arrangeGrob(grobs = panels, ncol = ncol)
  arrangeGrob(body, key, ncol = 1, heights = unit.c(unit(1, "null"), unit(1.6, "lines")))
}

# The four panels for one training size. Returns ggplot objects (no legends).
build_panels <- function(preds, metrics, Y_test, Z_test, n_label, n_draws,
                         cate_limits = NULL) {
  models <- c(Std = "BCF", Lin = "Semi-parametric BART")

  # --- per-draw curves ---
  upl <- list(); cal <- list()
  for (d in seq_len(n_draws)) {
    tp <- data.frame(Y = Y_test, Z = Z_test, Std_CATE_OOS = preds$Std_CATE[, d], Lin_CATE_OOS = preds$Lin_CATE[, d])
    for (m in names(models)) {
      col <- paste0(m, "_CATE_OOS")
      upl[[length(upl) + 1]] <- get_uplift(tp, col, models[[m]]) %>% mutate(draw = d)
      cal[[length(cal) + 1]] <- get_rolling_calibration(tp, col, models[[m]]) %>% mutate(draw = d)
    }
  }
  upl <- bind_rows(upl); cal <- bind_rows(cal)

  # Tied-score blocks can differ across draws (particularly after 1e-4 cache
  # quantization). Interpolate the same piecewise-linear curves used for AUUC
  # onto their union of knots before averaging; every draw contributes at each
  # point. Dotted limits are min/max across draws, NOT confidence intervals.
  upl_grid <- sort(unique(upl$frac))
  upl_common <- bind_rows(lapply(split(upl, interaction(upl$Model, upl$draw)), function(u) {
    data.frame(Model = u$Model[1], frac = upl_grid,
               uplift = approx(u$frac, u$uplift_normalized, xout = upl_grid)$y)
  }))
  upl_avg <- upl_common %>% group_by(Model, frac) %>%
    summarise(lo = min(uplift), hi = max(uplift), uplift = mean(uplift), .groups = "drop")

  # Calibration: interpolate each draw onto a common grid over the INTERSECTION of
  # the draws' predicted-CATE ranges (never extrapolate); band = sqrt(mean
  # within-draw SE^2 + between-draw variance). Within-draw part is NOT divided by
  # n_draws: it is the same 3004 patients every draw.
  cal_avg <- bind_rows(lapply(models, function(mn) {
    cm <- cal %>% filter(Model == mn)
    if (nrow(cm) == 0) return(NULL)
    rng <- cm %>% group_by(draw) %>% summarise(lo = min(Mean_Pred), hi = max(Mean_Pred), k = n(), .groups = "drop") %>% filter(k >= 2)
    if (nrow(rng) == 0) return(NULL)
    g_lo <- max(rng$lo); g_hi <- min(rng$hi)
    if (!(g_hi > g_lo)) return(NULL)          # a draw collapsed to homogeneity
    grid <- seq(g_lo, g_hi, length.out = 200)
    per <- lapply(rng$draw, function(d) {
      cd <- cm %>% filter(draw == d) %>% arrange(Mean_Pred)
      data.frame(x = grid, y = approx(cd$Mean_Pred, cd$Obs_Diff, xout = grid, ties = mean)$y,
                 se = approx(cd$Mean_Pred, cd$SE, xout = grid, ties = mean)$y)
    }) %>% bind_rows()
    per %>% group_by(x) %>%
      summarise(Obs_Diff = mean(y), band = sqrt(mean(se^2) + ifelse(n() > 1, var(y), 0)), .groups = "drop") %>%
      mutate(Model = mn, n_draws_used = nrow(rng))
  }))

  draw_label <- if (n_draws > 1) sprintf("Training N = %s; %d draws", n_label, n_draws) else sprintf("Training N = %s; one fit", n_label)

  # --- panel a: Y prediction, ONE fit (draw 1) ---
  p1_df <- data.frame(Y = Y_test, BCF = preds$Std_Y[, 1], `Semi-parametric BART` = preds$Lin_Y[, 1], check.names = FALSE) %>%
    pivot_longer(-Y, names_to = "Model", values_to = "Pred")
  if (nrow(p1_df) > 3000) { set.seed(1); p1_df <- sample_n(p1_df, 3000) }
  outcome_limits <- shared_plot_limits(c(p1_df$Pred, p1_df$Y))
  outcome_breaks <- shared_axis_breaks(outcome_limits)
  p1 <- ggplot(p1_df, aes(x = Pred, y = Y, colour = Model, shape = Model)) +
    geom_point(size = 0.5, stroke = 0.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.4) +
    scale_x_continuous(limits = outcome_limits, breaks = outcome_breaks) +
    scale_y_continuous(limits = outcome_limits, breaks = outcome_breaks) +
    coord_equal() +
    labs(title = "Outcome prediction", subtitle = draw_label,
         x = "Predicted change in HbA1c\n(percentage points)",
         y = "Observed change in HbA1c\n(percentage points)") +
    fig_theme() + scale_models(shape = TRUE, linetype = FALSE)

  # --- panel b: per-patient CATE averaged over draws; title = mean WITHIN-draw r ---
  sc <- data.frame(Std = rowMeans(preds$Std_CATE), Lin = rowMeans(preds$Lin_CATE))
  pca <- prcomp(~ Std + Lin, data = sc)
  tls_slope <- pca$rotation["Lin", "PC1"] / pca$rotation["Std", "PC1"]
  tls_int <- pca$center["Lin"] - tls_slope * pca$center["Std"]
  # mean WITHIN-draw correlation (what tab:model_metrics reports), from the predictions themselves
  r_within <- mean(sapply(seq_len(ncol(preds$Std_CATE)), function(d)
    suppressWarnings(cor(preds$Std_CATE[, d], preds$Lin_CATE[, d]))), na.rm = TRUE)
  if (is.null(cate_limits)) cate_limits <- shared_plot_limits(c(sc$Std, sc$Lin))
  cate_breaks <- shared_axis_breaks(cate_limits)
  p2 <- ggplot(sc, aes(x = Std, y = Lin)) + geom_point(size = 0.4, colour = "black") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.4) +
    geom_abline(slope = tls_slope, intercept = tls_int, colour = "black", linewidth = 0.6) +
    scale_x_continuous(limits = cate_limits, breaks = cate_breaks) +
    scale_y_continuous(limits = cate_limits, breaks = cate_breaks) +
    coord_equal() +
    labs(title = sprintf("CATE agreement (Pearson r = %.2f)", r_within),
         subtitle = draw_label,
         x = "BCF-predicted CATE\n(HbA1c percentage points)",
         y = "Semi-parametric BART-predicted CATE\n(HbA1c percentage points)") + fig_theme()

  # --- panel c: rolling calibration, averaged; band as dotted lines; support density as a line ---
  if (!is.null(cal_avg) && nrow(cal_avg) > 0) {
    min_y <- min(cal_avg$Obs_Diff - 1.96 * cal_avg$band); max_y <- max(cal_avg$Obs_Diff + 1.96 * cal_avg$band)
    y_range <- max_y - min_y
    dens <- bind_rows(lapply(names(models), function(m) {
      dd <- density(as.vector(preds[[paste0(m, "_CATE")]]), na.rm = TRUE)
      data.frame(CATE = dd$x, Density_Scaled = (dd$y / max(dd$y)) * (y_range * 0.2) + min_y - 0.05 * y_range, Model = models[[m]])
    }))
    p3 <- ggplot() +
      geom_line(data = dens, aes(x = CATE, y = Density_Scaled, colour = Model, linetype = Model), linewidth = 0.35) +
      geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.4) +
      geom_line(data = cal_avg, aes(x = x, y = Obs_Diff - 1.96 * band, colour = Model), linetype = "dotted", linewidth = 0.35) +
      geom_line(data = cal_avg, aes(x = x, y = Obs_Diff + 1.96 * band, colour = Model), linetype = "dotted", linewidth = 0.35) +
      geom_line(data = cal_avg, aes(x = x, y = Obs_Diff, colour = Model, linetype = Model), linewidth = 0.7) +
      labs(title = "Calibration", subtitle = draw_label,
           x = "Predicted CATE\n(HbA1c percentage points)",
           y = "Observed treatment effect\n(HbA1c percentage points)") + fig_theme() + scale_models()
  } else {
    p3 <- ggplot() + annotate("text", 0, 0, size = CFG$base_pt / .pt,
                              label = "Calibration curve undefined:\na model collapsed to a homogeneous CATE") + theme_void()
  }

  # --- panel d: uplift, averaged; min-max band as dotted lines; random-targeting chord ---
  ref <- upl_avg %>% group_by(Model) %>% summarise(end = uplift[which.max(frac)], .groups = "drop")
  p4 <- ggplot(upl_avg, aes(x = frac, colour = Model)) +
    geom_line(aes(y = lo), linetype = "dotted", linewidth = 0.35) +
    geom_line(aes(y = hi), linetype = "dotted", linewidth = 0.35) +
    geom_line(aes(y = uplift, linetype = Model), linewidth = 0.7) +
    geom_segment(data = ref, aes(x = 0, y = 0, xend = 1, yend = end), linetype = "dotted", linewidth = 0.4,
                 colour = "black", inherit.aes = FALSE) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       labels = function(x) paste0(round(100 * x), "%"),
                       expand = expansion(mult = c(0, 0.01))) +
    labs(title = "Uplift curve", subtitle = paste0(draw_label, "; solid: mean, dotted: range"),
         x = "Proportion targeted, ranked by predicted benefit",
         y = "Cumulative uplift per validation participant\n(HbA1c percentage points)") +
    fig_theme() + scale_models()

  list(outcome = p1, cate = p2, calibration = p3, uplift = p4)
}

# One training size -> the same two-pane comparison used in the manuscript:
# CATE agreement and uplift. The outcome and calibration plots remain available
# from build_panels() for diagnostics but are not added to the scenario figure.
make_dashboard <- function(preds, metrics, Y_test, Z_test, n_label, n_draws) {
  panels <- build_panels(preds, metrics, Y_test, Z_test, n_label, n_draws)
  panels <- panels[CFG$dashboard_panels]
  fig <- assemble_figure(panels, ncol = length(panels), title = sprintf("Out-of-sample performance, N_train = %s", n_label),
                         fig_number = CFG$fig_numbers[[paste0("N", n_label)]])
  save_fig(fig, paste0("OOS_N", n_label), fig_number = CFG$fig_numbers[[paste0("N", n_label)]],
           height = CFG$dashboard_height)
  invisible(panels)
}

# Several training sizes -> ONE compound figure: one row per N, selected panels
# as columns (FIG 4). Rebuilt from the saved .rds files, so it can be rerun
# offline with different sizes/panels without refitting.
make_compound_figure <- function(sizes, panels = c("cate", "uplift"), metrics_all, Y_test, Z_test, fig_number = NULL) {
  scenarios <- list()
  for (n in sizes) {
    f <- out_path("canvas_preds_", CFG$tag, "_N", n, ".rds")
    if (!file.exists(f)) { cat("  compound: no predictions for N =", n, "- skipped\n"); next }
    seeds <- CFG$seed_base + 1000 * n + seq_len(CFG$n_draws)
    validate_metrics(metrics_all, RUN_ID, length(Y_test), METRIC_COLUMNS, n, CFG$n_draws, seeds)
    preds <- load_preds(f, RUN_ID, n, length(Y_test), seeds)
    scenarios[[length(scenarios) + 1]] <- list(n = n, preds = preds)
  }
  if (length(scenarios) == 0) return(invisible(NULL))

  # Use one CATE range for both model axes and all displayed training sizes.
  # This makes changes in shrinkage/magnitude comparable across scenario rows.
  cate_values <- unlist(lapply(scenarios, function(s)
    c(rowMeans(s$preds$Std_CATE), rowMeans(s$preds$Lin_CATE))), use.names = FALSE)
  cate_limits <- shared_plot_limits(cate_values)
  plot_sets <- lapply(scenarios, function(s)
    build_panels(s$preds, metrics_for(metrics_all, s$n), Y_test, Z_test,
                 as.character(s$n), ncol(s$preds$Std_CATE), cate_limits = cate_limits))

  # The right-hand panels also keep one y scale across scenarios. Both models
  # already share an axis within each uplift panel because their curves overlay.
  if ("uplift" %in% panels) {
    uplift_values <- unlist(lapply(plot_sets, function(pl)
      c(pl$uplift$data$lo, pl$uplift$data$hi)), use.names = FALSE)
    uplift_limits <- shared_plot_limits(c(0, uplift_values))
    uplift_breaks <- shared_axis_breaks(uplift_limits)
    plot_sets <- lapply(plot_sets, function(pl) {
      pl$uplift <- pl$uplift + scale_y_continuous(limits = uplift_limits, breaks = uplift_breaks)
      pl
    })
  }
  rows <- unlist(lapply(plot_sets, function(pl) pl[panels]), recursive = FALSE)

  h <- CFG$fig_height * length(scenarios) / 2
  fig <- assemble_figure(rows, ncol = length(panels), title = sprintf("Out-of-sample performance, N_train = %s", paste(sizes, collapse = ", ")),
                         fig_number = fig_number)
  save_fig(fig, paste0("OOS_compound_N", paste(sizes, collapse = "-")),
           fig_number = fig_number, height = h)
  invisible(fig)
}


# ==============================================================================
# DATA
# ==============================================================================
cat("\n======================================================\n")
cat("External Holdout Validation (CANVAS-R -> CANVAS)\n")
cat("======================================================\n")

env_train <- new.env(); load(CFG$train_file, envir = env_train); data_train <- env_train$analsis; data_train$Dataset <- "Train"
env_test  <- new.env(); load(CFG$test_file,  envir = env_test);  data_test  <- env_test$analsis;  data_test$Dataset  <- "Test"

X_col <- c("BLBMI","AGE","HHFFL","BLSBP","SEX","DIADUR","HCEREBFL","BLA1C","BLGFR","INSULBL","GLP1BL","BLACR","RACE","METFORBL","RANGIBL")
needed_cols <- c("y", "ARM", "Dataset", X_col)
data_combined <- rbind(data_train[, needed_cols], data_test[, needed_cols])

# FIX 4: complete cases BEFORE the design matrix is built
n_before <- nrow(data_combined)
data_combined <- data_combined[complete.cases(data_combined), ]
cat(sprintf("Complete cases: %d of %d rows kept (%d dropped)\n", nrow(data_combined), n_before, n_before - nrow(data_combined)))
stopifnot(all(data_combined$ARM %in% c(0, 1)))

idx_train_all <- which(data_combined$Dataset == "Train")
idx_test_all  <- which(data_combined$Dataset == "Test")
cat(sprintf("train: n = %d   test: n = %d\n", length(idx_train_all), length(idx_test_all)))
if (length(idx_train_all) < length(idx_test_all)) warning("Training pool is smaller than the test set: are the files swapped?")

# Stack to align factor levels; scale with TRAIN moments only (FIX 5)
results_combined <- standardize_X_by_index(
  X_initial = data_combined[, X_col], process_data = TRUE,
  interaction_rule = "continuous_or_binary", cat_coding_method = "sum",
  scale_idx = idx_train_all)
X_final_combined <- results_combined$X_final
storage.mode(X_final_combined) <- "double"            # stochtree's C++ layer requires doubles
data_combined$y   <- as.numeric(data_combined$y)
data_combined$ARM <- as.numeric(data_combined$ARM)
stopifnot(nrow(X_final_combined) == nrow(data_combined))
print(results_combined$X_final_var_info[, c("col_name_final", "is_continuous", "is_binary", "is_categorical")])

X_train_full <- X_final_combined[idx_train_all, , drop = FALSE]
Y_train_full <- data_combined$y[idx_train_all]
Z_train_full <- data_combined$ARM[idx_train_all]
X_test <- X_final_combined[idx_test_all, , drop = FALSE]
Y_test <- data_combined$y[idx_test_all]
Z_test <- data_combined$ARM[idx_test_all]
n_test <- length(Y_test)

pi_train <- mean(Z_train_full == 1)      # same training-pool proportion for both models (FIX 7)
cat(sprintf("treated fraction: train %.3f, test %.3f\n", pi_train, mean(Z_test == 1)))
cat(sprintf("benefit_sign = %+d; empirical ATE train %+.3f, test %+.3f (both should be negative for HbA1c)\n",
            CFG$benefit_sign, mean(Y_train_full[Z_train_full == 1]) - mean(Y_train_full[Z_train_full == 0]),
            mean(Y_test[Z_test == 1]) - mean(Y_test[Z_test == 0])))

# ==============================================================================
# MONTE CARLO SWEEP
# ==============================================================================
all_metrics <- data.frame()
metrics_file <- out_path("External_MonteCarlo_metrics_", CFG$tag, ".csv")
METRIC_COLUMNS <- c("RMSE_Y_Std", "RMSE_Y_Lin", "CATE_Correlation", "CATE_Spearman",
  "SD_CATE_Std", "SD_CATE_Lin", "Calib_MAE_Std", "Calib_MAE_Lin",
  "AUUC_Normalized_Std", "AUUC_Normalized_Lin", "NetAUUC_Normalized_Std", "NetAUUC_Normalized_Lin",
  "AUUC_CountScaled_Std", "AUUC_CountScaled_Lin", "NetAUUC_CountScaled_Std", "NetAUUC_CountScaled_Lin",
  "AUUC_Std", "AUUC_Lin", "NetAUUC_Std", "NetAUUC_Lin", "ATE_Bias_Std", "ATE_Bias_Lin",
  "ATE_Cov_Std", "ATE_Cov_Lin", "ATE_CovAdj_Std", "ATE_CovAdj_Lin")
identity_helpers <- c("check_hchs_gibbs_package", "standardize_X_by_index", "gp_semiparam", "gp_bcf", "fit_both",
  "canvas_gp", "canvas_fit_chains", "canvas_chain_seeds", "canvas_diagnostics", "canvas_parameter_group", "canvas_diagnostic_groups", "canvas_mcmc_settings", "save_fit_diagnostics",
  "canvas_draw_matrix", "canvas_coefficient_matrix", "predict_linear_bcf_patched", "get_uplift", "get_auuc", "get_net_auuc", "get_calibration", "draw_metrics")
analysis_code <- setNames(lapply(identity_helpers, function(nm) function_signature(get(nm))), identity_helpers)
analysis_code$loaded_linear_fit <- function_signature(spbart::bcf_linear_probit)
analysis_code$loaded_bcf_fit <- function_signature(spbart::bcf)
input_hashes <- setNames(unname(tools::md5sum(c(CFG$train_file, CFG$test_file))), c("train", "validation"))
if (anyNA(input_hashes)) stop("Cannot fingerprint both restricted input files.")
RUN_IDENTITY <- make_run_identity(CFG,
  package_builds = lapply(c("spbart", "dplyr", "posterior"), package_fingerprint), input_hashes = input_hashes,
  design_hash = object_md5(list(X_train_full, Y_train_full, Z_train_full, X_test, Y_test, Z_test,
                               results_combined$X_final_var_info)), analysis_code = analysis_code)
manifest_file <- out_path("canvas_run_", CFG$tag, ".rds")
existing_outputs <- list.files(CFG$out_dir, full.names = TRUE)
existing_caches <- existing_outputs[startsWith(basename(existing_outputs), paste0("canvas_preds_", CFG$tag, "_"))]
RUN_ID <- ensure_run_manifest(RUN_IDENTITY, manifest_file, c(existing_caches, metrics_file))
cat("Run ID:", RUN_ID, "(Gibbs hc-hs; normalized AUUC; identified caches only)\n")
# all_metrics starts as a zero-column frame, where filter(N_Train ...) errors
metrics_for     <- function(df, nval) if ("N_Train" %in% names(df)) df %>% filter(N_Train == nval) else df
metrics_without <- function(df, nval) if ("N_Train" %in% names(df)) df %>% filter(N_Train != nval) else df
if (CFG$resume && file.exists(metrics_file)) {                           # FIX 11
  all_metrics <- read.csv(metrics_file, stringsAsFactors = FALSE, colClasses = c(Run_ID = "character"))
  validate_metrics(all_metrics, RUN_ID, n_test, METRIC_COLUMNS)
  cat(sprintf("Resume: %d metric rows loaded from %s\n", nrow(all_metrics), metrics_file))
}

for (n in CFG$n_targets) {
  if (n >= length(Y_train_full)) { cat("Skipping N =", n, "(>= training pool)\n"); next }
  cat("\n>>> N =", n, "<<<\n")
  rds_file <- out_path("canvas_preds_", CFG$tag, "_N", n, ".rds")
  seeds <- CFG$seed_base + 1000 * n + seq_len(CFG$n_draws)
  if (CFG$resume && file.exists(rds_file)) {                             # FIX 11
    cat("  Resume: predictions found, rebuilding the figure without refitting\n")
    validate_metrics(all_metrics, RUN_ID, n_test, METRIC_COLUMNS, n, CFG$n_draws, seeds)
    preds <- load_preds(rds_file, RUN_ID, n, n_test, seeds)
    make_dashboard(preds, metrics_for(all_metrics, n), Y_test, Z_test,
                   as.character(n), ncol(preds$Std_CATE))
    next
  }
  # refitting this size: discard any partial rows a crashed run left in the CSV
  all_metrics <- metrics_without(all_metrics, n)
  n_draws <- CFG$n_draws                                                # FIX 2
  preds <- list(Std_CATE = matrix(NA_real_, n_test, n_draws), Lin_CATE = matrix(NA_real_, n_test, n_draws),
                Std_Y    = matrix(NA_real_, n_test, n_draws), Lin_Y    = matrix(NA_real_, n_test, n_draws))

  for (draw in seq_len(n_draws)) {
    seed <- CFG$seed_base + 1000 * n + draw                              # FIX 3
    cat(sprintf("  Draw %d of %d (seed %d) ...\n", draw, n_draws, seed))
    set.seed(seed)
    # stratified by arm, WITHOUT replacement (a sample-size study wants distinct patients)
    n_treated <- round(n * pi_train); n_control <- n - n_treated
    train_idx <- c(sample(which(Z_train_full == 1), n_treated), sample(which(Z_train_full == 0), n_control))

    pr <- fit_both(X_train_full[train_idx, , drop = FALSE], Y_train_full[train_idx], Z_train_full[train_idx],
                   X_test, Z_test, pi_train, seed)
    save_fit_diagnostics(pr, n, draw, seed)
    preds$Std_CATE[, draw] <- pr$Std_CATE; preds$Lin_CATE[, draw] <- pr$Lin_CATE
    preds$Std_Y[, draw]    <- pr$Std_Y;    preds$Lin_Y[, draw]    <- pr$Lin_Y

    m <- cbind(data.frame(Run_ID = RUN_ID, N_Train = n, Draw = draw, Seed = seed), draw_metrics(pr, Y_test, Z_test))
    all_metrics <- bind_rows(all_metrics, m)
    write.csv(all_metrics, metrics_file, row.names = FALSE)              # FIX 10: flush every draw
  }
  save_preds(preds, rds_file, RUN_ID, n, seeds)

  cat("  Generating dashboard for N =", n, "...\n")
  make_dashboard(preds, metrics_for(all_metrics, n), Y_test, Z_test, as.character(n), n_draws)
}

# ==============================================================================
# FULL-POOL FIT (one fit on all of CANVAS-R; backs the "whole sample" sentence)
# ==============================================================================
if (CFG$fit_full) {
  cat("\n>>> Full CANVAS-R pool (n =", length(Y_train_full), ") <<<\n")
  rds_full <- out_path("canvas_preds_", CFG$tag, "_Nfull.rds")
  if (CFG$resume && file.exists(rds_full)) {                             # FIX 11
    cat("  Resume: predictions found, rebuilding the figure without refitting\n")
    validate_metrics(all_metrics, RUN_ID, n_test, METRIC_COLUMNS, length(Y_train_full), 1, CFG$seed_base + 999999)
    preds <- load_preds(rds_full, RUN_ID, "full", n_test, CFG$seed_base + 999999)
    make_dashboard(preds, metrics_for(all_metrics, length(Y_train_full)),
                   Y_test, Z_test, "full", 1)
    CFG$fit_full <- FALSE
  }
}
if (CFG$fit_full) {
  seed <- CFG$seed_base + 999999
  pr <- fit_both(X_train_full, Y_train_full, Z_train_full, X_test, Z_test, pi_train, seed)
  save_fit_diagnostics(pr, length(Y_train_full), 1L, seed)
  m <- cbind(data.frame(Run_ID = RUN_ID, N_Train = length(Y_train_full), Draw = 1, Seed = seed), draw_metrics(pr, Y_test, Z_test))
  all_metrics <- metrics_without(all_metrics, length(Y_train_full))
  all_metrics <- bind_rows(all_metrics, m)
  write.csv(all_metrics, metrics_file, row.names = FALSE)
  preds <- list(Std_CATE = cbind(pr$Std_CATE), Lin_CATE = cbind(pr$Lin_CATE), Std_Y = cbind(pr$Std_Y), Lin_Y = cbind(pr$Lin_Y))
  save_preds(preds, rds_full, RUN_ID, "full", seed)
  make_dashboard(preds, m, Y_test, Z_test, "full", 1)
  cat(sprintf("Full-pool OOS CATE correlation: Pearson %.3f, Spearman %.3f\n", m$CATE_Correlation, m$CATE_Spearman))
}

# ==============================================================================
# PUBLICATION TABLE. "mean Ã‚Â± sd" is the spread ACROSS DRAWS (not a standard
# error). Draws at large N overlap heavily (4000 of 5018 share >= 80% of
# patients), so the sd shrinks with N partly for that reason.
# ==============================================================================
cat("\nGenerating final aggregated table...\n")
fmt <- function(x, d = 3) {
  if (length(x) > 1) sprintf(paste0("%.", d, "f Ã‚Â± %.", d, "f"), mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
  else sprintf(paste0("%.", d, "f"), x)
}

summary_table <- all_metrics %>% group_by(N_Train) %>%
  summarise(
    Draws = n(),
    `Y RMSE (BCF)` = fmt(RMSE_Y_Std), `Y RMSE (Semi-parametric BART)` = fmt(RMSE_Y_Lin),
    `CATE Correlation` = fmt(CATE_Correlation), `CATE Spearman` = fmt(CATE_Spearman),
    `SD CATE (BCF)` = fmt(SD_CATE_Std), `SD CATE (Semi-parametric BART)` = fmt(SD_CATE_Lin),
    `Calib MAE (BCF)` = fmt(Calib_MAE_Std), `Calib MAE (Semi-parametric BART)` = fmt(Calib_MAE_Lin),
    `AUUC / validation n (BCF)` = fmt(AUUC_Normalized_Std),
    `AUUC / validation n (Semi-parametric BART)` = fmt(AUUC_Normalized_Lin),
    `Net AUUC / validation n (BCF)` = fmt(NetAUUC_Normalized_Std),
    `Net AUUC / validation n (Semi-parametric BART)` = fmt(NetAUUC_Normalized_Lin),
    `Legacy count-scaled AUUC (BCF)` = fmt(AUUC_CountScaled_Std, 1),
    `Legacy count-scaled AUUC (Semi-parametric BART)` = fmt(AUUC_CountScaled_Lin, 1),
    `Legacy count-scaled net AUUC (BCF)` = fmt(NetAUUC_CountScaled_Std, 1),
    `Legacy count-scaled net AUUC (Semi-parametric BART)` = fmt(NetAUUC_CountScaled_Lin, 1),
    `ATE Bias (BCF)` = fmt(ATE_Bias_Std), `ATE Bias (Semi-parametric BART)` = fmt(ATE_Bias_Lin),
    `ATE Cov raw/adj (BCF)` = sprintf("%.2f / %.2f", mean(ATE_Cov_Std), mean(ATE_CovAdj_Std)),
    `ATE Cov raw/adj (Semi-parametric BART)` = sprintf("%.2f / %.2f", mean(ATE_Cov_Lin), mean(ATE_CovAdj_Lin)),
    `NA correlations` = sum(is.na(CATE_Correlation)),
    .groups = "drop")

print(summary_table, width = Inf)
write.csv(summary_table, out_path("External_MonteCarlo_Summary_Table_", CFG$tag, ".csv"), row.names = FALSE)
png(out_path("my_table_final_", CFG$tag, ".png"), height = 50 * (nrow(summary_table) + 1), width = 160 * ncol(summary_table))
grid.table(summary_table); dev.off()

# ==============================================================================
# COMPOUND FIGURE: the sizes the manuscript shows side by side, as ONE file
# (journal: compound figures uploaded as one figure). Rebuilt from the .rds
# files, so `make_compound_figure()` can be rerun offline with other sizes or
# panels ("outcome", "cate", "calibration", "uplift") without refitting.
# ==============================================================================
cat("\nBuilding compound figure...\n")
make_compound_figure(CFG$compound_sizes, panels = CFG$compound_panels, metrics_all = all_metrics,
                     Y_test = Y_test, Z_test = Z_test, fig_number = CFG$fig_numbers[["compound"]])
cat("Done.\n")
