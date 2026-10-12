# ==============================================================================
# COMPLETE INTERPRETATION DIAGNOSTICS FOR THE CANVAS-R hc-hs FIT
# ==============================================================================
#
# This one-file script DOES NOT refit the model and DOES NOT modify the
# manuscript. It reads the files written by canvas_prepare.R / first (2).R and
# produces:
#
#   shapley_grouped_importance.png
#   shapley_distribution.png                (distribution plot from methodology-7)
#   shapley_grouped_importance.csv
#   shapley_design_column_importance.csv
#   posterior_coefficients.csv
#   posterior_BLA1C_and_RACE.csv
#   posterior_BLA1C_beta.png
#   posterior_BLA1C_beta.csv
#   posterior_BLA1C_contrast.png
#   posterior_BLA1C_contrast.csv
#   posterior_AGE_beta.png
#   posterior_AGE_beta.csv
#   posterior_AGE_contrast.png
#   posterior_AGE_contrast.csv
#   posterior_interactions.csv
#   mcmc_diagnostics_main.csv
#   mcmc_diagnostics_interactions.csv
#   posterior_trace_main_coefficients.png
#   all_covariate_contrasts.csv
#   all_covariate_contrasts.png
#   all_covariate_evidence.csv
#   classical_interaction_sensitivity.csv       (when analsis.RData is present)
#   race_level_main_effect_deviations.csv
#   race_pairwise_main_effect_contrasts.csv
#   race_unadjusted_treatment_effects.csv
#   race_counts.csv                         (when analysis files are present)
#   design_diagnostics.csv
#   diagnostic_notes.txt
#   interpretation_summary.txt
#
# The important distinction is:
#   * Shapley bar heights are descriptive importance measures, not tests.
#   * RACE1, ..., RACE(K-1) are sum-contrast coordinates, not race categories.
#   * RACE must be combined BEFORE taking absolute Shapley values.
#
# YODA example (PowerShell):
#   $env:BCF_WORKSPACE = "X:/Workspace"
#   $env:BCF_SHAPLEY_SRC = "X:/Workspace/shapley.R"
#   Rscript "X:/Workspace/interpretation_all_covariates.R"
#
# Required input files in BCF_WORKSPACE:
#   X_final.RData, nbcf_fit.RData, and shapley.R
# Recommended:
#   X_final_varinfo.RData, analsis.RData, analsis_2.RData
# ==============================================================================

# Shared sampler protocol; keep canvas_mcmc.R beside this script.
canvas_source_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(canvas_source_path)) {
  canvas_cli_path <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(canvas_cli_path)) canvas_source_path <- sub("^--file=", "", canvas_cli_path[1])
}
CANVAS_CODE_DIR <- if (is.null(canvas_source_path)) getwd() else dirname(normalizePath(canvas_source_path, mustWork = TRUE))
source(file.path(CANVAS_CODE_DIR, "canvas_mcmc.R"))
if (!requireNamespace("posterior", quietly = TRUE)) stop("Install the posterior R package before the CANVAS rerun")

INPUT_WORKSPACE <- Sys.getenv("BCF_WORKSPACE", unset = "X:/Workspace")
WORKSPACE_DIR <- file.path(canvas_run_root(INPUT_WORKSPACE), "interpretation")
LOW_SPACE <- tolower(Sys.getenv("BCF_LOW_SPACE", unset = "true")) %in%
  c("1", "true", "yes", "y")
SAVE_TRACE_PLOT <- tolower(Sys.getenv(
  "BCF_SAVE_TRACE_PLOT", unset = if (LOW_SPACE) "false" else "true"
)) %in% c("1", "true", "yes", "y")
PLOT_DPI <- if (LOW_SPACE) 150 else 300
EXPECTED_COVARIATES <- c(
  "BLBMI", "AGE", "HHFFL", "BLSBP", "SEX", "DIADUR", "HCEREBFL", "BLA1C",
  "BLGFR", "INSULBL", "GLP1BL", "BLACR", "RACE", "METFORBL", "RANGIBL"
)

script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(script_file_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[1]), mustWork = FALSE))
} else {
  getwd()
}

required_packages <- c("dplyr", "tidyr", "purrr", "ggplot2", "forcats", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the following R package(s) first: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

say <- function(...) {
  cat(sprintf(...))
  flush.console()
}

load_object <- function(path, object_name, required = TRUE) {
  if (!file.exists(path)) {
    if (required) stop("Cannot find ", path, call. = FALSE)
    return(NULL)
  }
  e <- new.env(parent = emptyenv())
  loaded <- load(path, envir = e)
  if (!object_name %in% loaded) {
    if (required) {
      stop(
        basename(path), " does not contain an object called '", object_name, "'.",
        call. = FALSE
      )
    }
    return(NULL)
  }
  e[[object_name]]
}

first_existing <- function(paths) {
  paths <- unique(paths[nzchar(paths)])
  hit <- paths[file.exists(paths)]
  if (length(hit)) hit[1] else NA_character_
}

posterior_matrix <- function(x, label) {
  d <- dim(x)
  if (length(d) == 3L) {
    # Expected layout: chain x draw x coefficient.
    return(do.call(rbind, lapply(seq_len(d[1]), function(chain) {
      matrix(x[chain, , , drop = FALSE], nrow = d[2], ncol = d[3])
    })))
  }
  if (length(d) == 2L) return(as.matrix(x))
  stop(label, " must be a 2D or 3D posterior array.", call. = FALSE)
}

subset_beta_columns <- function(beta, columns) {
  if (length(dim(beta)) == 3L) {
    beta[, , columns, drop = FALSE]
  } else {
    beta[, columns, drop = FALSE]
  }
}

posterior_summary <- function(draw_matrix, term_names, simultaneous = FALSE,
                              credible_level = 0.95) {
  stopifnot(ncol(draw_matrix) == length(term_names))
  ans <- data.frame(
    term = term_names,
    mean = colMeans(draw_matrix),
    median = apply(draw_matrix, 2, stats::median),
    sd = apply(draw_matrix, 2, stats::sd),
    lower_95 = apply(draw_matrix, 2, stats::quantile, probs = 0.025, names = FALSE),
    upper_95 = apply(draw_matrix, 2, stats::quantile, probs = 0.975, names = FALSE),
    probability_positive = colMeans(draw_matrix > 0),
    probability_negative = colMeans(draw_matrix < 0),
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(
      probability_same_direction = pmax(probability_positive, probability_negative),
      pointwise_95_interval_excludes_zero = lower_95 > 0 | upper_95 < 0
    )

  if (isTRUE(simultaneous)) {
    centre <- colMeans(draw_matrix)
    spread <- apply(draw_matrix, 2, stats::sd)
    active <- is.finite(spread) & spread > sqrt(.Machine$double.eps)
    if (any(active)) {
      z <- sweep(draw_matrix[, active, drop = FALSE], 2, centre[active], "-")
      z <- sweep(z, 2, spread[active], "/")
      max_z <- apply(abs(z), 1, max)
      critical <- unname(stats::quantile(max_z, credible_level, names = FALSE))
    } else {
      critical <- 0
    }
    ans$simultaneous_critical_value <- critical
    ans$simultaneous_lower_95 <- centre - critical * spread
    ans$simultaneous_upper_95 <- centre + critical * spread
    ans$simultaneous_95_interval_excludes_zero <-
      ans$simultaneous_lower_95 > 0 | ans$simultaneous_upper_95 < 0
  }
  ans
}

posterior_chains <- function(x, columns = NULL) {
  d <- dim(x)
  if (length(d) == 3L) {
    if (is.null(columns)) columns <- seq_len(d[3])
    return(lapply(seq_len(d[1]), function(chain) {
      matrix(x[chain, , columns, drop = FALSE], nrow = d[2], ncol = length(columns))
    }))
  }
  if (length(d) == 2L) {
    if (is.null(columns)) columns <- seq_len(d[2])
    return(list(as.matrix(x[, columns, drop = FALSE])))
  }
  stop("Posterior object must be a 2D or 3D array.", call. = FALSE)
}

mcmc_diagnostics <- function(x, term_names, columns = NULL) {
  chains <- posterior_chains(x, columns)
  stopifnot(ncol(chains[[1]]) == length(term_names), length(chains) == 4L)
  chains <- lapply(chains, function(value) { colnames(value) <- term_names; value })
  settings <- canvas_mcmc_settings(); settings$draws <- nrow(chains[[1]])
  result <- canvas_diagnostics(chains, settings)
  result$term <- result$variable
  result$rhat_assessable <- result$chains >= 2L
  result$diagnostic_flag <- !result$pass
  result
}

normal_credible_ellipsoid <- function(draw_matrix, level = 0.95) {
  draw_matrix <- as.matrix(draw_matrix)
  keep <- apply(draw_matrix, 2, stats::sd) > sqrt(.Machine$double.eps)
  draw_matrix <- draw_matrix[, keep, drop = FALSE]
  if (!ncol(draw_matrix)) {
    return(c(dimension = 0, mahalanobis_sq = 0, cutoff = 0,
             zero_outside = FALSE, tail_area = 1))
  }
  mu <- colMeans(draw_matrix)
  s <- stats::cov(draw_matrix)
  if (ncol(draw_matrix) == 1L) s <- matrix(s, 1, 1)
  eig <- eigen(s, symmetric = TRUE)
  tol <- max(eig$values) * max(dim(s)) * .Machine$double.eps^0.5
  positive <- eig$values > tol
  if (!any(positive)) {
    return(c(dimension = 0, mahalanobis_sq = NA, cutoff = NA,
             zero_outside = NA, tail_area = NA))
  }
  inv_s <- eig$vectors[, positive, drop = FALSE] %*%
    diag(1 / eig$values[positive], nrow = sum(positive)) %*%
    t(eig$vectors[, positive, drop = FALSE])
  d2 <- as.numeric(t(mu) %*% inv_s %*% mu)
  df <- sum(positive)
  cutoff <- stats::qchisq(level, df = df)
  c(dimension = df, mahalanobis_sq = d2, cutoff = cutoff,
    zero_outside = d2 > cutoff,
    tail_area = stats::pchisq(d2, df = df, lower.tail = FALSE))
}

write_csv <- function(x, filename) {
  path <- file.path(WORKSPACE_DIR, filename)
  utils::write.csv(x, path, row.names = FALSE, na = "")
  say("  -> %s\n", path)
}

plot_posterior <- function(draws, title, x_label, filename) {
  draws <- as.numeric(draws)
  if (length(draws) < 2L || any(!is.finite(draws))) {
    stop("Posterior plot requires at least two finite draws: ", filename,
         call. = FALSE)
  }
  qs <- unname(stats::quantile(draws, c(0.025, 0.5, 0.975)))
  p <- ggplot2::ggplot(data.frame(x = draws), ggplot2::aes(x = x))
  if (stats::sd(draws) > 0) {
    p <- p + ggplot2::geom_density(fill = "skyblue", alpha = 0.5, colour = "navy")
  } else {
    p <- p + ggplot2::geom_vline(xintercept = draws[1], colour = "navy")
  }
  p <- p +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    ggplot2::annotate("segment", x = qs[1], xend = qs[3], y = 0, yend = 0,
                      linewidth = 2, colour = "navy") +
    ggplot2::annotate("point", x = qs[2], y = 0, size = 4, colour = "firebrick") +
    ggplot2::labs(
      title = title,
      subtitle = sprintf(
        "mean %+.3f, median %+.3f\npointwise 95%% CrI [%+.3f, %+.3f], P(< 0) = %.3f",
        mean(draws), qs[2], qs[1], qs[3], mean(draws < 0)
      ),
      x = x_label, y = "Density"
    ) +
    ggplot2::theme_minimal(base_size = 13)
  plot_path <- file.path(WORKSPACE_DIR, filename)
  ggplot2::ggsave(plot_path, p, width = 7, height = 5, dpi = PLOT_DPI, bg = "white")
  say("  -> %s\n", plot_path)
}

say("============================================================\n")
say("Complete all-covariate interpretation diagnostics (no model refit)\n")
say("Workspace: %s\n", WORKSPACE_DIR)
say("Low-space mode: %s%s\n", LOW_SPACE,
    if (SAVE_TRACE_PLOT) "" else " (trace diagnostics remain in CSV; trace PNG skipped)")
say("============================================================\n\n")

# ---- Load the saved design and posterior -------------------------------------
X_final <- load_object(file.path(WORKSPACE_DIR, "X_final.RData"), "X_final")
my_list <- load_object(file.path(WORKSPACE_DIR, "nbcf_fit.RData"), "my_list")
if (!identical(my_list$mcmc$protocol, "canvas-independent-4chain-mu-propensity-8.3.5") ||
    !isTRUE(my_list$mcmc$independent_package_calls) ||
    !identical(my_list$model_params$propensity_covariate, "mu") ||
    !identical(my_list$model_params$propensity_seperate, "none") ||
    !identical(my_list$mcmc$settings$chains, 4L) ||
    !identical(my_list$mcmc$settings$draws, 4000L) ||
    !identical(my_list$mcmc$settings$gfr, 25L) ||
    !identical(my_list$mcmc$settings$burnin, 2000L) ||
    length(unique(my_list$mcmc$seeds)) != 4L)
  stop("Interpretation requires the new four-independent-chain CANVAS fit with propensity in the prognostic forest only; historical results cannot be reused.")
res_X <- load_object(
  file.path(WORKSPACE_DIR, "X_final_varinfo.RData"), "res_X", required = FALSE
)

if (is.null(colnames(X_final))) {
  if (!is.null(my_list$col_names) && length(my_list$col_names) == ncol(X_final)) {
    colnames(X_final) <- my_list$col_names
  } else {
    stop("X_final has no column names and my_list$col_names cannot restore them.", call. = FALSE)
  }
}
if (!is.null(my_list$col_names) && !identical(as.character(my_list$col_names), colnames(X_final))) {
  stop("my_list$col_names does not match colnames(X_final); these files are from different fits.",
       call. = FALSE)
}
if (is.null(my_list$beta) || is.null(my_list$beta_int) ||
    is.null(my_list$interaction_pairs)) {
  stop(
    "nbcf_fit.RData must contain my_list$beta, $beta_int, and $interaction_pairs. ",
    "Re-run the updated canvas_prepare.R if these fields are missing.",
    call. = FALSE
  )
}

B_all <- posterior_matrix(my_list$beta, "my_list$beta")
beta_offset <- ncol(B_all) - ncol(X_final)
if (!beta_offset %in% c(0L, 1L)) {
  stop(
    sprintf(
      "Cannot align %d beta columns with the %d design columns.",
      ncol(B_all), ncol(X_final)
    ),
    call. = FALSE
  )
}
beta_names <- if (beta_offset == 1L) c("(Intercept)", colnames(X_final)) else colnames(X_final)
colnames(B_all) <- beta_names
design_beta_columns <- seq.int(beta_offset + 1L, beta_offset + ncol(X_final))
B_design <- B_all[, design_beta_columns, drop = FALSE]
colnames(B_design) <- colnames(X_final)
say("Loaded %d posterior draws, %d design columns%s.\n",
    nrow(B_all), ncol(X_final), if (beta_offset == 1L) ", plus an intercept" else "")

# ---- Posterior coefficient summaries ----------------------------------------
# These are summaries of main-effect coefficients on the fitted model scale.
# Pointwise intervals are reported as diagnostics, not as multiplicity-adjusted
# declarations of significance.
coef_summary <- posterior_summary(B_all, beta_names, simultaneous = TRUE)
write_csv(coef_summary, "posterior_coefficients.csv")

focus_summary <- coef_summary |>
  dplyr::filter(term == "BLA1C" | grepl("^RACE[0-9]+$", term))
write_csv(focus_summary, "posterior_BLA1C_and_RACE.csv")

# Resolve named coefficients after the intercept offset and chains have been aligned.
# Keep the fitted model scale: legacy my_list exports do not save outcome_scale.
if ("BLA1C" %in% colnames(B_design)) {
  write_csv(dplyr::filter(coef_summary, term == "BLA1C"),
            "posterior_BLA1C_beta.csv")
  plot_posterior(
    B_design[, "BLA1C"], expression(paste("Posterior of ", beta[BLA1C])),
    "Baseline HbA1c main-effect coefficient (fitted model scale)",
    "posterior_BLA1C_beta.png"
  )
} else {
  warning("BLA1C is absent from the aligned beta/design columns; dedicated HbA1c beta outputs were not written.",
          call. = FALSE)
}

if ("AGE" %in% colnames(B_design)) {
  write_csv(dplyr::filter(coef_summary, term == "AGE"), "posterior_AGE_beta.csv")
  plot_posterior(
    B_design[, "AGE"], expression(paste("Posterior of ", beta[AGE])),
    "AGE main-effect coefficient (fitted model scale)", "posterior_AGE_beta.png"
  )
}

main_mcmc <- mcmc_diagnostics(my_list$beta, beta_names)
write_csv(main_mcmc, "mcmc_diagnostics_main.csv")

ipairs <- my_list$interaction_pairs
if (nrow(ipairs) != 2L && ncol(ipairs) == 2L) ipairs <- t(ipairs)
if (nrow(ipairs) != 2L) {
  stop("interaction_pairs must have two rows (or two columns).", call. = FALSE)
}
if (length(ipairs) &&
    (any(ipairs < 1L) || any(ipairs > ncol(X_final)))) {
  stop("interaction_pairs contains indices outside the X_final columns.", call. = FALSE)
}

B_interaction <- posterior_matrix(my_list$beta_int, "my_list$beta_int")
if (ncol(B_interaction) != ncol(ipairs)) {
  stop(
    sprintf(
      "beta_int has %d columns but interaction_pairs has %d pairs.",
      ncol(B_interaction), ncol(ipairs)
    ),
    call. = FALSE
  )
}
if (nrow(B_interaction) != nrow(B_design)) {
  stop("beta and beta_int contain different numbers of posterior draws.", call. = FALSE)
}
interaction_names <- if (ncol(ipairs)) {
  paste(colnames(X_final)[ipairs[1, ]], colnames(X_final)[ipairs[2, ]], sep = " x ")
} else {
  character(0)
}
colnames(B_interaction) <- interaction_names
if (length(interaction_names)) {
  interaction_summary <- posterior_summary(
    B_interaction, interaction_names, simultaneous = TRUE
  )
  write_csv(interaction_summary, "posterior_interactions.csv")
  interaction_mcmc <- mcmc_diagnostics(my_list$beta_int, interaction_names)
  write_csv(interaction_mcmc, "mcmc_diagnostics_interactions.csv")
}

# One compact trace figure for all stored main coefficients and all four chains.
# At most 2,000 draws per chain are plotted for file size; diagnostics use all draws.
main_chain_list <- posterior_chains(my_list$beta)
trace_data <- dplyr::bind_rows(lapply(seq_along(main_chain_list), function(ch) {
  z <- main_chain_list[[ch]]
  step <- max(1L, ceiling(nrow(z) / 2000L))
  keep <- seq.int(1L, nrow(z), by = step)
  as.data.frame(z[keep, , drop = FALSE]) |>
    stats::setNames(beta_names) |>
    dplyr::mutate(chain = factor(ch), iteration = keep, .before = 1) |>
    tidyr::pivot_longer(cols = -c(chain, iteration), names_to = "term", values_to = "draw")
}))
p_trace <- ggplot2::ggplot(
  trace_data, ggplot2::aes(x = iteration, y = draw, colour = chain)
) +
  ggplot2::geom_line(linewidth = 0.25, alpha = 0.75) +
  ggplot2::facet_wrap(~ term, scales = "free_y", ncol = 4) +
  ggplot2::labs(x = "Iteration", y = NULL, title = NULL, subtitle = NULL,
                colour = "Chain") +
  ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position = "bottom")
if (SAVE_TRACE_PLOT) {
  trace_path <- file.path(WORKSPACE_DIR, "posterior_trace_main_coefficients.png")
  ggplot2::ggsave(trace_path, p_trace, width = 10, height = 8,
                  dpi = if (LOW_SPACE) 120 else 220, bg = "white")
  say("  -> %s\n", trace_path)
} else {
  say("  -> trace PNG skipped to conserve space; numeric MCMC diagnostics were saved\n")
}

# ---- Source the same Shapley implementation used for the original plot ------
configured_shapley <- Sys.getenv("BCF_SHAPLEY_SRC", unset = "")
shapley_src <- first_existing(c(
  configured_shapley,
  file.path(CANVAS_CODE_DIR, "shapley_aux.R"),
  file.path(INPUT_WORKSPACE, "shapley.R"),
  file.path(SCRIPT_DIR, "shapley.R"),
  file.path(SCRIPT_DIR, "shapley_aux.R")
))
if (is.na(shapley_src)) {
  stop(
    "Cannot find shapley.R. Set BCF_SHAPLEY_SRC to the file used by canvas_prepare.R.",
    call. = FALSE
  )
}
say("Sourcing Shapley functions from %s\n", shapley_src)
source(shapley_src)
if (!exists("compute_shapley_all", mode = "function")) {
  stop("The Shapley source does not define compute_shapley_all().", call. = FALSE)
}

# compute_shapley_all() expects coefficients for the actual design columns. If
# the sampler stored a leading intercept, remove it here rather than silently
# assigning it to the first feature.
beta_for_shapley <- subset_beta_columns(my_list$beta, design_beta_columns)

shap_args <- list(
  X = X_final,
  beta_post = beta_for_shapley,
  beta_int_post = my_list$beta_int,
  ipairs = my_list$interaction_pairs,
  indices = seq_len(nrow(X_final))
)
shapley_formals <- names(formals(compute_shapley_all))
if ("boolean_interaction" %in% shapley_formals) {
  if (is.null(res_X) || is.null(res_X$X_final_var_info)) {
    stop(
      "This version of compute_shapley_all() needs boolean_interaction, but ",
      "X_final_varinfo.RData is missing.",
      call. = FALSE
    )
  }
  shap_args$boolean_interaction <- as.logical(
    res_X$X_final_var_info$is_continuous | res_X$X_final_var_info$is_binary
  )
}
extra_args <- setdiff(names(shap_args), shapley_formals)
if (length(extra_args)) shap_args[extra_args] <- NULL

say("Recomputing posterior-median Shapley values from the saved fit...\n")
shapley_results <- do.call(compute_shapley_all, shap_args)
needed_shapley <- c("obs", "feature", "shapley_total", "shapley_main", "shapley_int")
if (!all(needed_shapley %in% names(shapley_results))) {
  stop(
    "compute_shapley_all() did not return: ",
    paste(setdiff(needed_shapley, names(shapley_results)), collapse = ", "),
    call. = FALSE
  )
}

# Preserve the original design-column ranking for auditability.
column_importance <- shapley_results |>
  dplyr::group_by(feature) |>
  dplyr::summarise(
    main = mean(abs(shapley_main)),
    interaction = mean(abs(shapley_int)),
    total = main + interaction,
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(total))
write_csv(column_importance, "shapley_design_column_importance.csv")

# RACE1, RACE2, ... are coordinates of one categorical covariate. Sum their
# signed, observation-level contributions first; only then take absolute values.
shapley_grouped <- shapley_results |>
  dplyr::mutate(
    feature = ifelse(grepl("^RACE[0-9]+$", feature), "RACE (joint)", feature)
  ) |>
  dplyr::group_by(obs, feature) |>
  dplyr::summarise(
    shapley_total = sum(shapley_total),
    shapley_main = sum(shapley_main),
    shapley_int = sum(shapley_int),
    .groups = "drop"
  )

grouped_importance <- shapley_grouped |>
  dplyr::group_by(feature) |>
  dplyr::summarise(
    main = mean(abs(shapley_main)),
    interaction = mean(abs(shapley_int)),
    total = main + interaction,
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(total))
write_csv(grouped_importance, "shapley_grouped_importance.csv")

plot_data <- grouped_importance |>
  dplyr::mutate(feature = factor(feature, levels = rev(feature))) |>
  tidyr::pivot_longer(
    cols = c("main", "interaction"), names_to = "effect", values_to = "importance"
  ) |>
  dplyr::mutate(
    effect = factor(effect, levels = c("main", "interaction"),
                    labels = c("Main effect", "Interaction"))
  )

p_grouped <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = importance, y = feature, fill = effect)
) +
  ggplot2::geom_col(width = 0.78) +
  ggplot2::scale_fill_manual(
    values = c("Main effect" = "#5a5a5a", "Interaction" = "#d9d9d9"),
    name = NULL
  ) +
  ggplot2::labs(
    x = "Mean absolute Shapley contribution",
    y = NULL,
    title = NULL,
    subtitle = NULL
  ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(8, 12, 8, 8)
  )

grouped_plot_path <- file.path(WORKSPACE_DIR, "shapley_grouped_importance.png")
ggplot2::ggsave(
  grouped_plot_path, p_grouped, width = 8.2, height = 6.3,
  dpi = PLOT_DPI, bg = "white"
)
say("  -> %s\n", grouped_plot_path)

# ---- Individual Shapley distributions (methodology-7, Figure 3) -------------
# Keep the original design-column display, including separate RACE contrasts.
# Each point uses one individual's saved-fit Shapley summary; its colour is the
# observed feature value. This reuses the calculation above without a refit.
distribution_data <- shapley_results
if (!"feature_value" %in% names(distribution_data)) {
  feature_column <- match(
    as.character(distribution_data$feature),
    make.names(colnames(X_final), unique = TRUE)
  )
  if (anyNA(feature_column) || anyNA(distribution_data$obs) ||
      any(!distribution_data$obs %in% seq_len(nrow(X_final)))) {
    stop("Cannot align Shapley points with X_final for the distribution plot.",
         call. = FALSE)
  }
  distribution_data$feature_value <- X_final[
    cbind(distribution_data$obs, feature_column)
  ]
}
distribution_data <- distribution_data |>
  dplyr::group_by(feature) |>
  dplyr::mutate(scaled_value = scales::rescale(feature_value)) |>
  dplyr::ungroup()

p_distribution <- ggplot2::ggplot(
  distribution_data,
  ggplot2::aes(
    x = shapley_total,
    y = forcats::fct_reorder(feature, abs(shapley_total), .fun = mean)
  )
) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                    colour = "black", linewidth = 0.7) +
  ggplot2::geom_point(
    ggplot2::aes(colour = scaled_value), alpha = 0.7, size = 1.5,
    position = ggplot2::position_jitter(width = 0, height = 0.15, seed = 1848)
  ) +
  ggplot2::scale_colour_gradient2(
    low = "#2c7bb6", mid = "#ffffbf", high = "#d7191c", midpoint = 0.5,
    name = "Feature Value\n(Low to High)", breaks = c(0, 1),
    labels = c("Low", "High")
  ) +
  ggplot2::labs(
    title = "SHAP Summary: Feature Contribution & Distribution",
    x = "Total Shapley Value (Impact on Prediction)", y = "Feature"
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

distribution_plot_path <- file.path(WORKSPACE_DIR, "shapley_distribution.png")
ggplot2::ggsave(
  distribution_plot_path, p_distribution, width = 10, height = 8,
  dpi = PLOT_DPI, bg = "white"
)
say("  -> %s\n", distribution_plot_path)

# ---- RACE category counts and coding-invariant category contrasts ------------
train_path <- file.path(INPUT_WORKSPACE, "analsis.RData")
test_path <- file.path(INPUT_WORKSPACE, "analsis_2.RData")
train_data <- load_object(train_path, "analsis", required = FALSE)
test_data <- load_object(test_path, "analsis", required = FALSE)

count_one <- function(d, dataset_name) {
  if (is.null(d) || !"RACE" %in% names(d)) return(NULL)
  as.data.frame(table(RACE = as.character(d$RACE), useNA = "ifany"),
                stringsAsFactors = FALSE) |>
    dplyr::mutate(dataset = dataset_name, .before = 1)
}
race_counts <- dplyr::bind_rows(
  count_one(train_data, "CANVAS-R / interpretability fit"),
  count_one(test_data, "CANVAS / external validation")
)
if (nrow(race_counts)) write_csv(race_counts, "race_counts.csv")

# ---- Original-covariate map --------------------------------------------------
# res_X records which original input generated each model-matrix column. This
# makes RACE one feature with K-1 coordinates rather than K-1 opportunities to
# be selected. A name-based fallback is retained for older saved fits.
if (!is.null(res_X) && !is.null(res_X$X_final_var_info) &&
    nrow(res_X$X_final_var_info) == ncol(X_final)) {
  feature_map <- res_X$X_final_var_info
  feature_map$design_column <- colnames(X_final)
  feature_map$feature <- EXPECTED_COVARIATES[feature_map$original_idx]
  if (anyNA(feature_map$feature)) {
    stop("X_final_var_info contains an original_idx outside EXPECTED_COVARIATES.",
         call. = FALSE)
  }
} else {
  feature_map <- data.frame(
    original_idx = match(
      ifelse(grepl("^RACE[0-9]+$", colnames(X_final)), "RACE", colnames(X_final)),
      EXPECTED_COVARIATES
    ),
    is_continuous = !grepl("^(RACE[0-9]+|HHFFL|SEX|INSULBL|GLP1BL|METFORBL|RANGIBL)$",
                           colnames(X_final)),
    is_binary = grepl("^(HHFFL|SEX|INSULBL|GLP1BL|METFORBL|RANGIBL)$",
                      colnames(X_final)),
    is_categorical = grepl("^RACE[0-9]+$", colnames(X_final)),
    design_column = colnames(X_final),
    feature = ifelse(grepl("^RACE[0-9]+$", colnames(X_final)),
                     "RACE", colnames(X_final)),
    stringsAsFactors = FALSE
  )
  warning("X_final_varinfo.RData was unavailable; feature types were inferred from names.",
          call. = FALSE)
}

recover_factor_levels <- function(feature, k) {
  if (!is.null(train_data) && feature %in% names(train_data)) {
    lev <- levels(factor(train_data[[feature]]))
    if (length(lev) == k) return(lev)
  }
  paste(feature, "level", seq_len(k))
}

surface_contrast_draws <- function(X_a, X_b) {
  delta_main <- colMeans(X_a - X_b)
  main_part <- as.numeric(B_design %*% delta_main)
  if (!ncol(ipairs)) return(main_part)
  delta_int <- vapply(seq_len(ncol(ipairs)), function(j) {
    a <- ipairs[1, j]
    b <- ipairs[2, j]
    mean(X_a[, a] * X_a[, b] - X_b[, a] * X_b[, b])
  }, numeric(1))
  main_part + as.numeric(B_interaction %*% delta_int)
}

# ---- Multiplicity-aware model contrasts for every original covariate --------
# For continuous variables: average fitted CATE at the 90th versus 10th
# percentile, holding the empirical distribution of all other variables fixed.
# For binary variables: higher-coded versus lower-coded value.
# For categorical variables: every valid pair of levels. RACE is consequently
# tested as one factor, never as separate RACE1/RACE2/... variables.
contrast_draw_list <- list()
contrast_meta_list <- list()
add_contrast <- function(draws, feature, label, contrast_type) {
  id <- length(contrast_draw_list) + 1L
  contrast_draw_list[[id]] <<- as.numeric(draws)
  contrast_meta_list[[id]] <<- data.frame(
    contrast_id = sprintf("C%03d", id),
    feature = feature,
    contrast = label,
    contrast_type = contrast_type,
    stringsAsFactors = FALSE
  )
}

observed_features <- EXPECTED_COVARIATES[EXPECTED_COVARIATES %in% feature_map$feature]
for (feature in observed_features) {
  cols <- which(feature_map$feature == feature)
  info <- feature_map[cols, , drop = FALSE]
  X_a <- X_final
  X_b <- X_final

  if (all(info$is_continuous) && length(cols) == 1L) {
    q <- unname(stats::quantile(X_final[, cols], c(0.10, 0.90), names = FALSE))
    X_a[, cols] <- q[2]
    X_b[, cols] <- q[1]
    raw_label <- "90th percentile minus 10th percentile"
    if (!is.null(train_data) && feature %in% names(train_data) &&
        is.numeric(train_data[[feature]])) {
      qr <- unname(stats::quantile(train_data[[feature]], c(0.10, 0.90),
                                   na.rm = TRUE, names = FALSE))
      raw_label <- sprintf("90th minus 10th percentile (%.4g minus %.4g in raw units)",
                           qr[2], qr[1])
    }
    add_contrast(surface_contrast_draws(X_a, X_b), feature, raw_label,
                 "continuous partial-dependence contrast")
  } else if (all(info$is_binary) && length(cols) == 1L) {
    values <- sort(unique(X_final[, cols]))
    if (length(values) != 2L) next
    X_a[, cols] <- values[2]
    X_b[, cols] <- values[1]
    level_label <- "higher-coded level minus lower-coded level"
    if (!is.null(train_data) && feature %in% names(train_data)) {
      raw_values <- sort(unique(as.character(train_data[[feature]])))
      if (length(raw_values) == 2L) {
        level_label <- paste(raw_values[2], "minus", raw_values[1])
      }
    }
    add_contrast(surface_contrast_draws(X_a, X_b), feature, level_label,
                 "binary partial-dependence contrast")
  } else if (all(info$is_categorical)) {
    k <- length(cols) + 1L
    coding <- stats::contr.sum(k)
    level_names <- recover_factor_levels(feature, k)
    pairs <- utils::combn(seq_len(k), 2L)
    for (j in seq_len(ncol(pairs))) {
      a <- pairs[1, j]
      b <- pairs[2, j]
      X_a[, cols] <- matrix(coding[a, ], nrow(X_final), length(cols), byrow = TRUE)
      X_b[, cols] <- matrix(coding[b, ], nrow(X_final), length(cols), byrow = TRUE)
      add_contrast(
        surface_contrast_draws(X_a, X_b), feature,
        paste(level_names[a], "minus", level_names[b]),
        "categorical pairwise partial-dependence contrast"
      )
    }
  }
}

if (!length(contrast_draw_list)) {
  stop("No valid original-covariate contrasts could be constructed.", call. = FALSE)
}
all_contrast_draws <- do.call(cbind, contrast_draw_list)
contrast_meta <- dplyr::bind_rows(contrast_meta_list)
colnames(all_contrast_draws) <- contrast_meta$contrast_id
all_contrasts <- dplyr::bind_cols(
  contrast_meta,
  posterior_summary(
    all_contrast_draws, contrast_meta$contrast_id, simultaneous = TRUE
  ) |>
    dplyr::select(-term)
)
write_csv(all_contrasts, "all_covariate_contrasts.csv")

# These posteriors include interactions and are the estimands used for the
# feature-level classifications. Retain their global simultaneous intervals in
# the dedicated CSVs as well as the all-covariate table.
bla1c_contrast_rows <- which(contrast_meta$feature == "BLA1C")
if (length(bla1c_contrast_rows) == 1L) {
  write_csv(all_contrasts[bla1c_contrast_rows, , drop = FALSE],
            "posterior_BLA1C_contrast.csv")
  plot_posterior(
    all_contrast_draws[, bla1c_contrast_rows],
    "Posterior of the baseline HbA1c treatment-effect contrast",
    "CATE at baseline HbA1c's 90th minus 10th percentile (fitted model scale)",
    "posterior_BLA1C_contrast.png"
  )
} else {
  warning("Expected exactly one continuous BLA1C contrast; dedicated HbA1c contrast outputs were not written.",
          call. = FALSE)
}

age_contrast_rows <- which(contrast_meta$feature == "AGE")
if (length(age_contrast_rows) == 1L) {
  write_csv(all_contrasts[age_contrast_rows, , drop = FALSE],
            "posterior_AGE_contrast.csv")
  plot_posterior(
    all_contrast_draws[, age_contrast_rows], "Posterior of the AGE treatment-effect contrast",
    "CATE at age's 90th minus 10th percentile (fitted model scale)",
    "posterior_AGE_contrast.png"
  )
}

# Feature-level evidence table. The primary evidence flag is based on the 95%
# simultaneous credible intervals across ALL contrasts in the preceding table.
# It is deliberately more conservative than looking at many pointwise intervals.
feature_evidence <- all_contrasts |>
  dplyr::group_by(feature) |>
  dplyr::summarise(
    number_of_contrasts = dplyr::n(),
    any_pointwise_95_interval_excludes_zero =
      any(pointwise_95_interval_excludes_zero),
    any_global_simultaneous_95_interval_excludes_zero =
      any(simultaneous_95_interval_excludes_zero),
    strongest_posterior_direction_probability = max(probability_same_direction),
    evidence_classification = dplyr::case_when(
      any(simultaneous_95_interval_excludes_zero) ~
        "evidence after global simultaneous 95% adjustment",
      any(pointwise_95_interval_excludes_zero) ~
        "pointwise evidence only; not robust to global multiplicity",
      TRUE ~ "no 95% posterior contrast evidence"
    ),
    .groups = "drop"
  )

# Add approximate multivariate joint checks. These are secondary because the
# hc-hs posterior can be non-Gaussian; the contrast-based flag above is primary.
joint_checks <- dplyr::bind_rows(lapply(observed_features, function(feature) {
  main_cols <- which(feature_map$feature == feature)
  int_cols <- if (ncol(ipairs)) which(
    feature_map$feature[ipairs[1, ]] == feature |
      feature_map$feature[ipairs[2, ]] == feature
  ) else integer(0)
  main_check <- normal_credible_ellipsoid(B_design[, main_cols, drop = FALSE])
  surface_draws <- cbind(
    B_design[, main_cols, drop = FALSE],
    B_interaction[, int_cols, drop = FALSE]
  )
  surface_check <- normal_credible_ellipsoid(surface_draws)
  data.frame(
    feature = feature,
    main_parameter_dimension = unname(main_check["dimension"]),
    main_zero_outside_95_ellipsoid_approx = as.logical(main_check["zero_outside"]),
    main_ellipsoid_tail_area_approx = unname(main_check["tail_area"]),
    surface_parameter_dimension = unname(surface_check["dimension"]),
    surface_zero_outside_95_ellipsoid_approx = as.logical(surface_check["zero_outside"]),
    surface_ellipsoid_tail_area_approx = unname(surface_check["tail_area"]),
    stringsAsFactors = FALSE
  )
}))
feature_evidence <- dplyr::left_join(feature_evidence, joint_checks, by = "feature")
dropped_features <- setdiff(EXPECTED_COVARIATES, feature_evidence$feature)
if (length(dropped_features)) {
  feature_evidence <- dplyr::bind_rows(
    feature_evidence,
    data.frame(
      feature = dropped_features,
      number_of_contrasts = 0L,
      any_pointwise_95_interval_excludes_zero = NA,
      any_global_simultaneous_95_interval_excludes_zero = NA,
      strongest_posterior_direction_probability = NA,
      evidence_classification = "not estimable: covariate was absent or dropped from design",
      main_parameter_dimension = 0,
      main_zero_outside_95_ellipsoid_approx = NA,
      main_ellipsoid_tail_area_approx = NA,
      surface_parameter_dimension = 0,
      surface_zero_outside_95_ellipsoid_approx = NA,
      surface_ellipsoid_tail_area_approx = NA,
      stringsAsFactors = FALSE
    )
  )
}
feature_evidence$feature <- factor(
  feature_evidence$feature, levels = EXPECTED_COVARIATES
)
feature_evidence <- dplyr::arrange(feature_evidence, feature)
feature_evidence$feature <- as.character(feature_evidence$feature)
write_csv(feature_evidence, "all_covariate_evidence.csv")

contrast_plot_data <- all_contrasts |>
  dplyr::mutate(
    display = ifelse(feature == "RACE", paste0("RACE: ", contrast), feature),
    display = factor(display, levels = rev(unique(display)))
  )
p_contrasts <- ggplot2::ggplot(contrast_plot_data, ggplot2::aes(y = display)) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey45", linetype = "dashed") +
  ggplot2::geom_segment(
    ggplot2::aes(x = simultaneous_lower_95, xend = simultaneous_upper_95,
                 yend = display),
    linewidth = 0.7, colour = "grey65"
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(x = lower_95, xend = upper_95, yend = display),
    linewidth = 1.8, colour = "#2878B5"
  ) +
  ggplot2::geom_point(ggplot2::aes(x = mean), size = 1.8, colour = "black") +
  ggplot2::labs(x = "Posterior CATE contrast on fitted model scale", y = NULL,
                title = NULL, subtitle = NULL) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank())
contrast_plot_path <- file.path(WORKSPACE_DIR, "all_covariate_contrasts.png")
contrast_plot_height <- max(7, 1.5 + 0.26 * nrow(contrast_plot_data))
ggplot2::ggsave(contrast_plot_path, p_contrasts, width = 10,
                height = contrast_plot_height, dpi = PLOT_DPI, bg = "white")
say("  -> %s\n", contrast_plot_path)

# ---- Conventional randomized-trial interaction sensitivity analysis ---------
# This is not a replacement for the Bayesian model. It asks whether adding one
# ARM-by-covariate interaction to a common main-effects model improves fit, and
# reports both raw, BH, and Bonferroni p-values across the original covariates.
if (!is.null(train_data) && all(c("y", "ARM") %in% names(train_data))) {
  available_covariates <- EXPECTED_COVARIATES[
    EXPECTED_COVARIATES %in% names(train_data)
  ]
  classical_data <- train_data[, c("y", "ARM", available_covariates), drop = FALSE]
  classical_data <- classical_data[stats::complete.cases(classical_data), , drop = FALSE]
  base_formula <- stats::reformulate(c("ARM", available_covariates), response = "y")
  base_fit <- stats::lm(base_formula, data = classical_data)
  classical_tests <- dplyr::bind_rows(lapply(available_covariates, function(feature) {
    full_formula <- stats::reformulate(
      c("ARM", available_covariates, paste0("ARM:", feature)), response = "y"
    )
    full_fit <- try(stats::lm(full_formula, data = classical_data), silent = TRUE)
    if (inherits(full_fit, "try-error")) {
      return(data.frame(feature = feature, numerator_df = NA, F = NA,
                        p_value = NA, model_rank = NA))
    }
    cmp <- stats::anova(base_fit, full_fit)
    data.frame(
      feature = feature,
      numerator_df = cmp$Df[2],
      F = cmp$F[2],
      p_value = cmp$`Pr(>F)`[2],
      model_rank = full_fit$rank,
      stringsAsFactors = FALSE
    )
  })) |>
    dplyr::mutate(
      p_BH = stats::p.adjust(p_value, method = "BH"),
      p_Bonferroni = stats::p.adjust(p_value, method = "bonferroni"),
      significant_raw_0_05 = p_value < 0.05,
      significant_BH_0_05 = p_BH < 0.05,
      significant_Bonferroni_0_05 = p_Bonferroni < 0.05,
      analysis = "ARM interaction added to model with all baseline main effects"
    )
  write_csv(classical_tests, "classical_interaction_sensitivity.csv")

  if ("RACE" %in% names(classical_data)) {
    race_unadjusted <- classical_data |>
      dplyr::group_by(RACE) |>
      dplyr::summarise(
        n = dplyr::n(),
        treated_n = sum(ARM == 1),
        control_n = sum(ARM == 0),
        treated_fraction = mean(ARM),
        treated_mean = mean(y[ARM == 1]),
        control_mean = mean(y[ARM == 0]),
        treatment_difference = treated_mean - control_mean,
        standard_error = sqrt(
          stats::var(y[ARM == 1]) / treated_n +
            stats::var(y[ARM == 0]) / control_n
        ),
        lower_95_normal = treatment_difference - 1.96 * standard_error,
        upper_95_normal = treatment_difference + 1.96 * standard_error,
        .groups = "drop"
      )
    write_csv(race_unadjusted, "race_unadjusted_treatment_effects.csv")
  }
}

race_columns <- grep("^RACE[0-9]+$", colnames(X_final), value = TRUE)
if (length(race_columns)) {
  race_columns <- race_columns[order(as.integer(sub("^RACE", "", race_columns)))]
  k_race <- length(race_columns) + 1L
  race_levels <- recover_factor_levels("RACE", k_race)

  # Under contr.sum(K), each row of this K x (K-1) matrix maps the stored
  # contrast coefficients to the fitted category deviation from the grand mean.
  race_beta_draws <- B_all[, race_columns, drop = FALSE]
  race_coding <- stats::contr.sum(k_race)
  race_level_draws <- race_beta_draws %*% t(race_coding)
  race_level_summary <- posterior_summary(
    race_level_draws, race_levels, simultaneous = TRUE
  ) |>
    dplyr::rename(race_level = term) |>
    dplyr::mutate(
      estimand = "main-effect deviation from grand mean on fitted model scale",
      .after = race_level
    )
  write_csv(race_level_summary, "race_level_main_effect_deviations.csv")

  pairs <- utils::combn(seq_len(k_race), 2L)
  race_pair_draws <- do.call(cbind, lapply(seq_len(ncol(pairs)), function(j) {
    race_level_draws[, pairs[1, j]] - race_level_draws[, pairs[2, j]]
  }))
  race_pair_names <- vapply(seq_len(ncol(pairs)), function(j) {
    paste(race_levels[pairs[1, j]], "minus", race_levels[pairs[2, j]])
  }, character(1))
  pairwise_summary <- posterior_summary(
    race_pair_draws, race_pair_names, simultaneous = TRUE
  ) |>
    dplyr::mutate(
      level_1 = race_levels[pairs[1, ]],
      level_2 = race_levels[pairs[2, ]],
      estimand = "pairwise main-effect contrast on fitted model scale",
      .before = 1
    )
  write_csv(pairwise_summary, "race_pairwise_main_effect_contrasts.csv")
}

# ---- Rank check ---------------------------------------------------------------
rank_x <- qr(X_final)$rank
rank_with_intercept <- qr(cbind(`(Intercept)` = 1, X_final))$rank
condition_number <- tryCatch(kappa(cbind(`(Intercept)` = 1, X_final), exact = TRUE),
                             error = function(e) NA_real_)
design_diagnostics <- data.frame(
  rows = nrow(X_final),
  columns = ncol(X_final),
  rank = rank_x,
  rank_deficit = ncol(X_final) - rank_x,
  full_rank = rank_x == ncol(X_final),
  columns_with_intercept = ncol(X_final) + 1L,
  rank_with_intercept = rank_with_intercept,
  rank_deficit_with_intercept = ncol(X_final) + 1L - rank_with_intercept,
  full_rank_with_intercept = rank_with_intercept == ncol(X_final) + 1L,
  condition_number_with_intercept = condition_number
)
write_csv(design_diagnostics, "design_diagnostics.csv")

notes <- c(
  "HOW TO READ THESE OUTPUTS",
  "",
  "1. shapley_grouped_importance.png is the corrected descriptive ranking.",
  "   All RACE contrast-column contributions are summed within each person",
  "   before absolute values are taken.",
  "   shapley_distribution.png adds the coloured point plot from methodology-7.",
  "   It shows individual Shapley summaries, with colour indicating each",
  "   design column's observed value (blue = low, red = high). RACE contrast",
  "   columns are displayed separately, as in the original distribution plot.",
  "",
  "2. A large Shapley bar is not a declaration of statistical significance.",
  "   The Shapley calculation ranks contribution magnitude after fitting.",
  "",
  "3. all_covariate_evidence.csv is the main one-row-per-covariate evidence table.",
  "   Its primary flag uses simultaneous 95% posterior intervals across every",
  "   fitted CATE contrast in all_covariate_contrasts.csv. Continuous contrasts",
  "   compare the 90th with the 10th percentile; binary contrasts compare levels;",
  "   and RACE is assessed through valid pairwise category contrasts.",
  "",
  "4. Blue intervals in all_covariate_contrasts.png are pointwise 95% intervals;",
  "   wider grey intervals are simultaneous across all displayed contrasts.",
  "",
  "5. posterior_coefficients.csv and posterior_interactions.csv report both",
  "   pointwise and simultaneous intervals for model parameters. Main coefficients",
  "   alone are not complete modifier effects when interactions are present.",
  "   posterior_BLA1C_beta.png and posterior_AGE_beta.png show the main-effect",
  "   coefficient posteriors for baseline HbA1c and AGE, respectively.",
  "   Their corresponding contrast plots show 90th-minus-10th-percentile CATE",
  "   contrasts including interactions. The plots show pointwise intervals; their",
  "   companion CSVs retain the global simultaneous intervals.",
  "",
  "6. RACE1, RACE2, etc. are sum-contrast coordinates, not clinical categories.",
  "   Use race_level_main_effect_deviations.csv and the pairwise-contrast file",
  "   for category-level interpretation. They describe main effects at centered",
  "   values of other covariates and do not include the full interaction surface.",
  "",
  "7. Any RACE-related pattern is exploratory, conditional on recorded RACE,",
  "   model- and coding-sensitive, and is not evidence that race is a causal",
  "   biological treatment-effect modifier. Check counts and external stability.",
  "",
  "8. classical_interaction_sensitivity.csv is a separate conventional sensitivity",
  "   analysis. Raw, BH, and Bonferroni p-values must not be presented as though",
  "   they came from the Bayesian hc-hs model.",
  "",
  "9. Four independently initialized chains are assessed with rank-normalized",
  "   R-hat <= 1.01 and bulk and tail ESS >= 400. Inspect the trace plots and",
  "   Monte Carlo errors alongside these diagnostic flags.",
  "",
  "10. All posterior summaries are on the fitted model scale. Signs and zero",
  "    crossing are unaffected by positive rescaling."
)
notes_path <- file.path(WORKSPACE_DIR, "diagnostic_notes.txt")
writeLines(notes, notes_path)
say("  -> %s\n", notes_path)

feature_list <- function(x) {
  x <- as.character(x[!is.na(x) & nzchar(x)])
  if (length(x)) paste(x, collapse = ", ") else "None"
}
simultaneous_features <- feature_evidence$feature[
  !is.na(feature_evidence$any_global_simultaneous_95_interval_excludes_zero) &
    feature_evidence$any_global_simultaneous_95_interval_excludes_zero
]
pointwise_only_features <- feature_evidence$feature[
  !is.na(feature_evidence$any_pointwise_95_interval_excludes_zero) &
    feature_evidence$any_pointwise_95_interval_excludes_zero &
    !feature_evidence$any_global_simultaneous_95_interval_excludes_zero
]
no_evidence_features <- feature_evidence$feature[
  !is.na(feature_evidence$any_pointwise_95_interval_excludes_zero) &
    !feature_evidence$any_pointwise_95_interval_excludes_zero
]
not_estimable_features <- feature_evidence$feature[
  is.na(feature_evidence$any_pointwise_95_interval_excludes_zero)
]

automatic_summary <- c(
  "AUTOMATIC INTERPRETATION SUMMARY",
  "",
  "Primary Bayesian contrast screen:",
  paste0("  Evidence after the global simultaneous 95% adjustment: ",
         feature_list(simultaneous_features)),
  paste0("  Pointwise 95% evidence only: ", feature_list(pointwise_only_features)),
  paste0("  No pointwise 95% contrast evidence: ", feature_list(no_evidence_features)),
  paste0("  Not estimable / dropped: ", feature_list(not_estimable_features)),
  "",
  paste0("RACE classification: ",
         feature_evidence$evidence_classification[feature_evidence$feature == "RACE"]),
  "",
  sprintf("Design full rank with intercept: %s; condition number: %.3g",
          design_diagnostics$full_rank_with_intercept,
          design_diagnostics$condition_number_with_intercept),
  sprintf("Main-coefficient MCMC diagnostic flags: %d of %d",
          sum(main_mcmc$diagnostic_flag, na.rm = TRUE), nrow(main_mcmc)),
  if (all(!main_mcmc$rhat_assessable))
    "R-hat warning: unavailable because the saved fit has only one independent chain."
  else
    sprintf("Largest main-coefficient rank-normalized R-hat: %.4f",
            max(main_mcmc$rhat, na.rm = TRUE)),
  ""
)
if (exists("classical_tests")) {
  bh_features <- classical_tests$feature[
    !is.na(classical_tests$significant_BH_0_05) & classical_tests$significant_BH_0_05
  ]
  automatic_summary <- c(
    automatic_summary,
    "Separate conventional interaction sensitivity analysis:",
    paste0("  Covariates with BH-adjusted p < 0.05: ", feature_list(bh_features))
  )
  race_classical <- classical_tests[classical_tests$feature == "RACE", , drop = FALSE]
  if (nrow(race_classical)) {
    automatic_summary <- c(
      automatic_summary,
      sprintf("  RACE interaction: raw p = %.4g; BH p = %.4g; Bonferroni p = %.4g",
              race_classical$p_value, race_classical$p_BH,
              race_classical$p_Bonferroni)
    )
  }
  automatic_summary <- c(automatic_summary, "")
}
automatic_summary <- c(
  automatic_summary,
  "This is an automated diagnostic screen, not manuscript-ready causal language.",
  "Read diagnostic_notes.txt before interpreting RACE or the word 'evidence'."
)
summary_path <- file.path(WORKSPACE_DIR, "interpretation_summary.txt")
writeLines(automatic_summary, summary_path)
say("  -> %s\n", summary_path)

say("\nDone. No model was refitted and no input file was modified.\n")
