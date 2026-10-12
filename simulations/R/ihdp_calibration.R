# Reproduce the manuscript's OOS, 95% base-interval, multiplicative scale,
# leave-one-replicate-out calibration for the two Bayesian manuscript models.
# Read one fitted fold at a time; never combine its draw indices with another fold.

ihdp_calibration_stats <- function(draws, truth, base_level = .95) {
  if (!is.matrix(draws) || nrow(draws) != length(truth) || ncol(draws) < 2L ||
      any(!is.finite(draws)) || any(!is.finite(truth))) stop("Invalid IHDP calibration posterior")
  a <- (1 - base_level) / 2
  bounds <- t(vapply(seq_len(nrow(draws)), function(i)
    as.numeric(stats::quantile(draws[i, ], c(a, 1 - a), type = 8, names = FALSE)), numeric(2)))
  centre <- rowMeans(draws)
  # The historical scale construction assumes intervals expand monotonically
  # about their centre. Fail explicitly when an extreme skew violates that.
  if (any(centre < bounds[, 1] | centre > bounds[, 2]))
    stop("A posterior mean is outside its base interval; monotone interval scaling is not applicable")
  distance <- abs(truth - centre)
  radius <- ifelse(truth >= centre, bounds[, 2] - centre, centre - bounds[, 1])
  required <- ifelse(distance == 0, 0, ifelse(radius > 0, distance / radius, Inf))
  data.frame(centre = centre, base_lo = bounds[, 1], base_hi = bounds[, 2],
             truth = truth, required_scale = required)
}

ihdp_interval_metrics <- function(lower, upper, centre, truth) {
  if (any(!is.finite(c(lower, upper, centre, truth))) || any(lower > upper))
    stop("Invalid calibrated IHDP interval")
  direction <- ifelse(lower > 0, 1L, ifelse(upper < 0, -1L, 0L))
  excluded <- direction != 0L
  correct <- excluded & direction == sign(truth)
  fraction <- function(a, b) if (b > 0) a / b else NA_real_
  data.frame(n_used = length(truth), coverage = mean(lower <= truth & upper >= truth),
    power_excl0 = mean(excluded), power_correct_excl0 = mean(correct),
    sign_correct = fraction(sum(correct), sum(excluded)),
    type_s = fraction(sum(excluded & direction != sign(truth)), sum(excluded)),
    sens_benefit = fraction(sum(direction == 1L & truth > 0), sum(truth > 0)),
    sens_harm = fraction(sum(direction == -1L & truth < 0), sum(truth < 0)),
    mean_width = mean(upper - lower), median_width = median(upper - lower),
    bias = mean(centre - truth), mean_abs_error = mean(abs(centre - truth)))
}

ihdp_empirical_scale <- function(x, target = .95) {
  if (!length(x) || anyNA(x) || any(x < 0) || target <= 0 || target >= 1)
    stop("Invalid empirical calibration scales or target")
  value <- sort(x)[max(1L, ceiling(target * length(x)))]
  if (!is.finite(value)) stop("Target IHDP coverage is unattainable by finite interval scaling")
  value
}

ihdp_loo_calibrate <- function(statistics, replicates, target = .95, n_expected = 747L) {
  models <- c("SemiParam", "BCF")
  if (length(replicates) < 2L || anyDuplicated(replicates))
    stop("LOO IHDP calibration requires at least two distinct replicates")
  keys <- paste(statistics$replicate, statistics$model, statistics$index, sep = ":")
  if (anyDuplicated(keys) || any(!statistics$model %in% models) ||
      any(!statistics$replicate %in% replicates)) stop("Duplicate or unexpected calibration rows")
  for (replicate in replicates) for (model in models) {
    d <- statistics[statistics$replicate == replicate & statistics$model == model, ]
    if (nrow(d) != n_expected || !setequal(d$index, seq_len(n_expected)))
      stop("Incomplete paired IHDP calibration inputs: ", replicate, "/", model)
  }
  parameters <- list(); rows <- list(); global <- list()
  for (model in models) {
    d <- statistics[statistics$model == model, ]
    pooled <- ihdp_empirical_scale(d$required_scale, target)
    global[[model]] <- data.frame(model = model, calibration_value = pooled,
      apparent_coverage = mean(d$required_scale <= pooled), n_calibration = nrow(d),
      n_unattainable = sum(is.infinite(d$required_scale)))
    for (replicate in replicates) {
      training <- d[d$replicate != replicate, ]
      evaluation <- d[d$replicate == replicate, ]
      value <- ihdp_empirical_scale(training$required_scale, target)
      parameters[[length(parameters) + 1L]] <- data.frame(replicate = replicate, model = model,
        calibration_value = value, n_calibration = nrow(training),
        n_unattainable = sum(is.infinite(training$required_scale)),
        calibration_excludes_replicate = TRUE)
      lower <- evaluation$centre + value * (evaluation$base_lo - evaluation$centre)
      upper <- evaluation$centre + value * (evaluation$base_hi - evaluation$centre)
      row <- ihdp_interval_metrics(lower, upper, evaluation$centre, evaluation$truth)
      row$replicate <- replicate; row$model <- model; row$calibration_value <- value
      row$target <- target; row$sample <- "OOS"; row$method <- "scale"; row$calibration <- "loo"
      rows[[length(rows) + 1L]] <- row
    }
  }
  list(results = do.call(rbind, rows), parameters = do.call(rbind, parameters),
       global_parameters = do.call(rbind, global))
}

ihdp_mean_se <- function(x) {
  x <- x[is.finite(x)]
  c(mean = if (length(x)) mean(x) else NA_real_,
    se = if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_)
}

ihdp_summary <- function(data, groups, metrics) {
  keys <- do.call(paste, c(data[groups], sep = ":"))
  pieces <- lapply(split(data, keys), function(d) {
    row <- d[1, groups, drop = FALSE]; row$n_replicates <- nrow(d)
    for (metric in metrics) {
      value <- ihdp_mean_se(d[[metric]])
      row[[metric]] <- unname(value["mean"])
      row[[paste0(metric, "_se")]] <- unname(value["se"])
      row[[paste0(metric, "_n")]] <- sum(is.finite(d[[metric]]))
    }
    row
  })
  do.call(rbind, pieces)
}

ihdp_calibration_pairs <- function(data) {
  metrics <- c("coverage", "power_excl0", "power_correct_excl0", "sign_correct", "type_s",
               "sens_benefit", "sens_harm", "mean_width", "bias")
  a <- data[data$model == "SemiParam", c("replicate", metrics)]
  b <- data[data$model == "BCF", c("replicate", metrics)]
  paired <- merge(a, b, by = "replicate", suffixes = c("", ".ref"))
  if (nrow(paired) != nrow(a) || nrow(a) != nrow(b)) stop("Incomplete calibration pairs")
  do.call(rbind, lapply(metrics, function(metric) {
    delta <- paired[[metric]] - paired[[paste0(metric, ".ref")]]
    delta <- delta[is.finite(delta)]; n <- length(delta)
    value <- ihdp_mean_se(delta)
    critical <- if (n > 1L) stats::qt(.975, n - 1L) else NA_real_
    data.frame(model = "SemiParam", reference = "BCF", metric = metric, n_replicates = n,
      mean_difference = unname(value["mean"]), se = unname(value["se"]),
      ci_low = unname(value["mean"] - critical * value["se"]),
      ci_high = unname(value["mean"] + critical * value["se"]))
  }))
}

ihdp_read_calibration <- function(cfg) {
  jobs <- get_jobs("ihdp", cfg)
  jobs <- jobs[jobs$model %in% c("SemiParam", "BCF") & jobs$fold > 0L, , drop = FALSE]
  parts <- vector("list", nrow(jobs))
  for (k in seq_len(nrow(jobs))) {
    job <- jobs[k, , drop = FALSE]
    path <- result_path(job, cfg)
    if (!file.exists(path)) stop("Missing calibration result: ", path)
    result <- readRDS(path)
    if (!identical(result$job, job)) stop("Mismatched IHDP posterior job: ", path)
    draws <- result$tau_draws
    if (is.null(draws)) stop("Missing retained IHDP draws in ", path,
      "; keep_ihdp_posteriors=TRUE is required for interval calibration")
    individual <- result$individual
    value <- ihdp_calibration_stats(draws, individual$truth)
    if (max(abs(value$centre - individual$tau_mean)) > 1e-9 * max(1, max(abs(value$centre))))
      stop("IHDP posterior draws do not match saved summaries: ", path)
    value$replicate <- job$replicate; value$model <- job$model
    value$index <- individual$index; value$fold <- job$fold
    parts[[k]] <- value
  }
  do.call(rbind, parts)
}

ihdp_native_posterior_tables <- function(metrics, signal_groups) {
  posterior_metrics <- c("cate_coverage", "cate_width", "exclude_zero",
    "sign_correct", "sign_error_given_exclusion", "auc_sign")
  columns <- c("replicate", "model", "sample", "n", "sd_tau_true", posterior_metrics)
  if (!is.data.frame(metrics) || !all(columns %in% names(metrics)) ||
      !is.data.frame(signal_groups) ||
      !all(c("replicate", "sd_tau_true", "signal_quartile") %in% names(signal_groups)) ||
      !nrow(signal_groups) || anyDuplicated(signal_groups$replicate) ||
      anyNA(signal_groups$replicate) || any(!is.finite(signal_groups$sd_tau_true)) ||
      any(signal_groups$sd_tau_true < 0))
    stop("Incomplete IHDP native posterior table inputs")
  signal_groups <- signal_groups[order(signal_groups$sd_tau_true, signal_groups$replicate), ]
  expected_quartiles <- pmin(4L, ceiling(seq_len(nrow(signal_groups)) * 4 / nrow(signal_groups)))
  if (!identical(as.numeric(signal_groups$signal_quartile), as.numeric(expected_quartiles)))
    stop("IHDP native posterior quartiles differ from the paired point-performance groups")
  native <- metrics[metrics$model %in% c("SemiParam", "BCF"), columns, drop = FALSE]
  expected <- expand.grid(replicate = signal_groups$replicate,
    model = c("SemiParam", "BCF"), sample = c("IS", "OOS"))
  key <- function(x) paste(x$replicate, x$model, x$sample, sep = ":")
  if (anyDuplicated(key(native)) || !setequal(key(native), key(expected)))
    stop("Incomplete paired IHDP native posterior grid")
  if (any(!is.finite(native$n)) || any(native$n != 747L) ||
      any(!is.finite(native$sd_tau_true)) ||
      any(!is.finite(native$cate_width)) || any(native$cate_width < 0))
    stop("Invalid IHDP native posterior sizes, signal or interval widths")
  for (name in c("cate_coverage", "exclude_zero"))
    if (any(!is.finite(native[[name]])) || any(native[[name]] < 0 | native[[name]] > 1))
      stop("Invalid IHDP native posterior probability: ", name)
  for (name in c("sign_correct", "sign_error_given_exclusion", "auc_sign")) {
    values <- native[[name]]
    if (any(is.infinite(values)) || any(is.nan(values)) ||
        any(values < 0 | values > 1, na.rm = TRUE))
      stop("Invalid IHDP native posterior probability: ", name)
  }
  no_exclusion <- native$exclude_zero == 0
  if (any(is.na(native$sign_correct) != no_exclusion) ||
      any(is.na(native$sign_error_given_exclusion) != no_exclusion) ||
      any(abs(native$sign_correct + native$sign_error_given_exclusion - 1) > 1e-12,
          na.rm = TRUE))
    stop("IHDP native conditional sign metrics have inconsistent denominators")
  signal_index <- match(native$replicate, signal_groups$replicate)
  expected_signal <- signal_groups$sd_tau_true[signal_index]
  if (any(abs(native$sd_tau_true - expected_signal) > 1e-10 * pmax(1, abs(expected_signal))))
    stop("IHDP native models do not share the paired truth-signal groups")
  native <- native[order(native$replicate, native$model, native$sample), , drop = FALSE]
  summary <- ihdp_summary(native, c("model", "sample"), posterior_metrics)
  oos <- native[native$sample == "OOS", , drop = FALSE]
  oos$signal_quartile <- signal_groups$signal_quartile[match(oos$replicate, signal_groups$replicate)]
  bands <- ihdp_summary(oos, c("signal_quartile", "model", "sample"),
    c("sd_tau_true", posterior_metrics))
  list(by_replicate = native, summary = summary, signal_by_replicate = oos, signal_bands = bands)
}

export_ihdp_analysis <- function(combined, out, cfg) {
  metrics <- if (is.data.frame(combined)) combined else
    do.call(rbind, lapply(combined, function(x) x$metrics))
  replicates <- as.integer(cfg$ihdp_replicates)
  models <- c("SemiParam", "BCF", "DR", "R", "S")
  expected <- expand.grid(replicate = replicates, model = models, sample = c("IS", "OOS"))
  key <- function(x) paste(x$replicate, x$model, x$sample, sep = ":")
  if (anyDuplicated(key(metrics)) || !setequal(key(metrics), key(expected)))
    stop("IHDP table export requires the complete five-model replicate grid")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(metrics, file.path(out, "ihdp_per_replicate.csv"), row.names = FALSE)
  numerical <- setdiff(names(metrics)[vapply(metrics, is.numeric, logical(1))], c("replicate", "n"))
  summary <- ihdp_summary(metrics, c("model", "sample"), numerical)
  utils::write.csv(summary, file.path(out, "ihdp_native_summary.csv"), row.names = FALSE)
  oos <- metrics[metrics$sample == "OOS", ]
  ranks <- lapply(split(oos, oos$replicate), function(d) {
    if (any(!is.finite(d$cate_rmse))) stop("Cannot rank incomplete IHDP methods")
    d$rank <- rank(d$cate_rmse, ties.method = "average"); d
  })
  ranked <- do.call(rbind, ranks)
  points <- ihdp_summary(ranked, "model", c("cate_rmse", "cate_spearman", "calibration_slope", "rank"))
  utils::write.csv(points, file.path(out, "ihdp_point_performance.csv"), row.names = FALSE)

  sp <- oos[oos$model == "SemiParam", c("replicate", "sd_tau_true", "cate_rmse")]
  bcf <- oos[oos$model == "BCF", c("replicate", "cate_rmse")]
  pairs <- merge(sp, bcf, by = "replicate", suffixes = c("_sp", "_bcf"))
  pairs <- pairs[order(pairs$sd_tau_true, pairs$replicate), ]
  # Exactly 25 per quartile for the 100-replicate protocol; deterministic tie order.
  pairs$signal_quartile <- pmin(4L, ceiling(seq_len(nrow(pairs)) * 4 / nrow(pairs)))
  pairs$delta_rmse <- pairs$cate_rmse_sp - pairs$cate_rmse_bcf
  pairs$semiparam_win <- as.integer(pairs$delta_rmse < 0)
  band <- ihdp_summary(pairs, "signal_quartile", c("sd_tau_true", "delta_rmse", "semiparam_win"))
  band$semiparam_wins <- vapply(band$signal_quartile,
    function(q) sum(pairs$semiparam_win[pairs$signal_quartile == q]), integer(1))
  utils::write.csv(pairs, file.path(out, "ihdp_signal_by_replicate.csv"), row.names = FALSE)
  utils::write.csv(band, file.path(out, "ihdp_signal_bands.csv"), row.names = FALSE)

  # Native intervals and calibrated intervals are separate outputs. Reuse the
  # exact point-performance quartile assignment for both Bayesian models.
  posterior <- ihdp_native_posterior_tables(metrics, pairs)
  for (name in names(posterior)) utils::write.csv(posterior[[name]],
    file.path(out, paste0("ihdp_posterior_", name, ".csv")), row.names = FALSE)

  if (length(replicates) < 2L) stop("At least two IHDP replicates are needed for LOO calibration")
  statistics <- ihdp_read_calibration(cfg)
  calibrated <- ihdp_loo_calibrate(statistics, replicates)
  calibration_metrics <- c("coverage", "power_excl0", "power_correct_excl0", "sign_correct", "type_s",
    "sens_benefit", "sens_harm", "mean_width", "median_width", "bias", "mean_abs_error", "calibration_value")
  calibrated_summary <- ihdp_summary(calibrated$results, "model", calibration_metrics)
  calibrated_summary$target <- .95
  calibrated_summary$coverage_gap <- calibrated_summary$coverage - .95
  calibration_path <- file.path(out, "calibration_oos_scale_loo")
  dir.create(calibration_path, recursive = TRUE, showWarnings = FALSE)
  files <- list(calibration_parameters = calibrated$parameters,
    global_calibration_parameters = calibrated$global_parameters,
    matched_coverage_by_replicate = calibrated$results,
    matched_coverage_summary = calibrated_summary,
    matched_coverage_pairwise = ihdp_calibration_pairs(calibrated$results))
  for (name in names(files)) utils::write.csv(files[[name]],
    file.path(calibration_path, paste0(name, ".csv")), row.names = FALSE)
  writeLines(c("sample=OOS", "method=scale", "base_level=0.95", "target=0.95",
    "calibration=leave-one-replicate-out", "models=SemiParam,BCF", "quantile_type=8",
    "widen_only=false", "complete_replicates=true",
    paste0("replicates=", paste(replicates, collapse = ",")),
    "Every requested replicate is used; convergence failures are flagged by the parent analysis.",
    "Intervals calibrated to empirical coverage are not unchanged Bayesian credible intervals."),
    file.path(calibration_path, "calibration_config.txt"))
  invisible(list(native = summary, point = points, signal = band,
    posterior = posterior, calibration = calibrated_summary))
}
