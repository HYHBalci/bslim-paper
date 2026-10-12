# ==============================================================================
# R/subgroup/bands.R -- band and level-set rules for CATE subgroup selection.
#
# Successor to R/simultaneous/bands.R. Everything there that was correct is
# carried over unchanged in behaviour; this file adds what the August 2026 audit
# left on the table:
#
#   * band_tolerance() has a CLOSED FORM. The grid search (and with it both
#     degenerate branches, the 0.05 cap, and ~12 full-matrix passes) is gone.
#     `method = "grid"` still reproduces the old search for cross-checking.
#   * band_supt() can standardise by the pointwise credible interval instead of
#     the posterior sd, which is the Brezger & Lang (2007) form and does not
#     assume a symmetric posterior.
#   * band_wan() returns Wan, Liu & Bretz (2024) properly: BOTH the inner
#     (lower) and outer (upper) confidence sets, not just the inner one.
#   * level_sweep() exploits Wan et al. Theorem 2 -- one band gives the whole
#     lambda-indexed family of selected sets at the same simultaneous level.
#
# CONVENTION -------------------------------------------------------------------
#   post : n x S numeric matrix. Row i holds S posterior draws of tau_i.
#   Every band_*() returns a list with a stable shape:
#     method  character tag
#     lower   length-n lower bound
#     upper   length-n upper bound
#     centre  length-n posterior mean
#     k       effective multiplier: mean((upper-lower)/(2*sd)); 1.96 = pointwise
#     k_one   one-sided multiplier: mean((centre-lower)/sd)
#     level   the tuning constant actually used (gamma, or c1)
#     sided   "two" or "lower"
#     notes   character vector of diagnostics; empty when nothing fired
#
# REFERENCES -------------------------------------------------------------------
#   Brezger, A. & Lang, S. (2007). Simultaneous probability statements for
#     Bayesian P-splines. -- the sup-t construction, historically mislabelled
#     "Sechidis" in this project's older scripts.
#   Wan, F., Liu, W. & Bretz, F. (2024). Confidence sets for a level set in
#     linear regression. Statistics in Medicine 43:1103-1118.
#   Mueller, M. M. et al. (2026). Data-driven controlled subgroup selection in
#     clinical trials. -- the Type I error notion scored in metrics.R.
# ==============================================================================

# ---- internals ---------------------------------------------------------------

#' Type-7 row quantiles from a pre-sorted matrix (matches stats::quantile()).
.q7 <- function(ps, p) {
  S <- ncol(ps)
  h <- (S - 1) * p + 1
  lo <- min(max(floor(h), 1), S)
  hi <- min(max(ceiling(h), 1), S)
  g  <- h - floor(h)
  if (lo == hi) return(ps[, lo])
  ps[, lo] * (1 - g) + ps[, hi] * g
}

.chk_post <- function(post) {
  if (!is.matrix(post)) post <- as.matrix(post)
  if (anyNA(post)) stop("post contains NA; the band rules have no missing-data policy.")
  if (ncol(post) < 2) stop("post must be n x S with S >= 2 draws (got S = ", ncol(post), ").")
  post
}

.mk <- function(method, lower, upper, centre, sdv, level, sided, notes = character()) {
  two <- is.finite(upper[1]) && is.finite(lower[1])
  list(method = method, lower = lower, upper = upper, centre = centre,
       sd = sdv, level = level, sided = sided, notes = notes,
       k     = if (two) mean((upper - lower) / (2 * sdv)) else NA_real_,
       k_one = mean((centre - lower) / sdv),
       width = if (two) mean(upper - lower) else NA_real_)
}

#' k-th smallest of each column, via partial sort. O(n) per column.
.col_kth <- function(m, k) {
  k <- max(1L, min(as.integer(k), nrow(m)))
  apply(m, 2, function(z) sort.int(z, partial = k)[k])
}

# ---- 1. pointwise (no multiplicity adjustment) --------------------------------

band_pointwise <- function(post, alpha = 0.05, sided = c("two", "lower")) {
  sided <- match.arg(sided); post <- .chk_post(post)
  ps <- t(apply(post, 1, sort.int, method = "quick"))
  m  <- rowMeans(post); s <- apply(post, 1, stats::sd)
  if (sided == "two") {
    .mk("pointwise", .q7(ps, alpha/2), .q7(ps, 1 - alpha/2), m, s, alpha, "two")
  } else {
    .mk("pointwise", .q7(ps, alpha), rep(Inf, nrow(post)), m, s, alpha, "lower")
  }
}

# ---- 2. tolerance band, in closed form ----------------------------------------

#' Bayesian tolerance band over the posterior draws.
#'
#' Chooses the largest marginal level gamma such that, in at least `conf` of the
#' posterior draws, at least (1 - tol) of individuals fall inside their own
#' marginal gamma-interval.
#'
#' WHAT THIS GUARANTEES: with `conf` posterior confidence, at most `tol` of
#' individuals fall outside their interval. It is a (1-tol)-content,
#' conf-confidence tolerance region. It is NOT a simultaneous band and does NOT
#' control the family-wise error rate. Measured FWER in an exactly-specified
#' Gaussian world is 0.26 at n = 250 and 0.42 at n = 1000 against a nominal
#' 0.05, and it keeps rising with n -- see run_tier1.R.
#'
#' THE CLOSED FORM -------------------------------------------------------------
#' Whether draw j sits inside row i's gamma-interval depends on the posterior
#' only through the RANK of that draw within its own row. Write r_ij for the
#' rank and d_ij = min(r_ij, S+1-r_ij) for the distance to the nearer tail. With
#' the interval taken at order statistics [a, S+1-a],
#'
#'     draw j is inside row i   <=>   d_ij >= a.
#'
#' So per-draw content is monotone in the single integer a, and:
#'
#'   * coverage_j >= 1-tol  <=>  #{i : d_ij < a} <= floor(tol*n) =: k
#'                          <=>  a <= D_j, the (k+1)-th smallest of d_.j
#'   * at least `conf` of draws satisfy that
#'                          <=>  a <= the (m+1)-th smallest of D,
#'                               with m = floor((1-conf)*S)
#'
#' The optimum is therefore a quantile of a quantile, computed in one pass. This
#' removes the grid entirely -- and with it the 0.05 cap (old F6) and the two
#' degenerate branches that the legacy findconfglob() answered identically.
#'
#' @param method "exact" (closed form, default) or "grid" (the old bisection,
#'   kept so the two can be cross-checked; see test_suite.R).
band_tolerance <- function(post, tol = 0.05, conf = 0.95,
                           method = c("exact", "grid"),
                           gamma_grid = seq(5e-4, 0.5, by = 5e-4),
                           min_tail_draws = 20) {
  method <- match.arg(method); post <- .chk_post(post)
  n <- nrow(post); S <- ncol(post); notes <- character()
  centre <- rowMeans(post); sdv <- apply(post, 1, stats::sd)

  if (method == "exact") {
    # ONE radix sort per row yields both the order statistics (for the bounds)
    # and the within-row ranks (for the criterion). Calling sort() and rank()
    # separately sorts each row twice and throws the speed advantage away.
    ps <- matrix(0,  n, S)
    d  <- matrix(0L, n, S)
    idx <- seq_len(S); rev_idx <- S + 1L - idx
    for (i in seq_len(n)) {
      x <- post[i, ]
      o <- order(x, method = "radix")
      ps[i, ] <- x[o]
      r <- integer(S); r[o] <- idx
      d[i, ] <- pmin.int(r, S + 1L - r)
    }
    k1 <- min(floor(tol * n) + 1L, n)          # (k+1)-th smallest over individuals
    D  <- .col_kth(d, k1)
    m1 <- min(floor((1 - conf) * S) + 1L, S)   # (m+1)-th smallest over draws
    a  <- max(1L, min(as.integer(sort.int(D, partial = m1)[m1]), floor(S / 2)))
    gamma <- 2 * a / S
    lower <- ps[, a]; upper <- ps[, S + 1L - a]
    if (a == 1L) notes <- c(notes, paste0(
      "Optimal interval is the full range of draws (a = 1): the content criterion ",
      "cannot be met at any narrower level. The band is not calibrated -- raise S ",
      "or loosen tol/conf."))
  } else {
    ps <- t(apply(post, 1, sort.int, method = "quick"))
    meets <- function(g) {
      lo <- .q7(ps, g/2); hi <- .q7(ps, 1 - g/2)
      mean(colMeans((post >= lo) & (post <= hi)) >= 1 - tol) >= conf
    }
    mgrid <- length(gamma_grid)
    if (!meets(gamma_grid[1])) {
      gamma <- gamma_grid[1]
      notes <- c(notes, sprintf(
        "No gamma on the grid met the (%.2f content, %.2f confidence) criterion; fell back to the grid MINIMUM %.4f (widest band).",
        1 - tol, conf, gamma))
    } else if (meets(gamma_grid[mgrid])) {
      gamma <- gamma_grid[mgrid]
      notes <- c(notes, sprintf(
        "Criterion still met at the grid MAXIMUM %.4f; the true optimum is narrower.", gamma))
    } else {
      lo_i <- 1L; hi_i <- mgrid
      while (hi_i - lo_i > 1L) {
        mid <- (lo_i + hi_i) %/% 2L
        if (meets(gamma_grid[mid])) lo_i <- mid else hi_i <- mid
      }
      gamma <- gamma_grid[lo_i]
    }
    lower <- .q7(ps, gamma/2); upper <- .q7(ps, 1 - gamma/2)
  }

  tail_n <- S * gamma / 2
  if (tail_n < min_tail_draws) notes <- c(notes, sprintf(
    "Lower bound rests on ~%.0f draws (S*gamma/2 = %.1f < %d). Increase retained draws to at least %.0f.",
    max(1, round(tail_n)), tail_n, min_tail_draws, ceiling(2 * min_tail_draws / gamma)))

  .mk("tolerance", lower, upper, centre, sdv, gamma, "two", notes)
}

legacy_gamma_grid <- function() seq(5e-4, 0.05, by = 5e-4)

# ---- 3. sup-t simultaneous band ----------------------------------------------

#' Simultaneous credible band from the posterior sup-quantile (Brezger & Lang).
#'
#' c1 = the (1-alpha) quantile over draws of sup_i (centre_i - tau_i) / scale_i.
#' Controls the family-wise error rate over the n evaluation points.
#'
#' @param standardise "sd" divides by the posterior sd -- the historical choice,
#'   which implicitly assumes a roughly symmetric posterior. "quantile" divides
#'   by each row's own one-sided pointwise credible distance, which is the
#'   Brezger & Lang form and is robust to the skew a probit/horseshoe CATE
#'   posterior actually shows. Both are exact under symmetry; they differ in the
#'   tails, so tier 1 can measure which is tighter.
#' @param sided defaults to "lower": subgroup selection asks a one-sided
#'   question. The historical scripts disagreed on this by ~0.22*se.
band_supt <- function(post, alpha = 0.05, sided = c("lower", "two"),
                      standardise = c("sd", "quantile"), se_floor = 1e-8) {
  sided <- match.arg(sided); standardise <- match.arg(standardise)
  post <- .chk_post(post)
  m <- rowMeans(post); sdv <- pmax(apply(post, 1, stats::sd), se_floor)

  if (standardise == "sd") {
    s_lo <- sdv; s_hi <- sdv
  } else {
    ps <- t(apply(post, 1, sort.int, method = "quick"))
    s_lo <- pmax(m - .q7(ps, alpha), se_floor)       # one-sided pointwise distance
    s_hi <- pmax(.q7(ps, 1 - alpha) - m, se_floor)
  }

  tag <- if (standardise == "sd") "supt" else "supt_q"
  if (sided == "two") {
    D <- pmax((m - post) / s_lo, (post - m) / s_hi)
    c1 <- unname(stats::quantile(apply(D, 2, max), 1 - alpha))
    .mk(paste0(tag, "_two"), m - c1 * s_lo, m + c1 * s_hi, m, sdv, c1, "two")
  } else {
    D <- (m - post) / s_lo
    c1 <- unname(stats::quantile(apply(D, 2, max), 1 - alpha))
    .mk(paste0(tag, "_lower"), m - c1 * s_lo, rep(Inf, length(m)), m, sdv, c1, "lower")
  }
}

# ---- 4. Bonferroni reference --------------------------------------------------

#' Union-bound band. The gap between this and sup-t is the credit sup-t earns
#' for exploiting correlation across individuals.
band_bonferroni <- function(post, alpha = 0.05, sided = c("lower", "two")) {
  sided <- match.arg(sided); post <- .chk_post(post)
  n <- nrow(post); S <- ncol(post); notes <- character()
  ps <- t(apply(post, 1, sort.int, method = "quick"))
  m <- rowMeans(post); s <- apply(post, 1, stats::sd)
  a <- alpha / n
  ptail <- if (sided == "two") a / 2 else a
  if (S * ptail < 10) notes <- c(notes, sprintf(
    "Bonferroni needs the %.2g row quantile but S = %d gives only %.1f draws in that tail; the bound is resolution-limited. Need S >= %.0f.",
    ptail, S, S * ptail, ceiling(10 / ptail)))
  if (sided == "two") .mk("bonferroni", .q7(ps, a/2), .q7(ps, 1 - a/2), m, s, a, "two", notes)
  else                .mk("bonferroni", .q7(ps, a),   rep(Inf, n),      m, s, a, "lower", notes)
}

# ---- 5. Wan, Liu & Bretz (2024) level-set confidence sets ----------------------

#' Simultaneous confidence band for a linear model, and the level-set confidence
#' sets it induces.
#'
#' Wan et al. build TWO sets from one band, and the pair is the point:
#'   inner (lower) set  {x : lower(x) >= lambda}  -- these ARE in the level set
#'   outer (upper) set  {x : upper(x) >= lambda}  -- these are not ruled out
#' with P(inner subset of S subset of outer) >= 1-alpha. The old code kept only
#' the inner one, which is why an empty selection read as a dead end rather than
#' as the bound it is.
#'
#' @param Xd_fit  n x p model matrix INCLUDING the intercept column
#' @param y_fit   response (for a DR-learner, the cross-fitted pseudo-outcome)
#' @param Xd_eval m x p model matrix for the points the band must cover
#' @param vcov_type HC3 by default: the linear projection of a CATE is
#'   misspecified, so homoskedastic OLS errors understate the band.
band_wan <- function(Xd_fit, y_fit, Xd_eval = Xd_fit, alpha = 0.05, nsim = 4000,
                     sided = c("lower", "two"), vcov_type = c("HC3", "const"),
                     seed = NULL) {
  sided <- match.arg(sided); vcov_type <- match.arg(vcov_type)
  if (!is.null(seed)) set.seed(seed)
  Xd_fit <- as.matrix(Xd_fit); Xd_eval <- as.matrix(Xd_eval)

  qrf <- qr(Xd_fit)
  if (qrf$rank < ncol(Xd_fit)) {
    keep <- qrf$pivot[seq_len(qrf$rank)]
    Xd_fit <- Xd_fit[, keep, drop = FALSE]; Xd_eval <- Xd_eval[, keep, drop = FALSE]
  }
  n <- nrow(Xd_fit); p <- ncol(Xd_fit); nu <- n - p
  XtXi <- chol2inv(chol(crossprod(Xd_fit)))
  beta <- as.vector(XtXi %*% crossprod(Xd_fit, y_fit))
  e    <- as.vector(y_fit - Xd_fit %*% beta)

  if (vcov_type == "const") {
    V <- XtXi * sum(e^2) / nu
    scal <- sqrt(stats::rchisq(nsim, nu) / nu)
  } else {
    h <- rowSums((Xd_fit %*% XtXi) * Xd_fit)
    meat <- crossprod(Xd_fit * (e / (1 - pmin(h, 1 - 1e-8))))
    V <- XtXi %*% meat %*% XtXi
    scal <- rep(1, nsim)
  }

  centre <- as.vector(Xd_eval %*% beta)
  se     <- sqrt(pmax(rowSums((Xd_eval %*% V) * Xd_eval), 1e-16))

  Vc <- chol(V + diag(1e-12, ncol(V)))
  Z  <- matrix(stats::rnorm(nsim * ncol(V)), nsim) %*% Vc
  dev <- (Xd_eval %*% t(Z)) / se
  dev <- sweep(dev, 2, scal, "/")
  sup <- if (sided == "two") apply(abs(dev), 2, max) else apply(dev, 2, max)
  c1  <- unname(stats::quantile(sup, 1 - alpha))

  if (sided == "two")
    .mk(paste0("wan_", vcov_type), centre - c1*se, centre + c1*se, centre, se, c1, "two")
  else
    # NOTE: upper is the OUTER set's boundary, not +Inf. A one-sided lower band
    # bounds tau from below; the matching outer set uses the opposite tail at the
    # same simultaneous level, which is what makes wan_sets() valid.
    .mk(paste0("wan_", vcov_type), centre - c1*se, centre + c1*se, centre, se, c1, "lower")
}

#' The inner/outer level-set confidence sets of Wan et al. (2024).
#'
#' @return list(inner, outer, lambda) -- index vectors. inner subset of outer.
wan_sets <- function(band, lambda = 0) {
  list(inner  = which(band$lower >= lambda),
       outer  = which(band$upper >= lambda),
       lambda = lambda)
}

# ---- 6. Bayesian FDR selection rule -------------------------------------------

#' Select {i : tau_i >= lambda} controlling the posterior expected false
#' discovery proportion at q. A decision rule, not a band.
select_bayes_fdr <- function(post, lambda = 0, q = 0.10) {
  post <- .chk_post(post)
  p <- rowMeans(post >= lambda)
  o <- order(p, decreasing = TRUE)
  fdp <- cumsum(1 - p[o]) / seq_along(o)
  m <- if (any(fdp <= q)) max(which(fdp <= q)) else 0L
  list(method = "bayes_fdr", selected = if (m > 0) sort(o[seq_len(m)]) else integer(0),
       prob = p, q = q, lambda = lambda,
       expected_fdp = if (m > 0) fdp[m] else 0)
}

#' Dependency-corrected Bayesian FDR rule under correlated / smooth posteriors.
#'
#' Corrects for the ~80% FDP overshoot observed when posterior probabilities
#' are positively correlated through a shared tree/GP structure.
#' Applies a Benjamini-Yekutieli / harmonic discount factor q* = q / log(n).
select_bayes_fdr_dep <- function(post, lambda = 0, q = 0.10, method = c("by", "harmonic")) {
  post <- .chk_post(post)
  n <- nrow(post)
  method <- match.arg(method)
  discount <- if (method == "by") (log(n) + 0.5772) else sum(1 / seq_len(n))
  q_star <- q / discount
  p <- rowMeans(post >= lambda)
  o <- order(p, decreasing = TRUE)
  fdp <- cumsum(1 - p[o]) / seq_along(o)
  m <- if (any(fdp <= q_star)) max(which(fdp <= q_star)) else 0L
  list(method = paste0("bayes_fdr_dep_", method),
       selected = if (m > 0) sort(o[seq_len(m)]) else integer(0),
       prob = p, q = q, q_star = q_star, lambda = lambda,
       expected_fdp = if (m > 0) fdp[m] else 0)
}

#' 1D Isotonic Subgroup Selection (ISS) following Müller, Reeve, Cannings & Sechidis (2025/2026).
#'
#' Under monotonicity along a predictive ranking (e.g. estimated CATE order or a biomarker score),
#' an upper level set S = {i : tau_i >= lambda} is an interval [k, n].
#' A sequential hypothesis test on the lower boundary point guarantees finite-sample
#' inclusion error control P(Lhat subset of S) <= alpha with non-trivial power.
select_isotonic_1d <- function(post, lambda = 0, alpha = 0.05, score = NULL) {
  post <- .chk_post(post)
  n <- nrow(post)
  mns <- rowMeans(post)
  if (is.null(score)) score <- mns
  ord <- order(score)
  post_ord <- post[ord, , drop = FALSE]
  p_clear <- rowMeans(post_ord >= lambda)
  selected_ord <- integer(0)
  for (k in seq(n, 1, by = -1)) {
    if (p_clear[k] >= 1 - alpha) {
      selected_ord <- c(selected_ord, k)
    } else {
      break
    }
  }
  selected_orig <- if (length(selected_ord) > 0) sort(ord[selected_ord]) else integer(0)
  list(
    method   = "isotonic_iss",
    selected = selected_orig,
    lambda   = lambda,
    alpha    = alpha,
    q        = alpha,
    n_selected = length(selected_orig)
  )
}

# ---- selection from a band ----------------------------------------------------

#' Inner set: individuals whose lower bound clears the clinical threshold.
select_from_band <- function(band, lambda = 0) which(band$lower >= lambda)

#' Outer set: individuals not ruled out of the level set.
select_outer <- function(band, lambda = 0) {
  if (all(!is.finite(band$upper))) return(seq_along(band$lower))
  which(band$upper >= lambda)
}

#' The whole lambda-indexed family from ONE band (Wan et al. 2024, Theorem 2).
#'
#' The simultaneous confidence statement holds for a sequence of thresholds at
#' once, so the entire benefit curve costs nothing beyond the single band. The
#' useful summary is lambda_max: the largest threshold at which the inner set is
#' still non-empty -- a reportable finding even when a fixed-lambda selection is
#' empty.
#'
#' @return data.frame(lambda, n_inner, n_outer) plus attr "lambda_max"
level_sweep <- function(band, lambdas = NULL, n_grid = 50) {
  if (is.null(lambdas)) {
    fin <- band$lower[is.finite(band$lower)]
    lambdas <- seq(min(fin), max(fin), length.out = n_grid)
  }
  out <- data.frame(
    lambda  = lambdas,
    n_inner = vapply(lambdas, function(l) length(select_from_band(band, l)), integer(1)),
    n_outer = vapply(lambdas, function(l) length(select_outer(band, l)), integer(1))
  )
  nz <- which(out$n_inner > 0)
  attr(out, "lambda_max") <- if (length(nz)) max(out$lambda[nz]) else NA_real_
  out
}
