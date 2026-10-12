# ==============================================================================
# R/subgroup/compat.R -- drop-in replacements for the legacy band functions.
#
# WHY THIS EXISTS --------------------------------------------------------------
# The August 2026 consolidation fixed the band rules but never migrated the
# callers. As of commit 5eb272f, every one of these still defines its own
# findconfglob()/intfrompost1():
#
#   R/ACTG_complete.R      R/ACTG_fast.R        R/ACTG_new.R
#   R/complete_loop.R      R/crossval_ACTG.R    R/crossval_ACTG_New.R
#   stochtree/R/Shapley.R  (the local stochtree fork)
#
# So every number those scripts have produced comes from the superseded rule.
# This file makes migration a two-line edit: delete the local definitions and
#
#     source(file.path("R", "subgroup", "bands.R"))
#     source(file.path("R", "subgroup", "compat.R"))
#
# WHAT THE LEGACY CODE GOT WRONG -----------------------------------------------
# The old search was
#     grid <- seq(0.0005, 0.05, by = 0.0005)
#     wm   <- which(conftols - 0.95 < 0)[1] - 1
# which fails in three ways:
#
#  (1) The grid stops at 0.05, so the band could never be narrower than a
#      pointwise 95% interval however much the data supported.
#  (2) When the criterion holds across the whole grid, which() returns
#      integer(0) and wm becomes NA. The ACTG scripts then return 0.0005 -- the
#      grid MINIMUM, i.e. the WIDEST band -- when the correct answer is the
#      narrowest. Measured on a 300x6000 posterior at tol = 0.10: legacy gamma
#      0.0005 against a correct 0.0750, a band 1.9x too wide, selecting 5
#      patients where the fixed rule selects 38.
#  (3) stochtree/R/Shapley.R has no guard at all, so the same input returns NA
#      (and numeric(0) in the opposite branch), which propagates silently into
#      intfrompost1(post, conf = NA).
#
# The replacements below have a closed form, no grid, and no degenerate branch.
# ==============================================================================

if (!exists("band_tolerance")) stop("source bands.R before compat.R")

#' Legacy signature. Row summary c(lower, mean, upper) at marginal level `conf`.
postsum <- function(posters, conf = 0.05) {
  qs <- stats::quantile(posters, p = c(conf / 2, 1 - conf / 2))
  c(qs[1], mean(posters), qs[2])
}

#' Legacy signature, ACTG style: returns an n x 3 matrix.
intfrompost1 <- function(post, conf = 0.05) {
  t(apply(post, 1, postsum, conf = conf))
}

#' Legacy signature, stochtree/Shapley.R style: returns list(tot = <n x 3>).
#' Shapley.R calls intfrompost1(...)[[1]], which a bare matrix cannot serve.
intfrompost1_list <- function(post, conf = 0.05) {
  list(tot = intfrompost1(post, conf = conf))
}

#' Legacy signature. Fraction of draws whose content reaches 1 - tol.
conftol <- function(margconf = 0.05, post, tol = 0.05) {
  int  <- intfrompost1(post, conf = margconf)
  difl <- ((post - int[, 1]) >= 0) * ((post - int[, 3]) <= 0)
  mean(colMeans(difl) >= 1 - tol)
}

#' Legacy signature, correct implementation.
#'
#' Same call and same return type as the old findconfglob(); it just returns the
#' right number. Unlike the original it can return a gamma above 0.05, and it
#' never confuses "criterion met everywhere" with "criterion met nowhere".
#'
#' @param legacy_grid TRUE reproduces the pre-August-2026 search exactly,
#'   degenerate branches and 0.05 cap included. For regression checking only.
findconfglob <- function(post, tol = 0.05, conf = 0.95, legacy_grid = FALSE) {
  if (legacy_grid) {
    grid <- legacy_gamma_grid()
    ct <- vapply(grid, conftol, numeric(1), post = post, tol = tol)
    wm <- which(ct - conf < 0)[1] - 1
    if (is.na(wm) || wm < 1) return(grid[1])
    return(grid[wm])
  }
  band_tolerance(post, tol = tol, conf = conf, method = "exact")$level
}

#' The band the legacy pipeline was trying to build, in one call.
#'
#' Convenience for migrating plotting code: returns the same three objects the
#' old build_heterogeneity_plot() assembled by hand, plus the diagnostics the
#' old path discarded.
legacy_band_bundle <- function(post, tol = 0.05, conf = 0.95, alpha = 0.05) {
  tb <- band_tolerance(post, tol = tol, conf = conf)
  list(
    confglob = tb$level,
    cred_int = intfrompost1(post, conf = alpha),
    tol_int  = cbind(tb$lower, tb$centre, tb$upper),
    notes    = tb$notes
  )
}

# ---- the report's own question, with the real methods -------------------------

#' Drop-in replacement for check_significant_individuals() in
#' R/utility_functions.R, which competition.R uses to build the
#' Tolerance/Wan/Sechidis columns of papers/Full_Dataset_Report*.tex.
#'
#' The original computed, per individual:
#'   "Wan"      mu +/- t * sd/sqrt(S)          -- the MCSE of the posterior mean,
#'                                                ~55x too narrow at S = 3000;
#'                                                nothing to do with Wan (2024)
#'   "Sechidis" mu +/- 1.96 * sqrt(sd^2 + mu^2/S) -- a POINTWISE interval; no
#'                                                sup, no critical constant, no
#'                                                multiplicity adjustment
#'   "Tolerance" legacy findconfglob()          -- capped grid, degenerate branch
#'
#' This version answers the same question -- how many individuals have a CATE
#' that provably differs from `ref`? -- using bands that are what they are
#' called. The "Wan" column is dropped rather than renamed: Wan et al. (2024) is
#' a frequentist construction on (X, y) and cannot be computed from posterior
#' draws alone. Use subgroup_compare() in compare.R for that arm.
#'
#' @param het_draws n x S posterior draws of tau
#' @param ref       reference value; the original used 0
#' @return list with the same names the caller expects, plus what it was missing
check_significant_individuals <- function(het_draws, level = 0.95, tol = 0.05, ref = 0) {
  alpha <- 1 - level
  n <- nrow(het_draws)
  cnt <- function(b) sum(b$lower > ref | b$upper < ref, na.rm = TRUE)

  b_tol  <- band_tolerance(het_draws, tol = tol, conf = level)
  b_supt <- band_supt(het_draws, alpha = alpha, sided = "two", standardise = "sd")
  b_sq   <- band_supt(het_draws, alpha = alpha, sided = "two", standardise = "quantile")
  b_pw   <- band_pointwise(het_draws, alpha = alpha, sided = "two")
  b_bon  <- band_bonferroni(het_draws, alpha = alpha, sided = "two")

  list(
    tolerance_sig = cnt(b_tol),
    supt_sig      = cnt(b_supt),   # the simultaneous band the old "Sechidis" column claimed to be
    supt_q_sig    = cnt(b_sq),
    pointwise_sig = cnt(b_pw),     # what the old "Sechidis" column actually was
    bonferroni_sig= cnt(b_bon),
    wan_sig       = NA_integer_,   # needs (X, y); see subgroup_compare()
    sechidis_sig  = cnt(b_supt),   # legacy name, now pointing at a real sup-t band
    n_total       = n,
    k_one         = c(tolerance = b_tol$k_one, supt = b_supt$k_one, pointwise = b_pw$k_one),
    notes         = c(b_tol$notes, b_supt$notes, b_bon$notes)
  )
}
