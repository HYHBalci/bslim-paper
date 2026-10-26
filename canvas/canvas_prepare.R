# ==============================================================================
# CANVAS / CANVAS-R  --  ANALYSIS-SET PREPARATION AND THE INTERPRETABILITY FIT
#
# Replaces "R/first (2).R". Two jobs, in this order:
#
#   1. Build the analysis set from the raw ADaM/SDTM files and SAVE IT, for BOTH
#      trials, through ONE code path.
#   2. Fit BCF and the semi-parametric hc-hs model on the training trial and
#      produce the interpretability output: the Shapley figures, the posterior of
#      beta_BLA1C, and the CATE-against-BLA1C scatter.
#
# The .RData files written by step 1 are what new_OOS.R consumes, so a defect
# here corrupts BOTH case studies. Step 1 therefore runs first, saves first, and
# asserts loudly; step 2 cannot damage it.
#
# WHY THIS IS PARAMETERISED BY STUDY. "first (2).R" builds one trial. The second
# (analsis_2.RData) must come from a near-copy of it, and nothing in the repo
# writes that file -- so the two trials' recoding lives in two places that have
# to agree by hand. They do not have to agree for new_OOS.R to run: stacking the
# two frames aligns the factor LEVELS THAT ARE PRESENT, never the recoding
# DECISIONS. If one trial codes SEX as "F"/"M" and the other as
# "FEMALE"/"MALE", ifelse(SEX == "F", 0, 1) silently makes every patient in the
# second trial a 1, the design matrices disagree, and the external validation is
# quietly comparing two different encodings. build_analysis_set() below is
# called twice so that cannot happen.
#
# Defects fixed, each marked FIX Pn at its call site. The numbering is separate
# from new_OOS.R's FIX 1-11 so the two sets never get confused.
#
#   FIX P1  Recoding ran AFTER complete.cases(), and two of the recodes inject
#           NA: factor(RACE, levels = ...) silently maps any unlisted level to
#           NA, and as.numeric(BLACR) does the same for any non-numeric string.
#           Those NAs were never removed; model.matrix() then dropped the rows
#           and X_final came out SHORTER than y and Z. Recode first, filter
#           second, assert third.
#   FIX P2  ifelse(x == "N", 0, 1) maps every UNEXPECTED value to 1, so a flag
#           coded "NO" rather than "N" silently becomes "yes" for everyone.
#           Same shape of bug in the SEX and ARM recodes. Now a strict map that
#           errors on an unrecognised level.
#   FIX P3  One code path for both trials -- see above.
#   FIX P4  indices = 1:5018 was hard-coded in the Shapley call.
#   FIX P5  Posterior columns were addressed by number (Beta[1, , 8],
#           X_final[, 14]). Resolved by NAME now, with an explicit intercept
#           check. See beta_col().
#   FIX P6  Both models receive the same constant treatment proportion from
#           the analysis training pool as a prognostic-forest covariate only,
#           matching new_OOS.R. It is not added to the treatment design.
#   FIX P7  bcf() was unseeded; only the semi-parametric fit had random_seed.
#   FIX P8  VISIT == "WEEK 52" was an unchecked exact string match. Zero matches
#           produced no W52_ column and a confusing downstream failure.
#   FIX P9  Subjects with more than one week-52 record were silently averaged.
#           Still averaged -- but counted and reported.
#   FIX P10 The outcome switch had a dead `else { print('skip this') }` branch
#           that left y undefined.
#   FIX P11 save(analsis) happened AFTER the modelling, so a failure in the
#           Shapley step threw away the data preparation too.
#   FIX P12 compute_shapley_all() was called with a boolean_interaction
#           argument that the repo's shapley_aux.R does not accept. The call is
#           now built from the function's actual formals.
#
# KEPT FROM YOUR EDIT: sample_global_prior = "hc-hs" (your FIX 1). The published
# CANVAS numbers used "half-cauchy", so everything downstream of this script has
# to be regenerated, not merged with the old output.
#
# Usage:  Rscript canvas_prepare.R
# ==============================================================================

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(bslim)
})

# ==============================================================================
# CONFIG -- the only block you should need to touch.  On YODA, set
# BCF_WORKSPACE and (when needed) BCF_DATA_ROOT before calling Rscript instead
# of editing paths throughout the script.
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

WORKSPACE_DIR <- Sys.getenv("BCF_WORKSPACE", unset = "X:/Workspace")
DATA_ROOT <- Sys.getenv("BCF_DATA_ROOT", unset = "X:/Data")
CFG <- list(
  # One entry per trial. `label` names the object saved inside the .RData; both
  # files hold an object called `analsis`, which is what new_OOS.R expects.
  # FIX P3: both go through build_analysis_set(), so the recoding is identical
  # by construction rather than by inspection.
  studies = list(
    train = list(dir  = file.path(DATA_ROOT, "JnJ-28431754DIA4003xpt6/Data"),
                 out  = file.path(WORKSPACE_DIR, "analsis.RData"),
                 name = "CANVAS-R"),
    test  = list(dir  = file.path(DATA_ROOT, "JnJ-28431754DIA3008xpt6/Data"),   # <-- CHECK THIS PATH
                 out  = file.path(WORKSPACE_DIR, "analsis_2.RData"),
                 name = "CANVAS")
  ),
  # Set FALSE to prepare only the training trial (e.g. if the second study's
  # directory is not mounted). new_OOS.R needs both.
  prepare_test_trial = TRUE,
  # FALSE (the default): analsis.RData / analsis_2.RData ALREADY EXIST and are
  # ready -- do not touch the raw XPT files or the recoding logic, just load
  # what is on disk at `studies$*$out` and run the modelling/interpretability
  # step (BCF, semi-parametric hc-hs, Shapley, the beta_BLA1C posterior) on that.
  # This is the "rerun first (2).R's modelling on data that is already correct"
  # path. Set TRUE only when the analysis sets themselves need rebuilding from
  # the raw ADaM/SDTM files -- that path still goes through build_analysis_set()
  # and every FIX P1-P12 guard.
  rebuild_from_raw = FALSE,
  # Refuse to overwrite an analysis set that already exists. new_OOS.R has
  # already run successfully against the current files; this script has never
  # touched the real data. The first run therefore writes <name>_v2.RData beside
  # the original and prints a comparison, so a working pipeline cannot be broken
  # by an unverified one. Set TRUE only once the comparison looks right.
  overwrite_existing = FALSE,

  out_dir = file.path(canvas_run_root(WORKSPACE_DIR), "interpretation"),
  shapley_src = Sys.getenv("BCF_SHAPLEY_SRC", unset = file.path(CANVAS_CODE_DIR, "shapley_aux.R")),

  # --- outcome ------------------------------------------------------------
  # FIX P10: an explicit table, so an unknown name stops rather than leaving y
  # undefined. Y is the CHANGE from baseline, so lower is better for HbA1c --
  # the same convention new_OOS.R encodes as benefit_sign = -1.
  outcome = "HBA1C",
  outcome_defs = list(
    HBA1C = list(testcd = "HBA1C",    baseline = "BLA1C", label = "HbA1c"),
    GFR   = list(testcd = "GFRBSCRT", baseline = "BLGFR", label = "eGFR")
  ),
  visit = "WEEK 52",

  x_cols = c("BLBMI","AGE","HHFFL","BLSBP","SEX","DIADUR","HCEREBFL","BLA1C",
             "BLGFR","INSULBL","GLP1BL","BLACR","RACE","METFORBL","RANGIBL"),

  # --- modelling ----------------------------------------------------------
  fit_models = TRUE,
  seed = 1848L,
  mcmc = canvas_mcmc_settings(),

  # --- interpretability ---------------------------------------------------
  focus_var = "BLA1C",     # the effect modifier the manuscript reports
  trace_var = "AGE",
  n_tolerance_sample = 1000
)

RACE_LEVELS <- c("WHITE",
                 "BLACK OR AFRICAN AMERICAN",
                 "ASIAN",
                 "AMERICAN INDIAN OR ALASKA NATIVE",
                 "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER",
                 "OTHER/MULTIPLE")
RACE_MISSING  <- c("NOT REPORTED", "UNKNOWN", "")
RACE_COLLAPSE <- c("OTHER", "MULTIPLE")

YN_MAP  <- c(N = 0, Y = 1)
SEX_MAP <- c(F = 0, M = 1)
# ARM_MAP doubles as the inclusion criterion -- see build_analysis_set(). If a
# trial labels its treated arm by dose ("Canagliflozin 100 mg"), add those
# labels here; the script tells you which labels it dropped.
ARM_MAP <- c(Placebo = 0, Canagliflozin = 1)
YN_VARS <- c("HHFFL", "HCEREBFL", "INSULBL", "GLP1BL", "METFORBL", "RANGIBL")

say <- function(...) { cat(sprintf(...)); flush.console() }

# ==============================================================================
# SECTION 1 -- HELPERS
# ==============================================================================

# ------------------------------------------------------------------------------
# standardize_X_by_index: IDENTICAL to the copy in new_OOS.R, including its
# scale_idx argument and its two FIX 4 guards. It is duplicated rather than
# sourced because both files have to travel to the YODA VM on their own. If you
# change one, change the other -- the interpretability design matrix and the
# out-of-sample design matrix MUST be produced by the same code, or the Shapley
# decomposition describes a different model than the one that was validated.
# ------------------------------------------------------------------------------
standardize_X_by_index <- function(X_initial,
                                   process_data = TRUE,
                                   interaction_rule = c("continuous", "continuous_or_binary", "all"),
                                   cat_coding_method = c("sum", "difference"),
                                   scale_idx = NULL) {
  interaction_rule <- match.arg(interaction_rule)
  cat_coding_method <- match.arg(cat_coding_method)
  if (cat_coding_method == "difference" && !requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for difference coding (contr.sdif).", call. = FALSE)
  }
  if (!is.data.frame(X_initial) && !is.matrix(X_initial)) stop("X_initial must be a data.frame or a matrix.")
  if (anyNA(X_initial)) stop("X_initial contains NA; take complete cases first (model.matrix would drop rows).")

  original_colnames <- colnames(X_initial)
  if (is.null(original_colnames)) original_colnames <- paste0("V", 1:ncol(X_initial))
  X_df <- as.data.frame(X_initial); colnames(X_df) <- original_colnames
  if (is.null(scale_idx)) scale_idx <- seq_len(nrow(X_df))

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
        contrasts(col_as_factor) <- if (cat_coding_method == "difference")
          MASS::contr.sdif(nlevels(col_as_factor)) else contr.sum(nlevels(col_as_factor))
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
  stopifnot(nrow(X_final) == nrow(X_initial))

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
       non_continous_idx_cpp = non_continous_idx_cpp,
       X_final_var_info = X_final_var_info)
}

# ------------------------------------------------------------------------------
# FIX P2. ifelse(x == "N", 0, 1) is not a recode, it is a default-to-1. Anything
# that is not exactly "N" -- "NO", "n", " N", a trailing-space variant, a value
# that only appears in the other trial -- becomes 1, for every affected patient,
# with no message. This maps explicitly and refuses anything it does not know.
# ------------------------------------------------------------------------------
recode_strict <- function(x, map, varname, study) {
  chr <- trimws(as.character(x))
  chr[chr == ""] <- NA
  unknown <- setdiff(unique(chr[!is.na(chr)]), names(map))
  if (length(unknown))
    stop(sprintf("%s: %s has unexpected value(s) {%s}. Expected {%s}. Add them to the map deliberately -- do not let them default.",
                 study, varname, paste(unknown, collapse = ", "), paste(names(map), collapse = ", ")),
         call. = FALSE)
  out <- unname(map[chr])
  out[is.na(chr)] <- NA
  out
}

# ------------------------------------------------------------------------------
# FIX P1 (first half). factor(x, levels = L) maps every value not in L to NA,
# silently. In the original that call ran AFTER complete.cases(), so those NAs
# survived into the design matrix. Here anything unrecognised is an error, and
# the deliberate missing categories are counted and reported.
# ------------------------------------------------------------------------------
collapse_race <- function(x, study) {
  chr <- trimws(toupper(as.character(x)))
  n_missing_cat <- sum(chr %in% RACE_MISSING, na.rm = TRUE)
  chr[chr %in% RACE_MISSING] <- NA
  chr[chr %in% RACE_COLLAPSE] <- "OTHER/MULTIPLE"
  unknown <- setdiff(unique(chr[!is.na(chr)]), RACE_LEVELS)
  if (length(unknown))
    stop(sprintf("%s: RACE has unrecognised level(s) {%s}. factor() would have turned these into NA without saying so.",
                 study, paste(unknown, collapse = ", ")), call. = FALSE)
  if (n_missing_cat)
    say("    RACE: %d patient(s) with %s treated as missing (they are dropped by the complete-case filter)\n",
        n_missing_cat, paste(RACE_MISSING[RACE_MISSING != ""], collapse = "/"))
  factor(chr, levels = RACE_LEVELS)
}

# ------------------------------------------------------------------------------
# FIX P1 (second half). as.numeric("<0.5") is NA with only a warning. In the
# original this ran after complete.cases() too, so those NAs reached
# model.matrix() and shortened X_final.
# ------------------------------------------------------------------------------
to_numeric_strict <- function(x, varname, study) {
  if (is.numeric(x)) return(as.numeric(x))
  chr <- trimws(as.character(x))
  out <- suppressWarnings(as.numeric(chr))
  bad <- unique(chr[!is.na(chr) & is.na(out)])
  if (length(bad))
    say("    %s: %s has %d value(s) that are not numeric {%s}; they become NA and are dropped by the complete-case filter\n",
        study, varname, sum(!is.na(chr) & is.na(out)), paste(utils::head(bad, 5), collapse = ", "))
  out
}

# ------------------------------------------------------------------------------
# FIX P5. The original addressed posterior coefficients positionally:
# Beta[1, , 8] for BLA1C and Beta[1, , 2] for AGE. Two ways that goes wrong.
#
# (a) X_final[, 14] was read as METFORBL. It is not. RACE sits at position 13 of
#     x_cols with 6 levels, so contr.sum expands it to 5 columns and everything
#     after it shifts by 4: METFORBL is column 18, and column 14 is the second
#     RACE contrast.
#
# (b) Beta may or may not carry a leading intercept. predict_linear_bcf_patched()
#     explicitly handles both (`if (p_beta == p_linear + 1) X_target <- cbind(1,
#     X_linear)`), which is the tell. If it does here, Beta[1, , 8] is HCEREBFL
#     rather than BLA1C -- and the posterior figure in the manuscript would be
#     of the wrong variable, with nothing to reveal it.
#
# This resolves a name to an index and refuses if the offset is not 0 or 1.
# ------------------------------------------------------------------------------
beta_col <- function(fit, X_final, varname) {
  nm <- colnames(X_final)
  j <- match(varname, nm)
  if (is.na(j))
    stop(sprintf("'%s' is not a design column. Available: %s", varname, paste(nm, collapse = ", ")),
         call. = FALSE)
  p_beta <- dim(fit$Beta)[3]; p_x <- ncol(X_final)
  off <- p_beta - p_x
  if (!off %in% c(0L, 1L))
    stop(sprintf("Beta has %d columns for a %d-column design; cannot align (offset %d).",
                 p_beta, p_x, off), call. = FALSE)
  if (off == 1L)
    say("  note: Beta carries a leading intercept, so design column %d is Beta column %d.\n", j, j + 1L)
  j + off
}

fit_guarded <- function(label, expr) {
  r <- try(expr, silent = TRUE)
  if (inherits(r, "try-error")) {
    say("  [SKIP] %s failed: %s\n", label, trimws(conditionMessage(attr(r, "condition"))))
    return(NULL)
  }
  r
}

# ==============================================================================
# SECTION 2 -- BUILD ONE TRIAL'S ANALYSIS SET
#
# FIX P3: this is the ONLY place either trial is constructed, so the two cannot
# drift apart. Every filtering step is counted and printed, so the drop from
# randomised patients to analysed patients is auditable in one place rather than
# reconstructed afterwards.
# ==============================================================================
build_analysis_set <- function(study) {
  say("\n--- %s (%s) ---\n", study$name, study$dir)
  od <- CFG$outcome_defs[[CFG$outcome]]
  if (is.null(od))
    stop("Unknown outcome '", CFG$outcome, "'. Known: ",
         paste(names(CFG$outcome_defs), collapse = ", "), call. = FALSE)

  f_lb <- file.path(study$dir, "lb.xpt"); f_sl <- file.path(study$dir, "adsl.xpt")
  for (f in c(f_lb, f_sl)) if (!file.exists(f)) stop("cannot find ", f, call. = FALSE)
  lb <- haven::read_xpt(f_lb)
  adsl <- haven::read_xpt(f_sl)
  say("  adsl: %d rows | lb: %d rows\n", nrow(adsl), nrow(lb))

  # --- week-52 labs -------------------------------------------------------
  # FIX P8: an exact string match on VISIT with no check. If the visit label is
  # "Week 52" or carries a day suffix, this silently yields zero rows, no W52_
  # column is created, and the failure surfaces much later as an unrelated error.
  visits <- unique(as.character(lb$VISIT))
  hit <- lb$VISIT == CFG$visit & lb$LBTESTCD == od$testcd
  if (!any(hit, na.rm = TRUE))
    stop(sprintf("no rows with VISIT == '%s' and LBTESTCD == '%s'. VISIT values present: %s",
                 CFG$visit, od$testcd, paste(utils::head(sort(visits), 20), collapse = " | ")),
         call. = FALSE)

  lb_v <- lb %>%
    filter(VISIT == CFG$visit, LBTESTCD == od$testcd) %>%
    tidyr::drop_na(LBSTRESN)

  # FIX P9: the original averaged repeated week-52 records without saying so.
  # Averaging is still the rule -- it just gets counted now.
  dup <- lb_v %>% count(USUBJID) %>% filter(n > 1)
  if (nrow(dup))
    say("  %d subject(s) have >1 %s record at %s; the mean is used (max %d records)\n",
        nrow(dup), od$testcd, CFG$visit, max(dup$n))

  lb_wide <- lb_v %>%
    group_by(USUBJID) %>%
    summarise(W52 = mean(LBSTRESN, na.rm = TRUE), .groups = "drop")
  say("  %s at %s available for %d subject(s)\n", od$testcd, CFG$visit, nrow(lb_wide))

  d <- adsl %>% left_join(lb_wide, by = "USUBJID")
  if (!od$baseline %in% colnames(d))
    stop("adsl has no baseline column '", od$baseline, "'", call. = FALSE)
  d$y <- as.numeric(d$W52) - as.numeric(d[[od$baseline]])

  # --- recode BEFORE filtering -------------------------------------------
  # FIX P1: this is the whole point of the ordering. Every recode that can
  # introduce NA runs here; complete.cases() runs afterwards; the design matrix
  # is asserted to have the same number of rows as the frame. In the original
  # the recodes ran after the filter, so the NAs they injected were never
  # removed and model.matrix() silently dropped those rows from X but not from
  # y or Z.
  miss <- setdiff(c(CFG$x_cols, "ARM"), colnames(d))
  if (length(miss)) stop(study$name, ": adsl is missing ", paste(miss, collapse = ", "), call. = FALSE)

  # ARM_MAP is BOTH the recoding and the inclusion criterion: a trial can carry
  # arms that are not part of this comparison (CANVAS randomised two
  # canagliflozin doses plus placebo, and adsl can list other therapies), and the
  # original dropped them with filter(ARM %in% c("Placebo", "Canagliflozin")).
  # They are dropped here too -- but the labels are PRINTED, because the failure
  # mode that matters is a trial whose treated arm is labelled differently (by
  # dose, say). Silently dropping those would leave an all-placebo trial.
  arm_chr <- trimws(as.character(d$ARM))
  in_arms <- arm_chr %in% names(ARM_MAP)
  dropped_arms <- sort(unique(arm_chr[!in_arms & !is.na(arm_chr)]))
  if (length(dropped_arms))
    say("    arms present but not analysed: %s  (%d patient(s) dropped)\n",
        paste(dropped_arms, collapse = " | "), sum(!in_arms, na.rm = TRUE))
  d <- d[in_arms, , drop = FALSE]
  if (!nrow(d))
    stop(sprintf("%s: no patient is in an analysed arm. ARM_MAP expects {%s}; the trial contains {%s}.",
                 study$name, paste(names(ARM_MAP), collapse = ", "),
                 paste(dropped_arms, collapse = ", ")), call. = FALSE)
  d$ARM <- recode_strict(d$ARM, ARM_MAP, "ARM", study$name)
  if (length(unique(stats::na.omit(d$ARM))) < 2)
    stop(sprintf("%s: only one arm survives the ARM filter. ARM_MAP expects {%s}; the labels dropped were {%s}. If this trial names its treated arm differently (by dose, for instance), add those labels to ARM_MAP -- otherwise the trial is all-placebo and every downstream number is meaningless.",
                 study$name, paste(names(ARM_MAP), collapse = ", "),
                 paste(dropped_arms, collapse = ", ")), call. = FALSE)
  d$SEX <- recode_strict(d$SEX, SEX_MAP, "SEX", study$name)
  for (v in YN_VARS) d[[v]] <- recode_strict(d[[v]], YN_MAP, v, study$name)
  d$RACE  <- collapse_race(d$RACE, study$name)
  d$BLACR <- to_numeric_strict(d$BLACR, "BLACR", study$name)
  for (v in c("BLBMI","AGE","BLSBP","DIADUR","BLA1C","BLGFR"))
    d[[v]] <- to_numeric_strict(d[[v]], v, study$name)

  # --- filter -------------------------------------------------------------
  need <- c("y", "ARM", CFG$x_cols)
  n0 <- nrow(d)
  d <- d[!is.na(d$ARM), , drop = FALSE];             n1 <- nrow(d)
  d <- d[!is.na(d$y), , drop = FALSE];               n2 <- nrow(d)
  d <- d[complete.cases(d[, need, drop = FALSE]), , drop = FALSE]; n3 <- nrow(d)
  say("  patients: %d randomised -> %d in an analysed arm -> %d with an outcome -> %d complete on all covariates\n",
      n0, n1, n2, n3)
  if (n3 == 0) stop(study$name, ": no patients survive the filters", call. = FALSE)

  # USUBJID is carried through so two prepared versions can be compared patient
  # by patient. new_OOS.R subsets to its own needed_cols, so the extra column is
  # inert there, and the file never leaves the VM that already holds the raw data.
  keep_cols <- c(if ("USUBJID" %in% colnames(d)) "USUBJID", need)
  analsis <- as.data.frame(d[, keep_cols, drop = FALSE])
  analsis$y <- as.numeric(analsis$y)
  analsis$ARM <- as.numeric(analsis$ARM)
  stopifnot(all(analsis$ARM %in% c(0, 1)), !anyNA(analsis))

  a1 <- analsis$y[analsis$ARM == 1]; a0 <- analsis$y[analsis$ARM == 0]
  say("  n = %d, treated %.3f, ATE = %+.4f (se %.4f)   [%s change; negative = benefit]\n",
      nrow(analsis), mean(analsis$ARM), mean(a1) - mean(a0),
      sqrt(var(a1)/length(a1) + var(a0)/length(a0)), od$label)
  say("  RACE: %s\n", paste(sprintf("%s=%d", levels(analsis$RACE), as.integer(table(analsis$RACE))), collapse = ", "))
  analsis
}

# ==============================================================================
# SECTION 3 -- PREPARE AND SAVE BOTH TRIALS, BEFORE ANY MODELLING
#
# FIX P11: the original saved analsis at the very END, after the Shapley step,
# so any modelling failure discarded the data preparation as well. The .RData
# files are the handoff to new_OOS.R and are the more valuable artefact; they
# are written first and nothing below can touch them.
# ==============================================================================
CFG$out_dir <- canvas_output_directory(WORKSPACE_DIR, "interpretation")
say("\n======================================================\n")
say("CANVAS analysis-set preparation  (outcome: %s)\n", CFG$outcome)
say("======================================================\n")

to_do <- if (CFG$prepare_test_trial) c("train", "test") else "train"
prepared <- list()

# ------------------------------------------------------------------------------
# Load an existing analysis set from disk. Asserts the same contract
# build_analysis_set() guarantees, so the modelling code below cannot tell which
# path produced its input -- but it CANNOT catch a file that was prepared by a
# different, unaudited process (e.g. the original first (2).R). That file would
# load cleanly and still carry FIX P1/P2's silent-NA and default-to-1 defects.
# If in doubt, run once with rebuild_from_raw = TRUE and diff against it.
# ------------------------------------------------------------------------------
load_existing_set <- function(study) {
  if (!file.exists(study$out))
    stop(study$name, ": rebuild_from_raw = FALSE but ", study$out,
         " does not exist. Either point studies$*$out at the right file or set ",
         "rebuild_from_raw = TRUE to build it from the raw data.", call. = FALSE)
  e <- new.env(); load(study$out, envir = e)
  if (!"analsis" %in% ls(e))
    stop(study$out, " does not contain an object called 'analsis'.", call. = FALSE)
  d <- e$analsis
  need <- c("y", "ARM", CFG$x_cols)
  miss <- setdiff(need, colnames(d))
  if (length(miss)) stop(study$out, " is missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  if (anyNA(d[, need, drop = FALSE])) stop(study$out, " contains NA in a required column.", call. = FALSE)
  if (!all(d$ARM %in% c(0, 1))) stop(study$out, ": ARM is not coded 0/1.", call. = FALSE)
  if (!is.factor(d$RACE)) stop(study$out, ": RACE is not a factor.", call. = FALSE)
  a1 <- d$y[d$ARM == 1]; a0 <- d$y[d$ARM == 0]
  say("\n--- %s (loaded from %s) ---\n", study$name, study$out)
  say("  n = %d, treated %.3f, ATE = %+.4f (se %.4f)\n", nrow(d), mean(d$ARM),
      mean(a1) - mean(a0), sqrt(var(a1)/length(a1) + var(a0)/length(a0)))
  d
}

if (isTRUE(CFG$rebuild_from_raw)) {
  for (k in to_do) prepared[[k]] <- build_analysis_set(CFG$studies[[k]])
} else {
  say("\nrebuild_from_raw = FALSE: loading the existing analysis set(s) -- nothing\n")
  say("in the raw-data recoding runs, and nothing is re-saved.\n")
  for (k in to_do) prepared[[k]] <- load_existing_set(CFG$studies[[k]])
}

# Report any RACE level with no patients in this run. Not a fix -- see the
# comment at load_existing_set()/build_analysis_set() call site history: an
# unused level does not make the design rank-deficient, because
# standardize_X_by_index() calls factor(col_data) before assigning contr.sum
# contrasts, and factor() on an already-a-factor column drops unused levels on
# its own. This is informational only.
used <- unique(unlist(lapply(prepared, function(d) as.character(d$RACE))))
if (any(!RACE_LEVELS %in% used))
  say("\n  RACE level(s) with no patients in %s: %s (harmless -- dropped automatically downstream)\n",
      if (length(to_do) == 2) "either trial" else "this run",
      paste(setdiff(RACE_LEVELS, used), collapse = ", "))
keep_lev <- RACE_LEVELS[RACE_LEVELS %in% used]
if (length(keep_lev) < 2) stop("RACE has fewer than 2 observed levels.", call. = FALSE)

# Compare a freshly prepared set against one already on disk. Prints the things
# that would change downstream results; says nothing reassuring that it has not
# checked.
compare_sets <- function(old, new, label) {
  say("  comparison for %s (existing -> new)\n", label)
  say("    n: %d -> %d\n", nrow(old), nrow(new))
  co <- setdiff(names(old), names(new)); cn <- setdiff(names(new), names(old))
  if (length(co)) say("    columns only in existing: %s\n", paste(co, collapse = ", "))
  if (length(cn)) say("    columns only in new:      %s\n", paste(cn, collapse = ", "))
  if (!is.null(old$RACE) && !is.null(new$RACE)) {
    lo <- levels(old$RACE); ln <- levels(new$RACE)
    if (!identical(lo, ln))
      say("    RACE levels: %d -> %d  (dropped: %s)\n", length(lo), length(ln),
          paste(setdiff(lo, ln), collapse = ", "))
  }
  ate <- function(d) mean(d$y[d$ARM == 1]) - mean(d$y[d$ARM == 0])
  say("    treated fraction: %.4f -> %.4f\n", mean(old$ARM), mean(new$ARM))
  say("    ATE: %+.4f -> %+.4f\n", ate(old), ate(new))
  if (!is.null(old$USUBJID) && !is.null(new$USUBJID)) {
    say("    patients: %d only in existing, %d only in new, %d shared\n",
        length(setdiff(old$USUBJID, new$USUBJID)),
        length(setdiff(new$USUBJID, old$USUBJID)),
        length(intersect(old$USUBJID, new$USUBJID)))
  } else {
    say("    (the existing file has no USUBJID, so patient sets cannot be diffed)\n")
  }
}

if (isTRUE(CFG$rebuild_from_raw)) {
  for (k in to_do) {
    analsis <- prepared[[k]]
    stopifnot(!anyNA(analsis), nlevels(analsis$RACE) == length(keep_lev))
    target <- CFG$studies[[k]]$out
    if (file.exists(target) && !isTRUE(CFG$overwrite_existing)) {
      e <- new.env()
      ok <- !inherits(try(load(target, envir = e), silent = TRUE), "try-error")
      target <- sub("[.]RData$", "_v2.RData", target)
      say("\n  %s already exists and overwrite_existing = FALSE; writing %s instead\n",
          CFG$studies[[k]]$out, basename(target))
      if (ok && !is.null(e$analsis)) compare_sets(e$analsis, analsis, CFG$studies[[k]]$name)
    }
    save(analsis, file = target)
    say("  -> saved %s (%d rows, RACE with %d level(s))\n",
        target, nrow(analsis), nlevels(analsis$RACE))
  }
} else {
  say("\nrebuild_from_raw = FALSE: nothing re-saved -- %s\n",
      paste(vapply(to_do, function(k) CFG$studies[[k]]$out, ""), collapse = ", "))
}

if (length(to_do) == 2) {
  a <- prepared$train; b <- prepared$test
  say("\n  cross-trial check\n")
  if (nrow(a) < nrow(b))
    warning("The training trial has fewer patients than the test trial. new_OOS.R trains on the LARGER trial.",
            call. = FALSE, immediate. = TRUE)
  # FIX P3, verification half. Under rebuild_from_raw = TRUE, identical recoding
  # is guaranteed by construction; this proves it on the data actually in hand.
  # Under rebuild_from_raw = FALSE the two files were NOT necessarily produced
  # together, so this is a real check rather than a formality -- a level present
  # in one trial's saved levels() and absent from the other's shifts the
  # contrast columns, and new_OOS.R's rbind()-then-standardise step would
  # silently misalign the two designs.
  if (!identical(levels(a$RACE), levels(b$RACE)))
    stop(sprintf("RACE factor levels differ between the two saved files:\n  %s: %s\n  %s: %s\nThey were likely prepared by different code or at different times. Re-run with rebuild_from_raw = TRUE to regenerate both consistently.",
                 CFG$studies$train$name, paste(levels(a$RACE), collapse = ", "),
                 CFG$studies$test$name,  paste(levels(b$RACE), collapse = ", ")),
        call. = FALSE)
  for (v in c(YN_VARS, "SEX", "ARM")) {
    ua <- sort(unique(a[[v]])); ub <- sort(unique(b[[v]]))
    if (!identical(ua, ub))
      say("    note: %s takes {%s} in %s but {%s} in %s\n", v,
          paste(ua, collapse = ","), CFG$studies$train$name,
          paste(ub, collapse = ","), CFG$studies$test$name)
  }
  # A level carried by only ONE trial is real data and is kept -- but it is
  # perfectly confounded
  # with trial membership, so any effect attributed to it is really a
  # between-trial difference.
  ra <- table(a$RACE); rb <- table(b$RACE)
  one_sided <- names(ra)[(ra == 0) != (rb == 0)]
  if (length(one_sided))
    warning("RACE level(s) {", paste(one_sided, collapse = ", "), "} occur in only ",
            "one of the two trials, so they are perfectly confounded with trial ",
            "membership. Their coefficients and Shapley values are between-trial ",
            "differences, not effect modification.", call. = FALSE, immediate. = TRUE)
  say("    RACE levels match; %d train / %d test patients\n", nrow(a), nrow(b))
}
# The modelling below is wrapped in run_interpretability() rather than guarded
# by quit(), which would tear down an interactive RStudio session on the VM.
if (!isTRUE(CFG$fit_models))
  say("\nfit_models = FALSE; preparation only, nothing fitted.\n")

# ==============================================================================
# SECTION 4 -- THE INTERPRETABILITY FIT
#
# Fitted on the TRAINING trial only, on the whole sample. This is the fit the
# manuscript's Shapley figures and beta_BLA1C posterior come from.
# ==============================================================================
run_interpretability <- function() {

analsis <- prepared$train
X_mat <- analsis[, CFG$x_cols, drop = FALSE]
y <- as.numeric(analsis$y)
Z <- as.numeric(analsis$ARM)

res_X <- standardize_X_by_index(X_mat, process_data = TRUE,
                                interaction_rule = "continuous_or_binary",
                                cat_coding_method = "sum")
X_final <- res_X$X_final
storage.mode(X_final) <- "double"          # stochtree's C++ layer requires doubles
stopifnot(nrow(X_final) == length(y), length(y) == length(Z))
say("\n  design: %d x %d\n", nrow(X_final), ncol(X_final))
print(res_X$X_final_var_info[, c("col_name_final", "is_continuous", "is_binary", "is_categorical")],
      row.names = FALSE)

# FIX P5, demonstrated. RACE has 6 levels at position 13 of x_cols, so it
# expands to 5 contrast columns and every covariate after it shifts by 4. Anyone
# reading X_final[, 14] as METFORBL is reading the second RACE contrast.
say("\n  positional check -- what the old hard-coded indices actually pointed at:\n")
for (j in c(8, 14)) if (j <= ncol(X_final)) say("    X_final[, %d] = %s\n", j, colnames(X_final)[j])
say("    %s is column %d; METFORBL is column %d\n", CFG$focus_var,
    match(CFG$focus_var, colnames(X_final)), match("METFORBL", colnames(X_final)))

# FIX P6: pass the same constant analysis-pool treatment proportion to both
# models, only in the prognostic forest, as in new_OOS.R. The coefficient and
# interaction design remains X_final.
pi_hat <- rep(mean(Z), length(Z))
say("  randomisation propensity: %.4f\n", pi_hat[1])

# The same explicit Gibbs hierarchy and sampling budget are used in both analyses.
say("\n  fitting BCF: four independent chains...\n")
bcf_run <- canvas_fit_chains("BCF", X_final, y, Z, X_final, Z, pi_hat[1],
                             CFG$seed, CFG$mcmc)
tau_bcf <- bcf_run$tau_mean
write.csv(bcf_run$diagnostics, file.path(CFG$out_dir, "mcmc_diagnostics_BCF.csv"), row.names = FALSE)
write.csv(bcf_run$diagnostic_groups, file.path(CFG$out_dir, "mcmc_diagnostic_groups_BCF.csv"), row.names = FALSE)

say("  fitting semi-parametric (hc-hs): four independent chains...\n")
sp_run <- canvas_fit_chains("SemiParam", X_final, y, Z, X_final, Z, pi_hat[1],
                            CFG$seed, CFG$mcmc)
nbcf_fit <- sp_run$posterior
write.csv(sp_run$diagnostics, file.path(CFG$out_dir, "mcmc_diagnostics_SemiParam.csv"), row.names = FALSE)
write.csv(sp_run$diagnostic_groups, file.path(CFG$out_dir, "mcmc_diagnostic_groups_SemiParam.csv"), row.names = FALSE)
mcmc_manifest <- list(protocol = "canvas-independent-4chain-mu-propensity-8.3.5",
  BCF = bcf_run$mcmc, SemiParam = sp_run$mcmc,
  package = as.character(utils::packageVersion("bslim")),
  package_description = utils::packageDescription("bslim"),
  analysis_code_hashes = tools::md5sum(file.path(CANVAS_CODE_DIR, c("canvas_mcmc.R", "canvas_prepare.R"))),
  input_hashes = tools::md5sum(vapply(CFG$studies, `[[`, character(1), "out")))
saveRDS(mcmc_manifest, file.path(CFG$out_dir, "mcmc_manifest.rds"))

say("  ATE: BCF %+.4f | trial %+.4f\n",
    mean(tau_bcf, na.rm = TRUE), mean(y[Z == 1]) - mean(y[Z == 0]))

save(res_X, file = file.path(CFG$out_dir, "X_final_varinfo.RData"))
save(X_final, file = file.path(CFG$out_dir, "X_final.RData"))
my_list <- list(alpha = nbcf_fit$alpha, beta = nbcf_fit$Beta,
                beta_int = nbcf_fit$Beta_int, sigma = nbcf_fit$sigma2_global_samples,
                interaction_pairs = nbcf_fit$interaction_pairs,
                col_names = colnames(X_final), mcmc = sp_run$mcmc,
                model_params = nbcf_fit$model_params, scale_draws = nbcf_fit$scale_draws)
save(my_list, file = file.path(CFG$out_dir, "nbcf_fit.RData"), compress = "xz")
say("  saved posterior draws (with col_names, so nothing downstream needs positional indexing)\n")

# ==============================================================================
# SECTION 5 -- SHAPLEY DECOMPOSITION
# ==============================================================================
source(CFG$shapley_src)

# FIX P12: the original passed boolean_interaction =, which the repo's
# shapley_aux.R does not declare. Whether that call works depends on which
# variant of shapley.R the VM happens to have. Build the argument list from the
# function's actual formals instead of assuming.
shap_args <- list(X = X_final, beta_post = nbcf_fit$Beta,
                  beta_int_post = nbcf_fit$Beta_int,
                  ipairs = nbcf_fit$interaction_pairs,
                  indices = seq_len(nrow(X_final)))   # FIX P4: was 1:5018
fm <- names(formals(compute_shapley_all))
if ("boolean_interaction" %in% fm)
  shap_args$boolean_interaction <- as.logical(res_X$X_final_var_info$is_continuous |
                                              res_X$X_final_var_info$is_binary)
extra <- setdiff(names(shap_args), fm)
if (length(extra)) {
  say("  note: compute_shapley_all() does not accept {%s}; dropping\n", paste(extra, collapse = ", "))
  shap_args <- shap_args[intersect(names(shap_args), fm)]
}

shapley_results <- fit_guarded("Shapley", do.call(compute_shapley_all, shap_args))
if (!is.null(shapley_results)) {
  p1 <- plot_shapley_contribution(shapley_results)
  ggsave(file.path(CFG$out_dir, "contr_shap.png"), p1, width = 10, height = 8, dpi = 300, bg = "white")
  p2 <- plot_shapley_importance_breakdown(shapley_results)
  ggsave(file.path(CFG$out_dir, "break_shap.png"), p2, width = 10, height = 8, dpi = 300, bg = "white")
  say("  -> contr_shap.png, break_shap.png\n")

  set.seed(CFG$seed)
  idx <- sample(seq_len(nrow(X_final)), min(CFG$n_tolerance_sample, nrow(X_final)))
  post_focus <- fit_guarded("Shapley posterior", get_shapley_posterior_single(
    feature_name = CFG$focus_var, indices = idx, X = X_final,
    beta_post = nbcf_fit$Beta, beta_int_post = nbcf_fit$Beta_int,
    ipairs = nbcf_fit$interaction_pairs))
  if (!is.null(post_focus)) {
    p_tol <- plot_tolerance_intervals(post_matrix = post_focus,
                                      feature_values = X_final[idx, CFG$focus_var],
                                      feature_name = CFG$focus_var)
    ggsave(file.path(CFG$out_dir, sprintf("tolerance_%s.png", CFG$focus_var)),
           p_tol, width = 10, height = 6, dpi = 300, bg = "white")
    say("  -> tolerance_%s.png\n", CFG$focus_var)
  }
}

# ==============================================================================
# SECTION 6 -- POSTERIOR OF THE MAIN EFFECT MODIFIER
#
# FIX P5 in force: the coefficient is located by NAME and the intercept offset
# is checked, so this figure cannot silently become a different covariate.
# ==============================================================================
jb <- beta_col(nbcf_fit, X_final, CFG$focus_var)
draws <- as.numeric(nbcf_fit$Beta[, , jb])

# Back to the original HbA1c-per-unit scale. X_final's continuous columns were
# divided by their sd and the outcome was standardised internally, so the
# original-scale coefficient is beta * sd(y) / sd(x_raw).
# The outcome scale is taken from the FIT, not recomputed as sd(y). The model
# stores what it actually divided by (predict_linear_bcf_patched reads the same
# field), and assuming sd(y) would put every coefficient quoted in the paper off
# by a constant factor if stochtree ever standardises differently.
y_scale <- nbcf_fit$model_params$outcome_scale
if (is.null(y_scale) || !is.finite(y_scale)) {
  y_scale <- sd(y)
  say("  note: model_params$outcome_scale is unavailable; falling back to sd(y) = %.4f\n", y_scale)
} else if (abs(y_scale - sd(y)) > 1e-6 * max(1, sd(y))) {
  say("  note: outcome_scale = %.6f, sd(y) = %.6f -- using the model's value\n", y_scale, sd(y))
}
sd_raw <- sd(as.numeric(X_mat[[CFG$focus_var]]))
draws_orig <- draws * y_scale / sd_raw

qs <- quantile(draws_orig, c(0.025, 0.5, 0.975))
say("\n  beta_%s (original scale): mean %+.4f, median %+.4f, 95%% CrI [%+.4f, %+.4f], P(<0) = %.3f\n",
    CFG$focus_var, mean(draws_orig), qs[2], qs[1], qs[3], mean(draws_orig < 0))

p_post <- ggplot(data.frame(x = draws_orig), aes(x = x)) +
  geom_density(fill = "skyblue", alpha = 0.5, colour = "navy") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  annotate("segment", x = qs[1], xend = qs[3], y = 0, yend = 0, linewidth = 2, colour = "navy") +
  annotate("point", x = qs[2], y = 0, size = 4, colour = "firebrick") +
  labs(title = bquote(paste("Posterior of ", beta[.(CFG$focus_var)])),
       subtitle = sprintf("mean %+.3f, median %+.3f, 95%% CrI [%+.3f, %+.3f], P(< 0) = %.2f",
                          mean(draws_orig), qs[2], qs[1], qs[3], mean(draws_orig < 0)),
       x = sprintf("change in HbA1c per 1 unit of %s", CFG$focus_var), y = "Density") +
  theme_minimal(base_size = 13)
ggsave(file.path(CFG$out_dir, sprintf("posterior_%s_beta.png", CFG$focus_var)),
       p_post, width = 7, height = 5, dpi = 300, bg = "white")
say("  -> posterior_%s_beta.png\n", CFG$focus_var)

# Trace, also name-resolved
if (requireNamespace("coda", quietly = TRUE)) {
  jt <- beta_col(nbcf_fit, X_final, CFG$trace_var)
  tr <- coda::mcmc.list(lapply(seq_len(dim(nbcf_fit$Beta)[1]), function(chain)
    coda::as.mcmc(as.numeric(nbcf_fit$Beta[chain, , jt]) * y_scale / sd(as.numeric(X_mat[[CFG$trace_var]])))))
  png(file.path(CFG$out_dir, sprintf("trace_%s.png", CFG$trace_var)),
      width = 900, height = 500)
  coda::traceplot(tr, main = sprintf("beta_%s (original scale)", CFG$trace_var))
  dev.off()
  say("  -> trace_%s.png\n", CFG$trace_var)
}

# CATE against the effect modifier. The original plotted tau against
# X_final[, 14], which is a RACE contrast, not METFORBL (FIX P5).
if (!all(is.na(tau_bcf))) {
  pdf_df <- data.frame(x = X_mat[[CFG$focus_var]], cate = tau_bcf)
  p_cate <- ggplot(pdf_df, aes(x = x, y = cate)) +
    geom_point(alpha = 0.3, colour = "slategray") +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, colour = "firebrick") +
    labs(title = "Heterogeneous treatment effects (BCF)",
         subtitle = sprintf("CATE against baseline %s", CFG$focus_var),
         x = sprintf("Baseline %s", CFG$focus_var), y = "Estimated CATE") +
    theme_minimal(base_size = 13)
  ggsave(file.path(CFG$out_dir, sprintf("cate_vs_%s.png", CFG$focus_var)),
         p_cate, width = 7, height = 5, dpi = 300, bg = "white")
  say("  -> cate_vs_%s.png\n", CFG$focus_var)
}

say("\nDone. Analysis sets and figures are in %s\n", CFG$out_dir)
invisible(NULL)
}

if (isTRUE(CFG$fit_models)) run_interpretability()
