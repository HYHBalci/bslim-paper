# Edit BEFORE starting a new run; changing this file invalidates old checkpoints.
config <- list(
  protocol = "spbart-8.3.5-interactions-v5",
  n_values = c(250L, 500L, 750L, 1000L, 1500L, 3000L),
  replicates = 1:200,
  ihdp_replicates = 1:100, ihdp_folds = 10L,
  chains = 4L, gfr = 25L, burnin = 2000L, draws = 4000L,
  seed_base = 20260921L,
  rhat_limit = 1.01, ess_min = 400,
  # Full CATE draws are used for exact pooling, then compacted only after a
  # verified result is saved. Keep TRUE if you need arbitrary future summaries.
  keep_chain_posteriors = FALSE,
  keep_ihdp_posteriors = TRUE,
  # Retain sampled scales aligned with the already saved coefficient chains
  # for joint inspection, without retaining every synthetic CATE draw.
  keep_scale_posteriors = TRUE,
  # Stop starting new fits early; the shell enforces a separate hard timeout.
  worker_seconds = 33000, reserve_seconds = 5400,
  ihdp_data_dir = "data/ihdp"
)
