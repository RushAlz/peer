#!/usr/bin/env Rscript
# Correctness tests for OpenMP parallelization of PEER
# Uses PEER_setNThreads() to control thread count within a single process

library(peer)

cat("=== PEER OpenMP Correctness Tests ===\n\n")

pass_count <- 0
fail_count <- 0

check <- function(name, condition) {
  if (condition) {
    cat(sprintf("  PASS: %s\n", name))
    pass_count <<- pass_count + 1
  } else {
    cat(sprintf("  FAIL: %s\n", name))
    fail_count <<- fail_count + 1
  }
}

run_peer <- function(N, P, Nk, Niter, threads, covs_cols = 0, seed = 42) {
  set.seed(seed)
  pheno <- matrix(rnorm(N * P), N, P)
  covs <- if (covs_cols > 0) matrix(rnorm(N * covs_cols), N, covs_cols) else NULL
  model <- PEER()
  PEER_setPhenoMean(model, pheno)
  PEER_setNk(model, Nk)
  PEER_setNmax_iterations(model, Niter)
  PEER_setNThreads(model, threads)
  if (!is.null(covs)) PEER_setCovariates(model, covs)
  PEER_update(model)
  list(
    X = PEER_getX(model),
    W = PEER_getW(model),
    Alpha = PEER_getAlpha(model),
    Eps = PEER_getEps(model),
    residuals = PEER_getResiduals(model),
    bounds = PEER_getBounds(model)
  )
}

max_abs_diff <- function(a, b) max(abs(a - b))
max_rel_diff <- function(a, b) {
  denom <- pmax(abs(a), abs(b), 1e-10)
  max(abs(a - b) / denom)
}

tol <- 5e-4

# --- Test 1: Basic correctness (1 thread vs 4 threads) ---
cat("Test 1: Basic correctness (50x100, Nk=5)\n")
res1_st <- run_peer(50, 100, 5, 50, threads = 1)
res1_mt <- run_peer(50, 100, 5, 50, threads = 4)

check("X close", max_abs_diff(res1_st$X, res1_mt$X) < tol)
check("W close", max_abs_diff(res1_st$W, res1_mt$W) < tol)
check("Alpha close", max_rel_diff(res1_st$Alpha, res1_mt$Alpha) < tol)
check("Eps close", max_rel_diff(res1_st$Eps, res1_mt$Eps) < tol)
check("Residuals close", max_abs_diff(res1_st$residuals, res1_mt$residuals) < tol)

# --- Test 2: With covariates ---
cat("\nTest 2: With covariates (50x100, Nk=5, 3 covariates)\n")
res2_st <- run_peer(50, 100, 5, 50, threads = 1, covs_cols = 3)
res2_mt <- run_peer(50, 100, 5, 50, threads = 4, covs_cols = 3)

check("X close (covs)", max_abs_diff(res2_st$X, res2_mt$X) < tol)
check("W close (covs)", max_abs_diff(res2_st$W, res2_mt$W) < tol)
check("Residuals close (covs)", max_abs_diff(res2_st$residuals, res2_mt$residuals) < tol)

# --- Test 3: Larger problem ---
cat("\nTest 3: Larger problem (200x500, Nk=15)\n")
res3_st <- run_peer(200, 500, 15, 30, threads = 1)
res3_mt <- run_peer(200, 500, 15, 30, threads = 4)

check("X close (large)", max_abs_diff(res3_st$X, res3_mt$X) < tol)
check("W close (large)", max_abs_diff(res3_st$W, res3_mt$W) < tol)
check("Residuals close (large)", max_abs_diff(res3_st$residuals, res3_mt$residuals) < tol)

# --- Test 4: Determinism (single-thread, same seed = bitwise identical) ---
cat("\nTest 4: Determinism (single thread, same seed)\n")
res4a <- run_peer(50, 100, 5, 50, threads = 1, seed = 99)
res4b <- run_peer(50, 100, 5, 50, threads = 1, seed = 99)

check("X identical (1 thread)", identical(res4a$X, res4b$X))
check("W identical (1 thread)", identical(res4a$W, res4b$W))
check("Residuals identical (1 thread)", identical(res4a$residuals, res4b$residuals))

# --- Test 5: Convergence consistency ---
cat("\nTest 5: Convergence consistency\n")
check("Similar iteration count", abs(length(res1_st$bounds) - length(res1_mt$bounds)) <= 2)
if (length(res1_st$bounds) > 0 && length(res1_mt$bounds) > 0) {
  check("Final bounds close", abs(tail(res1_st$bounds, 1) - tail(res1_mt$bounds, 1)) < 1)
} else {
  check("Final bounds close", TRUE)
}

# --- Test 6: Output shapes ---
cat("\nTest 6: Output shapes\n")
check("X shape", all(dim(res1_st$X) == c(50, 5)))
check("W shape", all(dim(res1_st$W) == c(100, 5)))
check("Alpha shape", all(dim(res1_st$Alpha) == c(5, 1)))
check("Eps shape", all(dim(res1_st$Eps) == c(100, 1)))
check("Residuals shape", all(dim(res1_st$residuals) == c(50, 100)))

# --- Test 7: Convergence ---
cat("\nTest 7: Convergence\n")
check("Bounds increase (ST)", length(res1_st$bounds) <= 1 || all(diff(res1_st$bounds) >= -1e-2))
check("Bounds increase (MT)", length(res1_mt$bounds) <= 1 || all(diff(res1_mt$bounds) >= -1e-2))

# --- Test 8: PEER_setNThreads / PEER_getNThreads ---
cat("\nTest 8: Thread control API\n")
m <- PEER()
PEER_setNThreads(m, 3)
check("getNThreads returns 3", PEER_getNThreads(m) == 3)
PEER_setNThreads(m, 1)
check("getNThreads returns 1", PEER_getNThreads(m) == 1)

# --- Summary ---
cat(sprintf("\n=== Results: %d passed, %d failed ===\n", pass_count, fail_count))
if (fail_count > 0) quit(status = 1)
