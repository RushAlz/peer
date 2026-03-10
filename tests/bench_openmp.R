#!/usr/bin/env Rscript
# Performance scaling benchmark for OpenMP parallelization of PEER
# Uses PEER_setNThreads() to control thread count within a single process

library(peer)

cat("=== PEER OpenMP Performance Benchmark ===\n\n")

m_tmp <- PEER()
max_threads <- PEER_getNThreads(m_tmp)
thread_counts <- unique(c(1, 2, 4, 8, min(16, max_threads)))
thread_counts <- thread_counts[thread_counts <= max_threads]

cat(sprintf("Available cores: %d\n", max_threads))
cat(sprintf("Thread counts: %s\n\n", paste(thread_counts, collapse = ", ")))

run_bench <- function(N, P, Nk, Niter, threads, reps = 3) {
  times <- numeric(reps)
  for (r in seq_len(reps)) {
    set.seed(42)
    pheno <- matrix(rnorm(N * P), N, P)
    model <- PEER()
    PEER_setPhenoMean(model, pheno)
    PEER_setNk(model, Nk)
    PEER_setNmax_iterations(model, Niter)
    PEER_setNThreads(model, threads)
    times[r] <- system.time(PEER_update(model))["elapsed"]
  }
  median(times)
}

run_config <- function(N, P, Nk, Niter, label) {
  cat(sprintf("--- %s: %d samples x %d genes, Nk=%d, %d iterations ---\n", label, N, P, Nk, Niter))
  results <- data.frame(threads = integer(), time = numeric(), config = character(),
                        stringsAsFactors = FALSE)
  for (nt in thread_counts) {
    med_time <- run_bench(N, P, Nk, Niter, nt)
    cat(sprintf("  Threads=%d: %.2fs\n", nt, med_time))
    results <- rbind(results, data.frame(threads = nt, time = med_time, config = label))
  }
  results
}

res1 <- run_config(200, 2000, 20, 20, "Nk=20")
cat("\n")
res2 <- run_config(200, 2000, 50, 20, "Nk=50")

all_res <- rbind(res1, res2)

cat("\n| Config | Threads | Time (s) | Speedup | Efficiency |\n")
cat("|--------|---------|----------|---------|------------|\n")
for (cfg in unique(all_res$config)) {
  sub <- all_res[all_res$config == cfg, ]
  t1 <- sub$time[sub$threads == 1]
  for (j in seq_len(nrow(sub))) {
    speedup <- t1 / sub$time[j]
    efficiency <- speedup / sub$threads[j]
    cat(sprintf("| %s | %d | %.2f | %.2fx | %.0f%% |\n",
                cfg, sub$threads[j], sub$time[j], speedup, efficiency * 100))
  }
}

cat("\n--- Analysis ---\n")
for (cfg in unique(all_res$config)) {
  sub <- all_res[all_res$config == cfg, ]
  t1 <- sub$time[sub$threads == 1]
  if (any(sub$threads == 4)) {
    s4 <- t1 / sub$time[sub$threads == 4]
    if (s4 < 1.5) cat(sprintf("WARNING [%s]: Speedup at 4 threads = %.2fx (< 1.5x)\n", cfg, s4))
    else cat(sprintf("OK [%s]: Speedup at 4 threads = %.2fx\n", cfg, s4))
  }
}

tryCatch({
  write.csv(all_res, "tests/bench_results.csv", row.names = FALSE)
  cat("\nResults saved to tests/bench_results.csv\n")
}, error = function(e) cat("(Could not save CSV)\n"))
