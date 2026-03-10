# peer — Patched & Enhanced

R package for **PEER** (Probabilistic Estimation of Expression Residuals), a Bayesian factor analysis method for inferring hidden determinants from gene expression data.

This is a patched fork of the [original PEER R package](https://github.com/PMBio/peer) (v1.3, Stegle et al. 2012) that compiles with modern C++ compilers (g++ 11+, C++17) and adds OpenMP parallelization, reproducible seeding, and other improvements. **The underlying inference algorithm is unchanged**, but results will not be bitwise identical to the original due to the RNG change (see below).

## What changed from the original

### Compiler fixes

The original source bundled Eigen 2.92 (~2010), which fails to compile under C++11+. This fork **upgrades to Eigen 3.4.0** (the latest stable release) and applies the following additional fixes:

| Issue | Fix |
|-------|-----|
| Bundled Eigen 2.92 incompatible with C++11+ | Upgraded to Eigen 3.4.0 |
| Pre-compiled `libpeer.so` in source tree | Removed |
| Deprecated `.First.lib` loader | Modernized to `.onLoad` |

### Enhancements

| Change | Impact |
|--------|--------|
| **OpenMP parallelization** of 5 hot loops | Multi-core speedup (3-10x at 4-16 threads) |
| **Reproducible RNG**: `randn()` uses R's `unif_rand()` instead of C's `rand()` | `set.seed()` now controls PEER initialization. **This changes results vs the original** (which used an independent, non-seedable C RNG) |
| **`PEER_setNThreads(model, n)` / `PEER_getNThreads(model)`** | Control thread count from R without env vars |
| **Ctrl+C support** via `R_CheckUserInterrupt()` | Long runs can be cancelled |
| **Improved verbose output** | Shows dimensions, threads, residual variance, timing per iteration |

## Install

```r
# Install directly from GitHub
remotes::install_github("stasaki/peer")

# Or from source
R CMD INSTALL peer
```

No special compiler flags, apt sources, or legacy g++ versions needed.

## Verify

```r
library(peer)
model <- PEER()
```

## Usage

```r
library(peer)

model <- PEER()

# Set expression data (samples x genes matrix)
PEER_setPhenoMean(model, as.matrix(expr))

# Set number of hidden factors
PEER_setNk(model, 20)

# Optionally add known covariates
PEER_setCovariates(model, as.matrix(covs))

# Run inference
PEER_update(model)

# Extract results
factors    <- PEER_getX(model)      # hidden factors
weights    <- PEER_getW(model)      # factor weights
precision  <- PEER_getAlpha(model)  # factor relevance
residuals  <- PEER_getResiduals(model)
```

## Parallelization (OpenMP)

PEER uses OpenMP to parallelize the dominant per-phenotype loops. On multi-core machines this gives significant speedup, especially with larger numbers of factors (Nk).

### Controlling threads

```r
# Check current thread count
PEER_getNThreads(model)

# Set thread count (takes effect immediately, no restart needed)
PEER_setNThreads(model, 4)
```

Alternatively, set the `OMP_NUM_THREADS` environment variable before starting R:

```bash
OMP_NUM_THREADS=4 Rscript my_script.R
```

### Reproducibility

PEER's random initialization is controlled by R's `set.seed()`, so results are fully reproducible:

```r
set.seed(42)
PEER_update(model)  # same result every time (with same thread count)
```

With a single thread, results are bitwise identical across runs. With multiple threads, results may differ at the ~1e-10 level due to floating-point accumulation order, but are functionally equivalent.

### Typical speedup

Benchmarks on 200 samples x 2000 genes, 20 iterations:

| Nk | 1 thread | 4 threads | 8 threads | 16 threads |
|----|----------|-----------|-----------|------------|
| 20 | 1.6s     | 0.6s (2.9x) | 0.4s (4.1x) | 0.3s (4.9x) |
| 50 | 9.4s     | 2.6s (3.6x) | 1.5s (6.4x) | 0.9s (10.6x) |

Larger Nk (more factors) scales better because per-iteration K³ matrix inversions dominate.

### Verbose output

PEER reports progress during `PEER_update()` at two verbosity levels:

```r
# Level 1 (default): iteration count, residual variance, timing
#   PEER update: Nj=200, Np=2000, Nk=20, Nc=0, threads=4
#   iteration 0/20 | var(resid)=1.0012 | 0.03s/iter | 0.0s elapsed
#   Converged (bound) after 15 iterations in 0.5s

# Level 2: adds bound value and convergence deltas
peer:::setVerbose(2)
```

## Authors

- **Original PEER**: Oliver Stegle, Matias Piipari, Leopold Parts
- **This fork (v1.4.0)**: Shinya Tasaki (modern compiler fixes, OpenMP, enhancements)

## Citation

Stegle O, Parts L, Piipari M, Winn J, Durbin R. *Using probabilistic estimation of expression residuals (PEER) to obtain increased power and interpretability of gene expression analyses.* Nature Protocols 7, 500-507 (2012).

## License

CC BY-NC-SA 3.0 (original license from the PEER authors).
