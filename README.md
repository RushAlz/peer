# peer — Fast Multicore PEER

**PEER** (Probabilistic Estimation of Expression Residuals) is a Bayesian factor analysis method for inferring hidden determinants from gene expression data (Stegle et al. 2012). This fork adds **OpenMP parallelization** for 3-10x speedup on multicore machines while keeping the inference algorithm unchanged.

## Performance

Benchmarks on 200 samples x 2000 genes, 20 iterations:

| Factors (Nk) | 1 thread | 4 threads | 8 threads | 16 threads |
|---------------|----------|-----------|-----------|------------|
| 20 | 1.6s | 0.6s (2.9x) | 0.4s (4.1x) | 0.3s (4.9x) |
| 50 | 9.4s | 2.6s (3.6x) | 1.5s (6.4x) | 0.9s (10.6x) |

Speedup scales with the number of factors because per-iteration K³ matrix operations dominate at higher Nk.

## Install

```r
remotes::install_github("stasaki/peer")
```

Works out of the box with modern C++ compilers (g++ 11+, clang 14+, C++17). No special flags or legacy toolchains needed.

## Quick start

```r
library(peer)

model <- PEER()
PEER_setPhenoMean(model, as.matrix(expr))  # samples x genes
PEER_setNk(model, 20)                      # number of hidden factors
PEER_setNThreads(model, 8)                 # use 8 cores

set.seed(42)                               # reproducible results
PEER_update(model)

# Extract results
residuals <- PEER_getResiduals(model)
factors   <- PEER_getX(model)
weights   <- PEER_getW(model)
precision <- PEER_getAlpha(model)
```

## Controlling threads

```r
PEER_setNThreads(model, 4)   # set thread count (immediate, no restart)
PEER_getNThreads(model)      # check current count
```

Or set globally before starting R:

```bash
OMP_NUM_THREADS=4 Rscript my_script.R
```

Five hot loops in the variational update are parallelized with OpenMP. The default thread count follows `OMP_NUM_THREADS` (or all available cores if unset).

## Reproducibility

Random initialization uses R's `set.seed()`, so results are fully reproducible:

- **1 thread**: bitwise identical across runs.
- **Multiple threads**: differences at ~1e-10 due to floating-point accumulation order — functionally equivalent.

## Progress reporting

`PEER_update()` prints progress by default:

```
PEER update: Nj=200, Np=2000, Nk=20, Nc=0, threads=8
iteration 0/20 | var(resid)=1.0012 | 0.03s/iter | 0.0s elapsed
Converged (bound) after 15 iterations in 0.5s
```

Long-running jobs can be interrupted with Ctrl+C.

## Changes from the original

This is a fork of the [original PEER R package](https://github.com/PMBio/peer) (v1.3). The inference algorithm is identical, but results differ from the original due to the RNG change (R's `set.seed()` replaces C's non-seedable `rand()`).

| Change | Details |
|--------|---------|
| Eigen 2.92 → 3.4.0 | Fixes all C++11/17 compilation failures |
| OpenMP parallelization | 5 hot loops parallelized for 3-10x multicore speedup |
| `PEER_setNThreads` / `PEER_getNThreads` | Per-model thread control from R |
| Reproducible RNG | `randn()` backed by R's `unif_rand()` + `set.seed()` |
| Ctrl+C support | `R_CheckUserInterrupt()` between iterations |
| Verbose output | Dimensions, thread count, residual variance, timing |

## Citation

Stegle O, Parts L, Piipari M, Winn J, Durbin R. *Using probabilistic estimation of expression residuals (PEER) to obtain increased power and interpretability of gene expression analyses.* Nature Protocols 7, 500-507 (2012).

## Authors

- **Original PEER**: Oliver Stegle, Matias Piipari, Leopold Parts
- **This fork**: Shinya Tasaki

## License

CC BY-NC-SA 3.0 (original license from the PEER authors).
