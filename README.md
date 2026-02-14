# peer — Patched for Modern Compilers

R package for **PEER** (Probabilistic Estimation of Expression Residuals), a Bayesian factor analysis method for inferring hidden determinants from gene expression data.

This is a patched fork of the [original PEER R package](https://github.com/PMBio/peer) (v1.0, Stegle et al. 2011) that compiles with modern C++ compilers (g++ 11+, C++17) without requiring legacy toolchains.

## What was fixed

The original source bundles Eigen 2.92 (~2010), which fails to compile under C++11+ due to:

| Issue | File | Fix |
|-------|------|-----|
| `EIGEN_ASM_COMMENT` macro uses `"#"X` (invalid user-defined literal in C++11+) | `src/Eigen/src/Core/util/Macros.h` | Changed to `"#" #X` (proper token pasting) |
| `std::binder1st/2nd`, `register` keyword (removed in C++17) | Various Eigen headers | Suppressed via `-Wno-deprecated -Wno-register` in `Makevars` |
| Pre-compiled `libpeer.so` in source tree | `src/` | Removed |
| Deprecated `.First.lib` loader | `R/firstlib.R` | Modernized to `.onLoad` |

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

## Citation

Stegle O, Parts L, Piipari M, Winn J, Durbin R. *Using probabilistic estimation of expression residuals (PEER) to obtain increased power and interpretability of gene expression analyses.* Nature Protocols 7, 500-507 (2012).

## License

CC BY-NC-SA 3.0 (original license from the PEER authors).
