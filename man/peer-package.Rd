\name{peer-package}
\alias{peer-package}
\alias{peer}
\alias{PEER}
\alias{PEER_getNc}
\alias{PEER_getPriorAlphaA}
\alias{PEER_getSparsityPrior}
\alias{PEER_init_net_sparse}
\alias{PEER_setPhenoVar}
\alias{PEER_setSigmaOff}
\alias{PEER_calcBound}
\alias{PEER_getNj}
\alias{PEER_getPriorAlphaB}
\alias{PEER_getTolerance}
\alias{PEER_logprob}
\alias{PEER_setPriorAlpha}
\alias{PEER_setSparsityPrior}
\alias{PEER_getAdd_mean}
\alias{PEER_getNk}
\alias{PEER_getPriorEpsA}
\alias{PEER_getVarTolerance}
\alias{PEER_setAdd_mean}
\alias{PEER_setPriorAlphaA}
\alias{PEER_setTolerance}
\alias{PEER_getAlpha}
\alias{PEER_getNmax_iterations}
\alias{PEER_getPriorEpsB}
\alias{PEER_getW}
\alias{PEER_setCovariates}
\alias{PEER_setPriorAlphaB}
\alias{PEER_setVarTolerance}
\alias{PEER_getBounds}
\alias{PEER_getNp}
\alias{PEER_getResidualVars}
\alias{PEER_getX}
\alias{PEER_setNk}
\alias{PEER_setPriorEps}
\alias{PEER_update}
\alias{PEER_getCovariates}
\alias{PEER_getPhenoMean}
\alias{PEER_getResiduals}
\alias{PEER_getZ}
\alias{PEER_setNmax_iterations}
\alias{PEER_setPriorEpsA}
%% peer-\alias{package}
\alias{PEER_getEps}
\alias{PEER_getPhenoVar}
\alias{PEER_getSigmaOff}
\alias{PEER_init_net}
\alias{PEER_setPhenoMean}
\alias{PEER_setNThreads}
\alias{PEER_getNThreads}
\alias{PEER_plotModel}
\alias{PEER_plotQtl}
\docType{package}
\title{
PEER: Probabilistic Estimation of Expression Residuals
}
\description{

PEER is a collection of Bayesian approaches to infer hidden determinants and their effects from gene expression profiles using factor analysis methods. Applications of PEER have detected batch effects and experimental confounders
increased the number of expression QTL findings by threefold allowed inference of intermediate cellular traits, such as transcription factor or pathway activations. 

The PEER model, inference, and applications are described in ``A Bayesian framework to account for complex non-genetic factors in gene expression levels greatly increases power in eQTL studies.'' (PLoS Computational Biology, May 2010)
Joint Genetic Analysis of Gene Expression Data with Inferred Cellular Phenotypes. (PLoS Genetics, January 2011)
and several other projects have successfully used PEER.

This project offers an efficient and versatile C++ implementation of the underlying algorithms with user-friendly interfaces to R and python. The core loops are parallelized with OpenMP for multi-core speedup. To get started using PEER, download the source or binary versions, see the installation instructions, and take a look at the getting started tutorial.
}
\section{Differences from the original PEER package}{
This is a patched fork of the original PEER package (v1.3, Stegle et al.)
with the following changes. The underlying inference algorithm is unchanged,
but results will not be bitwise identical to the original due to the RNG
change below.

\describe{
  \item{Modern compiler support}{Fixed Eigen 2.92 compatibility issues with
    g++ 11+ and C++17 (macro syntax, deprecated keywords).}
  \item{OpenMP parallelization}{The 5 dominant per-phenotype/per-sample loops
    are parallelized with OpenMP. New functions \code{PEER_setNThreads(n)}
    and \code{PEER_getNThreads()} control thread count at runtime.}
  \item{Reproducible RNG}{The internal random number generator now uses R's
    RNG (\code{unif_rand()}) instead of C's \code{rand()}, so
    \code{set.seed()} controls PEER's random initialization. The original
    used an independent C RNG that was not seedable from R.}
  \item{Interruptible}{Long-running \code{PEER_update()} calls can be
    cancelled with Ctrl+C (via \code{R_CheckUserInterrupt()} between
    iterations).}
  \item{Improved logging}{Verbose output now shows data dimensions, thread
    count, residual variance, per-iteration timing, and total wall time.}
}
}
\details{
\tabular{ll}{
Package: \tab peer\cr
Type: \tab Package\cr
Version: \tab 1.4.0\cr
Date: \tab 2026-03-10\cr
License: \tab CC BY-NC-SA 3.0\cr
LazyLoad: \tab yes\cr
}
}
\author{
Oliver Stegle, Matias Piipari, Leopold Parts (original PEER)

Shinya Tasaki (modern compiler fixes, OpenMP, enhancements)

Maintainer: Shinya Tasaki <Shinya_Tasaki@rush.edu>
}
\references{

A Bayesian framework to account for complex non-genetic factors in gene expression levels greatly increases power in eQTL studies. (PLoS Computational Biology, May 2010)

Joint Genetic Analysis of Gene Expression Data with Inferred Cellular Phenotypes. (PLoS Genetics, January 2011)

}
\keyword{ package }
\section{Parallelization (OpenMP)}{
PEER uses OpenMP to parallelize the dominant per-phenotype loops in the
variational inference. On multi-core machines this gives significant speedup,
especially with larger numbers of factors (Nk).

Use \code{PEER_setNThreads(n)} to set the number of threads within an R
session, or set the \code{OMP_NUM_THREADS} environment variable before
starting R. Use \code{PEER_getNThreads()} to query the current setting.

PEER's random initialization is controlled by R's \code{set.seed()}, so
results are reproducible. With a single thread, results are bitwise
identical across runs. With multiple threads, results may differ at the
~1e-10 level due to floating-point accumulation order.
}
\examples{
model <- PEER()

# Set expression data (samples x genes)
PEER_setPhenoMean(model, matrix(rnorm(50 * 100), 50, 100))

# Set number of hidden factors
PEER_setNk(model, 10)

# Control parallelization
PEER_setNThreads(4)
PEER_getNThreads()

# Run inference (Ctrl+C to interrupt)
PEER_update(model)

# Extract results
factors   <- PEER_getX(model)
weights   <- PEER_getW(model)
precision <- PEER_getAlpha(model)
residuals <- PEER_getResiduals(model)
}
