\name{Trisfsvd}
\alias{Trisfsvd}
\title{Estimate Sparse Functional Components and Tri-Clusters}
\description{
Estimates sparse functional singular value decomposition components for
data-driven identification of biclustering and triclustering structure in
high-dimensional sparse multivariate functional data.
}
\usage{
Trisfsvd(
  X, k = 1, gamma_candidates, alpha_candidates, theta_candidates,
  lambda_candidates = NULL, max_outer_iter = 20, tol_outer = 1e-5,
  use_sgl = TRUE, verbose = TRUE, L_init = NULL, max_iter = 30000,
  eps = 1e-5, gamma_adaptive_power = 1, theta_adaptive_power = 1,
  lambda_adaptive_power = 1, extend_weight_domain, extend_weight_smooth,
  extend_weight_group, subj_overlap = FALSE, var_overlap = FALSE,
  subreg_overlap = FALSE, parallel = FALSE, n_cores = detectCores() - 1,
  extend_weight_u, backtracking
)
}
\arguments{
  \item{X}{
    A nested list of sparse multivariate functional observations. The outer
    list indexes samples. Each sample is a list of functional variables, and
    each variable is a numeric vector of domain measurements. Missing values
    should be represented with \code{NA}.
  }
  \item{k}{Number of rank-1 components to extract.}
  \item{gamma_candidates}{Candidate tuning values for sample-level sparsity.}
  \item{alpha_candidates}{Candidate tuning values for smoothness.}
  \item{theta_candidates}{Candidate tuning values for variable-level sparsity.}
  \item{lambda_candidates}{Candidate tuning values for domain-level sparsity.}
  \item{max_outer_iter}{Maximum number of outer iterations.}
  \item{tol_outer}{
    Convergence tolerance for the iterative algorithm used to estimate each
    rank-one regularized SVD component.
  }
  \item{use_sgl}{
    Logical indicating whether to use the sparse-group-lasso style update. If
    \code{TRUE}, sparse-group-lasso is used to obtain tri-cluster results. If
    \code{FALSE}, group-lasso is used to obtain bicluster results.
  }
  \item{verbose}{Logical indicating whether to print progress.}
  \item{L_init}{Initial step size parameter for the FISTA backtracking procedure.}
  \item{max_iter}{Maximum number of iterations for the sparse-group-lasso or group-lasso updates.}
  \item{eps}{
    Convergence tolerance for the FISTA iterations. The algorithm stops when
    the change between two successive iterates is smaller than \code{eps}.
  }
  \item{gamma_adaptive_power}{Adaptive weighting power for \code{gamma_candidates}.}
  \item{theta_adaptive_power}{Adaptive weighting power for variable-level sparsity tuning.}
  \item{lambda_adaptive_power}{Adaptive weighting power for domain-level sparsity tuning.}
  \item{extend_weight_domain}{Weights in the EBIC criterion for domain-level sparsity tuning selection.}
  \item{extend_weight_smooth}{Weights in the EBIC criterion for smoothing tuning selection.}
  \item{extend_weight_group}{Weights in the EBIC criterion for variable-level sparsity tuning selection.}
  \item{subj_overlap}{Logical indicating whether samples may overlap across estimated tri-clusters or biclusters.}
  \item{var_overlap}{Logical indicating whether variables may overlap across estimated tri-clusters or biclusters.}
  \item{subreg_overlap}{Logical indicating whether domain-level regions may overlap across estimated tri-clusters or biclusters.}
  \item{parallel}{Logical indicating whether to use parallel computation.}
  \item{n_cores}{Number of cores used when \code{parallel = TRUE}.}
  \item{extend_weight_u}{Weights in the EBIC criterion for sample-level sparsity tuning selection.}
  \item{backtracking}{
    Logical indicating whether to use backtracking in the FISTA algorithm. If
    \code{TRUE}, the backtracking FISTA algorithm is used. If \code{FALSE},
    FISTA is run without backtracking.
  }
}
\value{
A list with components:
\item{results_list}{
  A list of fitted rank-1 components. Each component contains estimated sample
  scores, variable/domain loading objects, selected tuning parameters, and
  related model-selection output.
}
\item{sample_cluster}{
  A list of sample-index vectors. Each element contains the sample indices
  whose fitted \code{u} values are nonzero for one extracted component.
}
\item{feature_cluster}{
  A list of feature-index vectors. Each element contains the feature indices
  whose fitted \code{varphi} functions are not identically zero for one
  extracted component.
}
\item{variance_explained}{
  Numeric vector giving the cumulative proportion of variation explained by the
  extracted components.
}
\item{BIC_score}{
  Numeric vector of BIC scores for the extracted components.
}
}
\details{
The function takes sparse multivariate functional observations stored as a
nested list and extracts rank-one sparse functional components. For each
component, the fitted sample scores, selected variables, and selected
time-domain subregions define an interpretable bicluster or tricluster
pattern, depending on the chosen penalty structure.

Typical use consists of three steps: prepare the nested-list input, fit
\code{Trisfsvd()}, and then inspect the returned sample clusters, feature
clusters, loading functions, and reconstruction summaries. When
\code{parallel = TRUE}, tuning-parameter path evaluations may use multiple
cores through the \pkg{parallel} package.
}
\examples{
\dontrun{
## ------------------------------------------------------------
## Example: generate data, then fit TriSfSVD
## ------------------------------------------------------------
## The package provides a simulation helper that generates a sparse
## multivariate functional data set in the nested-list format needed
## by Trisfsvd.

set.seed(2)

## Generate a synthetic data set using the package default settings.
## Here this means:
##   - 40 samples
##   - 4 latent rank-one components
##   - 10 functional features
##   - 30 domain points per feature
##   - Gaussian noise with standard deviation 0.2
##   - 50% missingness within each functional curve
sim <- genData(
  n_samples = 40,
  true_rank = 4,
  n_features = 10,
  domain_points = 30,
  noise_sd = 0.2,
  missing_rate = 0.5,
  signal_means = NULL,
  signal_sd = 0.02
)

## The returned object contains:
##   - sim$X_list: complete data without missing values
##   - sim$X_list_na: nested-list data with missing values
##   - sim$U, sim$V, sim$D: underlying latent signal objects
##   - sim$sample_group_sizes and sim$feature_group_sizes:
##     the true non-overlapping sample and feature allocations

## Candidate tuning grids. Larger grids may improve tuning resolution
## but increase computation time.
gamma_grid <- c(0, 0.05, 0.5)
theta_grid <- c(0, 0.01, 0.1)
lambda_grid <- c(0, 0.1, 1)
alpha_grid <- c(0, 1, 10)

## Fit the model. Setting use_sgl = TRUE uses sparse-group-lasso
## updates to estimate tri-cluster structure.
res <- Trisfsvd(
  X = sim$X_list_na,
  k = 2,
  gamma_candidates = gamma_grid,
  theta_candidates = theta_grid,
  lambda_candidates = lambda_grid,
  alpha_candidates = alpha_grid,
  max_outer_iter = 3,
  tol_outer = 2e-3,
  verbose = TRUE,
  use_sgl = TRUE,
  max_iter = 2000,
  eps = 1e-4,
  theta_adaptive_power = 0,
  lambda_adaptive_power = 0,
  gamma_adaptive_power = 0.5,
  extend_weight_domain = 0.5,
  extend_weight_u = 0.5,
  extend_weight_smooth = 0,
  extend_weight_group = 0,
  subj_overlap = TRUE,
  var_overlap = TRUE,
  subreg_overlap = TRUE,
  parallel = FALSE,
  backtracking = TRUE
)

res$sample_cluster
res$feature_cluster
}
}
