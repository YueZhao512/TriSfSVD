\name{plotResult}
\alias{plotResult}
\title{Plot Sample-Feature-Time Patterns for Paired Clusters}
\description{
Plots observed and reconstructed functional patterns for paired sample and
feature clusters obtained either directly from \code{\link{Trisfsvd}} or from
\code{\link{getClusterAssoc}}.
}
\usage{
plotResult(
  fit, X, sample_cluster, feature_cluster, association_indices = NULL,
  mean = TRUE, use_s = TRUE, show_variables = TRUE, random_seed = NULL,
  ranks = NULL, nrow_page = 4, ncol_page = 5, colours = NULL
)
}
\arguments{
  \item{fit}{A fitted object returned by \code{\link{Trisfsvd}}.}
  \item{X}{
    The original nested-list functional data. Each sample is a list of
    functional variables, and each variable is a numeric domain vector.
  }
  \item{sample_cluster}{
    A list of sample-index vectors, such as \code{fit$sample_cluster} or the
    \code{sample_cluster} output from \code{\link{getClusterAssoc}}.
  }
  \item{feature_cluster}{
    A list of feature-index vectors, such as \code{fit$feature_cluster} or the
    \code{feature_cluster} output from \code{\link{getClusterAssoc}}.
  }
  \item{association_indices}{
    Integer vector indicating which paired associations to plot. By default,
    all available pairs up to
    \code{min(length(sample_cluster), length(feature_cluster))} are shown.
  }
  \item{mean}{
    Logical indicating whether to plot cluster mean curves. If \code{TRUE},
    the mean observed curve and mean reconstructed curve are shown. If
    \code{FALSE}, all observed sample curves and all reconstructed sample
    curves are shown.
  }
  \item{use_s}{
    Logical indicating whether to include the fitted scalar values \code{s} in
    the reconstructed functional curves.
  }
  \item{show_variables}{
    Controls how many variables are displayed for each association. If
    \code{TRUE}, all variables in the feature cluster are shown. If a single
    integer is supplied, that many variables are randomly selected from the
    feature cluster.
  }
  \item{random_seed}{
    Optional random seed used when \code{show_variables} is a single integer.
  }
  \item{ranks}{
    Integer vector of fitted ranks to use in the reconstruction. By default,
    all fitted components in \code{fit$results_list} are used.
  }
  \item{nrow_page}{Number of panel rows per plotting page.}
  \item{ncol_page}{Number of panel columns per plotting page.}
  \item{colours}{Optional vector of colours used for the paired associations.}
}
\value{
Invisibly returns a list giving the feature indices plotted for each
association.
}
\details{
For each selected functional variable, the function reconstructs the fitted
sample-feature-time pattern by summing the chosen sparse functional components.
It then subsets the requested sample cluster and overlays the observed
functional trajectories with the corresponding reconstructed trajectories. This
provides a direct visual summary of how the fitted sparse functional
decomposition represents the paired sample and feature clusters.
}
\examples{
\dontrun{
set.seed(2)

sim <- genData(
  n_samples = 20,
  true_rank = 2,
  n_features = 6,
  domain_points = 20,
  noise_sd = 0.15,
  missing_rate = 0.3
)

fit <- Trisfsvd(
  X = sim$X_list_na,
  k = 2,
  gamma_candidates = c(0, 0.05),
  theta_candidates = c(0, 0.05),
  lambda_candidates = c(0, 0.1),
  alpha_candidates = c(0, 1),
  max_outer_iter = 2,
  tol_outer = 5e-3,
  verbose = FALSE,
  use_sgl = TRUE,
  max_iter = 500,
  eps = 1e-3,
  gamma_adaptive_power = 0.5,
  theta_adaptive_power = 0,
  lambda_adaptive_power = 0,
  extend_weight_domain = 0.5,
  extend_weight_smooth = 0,
  extend_weight_group = 0,
  extend_weight_u = 0.5,
  parallel = FALSE,
  backtracking = TRUE,
  subj_overlap = TRUE,
  var_overlap = TRUE,
  subreg_overlap = TRUE
)

plotResult(
  fit = fit,
  X = sim$X_list_na,
  sample_cluster = fit$sample_cluster,
  feature_cluster = fit$feature_cluster,
  mean = TRUE,
  show_variables = 3
)
}
}
