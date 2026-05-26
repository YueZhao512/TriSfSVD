\name{getClusterAssoc}
\alias{getClusterAssoc}
\title{Refine Sample Clusters and Summarize Cluster Associations}
\description{
Post-processes fitted \code{\link{Trisfsvd}} components to refine sample
clusters and quantify their association with component-specific feature
clusters and selected time-domain subregions.
}
\usage{
getClusterAssoc(
  fit, X, n_clusters, component_indices = NULL, scale_u = TRUE,
  km_iter_max = 10000, km_nstart = 1000, soft_rank = 1, soft_lambda = 0,
  dcor_R = 199, dcor_index = 1, bootstrapping = FALSE,
  bootstrap_B = 1000, bootstrap_seed = 2025, bootstrap_n_min = 5,
  keep_boot_array = FALSE,
  mask_align = c("strict", "truncate", "pad_false"),
  block_align = c("strict", "pad", "truncate"),
  block_pad_value = NA_real_
)
}
\arguments{
  \item{fit}{A fitted object returned by \code{\link{Trisfsvd}}.}
  \item{X}{
    The original nested-list functional data used for fitting. Missing values
    may be represented with \code{NA}.
  }
  \item{n_clusters}{Number of refined sample clusters to estimate by k-means.}
  \item{component_indices}{
    Integer vector specifying which fitted components to use. By default, the
    first \code{min(length(fit$results_list), n_clusters)} components are used.
  }
  \item{scale_u}{Logical indicating whether to standardize the selected sample scores before k-means clustering.}
  \item{km_iter_max}{Maximum number of iterations for k-means.}
  \item{km_nstart}{Number of random starts for k-means.}
  \item{soft_rank}{Rank used in the \pkg{softImpute} decomposition of the masked data matrices.}
  \item{soft_lambda}{Regularization parameter passed to \code{softImpute::softImpute()}.}
  \item{dcor_R}{Number of permutations used in \code{energy::dcor.test()}.}
  \item{dcor_index}{Index parameter passed to \code{energy::dcor.test()}.}
  \item{bootstrapping}{Logical indicating whether to compute bootstrap summaries for the distance correlation matrix.}
  \item{bootstrap_B}{Number of bootstrap resamples when \code{bootstrapping = TRUE}.}
  \item{bootstrap_seed}{Random seed used for the bootstrap resampling step.}
  \item{bootstrap_n_min}{Minimum number of finite paired observations required when computing bootstrap distance correlations.}
  \item{keep_boot_array}{Logical indicating whether to retain the full bootstrap array internally.}
  \item{mask_align}{Alignment rule used when masking selected subregions.}
  \item{block_align}{Alignment rule used when converting masked curves to a block matrix.}
  \item{block_pad_value}{Padding value used when \code{block_align = "pad"}.}
}
\value{
A list with components:
\item{sample_cluster}{
  A list of sample-index vectors defining the refined sample clusters.
}
\item{feature_cluster}{
  A list of feature-index vectors defining the component-specific feature
  clusters taken from the fitted \code{Trisfsvd} object.
}
\item{dcov_stat_matrix}{Distance covariance test statistic matrix.}
\item{dcor_hat}{Bootstrap target distance correlation matrix. Returned only when \code{bootstrapping = TRUE}.}
\item{ci_lo}{Lower 95\% bootstrap confidence bound matrix. Returned only when \code{bootstrapping = TRUE}.}
\item{ci_hi}{Upper 95\% bootstrap confidence bound matrix. Returned only when \code{bootstrapping = TRUE}.}
\item{se_hat}{Bootstrap standard error matrix. Returned only when \code{bootstrapping = TRUE}.}
}
\details{
The function is intended for settings where the fitted \code{u} vectors from
\code{\link{Trisfsvd}} are informative but not sufficiently distinct to define
final sample clusters on their own. It first forms a low-dimensional sample
representation by stacking selected \code{u} vectors, then applies k-means to
obtain refined sample clusters.

For each selected component, the fitted loading functions are converted into
feature-specific subregion masks. These masks are applied to the original
functional data, and a rank-reduced summary is obtained by \pkg{softImpute}.
Distance covariance and distance correlation summaries are then used to
quantify the association between the refined sample clusters and the
component-specific functional patterns.
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

assoc <- getClusterAssoc(
  fit = fit,
  X = sim$X_list_na,
  n_clusters = 2,
  component_indices = 1:2,
  bootstrapping = FALSE
)

assoc$sample_cluster
assoc$dcov_stat_matrix
}
}
