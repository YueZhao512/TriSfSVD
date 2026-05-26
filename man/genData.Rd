\name{genData}
\alias{genData}
\title{Generate Sparse Multivariate Functional Data}
\description{
Generates sparse multivariate functional data with known non-overlapping
sample clusters and feature clusters for illustrating and evaluating
\code{\link{Trisfsvd}}.
}
\usage{
genData(
  n_samples = 40,
  true_rank = 4,
  n_features = 10,
  domain_points = 30,
  noise_sd = 0.2,
  missing_rate = 0.5,
  signal_means = NULL,
  signal_sd = 0.02
)
}
\arguments{
  \item{n_samples}{Number of samples.}
  \item{true_rank}{Number of latent sparse functional components used to generate the data.}
  \item{n_features}{Number of functional variables.}
  \item{domain_points}{Number of observed domain points per functional variable.}
  \item{noise_sd}{Standard deviation of the additive Gaussian noise.}
  \item{missing_rate}{Proportion of within-curve observations replaced by \code{NA}.}
  \item{signal_means}{
    Optional vector controlling the component-specific signal strengths. If
    omitted, a decreasing default sequence is used.
  }
  \item{signal_sd}{Standard deviation used to generate nonzero sample scores.}
}
\value{
A list with components:
\item{X_list}{Complete nested-list data without missing values.}
\item{X_list_na}{Nested-list data with missing values inserted at random.}
\item{X_mat}{Underlying matrix representation of the simulated data.}
\item{U}{Matrix of true sample scores.}
\item{V}{Matrix of true functional loadings.}
\item{D}{Diagonal matrix of component strengths.}
\item{time_grid}{Common domain grid on \eqn{[0,1]}.}
\item{subject_groups}{List of true non-overlapping sample clusters.}
\item{feature_groups}{List of true component-specific feature clusters.}
\item{feature_group_sizes}{Numbers of active features assigned to each component.}
\item{sample_group_sizes}{Numbers of active samples assigned to each component.}
}
\details{
The generator is designed to mirror the structure targeted by
\code{\link{Trisfsvd}}. Samples are divided into non-overlapping latent groups,
and functional variables are assigned to component-specific active sets. The
output is returned in the nested-list format expected by the fitting and
post-processing functions in the package.
}
\examples{
set.seed(2)

sim <- genData()

str(sim$X_list_na, max.level = 2)
sim$feature_group_sizes
sim$sample_group_sizes

sim2 <- genData(
  n_samples = 60,
  true_rank = 4,
  n_features = 20,
  noise_sd = 0.3
)

sim2$feature_group_sizes
}
