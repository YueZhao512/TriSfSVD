# TriSfSVD

TriSfSVD implements sparse functional singular value decomposition methods for
identifying biclustering and triclustering structure in high-dimensional sparse
multivariate functional data.

The package estimates sparse functional components that describe coordinated
patterns across samples, variables, and time-domain subregions. It also includes
helpers for generating example data, refining sample-cluster associations, and
plotting observed versus reconstructed sample-feature-time patterns.

## Installation

TriSfSVD requires R 4.1.0 or later. It is not currently on CRAN.

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("YueZhao512/TriSfSVD")
```

Or install from a local source checkout:

```r
install.packages("/path/to/TriSfSVD", repos = NULL, type = "source")
```

## Dependencies

TriSfSVD imports the following R packages:

- `energy`
- `irlba`
- `Matrix`
- `parallel`
- `softImpute`
- `stats`
- `utils`

Install missing dependencies with:

```r
install.packages(c("energy", "irlba", "Matrix", "softImpute"))
```

## Main Functions

- `Trisfsvd()`: fit sparse functional SVD components and recover sample and
  feature clusters.
- `genData()`: generate synthetic sparse multivariate functional data in the
  nested-list format expected by `Trisfsvd()`.
- `getClusterAssoc()`: refine sample clusters and summarize their association
  with component-specific feature clusters.
- `plotResult()`: visualize observed and reconstructed functional patterns for
  paired sample and feature clusters.

## Data Format

Input data are stored as a nested list:

```r
X[[sample_index]][[feature_index]][domain_index]
```

The outer list indexes samples. Each sample contains a list of functional
variables. Each variable is a numeric vector observed over the functional
domain. Missing values should be represented as `NA`.

For example, `genData()` returns both complete and sparse data:

- `sim$X_list`: complete nested-list data without missing values.
- `sim$X_list_na`: nested-list data with randomly inserted `NA` values.
- `sim$U`, `sim$V`, `sim$D`: latent signal objects used to generate the data.
- `sim$subject_groups` and `sim$feature_groups`: true sample and feature groups.

## Basic Example

```r
library(TriSfSVD)

set.seed(2)

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

gamma_grid <- c(0, 0.05, 0.5)
theta_grid <- c(0, 0.01, 0.1)
lambda_grid <- c(0, 0.1, 1)
alpha_grid <- c(0, 1, 10)

fit <- Trisfsvd(
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

fit$sample_cluster
fit$feature_cluster
fit$variance_explained
fit$BIC_score
```

## Cluster Association and Plotting

After fitting the model, use `getClusterAssoc()` to refine sample clusters and
summarize their association with fitted feature clusters. Then use `plotResult()`
to inspect the fitted sample-feature-time patterns.

```r
assoc <- getClusterAssoc(
  fit = fit,
  X = sim$X_list_na,
  n_clusters = 2
)

assoc$sample_cluster
assoc$feature_cluster
assoc$dcov_stat_matrix

plotResult(
  fit = fit,
  X = sim$X_list_na,
  sample_cluster = assoc$sample_cluster,
  feature_cluster = assoc$feature_cluster,
  mean = TRUE,
  show_variables = 3
)
```

## Returned Objects

`Trisfsvd()` returns a list containing:

- `results_list`: fitted rank-1 components, including sample scores, functional
  loadings, selected tuning parameters, and model-selection summaries.
- `sample_cluster`: sample indices with nonzero fitted sample scores for each
  extracted component.
- `feature_cluster`: feature indices whose fitted loading functions are not
  identically zero for each extracted component.
- `variance_explained`: cumulative proportion of variation explained by the
  extracted components.
- `BIC_score`: BIC scores for the extracted components.

`getClusterAssoc()` returns refined sample clusters, fitted feature clusters,
and distance-covariance summaries. If `bootstrapping = TRUE`, it also returns
bootstrap confidence intervals and standard errors for distance correlations.

## Parallel Computation

Set `parallel = TRUE` in `Trisfsvd()` to evaluate tuning-parameter paths using
multiple cores:

```r
fit_parallel <- Trisfsvd(
  X = sim$X_list_na,
  k = 2,
  gamma_candidates = gamma_grid,
  theta_candidates = theta_grid,
  lambda_candidates = lambda_grid,
  alpha_candidates = alpha_grid,
  extend_weight_domain = 0.5,
  extend_weight_u = 0.5,
  extend_weight_smooth = 0,
  extend_weight_group = 0,
  subj_overlap = TRUE,
  var_overlap = TRUE,
  subreg_overlap = TRUE,
  parallel = TRUE,
  n_cores = 2,
  backtracking = TRUE
)
```

Choose `n_cores` based on the number of available CPU cores and the memory
required by your data set.

## Local Build and Check

From the package root:

```sh
R CMD build .
R CMD check TriSfSVD_0.1.0.tar.gz
```

From the parent directory:

```sh
R CMD build TriSfSVD
R CMD check TriSfSVD_0.1.0.tar.gz
```

## License

TriSfSVD is released under the MIT license.
