genData <- function(
    n_samples = 40,
    true_rank = 4,
    n_features = 10,
    domain_points = 30,
    noise_sd = 0.2,
    missing_rate = 0.5,
    signal_means = NULL,
    signal_sd = 0.02
) {
  if (true_rank < 1L || n_samples < 1L || n_features < 1L || domain_points < 2L) {
    stop("n_samples, true_rank, n_features, and domain_points must be positive.")
  }
  if (n_samples < true_rank) {
    stop("n_samples must be at least as large as true_rank.")
  }
  if (n_features < true_rank) {
    stop("n_features must be at least as large as true_rank so each rank has at least one active feature.")
  }
  if (!is.null(signal_means) && length(signal_means) != true_rank) {
    stop("signal_means must have length equal to true_rank.")
  }
  if (missing_rate < 0 || missing_rate > 1) {
    stop("missing_rate must be between 0 and 1.")
  }
  
  if (is.null(signal_means)) {
    signal_means <- seq(from = true_rank + 3, by = -1, length.out = true_rank)
  }
  
  time_grid <- seq(0, 1, length.out = domain_points)
  
  make_shape_library <- function(d) {
    t_full <- seq(0, 1, length.out = d)
    left_len <- floor(d / 2)
    right_len <- d - left_len
    
    t_left <- c(seq(0, 1, length.out = left_len), rep(0, right_len))
    t_right <- c(rep(0, left_len), seq(0, 1, length.out = right_len))
    
    list(
      sin_left = sin(pi * t_left),
      sin_full = sin(pi * t_full),
      sin_right = sin(2 * pi * t_right)
    )
  }
  
  allocate_counts <- function(total, weights, min_count = 1L) {
    weights <- as.numeric(weights)
    scaled <- total * weights / sum(weights)
    counts <- pmax(min_count, floor(scaled))
    excess <- sum(counts) - total
    
    if (excess > 0) {
      frac <- scaled - floor(scaled)
      order_idx <- order(frac, decreasing = FALSE)
      for (idx in order_idx) {
        removable <- counts[idx] - min_count
        if (removable <= 0) {
          next
        }
        step <- min(removable, excess)
        counts[idx] <- counts[idx] - step
        excess <- excess - step
        if (excess == 0) {
          break
        }
      }
    } else if (excess < 0) {
      deficit <- -excess
      frac <- scaled - floor(scaled)
      order_idx <- order(frac, decreasing = TRUE)
      pos <- 1L
      while (deficit > 0) {
        idx <- order_idx[pos]
        counts[idx] <- counts[idx] + 1L
        deficit <- deficit - 1L
        pos <- pos + 1L
        if (pos > length(order_idx)) {
          pos <- 1L
        }
      }
    }
    
    counts
  }
  
  matrix_to_nested_list <- function(X_mat, n_feat, d) {
    feature_blocks <- lapply(seq_len(n_feat), function(j) {
      cols <- ((j - 1) * d + 1):(j * d)
      X_mat[, cols, drop = FALSE]
    })
    
    lapply(seq_len(nrow(X_mat)), function(i) {
      lapply(feature_blocks, function(block) block[i, ])
    })
  }
  
  add_missingness <- function(X, rate) {
    lapply(X, function(subject_curves) {
      lapply(subject_curves, function(curve) {
        miss_size <- round(length(curve) * rate)
        if (miss_size > 0) {
          miss_id <- sample(seq_along(curve), size = miss_size)
          curve[miss_id] <- NA
        }
        curve
      })
    })
  }
  
  shape_lib <- make_shape_library(domain_points)
  shape_names <- names(shape_lib)
  
  sample_group_sizes <- rep(n_samples %/% true_rank, true_rank)
  sample_group_sizes[seq_len(n_samples %% true_rank)] <-
    sample_group_sizes[seq_len(n_samples %% true_rank)] + 1L
  sample_end <- cumsum(sample_group_sizes)
  sample_start <- c(1L, head(sample_end, -1L) + 1L)
  subject_groups <- Map(seq.int, sample_start, sample_end)
  
  feature_weights <- rep(c(3, 3, 2, 2), length.out = true_rank)
  feature_group_sizes <- allocate_counts(n_features, feature_weights, min_count = 1L)
  feature_end <- cumsum(feature_group_sizes)
  feature_start <- c(1L, head(feature_end, -1L) + 1L)
  feature_groups <- Map(seq.int, feature_start, feature_end)
  
  U <- matrix(0, nrow = n_samples, ncol = true_rank)
  for (k in seq_len(true_rank)) {
    U[subject_groups[[k]], k] <- stats::rnorm(
      length(subject_groups[[k]]),
      mean = 0.2,
      sd = signal_sd
    )
  }
  
  V <- matrix(0, nrow = true_rank, ncol = n_features * domain_points)
  for (k in seq_len(true_rank)) {
    active_features <- feature_groups[[k]]
    shape_ids <- rep(seq_along(shape_names), length.out = length(active_features))
    for (i in seq_along(active_features)) {
      feature_id <- active_features[i]
      cols <- ((feature_id - 1) * domain_points + 1):(feature_id * domain_points)
      V[k, cols] <- shape_lib[[shape_names[shape_ids[i]]]]
    }
  }
  
  D <- diag(signal_means, nrow = true_rank)
  X_mat <- U %*% D %*% V +
    matrix(stats::rnorm(n_samples * n_features * domain_points, sd = noise_sd),
           nrow = n_samples)
  
  X_list <- matrix_to_nested_list(X_mat, n_feat = n_features, d = domain_points)
  X_list_na <- add_missingness(X_list, missing_rate)
  
  list(
    X_list = X_list,
    X_list_na = X_list_na,
    X_mat = X_mat,
    U = U,
    V = V,
    D = D,
    time_grid = time_grid,
    subject_groups = subject_groups,
    feature_groups = feature_groups,
    feature_group_sizes = feature_group_sizes,
    sample_group_sizes = sample_group_sizes,
    true_rank = true_rank,
    n_samples = n_samples,
    n_features = n_features,
    domain_points = domain_points,
    noise_sd = noise_sd,
    missing_rate = missing_rate
  )
}
