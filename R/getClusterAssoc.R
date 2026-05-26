mask_X_by_bool_list <- function(
    X,
    bool_list,
    zero_value = 0,
    na_as_false = TRUE,
    align = c("strict", "truncate", "pad_false"),
    return_index = FALSE
) {
  stopifnot(is.list(X), is.list(bool_list))
  align <- match.arg(align)
  
  n <- length(X)
  
  b_list <- lapply(bool_list, function(b) {
    if (!is.logical(b)) {
      b <- as.logical(b)
    }
    if (na_as_false) {
      b[is.na(b)] <- FALSE
    }
    b
  })
  
  keep <- vapply(b_list, function(b) any(b), logical(1))
  keep_idx <- which(keep)
  drop_idx <- which(!keep)
  
  feat_names_src <- names(bool_list)
  feat_names_keep <- if (!is.null(feat_names_src)) feat_names_src[keep_idx] else NULL
  
  mask_vec <- function(v, b) {
    if (length(b) != length(v)) {
      if (align == "truncate") {
        m <- min(length(b), length(v))
        b <- b[seq_len(m)]
        v <- v[seq_len(m)]
      } else if (align == "pad_false") {
        if (length(b) < length(v)) {
          b <- c(b, rep(FALSE, length(v) - length(b)))
        } else if (length(b) > length(v)) {
          v <- c(v, rep(zero_value, length(b) - length(v)))
        }
      }
    }
    out <- v
    out[!b] <- zero_value
    out
  }
  
  X_new <- vector("list", n)
  names(X_new) <- names(X)
  for (i in seq_len(n)) {
    Xi <- X[[i]]
    Xi_new <- lapply(keep_idx, function(j) mask_vec(Xi[[j]], b_list[[j]]))
    if (!is.null(feat_names_keep)) {
      names(Xi_new) <- feat_names_keep
    } else if (!is.null(names(Xi))) {
      names(Xi_new) <- names(Xi)[keep_idx]
    }
    X_new[[i]] <- Xi_new
  }
  
  if (return_index) {
    attr(X_new, "keep_idx") <- keep_idx
    attr(X_new, "drop_idx") <- drop_idx
  }
  
  X_new
}

as_block_matrix <- function(
    X_new,
    align = c("strict", "pad", "truncate"),
    pad_value = NA_real_
) {
  stopifnot(is.list(X_new), length(X_new) > 0, is.list(X_new[[1]]))
  align <- match.arg(align)
  
  n <- length(X_new)
  p <- length(X_new[[1]])
  t_each <- vapply(seq_len(p), function(j) {
    max(vapply(X_new, function(xi) length(xi[[j]]), integer(1)))
  }, integer(1))
  
  if (align == "strict" && length(unique(t_each)) != 1L) {
    stop("All retained features must have the same length when align = 'strict'.")
  }
  
  t_target <- switch(
    align,
    strict = unique(t_each)[1],
    pad = max(t_each),
    truncate = min(t_each)
  )
  
  M <- matrix(NA_real_, nrow = n, ncol = p * t_target)
  fn <- names(X_new[[1]])
  if (is.null(fn)) {
    fn <- paste0("feat", seq_len(p))
  }
  colnames(M) <- unlist(lapply(seq_len(p), function(j) {
    paste0(fn[j], "_t", seq_len(t_target))
  }))
  rownames(M) <- names(X_new)
  
  norm_len <- function(v) {
    if (length(v) == t_target) {
      return(v)
    }
    if (length(v) > t_target) {
      return(v[seq_len(t_target)])
    }
    c(v, rep(pad_value, t_target - length(v)))
  }
  
  for (i in seq_len(n)) {
    xi <- X_new[[i]]
    for (j in seq_len(p)) {
      cols <- ((j - 1L) * t_target + 1L):(j * t_target)
      M[i, cols] <- norm_len(xi[[j]])
    }
  }
  
  M
}

svd_with_na_softImpute <- function(M, rank = 1, lambda = 0) {
  fit <- softImpute::softImpute(M, rank.max = rank, lambda = lambda, type = "als")
  U <- fit$u
  V <- fit$v
  if (is.null(dim(U))) {
    U <- matrix(U, ncol = 1L)
  }
  if (is.null(dim(V))) {
    V <- matrix(V, ncol = 1L)
  }
  list(
    U = U,
    d = fit$d,
    V = V,
    M_hat = softImpute::complete(M, fit),
    rank = length(fit$d)
  )
}

pair_dcov <- function(x, y, R = 199, index = 1) {
  ok <- is.finite(x) & is.finite(y)
  x2 <- x[ok]
  y2 <- y[ok]
  if (length(x2) < 3L) {
    return(list(stat = NA_real_, p = NA_real_))
  }
  tt <- energy::dcor.test(x2, y2, R = R, index = index)
  list(stat = unname(tt$statistic), p = tt$p.value)
}

compute_dcor_mat <- function(MEs, Tmat, n_min = 5) {
  K <- ncol(MEs)
  C <- ncol(Tmat)
  out <- matrix(NA_real_, K, C, dimnames = list(colnames(MEs), colnames(Tmat)))
  for (i in seq_len(K)) {
    x <- MEs[, i]
    if (!all(is.finite(x)) || stats::sd(x, na.rm = TRUE) == 0) {
      next
    }
    for (j in seq_len(C)) {
      y <- Tmat[, j]
      ok <- is.finite(x) & is.finite(y)
      x2 <- x[ok]
      y2 <- y[ok]
      if (length(x2) < n_min) {
        next
      }
      if (stats::sd(y2) == 0) {
        next
      }
      if (all(y2 %in% c(0, 1)) && length(unique(y2)) < 2) {
        next
      }
      out[i, j] <- energy::dcor(x2, y2)
    }
  }
  out
}

make_Tmat <- function(cluster_vec, levels_all) {
  T <- stats::model.matrix(~ 0 + factor(cluster_vec, levels = levels_all))
  colnames(T) <- paste0("C", levels_all)
  T
}

cluster_index_list <- function(cluster_vec, prefix = "C") {
  lev <- sort(unique(cluster_vec))
  out <- lapply(lev, function(cl) which(cluster_vec == cl))
  names(out) <- paste0(prefix, seq_along(out))
  out
}

boot_dcor_by_resampling_cluster <- function(
    MEs,
    cluster,
    B = 1000,
    n_min = 5,
    seed = 2025,
    keep_boot_array = FALSE
) {
  stopifnot(nrow(MEs) == length(cluster))
  set.seed(seed)
  n <- nrow(MEs)
  lev_all <- as.character(sort(unique(cluster)))
  Tmat0 <- make_Tmat(cluster, lev_all)
  dcor_hat <- compute_dcor_mat(MEs, Tmat0, n_min = n_min)
  K <- ncol(MEs)
  C <- length(lev_all)
  
  if (keep_boot_array) {
    boot_arr <- array(
      NA_real_,
      dim = c(K, C, B),
      dimnames = c(dimnames(dcor_hat), list(NULL))
    )
  } else {
    boot_arr <- array(NA_real_, dim = c(K, C, B))
  }
  
  for (b in seq_len(B)) {
    idx_b <- sample.int(n, n, replace = TRUE)
    MEs_b <- MEs[idx_b, , drop = FALSE]
    cl_b <- cluster[idx_b]
    Tmat_b <- make_Tmat(cl_b, lev_all)
    boot_arr[, , b] <- compute_dcor_mat(MEs_b, Tmat_b, n_min = n_min)
  }
  
  mean_mat <- apply(boot_arr, c(1, 2), function(v) mean(v, na.rm = TRUE))
  sd_mat <- apply(boot_arr, c(1, 2), function(v) stats::sd(v, na.rm = TRUE))
  ci_lo <- apply(boot_arr, c(1, 2), function(v) stats::quantile(v, 0.025, na.rm = TRUE, names = FALSE))
  ci_hi <- apply(boot_arr, c(1, 2), function(v) stats::quantile(v, 0.975, na.rm = TRUE, names = FALSE))
  
  if (!keep_boot_array) {
    boot_arr <- NULL
  }
  
  list(
    dcor_hat = dcor_hat,
    boot_arr = boot_arr,
    mean = mean_mat,
    sd = sd_mat,
    ci_lo = ci_lo,
    ci_hi = ci_hi
  )
}

getClusterAssoc <- function(
    fit,
    X,
    n_clusters,
    component_indices = NULL,
    scale_u = TRUE,
    km_iter_max = 10000,
    km_nstart = 1000,
    soft_rank = 1,
    soft_lambda = 0,
    dcor_R = 199,
    dcor_index = 1,
    bootstrapping = FALSE,
    bootstrap_B = 1000,
    bootstrap_seed = 2025,
    bootstrap_n_min = 5,
    keep_boot_array = FALSE,
    mask_align = c("strict", "truncate", "pad_false"),
    block_align = c("strict", "pad", "truncate"),
    block_pad_value = NA_real_
) {
  mask_align <- match.arg(mask_align)
  block_align <- match.arg(block_align)
  
  if (is.null(fit$results_list)) {
    stop("fit must be a Trisfsvd result containing results_list.")
  }
  if (!is.list(X) || length(X) == 0 || !is.list(X[[1]])) {
    stop("X must be a nested list in the format expected by Trisfsvd.")
  }
  if (n_clusters < 1L) {
    stop("n_clusters must be positive.")
  }
  
  n_comp_available <- length(fit$results_list)
  if (is.null(component_indices)) {
    component_indices <- seq_len(min(n_comp_available, n_clusters))
  }
  component_indices <- as.integer(component_indices)
  
  if (any(component_indices < 1L | component_indices > n_comp_available)) {
    stop("component_indices must refer to valid entries of fit$results_list.")
  }
  
  component_results <- fit$results_list[component_indices]
  u_matrix <- do.call(cbind, lapply(component_results, `[[`, "u"))
  if (scale_u) {
    u_for_kmeans <- scale(u_matrix)
  } else {
    u_for_kmeans <- u_matrix
  }
  
  km <- stats::kmeans(
    u_for_kmeans,
    centers = n_clusters,
    iter.max = km_iter_max,
    nstart = km_nstart
  )
  
  selected_subregion <- lapply(component_results, function(res) {
    lapply(res$varphi_list, function(v) {
      as.logical(v != 0)
    })
  })
  feature_cluster <- lapply(selected_subregion, function(mask_list) {
    which(vapply(mask_list, any, logical(1)))
  })
  names(feature_cluster) <- paste0("F", seq_along(feature_cluster))
  
  pc_scores <- NULL
  soft_svd_results <- vector("list", length(component_results))
  masked_data_list <- vector("list", length(component_results))
  masked_matrix_list <- vector("list", length(component_results))
  
  for (k in seq_along(component_results)) {
    X_new <- mask_X_by_bool_list(
      X = X,
      bool_list = selected_subregion[[k]],
      align = mask_align,
      return_index = TRUE
    )
    M <- as_block_matrix(
      X_new,
      align = block_align,
      pad_value = block_pad_value
    )
    svd_soft <- svd_with_na_softImpute(M, rank = soft_rank, lambda = soft_lambda)
    
    pc_scores <- cbind(pc_scores, svd_soft$U[, 1, drop = FALSE])
    soft_svd_results[[k]] <- svd_soft
    masked_data_list[[k]] <- X_new
    masked_matrix_list[[k]] <- M
  }
  
  colnames(pc_scores) <- paste0("Comp", component_indices)
  
  lev <- as.character(sort(unique(km$cluster)))
  cluster_indicator <- stats::model.matrix(~ 0 + factor(km$cluster, levels = lev))
  colnames(cluster_indicator) <- paste0("C", lev)
  
  p <- ncol(pc_scores)
  q <- ncol(cluster_indicator)
  
  dcov_stat_mat <- matrix(NA_real_, p, q, dimnames = list(colnames(pc_scores), colnames(cluster_indicator)))
  
  for (i in seq_len(p)) {
    for (j in seq_len(q)) {
      x <- pc_scores[, i]
      y <- cluster_indicator[, j]
      res <- pair_dcov(x, y, R = dcor_R, index = dcor_index)
      dcov_stat_mat[i, j] <- res$stat
    }
  }
  
  out <- list(
    sample_cluster = cluster_index_list(km$cluster, prefix = "S"),
    feature_cluster = feature_cluster,
    dcov_stat_matrix = dcov_stat_mat
  )
  
  if (bootstrapping) {
    res_bootstrapping <- boot_dcor_by_resampling_cluster(
      MEs = pc_scores,
      cluster = km$cluster,
      B = bootstrap_B,
      n_min = bootstrap_n_min,
      seed = bootstrap_seed,
      keep_boot_array = keep_boot_array
    )
    out$dcor_hat <- res_bootstrapping$dcor_hat
    out$ci_lo <- res_bootstrapping$ci_lo
    out$ci_hi <- res_bootstrapping$ci_hi
    out$se_hat <- res_bootstrapping$sd
  }
  
  out
}
