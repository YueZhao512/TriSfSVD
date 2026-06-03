###iterative in tuning selection###

# Package dependencies are declared in DESCRIPTION and NAMESPACE.

####sparse functional svd#######
Trisfsvd <- function(X, k = 1, gamma_candidates, alpha_candidates, theta_candidates, lambda_candidates = NULL,
                                   max_outer_iter = 20, tol_outer = 1e-5, use_sgl = TRUE, verbose = TRUE, L_init = NULL,
                                   max_iter = 30000, eps = 1e-5, gamma_adaptive_power = 1, theta_adaptive_power = 1,
                                   lambda_adaptive_power = 1, extend_weight_domain, extend_weight_smooth,extend_weight_group, subj_overlap = FALSE, var_overlap = FALSE, 
                                   subreg_overlap = FALSE, parallel = FALSE, n_cores = detectCores() - 1,extend_weight_u,backtracking) {
  view_sizes <- c(length(X[[1]]))
  results_list <- vector("list", k)
  X_res <- X
  ######integrative add######
  X_views <- split_X_into_views(X, view_sizes)
  index_view = get_non_na_indices_nested(X_views)
  view_var_index = make_view_indices(view_sizes)
  d_list = as.list(vapply(X[[1]], length, integer(1)))
  ###########################
  index = get_non_na_indices_nested(X)
  Sobj <- build_s_structure(X, index, d_list, use_sparse = TRUE)
  Sobj_views <- subset_s_structure_by_view(Sobj, view_var_index)
  
  u_index = c()
  varphi_index = zero_index_varphi = vector("list", length(X[[1]]))
  for (j in 1:length(X[[1]])) {
    varphi_index[[j]] = integer(0) 
  }
  # bic_gamma <- bic_alpha <- bic_lambda <- bic_theta <- list()
  
  iter <- 0
  
  if (parallel) {
    if (n_cores < 1 || n_cores > detectCores()) {
      stop("n_cores must be between 1 and the number of available cores.")
    }
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl))
    
    clusterEvalQ(cl, {
      library(Matrix)
      library(methods)
      NULL
    })
    
    clusterExport(
      cl,
      c(
        "fista_ga_lasso", "fista_sgl_lasso",
        "calc_bic", "calc_bic_smooth", "calc_bic_group",
        "compute_gradient_f",
        "prox_group_l2", "prox_sgl", "soft_threshold",
        "backtracking_update",
        "compute_f_value", "check_stopping", "compute_lipschitz",
        "calc_df_sgl_smooth",
        "normalize_design_matrix", "normalize_numeric_matrix", "normalize_numeric_vector",
        "compute_s_fast", "compute_s_fast_onevar",
        "alpha_worker_task", "lambda_worker_task", "theta_path_worker_task"
      ),
      envir = environment()
    )
    
    if (verbose) cat("Running in parallel mode with", n_cores, "cores.\n")
  } else {
    cl <- NULL
    if (verbose) cat("Running in serial mode.\n")
  }
  
  BIC_score = variance_explained = c()
  X_norm = frobenius_norm_sq(X)
  N = count_non_na_optimized(X)
  # N = length(X);
  # p = length(X[[1]])*length(X[[1]][[1]]);
  df_temp = 0
  d_view <- length(view_sizes)
  for (r in seq_len(k)) {
    iter <- r
    if (verbose) cat(sprintf("\n===== [Wrapper] Extracting rank-%d factor =====\n", r))
    # Call the appropriate version based on the 'parallel' parameter
    res <- bicluster_fsvd(X = X_res, d_view = d_view, gamma_candidates = gamma_candidates, alpha_candidates = alpha_candidates, 
                          theta_candidates = theta_candidates, lambda_candidates = lambda_candidates, 
                          max_outer_iter = max_outer_iter, tol_outer = tol_outer, use_sgl = use_sgl, 
                          verbose = verbose, L_init = L_init, max_iter = max_iter, eps = eps, 
                          gamma_adaptive_power = gamma_adaptive_power, theta_adaptive_power = theta_adaptive_power, 
                          lambda_adaptive_power = lambda_adaptive_power, extend_weight_domain = extend_weight_domain, extend_weight_smooth = extend_weight_smooth,extend_weight_group = extend_weight_group, 
                          u_index = u_index, varphi_index = varphi_index, parallel = parallel, cl = cl,extend_weight_u = extend_weight_u, index = index,backtracking = backtracking, 
                          index_view = index_view,
                          view_var_index = view_var_index,d_list = d_list,Sobj = Sobj,Sobj_views = Sobj_views)
    
    results_list[[r]] <- res
    # bic_gamma[[r]] = res$gamma_bic_result
    # bic_alpha[[r]] = res$alpha_bic_result
    # bic_lambda[[r]] = res$lambda_bic_result
    # bic_theta[[r]] = res$theta_bic_result
    
    if (!subj_overlap) {
      # u_index_temp <- which(res$u != 0 & abs(res$u) / ifelse(max(abs(res$u)) > 0, max(abs(res$u)), 1) > 0.06)
      u_index_temp <- which(res$u != 0 & abs(res$u) / ifelse(max(abs(res$u)) > 0, max(abs(res$u)), 1) > 0)
      u_index = unique(sort(c(u_index, u_index_temp)))
    }
    
    for (j in 1:length(X[[1]])) {
      zero_index_varphi[[j]] <- which(res$varphi_list[[j]] != 0 & 
                                        abs(res$varphi_list[[j]]) / ifelse(max(abs(res$varphi_list[[j]])) > 0, 
                                                                           max(abs(res$varphi_list[[j]])), 1) > 0.05)
      if (!var_overlap) {
        if (length(zero_index_varphi[[j]]) > length(X[[1]][[j]])/20) {
          varphi_index[[j]] = seq_len(length(X[[1]][[j]]))
        }
      } else {
        if (!subreg_overlap) {
          varphi_index[[j]] = unique(sort(c(varphi_index[[j]], zero_index_varphi[[j]])))
        }
      }
    }
    if ((!subreg_overlap || !var_overlap || !subj_overlap) && 
        (length(u_index) == length(X_res) || all(sapply(seq_along(varphi_index), function(i) {
          length(varphi_index[[i]]) == length(res$varphi_list[[i]])
        })))) {
      cat("All samples have been selected, terminating iterations early.\n")
      break
    }
    
    if (all(sapply(res$varphi_list, function(vec) all(vec == 0)))) {
      cat("varphi_list are all 0.\n")
      break
    }
    
    X_res <- subtract_rank1_from_X(X_res, res$u, res$varphi_list, res$s,index = index)
    X_res_norm = frobenius_norm_sq(X_res)
    nz_u = sum(res$u != 0);
    nz_phi = sum(unlist(res$varphi_list) != 0)
    df_temp = df_temp + (nz_u + nz_phi - 1)
    variance_explained[r] = 1 - X_res_norm/X_norm
    BIC_score[r] = N*log(X_res_norm/N)+df_temp*log(N)
  }
  
  results_list <- results_list[1:iter]
  sample_cluster <- lapply(results_list, function(res) {
    which(as.numeric(res$u) != 0)
  })
  names(sample_cluster) <- paste0("S", seq_along(sample_cluster))
  
  feature_cluster <- lapply(results_list, function(res) {
    which(!vapply(res$varphi_list, function(v) all(as.numeric(v) == 0), logical(1)))
  })
  names(feature_cluster) <- paste0("F", seq_along(feature_cluster))
  # bic_gamma <- bic_gamma[1:iter]
  # bic_alpha <- bic_alpha[1:iter]
  # bic_lambda <- bic_lambda[1:iter]
  # bic_theta <- bic_theta[1:iter]
  
  return(list(
    results_list = results_list, 
    sample_cluster = sample_cluster,
    feature_cluster = feature_cluster,
    variance_explained = variance_explained,
    BIC_score = BIC_score
    # bic_gamma = bic_gamma,
    # bic_alpha = bic_alpha,
    # bic_lambda = bic_lambda,
    # bic_theta = bic_theta
  ))
}

####SVD for ols#######
compute_svd_ols <- function(data_list,k = 1,fill_method = "zero",threshold = 0.5) {
  fill_method <- match.arg(fill_method)
  n <- length(data_list)
  p <- length(data_list[[1]])
  Tpts <- length(data_list[[1]][[1]])
  
  # 1) Build the original matrix M_orig (with NAs) and observation mask
  M_orig <- matrix(NA_real_, nrow = n, ncol = p * Tpts)
  for (i in seq_len(n)) {
    M_orig[i, ] <- unlist(data_list[[i]])
  }
  obs_mask <- !is.na(M_orig)
  
  # 2) Impute missing entries into M for SVD
  M <- M_orig
  if (fill_method == "mean") {
    col_means <- colMeans(M, na.rm = TRUE)
    na_idx <- which(is.na(M), arr.ind = TRUE)
    M[na_idx] <- col_means[na_idx[, 2]]
  } else {
    M[is.na(M)] <- 0
  }
  
  # 3) Perform SVD with threshold-based choice and error handling
  min_dim <- min(dim(M))
  if (k >= threshold * min_dim) {
    svd_res <- tryCatch({
      svd(M, nu = k, nv = k)
    }, error = function(e) {
      message("Error in svd(): ", e$message)
      NULL
    })
  } else {
    svd_res <- tryCatch({
      irlba::irlba(M, nu = k, nv = k)
    }, error = function(e) {
      message("Error in irlba(): ", e$message)
      NULL
    })
  }
  
  # 4) Extract u, v, d or fallback to zeros; truncate d to length k
  if (is.null(svd_res)) {
    u <- matrix(0, nrow = n, ncol = k)
    v <- matrix(0, nrow = ncol(M), ncol = k)
    d <- rep(0, k)
  } else {
    u <- svd_res$u
    v <- svd_res$v
    d <- head(svd_res$d, k)
  }
  
  # 5) Reconstruct low-rank approximation M_hat and compute squared Frobenius norm
  M_hat <- u %*% diag(d, nrow = k, ncol = k) %*% t(v)
  diff_mat <- (M - M_hat) * obs_mask
  frob_norm_sq <- sum(diff_mat^2)
  
  # 6) Compute residual_mat = M_orig - M_hat (preserves NAs in original positions)
  residual_mat <- M_orig - M_hat
  
  # 7) Build residual_list in same nested structure as data_list
  residual_list <- vector("list", n)
  for (i in seq_len(n)) {
    residual_list[[i]] <- vector("list", p)
    for (j in seq_len(p)) {
      start <- (j - 1) * Tpts + 1
      end   <- j * Tpts
      # extract the i-th row, columns start:end
      residual_list[[i]][[j]] <- residual_mat[i, start:end]
    }
  }
  
  # 8) Return all results
  list(
    u             = u,
    d             = d,
    v             = v,
    frob_norm_sq  = frob_norm_sq,
    residual_list = residual_list
  )
}


####frobenius norm of rank1 from X#######
frobenius_norm_sq <- function(X) {
  index <- get_non_na_indices_nested(X)
  
  n <- length(X)         
  p <- length(X[[1]]) 
  
  sum_sq <- 0
  for (i in seq_len(n)) {
    for (j in seq_len(p)) {
      obs_idx <- index[[i]][[j]]  
      for (t in obs_idx) {
        sum_sq <- sum_sq + X[[i]][[j]][t]^2
      }
    }
  }
  return(sum_sq)
}

####count non-missing data of rank1 from X#######
count_non_na_optimized <- function(X) {
  index <- get_non_na_indices_nested(X)
  total_count <- sum(unlist(lapply(index, function(idx) sapply(idx, length))))
  return(total_count)
}

####biculster rank1 estimation#######
bicluster_fsvd <- function(X,d_view, gamma_candidates, alpha_candidates, theta_candidates,
                           lambda_candidates = NULL, max_outer_iter = 20, tol_outer = 1e-5, use_sgl = FALSE, 
                           verbose = TRUE, L_init = NULL, max_iter = 30000, eps = 1e-5, 
                           gamma_adaptive_power = 1, theta_adaptive_power = 1, lambda_adaptive_power = 1, 
                           extend_weight_domain,extend_weight_smooth,extend_weight_group,  u_index, varphi_index, parallel = FALSE, cl = NULL,extend_weight_u, index,backtracking, index_view, view_var_index,d_list,
                           Sobj,Sobj_views) {
  n <- length(X)   
  p <- length(X[[1]])  
  alpha_vec <- rep(alpha_candidates[1], p)    
  if (use_sgl) {
    lambda_vec <- rep(lambda_candidates[1], p)
  } else {
    lambda_vec <- rep(0, p)
  }
  theta <- rep(theta_candidates[1],d_view)
  d <- length(X[[1]][[1]]) 
  # d_list = as.list(vapply(X[[1]], length, integer(1)))
  if (d>1) {
    Omega = get.pen(1:d)$OMEGA 
  } else{
    Omega = matrix(0) 
  }
  Omega <- as.matrix(Omega)
  storage.mode(Omega) <- "double"
  # Omega = get.pen(1:d)$OMEGA 
  # index = get_non_na_indices_nested(X)
  index_NA_list = get_na_indices_optimized(X)
  y_unlist <- lapply(seq_along(X), function(i) {
    lapply(seq_along(X[[i]]), function(j) {
      X[[i]][[j]][index[[i]][[j]]]
    })
  })
  y_list <- lapply(y_unlist, function(element) {
    unlist(element)
  })
  y_tilde_list_na <- combine_vectors_alternating(X)
  y_tilde_list <- lapply(y_tilde_list_na, function(vec) as.numeric(vec[!is.na(vec)]))
  u <- rnorm(n)
  u = u / norm(u, type = "2")
  initial_setting <- methodA_fill_svd_allvar(X, k = 1, fill_method = "zero", u_index = u_index, varphi_index = varphi_index)
  varphi_list <- vector("list", p);
  for (j in seq_len(p)) {
    varphi_list[[j]] <- initial_setting[[j]]$v
  }
  if (verbose) cat("Start bicluster_fsvd with n=", n, ", p=", p, "\n")
  gamma_bic_result = alpha_bic_result = lambda_bic_result = theta_bic_result = list()
  
  for (outer_iter in seq_len(max_outer_iter)) {
    c_val <- 0
    for (j in seq_len(p)) {
      ####new add####
      # c_val <- c_val + alpha_vec[j] * as.numeric(t(varphi_list[[j]]) %*% Omega %*% varphi_list[[j]])
      if (length(varphi_list[[j]]) > 1) {
        c_val <- c_val + alpha_vec[j] * as.numeric(t(varphi_list[[j]]) %*% Omega %*% varphi_list[[j]])
      }
    }
    old_u <- u
    old_varphi_list <- varphi_list
    varphi_list_step1 = lapply(index, function(index_set) {
      extracted_vectors <- lapply(seq_along(varphi_list), function(i) {
        varphi_list[[i]][index_set[[i]]]
      })
      unlist(extracted_vectors)
    })
    w1 <- 1 / (abs(mapply(function(v, y) {
      crossprod(v, y) / crossprod(v, v)
    }, varphi_list_step1, y_list)) + 1e-8) ^ gamma_adaptive_power
    w1 <- w1 / min(w1)
    res_u_tune <- update_u_with_gamma_tuning(
      y_list = y_list,
      varphi_list_step1 = varphi_list_step1,
      w1 = w1,
      c_val = c_val,
      gamma_candidates = gamma_candidates,
      eps_zero = 1e-14,
      verbose = verbose,
      u_index = u_index,
      cl = cl,
      extend_weight_u = extend_weight_u,
      varphi_list = varphi_list,
      Sobj = Sobj
    )

    u <- res_u_tune$u  
    gamma_bic_result[[outer_iter]] <- res_u_tune$bic_result
    
    if (norm(u, type = "2") > 0) {
      u = u / norm(u, type = "2")
    } else {
      set.seed(123)
      u = u + rnorm(length(u), 0, 10^-5)
      u = u / norm(u, type = "2")
    }
    
    # U <- filter_matrix_by_indices_list(create_vector_block_diag(u = u, d = d), index_NA_list)
    # d_list = as.list(vapply(X[[1]], length, integer(1)))
    U <- create_block_and_filter_by_na(u = u,d_list = d_list,index_NA_list = index_NA_list)
    U <- lapply(U, function(U_j) {
      if (ncol(U_j) == 1) {
        normalize_design_matrix(U_j)
      } else {
        Matrix(U_j, sparse = TRUE)
      }
    })
    
    varphi_init_list <- lapply(seq_along(U), function(j) {
      U_j <- U[[j]]
      y_j <- y_tilde_list[[j]]
      
      tryCatch({
        A <- as.matrix(crossprod(U_j))
        b <- as.vector(crossprod(U_j, y_j))
        as.vector(solve(A, b))
      }, error = function(e) {
        numeric(ncol(U_j))
      })
    })
    
    w2 <- sapply(varphi_init_list, function(varphi_j_init) {
      1 / ((sqrt(sum(varphi_j_init^2)) + 1e-14) ^ theta_adaptive_power)
    })
    w2 <- w2 / min(w2)
    
    w3 <- lapply(varphi_init_list, function(varphi_j_init) {
      w3_j <- 1 / ((abs(varphi_j_init) + 1e-14) ^ lambda_adaptive_power)
      w3_j / min(w3_j)
    })
    # Choose serial or parallel version based on the 'parallel' parameter
    if (parallel) {
      res_tune <- update_varphi_step_tuning_globaltheta_parallel(U = U, y_tilde_list = y_tilde_list, varphi_list = varphi_list,
                                                                 u = u, alpha_grid = alpha_candidates, theta_grid = theta_candidates,
                                                                 lambda_grid = lambda_candidates, w2 = w2, w3 = w3, use_sgl = use_sgl,
                                                                 Omega = Omega, L_init = L_init, alpha_vec_start = alpha_vec,
                                                                 lambda_vec_start = lambda_vec, theta_start = theta, tau = 1.01,
                                                                 eps = eps, max_iter = max_iter, backtracking = backtracking, verbose = verbose,
                                                                 extend_weight_domain = extend_weight_domain,extend_weight_smooth = extend_weight_smooth,extend_weight_group = extend_weight_group, varphi_index = varphi_index, 
                                                                 cl = cl,Sobj = Sobj, Sobj_views = Sobj_views,d_list = d_list,d_view = d_view,view_var_index = view_var_index,
                                                                 index_view = index_view)
    } else {
      res_tune <- update_varphi_step_tuning_globaltheta_serial(U = U, y_tilde_list = y_tilde_list, varphi_list = varphi_list,
                                                               u = u, alpha_grid = alpha_candidates, theta_grid = theta_candidates,
                                                               lambda_grid = lambda_candidates, w2 = w2, w3 = w3, use_sgl = use_sgl,
                                                               Omega = Omega, L_init = L_init, alpha_vec_start = alpha_vec,
                                                               lambda_vec_start = lambda_vec, theta_start = theta, tau = 1.01,
                                                               eps = eps, max_iter = max_iter, backtracking = backtracking, verbose = verbose,
                                                               extend_weight_domain = extend_weight_domain, extend_weight_smooth = extend_weight_smooth, extend_weight_group = extend_weight_group, varphi_index = varphi_index,
                                                               Sobj = Sobj, Sobj_views = Sobj_views, d_list = d_list, d_view = d_view, view_var_index = view_var_index,
                                                               index_view = index_view
      )    }
    # norms <- sapply(res_tune$varphi_list, function(vec) sqrt(sum(vec^2)))
    # varphi_list <- Map(function(vec, norm_val) {
    #   if (norm_val == 0) vec else vec / norm_val
    # }, res_tune$varphi_list, norms)
    global_scale <- sqrt(sum(unlist(res_tune$varphi_list)^2))
    norms <- rep(global_scale, length(res_tune$varphi_list))
    varphi_list <- Map(function(vec, norm_val) {
      if (norm_val == 0) vec else vec / norm_val
    }, res_tune$varphi_list, norms)
    # s = compute_s_masked(X, u, varphi_list, index)
    
    
    alpha_vec   <- res_tune$alpha_vec
    lambda_vec  <- res_tune$lambda_vec
    theta       <- res_tune$theta
    alpha_bic_result[[outer_iter]] <- res_tune$alpha_bic_result
    lambda_bic_result[[outer_iter]] <- res_tune$lambda_bic_result
    theta_bic_result[[outer_iter]] <- res_tune$theta_bic_result
    diff_u <- sqrt(mean((u - old_u)^2))
    diff_v <- 0
    for (j in seq_len(p)) {
      diff_v <- diff_v + mean((varphi_list[[j]] - old_varphi_list[[j]])^2)
    }
    diff_v <- sqrt(diff_v)
    cat(sprintf("Iter %d: diff_u=%.6f, diff_v=%.6f\n", outer_iter, diff_u, diff_v))
    
    if (max(diff_u, diff_v) < tol_outer) {
      if (verbose) cat("Converged.\n")
      break
    }
  }
  # print(1)
  # s = compute_s_masked(X, u, varphi_list, index)
  s <- compute_s_fast(Sobj, u, varphi_list)
  return(list(u = u, varphi_list = varphi_list, s = s, gamma_bic_result = gamma_bic_result, 
              alpha_bic_result = alpha_bic_result, lambda_bic_result = lambda_bic_result, 
              theta_bic_result = theta_bic_result
              # ,df_total_score = res_tune$bic_df,bic_total_rss = res_tune$bic_rss
  ))
}


#### Serial Version of update_varphi_step_tuning_globaltheta #######
update_varphi_step_tuning_globaltheta_serial <- function(
    U, y_tilde_list, varphi_list, u,
    alpha_grid, theta_grid, lambda_grid,
    w2, w3, use_sgl = FALSE, Omega,
    L_init = NULL, alpha_vec_start = NULL,
    lambda_vec_start = NULL, theta_start = NULL,
    tau = 1.01, eps = eps, max_iter,
    backtracking = TRUE, verbose = FALSE,
    extend_weight_domain, extend_weight_smooth, extend_weight_group,
    varphi_index,
    max_outer_iter = 1,
    Sobj, Sobj_views,
    d_list, d_view, view_var_index, index_view
) {
  
  # -----------------------
  # basic checks / setup
  # -----------------------
  p <- length(varphi_list)
  stopifnot(length(U) == p, length(y_tilde_list) == p)
  stopifnot(length(d_list) == p)
  stopifnot(length(view_var_index) == d_view)
  
  # map variable -> view
  var_to_view <- integer(p)
  for (d in seq_len(d_view)) {
    var_to_view[view_var_index[[d]]] <- d
  }
  stopifnot(all(var_to_view >= 1L & var_to_view <= d_view))
  
  # init starts
  if (is.null(alpha_vec_start)) {
    alpha_vec_start <- rep(alpha_grid[ceiling(length(alpha_grid) / 2)], p)
  }
  if (is.null(lambda_vec_start)) {
    lambda_vec_start <- rep(lambda_grid[ceiling(length(lambda_grid) / 2)], p)
  }
  if (is.null(theta_start)) {
    theta_start <- rep(theta_grid[ceiling(length(theta_grid) / 2)], d_view)
  }
  stopifnot(length(alpha_vec_start) == p, length(lambda_vec_start) == p)
  stopifnot(length(theta_start) == d_view)
  
  # initialize best parameters
  best_alpha <- alpha_vec_start
  best_lambda <- lambda_vec_start
  best_theta <- theta_start
  best_varphi_list <- varphi_list
  
  alpha_bic_result <- vector("list", p)
  lambda_bic_result <- vector("list", p)
  theta_bic_result <- vector("list", d_view)
  
  # -----------------------
  # outer loop
  # -----------------------
  for (iter in seq_len(max_outer_iter)) {
    if (verbose) cat("Starting iteration", iter, "\n")
    
    idx_keep <- which(d_list != 1)
    idx_drop <- which(d_list == 1)
    
    if (length(idx_drop) > 0) {
      best_alpha[idx_drop] <- 0
      best_lambda[idx_drop] <- 0
    }
    
    # -----------------------
    # Step 1: tune alpha
    # -----------------------
    if (length(idx_keep) > 0) {
      alpha_tasks <- lapply(idx_keep, function(j) {
        list(
          j = j,
          U_j = U[[j]],
          y_j = y_tilde_list[[j]],
          X_j = Sobj$X_list[[j]],
          M_j = Sobj$M_list[[j]],
          varphi_init_j = best_varphi_list[[j]],
          theta_j = best_theta[var_to_view[j]],
          lambda_j_now = best_lambda[j],
          w2_j = w2[j],
          w3_j = if (use_sgl) w3[[j]] else NULL,
          varphi_index_j = varphi_index[[j]],
          alpha_grid = alpha_grid
        )
      })
      
      alpha_results <- lapply(
        alpha_tasks,
        alpha_serial_task,
        use_sgl = use_sgl,
        u = u,
        Omega = Omega,
        L_init = L_init,
        tau = tau,
        eps = eps,
        max_iter = max_iter,
        backtracking = backtracking,
        extend_weight_smooth = extend_weight_smooth
      )
      
      for (res in alpha_results) {
        j <- res$j
        best_alpha[j] <- res$alpha_j
        best_varphi_list[[j]] <- res$varphi_j
        alpha_bic_result[[j]] <- res$bic
      }
    }
    
    if (verbose) cat("Iteration", iter, "- alpha selected:", best_alpha, "\n")
    
    # -----------------------
    # Step 2: tune lambda
    # -----------------------
    if (use_sgl) {
      if (length(idx_keep) > 0) {
        lambda_tasks <- lapply(idx_keep, function(j) {
          list(
            j = j,
            U_j = U[[j]],
            y_j = y_tilde_list[[j]],
            X_j = Sobj$X_list[[j]],
            M_j = Sobj$M_list[[j]],
            varphi_init_j = best_varphi_list[[j]],
            theta_j = best_theta[var_to_view[j]],
            alpha_j = best_alpha[j],
            w2_j = w2[j],
            w3_j = w3[[j]],
            varphi_index_j = varphi_index[[j]],
            lambda_grid = lambda_grid
          )
        })
        
        lambda_results <- lapply(
          lambda_tasks,
          lambda_serial_task,
          u = u,
          Omega = Omega,
          L_init = L_init,
          tau = tau,
          eps = eps,
          max_iter = max_iter,
          backtracking = backtracking,
          extend_weight_domain = extend_weight_domain
        )
        
        for (res in lambda_results) {
          j <- res$j
          best_lambda[j] <- res$lambda_j
          best_varphi_list[[j]] <- res$varphi_j
          lambda_bic_result[[j]] <- res$bic
        }
      }
      
      if (verbose) cat("Iteration", iter, "- lambda selected:", best_lambda, "\n")
    } else {
      best_lambda <- rep(0, length(best_lambda))
      lambda_bic_result <- vector("list", length(best_lambda))
    }
    
    # -----------------------
    # sanity check
    # -----------------------
    all_idx <- unlist(view_var_index, use.names = FALSE)
    stopifnot(length(all_idx) == length(unique(all_idx)))
    stopifnot(setequal(all_idx, seq_len(p)))
    
    # -----------------------
    # Step 3: per-view theta
    # -----------------------
    for (d in seq_len(d_view)) {
      idx_view <- view_var_index[[d]]
      stopifnot(length(idx_view) > 0)
      
      ncol_vec <- vapply(U[idx_view], ncol, integer(1))
      is_p1_view <- all(ncol_vec == 1L)
      stopifnot(is_p1_view || all(ncol_vec > 1L))
      
      U_view <- U[idx_view]
      U_view_bic <- if (is_p1_view) {
        U_view
      } else {
        lapply(U_view, normalize_numeric_matrix)
      }
      y_view <- y_tilde_list[idx_view]
      alpha_view <- best_alpha[idx_view]
      lambda_view <- best_lambda[idx_view]
      w2_view <- w2[idx_view]
      w3_view <- w3[idx_view]
      varphi_index_view <- varphi_index[idx_view]
      varphi_view_init <- best_varphi_list[idx_view]
      
      theta_tasks <- lapply(seq_along(idx_view), function(ii) {
        list(
          U_j = U_view[[ii]],
          y_j = y_view[[ii]],
          varphi_init_j = varphi_view_init[[ii]],
          alpha_j = alpha_view[ii],
          lambda_j = lambda_view[ii],
          w2_j = w2_view[ii],
          w3_j = if (use_sgl) w3_view[[ii]] else NULL,
          varphi_index_j = varphi_index_view[[ii]],
          theta_grid = theta_grid
        )
      })
      
      theta_path_by_var <- lapply(
        theta_tasks,
        theta_path_serial_task,
        use_sgl = use_sgl,
        u = u,
        Omega = Omega,
        L_init = L_init,
        tau = tau,
        eps = eps,
        max_iter = max_iter,
        backtracking = backtracking
      )
      
      best_bic_d <- Inf
      best_theta_d <- best_theta[d]
      bic_curve_d <- numeric(length(theta_grid))
      best_varphi_view <- varphi_view_init
      
      for (k in seq_along(theta_grid)) {
        t_cand <- theta_grid[k]
        
        varphi_temp_view_k <- lapply(theta_path_by_var, `[[`, k)
        
        if (length(d_view) == 1) {
          s_k = 1
        } else{
        s_k <- compute_s_fast(
          Sobj = Sobj_views[[d]],
          u = u,
          varphi_list = varphi_temp_view_k
        )}
        
        bic_view <- calc_bic_group(
          y_tilde_list = y_view,
          U_list = U_view_bic,
          varphi_list = varphi_temp_view_k,
          alpha_vec = alpha_view,
          u = u,
          Omega = Omega,
          theta = t_cand,
          w2 = w2_view,
          extend_weight_group = extend_weight_group,
          s = s_k
        )
        
        bic_curve_d[k] <- bic_view$bic_overall
        
        if (bic_view$bic_overall < best_bic_d) {
          best_bic_d <- bic_view$bic_overall
          best_theta_d <- t_cand
          best_varphi_view <- varphi_temp_view_k
        }
      }
      
      best_theta[d] <- best_theta_d
      theta_bic_result[[d]] <- bic_curve_d
      best_varphi_list[idx_view] <- best_varphi_view
      
      if (verbose) {
        cat(
          "Iteration", iter, "- view", d,
          if (is_p1_view) "(scalar)" else "(functional)",
          "theta selected:", best_theta_d, "\n"
        )
      }
    }
  }
  
  return(list(
    varphi_list = best_varphi_list,
    alpha_vec = best_alpha,
    lambda_vec = best_lambda,
    theta = best_theta,
    theta_bic_result = theta_bic_result,
    lambda_bic_result = lambda_bic_result,
    alpha_bic_result = alpha_bic_result
  ))
}

#### Parallel Version of update_varphi_step_tuning_globaltheta #######
update_varphi_step_tuning_globaltheta_parallel <- function(
    U, y_tilde_list, varphi_list, u,
    alpha_grid, theta_grid, lambda_grid,
    w2, w3, use_sgl = FALSE, Omega,
    L_init = NULL, alpha_vec_start = NULL,
    lambda_vec_start = NULL, theta_start = NULL,
    tau = 1.01, eps = eps, max_iter,
    backtracking = TRUE, verbose = FALSE,
    extend_weight_domain, extend_weight_smooth, extend_weight_group,
    varphi_index, cl,
    max_outer_iter = 1,
    Sobj, Sobj_views,
    d_list, d_view, view_var_index, index_view
) {
  
  # -----------------------
  # basic checks / setup
  # -----------------------
  p <- length(varphi_list)
  stopifnot(length(U) == p, length(y_tilde_list) == p)
  stopifnot(length(d_list) == p)
  stopifnot(length(view_var_index) == d_view)
  
  # map variable -> view
  var_to_view <- integer(p)
  for (d in seq_len(d_view)) {
    var_to_view[view_var_index[[d]]] <- d
  }
  stopifnot(all(var_to_view >= 1L & var_to_view <= d_view))
  
  # init starts
  if (is.null(alpha_vec_start)) {
    alpha_vec_start <- rep(alpha_grid[ceiling(length(alpha_grid) / 2)], p)
  }
  if (is.null(lambda_vec_start)) {
    lambda_vec_start <- rep(lambda_grid[ceiling(length(lambda_grid) / 2)], p)
  }
  if (is.null(theta_start)) {
    theta_start <- rep(theta_grid[ceiling(length(theta_grid) / 2)], d_view)
  }
  
  stopifnot(length(alpha_vec_start) == p, length(lambda_vec_start) == p)
  stopifnot(length(theta_start) == d_view)
  
  # initialize best parameters
  best_alpha <- alpha_vec_start
  best_lambda <- lambda_vec_start
  best_theta <- theta_start
  best_varphi_list <- varphi_list
  
  alpha_bic_result <- vector("list", p)
  lambda_bic_result <- vector("list", p)
  theta_bic_result <- vector("list", d_view)
  
  # ---------------------------------------------------------
  # export only shared small/common objects + functions ONCE
  # ---------------------------------------------------------
  use_sgl_global <- use_sgl
  u_global <- u
  Omega_global <- Omega
  L_init_global <- L_init
  tau_global <- tau
  eps_global <- eps
  max_iter_global <- max_iter
  backtracking_global <- backtracking
  extend_weight_domain_global <- extend_weight_domain
  extend_weight_smooth_global <- extend_weight_smooth
  extend_weight_group_global <- extend_weight_group
  
  
  clusterExport(
    cl,
    c(
      "use_sgl_global", "u_global", "Omega_global", "L_init_global",
      "tau_global", "eps_global", "max_iter_global", "backtracking_global",
      "extend_weight_domain_global", "extend_weight_smooth_global", "extend_weight_group_global",
      "alpha_worker_task", "lambda_worker_task", "theta_path_worker_task",
      "fista_ga_lasso", "fista_sgl_lasso",
      "calc_bic", "calc_bic_smooth", "calc_bic_group",
      "compute_s_fast", "compute_s_fast_onevar",
      "compute_lipschitz", "compute_gradient_f", "compute_f_value",
      "backtracking_update", "calc_df_sgl_smooth",
      "normalize_design_matrix", "normalize_numeric_matrix", "normalize_numeric_vector",
      "check_stopping", "prox_group_l2", "prox_sgl", "soft_threshold"
    ),
    envir = environment()
  )
  
  # -----------------------
  # outer loop
  # -----------------------
  for (iter in seq_len(max_outer_iter)) {
    if (verbose) cat("Starting iteration", iter, "\n")
    idx_keep <- which(d_list != 1)
    idx_drop <- which(d_list == 1)
    if (length(idx_drop) > 0) {
      best_alpha[idx_drop] <- 0
      best_lambda[idx_drop] <- 0
    }
    # -----------------------
    # Step 1: tune alpha
    # -----------------------
    if (length(idx_keep) > 0) {
      alpha_tasks <- lapply(idx_keep, function(j) {
        list(
          j = j,
          U_j = U[[j]],
          y_j = y_tilde_list[[j]],
          X_j = Sobj$X_list[[j]],
          M_j = Sobj$M_list[[j]],
          varphi_init_j = best_varphi_list[[j]],
          theta_j = best_theta[var_to_view[j]],
          lambda_j_now = best_lambda[j],
          w2_j = w2[j],
          w3_j = if (use_sgl) w3[[j]] else NULL,
          varphi_index_j = varphi_index[[j]],
          alpha_grid = alpha_grid
        )
      })
      
      alpha_results <- parLapply(cl, alpha_tasks, alpha_worker_task)
      
      for (res in alpha_results) {
        j <- res$j
        best_alpha[j] <- res$alpha_j
        best_varphi_list[[j]] <- res$varphi_j
        alpha_bic_result[[j]] <- res$bic
      }
    }
    if (verbose) cat("Iteration", iter, "- alpha selected:", best_alpha, "\n")
    
    # -----------------------
    # Step 2: tune lambda
    # -----------------------
    if (use_sgl) {
      if (length(idx_keep) > 0) {
        lambda_tasks <- lapply(idx_keep, function(j) {
          list(
            j = j,
            U_j = U[[j]],
            y_j = y_tilde_list[[j]],
            X_j = Sobj$X_list[[j]],
            M_j = Sobj$M_list[[j]],
            varphi_init_j = best_varphi_list[[j]],
            theta_j = best_theta[var_to_view[j]],
            alpha_j = best_alpha[j],
            w2_j = w2[j],
            w3_j = w3[[j]],
            varphi_index_j = varphi_index[[j]],
            lambda_grid = lambda_grid
          )
        })
        
        lambda_results <- parLapply(cl, lambda_tasks, lambda_worker_task)
        
        for (res in lambda_results) {
          j <- res$j
          best_lambda[j] <- res$lambda_j
          best_varphi_list[[j]] <- res$varphi_j
          lambda_bic_result[[j]] <- res$bic
        }
      }
      
      if (verbose) cat("Iteration", iter, "- lambda selected:", best_lambda, "\n")
    } else {
      best_lambda <- rep(0, length(best_lambda))
      lambda_bic_result <- vector("list", length(best_lambda))
    }
    
    # -----------------------
    # sanity check
    # -----------------------
    all_idx <- unlist(view_var_index, use.names = FALSE)
    stopifnot(length(all_idx) == length(unique(all_idx)))
stopifnot(setequal(all_idx, seq_len(p)))

# -----------------------
# Step 3: per-view theta
# -----------------------
for (d in seq_len(d_view)) {
  idx_view <- view_var_index[[d]]
  stopifnot(length(idx_view) > 0)
  
  ncol_vec <- vapply(U[idx_view], ncol, integer(1))
  is_p1_view <- all(ncol_vec == 1L)
  stopifnot(is_p1_view || all(ncol_vec > 1L))
  
  U_view <- U[idx_view]
  U_view_bic <- if (is_p1_view) {
    U_view
  } else {
    lapply(U_view, normalize_numeric_matrix)
  }
  y_view <- y_tilde_list[idx_view]
  alpha_view <- best_alpha[idx_view]
  lambda_view <- best_lambda[idx_view]
  w2_view <- w2[idx_view]
  w3_view <- w3[idx_view]
  varphi_index_view <- varphi_index[idx_view]
  varphi_view_init <- best_varphi_list[idx_view]
  
  theta_tasks <- lapply(seq_along(idx_view), function(ii) {
    list(
      U_j = U_view[[ii]],
      y_j = y_view[[ii]],
      varphi_init_j = varphi_view_init[[ii]],
      alpha_j = alpha_view[ii],
      lambda_j = lambda_view[ii],
      w2_j = w2_view[ii],
      w3_j = if (use_sgl) w3_view[[ii]] else NULL,
      varphi_index_j = varphi_index_view[[ii]],
      theta_grid = theta_grid
    )
  })
  
  theta_path_by_var <- parLapply(cl, theta_tasks, theta_path_worker_task)
  
  best_bic_d <- Inf
  best_theta_d <- best_theta[d]
  bic_curve_d <- numeric(length(theta_grid))
  best_varphi_view <- varphi_view_init
  
  for (k in seq_along(theta_grid)) {
    t_cand <- theta_grid[k]
    
    varphi_temp_view_k <- lapply(theta_path_by_var, `[[`, k)
    
    if (length(d_view) == 1) {
      s_k <- 1
    } else {
      s_k <- compute_s_fast(
        Sobj = Sobj_views[[d]],
        u = u,
        varphi_list = varphi_temp_view_k
      )
    }
    
    bic_view <- calc_bic_group(
      y_tilde_list = y_view,
      U_list = U_view_bic,
      varphi_list = varphi_temp_view_k,
      alpha_vec = alpha_view,
      u = u,
      Omega = Omega,
      theta = t_cand,
      w2 = w2_view,
      extend_weight_group = extend_weight_group,
      s = s_k
    )
    
    bic_curve_d[k] <- bic_view$bic_overall
    
    if (bic_view$bic_overall < best_bic_d) {
      best_bic_d <- bic_view$bic_overall
      best_theta_d <- t_cand
      best_varphi_view <- varphi_temp_view_k
    }
  }
  
  best_theta[d] <- best_theta_d
  theta_bic_result[[d]] <- bic_curve_d
  best_varphi_list[idx_view] <- best_varphi_view
  
  if (verbose) {
    cat(
      "Iteration", iter, "- view", d,
      if (is_p1_view) "(scalar)" else "(functional)",
      "theta selected:", best_theta_d, "\n"
    )
  }
}
  }
  
  return(list(
    varphi_list = best_varphi_list,
    alpha_vec = best_alpha,
    lambda_vec = best_lambda,
    theta = best_theta,
    theta_bic_result = theta_bic_result,
    lambda_bic_result = lambda_bic_result,
    alpha_bic_result = alpha_bic_result
  ))
}
################# functions to handle missing value ###################
get_non_na_indices_nested <- function(X) {
  fast_get_indices <- function(x) {
    if (is.list(x)) {
      # If the input is a nested list, apply the function recursively to each element
      lapply(x, fast_get_indices)
    } else {
      # If the input is a vector, return the indices of non-NA elements
      which(!is.na(x))
    }
  }
  fast_get_indices(X)
}

subtract_rank1_from_X <- function(X, u, varphi_list, s, index) {
  n <- length(X); p <- length(X[[1]])
  for (i in seq_len(n)) {
    ui <- u[i]
    for (j in seq_len(p)) {
      obs <- index[[i]][[j]]
      if (length(obs)) {
        X[[i]][[j]][obs] <- X[[i]][[j]][obs] - s * ui * varphi_list[[j]][obs]
      }
    }
  }
  X
}

combine_vectors_alternating <- function(y_unlist) {
  # Iterate over each position in the sublists
  y_tilde_list_fixed_length <- lapply(seq_along(y_unlist[[1]]), function(pos) {
    # Extract the vectors at position 'pos' from each sublist
    vectors <- lapply(y_unlist, `[[`, pos)
    
    # Get the length of the vectors (assuming all vectors have the same length)
    len <- length(vectors[[1]])
    
    # Initialize the result vector with NA
    result <- rep(NA, len * length(vectors))
    
    # Fill the result vector in alternating order for each vector
    for (i in seq_along(vectors)) {
      result[seq(i, length(result), by = length(vectors))] <- vectors[[i]]
    }
    
    return(result)  # Return the combined vector for the current position
  })
  
  return(y_tilde_list_fixed_length)  # Return the list of combined vectors
}

# Function to create a vector block diagonal matrix
create_vector_block_diag <- function(u, d) {
  # Calculate the number of rows in the resulting matrix
  n <- length(u)
  rows <- n * d
  
  # Initialize an empty matrix
  mat <- matrix(0, nrow = rows, ncol = d)
  
  # Fill the block diagonal
  for (i in seq_len(d)) {
    start_row <- (i - 1) * n + 1
    end_row <- i * n
    mat[start_row:end_row, i] <- u
  }
  
  return(mat)
}

get_na_indices_optimized <- function(X) {
  vectors_list <- lapply(seq_len(length(X[[1]])), function(i) lapply(X, `[[`, i))
  lapply(vectors_list, function(vectors) {
    combined_vector <- unlist(lapply(seq_along(vectors[[1]]), function(idx) {
      sapply(vectors, `[[`, idx)
    }))
    which(is.na(combined_vector))
  })
}

# Function to filter rows in a matrix based on a list of row indices to delete
filter_matrix_by_indices_list <- function(mat, row_indices_list) {
  # Iterate over the list of row indices
  filtered_matrices <- lapply(row_indices_list, function(row_indices) {
    # Ensure the row indices are within the valid range
    if (any(row_indices > nrow(mat) | row_indices < 1)) {
      stop("Row indices must be within the range of the matrix rows.")
    }
    
    # Create a logical vector to identify rows to keep
    rows_to_keep <- !(seq_len(nrow(mat)) %in% row_indices)
    
    # Subset the matrix to keep only the rows not in row_indices
    mat[rows_to_keep, , drop = FALSE]
  })
  
  return(filtered_matrices)  # Return a list of filtered matrices
}

# integrative version
create_block_and_filter_by_na <- function(u, d_list, index_NA_list) {
  if (length(d_list) != length(index_NA_list)) {
    stop("length(d_list) must equal length(index_NA_list).")
  }
  
  n <- length(u)
  
  out <- lapply(seq_along(index_NA_list), function(k) {
    d_k <- d_list[[k]]
    if (!is.numeric(d_k) || length(d_k) != 1 || is.na(d_k) || d_k < 1) {
      stop(sprintf("d_list[[%d]] must be a single positive number.", k))
    }
    d_k <- as.integer(d_k)
    
    # 1) build (n*d_k) x d_k block matrix
    mat_k <- matrix(0, nrow = n * d_k, ncol = d_k)
    for (j in seq_len(d_k)) {
      rr <- ((j - 1) * n + 1):(j * n)
      mat_k[rr, j] <- u
    }
    
    # 2) drop NA rows for this k
    row_indices <- index_NA_list[[k]]
    if (length(row_indices) == 0) return(mat_k)
    
    if (any(row_indices < 1L | row_indices > nrow(mat_k))) {
      stop(sprintf("index_NA_list[[%d]] has row indices out of range for nrow=%d.",
                   k, nrow(mat_k)))
    }
    
    rows_to_keep <- !(seq_len(nrow(mat_k)) %in% row_indices)
    mat_k[rows_to_keep, , drop = FALSE]
  })
  
  out
}
####################################

# Step(1) update u
solve_u <- function(y_list, varphi_list, gamma, w1, c_val, u_index){
  u <- numeric(length(y_list))
  for (i in 1:length(y_list)) {
    if (i %in% u_index) {
      u[i] <- 0
    } else{
      varphi_i <- varphi_list[[i]]
      y_i <- y_list[[i]]
      # Pre-compute reused values
      phiTy <- sum(varphi_i * y_i)
      phiTphi <- sum(varphi_i * varphi_i)
      if (is.na(w1[i]) ) {
        thresh_part <- abs(phiTy) - 0.5 * gamma 
      } else {
        thresh_part <- abs(phiTy) - 0.5 * gamma * w1[i]}
      # Calculate u[i] with early exit
      if (thresh_part <= 0 || phiTphi + c_val ==0) {
        u[i] <- 0
      } else {
        u[i] <- sign(phiTy) / (phiTphi + c_val) * thresh_part
      }
    }
  }
  return(u)
}

# Function to compute u[i] for a single index
compute_u_i <- function(i, y_list, varphi_list, gamma, w1, c_val, u_index) {
  if (i %in% u_index) {
    return(0)
  } else {
    varphi_i <- varphi_list[[i]]
    y_i <- y_list[[i]]
    phiTy <- sum(varphi_i * y_i)  # Compute dot product of varphi_i and y_i
    phiTphi <- sum(varphi_i * varphi_i)  # Compute squared norm of varphi_i
    if (is.na(w1[i])) {
      thresh_part <- abs(phiTy) - 0.5 * gamma  # Threshold without w1
    } else {
      thresh_part <- abs(phiTy) - 0.5 * gamma * w1[i]  # Threshold with w1
    }
    if (thresh_part <= 0 || phiTphi + c_val == 0) {
      return(0)  # Early exit if threshold condition fails
    } else {
      return(sign(phiTy) / (phiTphi + c_val) * thresh_part)  # Compute u[i]
    }
  }
}

compute_u_chunk <- function(idx, y_list, varphi_list, gamma, w1, c_val, u_index) {
  vapply(
    idx,
    compute_u_i,
    numeric(1),
    y_list = y_list,
    varphi_list = varphi_list,
    gamma = gamma,
    w1 = w1,
    c_val = c_val,
    u_index = u_index
  )
}

# Parallel version of solve_u with internal OS-based method selection
solve_u_parallel <- function(y_list, varphi_list, gamma, w1, c_val, u_index, cl) {
  n <- length(y_list)

  if (n <= 1L) {
    return(solve_u(y_list, varphi_list, gamma, w1, c_val, u_index))
  }

  n_workers <- length(cl)
  if (n <= 20L * n_workers) {
    return(solve_u(y_list, varphi_list, gamma, w1, c_val, u_index))
  }

  idx_chunks <- parallel::splitIndices(n, n_workers)
  u_chunks <- parLapply(
    cl,
    idx_chunks,
    compute_u_chunk,
    y_list = y_list,
    varphi_list = varphi_list,
    gamma = gamma,
    w1 = w1,
    c_val = c_val,
    u_index = u_index
  )

  unlist(u_chunks, use.names = FALSE)
}


normalize_numeric_vector <- function(x) {
  if (is.double(x) && is.atomic(x) && is.null(dim(x))) {
    return(x)
  }
  as.numeric(x)
}

normalize_design_matrix <- function(x) {
  if (is.matrix(x) || inherits(x, "Matrix")) {
    return(x)
  }
  as.matrix(x)
}

normalize_numeric_matrix <- function(x) {
  if (is.matrix(x) && is.double(x)) {
    return(x)
  }
  if (inherits(x, "Matrix") || !is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.double(x)) {
    storage.mode(x) <- "double"
  }
  x
}

# FISTA for group selection
fista_ga_lasso <- function(U, y,alpha, u, theta, w2_j, L_init = NULL, tau, eps, max_iter,
                           varphi_init = NULL,backtracking = TRUE,Omega,varphi_index) {
  U <- normalize_design_matrix(U)
  y <- normalize_numeric_vector(y)
  
  # (A) Initialization
  p <- ncol(U)
  
  if(is.null(varphi_init)) {
    tilde_varphi_prev <- rep(0, p)  # tilde_varphi^(0)
  } else {
    tilde_varphi_prev <- varphi_init
  }
  
  # scalar adaptive lasso:  min_v ||y - U v||^2 + theta*w*|v|
  if (p == 1) {
    a <- as.numeric(crossprod(U))       # U^T U
    b <- as.numeric(crossprod(U, y))    # U^T y
    w <- as.numeric(w2_j)
    
    if (!is.finite(a) || a <= 0) return(0)
    
    z   <- b / a
    thr <- (theta * w) / (2 * a)        # <-- divide by 2 to match "no 1/2" loss
    
    v_hat <- sign(z) * max(abs(z) - thr, 0)
    return(as.numeric(v_hat))
  }
  
  v_current <- tilde_varphi_prev   # v^(1) = tilde_varphi^(0)
  t_current <- 1                  # momentum
  # Omega = get.pen(1:ncol(U))$OMEGA
  if(is.null(L_init)) {
    L <- compute_lipschitz(U, alpha, u, Omega)
  } else {
    L <- L_init
  }
  # (B) Iteration
  for(k in seq_len(max_iter)) {
    # (1) Gradient step
    eta <- 1 / L
    grad_fk <- compute_gradient_f(v_current, U, y, alpha, u, Omega)
    z_k <- v_current - eta * grad_fk
    
    # (2) Prox for Group L2
    tilde_varphi_candidate <- prox_group_l2(z_k, theta, w2_j, eta, varphi_index = varphi_index)
    # (3) If backtracking => check condition
    if(backtracking) {
      # Here, sgl=FALSE => standard group-l2
      res <- backtracking_update(U, y, alpha, u, Omega,v_current, grad_fk,theta, w2_j,
                                 lambda = 0, w3_j = 0, eta, tau, tilde_varphi_candidate,sgl = FALSE,varphi_index = varphi_index)
      tilde_varphi_candidate <- res$beta_temp
      L <- res$L
    }
    
    tilde_varphi_current <- tilde_varphi_candidate
    
    # (4) Momentum
    t_next <- (1 + sqrt(1 + 4*t_current^2)) / 2
    v_next <- tilde_varphi_current + ((t_current - 1)/t_next)*
      (tilde_varphi_current - tilde_varphi_prev)
    
    # (5) Stop check
    if(check_stopping(tilde_varphi_current, tilde_varphi_prev, eps)) {
      # message("Algorithm 1 (Group L2) converged at iteration ", k)
      break
    }
    
    # Update
    tilde_varphi_prev <- tilde_varphi_current
    v_current <- v_next
    t_current <- t_next
  }
  
  return(tilde_varphi_current)
}

# FISTA for group and subregion selection
fista_sgl_lasso <- function(U, y, alpha, u, theta, w2_j, lambda, w3_j,
                            L_init = NULL, tau, eps, max_iter,
                            varphi_init = NULL, backtracking = TRUE,
                            Omega, varphi_index) {
  U <- normalize_design_matrix(U)
  y <- normalize_numeric_vector(y)
  
  p <- ncol(U)
  uu <- sum(u^2)
  
  if (is.null(varphi_init)) {
    tilde_varphi_prev <- rep(0, p)
  } else {
    tilde_varphi_prev <- varphi_init
  }
  
  if (p == 1) {
    a <- as.numeric(crossprod(U))
    b <- as.numeric(crossprod(U, y))
    w <- as.numeric(w2_j)
    
    if (!is.finite(a) || a <= 0) return(0)
    
    z <- b / a
    thr <- (theta * w) / (2 * a)
    v_hat <- sign(z) * max(abs(z) - thr, 0)
    return(as.numeric(v_hat))
  }
  
  v_current <- tilde_varphi_prev
  t_current <- 1
  
  if (is.null(L_init)) {
    L <- compute_lipschitz(U, alpha, uu, Omega)
  } else {
    L <- L_init
  }
  
  for (k in seq_len(max_iter)) {
    eta <- 1 / L
    grad_fk <- compute_gradient_f(v_current, U, y, alpha, uu, Omega)
    z_k <- v_current - eta * grad_fk
    
    tilde_varphi_candidate <- prox_sgl(
      z_k, theta, w2_j, lambda, w3_j, eta,
      varphi_index = varphi_index
    )
    
    if (backtracking) {
      res <- backtracking_update(
        U, y, alpha, uu, Omega,
        v_current, grad_fk, theta, w2_j,
        lambda, w3_j, eta, tau,
        tilde_varphi_candidate,
        sgl = TRUE,
        varphi_index = varphi_index
      )
      tilde_varphi_candidate <- res$beta_temp
      L <- res$L
    }
    
    tilde_varphi_current <- tilde_varphi_candidate
    
    t_next <- (1 + sqrt(1 + 4 * t_current^2)) / 2
    v_next <- tilde_varphi_current +
      ((t_current - 1) / t_next) * (tilde_varphi_current - tilde_varphi_prev)
    
    if (check_stopping(tilde_varphi_current, tilde_varphi_prev, eps)) {
      break
    }
    
    if (k == max(seq_len(max_iter))) {
      message("converge failed", k)
    }
    
    tilde_varphi_prev <- tilde_varphi_current
    v_current <- v_next
    t_current <- t_next
  }
  
  return(tilde_varphi_current)
}

backtracking_update <- function(U, y, alpha, uu, Omega,
                                v_current, grad_fk, theta, w2_j,
                                lambda, w3_j, eta, tau,
                                beta_candidate, sgl = FALSE,
                                varphi_index) {
  success <- FALSE
  L_local <- 1 / eta
  beta_temp <- beta_candidate
  
  while (!success) {
    f_temp <- compute_f_value(beta_temp, U, y, alpha, uu, Omega)
    f_vc   <- compute_f_value(v_current, U, y, alpha, uu, Omega)
    
    lhs <- f_temp
    rhs <- f_vc + sum(grad_fk * (beta_temp - v_current)) +
      (L_local / 2) * sum((beta_temp - v_current)^2)
    
    if (lhs > rhs) {
      L_local <- L_local * tau
      new_eta <- 1 / L_local
      z_k <- v_current - new_eta * grad_fk
      
      if (!sgl) {
        beta_temp <- prox_group_l2(
          z_k, theta, w2_j, new_eta,
          varphi_index = varphi_index
        )
      } else {
        beta_temp <- prox_sgl(
          z_k, theta, w2_j, lambda, w3_j, new_eta,
          varphi_index = varphi_index
        )
      }
      
      if (new_eta < 1e-16) {
        warning("Step size too small in backtracking_update.")
        break
      }
    } else {
      success <- TRUE
    }
  }
  
  return(list(beta_temp = beta_temp, L = L_local))
}

compute_lipschitz <- function(U, alpha, uu, Omega) {
  M <- crossprod(U) + alpha * uu * Omega
  M <- as.matrix(M)
  eigvals <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  L_value <- 2 * max(eigvals)
  return(L_value)
}

compute_gradient_f <- function(x, U, y, alpha, uu, Omega) {
  Ux <- as.vector(U %*% x)
  resid <- y - Ux
  part1 <- -2 * as.vector(crossprod(U, resid))
  part2 <-  2 * alpha * uu * as.vector(Omega %*% x)
  return(part1 + part2)
}

compute_f_value <- function(x, U, y, alpha, uu, Omega) {
  Ux <- as.vector(U %*% x)
  err <- y - Ux
  val_sq <- sum(err^2)
  val_smooth <- alpha * uu * sum(x * as.vector(Omega %*% x))
  return(val_sq + val_smooth)
}

check_stopping <- function(x_new, x_old, eps) {
  diff_vec <- x_new - x_old
  return(sqrt(sum(diff_vec^2)) <= eps)
}

prox_group_l2 <- function(z, theta, w2_j, eta,varphi_index = NULL) {
  norm_z <- sqrt(sum(z^2))
  threshold <- eta * theta * w2_j
  if(norm_z > threshold) {
    scale <- 1 - threshold / norm_z
    z_out <- scale * z
  } else {
    z_out <- rep(0, length(z))
  }
  
  #make required index to 0
  if(length(varphi_index) != 0) {
    z_out[varphi_index] <- 0
  }
  
  return(z_out)
}

prox_sgl <- function(z, theta, w2_j, lambda, w3_j, eta,varphi_index = NULL) {
  # (a) L1 soft-threshold
  z_soft <- soft_threshold(z, eta*lambda*w3_j)
  
  # (b) group-l2
  norm_zs <- sqrt(sum(z_soft^2))
  threshold <- eta * theta * w2_j
  if(norm_zs > threshold) {
    scale <- 1 - threshold / norm_zs
    z_out <- scale * z_soft
  } else {
    z_out <- rep(0, length(z_soft))
  }
  
  #make required index to 0
  if(length(varphi_index) != 0) {
    z_out[varphi_index] <- 0
  }
  return(z_out)
}

soft_threshold <- function(x, thr) {
  sign_x <- sign(x)
  mag <- pmax(abs(x) - thr, 0)
  return(sign_x * mag)
}

get.pen <- function(td) {
  m = length(td);
  h = td[2:m] - td[1:(m-1)]; 
  Q = matrix(0, m, m-1);
  R = matrix(0, m-1, m-1);
  for(k in 2:(m-1))
  {
    Q[k-1,k] = 1/h[k-1];
    Q[k,k] = -1/h[k-1] - 1/h[k];     
    Q[k+1,k] = 1/h[k]
  }
  for(j in 2:(m-2))
  {
    R[j,j] = 1/3 * (h[j-1] + h[j]);
    R[j,j+1] = 1/6 * h[j];
    R[j+1,j] = 1/6 * h[j]
  }
  R[m-1,m-1] = 1/3 * (h[m-2] + h[m-1]);
  s <- solve(R[2:(m-1), 2:(m-1)]) %*% t(Q[1:m, 2:(m-1)]);
  OMEGA = Q[1:m, 2:(m-1)] %*% s;
  EIG.O <- eigen(OMEGA);
  return(list(OMEGA=OMEGA, EIG.O=EIG.O, s=s))
}

####missing data svd for initial setting####
extract_var_matrix <- function(X, var_idx) {
  n <- length(X)  # number of subjects
  d <- length(X[[1]][[var_idx]])  # number of sample points for each variable
  
  out <- matrix(NA, nrow = n, ncol = d)
  for (i in seq_len(n)) {
    out[i, ] <- X[[i]][[var_idx]]
  }
  return(out)
}

fill_na_with_colmean <- function(X_mat) {
  X_filled <- X_mat
  d <- ncol(X_mat)
  for (col in seq_len(d)) {
    col_vals <- X_filled[, col]
    m <- mean(col_vals, na.rm = TRUE)
    X_filled[is.na(col_vals), col] <- m
  }
  return(X_filled)
}

methodA_fill_svd_allvar <- function(X, 
                                    k = 1, 
                                    fill_method = "zero",
                                    u_index = integer(0),     # indices of subjects to exclude from SVD computation
                                    varphi_index = vector("list", length(X[[1]])),  # for each variable, indices of columns to exclude
                                    threshold = 0.5  # threshold: when k exceeds a certain proportion of the allowed dimensions, use base svd
) {
  # Check the fill_method parameter
  if (!fill_method %in% c("zero", "mean")) {
    stop("fill_method must be either 'zero' or 'mean'.")
  }
  
  p <- length(X[[1]])  # number of variables
  n <- length(X)       # number of subjects
  
  # Create a global row mask: exclude subjects specified in u_index from SVD computation
  row_mask <- rep(TRUE, n)
  if (length(u_index) > 0) {
    row_mask[u_index] <- FALSE
  }
  
  # Lists to store processed submatrices, column masks, and original full column counts for each variable
  X_sub_list <- vector("list", p)
  col_masks <- vector("list", p)  # store column masks for each variable
  full_cols <- numeric(p)         # full number of columns for each variable (after filling NA)
  col_sizes <- numeric(p)         # number of allowed columns (TRUE in mask) for each variable
  
  # Process each variable individually
  for (j in seq_len(p)) {
    # 1) Extract the matrix for variable j (n x d_j)
    X_j <- extract_var_matrix(X, j)
    
    # 2) Fill missing values
    if (fill_method == "mean") {
      X_filled <- fill_na_with_colmean(X_j)
    } else {
      X_filled <- X_j
      X_filled[is.na(X_filled)] <- 0
    }
    
    # Save full number of columns for variable j
    full_cols[j] <- ncol(X_filled)
    
    # 3) Create a column mask: exclude columns specified in varphi_index[[j]]
    col_mask <- rep(TRUE, ncol(X_filled))
    if (length(varphi_index[[j]]) > 0) {
      col_mask[varphi_index[[j]]] <- FALSE
    }
    col_masks[[j]] <- col_mask  # save the column mask for later embedding
    
    # 4) Subset the matrix: keep only allowed rows and allowed columns
    X_sub <- X_filled[row_mask, col_mask, drop = FALSE]
    X_sub_list[[j]] <- X_sub
    col_sizes[j] <- ncol(X_sub)  # record the number of allowed columns for variable j
  }
  
  # Combine all variable submatrices horizontally (column bind) into one large matrix
  X_combined <- do.call(cbind, X_sub_list)
  
  # Compute SVD on the combined matrix; use base svd if k exceeds threshold * min(dim), otherwise use irlba
  min_dim <- min(dim(X_combined))
  if (k >= threshold * min_dim) {
    svd_res <- tryCatch({
      svd(X_combined, nu = k, nv = k)
    }, error = function(e) {
      message("Error caught in svd: ", e$message)
      NULL
    })
  } else {
    svd_res <- tryCatch({
      irlba::irlba(X_combined, nu = k, nv = k)
    }, error = function(e) {
      message("Error caught in irlba: ", e$message)
      NULL
    })
  }
  
  # If SVD computation fails, create zero matrices for u_sub and v_sub
  if (is.null(svd_res)) {
    u_sub <- matrix(0, nrow = nrow(X_combined), ncol = k)
    v_sub <- matrix(0, nrow = ncol(X_combined), ncol = k)
    d_combined <- rep(0, k)
  } else {
    u_sub <- svd_res$u
    v_sub <- svd_res$v
    d_combined <- svd_res$d
  }
  
  # Embed u_sub back into a full u matrix corresponding to all subjects (non-selected rows remain 0)
  u_full <- matrix(0, n, k)
  u_full[row_mask, ] <- u_sub
  
  # Partition and embed v_sub for each variable into a full v matrix of original dimension for that variable
  results_list <- vector("list", p)
  col_start <- 1
  for (j in seq_len(p)) {
    cols_allowed <- col_sizes[j]  # number of allowed columns for variable j in the combined matrix
    col_end <- col_start + cols_allowed - 1
    if (cols_allowed > 0) {
      # col_end <- col_start + cols_allowed - 1
      v_segment <- v_sub[col_start:col_end, , drop = FALSE]  # corresponding segment from v_sub
    } else{
      v_segment = 0
    }
    
    # Create a full v matrix with the original number of columns for variable j, filled with zeros
    d_full <- full_cols[j]
    v_full <- matrix(0, nrow = d_full, ncol = k)
    
    # Embed the computed v_segment into the positions indicated by the column mask for variable j
    current_mask <- col_masks[[j]]
    v_full[current_mask, ] <- v_segment
    
    results_list[[j]] <- list(
      u = u_full,      # global u matrix for all subjects (with zeros in non-allowed rows)
      d = d_combined,  # global singular values
      v = v_full       # full v matrix for variable j (non-allowed columns are zeros)
    )
    
    col_start <- col_end + 1
  }
  
  return(results_list)
}


calc_df_sgl_smooth <- function(y_tilde_j, U_j, alpha_j, u, Omega, theta, w2j, beta_hat,
                               tol_active = 1e-14) {
  U_j <- normalize_numeric_matrix(U_j)
  
  uu <- sum(u^2)
  beta_hat <- as.vector(beta_hat)
  
  Idx <- which(abs(beta_hat) > tol_active)
  if (length(Idx) == 0L) return(0)
  
  U_data_sub <- U_j[, Idx, drop = FALSE]
  
  if (alpha_j * uu != 0) {
    A_mat <- alpha_j * uu * Omega
    eps_mat <- 1e-10 * diag(nrow(Omega))
    Q <- chol(A_mat + eps_mat)
    Q_sub <- Q[, Idx, drop = FALSE]
    U_sub <- rbind(as.matrix(U_data_sub), Q_sub)
  } else {
    U_sub <- U_data_sub
  }
  
  xx_sub <- as.matrix(crossprod(U_sub))
  
  dA <- length(Idx)
  del <- matrix(0, dA, dA)
  
  v <- beta_hat[Idx]
  vnorm <- sqrt(sum(v^2))
  if (vnorm > 1e-14) {
    del_group <- diag(dA) - (v %o% v) / (vnorm^2)
    del_group <- (theta * w2j / vnorm) * del_group
    del <- del + del_group
  }
  
  mat_all <- xx_sub + del + 1e-10 * diag(dA)
  mat_inv <- solve(mat_all)
  
  XTX_data <- as.matrix(crossprod(U_data_sub))
  df_val <- sum(mat_inv * XTX_data)
  
  return(df_val)
}

calc_bic <- function(y_tilde_list, U_list, varphi_list, alpha_vec, u, Omega, theta, w2,
                     extend_weight_domain, s, j_var = NULL) {
  if (is.list(varphi_list)) {
    p <- length(varphi_list)
  } else {
    varphi_list <- list(varphi_list)
    p <- 1
  }
  if (!is.list(y_tilde_list)) y_tilde_list <- list(y_tilde_list)
  if (!is.list(U_list)) U_list <- list(U_list)
  
  bic_val <- 0
  
  for (j in seq_len(p)) {
    y_j <- y_tilde_list[[j]]
    U_j <- U_list[[j]]
    beta_hat <- as.vector(varphi_list[[j]])
    
    fit_j <- as.vector(U_j %*% beta_hat)
    resid_j <- y_j - s * fit_j
    RSS_j <- sum(resid_j^2)
    
    df_j <- sum(beta_hat != 0)
    n_j <- length(y_j)
    b_j <- length(beta_hat)
    
    bic_val <- bic_val +
      log(RSS_j / n_j) +
      (log(n_j) / n_j) * df_j +
      2 * extend_weight_domain * df_j * (log(b_j) / n_j)
  }
  
  return(list(
    bic_overall = bic_val
  ))
}

calc_bic_smooth <- function(y_tilde_list, U_list, varphi_list, alpha_vec, u, Omega, theta, w2,
                            extend_weight_smooth, s, j_var = NULL) {
  if (is.list(varphi_list)) {
    p <- length(varphi_list)
  } else {
    varphi_list <- list(varphi_list)
    p <- 1
  }
  if (!is.list(y_tilde_list)) y_tilde_list <- list(y_tilde_list)
  if (!is.list(U_list)) U_list <- list(U_list)
  
  bic_val <- 0
  
  for (j in seq_len(p)) {
    y_j <- y_tilde_list[[j]]
    U_j <- U_list[[j]]
    beta_hat <- as.vector(varphi_list[[j]])
    
    fit_j <- as.vector(U_j %*% beta_hat)
    resid_j <- y_j - s * fit_j
    RSS_j <- sum(resid_j^2)
    
    df_j <- calc_df_sgl_smooth(
      y_tilde_j = y_j, U_j = U_j, alpha_j = alpha_vec[j], u = u, Omega = Omega,
      theta = 0, w2j = w2[j], beta_hat = beta_hat
    )
    
    n_j <- length(y_j)
    b_j <- length(beta_hat)
    
    bic_val <- bic_val +
      log(RSS_j / n_j) +
      (log(n_j) / n_j) * df_j +
      2 * extend_weight_smooth * df_j * (log(b_j) / n_j)
  }
  
  return(list(
    bic_overall = bic_val
  ))
}

calc_bic_group <- function(y_tilde_list, U_list, varphi_list, alpha_vec, u, Omega, theta, w2,
                           extend_weight_group, s, j_var = NULL,
                           tol_active = 1e-12) {
  if (is.list(varphi_list)) {
    p <- length(varphi_list)
  } else {
    varphi_list <- list(varphi_list)
    p <- 1
  }
  if (!is.list(y_tilde_list)) y_tilde_list <- list(y_tilde_list)
  if (!is.list(U_list)) U_list <- list(U_list)
  
  is_p1_view <- all(vapply(U_list, ncol, integer(1)) == 1L)
  
  bic_val <- 0
  
  for (j in seq_len(p)) {
    y_j <- y_tilde_list[[j]]
    U_j <- U_list[[j]]
    beta_hat <- as.vector(varphi_list[[j]])
    
    fit_j <- as.vector(U_j %*% beta_hat)
    resid_j <- y_j - s * fit_j
    RSS_j <- sum(resid_j^2)
    
    n_j <- length(y_j)
    b_j <- length(beta_hat)
    
    if (is_p1_view) {
      df_j <- sum(abs(beta_hat) > tol_active)
    } else {
      df_j <- calc_df_sgl_smooth(
        y_tilde_j = y_j, U_j = U_j, alpha_j = 0, u = u, Omega = Omega,
        theta = theta, w2j = w2[j], beta_hat = beta_hat
      )
    }
    
    bic_val <- bic_val +
      log(RSS_j / n_j) +
      (log(n_j) / n_j) * df_j +
      2 * extend_weight_group * df_j * (log(b_j) / n_j)
  }
  
  list(bic_overall = bic_val)
}

compute_s_masked <- function(X, u, varphi_list, index, ridge = 0, j_var = NULL) {
  n <- length(X); p <- length(varphi_list)
  # stopifnot(length(u) == n, length(varphi_list) == p)
  phi2_list <- lapply(varphi_list, function(v) v * v)
  num <- 0.0
  den <- 0.0
  if (is.null(j_var)) {
    for (j in seq_len(p)) {
      phi_j  <- varphi_list[[j]]
      phi2_j <- phi2_list[[j]]
      for (i in seq_len(n)) {
        obs <- index[[i]][[j]]
        if (length(obs)) {
          ui  <- u[i]
          ui2 <- ui * ui
          xobs <- X[[i]][[j]][obs]
          num <- num + ui  * as.numeric(crossprod(xobs, phi_j[obs]))
          den <- den + ui2 * sum(phi2_j[obs])
        }
      }
    }
  } else{
    phi_j  <- varphi_list[[1]]
    phi2_j <- phi2_list[[1]]
    for (i in seq_len(n)) {
      obs <- index[[i]][[j_var]]
      if (length(obs)) {
        ui  <- u[i]
        ui2 <- ui * ui
        xobs <- X[[i]][[j_var]][obs]
        num <- num + ui  * as.numeric(crossprod(xobs, phi_j[obs]))
        den <- den + ui2 * sum(phi2_j[obs])
      }
    }
  }
  
  denom <- den + ridge
  if (denom > 0) num / denom else 0.0
}


update_u_with_gamma_tuning <- function(y_list, varphi_list_step1,w1,c_val,gamma_candidates,eps_zero = 1e-14,verbose,u_index,cl,extend_weight_u,varphi_list,Sobj){
  c_val <- 0
  # 1) Compute total sample size n_total
  n_total <- 0
  for (i in seq_along(y_list)) {
    n_total <- n_total + length(y_list[[i]])
  }
  
  best_bic   <- Inf
  best_u     <- NULL
  best_gamma <- NA
  best_df    <- NA
  bic_result <- c()
  count = 0
  # 2) Loop over gamma_candidates
  for (gamma_val in gamma_candidates) {
    count = count + 1
    # (a) Solve for u with the given gamma
    if (is.null(cl)) {
      u_test <- solve_u(y_list, varphi_list_step1, gamma_val, w1, c_val, u_index)
    } else {
      u_test <- solve_u_parallel(y_list, varphi_list_step1, gamma_val, w1, c_val, u_index, cl)
    }
    # (b) Normalize u_test
    norm_utest <- sqrt(sum(u_test^2))
    if (norm_utest > eps_zero) {
      u_test <- u_test / norm_utest
    }
    
    # (c) Compute df with the formula:
    df_test <- 0
    for (i in seq_along(u_test)) {
      if (abs(u_test[i]) > eps_zero) {
        varphi_i_temp <- varphi_list_step1[[i]]
        norm_varphi_i_sq <- sum(varphi_i_temp^2)
        df_test_i <- norm_varphi_i_sq / (norm_varphi_i_sq + c_val)
        # df_test_i <- norm_varphi_i_sq / (norm_varphi_i_sq)
        df_test <- df_test + df_test_i
      }
    }
    
    # s = compute_s_masked(X, u_test, varphi_list, index)
    s <- compute_s_fast(Sobj, u_test, varphi_list)
    
    RSS_test <- sum(mapply(function(y_i, phi_i, u_i) {
      r <- y_i - s*u_i * phi_i
      sum(r * r)
    }, y_list, varphi_list_step1, u_test, SIMPLIFY = TRUE))
    bic_result[count] <- log(RSS_test / n_total) + (log(n_total) / n_total) * df_test + 2*extend_weight_u* df_test*(log(length(u_test))/n_total)
    if (bic_result[count] < best_bic) {
      best_bic   <- bic_result[count]
      best_gamma <- gamma_val
      best_u     <- u_test
      best_df    <- df_test
    }
  }
  if (verbose) 
    cat("gamma selected", best_gamma,"\n")
  # 4) Return
  list(
    u          = best_u,
    best_gamma = best_gamma,
    best_bic   = best_bic,
    df         = best_df,
    bic_result = bic_result
  )
}


split_X_into_views <- function(X, view_sizes, return_indices = FALSE) {
  # ---- checks ----
  stopifnot(is.list(X), length(X) > 0)
  stopifnot(is.numeric(view_sizes), length(view_sizes) >= 1)
  if (any(view_sizes <= 0) || any(!is.finite(view_sizes))) {
    stop("view_sizes must be positive finite numbers.")
  }
  view_sizes <- as.integer(view_sizes)
  
  n <- length(X)
  p0 <- length(X[[1]])
  if (!all(vapply(X, length, integer(1)) == p0)) {
    stop("All X[[i]] must have the same length p.")
  }
  if (sum(view_sizes) != p0) {
    stop(sprintf("sum(view_sizes) must equal p (= %d). Got %d.",
                 p0, sum(view_sizes)))
  }
  
  # ---- build index blocks ----
  ends <- cumsum(view_sizes)
  starts <- c(1L, head(ends, -1L) + 1L)
  idx_list <- Map(seq.int, starts, ends)   # list of integer indices per view
  
  # ---- split ----
  views <- lapply(idx_list, function(idx) {
    lapply(X, function(xi) xi[idx])
  })
  
  if (return_indices) {
    attr(views, "view_indices") <- idx_list
  }
  views
}

make_view_indices <- function(view_sizes) {
  stopifnot(is.numeric(view_sizes), length(view_sizes) >= 1)
  if (any(view_sizes <= 0) || any(!is.finite(view_sizes))) {
    stop("view_sizes must be positive finite numbers.")
  }
  view_sizes <- as.integer(view_sizes)
  ends <- cumsum(view_sizes)
  starts <- c(1L, head(ends, -1L) + 1L)
  Map(seq.int, starts, ends)  # list of integer vectors
}

build_s_structure <- function(X, index, d_list, use_sparse = TRUE) {
  n <- length(X)
  p <- length(d_list)
  
  X_list <- vector("list", p)
  M_list <- vector("list", p)
  
  for (j in seq_len(p)) {
    d_j <- d_list[[j]]
    
    if (use_sparse) {
      ii <- integer(0)
      jj <- integer(0)
      xx <- numeric(0)
      
      mi <- integer(0)
      mj <- integer(0)
      mx <- numeric(0)
      
      for (i in seq_len(n)) {
        obs <- index[[i]][[j]]
        if (length(obs) > 0L) {
          vals <- X[[i]][[j]][obs]
          
          ii <- c(ii, rep.int(i, length(obs)))
          jj <- c(jj, obs)
          xx <- c(xx, vals)
          
          mi <- c(mi, rep.int(i, length(obs)))
          mj <- c(mj, obs)
          mx <- c(mx, rep.int(1, length(obs)))
        }
      }
      
      X_list[[j]] <- sparseMatrix(
        i = ii, j = jj, x = xx,
        dims = c(n, d_j)
      )
      
      M_list[[j]] <- sparseMatrix(
        i = mi, j = mj, x = mx,
        dims = c(n, d_j)
      )
    } else {
      Xmat_j <- matrix(0, nrow = n, ncol = d_j)
      Mmat_j <- matrix(0, nrow = n, ncol = d_j)
      
      for (i in seq_len(n)) {
        obs <- index[[i]][[j]]
        if (length(obs) > 0L) {
          Xmat_j[i, obs] <- X[[i]][[j]][obs]
          Mmat_j[i, obs] <- 1
        }
      }
      
      X_list[[j]] <- Xmat_j
      M_list[[j]] <- Mmat_j
    }
  }
  
  list(
    X_list = X_list,
    M_list = M_list
  )
}

subset_s_structure_by_view <- function(Sobj, view_var_index) {
  lapply(view_var_index, function(idx) {
    list(
      X_list = Sobj$X_list[idx],
      M_list = Sobj$M_list[idx]
    )
  })
}

compute_s_fast <- function(Sobj, u, varphi_list, ridge = 0) {
  p <- length(varphi_list)
  u2 <- u * u
  num <- 0
  den <- 0
  
  for (j in seq_len(p)) {
    phi_j <- as.vector(varphi_list[[j]])
    phi2_j <- phi_j * phi_j
    
    Xphi <- as.vector(Sobj$X_list[[j]] %*% phi_j)
    Mphi2 <- as.vector(Sobj$M_list[[j]] %*% phi2_j)
    
    num <- num + sum(u * Xphi)
    den <- den + sum(u2 * Mphi2)
  }
  
  denom <- den + ridge
  if (denom > 0) num / denom else 0
}

compute_s_fast_onevar <- function(X_j, M_j, u, phi_j, ridge = 0) {
  phi_j <- as.vector(phi_j)
  phi2_j <- phi_j * phi_j
  u2 <- u * u
  
  Xphi <- as.vector(X_j %*% phi_j)
  Mphi2 <- as.vector(M_j %*% phi2_j)
  
  num <- sum(u * Xphi)
  den <- sum(u2 * Mphi2)
  
  denom <- den + ridge
  if (denom > 0) num / denom else 0
}

alpha_worker_task <- function(task) {
  # tmp_j <- task$j
  # tmp_theta <- task$theta_j
  # tmp_lambda <- task$lambda_j_now
  # tmp_w2 <- task$w2_j
  # tmp_alpha_grid <- task$alpha_grid
  # tmp_varphi <- task$varphi_init_j
  # tmp_ncolU <- ncol(task$U_j)
  # tmp_leny <- length(task$y_j)
  
  j <- task$j
  U_j <- normalize_numeric_matrix(task$U_j)
  y_j <- normalize_numeric_vector(task$y_j)
  X_j <- normalize_numeric_matrix(task$X_j)
  M_j <- normalize_numeric_matrix(task$M_j)
  
  varphi_init_j <- normalize_numeric_vector(task$varphi_init_j)
  theta_j <- task$theta_j
  lambda_j_now <- task$lambda_j_now
  w2_j <- task$w2_j
  w3_j <- task$w3_j
  varphi_index_j <- task$varphi_index_j
  alpha_grid_local <- task$alpha_grid
  
  best_bic_alpha <- Inf
  best_alpha_j <- alpha_grid_local[1]
  best_varphi_j <- varphi_init_j
  alpha_bic <- numeric(length(alpha_grid_local))
  
  count <- 0
  for (a_cand in alpha_grid_local) {
    count <- count + 1
    
    if (!use_sgl_global) {
      varphi_temp <- fista_ga_lasso(
        U = U_j, y = y_j,
        alpha = a_cand, u = u_global,
        theta = theta_j, w2_j = w2_j,
        L_init = L_init_global, tau = tau_global,
        eps = eps_global, max_iter = max_iter_global,
        varphi_init = varphi_init_j,
        backtracking = backtracking_global,
        Omega = Omega_global,
        varphi_index = varphi_index_j
      )
    } else {
      varphi_temp <- fista_sgl_lasso(
        U = U_j, y = y_j,
        alpha = a_cand, u = u_global,
        theta = theta_j, w2_j = w2_j,
        lambda = lambda_j_now, w3_j = w3_j,
        L_init = L_init_global, tau = tau_global,
        eps = eps_global, max_iter = max_iter_global,
        varphi_init = varphi_init_j,
        backtracking = backtracking_global,
        Omega = Omega_global,
        varphi_index = varphi_index_j
      )
    }
    
    s_j <- compute_s_fast_onevar(
      X_j = X_j,
      M_j = M_j,
      u = u_global,
      phi_j = varphi_temp
    )
    bic_j_val <- calc_bic_smooth(
      y_tilde_list = y_j,
      U_list = U_j,
      alpha_vec = a_cand,
      u = u_global,
      Omega = Omega_global,
      theta = theta_j,
      w2 = w2_j,
      varphi_list = varphi_temp,
      extend_weight_smooth = extend_weight_smooth_global,
      s = s_j,
      j_var = j
    )
    
    alpha_bic[count] <- bic_j_val$bic_overall
    varphi_init_j <- varphi_temp
    
    if (bic_j_val$bic_overall < best_bic_alpha) {
      best_bic_alpha <- bic_j_val$bic_overall
      best_alpha_j <- a_cand
      best_varphi_j <- varphi_temp
    }
  }
  
  list(
    j = j,
    alpha_j = best_alpha_j,
    varphi_j = best_varphi_j,
    bic = alpha_bic
  )
}

lambda_worker_task <- function(task) {
  j <- task$j
  
  U_j <- normalize_numeric_matrix(task$U_j)
  y_j <- normalize_numeric_vector(task$y_j)
  X_j <- normalize_numeric_matrix(task$X_j)
  M_j <- normalize_numeric_matrix(task$M_j)
  
  varphi_init_j <- normalize_numeric_vector(task$varphi_init_j)
  theta_j <- task$theta_j
  alpha_j <- task$alpha_j
  w2_j <- task$w2_j
  w3_j <- task$w3_j
  varphi_index_j <- task$varphi_index_j
  lambda_grid_local <- task$lambda_grid
  
  best_bic_lambda <- Inf
  best_lambda_j <- lambda_grid_local[1]
  best_varphi_j <- varphi_init_j
  lambda_bic <- numeric(length(lambda_grid_local))
  
  count <- 0
  for (l_cand in lambda_grid_local) {
    count <- count + 1
    
    varphi_temp <- fista_sgl_lasso(
      U = U_j, y = y_j,
      alpha = alpha_j, u = u_global,
      theta = theta_j, w2_j = w2_j,
      lambda = l_cand, w3_j = w3_j,
      L_init = L_init_global, tau = tau_global,
      eps = eps_global, max_iter = max_iter_global,
      varphi_init = varphi_init_j,
      backtracking = backtracking_global,
      Omega = Omega_global,
      varphi_index = varphi_index_j
    )
    
    s_j <- compute_s_fast_onevar(
      X_j = X_j,
      M_j = M_j,
      u = u_global,
      phi_j = varphi_temp
    )
    
    bic_j_val <- calc_bic(
      y_tilde_list = y_j,
      U_list = U_j,
      alpha_vec = alpha_j,
      u = u_global,
      Omega = Omega_global,
      theta = theta_j,
      w2 = w2_j,
      varphi_list = varphi_temp,
      extend_weight_domain = extend_weight_domain_global,
      s = s_j,
      j_var = j
    )
    
    lambda_bic[count] <- bic_j_val$bic_overall
    varphi_init_j <- varphi_temp
    
    if (bic_j_val$bic_overall < best_bic_lambda) {
      best_bic_lambda <- bic_j_val$bic_overall
      best_lambda_j <- l_cand
      best_varphi_j <- varphi_temp
    }
  }
  
  list(
    j = j,
    lambda_j = best_lambda_j,
    varphi_j = best_varphi_j,
    bic = lambda_bic
  )
}

theta_path_worker_task <- function(task) {
  U_j <- normalize_numeric_matrix(task$U_j)
  y_j <- normalize_numeric_vector(task$y_j)
  varphi_init_j <- normalize_numeric_vector(task$varphi_init_j)
  alpha_j <- task$alpha_j
  lambda_j <- task$lambda_j
  w2_j <- task$w2_j
  w3_j <- task$w3_j
  varphi_index_j <- task$varphi_index_j
  theta_grid_local <- task$theta_grid
  
  phi_init_j <- varphi_init_j
  phi_path_j <- vector("list", length(theta_grid_local))
  
  for (k in seq_along(theta_grid_local)) {
    t_cand <- theta_grid_local[k]
    
    if (!use_sgl_global) {
      phi_temp <- fista_ga_lasso(
        U = U_j, y = y_j,
        alpha = alpha_j, u = u_global,
        theta = t_cand, w2_j = w2_j,
        L_init = L_init_global, tau = tau_global,
        eps = eps_global, max_iter = max_iter_global,
        varphi_init = phi_init_j,
        backtracking = backtracking_global,
        Omega = Omega_global,
        varphi_index = varphi_index_j
      )
    } else {
      phi_temp <- fista_sgl_lasso(
        U = U_j, y = y_j,
        alpha = alpha_j, u = u_global,
        theta = t_cand, w2_j = w2_j,
        lambda = lambda_j, w3_j = w3_j,
        L_init = L_init_global, tau = tau_global,
        eps = eps_global, max_iter = max_iter_global,
        varphi_init = phi_init_j,
        backtracking = backtracking_global,
        Omega = Omega_global,
        varphi_index = varphi_index_j
      )
    }
    
    phi_path_j[[k]] <- phi_temp
    phi_init_j <- phi_temp
  }
  
  phi_path_j
}

alpha_serial_task <- function(
    task,
    use_sgl, u, Omega, L_init, tau, eps, max_iter, backtracking,
    extend_weight_smooth
) {
  j <- task$j
  
  U_j <- normalize_numeric_matrix(task$U_j)
  y_j <- normalize_numeric_vector(task$y_j)
  X_j <- normalize_numeric_matrix(task$X_j)
  M_j <- normalize_numeric_matrix(task$M_j)
  
  varphi_init_j <- normalize_numeric_vector(task$varphi_init_j)
  theta_j <- task$theta_j
  lambda_j_now <- task$lambda_j_now
  w2_j <- task$w2_j
  w3_j <- task$w3_j
  varphi_index_j <- task$varphi_index_j
  alpha_grid_local <- task$alpha_grid
  
  best_bic_alpha <- Inf
  best_alpha_j <- alpha_grid_local[1]
  best_varphi_j <- varphi_init_j
  alpha_bic <- numeric(length(alpha_grid_local))
  
  count <- 0
  for (a_cand in alpha_grid_local) {
    count <- count + 1
    
    if (!use_sgl) {
      varphi_temp <- fista_ga_lasso(
        U = U_j, y = y_j,
        alpha = a_cand, u = u,
        theta = theta_j, w2_j = w2_j,
        L_init = L_init, tau = tau,
        eps = eps, max_iter = max_iter,
        varphi_init = varphi_init_j,
        backtracking = backtracking,
        Omega = Omega,
        varphi_index = varphi_index_j
      )
    } else {
      varphi_temp <- fista_sgl_lasso(
        U = U_j, y = y_j,
        alpha = a_cand, u = u,
        theta = theta_j, w2_j = w2_j,
        lambda = lambda_j_now, w3_j = w3_j,
        L_init = L_init, tau = tau,
        eps = eps, max_iter = max_iter,
        varphi_init = varphi_init_j,
        backtracking = backtracking,
        Omega = Omega,
        varphi_index = varphi_index_j
      )
    }
    
    s_j <- compute_s_fast_onevar(
      X_j = X_j,
      M_j = M_j,
      u = u,
      phi_j = varphi_temp
    )
    
    bic_j_val <- calc_bic_smooth(
      y_tilde_list = y_j,
      U_list = U_j,
      alpha_vec = a_cand,
      u = u,
      Omega = Omega,
      theta = theta_j,
      w2 = w2_j,
      varphi_list = varphi_temp,
      extend_weight_smooth = extend_weight_smooth,
      s = s_j,
      j_var = j
    )
    
    alpha_bic[count] <- bic_j_val$bic_overall
    varphi_init_j <- varphi_temp
    
    if (bic_j_val$bic_overall < best_bic_alpha) {
      best_bic_alpha <- bic_j_val$bic_overall
      best_alpha_j <- a_cand
      best_varphi_j <- varphi_temp
    }
  }
  
  list(
    j = j,
    alpha_j = best_alpha_j,
    varphi_j = best_varphi_j,
    bic = alpha_bic
  )
}

lambda_serial_task <- function(
    task,
    u, Omega, L_init, tau, eps, max_iter, backtracking,
    extend_weight_domain
) {
  j <- task$j
  
  U_j <- normalize_numeric_matrix(task$U_j)
  y_j <- normalize_numeric_vector(task$y_j)
  X_j <- normalize_numeric_matrix(task$X_j)
  M_j <- normalize_numeric_matrix(task$M_j)
  
  varphi_init_j <- normalize_numeric_vector(task$varphi_init_j)
  theta_j <- task$theta_j
  alpha_j <- task$alpha_j
  w2_j <- task$w2_j
  w3_j <- task$w3_j
  varphi_index_j <- task$varphi_index_j
  lambda_grid_local <- task$lambda_grid
  
  best_bic_lambda <- Inf
  best_lambda_j <- lambda_grid_local[1]
  best_varphi_j <- varphi_init_j
  lambda_bic <- numeric(length(lambda_grid_local))
  
  count <- 0
  for (l_cand in lambda_grid_local) {
    count <- count + 1
    
    varphi_temp <- fista_sgl_lasso(
      U = U_j, y = y_j,
      alpha = alpha_j, u = u,
      theta = theta_j, w2_j = w2_j,
      lambda = l_cand, w3_j = w3_j,
      L_init = L_init, tau = tau,
      eps = eps, max_iter = max_iter,
      varphi_init = varphi_init_j,
      backtracking = backtracking,
      Omega = Omega,
      varphi_index = varphi_index_j
    )
    
    s_j <- compute_s_fast_onevar(
      X_j = X_j,
      M_j = M_j,
      u = u,
      phi_j = varphi_temp
    )
    
    bic_j_val <- calc_bic(
      y_tilde_list = y_j,
      U_list = U_j,
      alpha_vec = alpha_j,
      u = u,
      Omega = Omega,
      theta = theta_j,
      w2 = w2_j,
      varphi_list = varphi_temp,
      extend_weight_domain = extend_weight_domain,
      s = s_j,
      j_var = j
    )
    
    lambda_bic[count] <- bic_j_val$bic_overall
    varphi_init_j <- varphi_temp
    
    if (bic_j_val$bic_overall < best_bic_lambda) {
      best_bic_lambda <- bic_j_val$bic_overall
      best_lambda_j <- l_cand
      best_varphi_j <- varphi_temp
    }
  }
  
  list(
    j = j,
    lambda_j = best_lambda_j,
    varphi_j = best_varphi_j,
    bic = lambda_bic
  )
}

theta_path_serial_task <- function(
    task,
    use_sgl, u, Omega, L_init, tau, eps, max_iter, backtracking
) {
  U_j <- normalize_numeric_matrix(task$U_j)
  y_j <- normalize_numeric_vector(task$y_j)
  varphi_init_j <- normalize_numeric_vector(task$varphi_init_j)
  alpha_j <- task$alpha_j
  lambda_j <- task$lambda_j
  w2_j <- task$w2_j
  w3_j <- task$w3_j
  varphi_index_j <- task$varphi_index_j
  theta_grid_local <- task$theta_grid
  
  phi_init_j <- varphi_init_j
  phi_path_j <- vector("list", length(theta_grid_local))
  
  for (k in seq_along(theta_grid_local)) {
    t_cand <- theta_grid_local[k]
    
    if (!use_sgl) {
      phi_temp <- fista_ga_lasso(
        U = U_j, y = y_j,
        alpha = alpha_j, u = u,
        theta = t_cand, w2_j = w2_j,
        L_init = L_init, tau = tau,
        eps = eps, max_iter = max_iter,
        varphi_init = phi_init_j,
        backtracking = backtracking,
        Omega = Omega,
        varphi_index = varphi_index_j
      )
    } else {
      phi_temp <- fista_sgl_lasso(
        U = U_j, y = y_j,
        alpha = alpha_j, u = u,
        theta = t_cand, w2_j = w2_j,
        lambda = lambda_j, w3_j = w3_j,
        L_init = L_init, tau = tau,
        eps = eps, max_iter = max_iter,
        varphi_init = phi_init_j,
        backtracking = backtracking,
        Omega = Omega,
        varphi_index = varphi_index_j
      )
    }
    
    phi_path_j[[k]] <- phi_temp
    phi_init_j <- phi_temp
  }
  
  phi_path_j
}
