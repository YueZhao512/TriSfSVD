reconstruct_rankK_var <- function(results_list, j, ranks = NULL, use_s = TRUE) {
  if (is.null(ranks)) {
    ranks <- seq_along(results_list)
  }
  ranks <- ranks[ranks <= length(results_list)]
  if (length(ranks) == 0) {
    stop("No valid ranks available in results_list.")
  }
  
  Xhat_j <- NULL
  for (h in ranks) {
    res_h <- results_list[[h]]
    u_h <- as.numeric(res_h$u)
    phi_h <- as.numeric(res_h$varphi_list[[j]])
    s_h <- if (use_s) as.numeric(res_h$s)[1] else 1
    if (!is.finite(s_h)) {
      s_h <- 1
    }
    part_h <- s_h * outer(u_h, phi_h)
    if (is.null(Xhat_j)) {
      Xhat_j <- part_h
    } else {
      Xhat_j <- Xhat_j + part_h
    }
  }
  Xhat_j
}

make_raw_mat_one_feature <- function(X, j) {
  do.call(rbind, lapply(seq_along(X), function(i) {
    as.numeric(X[[i]][[j]])
  }))
}

plot_association_one_feature <- function(
    j,
    sample_index,
    X,
    results_list,
    ranks,
    mean = TRUE,
    use_s = TRUE,
    col_use = "firebrick",
    raw_alpha = 0.35,
    recon_alpha = 0.8
) {
  raw_mat <- make_raw_mat_one_feature(X, j)
  recon_mat <- reconstruct_rankK_var(results_list, j = j, ranks = ranks, use_s = use_s)
  
  raw_sel <- raw_mat[sample_index, , drop = FALSE]
  recon_sel <- recon_mat[sample_index, , drop = FALSE]
  L <- ncol(raw_sel)
  
  if (mean) {
    raw_mean <- colMeans(raw_sel, na.rm = TRUE)
    recon_mean <- colMeans(recon_sel, na.rm = TRUE)
    ylim_use <- range(c(raw_mean, recon_mean), na.rm = TRUE)
  } else {
    ylim_use <- range(c(raw_sel, recon_sel), na.rm = TRUE)
  }
  
  if (!all(is.finite(ylim_use))) {
    ylim_use <- c(-1, 1)
  }
  d0 <- diff(ylim_use)
  if (!is.finite(d0) || d0 == 0) {
    d0 <- 0.2
  }
  pad <- 0.12 * d0
  ylim_use <- c(ylim_use[1] - pad, ylim_use[2] + pad)
  
  graphics::plot(
    seq_len(L),
    rep(NA_real_, L),
    type = "n",
    ylim = ylim_use,
    xlim = c(1, L),
    xlab = "Domain",
    ylab = "Value",
    main = paste0("Feature ", j)
  )
  graphics::abline(h = 0, col = "gray80", lty = 2)
  
  if (mean) {
    raw_col <- grDevices::adjustcolor(col_use, alpha.f = raw_alpha)
    graphics::points(seq_len(L), raw_mean, pch = 16, cex = 0.8, col = raw_col)
    graphics::lines(seq_len(L), recon_mean, col = col_use, lwd = 2.5)
  } else {
    raw_col <- grDevices::adjustcolor("gray40", alpha.f = raw_alpha)
    recon_col <- grDevices::adjustcolor(col_use, alpha.f = recon_alpha)
    for (ii in seq_len(nrow(raw_sel))) {
      graphics::lines(seq_len(L), raw_sel[ii, ], col = raw_col, lwd = 1)
    }
    for (ii in seq_len(nrow(recon_sel))) {
      graphics::lines(seq_len(L), recon_sel[ii, ], col = recon_col, lwd = 1.5)
    }
  }
  graphics::box()
}

plotResult <- function(
    fit,
    X,
    sample_cluster,
    feature_cluster,
    association_indices = NULL,
    mean = TRUE,
    use_s = TRUE,
    show_variables = TRUE,
    random_seed = NULL,
    ranks = NULL,
    nrow_page = 4,
    ncol_page = 5,
    colours = NULL
) {
  if (is.null(fit$results_list)) {
    stop("fit must be a Trisfsvd result containing results_list.")
  }
  if (!is.list(sample_cluster) || !is.list(feature_cluster)) {
    stop("sample_cluster and feature_cluster must both be lists of index vectors.")
  }
  if (is.null(association_indices)) {
    association_indices <- seq_len(min(length(sample_cluster), length(feature_cluster)))
  }
  association_indices <- as.integer(association_indices)
  if (any(association_indices < 1L)) {
    stop("association_indices must be positive.")
  }
  if (is.null(ranks)) {
    ranks <- seq_along(fit$results_list)
  }
  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }
  if (is.null(colours)) {
    colours <- grDevices::hcl.colors(length(association_indices), "Dark 3")
  }
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  
  selected_feature_output <- vector("list", length(association_indices))
  names(selected_feature_output) <- paste0("Association_", association_indices)
  
  for (idx_pos in seq_along(association_indices)) {
    k <- association_indices[idx_pos]
    s_idx <- sample_cluster[[k]]
    f_idx <- feature_cluster[[k]]
    
    if (isTRUE(show_variables)) {
      feat_show <- f_idx
    } else if (is.numeric(show_variables) && length(show_variables) == 1L) {
      n_show <- min(length(f_idx), as.integer(show_variables))
      feat_show <- sort(sample(f_idx, size = n_show, replace = FALSE))
    } else {
      stop("show_variables must be TRUE or a single integer.")
    }
    
    selected_feature_output[[idx_pos]] <- feat_show
    
    n_panels <- length(feat_show)
    n_per_page <- nrow_page * ncol_page
    start_idx <- seq(1L, n_panels, by = n_per_page)
    
    for (page_start in start_idx) {
      page_end <- min(page_start + n_per_page - 1L, n_panels)
      feat_page <- feat_show[page_start:page_end]
      graphics::par(mfrow = c(nrow_page, ncol_page), mar = c(3, 3, 2, 1), oma = c(0, 0, 2, 0))
      for (j in feat_page) {
        plot_association_one_feature(
          j = j,
          sample_index = s_idx,
          X = X,
          results_list = fit$results_list,
          ranks = ranks,
          mean = mean,
          use_s = use_s,
          col_use = colours[idx_pos]
        )
      }
      empty_panels <- n_per_page - length(feat_page)
      if (empty_panels > 0) {
        for (ii in seq_len(empty_panels)) {
          graphics::plot.new()
        }
      }
      graphics::mtext(
        paste0("Association ", k, ": Sample cluster ", k, " with Feature cluster ", k),
        outer = TRUE,
        cex = 1.1
      )
    }
  }
  
  invisible(selected_feature_output)
}
