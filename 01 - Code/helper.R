# ============================================================
# Gaussian Graphical Models for Partially Observed
# Multivariate Functional Data
#
# Helper functions and setup
# ============================================================

# ---------------------------
# 1. Package management
# ---------------------------

required_packages <- c("fda", "BDgraph", "Rcpp", "zoo", "pracma")

missing_packages <- required_packages[!vapply(
  required_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)]

if (length(missing_packages) > 0) {
  stop(
    "The following packages are required but not installed: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

# Load packages explicitly used in this script
library(fda)
library(BDgraph)
library(Rcpp)
library(zoo)

# Note:
# dplyr is not used in the current block of code.
# It is better not to attach unused packages in reproducible code.


# ---------------------------
# 2. Source project files
# ---------------------------

# A robust way to source local project files.
# Assumes this script is run from the project root or from a known subfolder.
source_local_code <- function(file) {
  if (!file.exists(file)) {
    stop("Cannot find required source file: ", file, call. = FALSE)
  }
  source(file = file)
}

source_local_cpp <- function(file) {
  if (!file.exists(file)) {
    stop("Cannot find required C++ file: ", file, call. = FALSE)
  }
  Rcpp::sourceCpp(file)
}

# Adjust these paths if needed for your project layout
source_local_code("../01 - Code/pofggm.R")
source_local_cpp("../01 - Code/pofggm.cpp")


# ---------------------------
# 3. Numerical integration weights
# ---------------------------

quadWeights <- function(argvals, method = c("trapezoidal", "midpoint")) {
  method <- match.arg(method)
  
  if (!is.numeric(argvals) || length(argvals) < 2L) {
    stop("'argvals' must be a numeric vector of length at least 2.", call. = FALSE)
  }
  
  if (is.unsorted(argvals, strictly = TRUE)) {
    stop("'argvals' must be strictly increasing.", call. = FALSE)
  }
  
  weights <- switch(
    method,
    trapezoidal = {
      D <- length(argvals)
      0.5 * c(
        argvals[2] - argvals[1],
        argvals[3:D] - argvals[1:(D - 2)],
        argvals[D] - argvals[D - 1]
      )
    },
    midpoint = c(0, diff(argvals))
  )
  
  return(weights)
}


# ---------------------------
# 4. Univariate FPCA from covariance surface
# ---------------------------

uFPCA <- function(G, argvals, pev = NULL, efunctions_sign = NULL) {
  # Computes eigenvalues/eigenfunctions of a covariance surface G
  # using numerical quadrature weights.
  #
  # Args:
  #   G: covariance matrix evaluated on a common grid
  #   argvals: grid points
  #   pev: proportion of explained variance used to select the number of PCs
  #   efunctions_sign: optional matrix used to align eigenfunction signs
  #
  # Returns:
  #   A list with eigenvalues, eigenfunctions, number of PCs, PVE, and weights.
  
  if (!is.matrix(G) || nrow(G) != ncol(G)) {
    stop("'G' must be a square matrix.", call. = FALSE)
  }
  
  if (length(argvals) != nrow(G)) {
    stop("Length of 'argvals' must match the dimension of 'G'.", call. = FALSE)
  }
  
  if (!is.null(pev) && (!is.numeric(pev) || length(pev) != 1L || pev <= 0 || pev > 1)) {
    stop("'pev' must be a single number in (0, 1].", call. = FALSE)
  }
  
  # Numerical integration weights
  w <- quadWeights(argvals, method = "trapezoidal")
  sqrt_w <- sqrt(w)
  
  # Weighted eigenproblem
  V <- G * outer(sqrt_w, sqrt_w)
  eig <- eigen(V, symmetric = TRUE)
  
  evalues <- eig$values
  evalues[evalues <= 0] <- 0
  
  positive_idx <- which(evalues > 0)
  if (length(positive_idx) == 0L) {
    stop("No positive eigenvalues found in covariance operator.", call. = FALSE)
  }
  
  cumpev <- cumsum(evalues[positive_idx]) / sum(evalues[positive_idx])
  
  npc <- length(positive_idx)
  if (!is.null(pev)) {
    npc <- which(cumpev >= pev)[1L]
  } else {
    pev <- cumpev[npc]
  }
  
  efunctions <- matrix(
    1 / sqrt_w,
    nrow = length(sqrt_w),
    ncol = npc
  ) * eig$vectors[, seq_len(npc), drop = FALSE]
  
  # Optional sign alignment
  if (!is.null(efunctions_sign)) {
    if (!is.matrix(efunctions_sign)) {
      stop("'efunctions_sign' must be a matrix.", call. = FALSE)
    }
    if (nrow(efunctions_sign) != nrow(efunctions)) {
      stop("Wrong number of rows in 'efunctions_sign'.", call. = FALSE)
    }
    if (ncol(efunctions_sign) < npc) {
      stop("Wrong number of columns in 'efunctions_sign'.", call. = FALSE)
    }
    
    ref_sign <- efunctions_sign[, seq_len(npc), drop = FALSE]
    sgn <- apply(efunctions * ref_sign, 2, function(x) sign(sum(w * x)))
    sgn[sgn == 0] <- 1
    efunctions <- sweep(efunctions, MARGIN = 2L, STATS = sgn, FUN = "*")
  }
  
  out <- list(
    evalues = evalues[seq_len(npc)],
    efunctions = efunctions,
    npc = npc,
    pev = pev,
    cumpev = cumpev,
    w = w
  )
  
  return(out)
}


# ---------------------------
# 5. Mean and covariance estimation with missing values
#    (Kraus-style empirical estimators)
# ---------------------------

meanKraus <- function(X_mat) {
  # Row-wise mean ignoring missing values
  if (!is.matrix(X_mat)) {
    stop("'X_mat' must be a matrix.", call. = FALSE)
  }
  rowMeans(X_mat, na.rm = TRUE)
}

covKraus <- function(X_mat) {
  # Empirical covariance estimator with missing values.
  # Rows = grid points / features, columns = curves / subjects.
  
  if (!is.matrix(X_mat)) {
    stop("'X_mat' must be a matrix.", call. = FALSE)
  }
  
  X_centered <- X_mat - meanKraus(X_mat)
  missing_idx <- is.na(X_centered)
  
  # Number of valid contributions for each covariance entry
  valid_counts <- tcrossprod(!missing_idx)
  
  # Replace NA by zero before cross-product
  X_centered[missing_idx] <- 0
  prod_mat <- tcrossprod(X_centered)
  
  cov_mat <- prod_mat / valid_counts
  cov_mat[valid_counts == 0] <- NA_real_
  
  return(cov_mat)
}


# ---------------------------
# 6. Curve reconstruction using Kraus-type approach
# ---------------------------

krauss_gia <- function(X_mat, alpha = NULL, Sgm = NULL) {
  # Reconstruct partially observed curves using conditional expectation
  # under a covariance-based regularization approach.
  #
  # Args:
  #   X_mat: matrix with rows = grid points, columns = curves
  #   alpha: optional regularization parameter; if NULL, selected by GCV
  #   Sgm: optional covariance matrix; if NULL, estimated from X_mat
  #
  # Returns:
  #   A list containing reconstructed curves and tuning information.
  
  if (!is.matrix(X_mat)) {
    stop("'X_mat' must be a matrix.", call. = FALSE)
  }
  
  if (is.null(Sgm)) {
    cov_mat <- covKraus(X_mat)
    mean_vec <- meanKraus(X_mat)
  } else {
    if (!is.matrix(Sgm) || nrow(Sgm) != ncol(Sgm)) {
      stop("'Sgm' must be a square matrix.", call. = FALSE)
    }
    if (nrow(Sgm) != nrow(X_mat)) {
      stop("nrow('Sgm') must match nrow('X_mat').", call. = FALSE)
    }
    cov_mat <- Sgm
    mean_vec <- rep(0, nrow(Sgm))
  }
  
  n_curves <- ncol(X_mat)
  
  # Indices of partially observed and fully observed curves
  reconst_fcts <- which(apply(is.na(X_mat), 2, any))
  NonNA_fcts <- setdiff(seq_len(n_curves), reconst_fcts)
  
  if (length(reconst_fcts) == 0L) {
    return(list(
      X_reconst = X_mat,
      alpha = numeric(0),
      df = numeric(0),
      id_obs = NonNA_fcts,
      id_pobs = integer(0)
    ))
  }
  
  X_reconst_mat <- X_mat[, reconst_fcts, drop = FALSE] - mean_vec
  X_Compl_mat <- X_mat[, NonNA_fcts, drop = FALSE]
  X_cent <- X_Compl_mat - mean_vec
  
  alpha_vec <- rep(NA_real_, length(reconst_fcts))
  df_vec <- rep(NA_real_, length(reconst_fcts))
  
  for (i in seq_along(reconst_fcts)) {
    X_tmp <- X_reconst_mat[, i]
    
    M_bool_vec <- is.na(X_tmp)
    O_bool_vec <- !M_bool_vec
    p_obs <- sum(O_bool_vec)
    
    covMO_mat <- cov_mat[M_bool_vec, O_bool_vec, drop = FALSE]
    covOO_mat <- cov_mat[O_bool_vec, O_bool_vec, drop = FALSE]
    
    eig <- eigen(covOO_mat, symmetric = TRUE)
    lambda <- pmax(eig$values, 0)
    Q <- eig$vectors
    
    if (is.null(alpha)) {
      temp <- stats::optimize(
        f = function(alpha_val) {
          positive_lambda <- lambda[lambda > 0]
          df <- sum(positive_lambda / (positive_lambda + alpha_val))
          
          covOO_a_mat_inv <- tcrossprod(
            Q * rep(sqrt(1 / (lambda + alpha_val)), each = p_obs)
          )
          B <- covMO_mat %*% covOO_a_mat_inv
          
          X_M_fit_cent_vec <- B %*% X_cent[O_bool_vec, , drop = FALSE]
          X_cent_tmp <- X_cent
          X_cent_tmp[M_bool_vec, ] <- X_M_fit_cent_vec
          
          X_fit <- X_cent_tmp + mean_vec
          rss <- sum((X_fit[M_bool_vec, ] - X_Compl_mat[M_bool_vec, ])^2)
          gcv <- rss / ((1 - df / ncol(X_Compl_mat))^2)
          
          attr(gcv, "df") <- df
          attr(gcv, "B") <- B
          gcv
        },
        interval = c(1e-12, sum(diag(cov_mat), na.rm = TRUE) * nrow(X_mat)),
        maximum = FALSE
      )
      
      alpha_vec[i] <- temp$minimum
      df_vec[i] <- attr(temp$objective, "df")
      B <- attr(temp$objective, "B")
    } else {
      alpha_vec[i] <- alpha
      df_vec[i] <- sum(lambda[lambda > 0] / (lambda[lambda > 0] + alpha_vec[i]))
      
      covOO_a_mat_inv <- tcrossprod(
        Q * rep(sqrt(1 / (lambda + alpha_vec[i])), each = p_obs)
      )
      B <- covMO_mat %*% covOO_a_mat_inv
    }
    
    X_tmp_M_cent_vec <- B %*% X_tmp[O_bool_vec]
    X_tmp[M_bool_vec] <- X_tmp_M_cent_vec
    X_reconst_mat[, i] <- X_tmp + mean_vec
  }
  
  X_out <- X_mat
  X_out[, reconst_fcts] <- X_reconst_mat
  
  return(list(
    X_reconst = X_out,
    alpha = alpha_vec,
    df = df_vec,
    id_obs = NonNA_fcts,
    id_pobs = reconst_fcts
  ))
}


# ---------------------------
# 7. Graph recovery metrics
# ---------------------------

auc_approx <- function(x, y) {
  # Trapezoidal approximation to AUC
  if (length(x) != length(y)) {
    stop("'x' and 'y' must have the same length.", call. = FALSE)
  }
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * zoo::rollmean(y, 2))
}

compute_metrics_graph <- function(Theta.true, Theta.hat) {
  # Compute ROC- and PR-related metrics across a list of estimated graphs.
  #
  # Args:
  #   Theta.true: true precision matrix
  #   Theta.hat: list of estimated precision matrices across thresholds
  #
  # Returns:
  #   A list of vectors and summary AUC values.
  
  if (!is.matrix(Theta.true) || nrow(Theta.true) != ncol(Theta.true)) {
    stop("'Theta.true' must be a square matrix.", call. = FALSE)
  }
  if (!is.list(Theta.hat) || length(Theta.hat) == 0L) {
    stop("'Theta.hat' must be a non-empty list of matrices.", call. = FALSE)
  }
  
  n_thresh <- length(Theta.hat)
  idx <- upper.tri(Theta.true, diag = FALSE)
  true_edges <- Theta.true[idx] != 0
  
  tpr <- numeric(n_thresh)
  fpr <- numeric(n_thresh)
  precision <- numeric(n_thresh)
  recall <- numeric(n_thresh)
  F1score <- numeric(n_thresh)
  
  for (k in seq_len(n_thresh)) {
    pred_edges <- Theta.hat[[k]][idx] != 0
    
    tp <- sum(true_edges & pred_edges)
    fp <- sum(!true_edges & pred_edges)
    fn <- sum(true_edges & !pred_edges)
    tn <- sum(!true_edges & !pred_edges)
    
    tpr[k] <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    fpr[k] <- if ((fp + tn) > 0) fp / (fp + tn) else 0
    precision[k] <- if ((tp + fp) > 0) tp / (tp + fp) else 0
    recall[k] <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    
    F1score[k] <- if ((precision[k] + recall[k]) > 0) {
      2 * precision[k] * recall[k] / (precision[k] + recall[k])
    } else {
      0
    }
  }
  
  return(list(
    TPR = tpr,
    FPR = fpr,
    precision = precision,
    recall = recall,
    F1score = F1score,
    AUC_ROC = auc_approx(fpr, tpr),
    AUC_PR = auc_approx(recall, precision)
  ))
}

relative_theta_error_layer <- function(theta_hat, theta_true) {
  # Normalized Frobenius error for one precision matrix.
  #
  # The normalization is the number of nonzero entries in theta_true.
  
  if (!is.matrix(theta_hat) || !is.matrix(theta_true)) {
    stop("'theta_hat' and 'theta_true' must be matrices.", call. = FALSE)
  }
  if (!all(dim(theta_hat) == dim(theta_true))) {
    stop("'theta_hat' and 'theta_true' must have the same dimensions.", call. = FALSE)
  }
  
  denom <- sum(theta_true != 0)
  if (denom == 0) {
    stop("The true precision matrix has no nonzero entries.", call. = FALSE)
  }
  
  norm(theta_hat - theta_true, type = "F") / denom
}

relative_theta_error_layer_fro <- function(theta_hat, theta_true) {
  # Relative Frobenius error for one precision matrix.
  
  if (!is.matrix(theta_hat) || !is.matrix(theta_true)) {
    stop("'theta_hat' and 'theta_true' must be matrices.", call. = FALSE)
  }
  if (!all(dim(theta_hat) == dim(theta_true))) {
    stop("'theta_hat' and 'theta_true' must have the same dimensions.", call. = FALSE)
  }
  
  denom <- norm(theta_true, type = "F")
  if (denom == 0) {
    stop("The Frobenius norm of 'theta_true' is zero.", call. = FALSE)
  }
  
  norm(theta_hat - theta_true, type = "F") / denom
}

mean_theta_error_array <- function(theta_hat_array,
                                   theta_true_array,
                                   normalize = c("fro", "nonzero")) {
  # Mean layer-wise error between estimated and true 3D arrays of precision matrices.
  #
  # Args:
  #   theta_hat_array: array p x p x K_hat
  #   theta_true_array: array p x p x K_true
  #   normalize: "nonzero" or "fro"
  #
  # Returns:
  #   Mean error across the K_true layers. If K_hat < K_true, missing layers
  #   are treated as zero matrices.
  
  normalize <- match.arg(normalize)
  
  if (length(dim(theta_hat_array)) != 3L || length(dim(theta_true_array)) != 3L) {
    stop("Both inputs must be 3D arrays.", call. = FALSE)
  }
  
  if (!all(dim(theta_hat_array)[1:2] == dim(theta_true_array)[1:2])) {
    stop("The first two dimensions of the arrays must match.", call. = FALSE)
  }
  
  p <- dim(theta_true_array)[1]
  K_true <- dim(theta_true_array)[3]
  K_hat <- dim(theta_hat_array)[3]
  
  error_fun <- switch(
    normalize,
    nonzero = relative_theta_error_layer,
    fro = relative_theta_error_layer_fro
  )
  
  errors <- numeric(K_true)
  
  for (k in seq_len(K_true)) {
    theta_true_k <- theta_true_array[, , k]
    
    theta_hat_k <- if (k <= K_hat) {
      theta_hat_array[, , k]
    } else {
      matrix(0, nrow = p, ncol = p)
    }
    
    errors[k] <- error_fun(theta_hat = theta_hat_k, theta_true = theta_true_k)
  }
  
  mean(errors)
}

compute_theta_error_path <- function(theta_hat_list,
                                     theta_true_array,
                                     normalize = c("fro", "nonzero")) {
  # Compute mean theta error along a path of estimates.
  
  normalize <- match.arg(normalize)
  
  if (!is.list(theta_hat_list) || length(theta_hat_list) == 0L) {
    stop("'theta_hat_list' must be a non-empty list.", call. = FALSE)
  }
  
  sapply(theta_hat_list, function(theta_hat_array) {
    mean_theta_error_array(
      theta_hat_array = theta_hat_array,
      theta_true_array = theta_true_array,
      normalize = normalize
    )
  })
}

compute_theta_error_path_pofggm <- function(fit_list,
                                            theta_true_array,
                                            normalize = c("fro", "nonzero")) {
  # Compute mean theta error along a path of pofggm fits.
  
  normalize <- match.arg(normalize)
  
  if (!is.list(fit_list) || length(fit_list) == 0L) {
    stop("'fit_list' must be a non-empty list.", call. = FALSE)
  }
  
  sapply(fit_list, function(fit) {
    mean_theta_error_array(
      theta_hat_array = fit$Tht,
      theta_true_array = theta_true_array,
      normalize = normalize
    )
  })
}

aggregate_theta_support <- function(theta_array, use_abs = TRUE) {
  # Aggregate a 3D array of precision matrices into a single matrix
  # by summing across layers.
  
  if (length(dim(theta_array)) != 3L) {
    stop("'theta_array' must be a 3D array.", call. = FALSE)
  }
  
  if (use_abs) {
    theta_array <- abs(theta_array)
  }
  
  Reduce("+", array2list(theta_array))
}

compute_auc_graph_path <- function(theta_true_array, theta_hat_list, use_abs = TRUE) {
  # Compute graph-recovery metrics after aggregating layer-specific precision matrices.
  
  theta_true_agg <- aggregate_theta_support(theta_true_array, use_abs = use_abs)
  
  theta_hat_agg_list <- lapply(theta_hat_list, function(theta_hat_array) {
    aggregate_theta_support(theta_hat_array, use_abs = use_abs)
  })
  
  compute_metrics_graph(
    Theta.true = theta_true_agg,
    Theta.hat = theta_hat_agg_list
  )
}

compute_auc_graph_path_pofggm <- function(theta_true_array, fit_list, use_abs = TRUE) {
  # Compute graph-recovery metrics for a list of pofggm fits.
  
  theta_true_agg <- aggregate_theta_support(theta_true_array, use_abs = use_abs)
  
  theta_hat_agg_list <- lapply(fit_list, function(fit) {
    aggregate_theta_support(fit$Tht, use_abs = use_abs)
  })
  
  compute_metrics_graph(
    Theta.true = theta_true_agg,
    Theta.hat = theta_hat_agg_list
  )
}

relative_curve_error <- function(X_hat, X_true, grid, relative = TRUE) {
  # Integrated reconstruction error for a set of curves on a common grid.
  #
  # Args:
  #   X_hat: estimated curves, matrix d x n
  #   X_true: true curves, matrix d x n
  #   grid: vector of grid points of length d
  #   relative: if TRUE, divide by the integrated squared norm of X_true
  #
  # Returns:
  #   Mean integrated error across columns.
  
  if (!is.matrix(X_hat) || !is.matrix(X_true)) {
    stop("'X_hat' and 'X_true' must be matrices.", call. = FALSE)
  }
  if (!all(dim(X_hat) == dim(X_true))) {
    stop("'X_hat' and 'X_true' must have the same dimensions.", call. = FALSE)
  }
  if (length(grid) != nrow(X_hat)) {
    stop("Length of 'grid' must match the number of rows of 'X_hat'.", call. = FALSE)
  }
  
  sq_err <- (X_hat - X_true)^2
  num <- 0.5 * crossprod(diff(grid), sq_err[-1, , drop = FALSE] + sq_err[-nrow(sq_err), , drop = FALSE])
  
  num <- as.numeric(num)
  
  if (!relative) {
    return(mean(num))
  }
  
  sq_true <- X_true^2
  den <- 0.5 * crossprod(diff(grid), sq_true[-1, , drop = FALSE] + sq_true[-nrow(sq_true), , drop = FALSE])
  den <- as.numeric(den)
  
  valid <- den > 0
  if (!any(valid)) {
    stop("All denominators are zero in relative curve error.", call. = FALSE)
  }
  
  mean(num[valid] / den[valid])
}

compute_imputation_error <- function(X_hat,
                                     X_true,
                                     X_obs,
                                     grid,
                                     relative = TRUE) {
  # Compute mean reconstruction error over partially observed curves only.
  #
  # Args:
  #   X_hat: estimated array d x n x p
  #   X_true: true array d x n x p
  #   X_obs: observed array d x n x p with NAs
  #   grid: time grid
  #   relative: if TRUE, compute relative integrated squared error
  #
  # Returns:
  #   Mean reconstruction error across variables and partially observed subjects.
  
  if (length(dim(X_hat)) != 3L || length(dim(X_true)) != 3L || length(dim(X_obs)) != 3L) {
    stop("'X_hat', 'X_true', and 'X_obs' must be 3D arrays.", call. = FALSE)
  }
  
  if (!all(dim(X_hat) == dim(X_true)) || !all(dim(X_hat) == dim(X_obs))) {
    stop("'X_hat', 'X_true', and 'X_obs' must have the same dimensions.", call. = FALSE)
  }
  
  p <- dim(X_hat)[3]
  errors <- numeric(p)
  
  for (j in seq_len(p)) {
    reconst_fcts <- which(apply(is.na(X_obs[, , j]), 2, any))
    
    if (length(reconst_fcts) == 0L) {
      errors[j] <- NA_real_
      next
    }
    
    X_hat_j <- X_hat[, reconst_fcts, j, drop = TRUE]
    X_true_j <- X_true[, reconst_fcts, j, drop = TRUE]
    
    # If only one partially observed curve is present, keep matrix structure
    if (is.vector(X_hat_j)) {
      X_hat_j <- matrix(X_hat_j, ncol = 1)
    }
    if (is.vector(X_true_j)) {
      X_true_j <- matrix(X_true_j, ncol = 1)
    }
    
    errors[j] <- relative_curve_error(
      X_hat = X_hat_j,
      X_true = X_true_j,
      grid = grid,
      relative = relative
    )
  }
  
  mean(errors, na.rm = TRUE)
}

# ------------------------------------------------------------
# Relative reconstruction error on the missing portion only
# ------------------------------------------------------------

relative_curve_error_missing <- function(x_hat,
                                         x_true,
                                         missing_idx,
                                         grid,
                                         relative = TRUE) {
  # Compute reconstruction error for a single curve on the missing region only.
  #
  # Args:
  #   x_hat: numeric vector of estimated values on the full grid
  #   x_true: numeric vector of true values on the full grid
  #   missing_idx: logical vector; TRUE where the curve is missing
  #   grid: numeric vector of grid points
  #   relative: if TRUE, divide by the squared norm of x_true on the missing region
  #
  # Returns:
  #   Integrated squared error (or relative integrated squared error)
  #   on the missing portion only.
  
  if (!is.numeric(x_hat) || !is.numeric(x_true)) {
    stop("'x_hat' and 'x_true' must be numeric vectors.", call. = FALSE)
  }
  if (length(x_hat) != length(x_true)) {
    stop("'x_hat' and 'x_true' must have the same length.", call. = FALSE)
  }
  if (!is.logical(missing_idx) || length(missing_idx) != length(x_hat)) {
    stop("'missing_idx' must be a logical vector with the same length as 'x_hat'.", call. = FALSE)
  }
  if (!is.numeric(grid) || length(grid) != length(x_hat)) {
    stop("'grid' must be a numeric vector with the same length as 'x_hat'.", call. = FALSE)
  }
  
  idx <- which(missing_idx)
  if (length(idx) == 0L) {
    return(NA_real_)
  }
  
  # If only one missing grid point is available, use pointwise squared error
  if (length(idx) == 1L) {
    num <- (x_hat[idx] - x_true[idx])^2
    
    if (!relative) {
      return(as.numeric(num))
    }
    
    den <- x_true[idx]^2
    if (den <= 0) {
      return(NA_real_)
    }
    
    return(as.numeric(num / den))
  }
  
  x_hat_mis <- x_hat[idx]
  x_true_mis <- x_true[idx]
  grid_mis <- grid[idx]
  
  sq_err <- (x_hat_mis - x_true_mis)^2
  num <- 0.5 * sum(diff(grid_mis) * (sq_err[-1] + sq_err[-length(sq_err)]))
  
  if (!relative) {
    return(as.numeric(num))
  }
  
  sq_true <- x_true_mis^2
  den <- 0.5 * sum(diff(grid_mis) * (sq_true[-1] + sq_true[-length(sq_true)]))
  
  if (den <= 0) {
    return(NA_real_)
  }
  
  as.numeric(num / den)
}


compute_imputation_error_missing <- function(X_hat,
                                             X_true,
                                             X_obs,
                                             grid,
                                             relative = TRUE) {
  # Compute mean reconstruction error over the missing portions only.
  #
  # Args:
  #   X_hat: estimated array d x n x p
  #   X_true: true array d x n x p
  #   X_obs: observed array d x n x p with NAs
  #   grid: time grid of length d
  #   relative: if TRUE, compute relative integrated squared error
  #
  # Returns:
  #   Mean reconstruction error over partially observed curves and variables.
  
  if (length(dim(X_hat)) != 3L || length(dim(X_true)) != 3L || length(dim(X_obs)) != 3L) {
    stop("'X_hat', 'X_true', and 'X_obs' must be 3D arrays.", call. = FALSE)
  }
  
  if (!all(dim(X_hat) == dim(X_true)) || !all(dim(X_hat) == dim(X_obs))) {
    stop("'X_hat', 'X_true', and 'X_obs' must have the same dimensions.", call. = FALSE)
  }
  
  d <- dim(X_hat)[1]
  n <- dim(X_hat)[2]
  p <- dim(X_hat)[3]
  
  if (length(grid) != d) {
    stop("Length of 'grid' must match the first dimension of 'X_hat'.", call. = FALSE)
  }
  
  error_values <- c()
  
  for (j in seq_len(p)) {
    for (i in seq_len(n)) {
      missing_idx <- is.na(X_obs[, i, j])
      
      if (!any(missing_idx)) {
        next
      }
      
      err_ij <- relative_curve_error_missing(
        x_hat = X_hat[, i, j],
        x_true = X_true[, i, j],
        missing_idx = missing_idx,
        grid = grid,
        relative = relative
      )
      
      error_values <- c(error_values, err_ij)
    }
  }
  
  mean(error_values, na.rm = TRUE)
}

prepare_theta_path_input <- function(X_array, Phi, grid) {
  # Prepare inputs for the precision-matrix estimation path.
  #
  # Args:
  #   X_array: array d x n x p containing curves
  #   Phi: matrix/array of basis functions used for projection
  #   grid: vector of time points
  #
  # Returns:
  #   A list with:
  #     Xi      = projected scores
  #     S       = empirical covariance array
  #     rho_max = maximum regularization value
  #     Tht_init = initial precision array obtained from solve(S_k)
  
  Xi <- try(integrate_cube(X_array, Phi, grid), silent = TRUE)
  if (inherits(Xi, "try-error")) {
    return(list(
      ok = FALSE,
      stage = "integrate_cube",
      error = Xi
    ))
  }
  
  tmp <- try(compute_S_and_rho_max(Xi), silent = TRUE)
  if (inherits(tmp, "try-error")) {
    return(list(
      ok = FALSE,
      stage = "compute_S_and_rho_max",
      error = tmp
    ))
  }
  
  S <- tmp$S
  rho_max <- tmp$rho.max
  
  K_hat <- dim(Xi)[3]
  p <- dim(S)[1]
  
  Tht_init <- array(0, dim = c(p, p, K_hat))
  for (k in seq_len(K_hat)) {
    Tht_init[, , k] <- 1 / diag(S[, , k])
  }
  
  return(list(
    ok = TRUE,
    Xi = Xi,
    S = S,
    rho_max = rho_max,
    Tht_init = Tht_init
  ))
}

run_theta_path <- function(S,
                           rho_max,
                           perc_rho,
                           wTht,
                           Tht_init = NULL,
                           alpha = 0,
                           pendiag = FALSE,
                           maxit = 1e5,
                           thr = 1e-6,
                           trace = 0,
                           verbose = FALSE) {
  # Compute a path of precision-matrix estimates over a rho grid.
  #
  # Args:
  #   S: array p x p x K of empirical covariance matrices
  #   rho_max: maximum value of rho
  #   perc_rho: vector of multipliers in [0, 1]
  #   wTht: weights used in the ADMM optimization
  #   Tht_init: optional initial array for Theta
  #   alpha, pendiag, maxit, thr, trace: ADMM tuning parameters
  #   verbose: logical; print progress if TRUE
  #
  # Returns:
  #   A list with:
  #     grid_rho  = rho grid
  #     theta_path = list of estimated precision arrays
  #     converged = logical
  #     last_fit  = last successful fit
  
  if (length(dim(S)) != 3L) {
    stop("'S' must be a 3D array.", call. = FALSE)
  }
  
  p <- dim(S)[1]
  K_hat <- dim(S)[3]
  
  if (is.null(Tht_init)) {
    Tht_init <- array(0, dim = c(p, p, K_hat))
    for (k in seq_len(K_hat)) {
      Tht_init[, , k] <- solve(S[, , k])
    }
  }
  
  grid_rho <- rho_max * perc_rho
  nrho <- length(grid_rho)
  
  theta_path <- vector("list", length = nrho)
  current_Tht <- Tht_init
  
  for (r in seq_len(nrho)) {
    if (verbose) {
      cat("\nRho iteration:", r, "of", nrho, "\n")
    }
    
    fit <- try(
      admm_tht_sub(
        p = p,
        N = K_hat,
        fk = rep(1 / K_hat, K_hat),
        S = S,
        wTht = wTht,
        pendiag = pendiag,
        rho = grid_rho[r],
        alpha = alpha,
        maxit = maxit,
        thr = thr,
        Tht = current_Tht,
        trace = trace
      ),
      silent = TRUE
    )
    
    if (inherits(fit, "try-error")) {
      return(list(
        ok = FALSE,
        grid_rho = grid_rho,
        theta_path = theta_path,
        converged = FALSE,
        failed_at = r,
        error = fit
      ))
    }
    
    current_Tht <- fit$Tht
    theta_path[[r]] <- current_Tht
  }
  
  return(list(
    ok = TRUE,
    grid_rho = grid_rho,
    theta_path = theta_path,
    converged = TRUE,
    last_fit = current_Tht
  ))
}

estimate_theta_path_from_curves <- function(X_array,
                                            Phi,
                                            grid,
                                            perc_rho,
                                            wTht,
                                            alpha = 0,
                                            pendiag = FALSE,
                                            maxit = 1e5,
                                            thr = 1e-6,
                                            trace = 0,
                                            verbose = FALSE) {
  # Full wrapper: project curves, build covariance objects, and estimate
  # a path of precision matrices over a rho grid.
  #
  # Returns:
  #   A list with Xi, S, rho_max, grid_rho, Tht_init, and theta_path.
  
  prep <- prepare_theta_path_input(
    X_array = X_array,
    Phi = Phi,
    grid = grid
  )
  
  if (!prep$ok) {
    return(prep)
  }
  
  fit_path <- run_theta_path(
    S = prep$S,
    rho_max = prep$rho_max,
    perc_rho = perc_rho,
    wTht = wTht,
    Tht_init = prep$Tht_init,
    alpha = alpha,
    pendiag = pendiag,
    maxit = maxit,
    thr = thr,
    trace = trace,
    verbose = verbose
  )
  
  if (!fit_path$ok) {
    return(c(prep, fit_path))
  }
  
  return(list(
    ok = TRUE,
    Xi = prep$Xi,
    S = prep$S,
    rho_max = prep$rho_max,
    Tht_init = prep$Tht_init,
    grid_rho = fit_path$grid_rho,
    theta_path = fit_path$theta_path
  ))
}

# ---------------------------
# 3. Helper: parse one configuration row
# ---------------------------

parse_sim_config <- function(config_row) {
  # Convert one row of 'configs' into a named list of parameters.
  
  n <- as.integer(config_row[["n"]])
  p_over_n <- as.numeric(config_row[["p_over_n"]])
  p <- as.integer(round(n * p_over_n))
  
  K_true <- as.integer(config_row[["K_true"]])
  perc_window <- as.numeric(config_row[["perc_window"]])
  perc_obs_curves <- as.numeric(config_row[["perc_obs_curves"]])
  graph_type <- as.character(config_row[["graph_type"]])
  pev <- as.numeric(config_row[["pev"]])
  perc_theta_share <- as.numeric(config_row[["perc_theta_share"]])
  scenario <- as.character(config_row[["scenario"]])
  
  d <- as.numeric(config_row[["d"]])
  
  # Width of the missing window on the observed grid
  window_size <- as.integer(round(d * perc_window))
  
  list(
    scenario = scenario,
    n = n,
    p = p,
    p_over_n = p_over_n,
    d = d,
    K_true = K_true,
    perc_window = perc_window,
    window_size = window_size,
    perc_obs_curves = perc_obs_curves,
    graph_type = graph_type,
    pev = pev,
    perc_theta_share = perc_theta_share
  )
}


# ---------------------------
# 4. Helper: initialize storage objects
# ---------------------------

initialize_simulation_storage <- function(d, n, p, K_true, n_sim, n_rho) {
  # Allocate arrays, matrices, and vectors used to store simulation outputs.
  
  X <- array(
    0,
    dim = c(d, n, p),
    dimnames = list(
      paste0("tp", seq_len(d)),
      seq_len(n),
      paste0("Var", seq_len(p))
    )
  )
  
  Xi <- array(
    0,
    dim = c(n, p, K_true),
    dimnames = list(
      seq_len(n),
      paste0("Score", seq_len(p)),
      paste0("k", seq_len(K_true))
    )
  )
  
  Xo_save <- Xpo_save <- array(
    0,
    dim = c(d, n, p, n_sim),
    dimnames = list(
      paste0("tp", seq_len(d)),
      seq_len(n),
      paste0("Var", seq_len(p)),
      seq_len(n_sim)
    )
  )
  
  Xi_save <- array(
    0,
    dim = c(n, p, K_true, n_sim),
    dimnames = list(
      seq_len(n),
      paste0("Score", seq_len(p)),
      paste0("k", seq_len(K_true)),
      seq_len(n_sim)
    )
  )
  
  # Precision-matrix estimation metrics
  theta_err_mat <- matrix(NA_real_, nrow = n_rho, ncol = n_sim)
  theta_err_kraus_mat <- matrix(NA_real_, nrow = n_rho, ncol = n_sim)
  theta_err_obs_mat <- matrix(NA_real_, nrow = n_rho, ncol = n_sim)
  
  # Graph-recovery metrics
  auc_theta_vec <- numeric(n_sim)
  auc_theta_kraus_vec <- numeric(n_sim)
  auc_theta_obs_vec <- numeric(n_sim)
  
  # Curve reconstruction metrics
  curve_err_mat <- matrix(NA_real_, nrow = n_rho, ncol = n_sim)
  curve_err_kraus_vec <- numeric(n_sim)
  
  # Computational time
  comp_time_vec <- numeric(n_sim)
  
  list(
    X = X,
    Xi = Xi,
    Xo_save = Xo_save,
    Xpo_save = Xpo_save,
    Xi_save = Xi_save,
    theta_err_mat = theta_err_mat,
    theta_err_kraus_mat = theta_err_kraus_mat,
    theta_err_obs_mat = theta_err_obs_mat,
    auc_theta_vec = auc_theta_vec,
    auc_theta_kraus_vec = auc_theta_kraus_vec,
    auc_theta_obs_vec = auc_theta_obs_vec,
    curve_err_mat = curve_err_mat,
    curve_err_kraus_vec = curve_err_kraus_vec,
    comp_time_vec = comp_time_vec
  )
}


# ---------------------------
# 5. Helper: create basis objects
# ---------------------------

create_basis_objects <- function(tp, K_true) {
  # Create the Fourier basis and evaluate it on the time grid.
  
  basis_obj <- create.fourier.basis(nbasis = K_true)
  Phi <- eval.basis(evalarg = tp, basisobj = basis_obj)
  rownames(Phi) <- paste0("tp", seq_along(tp))
  
  list(
    basis_obj = basis_obj,
    Phi = Phi
  )
}


# ---------------------------
# 6. Helper: generate true precision/covariance layers
# ---------------------------

# ============================================================
# Helper functions for simulation design
# ============================================================

# ---------------------------
# 1. Convert between lists and arrays
# ---------------------------

list2array <- function(x) {
  # Convert a list of matrices or vectors into an array/matrix.
  #
  # Args:
  #   x: list of matrices with common dimensions, or list of vectors
  #
  # Returns:
  #   If x is a list of matrices, returns a 3D array.
  #   If x is a list of vectors, returns a matrix.
  
  if (!is.list(x) || length(x) == 0L) {
    stop("'x' must be a non-empty list.", call. = FALSE)
  }
  
  first_elem <- x[[1]]
  
  if (is.matrix(first_elem)) {
    n_rows <- vapply(x, nrow, integer(1))
    n_cols <- vapply(x, ncol, integer(1))
    
    if (length(unique(n_cols)) != 1L) {
      stop("All matrices in 'x' must have the same number of columns.", call. = FALSE)
    }
    
    K <- length(x)
    p <- n_cols[1]
    
    out <- array(
      0,
      dim = c(max(n_rows), p, K),
      dimnames = list(
        NULL,
        colnames(first_elem),
        names(x)
      )
    )
    
    for (k in seq_len(K)) {
      out[seq_len(n_rows[k]), , k] <- x[[k]]
    }
    
    return(out)
  }
  
  if (is.atomic(first_elem)) {
    vec_lengths <- vapply(x, length, integer(1))
    if (length(unique(vec_lengths)) != 1L) {
      stop("All vectors in 'x' must have the same length.", call. = FALSE)
    }
    
    out <- matrix(
      unlist(x, use.names = FALSE),
      nrow = vec_lengths[1],
      ncol = length(x),
      dimnames = list(names(first_elem), names(x))
    )
    
    return(out)
  }
  
  stop("'x' must be a list of matrices or vectors.", call. = FALSE)
}


array2list <- function(x) {
  # Convert a 2D/3D array into a list.
  #
  # Args:
  #   x: matrix or array
  #
  # Returns:
  #   A list of vectors or matrices.
  
  if (!is.array(x)) {
    stop("'x' must be an array or matrix.", call. = FALSE)
  }
  
  dim_x <- dim(x)
  if (is.null(dim_x)) {
    stop("'x' must have dimensions.", call. = FALSE)
  }
  
  K <- tail(dim_x, 1)
  out <- vector(mode = "list", length = K)
  
  for (k in seq_len(K)) {
    out[[k]] <- if (length(dim_x) == 2L) {
      x[, k]
    } else {
      matrix(x[, , k], nrow = dim_x[1], ncol = dim_x[2])
    }
  }
  
  last_dimnames <- tail(dimnames(x), 1)[[1]]
  if (!is.null(last_dimnames)) {
    names(out) <- last_dimnames
  }
  
  return(out)
}


# ---------------------------
# 2. Shared support modification across layers
# ---------------------------

apply_shared_support <- function(Tht, perc_theta_share) {
  # Remove a fraction of nonzero upper-triangular entries from a precision matrix
  # to reduce support sharing across layers.
  #
  # Args:
  #   Tht: square symmetric precision matrix
  #   perc_theta_share: proportion of nonzero entries to keep
  #
  # Returns:
  #   Modified symmetric precision matrix.
  
  if (!is.matrix(Tht) || nrow(Tht) != ncol(Tht)) {
    stop("'Tht' must be a square matrix.", call. = FALSE)
  }
  
  if (!is.numeric(perc_theta_share) || length(perc_theta_share) != 1L ||
      perc_theta_share <= 0 || perc_theta_share > 1) {
    stop("'perc_theta_share' must be in (0, 1].", call. = FALSE)
  }
  
  if (perc_theta_share == 1.0) {
    return(Tht)
  }
  
  upper_idx <- upper.tri(Tht)
  nonzero_idx <- which(Tht[upper_idx] != 0)
  
  if (length(nonzero_idx) == 0L) {
    return(Tht)
  }
  
  n_remove <- round((1.0 - perc_theta_share) * length(nonzero_idx))
  if (n_remove == 0L) {
    return(Tht)
  }
  
  idx_remove <- sample(nonzero_idx, size = n_remove, replace = FALSE)
  
  tmp_upper <- Tht[upper_idx]
  tmp_upper[idx_remove] <- 0
  Tht[upper_idx] <- tmp_upper
  
  Tht <- t(Tht)
  tmp_upper <- Tht[upper_idx]
  tmp_upper[idx_remove] <- 0
  Tht[upper_idx] <- tmp_upper
  
  return(Tht)
}


# ---------------------------
# 3. Star-shaped precision matrix generator
# ---------------------------

rStars <- function(d, nstars = 1L, tht.min = 0.4, tht.max = 0.5) {
  # Generate a block-diagonal precision matrix with star-shaped blocks.
  #
  # Args:
  #   d: size of each block
  #   nstars: number of star blocks
  #   tht.min, tht.max: range for off-diagonal edge weights
  #
  # Returns:
  #   A block-diagonal precision matrix.
  
  if (!is.numeric(d) || length(d) != 1L || d < 2) {
    stop("'d' must be an integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(nstars) || length(nstars) != 1L || nstars < 1) {
    stop("'nstars' must be an integer >= 1.", call. = FALSE)
  }
  if (tht.min <= 0 || tht.max <= 0 || tht.min > tht.max) {
    stop("'tht.min' and 'tht.max' must be positive with tht.min <= tht.max.", call. = FALSE)
  }
  
  d <- as.integer(d)
  nstars <- as.integer(nstars)
  
  blocks <- vector(mode = "list", length = nstars)
  
  for (l in seq_len(nstars)) {
    Tht <- diag(d)
    
    edge_weights <- runif(d - 1L, min = tht.min, max = tht.max) *
      sample(c(-1, 1), size = d - 1L, replace = TRUE)
    
    K_mat <- outer(edge_weights, edge_weights)
    Tht_xx <- diag((1 + sqrt(1 + 4 * diag(K_mat))) / 2)
    
    Tht[1, 1] <- 1 + drop(edge_weights %*% solve(Tht_xx) %*% edge_weights)
    Tht[1, -1] <- edge_weights
    Tht[-1, 1] <- edge_weights
    Tht[-1, -1] <- Tht_xx
    
    blocks[[l]] <- Tht
  }
  
  return(do.call(pracma::blkdiag, blocks))
}


# ---------------------------
# 4. Band precision matrix generator
# ---------------------------

rTht <- function(p, id.diag, det.min = 0.5) {
  # Generate a banded precision matrix with entries rho^|i-j|
  # up to a given bandwidth.
  #
  # Args:
  #   p: matrix dimension
  #   id.diag: bandwidth
  #   det.min: target determinant used to calibrate rho
  #
  # Returns:
  #   A banded precision matrix.
  
  if (!is.numeric(p) || length(p) != 1L || p < 2) {
    stop("'p' must be an integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(id.diag) || length(id.diag) != 1L || id.diag < 0 || id.diag > (p - 1)) {
    stop("'id.diag' must be an integer between 0 and p - 1.", call. = FALSE)
  }
  if (!is.numeric(det.min) || length(det.min) != 1L || det.min <= 0) {
    stop("'det.min' must be a positive scalar.", call. = FALSE)
  }
  
  p <- as.integer(p)
  id.diag <- as.integer(id.diag)
  
  build_matrix <- function(rho) {
    outer(seq_len(p), seq_len(p), function(i, j) {
      dist_ij <- abs(i - j)
      ifelse(dist_ij <= id.diag, rho^dist_ij, 0)
    })
  }
  
  opt <- optimize(
    f = function(rho) {
      Tht <- build_matrix(rho)
      (det(Tht) - det.min)^2
    },
    interval = c(0, 1)
  )
  
  rho <- opt$minimum
  Tht <- build_matrix(rho)
  
  return(Tht)
}


rPar <- function(p, K, id.diag, s1 = 3, s2 = -1.8, perc.theta.share = 1.0) {
  # Generate layer-specific banded precision and covariance matrices.
  #
  # Args:
  #   p: dimension of each precision matrix
  #   K: number of layers
  #   id.diag: bandwidth of the banded matrices
  #   s1, s2: scaling parameters across layers
  #   perc.theta.share: proportion of nonzero entries shared across layers
  #
  # Returns:
  #   A list with arrays Tht and Sgm.
  
  if (!is.numeric(p) || length(p) != 1L || p < 2) {
    stop("'p' must be an integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(K) || length(K) != 1L || K < 1) {
    stop("'K' must be an integer >= 1.", call. = FALSE)
  }
  
  p <- as.integer(p)
  K <- as.integer(K)
  
  Tht <- array(
    0,
    dim = c(p, p, K),
    dimnames = list(
      paste0("X", seq_len(p)),
      paste0("X", seq_len(p)),
      paste0("K", seq_len(K))
    )
  )
  
  Sgm <- Tht
  
  for (k in seq_len(K)) {
    theta_k <- rTht(p = p, id.diag = id.diag)
    scale_factor <- s1 * k^s2
    theta_k <- theta_k / scale_factor
    
    if (k != 1L && perc.theta.share < 1.0) {
      theta_k <- apply_shared_support(theta_k, perc_theta_share = perc.theta.share)
    }
    
    Tht[, , k] <- theta_k
    Sgm[, , k] <- solve(theta_k)
  }
  
  return(list(Tht = Tht, Sgm = Sgm))
}


# ---------------------------
# 5. Random graph generator via BDgraph
# ---------------------------

simulTheta <- function(n,
                       p,
                       prob,
                       alpha,
                       vis = FALSE,
                       graph = c("random", "smallworld", "star", "circle",
                                 "cluster", "scale-free", "lattice", "hub",
                                 "AR(1)", "AR(2)"),
                       seed = NULL,
                       ...) {
  # Wrapper around BDgraph::bdgraph.sim to generate graph-based precision matrices.
  #
  # Args:
  #   n: sample size argument passed to bdgraph.sim
  #   p: number of variables
  #   prob: graph sparsity parameter
  #   alpha: graph density multiplier
  #   vis: logical, whether to visualize the graph
  #   graph: graph type
  #   seed: optional random seed
  #   ...: additional arguments passed to bdgraph.sim
  #
  # Returns:
  #   The object returned by BDgraph::bdgraph.sim().
  
  graph <- match.arg(graph)
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (!is.numeric(n) || length(n) != 1L || n < 1) {
    stop("'n' must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(p) || length(p) != 1L || p < 2) {
    stop("'p' must be an integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(prob) || length(prob) != 1L || prob <= 0 || prob > 1) {
    stop("'prob' must be in (0, 1].", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0) {
    stop("'alpha' must be a positive scalar.", call. = FALSE)
  }
  
  out <- BDgraph::bdgraph.sim(
    p = as.integer(p),
    n = as.integer(n),
    prob = prob * alpha,
    vis = vis,
    graph = graph,
    ...
  )
  
  return(out)
}


# ---------------------------
# 6. Full generator of true model
# ---------------------------

generate_true_model <- function(p,
                                K_true,
                                graph_type,
                                perc_theta_share,
                                tht_min,
                                tht_max,
                                s1,
                                s2,
                                seed_base = NULL) {
  # Generate the true layer-specific precision matrices and covariances.
  #
  # Args:
  #   p: number of variables
  #   K_true: number of layers
  #   graph_type: "star", "band", or "smallworld"
  #   perc_theta_share: proportion of nonzero entries shared across layers
  #   tht_min, tht_max: weight range for star graphs
  #   s1, s2: layer scaling parameters
  #   seed_base: optional base seed
  #
  # Returns:
  #   A list with:
  #     Tht = array p x p x K_true of precision matrices
  #     Sgm = array p x p x K_true of covariance matrices
  
  if (!is.null(seed_base)) {
    set.seed(seed_base)
  }
  
  if (!graph_type %in% c("star", "band", "smallworld")) {
    stop("'graph_type' must be one of: 'star', 'band', 'smallworld'.", call. = FALSE)
  }
  
  if (graph_type == "star") {
    n_stars <- 5L
    
    if (p %% n_stars != 0L) {
      stop("For 'star' graphs, 'p' must be divisible by 5.", call. = FALSE)
    }
    
    block_size <- as.integer(p / n_stars)
    
    out_list <- lapply(seq_len(K_true), function(k) {
      theta_k <- rStars(
        d = n_stars,
        nstars = block_size,
        tht.min = tht_min,
        tht.max = tht_max
      )
      
      scale_factor <- s1 * k^(-s2)
      theta_k <- theta_k / scale_factor
      
      if (k != 1L && perc_theta_share < 1.0) {
        theta_k <- apply_shared_support(theta_k, perc_theta_share = perc_theta_share)
      }
      
      sigma_k <- solve(theta_k)
      list(Tht = theta_k, Sgm = sigma_k)
    })
    
    return(list(
      Tht = list2array(lapply(out_list, function(x) x$Tht)),
      Sgm = list2array(lapply(out_list, function(x) x$Sgm))
    ))
  }
  
  if (graph_type == "band") {
    return(
      rPar(
        p = p,
        K = K_true,
        id.diag = 2,
        s1 = s1,
        s2 = -s2,
        perc.theta.share = perc_theta_share
      )
    )
  }
  
  if (graph_type == "smallworld") {
    out_list <- lapply(seq_len(K_true), function(k) {
      theta_k <- round(
        simulTheta(
          n = p,
          p = p,
          prob = 0.1,
          alpha = 1,
          vis = FALSE,
          graph = "smallworld",
          rewire = 0.5,
          seed = if (is.null(seed_base)) NULL else seed_base + k
        )$K,
        5
      )
      
      scale_factor <- s1 * k^(-s2)
      theta_k <- theta_k / scale_factor
      
      if (k != 1L && perc_theta_share < 1.0) {
        theta_k <- apply_shared_support(theta_k, perc_theta_share = perc_theta_share)
      }
      
      sigma_k <- solve(theta_k)
      list(Tht = theta_k, Sgm = sigma_k)
    })
    
    return(list(
      Tht = list2array(lapply(out_list, function(x) x$Tht)),
      Sgm = list2array(lapply(out_list, function(x) x$Sgm))
    ))
  }
}

create_full_matrix <- function(out_rPar, alpha = 0.5) {
  Tht_blocks <- out_rPar$Tht
  Sgm_blocks <- out_rPar$Sgm
  dimT <- dim(Tht_blocks)
  p <- dimT[1]
  K_true <- dimT[3]
  Tht_full <- matrix(0.0, nrow = p * K_true, ncol = p * K_true)
  for(k in 1:K_true){
    Tht_full[1:p + (k-1)*p, 1:p + (k-1)*p] <- Tht_blocks[,,k]
    if(k > 1) {
      Tht_full[1:p + (k-2)*p, 1:p + (k-1)*p] <- 0.005 * ((Tht_blocks[,,k-1] - diag(diag(Tht_blocks[,,k-1]))) + 
                                                         (Tht_blocks[,,k] - diag(diag(Tht_blocks[,,k]))))
      Tht_full[1:p + (k-1)*p, 1:p + (k-2)*p] <- t(Tht_full[1:p + (k-2)*p, 1:p + (k-1)*p])
    }
  }
  Sgm_full <- solve(Tht_full)
  out_rPar$Tht_full <- Tht_full
  out_rPar$Sgm_full <- Sgm_full
  return(out_rPar)
}

simulate_scores_and_curves <- function(Sgm_array, Phi, n, p, K_true) {
  # Simulate latent scores and reconstruct full curves.
  #
  # Args:
  #   Sgm_array: array p x p x K_true of covariance matrices
  #   Phi: basis matrix d x K_true
  #   n: sample size
  #   p: number of functional variables
  #   K_true: number of retained basis components
  #
  # Returns:
  #   A list with:
  #     Xi = array n x p x K_true of scores
  #     X  = array d x n x p of complete curves
  
  Xi <- list2array(
    apply(
      Sgm_array,
      3L,
      function(Sig) MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sig),
      simplify = FALSE
    )
  )
  
  X <- list2array(
    apply(
      Xi,
      2,
      function(xi) tcrossprod(Phi, xi),
      simplify = FALSE
    )
  )
  
  list(Xi = Xi, X = X)
}

generate_block_missing_mask <- function(d, n, p, id_pobs, window_size) {
  # Generate a block-missing mask for partially observed curves.
  #
  # Args:
  #   d: number of grid points
  #   n: number of subjects
  #   p: number of functional variables
  #   id_pobs: indices of partially observed subjects
  #   window_size: size of the missing block
  #
  # Returns:
  #   A logical array d x n x p. TRUE = missing.
  
  if (window_size > d) {
    stop("'window_size' cannot exceed 'd'.", call. = FALSE)
  }
  
  Mask <- array(FALSE, dim = c(d, n, p))
  
  for (i in id_pobs) {
    for (j in seq_len(p)) {
      start_idx <- sample.int(d - window_size + 1L, size = 1L)
      Mask[seq.int(from = start_idx, length.out = window_size), i, j] <- TRUE
    }
  }
  
  return(Mask)
}

estimate_empirical_basis <- function(X, tp, pev, smooth_k = 10) {
  # Estimate empirical basis functions from partially observed curves.
  #
  # Args:
  #   X: array d x n x p with missing values
  #   tp: time grid
  #   pev: proportion of explained variance for uFPCA truncation
  #   smooth_k: basis dimension for mgcv::te smoother
  #
  # Returns:
  #   A list with:
  #     mu_sample
  #     G_sample
  #     G_smooth
  #     H_sample
  #     ufpca
  #     Phi_emp
  
  d <- dim(X)[1]
  p <- dim(X)[3]
  
  mu_sample <- apply(X, 3, meanKraus)
  G_sample <- apply(X, 3, covKraus, simplify = FALSE)
  
  row_vec <- rep(tp, each = d)
  col_vec <- rep(tp, times = d)
  
  G_smooth <- vector("list", length = p)
  H_sample <- matrix(0, nrow = d, ncol = d)
  
  for (j in seq_len(p)) {
    G_smooth[[j]] <- matrix(
      mgcv::gam(
        as.vector(G_sample[[j]]) ~ te(row_vec, col_vec, k = smooth_k)
      )$fitted.values,
      nrow = d,
      ncol = d
    )
    
    G_smooth[[j]] <- (G_smooth[[j]] + t(G_smooth[[j]])) / 2
    H_sample <- H_sample + G_smooth[[j]]
  }
  
  H_sample <- H_sample / p
  
  ufpca_fit <- uFPCA(G = H_sample, argvals = tp, pev = pev)
  Phi_emp <- ufpca_fit$efunctions
  
  list(
    mu_sample = mu_sample,
    G_sample = G_sample,
    G_smooth = G_smooth,
    H_sample = H_sample,
    ufpca = ufpca_fit,
    Phi_emp = Phi_emp
  )
}

reconstruct_curves_kraus <- function(X) {
  # Reconstruct each functional variable separately using krauss_gia().
  #
  # Args:
  #   X: array d x n x p with missing values
  #
  # Returns:
  #   A list with:
  #     result_list
  #     X_reconstructed
  
  p <- dim(X)[3]
  
  result_list <- lapply(seq_len(p), function(j) {
    krauss_gia(X_mat = X[, , j])
  })
  
  X_reconstructed <- list2array(
    lapply(seq_len(p), function(j) result_list[[j]]$X_reconst)
  )
  
  list(
    result_list = result_list,
    X_reconstructed = X_reconstructed
  )
}

sample_partial_subjects <- function(n, perc_obs_curves) {
  # Sample partially observed subjects.
  n_pobs <- as.integer(round(n * perc_obs_curves))
  sort(sample.int(n, size = n_pobs, replace = FALSE))
}
