list2array <- function(x) {
  if(!is.list(x)) stop("x is not a list!")
  K <- length(x)
  if(is.matrix(x[[1]])) {
    n <- sapply(x, nrow)
    p <- ncol(x[[1]])
    y <- array(0, dim = c(max(n), p, K), dimnames = list(NULL, colnames(x[[1]]), names(x)))
    for(k in seq_len(K)){ y[seq_len(n[k]), , k] <- x[[k]] }
  } else {
    n <- length(x[[1]])
    y <- matrix(unlist(x), nrow = n, ncol = K, dimnames = list(names(x[[1]]), names(x)))
  }
  y
}
array2list <- function(x) {
  if(!is.array(x)) stop("x is not a list!")
  dimX <- dim(x)
  K <- tail(dimX, 1)
  y <- vector(mode = "list", length = K)
  for(k in seq_len(K)) y[[k]] <- if(length(dimX) == 2) x[, k] else na.omit(matrix(x[,,k], dimX[1], dimX[2]))
  names(y) <- tail(dimnames(x), 1)[[1]]
  y
}

# ------------------------------------------------------------
# Helper: regularized inverse from eigendecomposition
# ------------------------------------------------------------
compute_Ginv_R <- function(lambda, Q, alpha) {
  # Compute (G + alpha I)^(-1) from the eigendecomposition of G.
  #
  # Args:
  #   lambda: eigenvalues of G
  #   Q: eigenvectors of G
  #   alpha: regularization parameter
  #
  # Returns:
  #   A matrix equal to Q diag(1 / (lambda + alpha)) Q^T
  
  if (!is.numeric(lambda)) {
    stop("'lambda' must be a numeric vector.", call. = FALSE)
  }
  if (!is.matrix(Q)) {
    stop("'Q' must be a matrix.", call. = FALSE)
  }
  if (length(lambda) != ncol(Q)) {
    stop("Length of 'lambda' must match ncol('Q').", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha < 0) {
    stop("'alpha' must be a non-negative scalar.", call. = FALSE)
  }
  
  Q %*% diag(1 / (lambda + alpha), nrow = length(lambda)) %*% t(Q)
}
compute_A_fast_R <- function(alpha, setup) {
  # Compute A(alpha) = (Goo + alpha I)^(-1) Gom
  # using the eigendecomposition precomputed in setup.
  
  inv_diag <- 1 / (setup$lambda + alpha)
  weighted_Qt_Gom <- setup$Qt_Gom * inv_diag
  setup$Q %*% weighted_Qt_Gom
}
# ------------------------------------------------------------
# Helper: objective function for alpha selection
# ------------------------------------------------------------
alpha_objective_R <- function(alpha, lambda, Q, Gom, ZidOi, ZidMi, n_obs) {
  # GCV-like objective used to select alpha.
  #
  # Args:
  #   alpha: regularization parameter
  #   lambda, Q: eigendecomposition of Goo
  #   Gom: cross-covariance block G[Oi, Mi]
  #   ZidOi: observed part of Z for fully observed subjects
  #   ZidMi: missing part of Z for fully observed subjects
  #   n_obs: number of fully observed subjects
  #
  # Returns:
  #   Scalar objective value
  
  positive_lambda <- lambda[lambda > 0]
  df <- sum(positive_lambda / (positive_lambda + alpha))
  
  Ginv <- compute_Ginv_R(lambda = lambda, Q = Q, alpha = alpha)
  A <- Ginv %*% Gom
  
  # Predicted missing block from observed block
  ZidMi_hat <- ZidOi %*% A
  
  rss <- sum((ZidMi_hat - ZidMi)^2)
  gcv <- rss / ((1 - df / n_obs)^2)
  
  gcv
}
alpha_objective_fast <- function(alpha, setup) {
  # Fast GCV-like objective for alpha optimization.
  #
  # Args:
  #   alpha: scalar regularization parameter
  #   setup: object returned by fast_alpha_setup()
  #
  # Returns:
  #   Scalar objective value
  
  lambda <- setup$lambda
  inv_diag <- 1 / (lambda + alpha)
  
  # df(alpha) = sum_j lambda_j / (lambda_j + alpha)
  pos_lambda <- lambda[setup$pos_idx]
  df <- sum(pos_lambda / (pos_lambda + alpha))
  
  # Equivalent to:
  # A = Q %*% diag(inv_diag) %*% (Q^T Gom)
  # ZidMi_hat = ZidOi %*% A
  #
  # but without constructing Ginv or diag()
  weighted_Qt_Gom <- setup$Qt_Gom * inv_diag
  ZidMi_hat <- setup$ZidOi_Q %*% weighted_Qt_Gom
  
  rss <- sum((ZidMi_hat - setup$ZidMi)^2)
  rss / ((1 - df / setup$n_obs)^2)
}
gcv_fun_R <- function(alpha, Goo, Gom, ZidOi, ZidMi, n_obs) {
  R <- chol(Goo + diag(alpha, nrow(Goo)))
  A <- backsolve(R, forwardsolve(t(R), Gom))
  ZidMi_hat <- ZidOi %*% A
  rss <- sum((ZidMi_hat - ZidMi)^2) / n_obs # norm(ZidMi_hat - ZidMi, "F)^2 / n_obs
  return(rss)
  df <- norm(A, "F")^2 # sum(diag(tcrossprod(A)))
  # M_inv <- backsolve(R, forwardsolve(t(R), diag(nrow(Goo))))
  # df <- as.numeric(sum(diag(Goo %*% M_inv)))
  penalty <- 1 - df / n_obs
  rss / (penalty^2)
}

df_sparse_exact <- function(Goo, alpha) {
  # Exact df(alpha) = tr{Goo (Goo + alpha I)^(-1)}
  # using a sparse Cholesky factorization.
  
  if (!inherits(Goo, "Matrix")) {
    Goo <- Matrix::Matrix(Goo, sparse = TRUE)
  }
  
  m <- nrow(Goo)
  M_alpha <- Goo + Matrix::Diagonal(m, x = alpha)
  
  fac <- Matrix::Cholesky(M_alpha, LDL = FALSE, perm = TRUE)
  M_inv <- Matrix::solve(fac, Matrix::Diagonal(m))
  
  as.numeric(sum(diag(Goo %*% M_inv)))
}
alpha_objective_sparse <- function(alpha, Goo, Gom, ZidOi, ZidMi, n_obs) {
  # GCV-like objective based on sparse Cholesky solves.
  #
  # Args:
  #   alpha: regularization parameter
  #   Goo: sparse observed-observed block
  #   Gom: observed-missing block
  #   ZidOi: observed block for fully observed subjects
  #   ZidMi: missing block for fully observed subjects
  #   n_obs: number of fully observed subjects
  #
  # Returns:
  #   Scalar GCV criterion
  
  if (!inherits(Goo, "Matrix")) {
    Goo <- Matrix::Matrix(Goo, sparse = TRUE)
  }
  
  m <- nrow(Goo)
  M_alpha <- Goo + Matrix::Diagonal(m, x = alpha)
  
  fac <- Matrix::Cholesky(M_alpha, LDL = FALSE, perm = TRUE)
  
  # A(alpha) = solve(M_alpha, Gom)
  A <- Matrix::solve(fac, Gom)
  
  ZidMi_hat <- ZidOi %*% A
  rss <- sum((ZidMi_hat - ZidMi)^2) / n_obs
  return(rss)
  df <- df_sparse_exact(Goo = Goo, alpha = alpha)
  
  penalty <- 1 - df / n_obs
  if (abs(penalty) < 1e-8) {
    penalty <- sign(penalty + 1e-16) * 1e-8
  }
  
  rss / (penalty^2)
}

trace_inverse_hutchinson <- function(fac, m, n_probe = 20L) {
  # Approximate tr(M^{-1}) using Hutchinson's estimator.
  
  acc <- 0
  for (r in seq_len(n_probe)) {
    z <- sample(c(-1, 1), size = m, replace = TRUE)
    x <- Matrix::solve(fac, z)
    acc <- acc + sum(z * x)
  }
  acc / n_probe
}
df_sparse_hutchinson <- function(Goo, alpha, n_probe = 20L) {
  if (!inherits(Goo, "Matrix")) {
    Goo <- Matrix::Matrix(Goo, sparse = TRUE)
  }
  
  m <- nrow(Goo)
  M_alpha <- Goo + Matrix::Diagonal(m, x = alpha)
  fac <- Matrix::Cholesky(M_alpha, LDL = FALSE, perm = TRUE)
  
  tr_inv <- trace_inverse_hutchinson(fac, m = m, n_probe = n_probe)
  as.numeric(m - alpha * tr_inv)
}
alpha_objective_sparse_hutch <- function(alpha,
                                         Goo,
                                         Gom,
                                         ZidOi,
                                         ZidMi,
                                         n_obs,
                                         n_probe = 20L) {
  # Fast approximate GCV objective using sparse Cholesky + Hutchinson trace estimator.
  
  if (!inherits(Goo, "Matrix")) {
    Goo <- Matrix::Matrix(Goo, sparse = TRUE)
  }
  
  m <- nrow(Goo)
  M_alpha <- Goo + Matrix::Diagonal(m, x = alpha)
  fac <- Matrix::Cholesky(M_alpha, LDL = FALSE, perm = TRUE)
  
  A <- Matrix::solve(fac, Gom)
  ZidMi_hat <- ZidOi %*% A
  rss <- sum((ZidMi_hat - ZidMi)^2)
  
  tr_inv <- trace_inverse_hutchinson(fac, m = m, n_probe = n_probe)
  df <- as.numeric(m - alpha * tr_inv)
  
  penalty <- 1 - df / n_obs
  if (abs(penalty) < 1e-8) {
    penalty <- sign(penalty + 1e-16) * 1e-8
  }
  
  rss / (penalty^2)
}
# ------------------------------------------------------------
# Pure R version of reconsX_internal
# ------------------------------------------------------------
fast_alpha_setup <- function(Goo, Gom, ZidOi, ZidMi, n_obs) {
  # Precompute objects used repeatedly in the alpha optimization.
  #
  # Args:
  #   Goo: covariance block on observed coordinates
  #   Gom: cross-covariance block between observed and missing coordinates
  #   ZidOi: observed block of fully observed subjects
  #   ZidMi: missing block of fully observed subjects
  #   n_obs: number of fully observed subjects
  #
  # Returns:
  #   A list of precomputed quantities for fast evaluation.
  
  eig <- eigen(Goo, symmetric = TRUE)
  lambda <- pmax(eig$values, 0)
  Q <- eig$vectors
  
  Qt_Gom <- crossprod(Q, Gom)   # Q^T %*% Gom
  ZidOi_Q <- ZidOi %*% Q        # ZidOi %*% Q
  pos_idx <- lambda > 0
  
  list(
    lambda = lambda,
    Q = Q,
    Qt_Gom = Qt_Gom,
    ZidOi_Q = ZidOi_Q,
    ZidMi = ZidMi,
    n_obs = n_obs,
    pos_idx = pos_idx
  )
}
reconsX_internal_R <- function(O, Z, id_obs, G, i, alpha = 0.0) {
  # Reconstruct the missing part of row i of Z.
  #
  # Args:
  #   O: binary observation matrix (n x d), with 1 = observed and 0 = missing
  #   Z: data matrix (n x d)
  #   id_obs: indices of fully observed subjects (R-style, 1-based)
  #   G: covariance matrix (d x d)
  #   i: row index to reconstruct (R-style, 1-based)
  #
  # Returns:
  #   A 1 x d matrix containing the reconstructed i-th row of Z
  
  if (!is.matrix(O)) {
    stop("'O' must be a matrix.", call. = FALSE)
  }
  if (!is.matrix(Z)) {
    stop("'Z' must be a matrix.", call. = FALSE)
  }
  if (!is.matrix(G) || nrow(G) != ncol(G)) {
    stop("'G' must be a square matrix.", call. = FALSE)
  }
  if (!all(dim(O) == dim(Z))) {
    stop("'O' and 'Z' must have the same dimensions.", call. = FALSE)
  }
  if (ncol(Z) != nrow(G)) {
    stop("ncol('Z') must match nrow('G').", call. = FALSE)
  }
  if (!is.numeric(id_obs) || length(id_obs) == 0L) {
    stop("'id_obs' must be a non-empty numeric vector of row indices.", call. = FALSE)
  }
  if (!is.numeric(i) || length(i) != 1L || i < 1 || i > nrow(Z)) {
    stop("'i' must be a valid row index.", call. = FALSE)
  }
  
  # Observed and missing coordinates for subject i
  Oi <- which(O[i, ] == 1)
  Mi <- which(O[i, ] == 0)
  
  # Nothing to reconstruct
  if (length(Mi) == 0L) {
    return(Z[i, , drop = FALSE])
  }
  
  # No observed part: cannot reconstruct using this rule
  if (length(Oi) == 0L) {
    warning("Row ", i, " has no observed entries; returning original row.", call. = FALSE)
    return(Z[i, , drop = FALSE])
  }
  
  # Training blocks from fully observed subjects
  ZidOi <- Z[id_obs, Oi, drop = FALSE]
  ZidMi <- Z[id_obs, Mi, drop = FALSE]
  
  # Covariance sub-blocks
  Gom <- G[Oi, Mi, drop = FALSE]
  Goo <- G[Oi, Oi, drop = FALSE]
  
  
  if(alpha > 0.0) {
    alpha_opt <- alpha
    
    R <- chol(Goo + diag(alpha_opt, nrow(Goo)))
    A <- backsolve(R, forwardsolve(t(R), Gom))
    
    # gcv_fun_R(alpha_opt, Goo, Gom, ZidOi, ZidMi, length(id_obs))
    # 
    # alpha_vec <- c(100, 50, 25, 10, 1, .1, .01, .001, .0001, .00001, .000001)
    # system.time({
    #   gcv_val <- sapply(alpha_vec, gcv_fun_R, Goo = Goo, Gom = Gom, ZidOi = ZidOi,
    #                     ZidMi = ZidMi, n_obs = length(id_obs))
    # })
    # plot(alpha_vec, log(gcv_val), type = "b")
  } else{
    # setup_alpha <- fast_alpha_setup(
    #   Goo = Goo,
    #   Gom = Gom,
    #   ZidOi = ZidOi,
    #   ZidMi = ZidMi,
    #   n_obs = length(id_obs)
    # )  
    # 
    # opt <- optimize(
    #   f = alpha_objective_fast,
    #   interval = c(2.22e-16, 5),
    #   setup = setup_alpha
    # )
    # 
    # alpha_opt <- opt$minimum
    # 
    # A <- compute_A_fast_R(alpha = alpha_opt, setup = setup_alpha)
    
    # alpha_vec <- exp(seq(log(10), log(1E-6), l = 21))
    # gcv_val <- sapply(alpha_vec, alpha_objective_fast, setup = setup_alpha)
    # alpha_opt <- alpha_vec[which.min(gcv_val)]
    # 
    #alpha_vec <- exp(seq(log(1), log(1E-6), l = 11))
    #gcv_val <- sapply(alpha_vec, gcv_fun_R, Goo = Goo, Gom = Gom, ZidOi = ZidOi,
    #                  ZidMi = ZidMi, n_obs = length(id_obs))
    #plot(alpha_vec, log(gcv_val), type = "b")
    # alpha_opt <- alpha_vec[which.min(gcv_val)]
    # 
    opt <- optimize(
      f = gcv_fun_R,
      interval = c(2.22e-16, 10),
      Goo = Goo, Gom = Gom, ZidOi = ZidOi,
      ZidMi = ZidMi, n_obs = length(id_obs)
    )
    alpha_opt <- opt$minimum
    R <- chol(Goo + diag(alpha_opt, nrow(Goo)))
    A <- backsolve(R, forwardsolve(t(R), Gom))
  }
  
  # # Eigendecomposition of Goo
  # eig <- eigen(Goo, symmetric = TRUE)
  # lambda <- pmax(eig$values, 0)
  # Q <- eig$vectors
  # 
  # # Optimize alpha
  # opt <- optimize(
  #   f = alpha_objective_R,
  #   interval = c(2.22e-16, 100.0),
  #   lambda = lambda,
  #   Q = Q,
  #   Gom = Gom,
  #   ZidOi = ZidOi,
  #   ZidMi = ZidMi,
  #   n_obs = length(id_obs)
  # )
  # 
  # alpha_opt <- opt$minimum
  # if (isTRUE(all.equal(alpha_opt, 100.0))) {
  #   alpha_opt <- 1e-12
  # }
  # 
  # # Reconstruct missing part of row i
  # Ginv <- compute_Ginv_R(lambda = lambda, Q = Q, alpha = alpha_opt)
  # A <- Ginv %*% Gom
  
  Zi <- Z[i, , drop = FALSE]
  Zi[, Mi] <- Zi[, Oi, drop = FALSE] %*% A
  
  temp <- list(Zi = Zi, alpha_opt = alpha_opt)
  return(temp)
}
reconsX_internal_R_sparse <- function(O, Z, id_obs, G, i, alpha = 0.0) {
  # Reconstruct the missing part of row i of Z.
  #
  # Args:
  #   O: binary observation matrix (n x d), with 1 = observed and 0 = missing
  #   Z: data matrix (n x d)
  #   id_obs: indices of fully observed subjects (R-style, 1-based)
  #   G: covariance matrix (d x d)
  #   i: row index to reconstruct (R-style, 1-based)
  #
  # Returns:
  #   A 1 x d matrix containing the reconstructed i-th row of Z
  
  # if (!is.matrix(O)) {
  #   stop("'O' must be a matrix.", call. = FALSE)
  # }
  # if (!is.matrix(Z)) {
  #   stop("'Z' must be a matrix.", call. = FALSE)
  # }
  # if (!is.matrix(G) || nrow(G) != ncol(G)) {
  #   stop("'G' must be a square matrix.", call. = FALSE)
  # }
  if (!all(dim(O) == dim(Z))) {
    stop("'O' and 'Z' must have the same dimensions.", call. = FALSE)
  }
  if (ncol(Z) != nrow(G)) {
    stop("ncol('Z') must match nrow('G').", call. = FALSE)
  }
  if (!is.numeric(id_obs) || length(id_obs) == 0L) {
    stop("'id_obs' must be a non-empty numeric vector of row indices.", call. = FALSE)
  }
  if (!is.numeric(i) || length(i) != 1L || i < 1 || i > nrow(Z)) {
    stop("'i' must be a valid row index.", call. = FALSE)
  }
  
  # Observed and missing coordinates for subject i
  Oi <- which(O[i, ] == 1)
  Mi <- which(O[i, ] == 0)
  
  # Nothing to reconstruct
  if (length(Mi) == 0L) {
    return(Z[i, , drop = FALSE])
  }
  
  # No observed part: cannot reconstruct using this rule
  if (length(Oi) == 0L) {
    warning("Row ", i, " has no observed entries; returning original row.", call. = FALSE)
    return(Z[i, , drop = FALSE])
  }
  
  # Covariance sub-blocks
  Gom <- G[Oi, Mi, drop = FALSE]
  Goo <- G[Oi, Oi, drop = FALSE]
    
  if(alpha > 0) {
    alpha_opt <- alpha
  } else {
    # Training blocks from fully observed subjects
    ZidOi <- Z[id_obs, Oi, drop = FALSE]
    ZidMi <- Z[id_obs, Mi, drop = FALSE]
    
    opt <- optimize(
      f = alpha_objective_sparse,
      interval = c(2.22e-16, 10.0),
      Goo = Goo,
      Gom = Gom,
      ZidOi = ZidOi,
      ZidMi = ZidMi,
      n_obs = nrow(ZidOi)
    )
      
    alpha_opt <- opt$minimum
  }
    
  M_alpha <- Goo + Matrix::Diagonal(n = nrow(Goo), x = alpha_opt)
  fac <- Matrix::Cholesky(M_alpha, LDL = FALSE, perm = TRUE)
  A <- Matrix::solve(fac, Gom)
  
  Zi <- Z[i, , drop = FALSE]
  Zi[, Mi] <- Zi[, Oi, drop = FALSE] %*% A
  
  temp <- list(Zi = Zi, alpha_opt = alpha_opt)
  return(temp)
}

pofggm <- function(id_pobs, id_obs, X, Phi, tp, Sgm.hat = NULL, Tht.hat = NULL,
                   wTht = NULL, rho, gamma = 0.0, alpha = 0.0, pendiag = FALSE,
                   maxit.em = 100L, maxit.admm = 1E5, thr.em = 1E-5, thr.admm = 1e-6, 
                   # sparse = FALSE, 
                   ncores = 1L, verbose = FALSE, ...) {
  
  d <- dim(X)[1L]
  n <- dim(X)[2L]
  p <- dim(X)[3L]
  K <- dim(Phi)[2L]
  conv <- FALSE
  
  fk <- rep(1 / K, K)
  dt <- tp[2L] - tp[1L]
  
  # Costruzione Z e maschera di osservazione O
  Z <- t(apply(X, 2L, identity))
  O <- !is.na(Z)
  
  # if(sparse) {
  #   Z <- Matrix::Matrix(Z, sparse = TRUE)
  #   O <- Matrix::Matrix(O, sparse = TRUE)
  # }
  
  # Inizializzazione Xi.o (score osservati)
  Xi.o <- array(0, dim = c(n, p, K))
  X_na <- X
  X_na[is.na(X_na)] <- 0
  Xi.o <- integrate_cube(X_na, Phi, tp)
  
  # Stima iniziale Sgm se mancante
  if (is.null(Sgm.hat)) {
    Sgm.hat <- array(0, dim = c(p, p, K))
    if(n > p) {
      for (k in seq_len(K)) Sgm.hat[, , k] <- var(Xi.o[, , k]) * (n - 1) / n
    } else {
      for (k in seq_len(K)) diag(Sgm.hat[, , k]) <- 1 / (k * p)
    }
  }
  
  # Stima iniziale Theta se mancante
  if (is.null(Tht.hat)) {
    Tht.hat <- array(0, dim = c(p, p, K))
    for (k in seq_len(K)) {
      Tht.hat[, , k] <- solve(Sgm.hat[, , k])
    }
  }
  
  # Array dei pesi per la penalizzazione di Theta se mancante
  if(is.null(wTht)) {
    U <- 1.0 - diag(p)
    wTht <- list2array(lapply(seq_len(K), \(k) U))
  }
  
  # Iterazioni principali
  for (iter in seq_len(maxit.em)) {

    jj <- generate_jj(p, d) + 1
    G <- kron_sum(Sgm.hat, Phi)
    # if(sparse) G <- Matrix::Matrix(G, sparse = TRUE)
    Xhat <- X
    alpha_opt <- NULL
    count <- 0
    for (i in id_pobs) {
      # print(i)
      count <- count + 1
      # if(!sparse) {
        temp <- reconsX_internal_R(O, Z, id_obs, G, i, alpha = alpha[count])
        X_temp <- temp$Zi
      # } else {
      #   temp <- reconsX_internal_R_sparse(O, Z, id_obs, G, i, alpha = alpha[count]) 
      #   X_temp <- as.matrix(temp$Zi)
      # }
      alpha_opt <- c(alpha_opt, temp$alpha_opt)
      for (j in seq_len(p)) {
        id_j <- jj == j
        Xhat[, i, j] <- X_temp[, id_j]
      }
    }
    
    # system.time(reconsX_internal(O, Z, id_obs - 1, G, id_pobs[1] - 1))
    # system.time(reconsX_internal_R(O = O, Z = Z, id_obs = id_obs, G = G, i = id_pobs[1], alpha = 0.0))
    # system.time(reconsX_internal_fast(O, Z, id_obs - 1, G, id_pobs[1] - 1, alpha = 0.0))
    
    Xihat <- integrate_cube(Xhat, Phi, tp)
    tmp <- compute_S_and_rho_max(Xihat)
    tmp$X <- Xhat
    tmp$Xi <- Xihat

    S <- tmp$S
    tmp2 <- admm_tht_sub(p, K, fk, S, wTht, pendiag, rho, gamma,
                         as.integer(maxit.admm), thr.admm, Tht.hat, 0L)

    Tht.hat <- tmp2$Th
    df <- tmp2$df
    err <- 0
    Sgm <- Sgm.hat
    for(k in seq_len(K)) {
      Sgm.hat[, , k] <- chol2inv(chol(Tht.hat[, , k]))
      err <- err + norm(Sgm.hat[, , k] - Sgm[, , k], type = "F") / (df[k] + p)
    }
    tmp2$Sgm <- Sgm.hat
    tmp$Theta <- tmp2
    tmp$err <- err / K

    out <- tmp
    
    # out <- reconsX(O = O, Z = Z, id_obs = id_obs - 1, id_pobs = id_pobs - 1, Phi = Phi, Sgm = Sgm.hat,
    #                X = X, Xi = Xi.o, tp = tp, fk = fk, wTht = wTht, pendiag = pendiag, rho = rho, gamma = gamma,
    #                alpha = alpha, maxit = maxit.admm, thr = thr.admm, Tht = Tht.hat, n_threads = ncores)
    
    # Verifica convergenza
    err <- out$err
    if (verbose) cat("iter:", iter, "error =", err, "\n")
    if (err < thr.em) {
      conv <- TRUE
      break
    }
    
    # Reinizializzo Sgm e Tht se non è arrivato a convergenza
    Sgm.hat <- out$Theta$Sgm
    Tht.hat <- out$Theta$Tht
  }
  
  # Output list
  list(Sgm.hat = Sgm.hat,
       Tht.hat = Tht.hat,
       S = out$S,
       Xi = out$Xi,
       Phi = Phi,
       Ximputed = out$X, #compute_Ximputed(Phi, out$Xi$Xi),
       wTht = wTht,
       df = out$Theta$df,
       rho.max = out$rho.max,
       rho = rho,
       gamma = gamma,
       alpha_opt = alpha_opt,
       nstep = iter,
       conv = conv)
}

fit_pofggm_path <- function(X,
                            Phi_emp,
                            tp,
                            id_pobs,
                            id_obs,
                            perc_rho,
                            # sparse = FALSE,
                            ncores,
                            verbose = FALSE,
                            alpha = 0.0,
                            gamma = 0.0,
                            maxit_admm = 1e5,
                            thr_em = 1e-5,
                            thr_admm = 1e-6) {
  # Fit the proposed method along a rho path.
  #
  # Returns:
  #   A list with:
  #     ok
  #     fit_list
  #     rho_max
  #     grid_rho
  #     elapsed
  
  n_rho <- length(perc_rho)
  fit_list <- vector("list", length = n_rho)

  alpha_vec <- alpha
  if(length(alpha) == 1) alpha_vec <- rep(alpha, length(id_pobs))
  
  elapsed <- system.time({
    if (verbose) {
      cat("\nRho iteration:", 1L, "of", n_rho, "\n")
    }
    fit_list[[1]] <- pofggm(
      id_pobs = id_pobs,
      id_obs = id_obs,
      X = X,
      Phi = Phi_emp,
      tp = tp,
      Sgm.hat = NULL,
      Tht.hat = NULL,
      wTht = NULL,
      maxit.admm = maxit_admm,
      rho = .Machine$double.xmax,
      gamma = gamma,
      alpha = alpha_vec,
      # sparse = sparse,
      ncores = ncores,
      verbose = verbose,
      thr.em = thr_em,
      thr.admm = thr_admm
    )
    
    rho_max <- fit_list[[1]]$rho.max
    grid_rho <- rho_max * perc_rho
    alpha_vec <- fit_list[[1]]$alpha_opt
    
    if (n_rho >= 2L) {
      for (r in 2:n_rho) {
        if (verbose) {
          cat("\nRho iteration:", r, "of", n_rho, "\n")
        }
        
        fit_list[[r]] <- pofggm(
          id_pobs = id_pobs,
          id_obs = id_obs,
          X = X,
          Phi = Phi_emp,
          tp = tp,
          Sgm.hat = fit_list[[r - 1]]$Sgm.hat,
          Tht.hat = fit_list[[r - 1]]$Tht.hat,
          wTht = NULL,
          maxit.admm = maxit_admm,
          rho = grid_rho[r],
          gamma = gamma,
          alpha = alpha_vec,
          # sparse = sparse,
          ncores = ncores,
          verbose = verbose,
          thr.em = thr_em,
          thr.admm = thr_admm
        )
      }
    }
  })["elapsed"]
  
  list(
    ok = TRUE,
    fit_list = fit_list,
    rho_max = rho_max,
    grid_rho = grid_rho,
    elapsed = unname(elapsed)
  )
}

