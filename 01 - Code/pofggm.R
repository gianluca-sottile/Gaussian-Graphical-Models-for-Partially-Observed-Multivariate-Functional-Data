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

pofggm <- function(id_pobs, id_obs, X, Phi, tp, Sgm.hat = NULL, Tht.hat = NULL,
                   wTht = NULL, rho, gamma = 0.0, pendiag = FALSE,
                   maxit.em = 100L, maxit.admm = 1E5, thr.em = 1E-5, thr.admm = 1e-6, 
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
  
  # Inizializzazione Xi.o (score osservati)
  Xi.o <- array(0, dim = c(n, p, K))
  X_na <- X
  X_na[is.na(X_na)] <- 0
  Xi.o <- integrate_cube(X_na, Phi, tp)
  
  # Stima iniziale Sgm se mancante
  if (is.null(Sgm.hat)) {
    Sgm.hat <- array(0, dim = c(p, p, K))
    for (k in seq_len(K)) Sgm.hat[, , k] <- var(Xi.o[, , k]) * (n - 1) / n
    # for (k in seq_len(K)) diag(Sgm.hat[, , k]) <- 1 / (k * p)
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

    # jj <- generate_jj(p, d) + 1
    # G <- kron_sum(Sgm.hat, Phi)
    # Xhat <- X
    # for (i in id_pobs - 1) {
    #   print(i + 1)
    #   X_temp <- reconsX_internal(O, Z, id_obs - 1, G, i)
    #   for (j in seq_len(p)) {
    #     id_j <- jj == j
    #     Xhat[, i + 1, j] <- X_temp[, id_j]
    #   }
    # }
    # 
    # Xihat <- integrate_cube(Xhat, Phi, tp)
    # tmp <- compute_S_and_rho_max(Xihat)
    # tmp$X <- Xhat
    # tmp$Xi <- Xihat
    # 
    # S <- tmp$S
    # tmp2 <- admm_tht_sub(p, K, fk, S, wTht, pendiag, rho, gamma,
    #                      as.integer(maxit.admm), thr.admm, Tht.hat, 9L)
    # 
    # Tht.hat <- tmp2$Th
    # df <- tmp2$df
    # err <- 0
    # Sgm <- Sgm.hat
    # for(k in seq_len(K)) {
    #   Sgm.hat[, , k] <- chol2inv(chol(Tht.hat[, , k]))
    #   err <- err + norm(Sgm.hat[, , k] - Sgm[, , k], type = "F") / (df[k] + p)
    # }
    # tmp2$Sgm <- Sgm.hat
    # tmp$Theta <- tmp2
    # tmp$err <- err / K
    # 
    # out <- tmp
    
    out <- reconsX(O = O, Z = Z, id_obs = id_obs - 1, id_pobs = id_pobs - 1, Phi = Phi, Sgm = Sgm.hat, 
                   X = X, Xi = Xi.o, tp = tp, fk = fk, wTht = wTht, pendiag = pendiag, rho = rho, alpha = gamma, 
                   maxit = maxit.admm, thr = thr.admm, Tht = Tht.hat, n_threads = ncores)
    
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
       nstep = iter,
       conv = conv)
}