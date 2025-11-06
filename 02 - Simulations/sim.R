#--------- Library -----------------------
library(fda)
library(BDgraph)

#---- set path ------------
source (file = "../01 - Code/pofggm.R")
Rcpp::sourceCpp("../01 - Code/pofggm.cpp")

uFPCA <- function(G, argvals, pev, efunctions_sign) {
  # Numerical integration for calculation of 
  # eigenvalues / eigenfunctions (see Ramsay & Silverman, Ch. 8.4.3)
  w  <- quadWeights(argvals, method = "trapezoidal")
  sqrt.w <- sqrt(w)
  V <- G * outer(sqrt.w, sqrt.w)
  out.eigen <- eigen(V, symmetric = TRUE)
  evalues <- replace(out.eigen$values, which(out.eigen$values <= 0), 0)
  npc <- length(evalues[evalues > 0])
  cumpev <- cumsum(evalues[evalues > 0]) / sum(evalues[evalues > 0])
  if (!missing(pev)) {
    npc <- which(cumpev >= pev)[1L]
  } else {
    pev <- cumpev[npc]
  }
  efunctions <- matrix(1 / sqrt.w, nrow = length(sqrt.w), ncol = npc) * out.eigen$vectors[, seq_len(npc), drop = FALSE]
  if (!missing(efunctions_sign)) {
    if (nrow(efunctions_sign) != nrow(efunctions)) {
      stop("wrong number of rows in efunctions_sign")
    }
    if (ncol(efunctions_sign) < npc) {
      stop("wrong number of columns in efunctions_sign")
    }
    efunctions_sign <- efunctions_sign[, seq_len(npc), drop = FALSE]
    sgn <- apply(efunctions * efunctions_sign, 2, \(x) sign(sum(w * x)))
    efunctions <- sweep(efunctions, MARGIN = 2L, STATS = sgn, FUN = "*")
  }
  # # use correct matrix for eigenvalue problem
  evalues <- eigen(V, symmetric = TRUE, only.values = TRUE)$values[1:npc]
  out <- list(evalues = evalues, efunctions = efunctions, npc = npc, pev = pev, 
              cumpev = cumpev, w = w)
  return(out)
}
quadWeights <- function(argvals, method = "trapezoidal"){
  ret <- switch(method, 
                trapezoidal = {
                  D <- length(argvals)
                  0.5 * c(argvals[2] - argvals[1], argvals[3:D] - argvals[1:(D - 2)], argvals[D] - argvals[D - 1])
                }, 
                midpoint = c(0, diff(argvals)), 
                stop("function quadWeights: choose either trapezoidal or midpoint quadrature rule")
  )
  return(ret)
}
meanKraus <- function(X_mat){
  rowMeans(X_mat, na.rm = TRUE)
}
covKraus <- function(X_mat) {
  p <- nrow(X_mat)
  n <- ncol(X_mat)
  
  # Centering
  X_cent_mat <- X_mat - meanKraus(X_mat)
  idMiss <- is.na(X_cent_mat)
  
  # Conta quanti valori non mancanti contribuiscono a ciascuna cella
  valid_counts <- tcrossprod(!idMiss)
  
  # Per ogni coppia (s,t) vogliamo la media degli s * t ignorando NA
  # Uso tcrossprod per calcolare tutte le cross-produttorie insieme
  X_cent_mat[idMiss] <- 0
  prod_mat <- tcrossprod(X_cent_mat)  # p x p con somma dei prodotti (inclusi NA)
  
  # Stima della covarianza
  covKraus_mat <- prod_mat / valid_counts
  
  # Dove non c'è nessun contributo valido, mettiamo NA
  covKraus_mat[valid_counts == 0] <- NA
  
  return(covKraus_mat)
}
krauss_gia <- function(X_mat, alpha = NULL, Sgm = NULL){
  meanKraus <- function(X_mat){
    rowMeans(X_mat, na.rm = TRUE)
  }
  ##
  covKraus <- function(X_mat) {
    p <- nrow(X_mat)
    n <- ncol(X_mat)
    
    # Centering
    X_cent_mat <- X_mat - meanKraus(X_mat)
    idMiss <- is.na(X_cent_mat)
    
    # Conta quanti valori non mancanti contribuiscono a ciascuna cella
    valid_counts <- tcrossprod(!idMiss)
    
    # Per ogni coppia (s,t) vogliamo la media degli s * t ignorando NA
    # Uso tcrossprod per calcolare tutte le cross-produttorie insieme
    X_cent_mat[idMiss] <- 0
    prod_mat <- tcrossprod(X_cent_mat)  # p x p con somma dei prodotti (inclusi NA)
    
    # Stima della covarianza
    covKraus_mat <- prod_mat / valid_counts
    
    # Dove non c'è nessun contributo valido, mettiamo NA
    covKraus_mat[valid_counts == 0] <- NA
    
    return(covKraus_mat)
  }
  ##
  if(is.null(Sgm)) {
    cov_mat       <- covKraus(X_mat)
    mean_vec      <- meanKraus(X_mat)
  } else {
    cov_mat <- Sgm
    mean_vec <- rep(0, nrow(Sgm))
  }
  n             <- ncol(X_mat)
  ##
  reconst_fcts <- which(apply(is.na(X_mat), 2, any))
  NonNA_fcts    <- setdiff(1:n, reconst_fcts)
  ##
  X_reconst_mat <- X_mat[, reconst_fcts, drop = FALSE] - mean_vec
  X_Compl_mat   <- X_mat[, NonNA_fcts, drop = FALSE]
  X_cent <- X_Compl_mat - mean_vec
  ##
  alpha_vec     <- rep(NA, length(reconst_fcts)) 
  df_vec        <- rep(NA, length(reconst_fcts)) 
  ##
  for(i in seq_along(reconst_fcts)){
    X_tmp <- X_reconst_mat[, i]
    ##
    M_bool_vec <- is.na(X_tmp)
    O_bool_vec <- !M_bool_vec
    p_obs <- sum(O_bool_vec)
    ##
    covMO_mat  <- cov_mat[M_bool_vec, O_bool_vec] 
    covOO_mat  <- cov_mat[O_bool_vec, O_bool_vec]
    ##
    tmp <- base::eigen(covOO_mat, symmetric = TRUE)
    lambda <- pmax(tmp$values, 0)
    Q <- tmp$vectors
    ##
    if(is.null(alpha)){
      temp <- stats::optimize(
        f = function(alpha) {
          ##
          df <- sum(lambda[lambda > 0] / (lambda[lambda > 0] + alpha))
          covOO_a_mat_inv <- tcrossprod(Q * rep(sqrt(1 / (lambda + alpha)), each = p_obs))
          B <- covMO_mat %*% covOO_a_mat_inv
          ##
          X_M_fit_cent_vec <- B %*% X_cent[O_bool_vec, ]
          X_cent[M_bool_vec, ] <- X_M_fit_cent_vec
          X_fit <- X_cent + mean_vec
          ##
          rss <- sum((X_fit[M_bool_vec, ] - X_Compl_mat[M_bool_vec, ])^2)
          gcv <- rss / ((1 - df / ncol(X_Compl_mat))^2)
          attr(gcv, "df") <- df
          attr(gcv, "B") <- B
          gcv
        },
        interval = c(1E-12, sum(diag(cov_mat))*n), maximum = FALSE)
      alpha_vec[i] <- temp$minimum
      df_vec[i] <- attr(temp$objective, "df")
      B <- attr(temp$objective, "B")
    } else {
      alpha_vec[i] <- alpha
      ##
      df_vec[i] <- sum(lambda[lambda > 0] / (lambda[lambda > 0] + alpha_vec[i]))
      covOO_a_mat_inv <- tcrossprod(Q * rep(sqrt(1 / (lambda + alpha_vec[i])), each = p_obs))
      B <- covMO_mat %*% covOO_a_mat_inv
    }
    X_tmp_M_cent_vec <- B %*% X_tmp[O_bool_vec]
    X_tmp[M_bool_vec] <- X_tmp_M_cent_vec
    X_reconst_mat[, i] <- X_tmp + mean_vec
  }
  
  X_mat[, reconst_fcts] <- X_reconst_mat
  
  return(list("X_reconst" = X_mat,
              "alpha"     = alpha_vec, 
              "df"        = df_vec,
              "id_obs"    = NonNA_fcts,
              "id_pobs"   = reconst_fcts))
}
compute_metrics_graph <- function(Theta.true, Theta.hat) {
  n_thresh <- length(Theta.hat)
  
  # Considera solo la parte superiore (simmetrica, senza diagonale)
  idx <- upper.tri(Theta.true, diag = FALSE)
  
  true_edges <- Theta.true[idx] != 0
  
  # Soglie decrescenti da max a min
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
    AUC_PR  = auc_approx(recall, precision)
  ))
}
auc_approx <- function(x, y) {
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * zoo::rollmean(y, 2))
}
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
rStars <- function(d, nstars = 1, tht.min = 0.4, tht.max = 0.5) {
  # Arguments
  # d: number of outer vertices 
  # nstars: number of stars
  # tht.min, tht.max:
  require("pracma")
  out <- vector(mode = "list", length = nstars)
  for(l in seq_len(nstars)) {
    Tht <- diag(d)
    Tht.yx <- runif(d - 1, tht.min, tht.max) * sample(c(-1, 1), d - 1, TRUE)
    K <- outer(Tht.yx, Tht.yx)
    Tht.xx <- diag((1 + sqrt(1 + 4 * diag(K))) / 2)
    Tht[1, 1] <- 1 + drop(Tht.yx %*% solve(Tht.xx) %*% Tht.yx)
    Tht[1, -1] <- Tht[-1, 1] <- Tht.yx
    Tht[-1, -1] <- Tht.xx
    out[[l]] <- Tht
  }
  return(do.call(blkdiag, out))
}
rPar <- function(p, K, id.diag, s1 = 3, s2 = -1.8) {
  # Argomenti
  # p: dimesione del vettore degli scores
  # K: numero di termini utilizzati per l'espansione multivariate di KL
  # id.diag = 
  # s1, s2: fattore di scale
  
  #--- Setting Output
  Sgm <- Tht <- array(0, dim = c(p, p, K), 
                      dimnames = list(paste0("X", seq_len(p)),
                                      paste0("X", seq_len(p)),
                                      paste0("K", seq_len(K))))
  for (k in seq_len(K)) {
    Tht[, , k] <- rTht(p = p, id.diag = id.diag)
    sc <- s1 * k^s2 
    Tht[, , k] <- Tht[, , k] / sc
    Sgm[, , k] <- solve(Tht[, , k])
  }
  return(list(Tht = Tht, Sgm = Sgm))
}
simulTheta <- function(n, p, prob, alpha, seed = 1234, vis = TRUE, 
                       graph = c("random", "smallworld", "star", "circle", "cluster", "scale-free", 
                                 "lattice", "hub", "AR(1)", "AR(2)"), ...){
  graph <- match.arg(graph)
  # full.size <- p * (p-1) / 2
  # size <- round(alpha * prob * full.size)
  set.seed(seed)
  data.sim <- bdgraph.sim(p = p, n = n, prob = prob * alpha, vis = vis, graph = graph, ...)
  data.sim
}

configs <- expand.grid(
  perc_window = c(.25, .5, .75),
  perc_obs_curves = c(.25, .5, .75),
  n = c(100), 
  p = c(15),
  method = c("star", "band", "smallworld")
)

for(iconfig in 1:nrow(configs)) {
  # iconfig <- 1
  
  n.sim <- 250                                      # numero di simulazioni
  d <- 50                                          # numero di punti temporali
  K <- 3                                           # numero di basi utilizzate per il troncamento
  p <- configs[iconfig, ][["p"]]                   # number of variables c(5, 10, 20)
  n <- configs[iconfig, ][["n"]]                   # sample size n = 50 / 100 / 200
  w <- d * configs[iconfig, ][["perc_window"]]     # ampiezza della finestra non osservata 20%, 30%, 40% (rispetto a d) c(.2, .3, .4) * d
  p_obs <- configs[iconfig, ][["perc_obs_curves"]] # percentuale di curve non osservate (rispetto a n) c(.25, .50, .75) * n
  method <- configs[iconfig, ][["method"]]         # specific method used for the Theta structure
  ncores <- 4L
  
  perc.rho <- seq(1, 0, by = -.05)
  nrho <- length(perc.rho)
  pofggm.rho <- vector(mode = "list", length = nrho)
  
  # argomenti rPar
  tht.min <- .4
  tht.max <- .5
  s1 <- 3
  s2 <- 1.8
  
  # 1. definizione output
  tp <- seq(from = 0, to = 1, length = d)
  X <- array(0, dim = c(d, n, p), 
             dimnames = list(paste0("tp", seq_len(d)), 
                             seq_len(n), 
                             paste0("Var", seq_len(p))))
  Xi <- array(0, dim = c(n, p, K),
              dimnames = list(seq_len(n), 
                              paste0("Score", seq_len(p)),
                              paste0("k", seq_len(K))))
  
  Xo.save <- Xpo.save <- array(0, dim = c(d, n, p, n.sim),
                               dimnames = list(paste0("tp", seq_len(d)), 
                                               seq_len(n), 
                                               paste0("Var", seq_len(p)),
                                               seq_len(n.sim)))
  Xi.save <- array(0, dim = c(n, p, K, n.sim),
                   dimnames = list(seq_len(n), 
                                   paste0("Score", seq_len(p)),
                                   paste0("k", seq_len(K)),
                                   seq_len(n.sim)))
  
  MSE_Tht_mat2 <- matrix(NA, nrow = nrho, ncol = n.sim)
  MSE_Tht_mat <- matrix(NA, nrow = nrho, ncol = n.sim)
  AUC_Tht_vec2 <- vector(mode = "double", length = n.sim)
  AUC_Tht_vec <- vector(mode = "double", length = n.sim)
  MSE_Xi_mat <- matrix(NA, nrow = nrho, ncol = n.sim)
  MSE_X_mat <- matrix(NA, nrow = nrho, ncol = n.sim)
  MSE_X_Krauss_vec <- vector(mode = "double", length = n.sim)
  
  # 2. creazione basi
  phi.f <- create.fourier.basis(nbasis = K) #Create an Fourier basis object defining a set of Fourier functions with specified period
  Phi <- eval.basis(evalarg = tp, basisobj = phi.f) #creo un vettore formato dalla base valuata nella griglia di punti temporali
  rownames(Phi) <- paste("tp =",  1:d)
  
  # 3. simulazione score sotto ipotesi di parziale separabilità
  if(method == "star") {
    out.rPar <- lapply(seq_len(K), \(k) {
      Tht <- rStars(p / 3, 3, tht.min, tht.max)
      sc <- s1 * k^(-s2)
      Tht <- Tht / sc
      Sgm <- solve(Tht)
      list(Tht = Tht, Sgm = Sgm)
    })
    out.rPar <- list(Tht = list2array(lapply(out.rPar, \(out) out$Tht)),
                     Sgm = list2array(lapply(out.rPar, \(out) out$Sgm)))
  } else {
    if(method == "band") {
      out.rPar <- rPar(p = p, K = K, id.diag = 2, s1 = s1, s2 = -s2) 
    } else {
      out.rPar <- lapply(seq_len(K), \(k) {
        Tht <- round(simulTheta(n, p, prob = .1, alpha = 1, seed = 1234, vis = FALSE, 
                                graph = "smallworld", rewire = .5)$K, 5)
        sc <- s1 * k^(-s2)
        Tht <- Tht / sc
        Sgm <- solve(Tht)
        list(Tht = Tht, Sgm = Sgm)
      })
      out.rPar <- list(Tht = list2array(lapply(out.rPar, \(out) out$Tht)),
                       Sgm = list2array(lapply(out.rPar, \(out) out$Sgm)))
    }
  }
  
  apply(out.rPar$Sgm, 3, \(x) sum(diag(x)))
  image(out.rPar$Tht[, , 1] != 0)
  BDgraph::plot.graph(out.rPar$Tht[, , 1] != 0, main = "Graph structure")
  
  #-------------------------------------------------------------------------------
  # 5. Simulazione
  set.seed(1)
  id_pobs <- sort(sample(seq_len(n), size = n * p_obs, replace = FALSE))
  id_obs <- setdiff(1:n, id_pobs) # insieme C nel paper di Kraus
  
  plot.it <- TRUE
  verbose <- TRUE
  if(!verbose) pb <- txtProgressBar(max = n.sim, style = 3L)
  
  for(isim in seq_len(n.sim)){
    # isim <- 1
    
    set.seed(1234 + isim)
    if(verbose) {
      cat("\nSimulazione", isim, "\n")
    } else {
      setTxtProgressBar(pb = pb, value = isim)
    }
    
    #####################################################
    # Section 1: Simulazione Curve Parzialmente Osservate
    
    Xi <- list2array(apply(out.rPar$Sgm, 3L, \(Sig)
                           MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sig),
                           simplify = FALSE))
    Xi.save[, , , isim] <- Xi
    
    X <- list2array(apply(Xi, 2, \(xi) tcrossprod(Phi, xi), simplify = FALSE))
    Xo.save[, , , isim] <- X
    
    Mask <- array(FALSE, dim = c(d, n, p))
    for(i in id_pobs){
      for (j in seq_len(p)) {
        ii <- sample(d - w + 1, 1)
        Mask[seq(from = ii, length = w), i, j] <- TRUE
      }
    }
    X[Mask] <- NA
    Xpo.save[, , , isim] <- X
    
    if(plot.it) {
      par(mfrow = c(1, 1))
      matplot(Xo.save[, , 1, isim], type = "l", col = 1, lty = 1)
      invisible(sapply(1:n, \(i) lines(which(Mask[, i, 1]),
                                       Xo.save[Mask[, i, 1], i, 1, isim], col = "red", lwd = 2)))
    }
    
    #####################################################
    # Section 2: Stima Basi
    
    Xsmooth <- X
    for(i in 1:n){
      for(j in 1:p){
        curve_tmp <- smooth.spline(x = tp[!is.na(X[,i,j])], y = X[!is.na(X[,i,j]),i,j])$y
        Xsmooth[!is.na(X[,i,j]),i,j] <- curve_tmp
      }
    }
    
    mu.sample <- apply(Xsmooth, 3, meanKraus)
    G.sample <- apply(Xsmooth, 3, covKraus, simplify = FALSE)
    
    argvals <- tp
    row.vec <- rep(argvals, each = d)
    col.vec <- rep(argvals, d)
    G.smooth <- list()
    H.sample <- matrix(0, nrow = d, ncol = d)
    for(i in 1:p){
      G.smooth[[i]] <- matrix(mgcv::gam(as.vector(G.sample[[i]]) ~ te(row.vec, col.vec, k = 10),
                                        weights = as.vector(attr(G.sample[[i]], "cov.count")))$fitted.values, d, d)
      G.smooth[[i]] <- (G.smooth[[i]] + t(G.smooth[[i]])) / 2
      H.sample <- H.sample + G.smooth[[i]]
    }
    H.sample <- H.sample / p
    
    temp_ufpca <- uFPCA(G = H.sample, argvals = tp, pev = .9999)
    temp_ufpca$npc
    temp_ufpca$cumpev
    Phi.emp <- temp_ufpca$efunctions
    
    #####################################################
    # Section 3: Stima Modello Proposto
    mod1 <- pofggm(id_pobs = id_pobs, id_obs = id_obs, X = X, Phi = Phi.emp, tp = tp,
                       Sgm.hat = NULL, Tht.hat = NULL, wTht = NULL, maxit.admm = 1E5, 
                       rho = 100.0, ncores = ncores, verbose = verbose)
    pofggm.rho[[1]] <- mod1
    
    #############################################################################
    # Section 3 bis: Stima Modello Proposto Path Rho
    
    rho.max <- pofggm.rho[[1]]$rho.max
    grid.rho <- rho.max * perc.rho #seq(rho.max, 1E-6, l = nrho)
    
    # fitting path for rho
    for (r in 2:nrho) {
      if(verbose) cat("\nRho iter: ", r, "\n")
      pofggm.rho[[r]] <- pofggm(id_pobs = id_pobs, id_obs = id_obs, X = X, Phi = Phi.emp, tp = tp,
                                    Sgm.hat = pofggm.rho[[r - 1]]$Sgm.hat, 
                                    Tht.hat = pofggm.rho[[r - 1]]$Tht.hat, wTht = NULL, maxit.admm = 1E5, 
                                    rho = grid.rho[r], ncores = ncores, verbose = verbose)
    }
    
    ### KRAUS
    
    result <- lapply(seq_len(p), \(j) krauss_gia(X_mat = X[, , j]))
    XreconsKrauss <- list2array(lapply(seq_len(p), \(j) result[[j]]$X_reconst))

    XiKrauss <- try(integrate_cube(XreconsKrauss, Phi.emp, tp), 
                    silent = TRUE)
    if(inherits(XiKrauss, "try-error")) break
    tmp <- try(compute_S_and_rho_max(XiKrauss), 
               silent = TRUE)
    if(inherits(tmp, "try-error")) break
    rho.max.Krauss <- tmp$rho.max
    grid.rho.Krauss <- rho.max.Krauss * perc.rho
    S_krauss <- tmp$S
    Tht_krauss <- array(0, dim = c(p, p, dim(XiKrauss)[3]))
    for (k in seq_len(dim(XiKrauss)[3])) {
      Tht_krauss[, , k] <- solve(S_krauss[, , k])
    }
    
    THT.K <- vector(mode = "list", length = nrho)
    for(r in seq_len(nrho)){
      if(verbose) cat("\nRho iter: ", r, "\n")
      tmp2 <- try(admm_tht_sub(p = p, N = dim(XiKrauss)[3], fk = rep(1 / dim(XiKrauss)[3], dim(XiKrauss)[3]), 
                               S = S_krauss, wTht = pofggm.rho[[1]]$wTht, 
                               pendiag = FALSE, rho = grid.rho.Krauss[r], alpha = 0.0,
                               maxit = 1E5, thr = 1e-6, Tht = Tht_krauss, trace = 0), 
                  silent = TRUE)
      if(inherits(tmp2, "try-error")) break
      THT.K[[r]] <- Tht_krauss <- tmp2$Tht
    }
    
    #####
    
    if(plot.it) {
      par(mfrow = c(1, 1))
      id <- id_pobs[1:2]
      matplot(tp, Xo.save[, id, 1, isim], type = "l", lty = 1, col = "pink2", lwd = 1.5)
      matpoints(tp, X[, id, 1], pch = 16, col = "black")
      
      matlines(tp, pofggm.rho[[1]]$Ximputed[, id, 1], col = "red2", lty = 2)
      matlines(tp, pofggm.rho[[6]]$Ximputed[, id, 1], col = "blue2", lty = 2)
      matlines(tp, pofggm.rho[[11]]$Ximputed[, id, 1], col = "green2", lty = 2)
      matlines(tp, pofggm.rho[[16]]$Ximputed[, id, 1], col = "purple", lty = 2)
      matlines(tp, pofggm.rho[[21]]$Ximputed[, id, 1], col = "yellow2", lty = 2)
      
      matlines(tp, result[[1]]$X_reconst[, id], col = "cyan2", lty = 2)
    }
    
    #############################################################################

    MSE_Tht <- sapply(pofggm.rho, 
                      \(out) {
                        mean(sapply(seq_len(K), 
                                    \(k) norm(out$Tht[, , k] - out.rPar$Tht[,,k], type = "F") / sum(out.rPar$Tht[,,k] != 0)))
                      })
    (MSE_Tht_mat[, isim] <- MSE_Tht)
    
    MSE_Tht2 <- sapply(THT.K, 
                       \(out) {
                         mean(sapply(seq_len(K), 
                                     \(k) norm(out[, , k] - out.rPar$Tht[,,k], type = "F") / sum(out.rPar$Tht[,,k] != 0)))
                       })
    (MSE_Tht_mat2[, isim] <- MSE_Tht2)
    
    if(plot.it) {
      par(mfrow = c(1, 2))
      plot(perc.rho, MSE_Tht, type = "l", main = "MSE Tht")
      lines(perc.rho, MSE_Tht2, col = 2, lty = 2)
    }
    
    res <- compute_metrics_graph(Theta.true = Reduce("+", array2list(abs(out.rPar$Tht))),
                                 Theta.hat = lapply(pofggm.rho, \(out) Reduce("+", array2list(abs(out$Tht)))))
    
    AUC_Tht <- res$AUC_ROC
    (AUC_Tht_vec[isim] <- AUC_Tht)
    
    res2 <- compute_metrics_graph(Theta.true = Reduce("+", array2list(abs(out.rPar$Tht))),
                                  Theta.hat = lapply(THT.K, \(out) Reduce("+", array2list(abs(out)))))
    
    AUC_Tht2 <- res2$AUC_ROC
    (AUC_Tht_vec2[isim] <- AUC_Tht2)
    
    MSE_X <- rowMeans(
      sapply(seq_len(p), \(j) {
        reconst_fcts <- which(apply(is.na(X[, , j]), 2, any))
        id <- Mask[, reconst_fcts, j]
        colMeans(
          sapply(seq_len(nrho), \(r) {
            sapply(1:length(reconst_fcts), \(k) {
              temp <- (pofggm.rho[[r]]$Ximputed[id[, k], reconst_fcts[k], j] - Xo.save[id[, k], reconst_fcts[k], j, isim])^2
              .5 * diff(tp[id[, k]]) %*% (temp[-1] + temp[-length(temp)])
            })
          })
        )
      })
    )
    (MSE_X_mat[, isim] <- MSE_X)
    
    MSE_X_Krauss <- mean(
      sapply(seq_len(p), \(j) {
        reconst_fcts <- which(apply(is.na(X[, , j]), 2, any))
        id <- Mask[, reconst_fcts, j]
        mean(
          sapply(1:length(reconst_fcts), \(k) {
            temp <- (XreconsKrauss[id[, k], reconst_fcts[k], j] - Xo.save[id[, k], reconst_fcts[k], j, isim])^2
            .5 * diff(tp[id[, k]]) %*% (temp[-1] + temp[-length(temp)])
          })
        )
      })
    )
    (MSE_X_Krauss_vec[isim] <- MSE_X_Krauss)
    
    if(isim %% 25 == 0) {
      temp <- rbind("MSE_THT" = rowMeans(MSE_Tht_mat[, seq_len(isim)]),
                    "MSE_THT_Krauss" = rowMeans(MSE_Tht_mat2[, seq_len(isim)]),
                    "AUC_THT" = mean(AUC_Tht_vec[seq_len(isim)]),
                    "AUC_THT_Krauss" = mean(AUC_Tht_vec2[seq_len(isim)]),
                    "MSE_X" = rowMeans(MSE_X_mat[, seq_len(isim)]),
                    "MSE_X_Kraussvirtualgravity" = mean(MSE_X_Krauss_vec[seq_len(isim)]))
      print(temp)
    }
  }
}
