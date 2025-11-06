library(ggplot2)
library(reshape2)
library(BDgraph)

array2list <- function(x) {
  if(!is.array(x)) stop("x is not a list!")
  dimX <- dim(x)
  K <- tail(dimX, 1)
  y <- vector(mode = "list", length = K)
  for(k in seq_len(K)) y[[k]] <- if(length(dimX) == 2) x[, k] else na.omit(matrix(x[,,k], dimX[1], dimX[2]))
  names(y) <- tail(dimnames(x), 1)[[1]]
  y
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
rTht <- function(p, id.diag, det.min = 0.5) {
  # argomenti
  # p = dimensione della matrice
  # id.diag = 
  # det.min = valore del determinate della matrice di precisione generata
  if (!(id.diag %in% 0:(p - 1)))
    stop("prova")
  rho <- optimise(\(rho) {
    outer(1:p, 1:p, function(i, j) {
      d <- abs(i - j)
      ifelse(d <= id.diag, rho^abs(i - j), 0)
    }) -> Tht
    (det(Tht) - det.min)^2
  }, c(0, 1))$min
  
  rho <- 0.5
  
  outer(1:p, 1:p, function(i, j) {
    d <- abs(i - j)
    ifelse(d <= id.diag, rho^abs(i - j), 0)
  }) -> Tht
  return(Tht)
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

iconfig_vec <- c(1, 10, 19)
out.rPar <- vector(mode = "list", length = length(iconfig_vec))

for(i in seq_along(iconfig_vec)) {
  iconfig <- iconfig_vec[i]
  d <- 50                                          # numero di punti temporali
  K <- 3                                           # numero di basi utilizzate per il troncamento
  p <- configs[iconfig, ][["p"]]                   # number of variables c(5, 10, 20)
  n <- configs[iconfig, ][["n"]]                   # sample size n = 50 / 100 / 200
  w <- d * configs[iconfig, ][["perc_window"]]     # ampiezza della finestra non osservata 20%, 30%, 40% (rispetto a d) c(.2, .3, .4) * d
  p_obs <- configs[iconfig, ][["perc_obs_curves"]] # percentuale di curve non osservate (rispetto a n) c(.25, .50, .75) * n
  method <- configs[iconfig, ][["method"]]         # specific method used for the Theta structure
  
  # argomenti rPar
  tht.min <- .4
  tht.max <- .5
  s1 <- 1
  s2 <- .2
  
  # 3. simulazione score sotto ipotesi di parziale separabilità
  if(method == "star") {
    out.rPar[[i]] <- lapply(seq_len(K), \(k) {
      Tht <- rStars(p / 3, 3, tht.min, tht.max)
      sc <- s1 * k^(-s2)
      Tht <- Tht / sc
      Sgm <- solve(Tht)
      list(Tht = Tht, Sgm = Sgm)
    })
    out.rPar[[i]] <- list(Tht = list2array(lapply(out.rPar[[i]], \(out) out$Tht)),
                          Sgm = list2array(lapply(out.rPar[[i]], \(out) out$Sgm)))
  } else {
    if(method == "band") {
      out.rPar[[i]] <- rPar(p = p, K = K, id.diag = 2, s1 = s1, s2 = -s2) 
    } else {
      out.rPar[[i]] <- lapply(seq_len(K), \(k) {
        Tht <- round(simulTheta(n, p, prob = .1, alpha = 1, seed = 1234, vis = FALSE, 
                                graph = "smallworld", rewire = .5)$K, 5)
        sc <- s1 * k^(-s2)
        Tht <- Tht / sc
        Sgm <- solve(Tht)
        list(Tht = Tht, Sgm = Sgm)
      })
      out.rPar[[i]] <- list(Tht = list2array(lapply(out.rPar[[i]], \(out) out$Tht)),
                            Sgm = list2array(lapply(out.rPar[[i]], \(out) out$Sgm)))
    }
  }
}

Tht <- lapply(out.rPar, \(x) Reduce("+", array2list(abs(x$Tht))))
Tht <- 1*(list2array(Tht) != 0)

# Portiamo in formato long
df <- melt(Tht)
colnames(df) <- c("from", "to", "time", "value")
df$time <- factor(df$time, labels = c("Stars", "Bands", "Smallworld"))

# Creiamo il grafico
p <- ggplot(df, aes(x = from, y = to, fill = factor(value))) +
  geom_tile(color = "grey80") +
  scale_fill_manual(values = c("0" = "white", "1" = "steelblue")) +
  coord_fixed() +
  facet_wrap(~ time, ncol = 3) +
  labs(title = "Adjacency structures",
       x = "",
       y = "",
       fill = "") +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = "bottom")
p
