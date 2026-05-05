# ============================================================
# Precision-matrix structures for simulation study
# ============================================================
# This script generates three representative precision-matrix
# structures used in the simulation study:
#   1. star
#   2. band
#   3. smallworld
# and produces a compact visualization of their adjacency
# patterns aggregated over the K latent layers.
# ============================================================

# ============================================================
# 0. Setup
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(BDgraph)
  library(pracma)
})

set.seed(123)

# ============================================================
# 1. Utility functions
# ============================================================

array2list <- function(x) {
  if (!is.array(x)) stop("`x` must be an array.")
  
  dim_x <- dim(x)
  K <- tail(dim_x, 1)
  out <- vector("list", length = K)
  
  for (k in seq_len(K)) {
    out[[k]] <- if (length(dim_x) == 2) {
      x[, k]
    } else {
      na.omit(matrix(x[, , k], dim_x[1], dim_x[2]))
    }
  }
  
  names(out) <- tail(dimnames(x), 1)[[1]]
  out
}

list2array <- function(x) {
  if (!is.list(x)) stop("`x` must be a list.")
  
  K <- length(x)
  
  if (is.matrix(x[[1]])) {
    n <- sapply(x, nrow)
    p <- ncol(x[[1]])
    
    out <- array(
      0,
      dim = c(max(n), p, K),
      dimnames = list(NULL, colnames(x[[1]]), names(x))
    )
    
    for (k in seq_len(K)) {
      out[seq_len(n[k]), , k] <- x[[k]]
    }
  } else {
    n <- length(x[[1]])
    out <- matrix(
      unlist(x),
      nrow = n,
      ncol = K,
      dimnames = list(names(x[[1]]), names(x))
    )
  }
  
  out
}

# ============================================================
# 2. Precision-matrix generators
# ============================================================

rTht_band <- function(p, id_diag, det_min = 0.5, fix_rho = 0.5) {
  if (!(id_diag %in% 0:(p - 1))) {
    stop("`id_diag` must be between 0 and p - 1.")
  }
  
  # Optional optimization kept for future use
  # rho <- optimise(function(rho) {
  #   Tht <- outer(1:p, 1:p, function(i, j) {
  #     d <- abs(i - j)
  #     ifelse(d <= id_diag, rho^d, 0)
  #   })
  #   (det(Tht) - det_min)^2
  # }, c(0, 1))$minimum
  
  rho <- fix_rho
  
  outer(1:p, 1:p, function(i, j) {
    d <- abs(i - j)
    ifelse(d <= id_diag, rho^d, 0)
  })
}

rPar_band <- function(p, K, id_diag, s1 = 3, s2 = -1.8) {
  Sgm <- Tht <- array(
    0,
    dim = c(p, p, K),
    dimnames = list(
      paste0("X", seq_len(p)),
      paste0("X", seq_len(p)),
      paste0("K", seq_len(K))
    )
  )
  
  for (k in seq_len(K)) {
    Tht[, , k] <- rTht_band(p = p, id_diag = id_diag)
    scale_k <- s1 * k^s2
    Tht[, , k] <- Tht[, , k] / scale_k
    Sgm[, , k] <- solve(Tht[, , k])
  }
  
  list(Tht = Tht, Sgm = Sgm)
}

rStars <- function(d, nstars = 1, tht_min = 0.4, tht_max = 0.5) {
  out <- vector("list", length = nstars)
  
  for (l in seq_len(nstars)) {
    Tht <- diag(d)
    Tht_yx <- runif(d - 1, tht_min, tht_max) * sample(c(-1, 1), d - 1, replace = TRUE)
    
    Kmat <- outer(Tht_yx, Tht_yx)
    Tht_xx <- diag((1 + sqrt(1 + 4 * diag(Kmat))) / 2)
    
    Tht[1, 1] <- 1 + drop(Tht_yx %*% solve(Tht_xx) %*% Tht_yx)
    Tht[1, -1] <- Tht[-1, 1] <- Tht_yx
    Tht[-1, -1] <- Tht_xx
    
    out[[l]] <- Tht
  }
  
  do.call(blkdiag, out)
}

simulTheta <- function(n, p, prob, alpha, seed = 1234, vis = FALSE,
                       graph = c("random", "smallworld", "star", "circle",
                                 "cluster", "scale-free", "lattice",
                                 "hub", "AR(1)", "AR(2)"), ...) {
  graph <- match.arg(graph)
  
  set.seed(seed)
  bdgraph.sim(
    p = p,
    n = n,
    prob = prob * alpha,
    vis = vis,
    graph = graph,
    ...
  )
}

generate_precision_structure <- function(method, p, K, n,
                                         tht_min = 0.4,
                                         tht_max = 0.5,
                                         s1 = 1,
                                         s2 = 0.2,
                                         smallworld_prob = 0.1,
                                         smallworld_rewire = 0.5,
                                         seed = 1234) {
  if (method == "star") {
    out <- lapply(seq_len(K), function(k) {
      Tht <- rStars(d = p / 4, nstars = 4, tht_min = tht_min, tht_max = tht_max)
      scale_k <- s1 * k^(-s2)
      Tht <- Tht / scale_k
      Sgm <- solve(Tht)
      list(Tht = Tht, Sgm = Sgm)
    })
    
    return(list(
      Tht = list2array(lapply(out, `[[`, "Tht")),
      Sgm = list2array(lapply(out, `[[`, "Sgm"))
    ))
  }
  
  if (method == "band") {
    return(rPar_band(
      p = p,
      K = K,
      id_diag = 2,
      s1 = s1,
      s2 = -s2
    ))
  }
  
  if (method == "smallworld") {
    out <- lapply(seq_len(K), function(k) {
      Tht <- round(
        simulTheta(
          n = n,
          p = p,
          prob = smallworld_prob,
          alpha = 1,
          seed = seed,
          vis = FALSE,
          graph = "smallworld",
          rewire = smallworld_rewire
        )$K,
        5
      )
      scale_k <- s1 * k^(-s2)
      Tht <- Tht / scale_k
      Sgm <- solve(Tht)
      list(Tht = Tht, Sgm = Sgm)
    })
    
    return(list(
      Tht = list2array(lapply(out, `[[`, "Tht")),
      Sgm = list2array(lapply(out, `[[`, "Sgm"))
    ))
  }
  
  stop("Unknown method.")
}

# ============================================================
# 3. Simulation design
# ============================================================

configs <- expand.grid(
  perc_window = 0.50,
  perc_obs_curves = 0.50,
  n = 100L,
  p = 20L,
  method = c("star", "band", "smallworld"),
  stringsAsFactors = FALSE
)

# Indices chosen to display one representative configuration
# for each graph type
iconfig_vec <- 1:3

# Common parameters
d <- 50
K <- 5
tht_min <- 0.4
tht_max <- 0.5
s1 <- 3
s2 <- 1.8

# ============================================================
# 4. Generate structures
# ============================================================

out_rPar <- vector("list", length(iconfig_vec))

for (i in seq_along(iconfig_vec)) {
  iconfig <- iconfig_vec[i]
  cfg <- configs[iconfig, ]
  
  out_rPar[[i]] <- generate_precision_structure(
    method = cfg$method,
    p = cfg$p,
    K = K,
    n = cfg$n,
    tht_min = tht_min,
    tht_max = tht_max,
    s1 = s1,
    s2 = s2,
    smallworld_prob = 0.1,
    smallworld_rewire = 0.5,
    seed = 1234
  )
}

names(out_rPar) <- c("Stars", "Bands", "Smallworld")

# ============================================================
# 5. Aggregate adjacency patterns across layers
# ============================================================

adjacency_list <- lapply(out_rPar, function(x) {
  Reduce("+", array2list(abs(x$Tht)))
})

adjacency_array <- 1 * (list2array(adjacency_list) != 0)

df_plot <- melt(adjacency_array)
colnames(df_plot) <- c("from", "to", "structure", "value")
df_plot$structure <- factor(df_plot$structure, levels = c("Stars", "Bands", "Smallworld"))

# ============================================================
# 6. Plot
# ============================================================

p <- ggplot(df_plot, aes(x = from, y = to, fill = factor(value))) +
  geom_tile(color = "grey80") +
  scale_fill_manual(values = c("0" = "white", "1" = "gray40")) +
  coord_fixed() +
  facet_wrap(~structure, ncol = 3) +
  labs(
    title = NULL,
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

print(p)

ggsave("output/figures/adj_struct.pdf", plot = p, 
       device = "pdf", scale = .7, width = 12, height = 6, units = "in", dpi = 300)
