# source("01 - Code/helper.R")

library(ggplot2)

# ============================================================
# Simulation design and configuration setup
# ============================================================

safe_logdet <- function(Theta, jitter = 1e-10, max_jitter = 1e-3) {
  stopifnot(is.matrix(Theta), nrow(Theta) == ncol(Theta))
  Theta <- 0.5 * (Theta + t(Theta))
  
  try_det <- try(determinant(Theta, logarithm = TRUE), silent = TRUE)
  if (!inherits(try_det, "try-error")) {
    return(as.numeric(try_det$modulus))
  }
  
  jj <- jitter
  repeat {
    Theta_j <- Theta
    diag(Theta_j) <- diag(Theta_j) + jj
    try_det <- try(determinant(Theta_j, logarithm = TRUE), silent = TRUE)
    if (!inherits(try_det, "try-error") && is.finite(as.numeric(try_det$modulus))) {
      return(as.numeric(try_det$modulus))
    }
    jj <- jj * 10
    if (jj > max_jitter) break
  }
  
  NA_real_
}

trace_prod_sym <- function(S, Theta) {
  sum(S * Theta)
}

df_count <- function(Theta, rule = c("edges", "offdiag", "entries")) {
  rule <- match.arg(rule)
  p <- nrow(Theta)
  nz_off_upper <- sum(Theta[upper.tri(Theta)] != 0)
  if (rule == "edges") return(p + nz_off_upper)
  if (rule == "offdiag") return(nz_off_upper)
  if (rule == "entries") return(sum(Theta != 0))
}

compute_ic <- function(n, S, Theta, type = c("AIC", "BIC", "EBIC"),
                       gamma_ebic = 0.5,
                       df_rule = c("edges", "offdiag", "entries"),
                       return_per_slice = FALSE) {
  type <- match.arg(type)
  df_rule <- match.arg(df_rule)
  
  stopifnot(
    length(dim(S)) == 3L,
    length(dim(Theta)) == 3L,
    all(dim(S)[1:2] == dim(Theta)[1:2]),
    dim(S)[3] == dim(Theta)[3]
  )
  
  p <- dim(Theta)[1]
  K <- dim(Theta)[3]
  vals <- numeric(K)
  
  for (k in seq_len(K)) {
    Sk <- S[, , k]
    Tk <- Theta[, , k]
    trST <- trace_prod_sym(Sk, Tk)
    logdetT <- safe_logdet(Tk)
    
    if (!is.finite(logdetT)) {
      vals[k] <- NA_real_
      next
    }
    
    neg2loglik <- n * (trST - logdetT)
    df_k <- df_count(Tk, rule = df_rule)
    
    vals[k] <- switch(
      type,
      AIC = neg2loglik + 2 * df_k,
      BIC = neg2loglik + df_k * log(n),
      EBIC = {
        E_k <- df_count(Tk, rule = "offdiag")
        neg2loglik + df_k * log(n) + 4 * gamma_ebic * E_k * log(p)
      }
    )
  }
  
  if (isTRUE(return_per_slice)) return(vals)
  sum(vals, na.rm = TRUE)
}

compute_aic <- function(n, S, Theta, df_rule = c("edges", "offdiag", "entries"),
                        return_per_slice = FALSE) {
  compute_ic(n, S, Theta, type = "AIC", df_rule = df_rule,
             return_per_slice = return_per_slice)
}

compute_bic <- function(n, S, Theta, df_rule = c("edges", "offdiag", "entries"),
                        return_per_slice = FALSE) {
  compute_ic(n, S, Theta, type = "BIC", df_rule = df_rule,
             return_per_slice = return_per_slice)
}

compute_ebic <- function(gamma_ebic, n, S, Theta,
                         df_rule = c("edges", "offdiag", "entries"),
                         return_per_slice = FALSE) {
  compute_ic(n, S, Theta, type = "EBIC", gamma_ebic = gamma_ebic,
             df_rule = df_rule, return_per_slice = return_per_slice)
}


# ---------------------------
# 1. Global simulation settings
# ---------------------------

n_sim <- 100L                # number of Monte Carlo replications
ncores <- 4L                 # number of cores (if parallelization is used)

plot_it <- FALSE
verbose <- FALSE

# Regularization path for rho
perc_rho <- seq(1, 0, by = -0.05)
n_rho <- length(perc_rho)

# Parameters controlling the true precision matrices
tht_min <- 0.4
tht_max <- 0.5
s1 <- 3
s2 <- 1.8

# ---------------------------
# 2. Simulation scenarios
# ---------------------------

configs <- expand.grid(
  scenario = c("gamma"),
  n = c(100L),
  p_over_n = c(0.2),
  K_true = c(5L),
  perc_window = c(0.5),
  perc_obs_curves = c(0.5),
  pev = c(0.99),
  perc_theta_share = c(1.0),
  d = c(50L),
  graph_type = c("star", "band", "smallworld"),
  stringsAsFactors = FALSE
)

for (iconfig in seq_len(nrow(configs))) {
  
  config <- parse_sim_config(configs[iconfig, ])
  
  n <- config$n
  p <- config$p
  K_true <- config$K_true
  w <- config$window_size
  p_obs <- config$perc_obs_curves
  graph_type <- config$graph_type
  pev <- config$pev
  perc_theta_share <- config$perc_theta_share
  
  d <- config$d                     # number of grid points
  
  # Common time grid
  tp <- seq(from = 0, to = 1, length.out = d)
  
  storage <- initialize_simulation_storage(
    d = d,
    n = n,
    p = p,
    K_true = K_true,
    n_sim = n_sim,
    n_rho = n_rho
  )
  storage$theta_err_mat_2 <- storage$theta_err_mat
  storage$auc_theta_vec_2 <- storage$auc_theta_vec
  storage$curve_err_mat_2 <- storage$curve_err_mat
  storage$comp_time_vec_2 <- storage$comp_time_vec
  
  basis_info <- create_basis_objects(tp = tp, K_true = K_true)
  Phi <- basis_info$Phi
  
  out_rPar <- generate_true_model(
    p = p,
    K_true = K_true,
    graph_type = graph_type,
    perc_theta_share = perc_theta_share,
    tht_min = tht_min,
    tht_max = tht_max,
    s1 = s1,
    s2 = s2, 
    seed_base = 1235 + seq_len(K_true)
  )
  
  if (!verbose) {
    pb <- txtProgressBar(max = n_sim, style = 3L)
  }
  
  for (isim in seq_len(n_sim)) {
    
    set.seed(1234 + isim)
    
    if (verbose) {
      cat("\nSimulation", isim, "\n")
    } else {
      setTxtProgressBar(pb = pb, value = isim)
    }
    
    # ----------------------------------------------------------
    # 1. Define observed and partially observed subjects
    # ----------------------------------------------------------
    
    id_pobs <- sample_partial_subjects(
      n = n,
      perc_obs_curves = p_obs
    )
    id_obs <- setdiff(seq_len(n), id_pobs)
    
    # ----------------------------------------------------------
    # 2. Simulate latent scores and full curves
    # ----------------------------------------------------------
    
    sim_data <- simulate_scores_and_curves(
      Sgm_array = out_rPar$Sgm,
      Phi = Phi,
      n = n,
      p = p,
      K_true = K_true
    )
    
    Xi <- sim_data$Xi
    X_full <- sim_data$X
    
    storage$Xi_save[, , , isim] <- Xi
    storage$Xo_save[, , , isim] <- X_full
    
    # ----------------------------------------------------------
    # 3. Introduce block missingness
    # ----------------------------------------------------------
    
    Mask <- generate_block_missing_mask(
      d = d,
      n = n,
      p = p,
      id_pobs = id_pobs,
      window_size = w
    )
    
    X_obs <- X_full
    X_obs[Mask] <- NA
    storage$Xpo_save[, , , isim] <- X_obs
    
    if (plot_it) {
      par(mfrow = c(1, 1))
      matplot(X_full[, , 1], type = "l", col = 1, lty = 1)
      invisible(
        lapply(seq_len(n), function(i) {
          lines(which(Mask[, i, 1]), X_full[Mask[, i, 1], i, 1], col = "red", lwd = 2)
        })
      )
    }
    
    # ----------------------------------------------------------
    # 4. Estimate empirical basis from partially observed data
    # ----------------------------------------------------------
    
    basis_fit <- try(
      estimate_empirical_basis(
        X = X_obs,
        tp = tp,
        pev = pev,
        smooth_k = 12
      ),
      silent = TRUE
    )
    
    if (inherits(basis_fit, "try-error")) {
      warning("Basis estimation failed at simulation ", isim)
      next
    }
    
    Phi_emp <- basis_fit$Phi_emp
    
    if (verbose) {
      cat("Estimated number of PCs:", basis_fit$ufpca$npc, "\n")
      print(basis_fit$ufpca$cumpev)
    }
    
    # ----------------------------------------------------------
    # 5. Fit proposed method along the rho path
    # ----------------------------------------------------------
    
    fit_pofggm <- vector(mode = "list", length = 5L)
    fit_pofggm_mle <- vector(mode = "list", length = 5L)
    df_ic <- vector(mode = "list", length = 5L)
    count <- 0
    for(gamma in c(0, .25, .5, .75, 1)) {
      count <- count + 1L
      
      fit_pofggm[[count]] <- try(
        fit_pofggm_path(
          X = X_obs,
          Phi_emp = Phi_emp,
          tp = tp,
          id_pobs = id_pobs,
          id_obs = id_obs,
          perc_rho = perc_rho,
          gamma = gamma,
          alpha = 0.0,
          ncores = ncores,
          verbose = verbose,
          maxit_admm = 1e5,
          thr_em = 1e-5,
          thr_admm = 1e-6
        ),
        silent = TRUE
      )
      
      fitList <- fit_pofggm[[count]]$fit_list
      out.pfggm.rho.mle <- vector("list", length = length(fitList))
      for (kk in seq_along(fitList)) {
        out.pfggm.rho <- fitList[[kk]]
        out.pfggm.rho.mle[[kk]] <- pofggm(
          id_pobs = id_pobs,
          id_obs = id_obs,
          X = X_obs,
          Phi = Phi_emp,
          tp = tp,
          Sgm.hat = out.pfggm.rho$Sgm.hat,
          Tht.hat = out.pfggm.rho$Tht.hat,
          wTht = (out.pfggm.rho$Tht.hat == 0) * .Machine$double.xmax,
          gamma = gamma,
          alpha = out.pfggm.rho$alpha_opt,
          rho = 1e-12,
          ncores = 4L,
          verbose = verbose
        )
      }
      fit_pofggm_mle[[count]] <- out.pfggm.rho.mle
      
      idx <- seq_along(fitList)
      gamma_ebic <- 0.5
      
      aic_vector <- vapply(idx, function(i) compute_aic(n, S = fitList[[i]]$S, Theta = fitList[[i]]$Tht.hat), numeric(1))
      aic_vector.mle <- vapply(idx, function(i) compute_aic(n, S = out.pfggm.rho.mle[[i]]$S, Theta = out.pfggm.rho.mle[[i]]$Tht.hat), numeric(1))
      
      bic_vector <- vapply(idx, function(i) compute_bic(n, S = fitList[[i]]$S, Theta = fitList[[i]]$Tht.hat), numeric(1))
      bic_vector.mle <- vapply(idx, function(i) compute_bic(n, S = out.pfggm.rho.mle[[i]]$S, Theta = out.pfggm.rho.mle[[i]]$Tht.hat), numeric(1))
      
      ebic_vector <- vapply(idx, function(i) compute_ebic(gamma_ebic, n, S = fitList[[i]]$S, Theta = fitList[[i]]$Tht.hat), numeric(1))
      ebic_vector.mle <- vapply(idx, function(i) compute_ebic(gamma_ebic, n, S = out.pfggm.rho.mle[[i]]$S, Theta = out.pfggm.rho.mle[[i]]$Tht.hat), numeric(1))
      
      df_ic[[count]] <- data.frame(
        gamma1 = fit_pofggm[[count]]$grid_rho,
        AIC = aic_vector,
        AIC_MLE = aic_vector.mle,
        BIC = bic_vector,
        BIC_MLE = bic_vector.mle,
        eBIC = ebic_vector,
        eBIC_MLE = ebic_vector.mle
      )
      
      # id.opt <- which.min(ebic_vector.mle)
      # selected_fit <- out.pfggm.rho.mle[[id.opt]]
      # 
      # gamma1_min <- df_ic[[count]]$gamma1[id.opt]
      # 
      # ggplot(df_ic[[count]], aes(x = gamma1, y = eBIC_MLE)) +
      #   geom_line(linewidth = 0.8, color = "gray") +
      #   geom_point(size = 2, color = "gray") +
      #   geom_vline(xintercept = gamma1_min, linetype = "dashed", color = "gray20", linewidth = 1) +
      #   labs(x = expression(gamma[1]), y = "eBIC") +
      #   theme_bw(base_size = 12) +
      #   theme(
      #     panel.grid.minor = element_blank(),
      #     panel.grid.major.x = element_blank(),
      #     axis.title = element_text(face = "bold"),
      #     axis.text = element_text(color = "black")
      #   )
    }
    
    df_ic_full <- do.call("rbind", df_ic) |> 
      dplyr::mutate(gamma2 = rep(c(0, .25, .5, .75, 1), each = length(perc_rho)))
    id_opt <- which.min(df_ic_full$eBIC_MLE)
    min_point <- df_ic_full[id_opt, ]
    
    ggplot(df_ic_full, aes(x = gamma1, y = gamma2, z = eBIC_MLE)) +
      geom_contour_filled() +
      geom_point(
        data = min_point,
        aes(x = gamma1, y = gamma2),
        color = "red",
        shape = 17,
        size = 4
      ) +
      geom_text(
        data = min_point,
        aes(x = gamma1, y = gamma2, label = paste0("min = ", round(eBIC_MLE, 3))),
        color = "red",
        hjust = -0.1,
        vjust = -0.5
      ) +
      labs(x = expression(gamma[1]), y = expression(gamma[2]), title = "Curve di livello") +
      theme_bw()
    
    selected_models_id <- c(1, which.min(sapply(df_ic, \(z) min(z$eBIC_MLE))))
    selected_fit_it <- sapply(df_ic, \(z) which.min(z$eBIC_MLE))[selected_models_id]
    
    pofggm_gamma2_0 <- fit_pofggm[[selected_models_id[[1]]]]$fit_list#[[selected_fit_it[[1]]]]
    pofggm_best <- fit_pofggm[[selected_models_id[[2]]]]$fit_list#[[selected_fit_it[[2]]]]
    
    if (plot_it) {
      par(mfrow = c(1, 1))
      id_plot <- id_pobs[seq_len(min(2, length(id_pobs)))]
      
      matplot(tp, X_full[, id_plot, 1], type = "l", lty = 1, col = "pink2", lwd = 1.5)
      matpoints(tp, X_obs[, id_plot, 1], pch = 16, col = "black")
      
      matlines(tp, pofggm_gamma2_0[[selected_fit_it[[1]]]]$Ximputed[, id_plot, 1], col = "red2", lty = 2)
      matlines(tp, pofggm_best[[selected_fit_it[[2]]]]$Ximputed[, id_plot, 1], col = "blue2", lty = 2)
    }
    
    # ----------------------------------------------------------
    # 9. Precision-matrix estimation metrics
    # ----------------------------------------------------------
    
    # theta_err_pofggm_gamma2_0 <- compute_theta_error_path_pofggm(
    #   fit_list = pofggm_gamma2_0,
    #   theta_true_array = out_rPar$Tht,
    #   normalize = "fro"
    # )
    # storage$theta_err_mat[, isim] <- theta_err_pofggm_gamma2_0
    
    theta_err_pofggm_gamma2_0 <- mean_theta_error_array(
      theta_hat_array = pofggm_gamma2_0[[selected_fit_it[1]]]$Tht,
      theta_true_array = out_rPar$Tht,
      normalize = "fro"
    )
    storage$theta_err_mat[, isim] <- theta_err_pofggm_gamma2_0
    
    # theta_err_pofggm_best <- compute_theta_error_path_pofggm(
    #   fit_list = pofggm_best,
    #   theta_true_array = out_rPar$Tht,
    #   normalize = "fro"
    # )
    # storage$theta_err_mat_2[, isim] <- theta_err_pofggm_best
    
    theta_err_pofggm_best <- mean_theta_error_array(
      theta_hat_array = pofggm_best[[selected_fit_it[2]]]$Tht,
      theta_true_array = out_rPar$Tht,
      normalize = "fro"
    )
    storage$theta_err_mat_2[, isim] <- theta_err_pofggm_best

    # ----------------------------------------------------------
    # 10. Graph-recovery metrics
    # ----------------------------------------------------------
    
    res_pofggm_gamma2_0 <- compute_auc_graph_path_pofggm(
      theta_true_array = out_rPar$Tht,
      fit_list = pofggm_gamma2_0
    )
    storage$auc_theta_vec[isim] <- res_pofggm_gamma2_0$AUC_ROC
    
    res_pofggm_best <- compute_auc_graph_path_pofggm(
      theta_true_array = out_rPar$Tht,
      fit_list = pofggm_best
    )
    storage$auc_theta_vec_2[isim] <- res_pofggm_best$AUC_ROC
    
    # ----------------------------------------------------------
    # 11. Curve reconstruction metrics
    # ----------------------------------------------------------
    
    # curve_err_pofggm_gamma2_0 <- sapply(pofggm_gamma2_0, function(fit) {
    #   compute_imputation_error_missing(
    #     X_hat = fit$Ximputed,
    #     X_true = X_full,
    #     X_obs = X_obs,
    #     grid = tp,
    #     relative = TRUE
    #   )
    # })
    # storage$curve_err_mat[, isim] <- curve_err_pofggm_gamma2_0
    
    curve_err_pofggm_gamma2_0 <- compute_imputation_error_missing(
      X_hat = pofggm_gamma2_0[[selected_fit_it[[1]]]]$Ximputed,
      X_true = X_full,
      X_obs = X_obs,
      grid = tp,
      relative = TRUE
    )
    storage$curve_err_mat[, isim] <- curve_err_pofggm_gamma2_0
    
    # curve_err_pofggm_best <- sapply(pofggm_best, function(fit) {
    #   compute_imputation_error_missing(
    #     X_hat = fit$Ximputed,
    #     X_true = X_full,
    #     X_obs = X_obs,
    #     grid = tp,
    #     relative = TRUE
    #   )
    # })
    # storage$curve_err_mat_2[, isim] <- curve_err_pofggm_best
    
    curve_err_pofggm_best <- compute_imputation_error_missing(
      X_hat = pofggm_best[[selected_fit_it[[2]]]]$Ximputed,
      X_true = X_full,
      X_obs = X_obs,
      grid = tp,
      relative = TRUE
    )
    storage$curve_err_mat_2[, isim] <- curve_err_pofggm_best

    # ----------------------------------------------------------
    # 12. Intermediate progress report
    # ----------------------------------------------------------
    
    if (isim %% 10L == 0L) {
      temp <- rbind(
        "ThetaError_POFGGM" = rowMeans(storage$theta_err_mat[, seq_len(isim), drop = FALSE], na.rm = TRUE),
        "ThetaError_POFGGM_2" = rowMeans(storage$theta_err_mat_2[, seq_len(isim), drop = FALSE], na.rm = TRUE),
        "AUC_POFGGM" = mean(storage$auc_theta_vec[seq_len(isim)], na.rm = TRUE),
        "AUC_POFGGM_2" = mean(storage$auc_theta_vec_2[seq_len(isim)], na.rm = TRUE),
        "CurveError_POFGGM" = rowMeans(storage$curve_err_mat[, seq_len(isim), drop = FALSE], na.rm = TRUE),
        "CurveError_POFGGM_2" = rowMeans(storage$curve_err_mat_2[, seq_len(isim), drop = FALSE], na.rm = TRUE),
        "ComputationalTime_POFGGM" = mean(storage$comp_time_vec[seq_len(isim)], na.rm = TRUE),
        "ComputationalTime_POFGGM_2" = mean(storage$comp_time_vec_2[seq_len(isim)], na.rm = TRUE)
      )
      print(round(temp[, seq(1, n_rho, by = 2L)], digits = 4L))
      save.image(paste0("~/Downloads/jcgs_simul_gamma_config", iconfig, ".Rdata"))
    }
  }
}
