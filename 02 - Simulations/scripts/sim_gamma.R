# source("01 - Code/helper.R")

# ============================================================
# Simulation design and configuration setup
# ============================================================

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
  perc_theta_share = c(1.0, 0.5),
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
  storage$theta_err_mat_3 <- storage$theta_err_mat_2 <- storage$theta_err_mat
  storage$auc_theta_vec_3 <- storage$auc_theta_vec_2 <- storage$auc_theta_vec
  storage$curve_err_mat_3 <- storage$curve_err_mat_2 <- storage$curve_err_mat
  storage$comp_time_vec_3 <- storage$comp_time_vec_2 <- storage$comp_time_vec
  
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
    s2 = s2
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
    
    fit_pofggm <- try(
      fit_pofggm_path(
        X = X_obs,
        Phi_emp = Phi_emp,
        tp = tp,
        id_pobs = id_pobs,
        id_obs = id_obs,
        perc_rho = perc_rho,
        gamma = 0.0,
        alpha = 0.0,
        ncores = ncores,
        verbose = verbose,
        maxit_admm = 1e5,
        thr_em = 1e-5,
        thr_admm = 1e-6
      ),
      silent = TRUE
    )
    
    if (inherits(fit_pofggm, "try-error")) {
      warning("POFGGM path failed at simulation ", isim)
      next
    }
    
    pofggm_rho <- fit_pofggm$fit_list
    rho_max <- fit_pofggm$rho_max
    grid_rho <- fit_pofggm$grid_rho
    storage$comp_time_vec[isim] <- fit_pofggm$elapsed
    
    # ----------------------------------------------------------
    
    fit_pofggm_2 <- try(
      fit_pofggm_path(
        X = X_obs,
        Phi_emp = Phi_emp,
        tp = tp,
        id_pobs = id_pobs,
        id_obs = id_obs,
        perc_rho = perc_rho,
        gamma = 0.5,
        alpha = 0.0,
        ncores = ncores,
        verbose = verbose,
        maxit_admm = 1e5,
        thr_em = 1e-5,
        thr_admm = 1e-6
      ),
      silent = TRUE
    )
    
    if (inherits(fit_pofggm_2, "try-error")) {
      warning("POFGGM path failed at simulation ", isim)
      next
    }
    
    pofggm_rho_2 <- fit_pofggm_2$fit_list
    rho_max_2 <- fit_pofggm_2$rho_max
    grid_rho_2 <- fit_pofggm_2$grid_rho
    storage$comp_time_vec_2[isim] <- fit_pofggm_2$elapsed
    
    # ----------------------------------------------------------
    
    fit_pofggm_3 <- try(
      fit_pofggm_path(
        X = X_obs,
        Phi_emp = Phi_emp,
        tp = tp,
        id_pobs = id_pobs,
        id_obs = id_obs,
        perc_rho = perc_rho,
        gamma = 1.0,
        alpha = 0.00,
        ncores = ncores,
        verbose = verbose,
        maxit_admm = 1e5,
        thr_em = 1e-5,
        thr_admm = 1e-6
      ),
      silent = TRUE
    )
    
    if (inherits(fit_pofggm_3, "try-error")) {
      warning("POFGGM path failed at simulation ", isim)
      next
    }
    
    pofggm_rho_3 <- fit_pofggm_3$fit_list
    rho_max_3 <- fit_pofggm_3$rho_max
    grid_rho_3 <- fit_pofggm_3$grid_rho
    storage$comp_time_vec_3[isim] <- fit_pofggm_3$elapsed
    
    # ----------------------------------------------------------
    # 6. Oracle competitor: graph estimation from full curves
    # ----------------------------------------------------------
    
    fit_oracle <- estimate_theta_path_from_curves(
      X_array = X_full,
      Phi = Phi_emp,
      grid = tp,
      perc_rho = perc_rho,
      wTht = pofggm_rho[[1]]$wTht,
      alpha = 0,
      pendiag = FALSE,
      maxit = 1e5,
      thr = 1e-6,
      trace = 0,
      verbose = verbose
    )
    
    if (!fit_oracle$ok) {
      warning("Oracle path failed at simulation ", isim)
      next
    }
    
    THT_O <- fit_oracle$theta_path
    
    # ----------------------------------------------------------
    # 7. Kraus competitor: impute curves, then estimate graph
    # ----------------------------------------------------------
    
    kraus_fit <- try(
      reconstruct_curves_kraus(X_obs),
      silent = TRUE
    )
    
    if (inherits(kraus_fit, "try-error")) {
      warning("Kraus curve reconstruction failed at simulation ", isim)
      next
    }
    
    result_kraus <- kraus_fit$result_list
    XreconsKrauss <- kraus_fit$X_reconstructed
    
    fit_kraus <- estimate_theta_path_from_curves(
      X_array = XreconsKrauss,
      Phi = Phi_emp,
      grid = tp,
      perc_rho = perc_rho,
      wTht = pofggm_rho[[1]]$wTht,
      alpha = 0,
      pendiag = FALSE,
      maxit = 1e5,
      thr = 1e-6,
      trace = 0,
      verbose = verbose
    )
    
    if (!fit_kraus$ok) {
      warning("Kraus path failed at simulation ", isim)
      next
    }
    
    THT_K <- fit_kraus$theta_path
    
    # ----------------------------------------------------------
    # 8. Optional plotting
    # ----------------------------------------------------------
    
    if (plot_it) {
      par(mfrow = c(1, 1))
      id_plot <- id_pobs[seq_len(min(2, length(id_pobs)))]
      
      matplot(tp, X_full[, id_plot, 1], type = "l", lty = 1, col = "pink2", lwd = 1.5)
      matpoints(tp, X_obs[, id_plot, 1], pch = 16, col = "black")
      
      matlines(tp, pofggm_rho[[1]]$Ximputed[, id_plot, 1], col = "red2", lty = 2)
      if (length(pofggm_rho) >= 6)  matlines(tp, pofggm_rho[[6]]$Ximputed[, id_plot, 1], col = "blue2", lty = 2)
      if (length(pofggm_rho) >= 11) matlines(tp, pofggm_rho[[11]]$Ximputed[, id_plot, 1], col = "green2", lty = 2)
      if (length(pofggm_rho) >= 16) matlines(tp, pofggm_rho[[16]]$Ximputed[, id_plot, 1], col = "purple", lty = 2)
      if (length(pofggm_rho) >= 21) matlines(tp, pofggm_rho[[21]]$Ximputed[, id_plot, 1], col = "yellow2", lty = 2)
      
      matlines(tp, result_kraus[[1]]$X_reconst[, id_plot], col = "cyan2", lty = 2)
    }
    
    if (plot_it) {
      par(mfrow = c(1, 1))
      id_plot <- id_pobs[seq_len(min(2, length(id_pobs)))]
      
      matplot(tp, X_full[, id_plot, 1], type = "l", lty = 1, col = "pink2", lwd = 1.5)
      matpoints(tp, X_obs[, id_plot, 1], pch = 16, col = "black")
      
      matlines(tp, pofggm_rho_2[[1]]$Ximputed[, id_plot, 1], col = "red2", lty = 2)
      if (length(pofggm_rho) >= 6)  matlines(tp, pofggm_rho_2[[6]]$Ximputed[, id_plot, 1], col = "blue2", lty = 2)
      if (length(pofggm_rho) >= 11) matlines(tp, pofggm_rho_2[[11]]$Ximputed[, id_plot, 1], col = "green2", lty = 2)
      if (length(pofggm_rho) >= 16) matlines(tp, pofggm_rho_2[[16]]$Ximputed[, id_plot, 1], col = "purple", lty = 2)
      if (length(pofggm_rho) >= 21) matlines(tp, pofggm_rho_2[[21]]$Ximputed[, id_plot, 1], col = "yellow2", lty = 2)
      
      matlines(tp, result_kraus[[1]]$X_reconst[, id_plot], col = "cyan2", lty = 2)
    }
    
    if (plot_it) {
      par(mfrow = c(1, 1))
      id_plot <- id_pobs[seq_len(min(2, length(id_pobs)))]
      
      matplot(tp, X_full[, id_plot, 1], type = "l", lty = 1, col = "pink2", lwd = 1.5)
      matpoints(tp, X_obs[, id_plot, 1], pch = 16, col = "black")
      
      matlines(tp, pofggm_rho_3[[1]]$Ximputed[, id_plot, 1], col = "red2", lty = 2)
      if (length(pofggm_rho) >= 6)  matlines(tp, pofggm_rho_3[[6]]$Ximputed[, id_plot, 1], col = "blue2", lty = 2)
      if (length(pofggm_rho) >= 11) matlines(tp, pofggm_rho_3[[11]]$Ximputed[, id_plot, 1], col = "green2", lty = 2)
      if (length(pofggm_rho) >= 16) matlines(tp, pofggm_rho_3[[16]]$Ximputed[, id_plot, 1], col = "purple", lty = 2)
      if (length(pofggm_rho) >= 21) matlines(tp, pofggm_rho_3[[21]]$Ximputed[, id_plot, 1], col = "yellow2", lty = 2)
      
      matlines(tp, result_kraus[[1]]$X_reconst[, id_plot], col = "cyan2", lty = 2)
    }
    
    # ----------------------------------------------------------
    # 9. Precision-matrix estimation metrics
    # ----------------------------------------------------------
    
    theta_err_pofggm <- compute_theta_error_path_pofggm(
      fit_list = pofggm_rho,
      theta_true_array = out_rPar$Tht,
      normalize = "fro"
    )
    storage$theta_err_mat[, isim] <- theta_err_pofggm
    
    theta_err_pofggm_2 <- compute_theta_error_path_pofggm(
      fit_list = pofggm_rho_2,
      theta_true_array = out_rPar$Tht,
      normalize = "fro"
    )
    storage$theta_err_mat_2[, isim] <- theta_err_pofggm_2
    
    theta_err_pofggm_3 <- compute_theta_error_path_pofggm(
      fit_list = pofggm_rho_3,
      theta_true_array = out_rPar$Tht,
      normalize = "fro"
    )
    storage$theta_err_mat_3[, isim] <- theta_err_pofggm_3
    
    theta_err_kraus <- compute_theta_error_path(
      theta_hat_list = THT_K,
      theta_true_array = out_rPar$Tht,
      normalize = "fro"
    )
    storage$theta_err_kraus_mat[, isim] <- theta_err_kraus
    
    theta_err_obs <- compute_theta_error_path(
      theta_hat_list = THT_O,
      theta_true_array = out_rPar$Tht,
      normalize = "fro"
    )
    storage$theta_err_obs_mat[, isim] <- theta_err_obs
    
    # ----------------------------------------------------------
    # 10. Graph-recovery metrics
    # ----------------------------------------------------------
    
    res_pofggm <- compute_auc_graph_path_pofggm(
      theta_true_array = out_rPar$Tht,
      fit_list = pofggm_rho
    )
    storage$auc_theta_vec[isim] <- res_pofggm$AUC_ROC
    
    res_pofggm_2 <- compute_auc_graph_path_pofggm(
      theta_true_array = out_rPar$Tht,
      fit_list = pofggm_rho_2
    )
    storage$auc_theta_vec_2[isim] <- res_pofggm_2$AUC_ROC
    
    res_pofggm_3 <- compute_auc_graph_path_pofggm(
      theta_true_array = out_rPar$Tht,
      fit_list = pofggm_rho_3
    )
    storage$auc_theta_vec_3[isim] <- res_pofggm_3$AUC_ROC
    
    res_kraus <- compute_auc_graph_path(
      theta_true_array = out_rPar$Tht,
      theta_hat_list = THT_K
    )
    storage$auc_theta_kraus_vec[isim] <- res_kraus$AUC_ROC
    
    res_obs <- compute_auc_graph_path(
      theta_true_array = out_rPar$Tht,
      theta_hat_list = THT_O
    )
    storage$auc_theta_obs_vec[isim] <- res_obs$AUC_ROC
    
    # ----------------------------------------------------------
    # 11. Curve reconstruction metrics
    # ----------------------------------------------------------
    
    curve_err_pofggm <- sapply(pofggm_rho, function(fit) {
      compute_imputation_error_missing(
        X_hat = fit$Ximputed,
        X_true = X_full,
        X_obs = X_obs,
        grid = tp,
        relative = TRUE
      )
    })
    storage$curve_err_mat[, isim] <- curve_err_pofggm
    
    curve_err_pofggm_2 <- sapply(pofggm_rho_2, function(fit) {
      compute_imputation_error_missing(
        X_hat = fit$Ximputed,
        X_true = X_full,
        X_obs = X_obs,
        grid = tp,
        relative = TRUE
      )
    })
    storage$curve_err_mat_2[, isim] <- curve_err_pofggm_2
    
    curve_err_pofggm_3 <- sapply(pofggm_rho_3, function(fit) {
      compute_imputation_error_missing(
        X_hat = fit$Ximputed,
        X_true = X_full,
        X_obs = X_obs,
        grid = tp,
        relative = TRUE
      )
    })
    storage$curve_err_mat_3[, isim] <- curve_err_pofggm_3
    
    curve_err_kraus <- compute_imputation_error_missing(
      X_hat = XreconsKrauss,
      X_true = X_full,
      X_obs = X_obs,
      grid = tp,
      relative = TRUE
    )
    storage$curve_err_kraus_vec[isim] <- curve_err_kraus
    
    # ----------------------------------------------------------
    # 12. Intermediate progress report
    # ----------------------------------------------------------
    
    if (isim %% 10L == 0L) {
      temp <- rbind(
        "ThetaError_POFGGM" = rowMeans(storage$theta_err_mat[, seq_len(isim), drop = FALSE], na.rm = TRUE),
        "ThetaError_POFGGM_2" = rowMeans(storage$theta_err_mat_2[, seq_len(isim), drop = FALSE], na.rm = TRUE),
        "ThetaError_POFGGM_3" = rowMeans(storage$theta_err_mat_3[, seq_len(isim), drop = FALSE], na.rm = TRUE),
        "ThetaError_Kraus" = rowMeans(storage$theta_err_kraus_mat[, seq_len(isim), drop = FALSE], na.rm = TRUE),
        "ThetaError_Oracle" = rowMeans(storage$theta_err_obs_mat[, seq_len(isim), drop = FALSE], na.rm = TRUE),
        "AUC_POFGGM" = mean(storage$auc_theta_vec[seq_len(isim)], na.rm = TRUE),
        "AUC_POFGGM_2" = mean(storage$auc_theta_vec_2[seq_len(isim)], na.rm = TRUE),
        "AUC_POFGGM_3" = mean(storage$auc_theta_vec_3[seq_len(isim)], na.rm = TRUE),
        "AUC_Kraus" = mean(storage$auc_theta_kraus_vec[seq_len(isim)], na.rm = TRUE),
        "AUC_Oracle" = mean(storage$auc_theta_obs_vec[seq_len(isim)], na.rm = TRUE),
        "CurveError_POFGGM" = rowMeans(storage$curve_err_mat[, seq_len(isim), drop = FALSE], na.rm = TRUE),
        "CurveError_POFGGM_2" = rowMeans(storage$curve_err_mat_2[, seq_len(isim), drop = FALSE], na.rm = TRUE),
        "CurveError_POFGGM_3" = rowMeans(storage$curve_err_mat_3[, seq_len(isim), drop = FALSE], na.rm = TRUE),
        "CurveError_Kraus" = mean(storage$curve_err_kraus_vec[seq_len(isim)], na.rm = TRUE),
        "ComputationalTime_POFGGM" = mean(storage$comp_time_vec[seq_len(isim)], na.rm = TRUE),
        "ComputationalTime_POFGGM_2" = mean(storage$comp_time_vec_2[seq_len(isim)], na.rm = TRUE),
        "ComputationalTime_POFGGM_3" = mean(storage$comp_time_vec_3[seq_len(isim)], na.rm = TRUE)
      )
      print(round(temp[, seq(1, n_rho, by = 2L)], digits = 4L))
      save.image(paste0("~/Downloads/jcgs_simul_gamma_config", iconfig, ".Rdata"))
    }
  }
}
