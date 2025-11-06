library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(fields)
library(patchwork)

source (file = "../01 - Code/pofggm.R")
Rcpp::sourceCpp("../01 - Code/pofggm.cpp")

build_weekly_matrices <- function(
    electr_file,
    airtemp_file,
    holidays_file,
    sep_holidays = ";",
    date_format_electr = "%d.%m.%Y",  # es: "31.12.2009"
    date_format_holidays = "%d.%m.%Y",
    date_col_airtemp_fmt = "%Y%m%d",  # es: 20120331
    remove_sundays = TRUE,
    remove_holidays = TRUE,
    drop_all_na_columns = TRUE,
    drop_first_last_weeks = TRUE,
    expect_24_hours = TRUE
) {
  library(data.table)
  
  ## ---- Lettura dati ----
  electr <- data.table::fread(electr_file, colClasses = list(character = "Date"))
  stopifnot(all(c("Date", "Price", "Demand", "Wind", "Netimport") %in% names(electr)))
  
  electr[, Date := as.Date(Date, format = date_format_electr)]
  data.table::setorder(electr, Date)
  
  ## Ora per ogni giorno creo l'indice ora 1..24 (assumo 24 righe/giorno)
  electr[, hour := seq_len(.N), by = Date]
  if (expect_24_hours) {
    mx <- electr[, max(hour)]
    if (!identical(mx, 24L)) {
      warning(sprintf("Trovate %d righe nel giorno più lungo, non 24. Imposta expect_24_hours=FALSE se voluto.", mx))
    }
  }
  
  ## Inizio settimana (lunedì) e giorno della settimana indipendente dalla locale
  electr[, week_start := as.Date(cut(Date, "week", start.on.monday = TRUE))]
  electr[, wday := as.integer(strftime(Date, "%u"))]  # 1=Mon ... 7=Sun
  
  ## ---- Festivi ----
  holidays <- data.table::fread(holidays_file, sep = sep_holidays)
  ## Primo campo: "Holidays"
  hol_vec <- as.Date(as.character(holidays[[1L]]), format = date_format_holidays)
  
  if (remove_holidays) {
    ## Mantengo la struttura ma nollo i valori (così gli indici restano coerenti)
    electr[Date %in% hol_vec, c("Price","Demand","Wind","Netimport") := NA_real_]
  }
  
  ## ---- Domeniche ----
  if (remove_sundays) {
    electr <- electr[wday <= 6L]  # rimuovo completamente le domeniche
  }
  
  ## ---- Temperatura: join by Date ----
  air <- data.table::fread(airtemp_file, colClasses = list(character = 1))
  if (ncol(air) < 2) stop("Il file temperatura deve avere almeno 2 colonne (data, temperatura).")
  data.table::setnames(air, old = names(air)[1:2], new = c("date_raw","tempr"))
  air[, date := as.Date(date_raw, format = date_col_airtemp_fmt)]
  air[, date_raw := NULL]
  
  ## Join: assegno temperatura per data (replicata per le 24 ore)
  electr <- air[electr, on = .(date = Date)]
  data.table::setnames(electr, "tempr", "Z")
  
  ## ---- Pivot verso matrici [m, n] ----
  weeks <- sort(unique(electr$week_start))
  n <- length(weeks)
  m <- if (remove_sundays) 24L * 6L else 24L * 7L
  
  ## Indice di riga (1..m) e colonna (1..n)
  electr[, row_idx := (wday - 1L) * 24L + hour]
  if (remove_sundays && electr[, any(wday == 7L)]) {
    stop("Sono presenti domeniche non rimosse: controlla remove_sundays.")
  }
  
  electr[, col_idx := match(week_start, weeks)]
  
  ## Matrici preallocate
  Y.mat <- matrix(NA_real_, nrow = m, ncol = n)
  D.mat <- matrix(NA_real_, nrow = m, ncol = n)
  W.mat <- matrix(NA_real_, nrow = m, ncol = n)
  N.mat <- matrix(NA_real_, nrow = m, ncol = n)
  Z.mat <- matrix(NA_real_, nrow = m, ncol = n)
  
  keep <- electr$row_idx >= 1L & electr$row_idx <= m
  rr <- electr$row_idx[keep]
  cc <- electr$col_idx[keep]
  
  Y.mat[cbind(rr, cc)] <- electr$Price[keep]
  D.mat[cbind(rr, cc)] <- electr$Demand[keep]
  W.mat[cbind(rr, cc)] <- electr$Wind[keep]
  N.mat[cbind(rr, cc)] <- electr$Netimport[keep]
  Z.mat[cbind(rr, cc)] <- electr$Z[keep]
  
  ## ---- Drop colonne completamente non osservate ----
  if (drop_all_na_columns) {
    keep_cols <- which(colSums(!is.na(Y.mat)) > 0L)
    if (length(keep_cols) < n) {
      Y.mat <- Y.mat[, keep_cols, drop = FALSE]
      D.mat <- D.mat[, keep_cols, drop = FALSE]
      W.mat <- W.mat[, keep_cols, drop = FALSE]
      N.mat <- N.mat[, keep_cols, drop = FALSE]
      Z.mat <- Z.mat[, keep_cols, drop = FALSE]
      weeks <- weeks[keep_cols]
      n <- length(weeks)
    }
  }
  
  ## ---- (Opzionale) rimuovi prima/ultima settimana (parziali) ----
  if (drop_first_last_weeks && n >= 2) {
    Y.mat <- Y.mat[, -c(1L, ncol(Y.mat)), drop = FALSE]
    D.mat <- D.mat[, -c(1L, ncol(D.mat)), drop = FALSE]
    W.mat <- W.mat[, -c(1L, ncol(W.mat)), drop = FALSE]
    N.mat <- N.mat[, -c(1L, ncol(N.mat)), drop = FALSE]
    Z.mat <- Z.mat[, -c(1L, ncol(Z.mat)), drop = FALSE]
    weeks <- weeks[-c(1L, length(weeks))]
    n <- length(weeks)
  }
  
  ## ---- Sanity checks e report ----
  msg_na <- function(M) sprintf("%d NA (%.2f%%)", sum(is.na(M)), 100*mean(is.na(M)))
  message(sprintf("Y: %s | D: %s | W: %s | N: %s | Z: %s",
                  msg_na(Y.mat), msg_na(D.mat), msg_na(W.mat), msg_na(N.mat), msg_na(Z.mat)))
  
  ## Ritorno
  list(
    Y.mat = Y.mat,
    W.mat = W.mat,
    D.mat = D.mat,
    N.mat = N.mat,
    Z.mat = Z.mat,
    weeks = weeks,
    m = nrow(Y.mat),
    n = ncol(Y.mat)
  )
}

root <- "../03 - Analysis/"

res <- build_weekly_matrices(
  electr_file = file.path(root, "DATA", "DATA_electricity.csv"),
  airtemp_file = file.path(root, "DATA", "DATA_airtemp.csv"),
  holidays_file = file.path(root, "DATA", "DATA_holidays.csv"),
  sep_holidays = ";",
  remove_sundays = TRUE,
  remove_holidays = TRUE,
  drop_all_na_columns = TRUE,
  drop_first_last_weeks = TRUE
)

Y.mat <- res$Y.mat; W.mat <- res$W.mat; D.mat <- res$D.mat; N.mat <- res$N.mat; Z.mat <- res$Z.mat
n <- res$n; m <- res$m
weeks <- res$weeks

## ------------------------------------------------------------
## Building array X = (m x n x p)
## ------------------------------------------------------------

make_X <- function(Y.mat, D.mat, W.mat, N.mat) {
  stopifnot(identical(dim(Y.mat), dim(D.mat)),
            identical(dim(Y.mat), dim(W.mat)),
            identical(dim(Y.mat), dim(N.mat)))
  m <- nrow(Y.mat); n <- ncol(Y.mat); p <- 4L
  X <- array(NA_real_, dim = c(m, n, p),
             dimnames = list(NULL, NULL, c("Price", "Demand", "Wind", "Netimport")))
  X[,,1] <- Y.mat
  X[,,2] <- D.mat
  X[,,3] <- W.mat
  X[,,4] <- N.mat
  X
}

X <- make_X(Y.mat, D.mat, W.mat, N.mat)
d <- dim(X)[1]; n <- dim(X)[2]; p <- dim(X)[3]
tp <- seq(0, 1, length.out = d)

## ------------------------------------------------------------
## Smoothing of X along time
## ------------------------------------------------------------

smooth_X <- function(X,
                     tp = seq(0, 1, length.out = dim(X)[1]),
                     spar = NULL, df = NULL, cv = FALSE,
                     fill_missing = FALSE,           # TRUE => predice su tutta la griglia
                     min_points = 6L,                # minimo di punti osservati per spline
                     fallback = c("approx", "none"), # fallback se pochi punti
                     ncores = 1L) {
  fallback <- match.arg(fallback)
  stopifnot(is.array(X), length(dim(X)) == 3L)
  m <- dim(X)[1]; n <- dim(X)[2]; p <- dim(X)[3]
  stopifnot(length(tp) == m)
  
  ## Helper: smoothing di una singola curva y (lunghezza m)
  smooth_one <- function(y) {
    idx <- !is.na(y)
    if (sum(idx) < max(4L, min_points) || length(unique(tp[idx])) < 4L) {
      ## Non abbastanza informazione per spline
      if (fallback == "approx") {
        if (sum(idx) >= 2L) {
          if (fill_missing) {
            y_hat <- approx(x = tp[idx], y = y[idx], xout = tp, ties = "ordered", rule = 2)$y
            return(y_hat)
          } else {
            y_out <- y
            y_out[idx] <- approx(x = tp[idx], y = y[idx], xout = tp[idx], ties = "ordered", rule = 2)$y
            return(y_out)
          }
        } else {
          ## 0 o 1 punto → ritorno l'originale
          return(y)
        }
      } else {
        return(y)
      }
    }
    
    ## Costruisco la chiamata a smooth.spline senza passare argomenti NULL
    ss <- tryCatch({
      if (!is.null(spar)) {
        stats::smooth.spline(x = tp[idx], y = y[idx], spar = spar, cv = cv)
      } else if (!is.null(df)) {
        stats::smooth.spline(x = tp[idx], y = y[idx], df = df, cv = cv)
      } else {
        stats::smooth.spline(x = tp[idx], y = y[idx], cv = cv)
      }
    }, error = function(e) NULL)
    
    if (is.null(ss)) {
      ## fallback in caso di errore runtime di smooth.spline
      if (fallback == "approx" && sum(idx) >= 2L) {
        if (fill_missing) {
          return(approx(x = tp[idx], y = y[idx], xout = tp, ties = "ordered", rule = 2)$y)
        } else {
          y_out <- y
          y_out[idx] <- approx(x = tp[idx], y = y[idx], xout = tp[idx], ties = "ordered", rule = 2)$y
          return(y_out)
        }
      } else {
        return(y)
      }
    }
    
    if (fill_missing) {
      stats::predict(ss, x = tp)$y
    } else {
      y_out <- y
      y_out[idx] <- stats::predict(ss, x = tp[idx])$y
      y_out
    }
  }
  
  ## Applico per ogni variabile (terza dimensione) e per ogni colonna (settimana)
  Xsmooth <- array(NA_real_, dim = c(m, n, p), dimnames = dimnames(X))
  
  ## --- Strategia di parallelizzazione leggera (opzionale) ---
  use_mc <- ncores > 1L && .Platform$OS.type != "windows" && requireNamespace("parallel", quietly = TRUE)
  if (use_mc) {
    parallel <- getNamespace("parallel")
  }
  
  for (j in seq_len(p)) {
    M <- X[,,j, drop = TRUE]  # m x n
    
    sm_cols <- if (!use_mc) {
      ## vapply su ogni colonna: veloce e senza allocazioni inutili
      vapply(seq_len(n), function(col) smooth_one(M[, col]), numeric(m))
    } else {
      ## mclapply su colonne, poi cbind
      cols_list <- parallel::mclapply(seq_len(n), function(col) smooth_one(M[, col]),
                                      mc.cores = ncores)
      do.call(cbind, cols_list)
    }
    
    ## vapply/ mclapply sopra ritorna m x n (colonne)
    Xsmooth[,,j] <- sm_cols
  }
  
  Xsmooth
}

Xsmooth <- smooth_X(
  X,
  tp = tp,
  spar = NULL, df = NULL, cv = FALSE,  # lascia scegliere a GCV; oppure imposta tu spar/df
  fill_missing = FALSE,                 # come nel tuo codice: non imputo i missing
  min_points = 6L,
  fallback = "approx",                  # in caso di pochi punti: lineare
  ncores = 1L                           # aumenta su Linux/macOS se vuoi
)

op <- par(mfrow = c(2,2), mar = c(3,3,2,1), mgp = c(1.6,0.6,0))
on.exit(par(op), add = TRUE)

plot(tp, X[,1,1], type = "l", main = "Price (settimana 1)")
lines(tp, Xsmooth[,1,1], col = "red")

plot(tp, X[,1,2], type = "l", main = "Demand (settimana 1)")
lines(tp, Xsmooth[,1,2], col = "red")

plot(tp, X[,1,3], type = "l", main = "Wind (settimana 1)")
lines(tp, Xsmooth[,1,3], col = "red")

plot(tp, X[,1,4], type = "l", main = "Netimport (settimana 1)")
lines(tp, Xsmooth[,1,4], col = "red")

## ------------------------------------------------------------
## Media di Kraus (per riga), ignora NA
## ------------------------------------------------------------
meanKraus <- function(X_mat) {
  # X_mat: (d x n)
  rowMeans(X_mat, na.rm = TRUE)
}

## ------------------------------------------------------------
## Covarianza di Kraus vettorializzata
## ------------------------------------------------------------
covKraus <- function(X_mat, mu = NULL) {
  # X_mat: (d x n)  con NA
  # mu:    (d)      media per riga; se NULL usa meanKraus(X_mat)
  d <- nrow(X_mat)
  
  if (is.null(mu)) mu <- meanKraus(X_mat)
  stopifnot(length(mu) == d)
  
  # Centra, marca i NA, conta coperture valide e accumula prodotti
  X_cent <- sweep(X_mat, 1L, mu, FUN = "-")
  idMiss <- is.na(X_cent)
  
  # valid_counts[s,t] = #colonne in cui entrambe (s,t) sono osservate
  valid_counts <- tcrossprod(!idMiss)
  
  # Sostituisci NA con 0: la crossprod somma solo dove non NA
  X_cent[idMiss] <- 0
  prod_mat <- tcrossprod(X_cent)  # somma dei prodotti
  
  # Stima media dei prodotti (covarianza) ignorando NA
  covKraus_mat <- prod_mat / valid_counts
  covKraus_mat[valid_counts == 0] <- NA_real_
  
  attr(covKraus_mat, "cov.count") <- valid_counts
  covKraus_mat
}

## ------------------------------------------------------------
## Smoothing della media campionaria per canale
## ------------------------------------------------------------
smooth_means <- function(mu_sample, tp,
                         spar = NULL, df = NULL, cv = FALSE,
                         min_points = 6L, fallback = c("approx","none")) {
  # mu_sample: (d x p) media grezza per canale
  fallback <- match.arg(fallback)
  d <- nrow(mu_sample); p <- ncol(mu_sample)
  mu_smooth <- matrix(NA_real_, nrow = d, ncol = p)
  
  smooth_one <- function(y) {
    ok <- is.finite(y)
    if (sum(ok) < max(4L, min_points)) {
      if (fallback == "approx" && sum(ok) >= 2L) {
        return(approx(tp[ok], y[ok], xout = tp, rule = 2)$y)
      } else {
        return(y)  # non tocco
      }
    }
    ss <- tryCatch({
      if (!is.null(spar)) stats::smooth.spline(tp[ok], y[ok], spar = spar, cv = cv)
      else if (!is.null(df)) stats::smooth.spline(tp[ok], y[ok], df = df, cv = cv)
      else stats::smooth.spline(tp[ok], y[ok], cv = cv)
    }, error = function(e) NULL)
    if (is.null(ss)) {
      if (fallback == "approx" && sum(ok) >= 2L) {
        approx(tp[ok], y[ok], xout = tp, rule = 2)$y
      } else y
    } else {
      stats::predict(ss, x = tp)$y
    }
  }
  
  for (j in seq_len(p)) mu_smooth[, j] <- smooth_one(mu_sample[, j])
  mu_smooth
}

## ------------------------------------------------------------
## Proiezione PSD: clip autovalori negativi
## ------------------------------------------------------------
make_psd <- function(M, eps = 1e-10) {
  M <- 0.5 * (M + t(M))  # simmetrizza
  ei <- eigen(M, symmetric = TRUE)
  lam <- pmax(ei$values, eps)
  ei$vectors %*% (lam * t(ei$vectors))
}

## ------------------------------------------------------------
## Smoothing G_i(s,t) con GAM tensoriale + pesi
## ------------------------------------------------------------
smooth_G_channel <- function(Gi, counts, tp, k = 10, method = c("REML","GCV.Cp"),
                             use_bam = FALSE, ensure_psd = TRUE) {
  method <- match.arg(method)
  d <- length(tp)
  stopifnot(identical(dim(Gi), c(d, d)), identical(dim(counts), c(d, d)))
  
  # Griglia
  row.vec <- rep(tp, each = d)
  col.vec <- rep(tp, times = d)
  
  # Dati vettorializzati e filtrati (niente NA, niente pesi 0)
  y <- as.vector(Gi)
  w <- as.vector(counts)
  keep <- is.finite(y) & is.finite(w) & (w > 0)
  df <- data.frame(y = y[keep], row = row.vec[keep], col = col.vec[keep], w = w[keep])
  
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Serve il pacchetto 'mgcv'. Esegui install.packages('mgcv').")
  }
  
  # Normalizzo i pesi per stabilità numerica (facoltativo)
  df$w <- df$w / mean(df$w)
  
  form <- y ~ te(row, col, k = c(k, k))
  
  fit <- if (use_bam) {
    mgcv::bam(form, data = df, weights = w, method = method, discrete = TRUE)
  } else {
    mgcv::gam(form, data = df, weights = w, method = method)
  }
  # browser()
  G_hat <- matrix(NA_real_, nrow = d, ncol = d)
  # Predico su tutta la griglia (in ordine coerente con row.vec/col.vec)
  pred <- stats::predict(fit, newdata = data.frame(row = row.vec, col = col.vec))
  G_hat[,] <- matrix(pred, d, d)
  
  # Simmetrizza
  G_hat <- 0.5 * (G_hat + t(G_hat))
  
  # Garantisco PSD se richiesto
  if (ensure_psd) G_hat <- make_psd(G_hat)
  
  G_hat
}

## ------------------------------------------------------------
## Pipeline: mu, G_i, G_i^smooth, H
## ------------------------------------------------------------
estimate_mean_and_cov_surfaces <- function(Xsmooth, tp,
                                           mean_spar = NULL, mean_df = NULL, mean_cv = FALSE,
                                           cov_k = 10, cov_method = "REML", use_bam = FALSE,
                                           ensure_psd_each = TRUE, ensure_psd_H = TRUE) {
  d <- dim(Xsmooth)[1]; n <- dim(Xsmooth)[2]; p <- dim(Xsmooth)[3]
  stopifnot(length(tp) == d)
  
  ## --- Media campionaria per canale (d x p)
  mu_sample <- apply(Xsmooth, 3, meanKraus)        # (d x p)
  mu_smooth <- smooth_means(mu_sample, tp, spar = mean_spar, df = mean_df, cv = mean_cv)
  
  ## --- Covarianze per canale (lista) e smoothing
  G_sample  <- vector("list", p)
  G_smooth  <- vector("list", p)
  
  for (i in seq_len(p)) {
    Gi <- covKraus(Xsmooth[,, i, drop = TRUE], mu = mu_smooth[, i])
    G_sample[[i]] <- Gi
    counts_i <- attr(Gi, "cov.count")
    G_smooth[[i]] <- smooth_G_channel(Gi, counts_i, tp,
                                      k = cov_k, method = cov_method,
                                      use_bam = use_bam, ensure_psd = ensure_psd_each)
  }
  
  ## --- Media delle covarianze lisciate
  H <- Reduce(`+`, G_smooth) / p
  if (ensure_psd_H) H <- make_psd(H)
  
  list(
    mu_sample = mu_sample,
    mu_smooth = mu_smooth,
    G_sample  = G_sample,
    G_smooth  = G_smooth,
    H_sample  = H
  )
}

# Parametri (puoi lasciarli così per iniziare)
mean_df   <- NULL     # oppure imposta df=10, p.es.
mean_spar <- NULL     # oppure spar=0.4
cov_k     <- 30       # dimensione base di ciascun margine del 'te'
cov_meth  <- "GCV.Cp"   # oppure "GCV.Cp" (REML)
use_bam   <- FALSE    # TRUE se d^2 è grande o se vuoi velocità extra

res_cov <- estimate_mean_and_cov_surfaces(
  Xsmooth, tp,
  mean_spar = mean_spar, mean_df = mean_df, mean_cv = FALSE,
  cov_k = cov_k, cov_method = cov_meth, use_bam = use_bam,
  ensure_psd_each = TRUE, ensure_psd_H = TRUE
)

mu.sample <- res_cov$mu_sample
mu.smooth <- res_cov$mu_smooth

G.sample  <- res_cov$G_sample
G.smooth  <- res_cov$G_smooth

H.sample <- Reduce(`+`, G.sample) / p
H.smooth  <- res_cov$H_sample

d <- dim(Xsmooth)[1]; n <- dim(Xsmooth)[2]; p <- dim(Xsmooth)[3]
stopifnot(identical(dim(mu.smooth), c(d, p)))

# Costruisco una "copia" di mu.smooth con n colonne (una per settimana)
mu3 <- array(mu.sample, dim = c(d, 1, p))
mu3 <- mu3[, rep(1, n), , drop = FALSE]   # diventa d x n x p

# Centro: per ogni canale j sottraggo mu.smooth[,j] a tutte le colonne
Xcentered <- Xsmooth - mu3

# Costruisco una "copia" di sd.smooth con n colonne (una per settimana)
sd.smooth <- sqrt(sapply(G.sample, diag))
sd3 <- array(sd.smooth, dim = c(d, 1, p))
sd3 <- sd3[, rep(1, n), , drop = FALSE]   # diventa d x n x p

# Scalo: per ogni canale j divido sd.smooth[,j] a tutte le colonne
Xscaled <- Xcentered / sd3

# Media residua per canale dopo il centraggio (dovrebbe essere “piccola” se mu.smooth è vicino a mu.sample)
res_rowmeans <- apply(Xscaled, 3, function(M) rowMeans(M, na.rm = TRUE))
sapply(1:dim(res_rowmeans)[2], function(j) max(abs(res_rowmeans[, j]), na.rm = TRUE))

# Nomi delle variabili
var_names <- c("Electricity price", "Electricity demand", "Wind power", "Electricity imports")

# Convertiamo Xscaled in formato long
df_long <- map_dfr(1:dim(Xscaled)[3], function(j) {
  as.data.frame(Xscaled[,,j]) |>
    mutate(Time = tp, Variable = var_names[j])
})

# Portiamo in formato tidy (escludendo Time e Variable dal pivot)
df_tidy <- df_long |>
  pivot_longer(
    cols = -c(Time, Variable),
    names_to = "Sim",
    values_to = "Value"
  ) |>
  mutate(Highlight = ifelse(Sim == "V1", "Main", "Background"),
         Variable = factor(Variable, levels = var_names))

# Grafico ggplot
# Grafico ggplot con trasparenza e ordine corretto
p <- ggplot() +
  # Linee grigie di background (prima, con alpha)
  geom_line(
    data = df_tidy %>% filter(Highlight == "Background"),
    aes(x = Time, y = Value, group = Sim),
    color = "gray70", linewidth = 0.4, alpha = 0.5
  ) +
  # Linee rosse in evidenza (poi)
  geom_line(
    data = df_tidy %>% filter(Highlight == "Main"),
    aes(x = Time, y = Value, group = Sim),
    color = "gray40", linewidth = 1.2
  ) +
  facet_wrap(~Variable, scales = "free_y") +
  labs(x = "Time", y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold")
  )
p                                    

## ------------------------------------------------------------
## quadratic weights
## ------------------------------------------------------------
quadWeights <- function(argvals, method = c("trapezoidal", "midpoint")) {
  method <- match.arg(method)
  if (any(!is.finite(argvals))) stop("'argvals' contiene NA/Inf.")
  if (is.unsorted(argvals, strictly = TRUE)) {
    stop("'argvals' deve essere strettamente crescente (niente duplicati).")
  }
  D <- length(argvals)
  if (D < 2L) stop("Servono almeno 2 punti di griglia.")
  
  if (method == "trapezoidal") {
    # 0.5 * (Δ1, Δ2+... , Δ_last)
    d <- diff(argvals)
    w <- numeric(D)
    if (D == 2L) {
      w[] <- 0.5 * d
    } else {
      w[1L] <- 0.5 * d[1L]
      w[D]  <- 0.5 * d[D - 1L]
      w[2:(D-1L)] <- 0.5 * (d[1:(D-2L)] + d[2:(D-1L)])
    }
  } else {
    # midpoint: [0, diff]
    w <- c(0, diff(argvals))
  }
  
  if (any(w <= 0)) {
    stop("Pesi non positivi; controlla 'argvals' (duplicati o non monotono).")
  }
  w
}

## ------------------------------------------------------------
## univariate functional principal component analysis
## ------------------------------------------------------------
uFPCA <- function(G, argvals, pev, efunctions_sign = NULL, 
                  method = c("trapezoidal", "midpoint"), tol = NULL, 
                  max_comp = NULL, check_orthonorm = FALSE, ...) {
  method <- match.arg(method)
  if (!is.matrix(G) || nrow(G) != ncol(G)) stop("'G' deve essere matrice quadrata.")
  d <- nrow(G)
  if (length(argvals) != d) stop("Lunghezza 'argvals' diversa da nrow(G).")
  if (!isSymmetric(G, tol = 1e-10)) {
    # Simmetrizza leggermente per sicurezza numerica
    G <- 0.5 * (G + t(G))
  }
  
  # Pesi e matrice 'V' del problema agli autovalori
  w      <- quadWeights(argvals, method = method)
  sqrt.w <- sqrt(w)
  
  # Costruzione efficiente: V = D_{sqrt.w} G D_{sqrt.w}
  # evito diag(d) usando broadcasting vettoriale
  V <- (sqrt.w * G) * rep(sqrt.w, each = d)
  
  # Unico eigen-decomposition (ordinata in decrescente per matrici simmetriche)
  eig <- eigen(V, symmetric = TRUE)
  vals <- eig$values
  vecs <- eig$vectors
  
  # Soglia per negativi numerici
  if (is.null(tol)) tol <- max(1e-12, max(vals, na.rm = TRUE) * 1e-12)
  
  # Seleziona solo gli autovalori positivi (o > tol)
  pos <- which(vals > tol)
  if (length(pos) == 0L) {
    warning("Tutti gli autovalori sono <= tol; ritorno struttura vuota.")
    return(list(evalues = numeric(), efunctions = matrix(numeric(0), d, 0),
                npc = 0L, pev = 0, cumpev = numeric(0), w = w))
  }
  vals_pos <- vals[pos]
  vecs_pos <- vecs[, pos, drop = FALSE]
  
  # Varianza spiegata cumulata sui positivi
  cumpev <- cumsum(vals_pos) / sum(vals_pos)
  
  # Scelta npc: o per pev richiesto o tutto il positivo
  if (!missing(pev) && !is.null(pev)) {
    if (pev <= 0 || pev > 1) stop("'pev' deve essere in (0,1].")
    npc <- which(cumpev >= pev)[1L]
  } else {
    npc <- length(vals_pos)
    pev <- cumpev[npc]
  }
  if (!is.null(max_comp)) npc <- min(npc, as.integer(max_comp))
  
  # Autovalori finali (primi npc)
  evalues <- vals_pos[seq_len(npc)]
  
  # Autovettori/efunzioni su G: φ = D_{-1/2} * v
  # (normalizzati in L2(w): t(phi) %*% diag(w) %*% phi = I)
  efunctions <- vecs_pos[, seq_len(npc), drop = FALSE] / sqrt.w
  
  # Allineamento segni (opzionale)
  if (!is.null(efunctions_sign)) {
    if (!is.matrix(efunctions_sign) || nrow(efunctions_sign) != d) {
      stop("'efunctions_sign' deve essere una matrice d x >= npc.")
    }
    if (ncol(efunctions_sign) < npc) {
      stop("'efunctions_sign' ha meno colonne di 'npc'.")
    }
    ref <- efunctions_sign[, seq_len(npc), drop = FALSE]
    # criterio: segno del prodotto interno pesato
    sgn <- sign(colSums(ref * efunctions * w))
    sgn[sgn == 0] <- 1
    efunctions <- sweep(efunctions, 2L, sgn, `*`)
  }
  
  # Verifica ortonormalità (facoltativa)
  if (isTRUE(check_orthonorm)) {
    # t(phi) W phi ~ I
    Wphi <- efunctions * w
    Gcheck <- crossprod(efunctions, Wphi)
    ortho_err <- norm(Gcheck - diag(npc), "F")
  } else {
    ortho_err <- NULL
  }
  
  list(
    evalues    = evalues,          # autovalori > 0 (o > tol), primi npc
    efunctions = efunctions,       # d x npc, L2(w)-ortonormali
    npc        = npc,
    pev        = pev,
    cumpev     = cumpev,           # solo sui positivi
    w          = w,
    tol        = tol,
    ortho_err  = ortho_err
  )
}

temp_ufpca <- uFPCA(G = H.smooth, argvals = tp, pev = 0.9999, check_orthonorm = TRUE)

temp_ufpca$npc
temp_ufpca$cumpev

## K e funzioni proprie dalla uFPCA su H
K <- temp_ufpca$npc
Phi.sample <- temp_ufpca$efunctions   # matrice d x K
w <- temp_ufpca$w                     # pesi di quadratura (ti servono per i punteggi)

id_pobs <- which(apply(X, 2, function(slice) anyNA(slice)))
id_obs  <- setdiff(seq_len(n), id_pobs)

## Verifica/diagnostica
length(id_pobs)
length(id_obs)

na_by_week_channel <- apply(is.na(X), c(2, 3), any)  # n x p (TRUE se la variabile j in settimana i ha NA)
id_pobs <- which(rowSums(na_by_week_channel) > 0)
id_obs  <- which(rowSums(na_by_week_channel) == 0)

par(mfrow = c(2, 2))
matplot(tp, Xscaled[,,1], type = "l", lty = 1, col = gray(.5), lwd = .75, 
        xlab = "Time", ylab = "Price")
lines(tp, Xscaled[,1,1], lty = 1, col = "red", lwd = 1.5)
matplot(tp, Xscaled[,,2], type = "l", lty = 1, col = gray(.5), lwd = .75, 
        xlab = "Time", ylab = "Demand")
lines(tp, Xscaled[,1,2], lty = 1, col = "red", lwd = 1.5)
matplot(tp, Xscaled[,,3], type = "l", lty = 1, col = gray(.5), lwd = .75, 
        xlab = "Time", ylab = "Wind")
lines(tp, Xscaled[,1,3], lty = 1, col = "red", lwd = 1.5)
matplot(tp, Xscaled[,,4], type = "l", lty = 1, col = gray(.5), lwd = .75, 
        xlab = "Time", ylab = "Import")
lines(tp, Xscaled[,1,4], lty = 1, col = "red", lwd = 1.5)

first.iteration <- pofggm(id_pobs = id_pobs, id_obs = id_obs, X = Xscaled, 
                          Phi = Phi.sample, tp = tp, 
                          Sgm.hat = NULL, Tht.hat = NULL, wTht = NULL,
                          rho = 100, ncores = 1L, verbose = TRUE)

rho.max <- first.iteration$rho.max
grid.rho <- seq(rho.max, 0.0, length = 51)

out.pfggm.rho <- vector(mode = "list")
out.pfggm.rho[[1]] <- first.iteration

# fitting path for rho
for (kk in 2:length(grid.rho)) {
  print(kk)
  out.pfggm.rho[[kk]] <- pofggm(id_pobs = id_pobs, id_obs = id_obs, X = Xscaled, 
                                Phi = Phi.sample, tp = tp, 
                                Sgm.hat = out.pfggm.rho[[kk - 1]]$Sgm.hat, 
                                Tht.hat = out.pfggm.rho[[kk - 1]]$Tht.hat, wTht = NULL,
                                rho = grid.rho[kk], ncores = 1L, verbose = TRUE)
}

out.pfggm.rho.mle <- vector(mode = "list")

# fitting path for rho
for (kk in 1:length(grid.rho)) {
  print(kk)
  out.pfggm.rho.mle[[kk]] <- pofggm(id_pobs = id_pobs, id_obs = id_obs, X = Xscaled, 
                                    Phi = Phi.sample, tp = tp, 
                                    Sgm.hat = out.pfggm.rho[[kk]]$Sgm.hat, 
                                    Tht.hat = out.pfggm.rho[[kk]]$Tht.hat, 
                                    wTht = (out.pfggm.rho[[kk]]$Tht.hat == 0) * .Machine$double.xmax,
                                    rho = 1e-12, ncores = 4L, verbose = TRUE)
}


## --------- Utilità numeriche ----------------------------------------------

# log|Theta| stabile, con jitter diag se necessario
safe_logdet <- function(Theta, jitter = 1e-10, max_jitter = 1e-3) {
  stopifnot(is.matrix(Theta), nrow(Theta) == ncol(Theta))
  Theta <- 0.5 * (Theta + t(Theta)) # simmetrizza numericamente
  try_det <- try(determinant(Theta, logarithm = TRUE), silent = TRUE)
  if (!inherits(try_det, "try-error")) {
    return(as.numeric(try_det$modulus))
  }
  # aggiunge jitter crescente sulla diagonale
  jj <- jitter
  repeat {
    Theta_j <- Theta; diag(Theta_j) <- diag(Theta_j) + jj
    try_det <- try(determinant(Theta_j, logarithm = TRUE), silent = TRUE)
    if (!inherits(try_det, "try-error") && is.finite(as.numeric(try_det$modulus))) {
      return(as.numeric(try_det$modulus))
    }
    jj <- jj * 10
    if (jj > max_jitter) break
  }
  return(NA_real_) # fallito
}

# tr(S %*% Theta) veloce (Theta simmetrica)
trace_prod_sym <- function(S, Theta) {
  # tr(S Theta) = sum_{ij} S_ij * Theta_ji = sum(S * Theta') ; se Theta simmetrica, uguale a sum(S * Theta)
  sum(S * Theta)
}

# df counting: "edges" (default), "offdiag", "entries"
df_count <- function(Theta, rule = c("edges","offdiag","entries")) {
  rule <- match.arg(rule)
  p <- nrow(Theta)
  nz_off_upper <- sum(Theta[upper.tri(Theta)] != 0)
  if (rule == "edges")   return(p + nz_off_upper)       # p diag + #archi (upper)
  if (rule == "offdiag") return(nz_off_upper)           # solo archi
  if (rule == "entries") return(sum(Theta != 0))        # tutte le celle (diag + 2*upper)
}

## --------- IC per K slice --------------------------------------------------

# Ritorna somma su K (e opzionalmente per-slice)
compute_ic <- function(n, S, Theta, type = c("AIC","BIC","EBIC"),
                       gamma_ebic = 0.5,
                       df_rule = c("edges","offdiag","entries"),
                       return_per_slice = FALSE) {
  type   <- match.arg(type)
  df_rule <- match.arg(df_rule)
  
  stopifnot(length(dim(S))    == 3L,
            length(dim(Theta))== 3L,
            S[,,1] |> is.matrix(),
            all(dim(S)[1:2] == dim(Theta)[1:2]),
            dim(S)[3] == dim(Theta)[3])
  
  p <- dim(Theta)[1]; K <- dim(Theta)[3]
  vals <- numeric(K)
  
  for (k in seq_len(K)) {
    Sk  <- S[,,k]; Tk <- Theta[,,k]
    trST    <- trace_prod_sym(Sk, Tk)
    logdetT <- safe_logdet(Tk)
    if (!is.finite(logdetT)) {
      vals[k] <- NA_real_
      next
    }
    neg2loglik <- n * (trST - logdetT)
    
    # gradi di liberta'
    df_k  <- df_count(Tk, rule = df_rule)
    
    if (type == "AIC") {
      vals[k] <- neg2loglik + 2 * df_k
    } else if (type == "BIC") {
      vals[k] <- neg2loglik + df_k * log(n)
    } else { # EBIC
      # |E| = #archi = nz_off_upper
      E_k <- df_count(Tk, rule = "offdiag")
      vals[k] <- neg2loglik + df_k * log(n) + 4 * gamma_ebic * E_k * log(p)
    }
  }
  
  if (isTRUE(return_per_slice)) return(vals)
  sum(vals, na.rm = TRUE)
}

# Wrapper comodi
compute_aic  <- function(n, S, Theta, df_rule = c("edges","offdiag","entries"),
                         return_per_slice = FALSE) {
  compute_ic(n, S, Theta, type = "AIC", df_rule = df_rule,
             return_per_slice = return_per_slice)
}
compute_bic  <- function(n, S, Theta, df_rule = c("edges","offdiag","entries"),
                         return_per_slice = FALSE) {
  compute_ic(n, S, Theta, type = "BIC", df_rule = df_rule,
             return_per_slice = return_per_slice)
}
compute_ebic <- function(gamma_ebic, n, S, Theta, df_rule = c("edges","offdiag","entries"),
                         return_per_slice = FALSE) {
  compute_ic(n, S, Theta, type = "EBIC", gamma_ebic = gamma_ebic,
             df_rule = df_rule, return_per_slice = return_per_slice)
}

idx <- seq_along(grid.rho)

# AIC
aic_vector      <- vapply(idx, function(i) compute_aic(n, S = out.pfggm.rho[[i]][["S"]],
                                                       Theta = out.pfggm.rho[[i]][["Tht.hat"]]),
                          numeric(1))
aic_vector.mle  <- vapply(idx, function(i) compute_aic(n, S = out.pfggm.rho.mle[[i]][["S"]],
                                                       Theta = out.pfggm.rho.mle[[i]][["Tht.hat"]]),
                          numeric(1))

# BIC
bic_vector      <- vapply(idx, function(i) compute_bic(n, S = out.pfggm.rho[[i]][["S"]],
                                                       Theta = out.pfggm.rho[[i]][["Tht.hat"]]),
                          numeric(1))
bic_vector.mle  <- vapply(idx, function(i) compute_bic(n, S = out.pfggm.rho.mle[[i]][["S"]],
                                                       Theta = out.pfggm.rho.mle[[i]][["Tht.hat"]]),
                          numeric(1))

# EBIC (gamma tipico 0.0–1.0; 0.5–1 più severo in alta dimensione)
gamma_ebic <- 0.5
ebic_vector     <- vapply(idx, function(i) compute_ebic(gamma_ebic, n,
                                                        S = out.pfggm.rho[[i]][["S"]],
                                                        Theta = out.pfggm.rho[[i]][["Tht.hat"]]),
                          numeric(1))
ebic_vector.mle <- vapply(idx, function(i) compute_ebic(gamma_ebic, n,
                                                        S = out.pfggm.rho.mle[[i]][["S"]],
                                                        Theta = out.pfggm.rho.mle[[i]][["Tht.hat"]]),
                          numeric(1))

# Costruiamo il data frame
df_ebic <- data.frame(
  rho = grid.rho,
  eBIC = ebic_vector.mle
)

# Calcolo della posizione della linea verticale
rho_min <- grid.rho[which.min(bic_vector.mle)]

# Grafico ggplot
p2 <- ggplot(df_ebic, aes(x = rho, y = eBIC)) +
  geom_line(linewidth = 0.8, color = "gray") +
  geom_point(size = 2, color = "gray") +
  geom_vline(xintercept = rho_min, linetype = "dashed", color = "gray20", linewidth = 1) +
  labs(x = expression(gamma[1]), y = "eBIC") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )
p2

id.opt <- which.min(bic_vector.mle)
THT.graph <- 1*(Reduce(`+`, array2list(out.pfggm.rho.mle[[id.opt]]$Tht.hat)) != 0)
dimnames(THT.graph) <- list(var_names, var_names)

par(mfrow = c(1, 1))
BDgraph::plot.graph(THT.graph, main = "",
                    vertex.size = 5, vertex.color = "black", vertex.frame.color = "black", vertex.label.dist = 1.25,
                    vertex.label.color = "black", edge.color = "gray")

df_long2 <- map_dfr(1:dim(out.pfggm.rho.mle[[id.opt]]$Ximputed)[3], function(j) {
  as.data.frame(compute_Ximputed(Phi.sample, out.pfggm.rho.mle[[id.opt]]$Xi)[,,j]) |>
    mutate(Time = tp, Variable = var_names[j])
})

df_tidy2 <- df_long2 |>
  pivot_longer(
    cols = -c(Time, Variable),
    names_to = "Sim",
    values_to = "Value"
  ) |>
  mutate(Highlight = ifelse(Sim == "V1", "Main", "Background"),
         Variable = factor(Variable, levels = var_names),
         isNA = factor(1*is.na(df_tidy$Value)))

p.2 <- ggplot() +
  # Linee grigie di background (prima, con alpha)
  geom_line(
    data = df_tidy2 %>% filter(Highlight == "Background"),
    aes(x = Time, y = Value, group = Sim),
    color = "grey80", linewidth = 0.4, alpha = 0.5
  ) +
  # Linee rosse in evidenza (poi)
  geom_line(
    data = df_tidy2 %>% filter(Highlight == "Main"),
    aes(x = Time, y = Value, group = Sim, col = isNA),
    linewidth = 1.2
  ) +
  scale_color_manual(values = c("0" = "gray60", "1" = "gray20")) +
  facet_wrap(~Variable, scales = "free_y") +
  labs(x = "Time", y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold")
  )
p.2                                    

# Estrae la stima Tht.hat
Tht <- out.pfggm.rho.mle[[id.opt]]$Tht.hat
Sgm <- out.pfggm.rho.mle[[id.opt]]$Sgm.hat

# Calcolo vettorizzato delle tilde sigma
denom_WI <- Tht[3,3,] * Tht[1,1,] - Tht[3,1,]^2
Tilde.sigma_Wind_Price <- - Tht[3,1,] / denom_WI
Sgm[3, 1, ]

denom_WI <- Tht[3,3,] * Tht[4,4,] - Tht[3,4,]^2
Tilde.sigma_Wind_Import <- - Tht[3,4,] / denom_WI
Sgm[3, 4, ]

denom_DI <- Tht[2,2,] * Tht[4,4,] - Tht[2,4,]^2
Tilde.sigma_Demand_Import <- - Tht[2,4,] / denom_DI
Sgm[2, 4, ]

# Funzione per costruire la matrice C in modo vettorizzato
make_C_matrix <- function(sigma_vec, Phi) {
  # moltiplica ogni colonna di Phi per il corrispondente sigma
  weighted_Phi <- sweep(Phi, 2, sigma_vec, FUN = "*")
  # calcola Phi %*% diag(sigma) %*% t(Phi) in un solo passaggio
  C <- Phi %*% t(weighted_Phi)
  return(C)
}

# Calcolo delle due matrici C
C_Wind_Price    <- make_C_matrix(Sgm[3, 1, ], Phi.sample)
C_Wind_Import   <- make_C_matrix(Sgm[3, 4, ], Phi.sample)
C_Demand_Import <- make_C_matrix(Sgm[2, 4, ], Phi.sample)

# Definisci i tick ogni 24 punti e le etichette corrispondenti
tick_positions <- seq(24, 144, by = 24)
day_labels <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat")

mesh <- expand.grid(x = 1:d, y = 1:d)
mesh$z <- c(C_Wind_Price)
p4 <- ggplot(mesh, aes(x = x, y = y, fill = z)) +
  geom_raster(interpolate = FALSE) +
  scale_fill_gradient2(
    low = "blue", mid = "green", high = "red",
    midpoint = 0, limits = c(min(min(C_Wind_Price),min(C_Wind_Import)),
                             max(max(C_Wind_Price),max(C_Wind_Import))), 
    name = NULL
  ) +
  coord_fixed() +  # mantiene le celle quadrate
  scale_x_continuous(
    breaks = tick_positions,
    labels = day_labels,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = tick_positions,
    labels = day_labels,
    expand = c(0, 0)
  ) +
  labs(
    x = "Wind power",
    y = "Electricity price"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    axis.ticks = element_line(color = "grey40"),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    legend.key.height = unit(0.6, "cm"),  # aumenta la lunghezza verticale
    legend.key.width  = unit(2.5, "cm"),  # aumenta lo spessore orizzontale
    plot.margin = margin(5, 10, 5, 5)
  ) +
  guides(
    fill = guide_colorbar(
      barwidth = 20,   # più largo orizzontalmente
      barheight = 1.2,   # più lungo verticalmente
      ticks.colour = "black",
      frame.colour = "black"
    )
  )
p4  

mesh2 <- expand.grid(x = 1:d, y = 1:d)
mesh2$z <- c(C_Wind_Import)
p5 <- ggplot(mesh2, aes(x = x, y = y, fill = z)) +
  geom_raster(interpolate = FALSE) +
  scale_fill_gradient2(
    low = "blue", mid = "green", high = "red",
    midpoint = 0, limits = c(min(min(C_Wind_Price),min(C_Wind_Import)),
                             max(max(C_Wind_Price),max(C_Wind_Import))), 
    name = NULL
  ) +
  coord_fixed() +  # mantiene le celle quadrate
  scale_x_continuous(
    breaks = tick_positions,
    labels = day_labels,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = tick_positions,
    labels = day_labels,
    expand = c(0, 0)
  ) +
  labs(
    x = "Wind power",
    y = "Electricity imports"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    axis.ticks = element_line(color = "grey40"),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    legend.key.height = unit(0.6, "cm"),  # aumenta la lunghezza verticale
    legend.key.width  = unit(2.5, "cm"),  # aumenta lo spessore orizzontale
    plot.margin = margin(5, 10, 5, 5)
  ) +
  guides(
    fill = guide_colorbar(
      barwidth = 20,   # più largo orizzontalmente
      barheight = 1.2,   # più lungo verticalmente
      ticks.colour = "black",
      frame.colour = "black"
    )
  )
p5

p6 <- (p4 + p5 + plot_layout(guides = "collect")) &
  theme(legend.position = "bottom",
        legend.margin = margin(t = -5, unit = "pt"))
p6

mean(C_Demand_Import)
