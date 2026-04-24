# ============================================================
# Sovereign ESG application pipeline
# Partially observed functional graphical model (PFGGM)
# ============================================================
# This script prepares the World Bank Sovereign ESG data,
# computes missingness summaries, estimates the empirical basis,
# fits the penalized graph path, selects the optimal model,
# and generates the figures used in the paper.
#
# The script is organized into sections so it can be run end-to-end
# or section by section.
#
# Required external helper functions (sourced below):
# - estimate_empirical_basis()
# - integrate_cube()
# - pofggm()
# - compute_Ximputed()
#
# ============================================================

# ============================================================
# 0. Setup
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(here)
  library(patchwork)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(BDgraph)
  library(grid)
})

set.seed(123)
options(stringsAsFactors = FALSE)

# ---- Paths ----
RAW_XLSX <- here("data", "raw", "esgdata_download-2026-01-09.xlsx")
HELPER_R <- "../01 - Code/helper.R"

DERIVED_DIR <- here("data", "derived")
FIGURES_DIR <- here("output", "figures")
TABLES_DIR  <- here("output", "tables")

for (path in c(DERIVED_DIR, FIGURES_DIR, TABLES_DIR)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

if (file.exists(HELPER_R)) {
  source(HELPER_R)
} else {
  message("helper.R not found at: ", HELPER_R)
  message("Functions such as estimate_empirical_basis() and pofggm() must be available in the environment.")
}

# ============================================================
# 1. Configuration
# ============================================================

years_keep <- as.character(1990:2020)
tp <- 1990:2020
sheet_id <- 4

# Core ESG framework indicators (78)
esg_indicators <- c(
  # Environment
  # Emissions & pollution
  "EN.ATM.PM25.MC.M3", "EN.GHG.CO2.MT.CE.AR5", "EN.GHG.CO2.PC.CE.AR5", "EN.GHG.CO2.LU.MT.CE.AR5",
  "EN.GHG.CH4.MT.CE.AR5", "EN.GHG.N2O.MT.CE.AR5", "EN.GHG.ALL.MT.CE.AR5", "EN.GHG.ALL.PC.CE.AR5",
  
  # Energy use & security
  "EG.ELC.COAL.ZS", "EG.IMP.CONS.ZS", "EG.EGY.PRIM.PP.KD", "EG.USE.PCAP.KG.OE", "EG.USE.COMM.FO.ZS",
  "EG.ELC.RNEW.ZS", "EG.FEC.RNEW.ZS",
  
  # Climate risk & resilience
  "EN.CLC.CSTP.ZS", "EN.CLC.CDDY.XD", "EN.CLC.HEAT.XD", "EN.CLC.HDDY.XD", "EN.LND.LTMP.DC",
  "ER.H2O.FWST.ZS", "EN.POP.DNST", "EN.H2O.BDYS.ZS", "EN.CLC.SPEI.XD",
  
  # Food Security
  "AG.LND.AGRI.ZS", "NV.AGR.TOTL.ZS", "AG.PRD.FOOD.XD",
  
  # Natural capital endowment & management
  "NY.ADJ.DRES.GN.ZS", "NY.ADJ.DFOR.GN.ZS", "ER.H2O.FWTL.ZS", "AG.LND.FRST.ZS", "EN.MAM.THRD.NO",
  "AG.LND.PFLS.HA", "ER.PTD.TOTL.ZS", "AG.LND.FRLS.HA",
  
  # Social
  # Access to Services
  "EG.CFT.ACCS.ZS", "EG.ELC.ACCS.ZS", "SH.H2O.SMDW.ZS", "SH.STA.SMSS.ZS",
  
  # Demography
  "SP.DYN.TFRT.IN", "SP.DYN.LE00.IN", "SP.POP.65UP.TO.ZS",
  
  # Education & skills
  "SE.XPD.TOTL.GB.ZS", "SE.ADT.LITR.ZS", "SE.PRM.ENRR",
  
  # Employment
  "SL.TLF.0714.ZS", "SL.TLF.ACTI.ZS", "SL.UEM.NEET.ME.ZS", "SL.UEM.TOTL.ZS", "SL.EMP.WORK.ZS",
  
  # Health & Nutrition
  "SH.DTH.COMM.ZS", "SH.MED.BEDS.ZS", "SH.DYN.MORT", "SH.STA.OWAD.ZS", "SN.ITK.DEFC.ZS",
  
  # Poverty & Inequality
  "SI.POV.GINI", "SI.DST.FRST.20", "SI.POV.DDAY", "SI.POV.UMIC", "SI.POV.NAHC", "SI.SPR.PGAP",
  
  # Governance
  # Economic Environment
  "NY.GDP.MKTP.KD.ZG", "IT.NET.USER.ZS",
  
  # Gender
  "SG.GEN.PARL.ZS", "SL.TLF.CACT.FM.ZS", "SE.ENR.PRSC.FM.ZS", "SP.UWT.TFRT",
  
  # Government Effectiveness
  "GE.EST", "RQ.EST",
  
  # Human Rights
  "SD.ESR.PERF.XQ", "VA.EST",
  
  # Innovation
  "IP.PAT.RESD", "GB.XPD.RSDV.GD.ZS", "IP.JRN.ARTC.SC",
  
  # Stability & Rule of Law
  "CC.EST", "SM.POP.NETM", "PV.EST", "RL.EST"
)

# Indicator-specific preprocessing
log1p_vars <- c(
  "IP.PAT.RESD",
  "ER.H2O.FWTL.ZS",
  "ER.H2O.FWST.ZS",
  "EG.USE.PCAP.KG.OE",
  "NY.ADJ.DRES.GN.ZS",
  "SI.SPR.PGAP",
  "AG.PRD.FOOD.XD",
  "AG.LND.FRLS.HA",
  "SH.DYN.MORT",
  "IP.JRN.ARTC.SC",
  "EG.EGY.PRIM.PP.KD",
  "SH.MED.BEDS.ZS",
  "EN.ATM.PM25.MC.M3",
  "EN.GHG.CO2.MT.CE.AR5",
  "EN.GHG.CO2.PC.CE.AR5",
  "EN.GHG.CH4.MT.CE.AR5",
  "EN.GHG.N2O.MT.CE.AR5",
  "EN.GHG.ALL.MT.CE.AR5",
  "EN.GHG.ALL.PC.CE.AR5",
  "EN.MAM.THRD.NO",
  "EN.POP.DNST"
)

# Selected indicator display names used in figures
short_names <- c(
  "EN.GHG.CO2.MT.CE.AR5" = "CO2 total",
  "EN.GHG.CO2.PC.CE.AR5" = "CO2 pc",
  "EN.GHG.CH4.MT.CE.AR5" = "CH4 total",
  "EN.GHG.N2O.MT.CE.AR5" = "N2O total",
  "EN.GHG.ALL.MT.CE.AR5" = "GHG total",
  "EN.GHG.ALL.PC.CE.AR5" = "GHG pc",
  "EG.USE.PCAP.KG.OE"    = "Energy pc",
  "EG.FEC.RNEW.ZS"       = "Renewables",
  "EN.CLC.CDDY.XD"       = "CDD",
  "EN.CLC.HEAT.XD"       = "Heat",
  "EN.CLC.HDDY.XD"       = "HDD",
  "ER.H2O.FWST.ZS"       = "Water stress",
  "ER.H2O.FWTL.ZS"       = "Water withdr.",
  "AG.LND.AGRI.ZS"       = "Agri land",
  "AG.LND.FRST.ZS"       = "Forest area",
  "NV.AGR.TOTL.ZS"       = "Agri VA",
  "SP.DYN.TFRT.IN"       = "Fertility",
  "SP.DYN.LE00.IN"       = "Life exp.",
  "SP.POP.65UP.TO.ZS"    = "Age 65+",
  "SH.DYN.MORT"          = "U5 mortality",
  "SH.STA.OWAD.ZS"       = "Overweight"
)

indicator_names_missing <- c(
  "EG.ELC.COAL.ZS"       = "Coal electricity",
  "EG.USE.PCAP.KG.OE"    = "Energy use per capita",
  "EG.ELC.RNEW.ZS"       = "Renewable electricity",
  "EG.FEC.RNEW.ZS"       = "Renewable energy use",
  "ER.H2O.FWST.ZS"       = "Water stress",
  "EN.POP.DNST"          = "Population density",
  "AG.LND.AGRI.ZS"       = "Agricultural land",
  "NV.AGR.TOTL.ZS"       = "Agriculture value added",
  "AG.PRD.FOOD.XD"       = "Food production",
  "NY.ADJ.DRES.GN.ZS"    = "Natural resource depletion",
  "NY.ADJ.DFOR.GN.ZS"    = "Forest depletion",
  "ER.H2O.FWTL.ZS"       = "Freshwater withdrawals",
  "AG.LND.FRST.ZS"       = "Forest area",
  "NY.GDP.MKTP.KD.ZG"    = "GDP growth"
)

# Graph coloring
node_pillar <- c(
  "EN.GHG.CO2.MT.CE.AR5" = "Environment",
  "EN.GHG.CO2.PC.CE.AR5" = "Environment",
  "EN.GHG.CH4.MT.CE.AR5" = "Environment",
  "EN.GHG.N2O.MT.CE.AR5" = "Environment",
  "EN.GHG.ALL.MT.CE.AR5" = "Environment",
  "EN.GHG.ALL.PC.CE.AR5" = "Environment",
  "EG.USE.PCAP.KG.OE"    = "Environment",
  "EG.FEC.RNEW.ZS"       = "Environment",
  "EN.CLC.CDDY.XD"       = "Environment",
  "EN.CLC.HEAT.XD"       = "Environment",
  "EN.CLC.HDDY.XD"       = "Environment",
  "ER.H2O.FWST.ZS"       = "Environment",
  "ER.H2O.FWTL.ZS"       = "Environment",
  "AG.LND.AGRI.ZS"       = "Environment",
  "AG.LND.FRST.ZS"       = "Environment",
  "NV.AGR.TOTL.ZS"       = "Environment",
  "SP.DYN.TFRT.IN"       = "Social",
  "SP.DYN.LE00.IN"       = "Social",
  "SP.POP.65UP.TO.ZS"    = "Social",
  "SH.DYN.MORT"          = "Social",
  "SH.STA.OWAD.ZS"       = "Social"
)

node_group <- c(
  "EN.GHG.CO2.MT.CE.AR5" = "Emissions & pollution",
  "EN.GHG.CO2.PC.CE.AR5" = "Emissions & pollution",
  "EN.GHG.CH4.MT.CE.AR5" = "Emissions & pollution",
  "EN.GHG.N2O.MT.CE.AR5" = "Emissions & pollution",
  "EN.GHG.ALL.MT.CE.AR5" = "Emissions & pollution",
  "EN.GHG.ALL.PC.CE.AR5" = "Emissions & pollution",
  "EG.USE.PCAP.KG.OE"    = "Energy use & security",
  "EG.FEC.RNEW.ZS"       = "Energy use & security",
  "EN.CLC.CDDY.XD"       = "Climate risk & resilience",
  "EN.CLC.HEAT.XD"       = "Climate risk & resilience",
  "EN.CLC.HDDY.XD"       = "Climate risk & resilience",
  "ER.H2O.FWST.ZS"       = "Climate risk & resilience",
  "ER.H2O.FWTL.ZS"       = "Natural capital use & management",
  "AG.LND.AGRI.ZS"       = "Food security",
  "NV.AGR.TOTL.ZS"       = "Food security",
  "AG.LND.FRST.ZS"       = "Natural capital use & management",
  "SP.DYN.TFRT.IN"       = "Demography",
  "SP.DYN.LE00.IN"       = "Demography",
  "SP.POP.65UP.TO.ZS"    = "Demography",
  "SH.DYN.MORT"          = "Health & nutrition",
  "SH.STA.OWAD.ZS"       = "Health & nutrition"
)

group_cols <- c(
  "Emissions & pollution"            = "#8E1B1B",
  "Energy use & security"            = "#C73E1D",
  "Climate risk & resilience"        = "#E76F51",
  "Food security"                    = "#F4A261",
  "Natural capital use & management" = "#F6D6AD",
  "Demography"                       = "#2166ac",
  "Employment"                       = "#4393c3",
  "Health & nutrition"               = "#92c5de",
  "Economic environment"             = "#1b7837",
  "Gender"                           = "#5aae61",
  "Stability & rule of law"          = "#a6dba0"
)

pillar_cols <- c(
  "Environment" = "#d73027",
  "Social"      = "#4575b4",
  "Governance"  = "#1a9850"
)

# ============================================================
# 2. Utility functions
# ============================================================

build_data_array <- function(data_long_norm, years_keep) {
  data_list <- split(data_long_norm, data_long_norm$`Indicator code`)
  countries <- sort(unique(data_long_norm$Economy))
  d <- length(years_keep)
  p <- length(data_list)
  
  arr <- array(
    NA_real_,
    dim = c(length(countries), d, p),
    dimnames = list(countries, years_keep, names(data_list))
  )
  
  for (j in seq_len(p)) {
    arr[, , j] <- data_list[[j]] %>%
      dplyr::select(Economy, year, value_norm) %>%
      pivot_wider(names_from = year, values_from = value_norm) %>%
      arrange(Economy) %>%
      dplyr::select(-Economy) %>%
      as.matrix()
  }
  
  arr
}

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

W_from_precision <- function(Theta) {
  p <- nrow(Theta)
  W <- matrix(0, p, p)
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      if (i != j) {
        W[i, j] <- -Theta[i, j] / (Theta[i, i] * Theta[j, j] - Theta[i, j]^2)
      }
    }
  }
  W
}

PHI_smooth <- function(tp, phi, n = 100) {
  tp_fine <- seq(tp[1], tp[length(tp)], length.out = n)
  phi_smooth_func <- smooth.spline(x = tp, y = phi)
  phi_fine <- predict(phi_smooth_func, x = tp_fine)$y
  tcrossprod(phi_fine)
}

# ============================================================
# 3. Read raw data and build long panel
# ============================================================

raw_data <- read_excel(RAW_XLSX, sheet = sheet_id)

esg_data <- raw_data %>%
  filter(`Indicator code` %in% esg_indicators) %>%
  mutate(`Indicator code` = factor(`Indicator code`, levels = esg_indicators)) %>%
  distinct() %>%
  dplyr::select(`ISO3 code`, Economy, `Indicator code`, `Indicator name`, all_of(years_keep))

esg_long <- esg_data %>%
  dplyr::select(`ISO3 code`, Economy, `Indicator code`, `Indicator name`, all_of(years_keep)) %>%
  pivot_longer(
    cols = all_of(years_keep),
    names_to = "year",
    values_to = "value"
  ) %>%
  mutate(
    year = as.integer(year),
    value = as.numeric(value)
  )

saveRDS(esg_long, file.path(DERIVED_DIR, "esg_long_raw.rds"))

# ============================================================
# 4. Indicator and country screening
# ============================================================

indicator_counts <- esg_data %>%
  count(`Indicator code`, Economy) %>%
  count(`Indicator code`, name = "n_countries")

to_keep <- indicator_counts %>%
  filter(n_countries >= 170) %>%
  pull(`Indicator code`) %>%
  as.character()

esg_data_1 <- esg_data %>% filter(`Indicator code` %in% to_keep)
esg_long_1 <- esg_long %>% filter(`Indicator code` %in% to_keep)

country_counts <- esg_data_1 %>%
  count(Economy, `Indicator code`) %>%
  count(Economy, name = "n_indicators")

to_keep_countries <- country_counts %>%
  filter(n_indicators >= length(to_keep)) %>%
  pull(Economy)

esg_data_2 <- esg_data_1 %>%
  filter(Economy %in% to_keep_countries) %>%
  mutate(`Indicator code` = factor(`Indicator code`, levels = to_keep))

esg_long_2 <- esg_long_1 %>%
  filter(Economy %in% to_keep_countries) %>%
  mutate(`Indicator code` = factor(`Indicator code`, levels = to_keep))

# ============================================================
# 5. Indicator-specific transformation and standardization
# ============================================================

esg_long_norm <- esg_long_2 %>%
  mutate(
    value_tr = case_when(
      `Indicator code` %in% log1p_vars ~ log1p(value),
      TRUE ~ value
    )
  ) %>%
  group_by(`Indicator code`, year) %>%
  mutate(
    mu = mean(value_tr, na.rm = TRUE),
    sdv = sd(value_tr, na.rm = TRUE),
    value_norm = (value_tr - mu) / ifelse(is.na(sdv) | sdv == 0, 1, sdv)
  ) %>%
  ungroup()

saveRDS(esg_long_norm, file.path(DERIVED_DIR, "esg_long_norm.rds"))

# ============================================================
# 6. Build array and screen by missingness
# ============================================================

esg_array <- build_data_array(esg_long_norm, years_keep)

missing_values <- sapply(seq_len(dim(esg_array)[3]), function(j) {
  is.na(esg_array[, , j]) %>% colMeans()
})

dimnames(missing_values) <- list(tp, levels(esg_data_2$`Indicator code`))
missing_values <- t(missing_values)

code_filter <- which(!apply(missing_values, 1, function(r) any(r >= 0.5)))

# Manual exclusion used in the paper
code_filter <- code_filter[-10]

esg_array_filter <- esg_array[, , code_filter, drop = FALSE]
esg_array_filter <- aperm(esg_array_filter, c(2, 1, 3))

saveRDS(esg_array_filter, file.path(DERIVED_DIR, "esg_array_filter.rds"))

# ============================================================
# 7. Missingness summaries and plots
# ============================================================

var_names <- dimnames(esg_array_filter)[[3]]

tp_filter <- seq_along(tp)

pi_po <- setNames(
  sapply(seq_len(dim(esg_array_filter)[3]), function(i) {
    mean(colMeans(is.na(esg_array_filter[, , i])) > 0)
  }),
  nm = var_names
)

pi_w <- sapply(which(pi_po > 0), function(i) {
  temp <- colMeans(is.na(esg_array_filter[, , i]))
  c(
    pi_w = mean(temp[temp > 0]),
    pi_w_min = min(temp[temp > 0]),
    pi_w_max = max(temp[temp > 0])
  )
})

missingness_table <- cbind(`pi_po` = pi_po[pi_po > 0], t(pi_w)) %>% round(3)
write.csv(missingness_table, file.path(TABLES_DIR, "missingness_table.csv"))
saveRDS(list(pi_po = pi_po, pi_w = pi_w), file.path(DERIVED_DIR, "missingness_stats.rds"))

# Histogram of pi_po
df_pi_po <- tibble(
  indicator = names(pi_po),
  pi_po = as.numeric(pi_po)
)

p_pi_po <- ggplot(df_pi_po, aes(x = pi_po)) +
  geom_histogram(bins = 20, color = "black", fill = "gray70") +
  labs(x = expression(pi[po]), y = "Count") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

# Dot plot of pi_w
# pi_w has rows pi_w, pi_w_min, pi_w_max and columns indicator codes
df_pi_w <- data.frame(
  `Indicator code` = colnames(pi_w),
  pi_w = as.numeric(pi_w[1, ]),
  Indicator = unname(indicator_names_missing[colnames(pi_w)])
) %>%
  mutate(Indicator = ifelse(is.na(Indicator), `Indicator code`, Indicator)) %>%
  arrange(pi_w) %>%
  mutate(Indicator = factor(Indicator, levels = Indicator))

p_pi_w <- ggplot(df_pi_w, aes(x = pi_w, y = Indicator)) +
  geom_segment(aes(x = 0, xend = pi_w, y = Indicator, yend = Indicator), color = "gray70") +
  geom_point(size = 2.8, color = "black") +
  labs(x = expression(pi[w]), y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

combined_missingness <- p_pi_po + p_pi_w + plot_layout(ncol = 2)

ggsave(
  filename = file.path(FIGURES_DIR, "missingness_summary.pdf"),
  plot = combined_missingness,
  width = 12,
  height = 5,
  device = "pdf"
)

# ============================================================
# 8. Estimate empirical basis
# ============================================================

basis_fit <- try(
  estimate_empirical_basis(
    X = esg_array_filter,
    tp = tp[tp_filter],
    pev = 0.9999,
    smooth_k = 15
  ),
  silent = TRUE
)

if (inherits(basis_fit, "try-error")) {
  stop("estimate_empirical_basis() failed.")
}

saveRDS(basis_fit, file.path(DERIVED_DIR, "basis_fit.rds"))

Phi_emp <- basis_fit$Phi_emp
XiEst <- integrate_cube(esg_array_filter, Phi_emp, tp[tp_filter])

# ============================================================
# 9. Fit penalized path
# ============================================================

id_pobs <- which(apply(apply(is.na(esg_array_filter), 3, colSums) > 0, 1, any))
id_obs <- setdiff(seq_len(dim(esg_array_filter)[2]), id_pobs)

fit0 <- pofggm(
  id_pobs = id_pobs,
  id_obs = id_obs,
  X = esg_array_filter,
  Phi = Phi_emp,
  tp = tp[tp_filter],
  Sgm.hat = NULL,
  Tht.hat = NULL,
  wTht = NULL,
  maxit.admm = 1e5,
  rho = .Machine$double.xmax,
  gamma = 0.0,
  alpha = rep(0.0, length(id_pobs)),
  ncores = 4L,
  verbose = TRUE,
  thr.em = 1e-5,
  thr.admm = 1e-6
)

gamma1.vec <- fit0$rho.max * exp(seq(log(1), log(0.01), length.out = 21))

fitList <- vector("list", length = length(gamma1.vec))
for (i in seq_along(gamma1.vec)) {
  fitList[[i]] <- pofggm(
    id_pobs = id_pobs,
    id_obs = id_obs,
    X = esg_array_filter,
    Phi = Phi_emp,
    tp = tp[tp_filter],
    Sgm.hat = NULL,
    Tht.hat = NULL,
    wTht = NULL,
    maxit.admm = 1e5,
    rho = gamma1.vec[i],
    gamma = 0.0,
    alpha = fit0$alpha_opt,
    verbose = TRUE,
    thr.em = 1e-5,
    thr.admm = 1e-6
  )
}

saveRDS(list(fit0 = fit0, gamma1.vec = gamma1.vec, fitList = fitList),
        file.path(DERIVED_DIR, "fit_path.rds"))

# ============================================================
# 10. Refit and model selection
# ============================================================

out.pfggm.rho.mle <- vector("list", length = length(fitList))
for (kk in seq_along(fitList)) {
  out.pfggm.rho <- fitList[[kk]]
  out.pfggm.rho.mle[[kk]] <- pofggm(
    id_pobs = id_pobs,
    id_obs = id_obs,
    X = esg_array_filter,
    Phi = Phi_emp,
    tp = tp[tp_filter],
    Sgm.hat = out.pfggm.rho$Sgm.hat,
    Tht.hat = out.pfggm.rho$Tht.hat,
    wTht = (out.pfggm.rho$Tht.hat == 0) * .Machine$double.xmax,
    alpha = out.pfggm.rho$alpha_opt,
    rho = 1e-12,
    ncores = 4L,
    verbose = TRUE
  )
}

idx <- seq_along(fitList)
gamma_ebic <- 0.5
n_countries <- dim(esg_array_filter)[2]

aic_vector <- vapply(idx, function(i) compute_aic(n_countries, S = fitList[[i]]$S, Theta = fitList[[i]]$Tht.hat), numeric(1))
aic_vector.mle <- vapply(idx, function(i) compute_aic(n_countries, S = out.pfggm.rho.mle[[i]]$S, Theta = out.pfggm.rho.mle[[i]]$Tht.hat), numeric(1))

bic_vector <- vapply(idx, function(i) compute_bic(n_countries, S = fitList[[i]]$S, Theta = fitList[[i]]$Tht.hat), numeric(1))
bic_vector.mle <- vapply(idx, function(i) compute_bic(n_countries, S = out.pfggm.rho.mle[[i]]$S, Theta = out.pfggm.rho.mle[[i]]$Tht.hat), numeric(1))

ebic_vector <- vapply(idx, function(i) compute_ebic(gamma_ebic, n_countries, S = fitList[[i]]$S, Theta = fitList[[i]]$Tht.hat), numeric(1))
ebic_vector.mle <- vapply(idx, function(i) compute_ebic(gamma_ebic, n_countries, S = out.pfggm.rho.mle[[i]]$S, Theta = out.pfggm.rho.mle[[i]]$Tht.hat), numeric(1))

df_ic <- data.frame(
  gamma1 = gamma1.vec,
  AIC = aic_vector,
  AIC_MLE = aic_vector.mle,
  BIC = bic_vector,
  BIC_MLE = bic_vector.mle,
  eBIC = ebic_vector,
  eBIC_MLE = ebic_vector.mle
)

saveRDS(list(refit = out.pfggm.rho.mle, ic_table = df_ic), file.path(DERIVED_DIR, "fit_refit.rds"))

id.opt <- which.min(ebic_vector.mle)
selected_fit <- out.pfggm.rho.mle[[id.opt]]
saveRDS(selected_fit, file.path(DERIVED_DIR, "selected_model.rds"))

# ============================================================
# 11. eBIC plot and network figure
# ============================================================

gamma1_min <- gamma1.vec[id.opt]

p_ebic <- ggplot(df_ic, aes(x = gamma1, y = eBIC_MLE)) +
  geom_line(linewidth = 0.8, color = "gray") +
  geom_point(size = 2, color = "gray") +
  geom_vline(xintercept = gamma1_min, linetype = "dashed", color = "gray20", linewidth = 1) +
  labs(x = expression(gamma[1]), y = "eBIC") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

THT.graph <- 1 * (Reduce(`+`, array2list(selected_fit$Tht.hat)) != 0)
dimnames(THT.graph) <- list(var_names, var_names)

graph <- BDgraph::get_graph(THT.graph, cut = 0.5)
graph_ig <- graph_from_adjacency_matrix(graph, mode = "undirected", diag = FALSE)
g2 <- delete_vertices(graph_ig, V(graph_ig)[degree(graph_ig) == 0])

V(g2)$label <- short_names[V(g2)$name]
V(g2)$pillar <- node_pillar[V(g2)$name]
V(g2)$group <- node_group[V(g2)$name]
V(g2)$deg <- degree(g2)

g_tbl <- as_tbl_graph(g2) %>%
  activate(nodes) %>%
  mutate(
    label = short_names[name],
    pillar = factor(node_pillar[name], levels = c("Environment", "Social", "Governance")),
    group = factor(
      node_group[name],
      levels = c(
        "Emissions & pollution", "Energy use & security", "Climate risk & resilience",
        "Food security", "Natural capital use & management",
        "Demography", "Employment", "Health & nutrition",
        "Economic environment", "Gender", "Stability & rule of law"
      )
    ),
    deg = centrality_degree()
  )

p_net <- ggraph(g_tbl, layout = "kk") +
  geom_edge_link(color = "gray70", width = 0.8) +
  geom_node_point(
    aes(size = deg, fill = group, color = pillar),
    shape = 21, stroke = 1.2
  ) +
  geom_node_text(
    aes(label = label),
    repel = TRUE,
    size = 4,
    force = 3,
    box.padding = 0.6,
    point.padding = 0.4,
    max.overlaps = Inf
  ) +
  scale_size(range = c(4, 10), guide = "none") +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = pillar_cols) +
  guides(
    fill = guide_legend(
      title = "Group",
      override.aes = list(size = 6, shape = 21, color = "black")
    ),
    color = guide_legend(
      title = "Pillar",
      override.aes = list(size = 6, shape = 21, fill = "white", stroke = 1.5)
    )
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12)
  )

ggsave(
  filename = file.path(FIGURES_DIR, "esg_network.pdf"),
  plot = p_net,
  width = 11,
  height = 8.5,
  device = "pdf"
)

combined_network <- p_ebic + p_net + plot_layout(
  guides = "collect",
  design = c(
    patchwork::area(1, 1, 1, 1),
    patchwork::area(1, 2, 1, 3)
  )
)

ggsave(
  filename = file.path(FIGURES_DIR, "ebic_esg_network.pdf"),
  plot = combined_network,
  width = 14,
  height = 6,
  device = "pdf"
)

# ============================================================
# 12. Cross-covariance surfaces
# ============================================================

edge_tbl <- as_data_frame(g2, what = "edges") %>%
  mutate(
    from_label = short_names[from],
    to_label = short_names[to]
  )

W_list <- lapply(array2list(selected_fit$Tht.hat), W_from_precision)

idx_3_26 <- match(c("EN.GHG.CO2.PC.CE.AR5", "SP.DYN.LE00.IN"), var_names)
idx_9_29 <- match(c("EG.USE.PCAP.KG.OE", "SH.DYN.MORT"), var_names)
idx_11_19 <- match(c("EG.FEC.RNEW.ZS", "NV.AGR.TOTL.ZS"), var_names)

W_3_26 <- sapply(W_list, function(W) W[idx_3_26[1], idx_3_26[2]])
W_9_29 <- sapply(W_list, function(W) W[idx_9_29[1], idx_9_29[2]])
W_11_19 <- sapply(W_list, function(W) W[idx_11_19[1], idx_11_19[2]])

SIGMA_3_26 <- Reduce("+", lapply(seq_len(length(W_3_26)), function(i) W_3_26[i] * PHI_smooth(tp[tp_filter], Phi_emp[, i], n = 301)))
SIGMA_9_29 <- Reduce("+", lapply(seq_len(length(W_9_29)), function(i) W_9_29[i] * PHI_smooth(tp[tp_filter], Phi_emp[, i], n = 301)))
SIGMA_11_19 <- Reduce("+", lapply(seq_len(length(W_11_19)), function(i) W_11_19[i] * PHI_smooth(tp[tp_filter], Phi_emp[, i], n = 301)))

cross_cov_summary <- tibble(
  pair = c("CO2pc-LifeExp", "Energy-U5Mort", "Renewables-AgriVA"),
  mean = c(mean(SIGMA_3_26), mean(SIGMA_9_29), mean(SIGMA_11_19)),
  min = c(min(SIGMA_3_26), min(SIGMA_9_29), min(SIGMA_11_19)),
  max = c(max(SIGMA_3_26), max(SIGMA_9_29), max(SIGMA_11_19))
)
write.csv(cross_cov_summary, file.path(TABLES_DIR, "cross_cov_summary.csv"), row.names = FALSE)

tp_fine <- seq(tp[1], tp[length(tp)], length.out = 301)
tick_positions <- seq(1, 301, by = 50)
year_labels <- as.character(round(tp_fine[tick_positions]))

make_crosscov_plot <- function(Sigma_mat, midpoint, xlab, ylab) {
  mesh <- expand.grid(x = seq_len(nrow(Sigma_mat)), y = seq_len(ncol(Sigma_mat)))
  mesh$z <- c(Sigma_mat)
  
  ggplot(mesh, aes(x = x, y = y, fill = z)) +
    geom_raster(interpolate = FALSE) +
    scale_fill_gradient2(
      low = "blue", mid = "green", high = "red",
      midpoint = midpoint,
      limits = range(Sigma_mat),
      name = NULL
    ) +
    coord_fixed() +
    scale_x_continuous(breaks = tick_positions, labels = year_labels, expand = c(0, 0)) +
    scale_y_continuous(breaks = tick_positions, labels = year_labels, expand = c(0, 0)) +
    labs(x = xlab, y = ylab) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13),
      axis.ticks = element_line(color = "grey40"),
      legend.text = element_text(size = 12),
      legend.position = "bottom",
      legend.key.height = unit(0.6, "cm"),
      legend.key.width = unit(2.5, "cm"),
      plot.margin = margin(5, 10, 5, 5)
    ) +
    guides(
      fill = guide_colorbar(
        barwidth = 20,
        barheight = 1.2,
        ticks.colour = "black",
        frame.colour = "black"
      )
    )
}

p4 <- make_crosscov_plot(
  Sigma_mat = SIGMA_3_26,
  midpoint = min(SIGMA_3_26),
  xlab = expression(CO[2] ~ emissions ~ per ~ capita),
  ylab = "Life expectancy at birth"
)

p5 <- make_crosscov_plot(
  Sigma_mat = SIGMA_9_29,
  midpoint = max(SIGMA_9_29),
  xlab = "Energy use per capita",
  ylab = "Under-5 mortality rate"
)

p6 <- make_crosscov_plot(
  Sigma_mat = SIGMA_11_19,
  midpoint = min(SIGMA_11_19),
  xlab = "Renewable energy consumption",
  ylab = "Agriculture, forestry and fishing value added"
)

combined_crosscov <- p4 + p5 + p6 + plot_layout(
  design = c(
    patchwork::area(1, 1, 1, 1),
    patchwork::area(1, 2, 1, 2),
    patchwork::area(1, 3, 1, 3)
  )
)

ggsave(
  filename = file.path(FIGURES_DIR, "cross_cov.pdf"),
  plot = combined_crosscov,
  width = 25,
  height = 10,
  device = "pdf",
  scale = 0.6,
  dpi = 300,
  units = "in"
)

# ============================================================
# 13. Example raw and reconstructed curves
# ============================================================

# Example country used in the paper plots
main_country <- "Equatorial Guinea"
plot_vars <- c("EG.ELC.COAL.ZS", "EG.ELC.RNEW.ZS", "NV.AGR.TOTL.ZS")
plot_var_labels <- c(
  "EG.ELC.COAL.ZS" = "Electricity production from coal sources",
  "EG.ELC.RNEW.ZS" = "Renewable electricity output",
  "NV.AGR.TOTL.ZS" = "Agriculture, forestry, and fishing, value added"
)

# Raw standardized curves
raw_long_plot <- purrr::map_dfr(seq_len(dim(esg_array_filter)[3]), function(j) {
  as.data.frame(esg_array_filter[, , j]) %>%
    mutate(Time = tp, Variable = var_names[j])
})

raw_tidy_plot <- raw_long_plot %>%
  pivot_longer(
    cols = -c(Time, Variable),
    names_to = "Country",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight = ifelse(Country == main_country, "Main", "Background"),
    Variable = factor(Variable, levels = var_names)
  )

p_raw <- ggplot() +
  geom_line(
    data = raw_tidy_plot %>% filter(Highlight == "Background", Variable %in% plot_vars),
    aes(x = Time, y = Value, group = Country),
    color = "gray70", linewidth = 0.4, alpha = 0.5
  ) +
  geom_line(
    data = raw_tidy_plot %>% filter(Highlight == "Main", Variable %in% plot_vars),
    aes(x = Time, y = Value, group = Country),
    color = "gray10", linewidth = 1.2
  ) +
  facet_grid(
    cols = vars(Variable),
    scales = "free_y",
    labeller = labeller(Variable = as_labeller(plot_var_labels))
  ) +
  labs(x = "Time", y = "Raw curves") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Reconstructed curves
Ximputed_array <- compute_Ximputed(Phi_emp, selected_fit$Xi)

recon_long_plot <- purrr::map_dfr(seq_len(dim(Ximputed_array)[3]), function(j) {
  temp_data <- Ximputed_array[, , j]
  colnames(temp_data) <- colnames(esg_array_filter[, , j])
  as.data.frame(temp_data) %>%
    mutate(Time = tp[tp_filter], Variable = var_names[j])
})

recon_tidy_plot <- recon_long_plot %>%
  pivot_longer(
    cols = -c(Time, Variable),
    names_to = "Country",
    values_to = "Value"
  ) %>%
  mutate(
    Highlight = ifelse(Country == main_country, "Main", "Background"),
    Variable = factor(Variable, levels = var_names)
  )

main_country_na <- raw_tidy_plot %>%
  filter(Country == main_country) %>%
  mutate(isNA = factor(1 * is.na(Value))) %>%
  dplyr::select(Time, Variable, Country, isNA)

recon_tidy_plot <- recon_tidy_plot %>%
  left_join(main_country_na, by = c("Time", "Variable", "Country"))

p_recon <- ggplot() +
  geom_line(
    data = recon_tidy_plot %>% filter(Highlight == "Background", Variable %in% plot_vars),
    aes(x = Time, y = Value, group = Country),
    color = "grey80", linewidth = 0.4, alpha = 0.5
  ) +
  geom_line(
    data = recon_tidy_plot %>% filter(Highlight == "Main", Variable %in% plot_vars),
    aes(x = Time, y = Value, group = Country, color = isNA),
    linewidth = 1.2
  ) +
  scale_color_manual(values = c("0" = "gray10", "1" = "gray50")) +
  facet_grid(
    cols = vars(Variable),
    scales = "free_y",
    labeller = labeller(Variable = as_labeller(plot_var_labels))
  ) +
  labs(x = "Time", y = "Reconstructed curves") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

combined_example <- (p_raw + theme(axis.text.x = element_blank(), axis.title.x = element_blank())) /
  (p_recon + theme(strip.text.x = element_blank(), strip.background.x = element_blank()))

ggsave(
  filename = file.path(FIGURES_DIR, "example_esg.pdf"),
  plot = combined_example,
  device = "pdf",
  scale = 0.7,
  width = 18,
  height = 11,
  units = "in",
  dpi = 300
)

# ============================================================
# 14. Session outputs
# ============================================================

message("Done. Outputs saved in:")
message("- ", DERIVED_DIR)
message("- ", FIGURES_DIR)
message("- ", TABLES_DIR)
