# ============================================================
# Simulation results pipeline
# ============================================================
# Purpose:
#   - collect simulation outputs saved as .RData files
#   - assemble a tidy metrics table
#   - generate summary tables and figures for the paper
#
# Main outputs:
#   - output/tables/*.csv
#   - output/figures/*.pdf
#
# Notes:
#   - This script replaces the original monolithic plotting script
#     with a structured, repository-friendly workflow.
#   - It assumes that files named
#       jcgs_simul_<scenario>_config<id>.RData
#     are available in RESULTS_DIR.
#   - Each .RData file is expected to contain an object named `storage`.
# ============================================================

# ============================================================
# 0. Setup
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
set.seed(123)

RESULTS_DIR <- "~/Downloads/JCGS SIMUL/"
FIGURES_DIR <- here("output", "figures")
TABLES_DIR  <- here("output", "tables")

for (path in c(FIGURES_DIR, TABLES_DIR)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# ============================================================
# 1. Simulation design
# ============================================================

make_config_grid <- function(
    scenario,
    n,
    p_over_n,
    K_true,
    perc_window,
    perc_obs_curves,
    pev,
    perc_theta_share,
    d,
    graph_type
) {
  expand.grid(
    scenario = scenario,
    n = n,
    p_over_n = p_over_n,
    K_true = K_true,
    perc_window = perc_window,
    perc_obs_curves = perc_obs_curves,
    pev = pev,
    perc_theta_share = perc_theta_share,
    d = d,
    graph_type = graph_type,
    stringsAsFactors = FALSE
  )
}

CONFIGS_list <- list(
  make_config_grid(
    scenario = "p_and_n",
    n = c(50L, 100L),
    p_over_n = c(0.2, 1.2),
    K_true = 5L,
    perc_window = 0.5,
    perc_obs_curves = 0.5,
    pev = 0.99,
    perc_theta_share = 1.0,
    d = 50L,
    graph_type = c("star", "band", "smallworld")
  ),
  make_config_grid(
    scenario = "missing",
    n = 100L,
    p_over_n = 0.2,
    K_true = 5L,
    perc_window = c(0.3, 0.5, 0.7),
    perc_obs_curves = c(0.3, 0.5, 0.7),
    pev = 0.99,
    perc_theta_share = 1.0,
    d = 50L,
    graph_type = c("star", "band", "smallworld")
  ),
  make_config_grid(
    scenario = "pev",
    n = 100L,
    p_over_n = 0.2,
    K_true = 5L,
    perc_window = 0.5,
    perc_obs_curves = 0.5,
    pev = c(0.90, 0.95, 0.99),
    perc_theta_share = 1.0,
    d = 50L,
    graph_type = c("star", "band", "smallworld")
  ),
  make_config_grid(
    scenario = "L",
    n = 100L,
    p_over_n = 0.2,
    K_true = c(5L, 9L),
    perc_window = 0.5,
    perc_obs_curves = 0.5,
    pev = 0.99,
    perc_theta_share = 1.0,
    d = 50L,
    graph_type = c("star", "band", "smallworld")
  ),
  make_config_grid(
    scenario = "part_sep",
    n = 100L,
    p_over_n = 0.2,
    K_true = 5L,
    perc_window = 0.5,
    perc_obs_curves = 0.5,
    pev = 0.99,
    perc_theta_share = 1.0,
    d = 50L,
    graph_type = c("star", "band", "smallworld")
  ),
  make_config_grid(
    scenario = "d",
    n = 100L,
    p_over_n = 0.2,
    K_true = 5L,
    perc_window = 0.5,
    perc_obs_curves = 0.5,
    pev = 0.99,
    perc_theta_share = 1.0,
    d = c(50L, 100L, 200L),
    graph_type = c("star", "band", "smallworld")
  ),
  make_config_grid(
    scenario = "gamma",
    n = 100L,
    p_over_n = 0.2,
    K_true = 5L,
    perc_window = 0.5,
    perc_obs_curves = 0.5,
    pev = 0.99,
    perc_theta_share = c(1.0, 0.5),
    d = 50L,
    graph_type = c("star", "band", "smallworld")
  ),
  make_config_grid(
    scenario = "alpha",
    n = 100L,
    p_over_n = 0.2,
    K_true = 5L,
    perc_window = 0.5,
    perc_obs_curves = 0.5,
    pev = 0.99,
    perc_theta_share = 1.0,
    d = 50L,
    graph_type = c("star", "band", "smallworld")
  )
)

# ============================================================
# 2. Helpers for loading and tidying metrics
# ============================================================

new_metric_row <- function(method, metric, val, cfg, gamma = "0.0", alpha = "best") {
  tibble(
    method = method,
    metric = metric,
    val = val,
    scenario = cfg$scenario,
    n = cfg$n,
    p_over_n = cfg$p_over_n,
    L = cfg$K_true,
    pi_po = cfg$perc_obs_curves,
    pi_d = cfg$perc_window,
    d = cfg$d,
    pev = cfg$pev,
    gamma = gamma,
    perc_theta_share = cfg$perc_theta_share,
    alpha = alpha,
    graph_type = cfg$graph_type
  )
}

metric_spec_default <- tribble(
  ~method,   ~metric,              ~storage_name,             ~aggregator, ~gamma, ~alpha,
  "poFGGM", "ThetaError",         "theta_err_mat",          "min_col",  "0.0", "best",
  "Kraus",  "ThetaError",         "theta_err_kraus_mat",    "min_col",  "0.0", "best",
  "Oracle", "ThetaError",         "theta_err_obs_mat",      "min_col",  "0.0", "none",
  "poFGGM", "AUC",                "auc_theta_vec",          "identity", "0.0", "best",
  "Kraus",  "AUC",                "auc_theta_kraus_vec",    "identity", "0.0", "best",
  "Oracle", "AUC",                "auc_theta_obs_vec",      "identity", "0.0", "none",
  "poFGGM", "CurveError",         "curve_err_mat",          "min_col",  "0.0", "best",
  "Kraus",  "CurveError",         "curve_err_kraus_vec",    "identity", "0.0", "best",
  "poFGGM", "ComputationalTime",  "comp_time_vec",          "identity", "0.0", "best"
)

metric_spec_gamma <- tribble(
  ~method,   ~metric,              ~storage_name,           ~aggregator, ~gamma, ~alpha,
  "poFGGM", "ThetaError",         "theta_err_mat",        "min_col",  "0.0", "best",
  "poFGGM", "ThetaError",         "theta_err_mat_2",      "min_col",  "0.5", "best",
  "poFGGM", "ThetaError",         "theta_err_mat_3",      "min_col",  "1.0", "best",
  "Kraus",  "ThetaError",         "theta_err_kraus_mat",  "min_col",  "0.0", "best",
  "Oracle", "ThetaError",         "theta_err_obs_mat",    "min_col",  "0.0", "none",
  "poFGGM", "AUC",                "auc_theta_vec",        "identity", "0.0", "best",
  "poFGGM", "AUC",                "auc_theta_vec_2",      "identity", "0.5", "best",
  "poFGGM", "AUC",                "auc_theta_vec_3",      "identity", "1.0", "best",
  "Kraus",  "AUC",                "auc_theta_kraus_vec",  "identity", "0.0", "best",
  "Oracle", "AUC",                "auc_theta_obs_vec",    "identity", "0.0", "none",
  "poFGGM", "CurveError",         "curve_err_mat",        "min_col",  "0.0", "best",
  "poFGGM", "CurveError",         "curve_err_mat_2",      "min_col",  "0.5", "best",
  "poFGGM", "CurveError",         "curve_err_mat_3",      "min_col",  "1.0", "best",
  "Kraus",  "CurveError",         "curve_err_kraus_vec",  "identity", "0.0", "best",
  "poFGGM", "ComputationalTime",  "comp_time_vec",        "identity", "0.0", "best",
  "poFGGM", "ComputationalTime",  "comp_time_vec_2",      "identity", "0.5", "best",
  "poFGGM", "ComputationalTime",  "comp_time_vec_3",      "identity", "1.0", "best"
)

metric_spec_alpha <- tribble(
  ~method,   ~metric,              ~storage_name,           ~aggregator, ~gamma, ~alpha,
  "poFGGM", "ThetaError",         "theta_err_mat",        "min_col",  "0.0", "best",
  "poFGGM", "ThetaError",         "theta_err_mat_2",      "min_col",  "0.0", "10^0*mean",
  "poFGGM", "ThetaError",         "theta_err_mat_3",      "min_col",  "0.0", "10^2*mean",
  "poFGGM", "ThetaError",         "theta_err_mat_4",      "min_col",  "0.0", "10^4*mean",
  "poFGGM", "ThetaError",         "theta_err_mat_5",      "min_col",  "0.0", "10^6*mean",
  "Kraus",  "ThetaError",         "theta_err_kraus_mat",  "min_col",  "0.0", "best",
  "Oracle", "ThetaError",         "theta_err_obs_mat",    "min_col",  "0.0", "none",
  "poFGGM", "AUC",                "auc_theta_vec",        "identity", "0.0", "best",
  "poFGGM", "AUC",                "auc_theta_vec_2",      "identity", "0.0", "10^0*mean",
  "poFGGM", "AUC",                "auc_theta_vec_3",      "identity", "0.0", "10^2*mean",
  "poFGGM", "AUC",                "auc_theta_vec_4",      "identity", "0.0", "10^4*mean",
  "poFGGM", "AUC",                "auc_theta_vec_5",      "identity", "0.0", "10^6*mean",
  "Kraus",  "AUC",                "auc_theta_kraus_vec",  "identity", "0.0", "best",
  "Oracle", "AUC",                "auc_theta_obs_vec",    "identity", "0.0", "none",
  "poFGGM", "CurveError",         "curve_err_mat",        "min_col",  "0.0", "best",
  "poFGGM", "CurveError",         "curve_err_mat_2",      "min_col",  "0.0", "10^0*mean",
  "poFGGM", "CurveError",         "curve_err_mat_3",      "min_col",  "0.0", "10^2*mean",
  "poFGGM", "CurveError",         "curve_err_mat_4",      "min_col",  "0.0", "10^4*mean",
  "poFGGM", "CurveError",         "curve_err_mat_5",      "min_col",  "0.0", "10^6*mean",
  "Kraus",  "CurveError",         "curve_err_kraus_vec",  "identity", "0.0", "best",
  "poFGGM", "ComputationalTime",  "comp_time_vec",        "identity", "0.0", "best",
  "poFGGM", "ComputationalTime",  "comp_time_vec_2",      "identity", "0.0", "10^0*mean",
  "poFGGM", "ComputationalTime",  "comp_time_vec_3",      "identity", "0.0", "10^2*mean",
  "poFGGM", "ComputationalTime",  "comp_time_vec_4",      "identity", "0.0", "10^4*mean",
  "poFGGM", "ComputationalTime",  "comp_time_vec_5",      "identity", "0.0", "10^6*mean"
)

apply_aggregator <- function(x, aggregator) {
  if (aggregator == "min_col") return(apply(x, 2, min))
  if (aggregator == "identity") return(x)
  stop("Unknown aggregator: ", aggregator)
}

extract_metrics_from_storage <- function(storage, cfg) {
  spec <- switch(
    cfg$scenario,
    gamma = metric_spec_gamma,
    alpha = metric_spec_alpha,
    metric_spec_default
  )
  
  purrr::map_dfr(seq_len(nrow(spec)), function(i) {
    obj_name <- spec$storage_name[i]
    if (!obj_name %in% names(storage)) {
      warning("Missing object in storage: ", obj_name)
      return(NULL)
    }
    
    values <- apply_aggregator(storage[[obj_name]], spec$aggregator[i])
    new_metric_row(
      method = spec$method[i],
      metric = spec$metric[i],
      val = values,
      cfg = cfg,
      gamma = spec$gamma[i],
      alpha = spec$alpha[i]
    )
  })
}

load_storage_file <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (!exists("storage", envir = e, inherits = FALSE)) {
    stop("Object `storage` not found in file: ", path)
  }
  get("storage", envir = e, inherits = FALSE)
}

# ============================================================
# 3. Build METRICS table
# ============================================================

METRICS <- purrr::map_dfr(seq_along(CONFIGS_list), function(i_cfg_set) {
  CONFIGS <- CONFIGS_list[[i_cfg_set]]
  
  purrr::map_dfr(seq_len(nrow(CONFIGS)), function(iconfig) {
    cfg <- CONFIGS[iconfig, ]
    file_name <- sprintf("jcgs_simul_%s_config%s.RData", cfg$scenario, iconfig)
    file_path <- file.path(RESULTS_DIR, file_name)
    
    if (!file.exists(file_path)) {
      warning("Missing results file: ", file_path)
      return(NULL)
    }
    
    message("Loading: ", basename(file_path))
    storage <- load_storage_file(file_path)
    extract_metrics_from_storage(storage, cfg)
  })
})

write.csv(METRICS, file.path(TABLES_DIR, "simulation_metrics_raw.csv"), row.names = FALSE)

# ============================================================
# 4. Factor coding and summaries
# ============================================================

METRICS <- METRICS %>%
  mutate(
    method = factor(method, levels = c("Oracle", "Kraus", "poFGGM")),
    metric = factor(metric, levels = c("ThetaError", "AUC", "CurveError", "ComputationalTime")),
    scenario = factor(scenario, levels = c("p_and_n", "missing", "pev", "L", "part_sep", "d", "gamma", "alpha")),
    n = factor(as.character(n), levels = c("50", "100")),
    p_over_n = factor(as.character(p_over_n), levels = c("0.2", "1.2")),
    L = factor(as.character(L), levels = c("5", "9")),
    pi_po = factor(as.character(pi_po), levels = c("0.3", "0.5", "0.7")),
    pi_d = factor(as.character(pi_d), levels = c("0.3", "0.5", "0.7")),
    d = factor(as.character(d), levels = c("50", "100", "200")),
    pev = factor(as.character(pev), levels = c("0.9", "0.95", "0.99"), labels = c("0.90", "0.95", "0.99")),
    gamma = factor(as.character(gamma), levels = c("0.0", "0.5", "1.0")),
    perc_theta_share = factor(as.character(perc_theta_share), levels = c("1", "0.5"), labels = c("1.0", "0.5")),
    alpha = factor(alpha, levels = c("best", "10^0*mean", "10^2*mean", "10^4*mean", "10^6*mean", "none")),
    graph_type = factor(graph_type, levels = c("star", "band", "smallworld"))
  )

write.csv(METRICS, file.path(TABLES_DIR, "simulation_metrics_tidy.csv"), row.names = FALSE)

# ============================================================
# 5. Generic summaries and plotting helpers
# ============================================================

summarise_metric_table <- function(data, time_in_minutes = FALSE) {
  data %>%
    group_by(n, p_over_n, pev, L, pi_po, pi_d, d, gamma, alpha, scenario, metric, method, graph_type, perc_theta_share) %>%
    summarise(
      mean = mean(if (time_in_minutes) val / 60 else val, na.rm = TRUE),
      sd = sd(if (time_in_minutes) val / 60 else val, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(summary = sprintf("%.3f (%.3f)", mean, sd))
}

metric_labeller <- as_labeller(c(
  ThetaError = "min[gamma[1]]~Err[Theta]",
  AUC = "AUC[Theta]",
  CurveError = "min[gamma[1]]~Err[X]"
), label_parsed)

method_fill <- c("Oracle" = "gray10", "Kraus" = "gray50", "poFGGM" = "gray90")

theme_boxpaper <- function() {
  theme_bw() +
    theme(
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white")
    )
}

plot_method_boxes <- function(data, x_var, facet_rows, facet_cols = NULL,
                              x_lab = NULL, y_lab = "",
                              fill_values = method_fill,
                              title = NULL,
                              scales = "free_y",
                              x_text_blank = FALSE,
                              parsed_col_labeller = NULL) {
  p <- ggplot(data, aes_string(x = x_var, y = "val", fill = "method")) +
    geom_boxplot(outlier.shape = NA) +
    scale_fill_manual(values = fill_values) +
    theme_boxpaper() +
    labs(x = x_lab, y = y_lab, title = title)
  
  row_labellers <- setNames(list(metric_labeller), facet_rows)
  if (!is.null(facet_cols) && !is.null(parsed_col_labeller)) {
    all_labellers <- c(row_labellers, setNames(list(parsed_col_labeller), facet_cols))
    p <- p + facet_grid(
      rows = vars(!!rlang::sym(facet_rows)),
      cols = vars(!!rlang::sym(facet_cols)),
      scales = scales,
      labeller = do.call(labeller, all_labellers)
    )
  } else if (!is.null(facet_cols)) {
    p <- p + facet_grid(
      rows = vars(!!rlang::sym(facet_rows)),
      cols = vars(!!rlang::sym(facet_cols)),
      scales = scales,
      labeller = labeller(!!facet_rows := metric_labeller)
    )
  } else {
    p <- p + facet_grid(
      rows = vars(!!rlang::sym(facet_rows)),
      scales = scales,
      labeller = labeller(!!facet_rows := metric_labeller)
    )
  }
  
  if (x_text_blank) {
    p <- p + theme(axis.text.x = element_blank())
  }
  
  p
}

plot_alpha_sensitivity <- function(data, title = NULL) {
  alpha_fill <- c(
    "best" = "gray5",
    "10^0*mean" = "gray25",
    "10^2*mean" = "gray50",
    "10^4*mean" = "gray75",
    "10^6*mean" = "gray95"
  )
  
  ggplot(data, aes(x = alpha, y = val, fill = alpha)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(rows = vars(metric), scales = "free_y", labeller = labeller(metric = metric_labeller)) +
    scale_fill_manual(
      values = alpha_fill,
      breaks = names(alpha_fill),
      labels = c(
        expression(hat(alpha)),
        expression(10^0 * bar(alpha)),
        expression(10^2 * bar(alpha)),
        expression(10^4 * bar(alpha)),
        expression(10^6 * bar(alpha))
      )
    ) +
    theme_boxpaper() +
    labs(x = NULL, y = NULL, title = title) +
    theme(
      axis.text.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14, face = "italic")
    )
}

plot_gamma_sensitivity <- function(data, title = NULL) {
  gamma_col_labeller <- as_labeller(c(
    "1.0" = "'100%'",
    "0.5" = "'50%'"
  ), label_parsed)
  
  ggplot(data, aes(x = gamma, y = val, fill = gamma)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(
      rows = vars(metric), cols = vars(perc_theta_share), scales = "free_y",
      labeller = labeller(metric = metric_labeller, perc_theta_share = gamma_col_labeller)
    ) +
    scale_fill_manual(
      values = c("0.0" = "gray10", "0.5" = "gray50", "1.0" = "gray90"),
      labels = c(
        expression(gamma[2] == 0.0),
        expression(gamma[2] == 0.5),
        expression(gamma[2] == 1.0)
      )
    ) +
    theme_boxpaper() +
    labs(x = NULL, y = NULL, title = title) +
    theme(
      axis.text.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14, face = "italic")
    )
}

# ============================================================
# 6. Tables
# ============================================================

# Main table: n and p_over_n, small-world, n = 100
main_np_table <- bind_rows(
  summarise_metric_table(
    METRICS %>%
      filter(scenario == "p_and_n", graph_type == "smallworld", metric %in% c("ThetaError", "AUC", "CurveError"))
  ),
  summarise_metric_table(
    METRICS %>%
      filter(scenario == "p_and_n", graph_type == "smallworld", metric == "ComputationalTime"),
    time_in_minutes = TRUE
  )
) %>%
  filter(n == "100") %>%
  arrange(metric, n, p_over_n, method) %>%
  select(metric, method, n, p_over_n, summary)

write.csv(main_np_table, file.path(TABLES_DIR, "table_np_smallworld_n100.csv"), row.names = FALSE)

# Computational cost tables by graph type
for (gtype in c("smallworld", "star", "band")) {
  comp_tab <- summarise_metric_table(
    METRICS %>% filter(method == "poFGGM", metric == "ComputationalTime", graph_type == gtype),
    time_in_minutes = TRUE
  ) %>%
    arrange(graph_type, scenario) %>%
    select(graph_type, scenario, n, p_over_n, pi_po, pi_d, pev, L, d, gamma, alpha, summary)
  
  write.csv(comp_tab, file.path(TABLES_DIR, paste0("computational_cost_", gtype, ".csv")), row.names = FALSE)
}

# ============================================================
# 7. Main-paper figures (small-world)
# ============================================================

# Missingness figure
missing_col_labeller <- as_labeller(c(
  `0.3` = "pi[po] == 0.3",
  `0.5` = "pi[po] == 0.5",
  `0.7` = "pi[po] == 0.7"
), label_parsed)

p_missing_smallworld <- METRICS %>%
  filter(scenario == "missing", graph_type == "smallworld", metric %in% c("ThetaError", "AUC", "CurveError")) %>%
  ggplot(aes(x = pi_d, y = val, fill = method)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(
    rows = vars(metric), cols = vars(pi_po), scales = "free_y",
    labeller = labeller(metric = metric_labeller, pi_po = missing_col_labeller)
  ) +
  scale_fill_manual(values = method_fill) +
  theme_boxpaper() +
  labs(x = expression(pi[w]), y = "")

ggsave(file.path(FIGURES_DIR, "missing.pdf"), p_missing_smallworld,
       device = "pdf", width = 10, height = 9, units = "in", dpi = 300)

# PEV and L figure
pev_col_labeller <- as_labeller(c(
  `0.90` = "PEV == 0.90",
  `0.95` = "PEV == 0.95",
  `0.99` = "PEV == 0.99"
), label_parsed)

L_col_labeller <- as_labeller(c(
  `5` = "L == 5",
  `9` = "L == 9"
), label_parsed)

p_pev <- METRICS %>%
  filter(scenario == "pev", graph_type == "smallworld", metric %in% c("ThetaError", "AUC", "CurveError")) %>%
  ggplot(aes(x = method, y = val, fill = method)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(
    rows = vars(metric), cols = vars(pev), scales = "free_y",
    labeller = labeller(metric = metric_labeller, pev = pev_col_labeller)
  ) +
  scale_fill_manual(values = method_fill) +
  theme_boxpaper() +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank())

p_L <- METRICS %>%
  filter(scenario == "L", graph_type == "smallworld", metric %in% c("ThetaError", "AUC", "CurveError")) %>%
  ggplot(aes(x = method, y = val, fill = method)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(
    rows = vars(metric), cols = vars(L), scales = "free_y",
    labeller = labeller(metric = metric_labeller, L = L_col_labeller)
  ) +
  scale_fill_manual(values = method_fill) +
  theme_boxpaper() +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank())

combined_pev_L <- (p_pev + theme(strip.text.y = element_blank(), strip.background.y = element_blank())) +
  p_L +
  plot_layout(widths = c(1.1, 1), guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(FIGURES_DIR, "pev_l.pdf"), combined_pev_L,
       device = "pdf", width = 12, height = 8, units = "in", dpi = 300)

# Sensitivity figure: small-world
p_alpha_sw <- METRICS %>%
  filter(
    scenario == "alpha", method == "poFGGM",
    alpha %in% c("best", "10^0*mean", "10^2*mean", "10^4*mean", "10^6*mean"),
    graph_type == "smallworld", metric %in% c("ThetaError", "AUC", "CurveError")
  ) %>%
  plot_alpha_sensitivity(title = "(i)")

p_gamma_sw <- METRICS %>%
  filter(scenario == "gamma", method == "poFGGM", graph_type == "smallworld", metric %in% c("ThetaError", "AUC", "CurveError")) %>%
  plot_gamma_sensitivity(title = "(ii)")

p_d_sw <- METRICS %>%
  filter(scenario == "d", graph_type == "smallworld", metric %in% c("ThetaError", "AUC", "CurveError")) %>%
  ggplot(aes(x = method, y = val, fill = method)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(rows = vars(metric), cols = vars(d), scales = "free_y", labeller = labeller(metric = metric_labeller)) +
  scale_fill_manual(values = method_fill) +
  theme_boxpaper() +
  labs(x = NULL, y = NULL, title = "(iii)") +
  theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 14, face = "italic"))

p_part_sep_sw <- METRICS %>%
  filter(scenario == "part_sep", graph_type == "smallworld", metric %in% c("ThetaError", "AUC", "CurveError")) %>%
  ggplot(aes(x = method, y = val, fill = method)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(rows = vars(metric), scales = "free_y", labeller = labeller(metric = metric_labeller)) +
  scale_fill_manual(values = method_fill) +
  theme_boxpaper() +
  labs(x = NULL, y = NULL, title = "(iv)") +
  theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 14, face = "italic"))

combined_sens_sw <- ((p_alpha_sw + theme(strip.text.y = element_blank(), strip.background.y = element_blank())) + p_gamma_sw) +
  ((p_d_sw + theme(strip.text.y = element_blank(), strip.background.y = element_blank())) +
     p_part_sep_sw +
     plot_layout(guides = "collect", design = c(area(1, 1, 1, 2), area(1, 3, 1, 3))) &
     theme(legend.position = "bottom")) +
  plot_layout(
    design = c(
      area(1, 1, 1, 1),
      area(1, 2, 1, 3),
      area(2, 1, 2, 3)
    )
  )

ggsave(file.path(FIGURES_DIR, "sens_an.pdf"), combined_sens_sw,
       device = "pdf", width = 14, height = 14, units = "in", dpi = 300)

# ============================================================
# 8. Supplementary figures by graph type
# ============================================================

# Helper to save patchwork object
save_pdf <- function(plot_obj, filename, width, height) {
  ggsave(file.path(FIGURES_DIR, filename), plot_obj,
         device = "pdf", width = width, height = height,
         units = "in", dpi = 300, scale = .7)
}

p_n50 <- p_n100 <- setNames(vector(mode = "list", length = 2L), nm = c("star", "band"))

# p_and_n: star and band
for (gtype in c("star", "band")) {
  p_n50[[gtype]] <- METRICS %>%
    filter(scenario == "p_and_n", graph_type == gtype, metric %in% c("ThetaError", "AUC", "CurveError"), n == "50") %>%
    ggplot(aes(x = method, y = val, fill = method)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(
      rows = vars(metric), cols = vars(p_over_n), scales = "free_y",
      labeller = labeller(
        metric = metric_labeller,
        p_over_n = as_labeller(c(`0.2` = "p/n == 0.2", `1.2` = "p/n == 1.2"), label_parsed)
      )
    ) +
    scale_fill_manual(values = method_fill) +
    theme_boxpaper() +
    labs(x = NULL, y = ifelse(gtype == "star", "Star graph structure", "Banded graph structure"), title = expression(n == 50)) +
    theme(
      axis.text.x = element_blank(),
      axis.title = element_text(size = 14, face = "italic"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "italic")
    )
  
  p_n100[[gtype]] <- METRICS %>%
    filter(scenario == "p_and_n", graph_type == gtype, metric %in% c("ThetaError", "AUC", "CurveError"), n == "100", !is.na(val)) %>%
    ggplot(aes(x = method, y = val, fill = method)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(
      rows = vars(metric), cols = vars(p_over_n), scales = "free_y",
      labeller = labeller(
        metric = metric_labeller,
        p_over_n = as_labeller(c(`0.2` = "p/n == 0.2", `1.2` = "p/n == 1.2"), label_parsed)
      )
    ) +
    scale_fill_manual(values = method_fill) +
    theme_boxpaper() +
    labs(x = NULL, y = NULL, title = expression(n == 100)) +
    theme(
      axis.text.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14, face = "italic")
    )
  
  combined_np <- (p_n50[[gtype]] + theme(strip.text.y = element_blank(), strip.background.y = element_blank())) +
    p_n100[[gtype]] +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  save_pdf(combined_np, paste0("supp_np_", gtype, ".pdf"), width = 14, height = 8)
}

combined_np_graph <- (p_n50[["star"]] + theme(strip.text.y = element_blank(), strip.background.y = element_blank())) + 
  p_n100[["star"]] + (p_n50[["band"]] + theme(strip.text.y = element_blank(), strip.background.y = element_blank(),
                   strip.text.x = element_blank(), strip.background.x = element_blank())) + 
  (p_n100[["band"]] + theme(strip.text.x = element_blank(), strip.background.x = element_blank())) + 
  plot_layout(
    guides = "collect",
    design = c(
      area(1, 1, 1, 1),
      area(1, 2, 1, 2),
      area(2, 1, 2, 1),
      area(2, 2, 2, 2)
    )
  ) &
  theme(legend.position = "bottom")

combined_np_graph

save_pdf(combined_np_graph, "supp_np_graphs.pdf", width = 14, height = 14)


p_missing <- setNames(vector(mode = "list", length = 2L), nm = c("star", "band"))

# missingness: star and band
for (gtype in c("star", "band")) {
  p_missing[[gtype]] <- METRICS %>%
    filter(scenario == "missing", graph_type == gtype, metric %in% c("ThetaError", "AUC", "CurveError")) %>%
    ggplot(aes(x = pi_d, y = val, fill = method)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(
      rows = vars(metric), cols = vars(pi_po), scales = "free_y",
      labeller = labeller(metric = metric_labeller, pi_po = missing_col_labeller)
    ) +
    scale_fill_manual(values = method_fill) +
    theme_boxpaper() +
    labs(x = expression(pi[w]), y = ifelse(gtype == "star", "Star graph structure", "Banded graph structure"))
  
  save_pdf(p_missing[[gtype]], paste0("supp_missing_", gtype, ".pdf"), width = 10, height = 9)
}

combined_missing_graph <- (p_missing[["star"]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank())) +
  (p_missing[["star"]] + theme(strip.text.x = element_blank(), strip.background.x = element_blank())) + 
  plot_layout(
    guides = "collect",
    design = c(
      area(1, 1, 1, 1),
      area(2, 1, 2, 1)
    )
  ) &
  theme(legend.position = "bottom")

combined_missing_graph

save_pdf(combined_missing_graph, "supp_missing_graphs.pdf", width = 10, height = 11.5)


p_pev_g <- p_L_g <- setNames(vector(mode = "list", length = 2L), nm = c("star", "band"))

# pev and L: star and band
for (gtype in c("star", "band")) {
  y_lab <- ifelse(gtype == "star", "Star graph structure", "Banded graph structure")
  
  p_pev_g[[gtype]] <- METRICS %>%
    filter(scenario == "pev", graph_type == gtype, metric %in% c("ThetaError", "AUC", "CurveError")) %>%
    ggplot(aes(x = method, y = val, fill = method)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(rows = vars(metric), cols = vars(pev), scales = "free_y",
               labeller = labeller(metric = metric_labeller, pev = pev_col_labeller)) +
    scale_fill_manual(values = method_fill) +
    theme_boxpaper() +
    labs(x = NULL, y = y_lab) +
    theme(axis.text.x = element_blank())
  
  p_L_g[[gtype]] <- METRICS %>%
    filter(scenario == "L", graph_type == gtype, metric %in% c("ThetaError", "AUC", "CurveError")) %>%
    ggplot(aes(x = method, y = val, fill = method)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(rows = vars(metric), cols = vars(L), scales = "free_y",
               labeller = labeller(metric = metric_labeller, L = L_col_labeller)) +
    scale_fill_manual(values = method_fill) +
    theme_boxpaper() +
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank())
  
  combined_pevL_g <- (p_pev_g[[gtype]] + theme(strip.text.y = element_blank(), strip.background.y = element_blank())) +
    p_L_g[[gtype]] +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  save_pdf(combined_pevL_g, paste0("supp_pev_l_", gtype, ".pdf"), width = 14, height = 8)
}

combined_pevL_graph <- (p_pev_g[["star"]] + theme(strip.text.y = element_blank(), strip.background.y = element_blank())) + 
  p_L_g[["star"]] + (p_pev_g[["band"]] + theme(strip.text.y = element_blank(), strip.background.y = element_blank(),
                                               strip.text.x = element_blank(), strip.background.x = element_blank(),
                                               axis.text.x = element_blank(), axis.title.x = element_blank())) + 
  (p_L_g[["band"]] + theme(strip.text.x = element_blank(), strip.background.x = element_blank(),
                           axis.text.x = element_blank(), axis.title.x = element_blank())) + 
  plot_layout(
    guides = "collect",
    design = c(
      area(1, 1, 1, 1),
      area(1, 2, 1, 2),
      area(2, 1, 2, 1),
      area(2, 2, 2, 2)
    )
  ) &
  theme(legend.position = "bottom")

save_pdf(combined_pevL_graph, "supp_pev_l_graphs.pdf", width = 14, height = 14)


# sensitivity: star and band
for (gtype in c("star", "band")) {
  p_alpha_g <- METRICS %>%
    filter(
      scenario == "alpha", method == "poFGGM",
      alpha %in% c("best", "10^0*mean", "10^2*mean", "10^4*mean", "10^6*mean"),
      graph_type == gtype, metric %in% c("ThetaError", "AUC", "CurveError")
    ) %>%
    plot_alpha_sensitivity(title = "(i)")
  
  p_gamma_g <- METRICS %>%
    filter(scenario == "gamma", method == "poFGGM", graph_type == gtype, metric %in% c("ThetaError", "AUC", "CurveError")) %>%
    plot_gamma_sensitivity(title = "(ii)")
  
  p_d_g <- METRICS %>%
    filter(scenario == "d", graph_type == gtype, metric %in% c("ThetaError", "AUC", "CurveError")) %>%
    ggplot(aes(x = method, y = val, fill = method)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(rows = vars(metric), cols = vars(d), scales = "free_y", labeller = labeller(metric = metric_labeller)) +
    scale_fill_manual(values = method_fill) +
    theme_boxpaper() +
    labs(x = NULL, y = NULL, title = "(iii)") +
    theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 14, face = "italic"))
  
  p_part_sep_g <- METRICS %>%
    filter(scenario == "part_sep", graph_type == gtype, metric %in% c("ThetaError", "AUC", "CurveError")) %>%
    ggplot(aes(x = method, y = val, fill = method)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(rows = vars(metric), scales = "free_y", labeller = labeller(metric = metric_labeller)) +
    scale_fill_manual(values = method_fill) +
    theme_boxpaper() +
    labs(x = NULL, y = NULL, title = "(iv)") +
    theme(axis.text.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 14, face = "italic"))
  
  combined_sens_g <- ((p_alpha_g + theme(strip.text.y = element_blank(), strip.background.y = element_blank())) + p_gamma_g) +
    ((p_d_g + theme(strip.text.y = element_blank(), strip.background.y = element_blank())) +
       p_part_sep_g +
       plot_layout(guides = "collect", design = c(area(1, 1, 1, 2), area(1, 3, 1, 3))) &
       theme(legend.position = "bottom")) +
    plot_layout(
      design = c(
        area(1, 1, 1, 1),
        area(1, 2, 1, 3),
        area(2, 1, 2, 3)
      )
    ) +
    plot_annotation(title = ifelse(gtype == "star", "Star graph structure", "Banded graph structure")) &
    theme(plot.title = element_text(hjust = 0.5, face = "italic"))
  
  save_pdf(combined_sens_g, paste0("supp_sens_", gtype, ".pdf"), width = 14, height = 14)
}

# ============================================================
# 9. Console summaries
# ============================================================

message("Summary tables and figures created.")
message("Number of rows in METRICS: ", nrow(METRICS))
message("Saved tidy metrics to: ", file.path(TABLES_DIR, "simulation_metrics_tidy.csv"))
