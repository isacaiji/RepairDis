############################################################
# Scatter + fitted-line drug sensitivity plots
# Each dot is one TCGA sample.
# X axis: continuous MO-DDRscore or DPRS
# Y axis: predicted IC50 from oncoPredict
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260531)

pkgs <- c("data.table", "dplyr", "tidyr", "ggplot2", "scales")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

PROJECT_DIR <- file.path("D:/R_workspace", intToUtf8(c(0x8bc4, 0x5206)), "AD_DDR_project")
DRUG_DIR <- file.path(PROJECT_DIR, "drug")
TABLE_DIR <- file.path(DRUG_DIR, "03-res", "tables")
DATA_DIR <- file.path(DRUG_DIR, "02-data")
FIG_DIR <- file.path(DRUG_DIR, "03-res", "figures", "publication_ready_v3_scatter_fit")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

############################################################
# Helpers
############################################################

save_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file)
}

p_star <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

theme_scatter <- function(base_size = 8) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.major = element_line(color = "#E8EDF3", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "#2C3440", fill = NA, linewidth = 0.45),
      axis.text = element_text(color = "#111827"),
      axis.title = element_text(color = "#111827", face = "bold"),
      plot.title = element_text(color = "#111827", face = "bold", hjust = 0.5, size = rel(1.08)),
      strip.background = element_rect(fill = "#F3F6FA", color = "#AEB9CA", linewidth = 0.35),
      strip.text = element_text(color = "#111827", face = "bold"),
      legend.title = element_text(color = "#111827", face = "bold"),
      legend.text = element_text(color = "#111827"),
      plot.margin = margin(7, 8, 7, 8)
    )
}

save_plot <- function(p, name, width, height) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = width, height = height, dpi = 600, bg = "white")
}

select_drugs_by_name <- function(tab, wanted) {
  lapply(wanted, function(nm) {
    hit <- tab %>% filter(toupper(DrugName) == toupper(nm)) %>% arrange(Spearman_FDR, Spearman_P)
    if (nrow(hit) == 0) return(NULL)
    hit[1, , drop = FALSE]
  }) %>%
    bind_rows() %>%
    distinct(DrugName, .keep_all = TRUE)
}

make_long <- function(pred, anno, selected, score_col, group_col, score_label) {
  pred %>%
    select(Sample, all_of(selected$DrugKey)) %>%
    pivot_longer(-Sample, names_to = "DrugKey", values_to = "Predicted_IC50") %>%
    left_join(anno, by = "Sample") %>%
    left_join(
      selected %>% select(DrugKey, DrugName, PathwayName, Spearman_Rho, Spearman_P, Spearman_FDR),
      by = "DrugKey"
    ) %>%
    mutate(
      Score = as.numeric(.data[[score_col]]),
      Group = factor(.data[[group_col]], levels = c("Low", "High")),
      Log10_IC50 = log10(Predicted_IC50 + 1),
      DrugName = factor(DrugName, levels = selected$DrugName),
      ScoreLabel = score_label
    ) %>%
    filter(is.finite(Score), is.finite(Log10_IC50), !is.na(Group))
}

label_df <- function(plot_df) {
  plot_df %>%
    group_by(DrugName) %>%
    summarise(
      x = quantile(Score, 0.04, na.rm = TRUE),
      y = quantile(Log10_IC50, 0.96, na.rm = TRUE),
      Spearman_Rho = dplyr::first(Spearman_Rho),
      Spearman_P = dplyr::first(Spearman_P),
      Spearman_FDR = dplyr::first(Spearman_FDR),
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0(
        "rho = ", sprintf("%.2f", Spearman_Rho),
        "\nFDR = ", formatC(Spearman_FDR, format = "e", digits = 1),
        " ", p_star(Spearman_FDR)
      )
    )
}

scatter_fit_panel <- function(plot_df, title, x_lab, color_title, ncol = 4) {
  labs <- label_df(plot_df)

  ggplot(plot_df, aes(Score, Log10_IC50)) +
    geom_point(aes(color = Group), size = 0.72, alpha = 0.58) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      linewidth = 0.62,
      color = "#111827",
      fill = "#8FA7C8",
      alpha = 0.20
    ) +
    geom_text(
      data = labs,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.55,
      lineheight = 0.92,
      color = "#111827"
    ) +
    facet_wrap(~DrugName, scales = "free_y", ncol = ncol) +
    scale_color_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
    labs(
      x = x_lab,
      y = "Predicted IC50 (log10)",
      title = title,
      color = color_title
    ) +
    theme_scatter(8) +
    theme(
      legend.position = "top",
      panel.spacing = unit(1.0, "lines")
    )
}

single_scatter <- function(plot_df, drug_name, title, x_lab, color_title) {
  df <- plot_df %>% filter(as.character(DrugName) == drug_name)
  labs <- label_df(df)

  ggplot(df, aes(Score, Log10_IC50)) +
    geom_point(aes(color = Group), size = 1.05, alpha = 0.62) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      linewidth = 0.72,
      color = "#111827",
      fill = "#8FA7C8",
      alpha = 0.22
    ) +
    geom_text(
      data = labs,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 3.1,
      lineheight = 0.95,
      color = "#111827"
    ) +
    scale_color_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
    labs(
      x = x_lab,
      y = "Predicted IC50 (log10)",
      title = title,
      color = color_title
    ) +
    theme_scatter(9) +
    theme(legend.position = "top")
}

############################################################
# Load data
############################################################

pred <- fread(file.path(TABLE_DIR, "Predicted_IC50_TCGA_all_drugs.csv"), data.table = FALSE, check.names = FALSE)
anno <- fread(file.path(DATA_DIR, "Drug_analysis_sample_annotation.csv"), data.table = FALSE, check.names = FALSE)
mo_cor <- fread(file.path(TABLE_DIR, "Drug_sensitivity_MO_DDRscore_spearman.csv"), data.table = FALSE, check.names = FALSE)
dprs_cor <- fread(file.path(TABLE_DIR, "Drug_sensitivity_DPRS_spearman.csv"), data.table = FALSE, check.names = FALSE)

############################################################
# Selected drugs
############################################################

mo_selected <- select_drugs_by_name(
  mo_cor,
  c("AZD7762", "MK-1775", "Wee1 Inhibitor", "AZD6738", "Olaparib", "Talazoparib", "Cisplatin", "Paclitaxel", "Pictilisib", "AZD8055", "Selumetinib", "Temsirolimus")
)

dprs_selected <- select_drugs_by_name(
  dprs_cor,
  c("AZD7762", "AZD6738", "MK-1775", "Wee1 Inhibitor", "Docetaxel", "Paclitaxel", "Gemcitabine", "Erlotinib", "Gefitinib", "Selumetinib", "VE-822", "Cisplatin")
)

mo_top <- mo_cor %>%
  arrange(Spearman_FDR, Spearman_P) %>%
  group_by(DrugName) %>%
  slice(1) %>%
  ungroup() %>%
  filter(Spearman_FDR < 0.05) %>%
  slice_head(n = 12)

dprs_top <- dprs_cor %>%
  arrange(Spearman_FDR, Spearman_P) %>%
  group_by(DrugName) %>%
  slice(1) %>%
  ungroup() %>%
  filter(Spearman_FDR < 0.05) %>%
  slice_head(n = 12)

############################################################
# Plot data
############################################################

mo_df <- make_long(
  pred, anno, mo_selected,
  score_col = "MO_DDRscore_raw",
  group_col = "MO_DDRscore_group",
  score_label = "MO-DDRscore"
)

dprs_df <- make_long(
  pred, anno, dprs_selected,
  score_col = "DPRS",
  group_col = "DPRS_RiskGroup",
  score_label = "DPRS"
)

mo_top_df <- make_long(
  pred, anno, mo_top,
  score_col = "MO_DDRscore_raw",
  group_col = "MO_DDRscore_group",
  score_label = "MO-DDRscore"
)

dprs_top_df <- make_long(
  pred, anno, dprs_top,
  score_col = "DPRS",
  group_col = "DPRS_RiskGroup",
  score_label = "DPRS"
)

save_csv(mo_df, file.path(FIG_DIR, "Drug_A_MO_DDRscore_selected_scatter_data.csv"))
save_csv(dprs_df, file.path(FIG_DIR, "Drug_B_DPRS_selected_scatter_data.csv"))
save_csv(mo_top_df, file.path(FIG_DIR, "Drug_C_MO_DDRscore_top_scatter_data.csv"))
save_csv(dprs_top_df, file.path(FIG_DIR, "Drug_D_DPRS_top_scatter_data.csv"))

############################################################
# Multi-panel scatter-fit figures
############################################################

p_mo <- scatter_fit_panel(
  mo_df,
  title = "MO-DDRscore and predicted drug sensitivity",
  x_lab = "MO-DDRscore",
  color_title = "MO-DDRscore group",
  ncol = 4
)
p_dprs <- scatter_fit_panel(
  dprs_df,
  title = "DPRS and predicted drug sensitivity",
  x_lab = "DPRS",
  color_title = "DPRS risk group",
  ncol = 4
)
p_mo_top <- scatter_fit_panel(
  mo_top_df,
  title = "Top MO-DDRscore-drug IC50 correlations",
  x_lab = "MO-DDRscore",
  color_title = "MO-DDRscore group",
  ncol = 4
)
p_dprs_top <- scatter_fit_panel(
  dprs_top_df,
  title = "Top DPRS-drug IC50 correlations",
  x_lab = "DPRS",
  color_title = "DPRS risk group",
  ncol = 4
)

save_plot(p_mo, "Drug_A_scatterfit_MO_DDRscore_selected_drugs", 9.6, 7.2)
save_plot(p_dprs, "Drug_B_scatterfit_DPRS_selected_drugs", 9.6, 7.2)
save_plot(p_mo_top, "Drug_C_scatterfit_MO_DDRscore_top_correlated_drugs", 9.6, 7.2)
save_plot(p_dprs_top, "Drug_D_scatterfit_DPRS_top_correlated_drugs", 9.6, 7.2)

############################################################
# Single-drug examples for main text
############################################################

for (drug in c("AZD7762", "Cisplatin", "Pictilisib", "Olaparib")) {
  if (drug %in% as.character(mo_df$DrugName)) {
    p <- single_scatter(
      mo_df,
      drug,
      title = paste0("MO-DDRscore vs ", drug, " predicted IC50"),
      x_lab = "MO-DDRscore",
      color_title = "MO-DDRscore group"
    )
    save_plot(p, paste0("Drug_single_MO_DDRscore_", gsub("[^A-Za-z0-9]+", "_", drug)), 4.2, 3.8)
  }
}

for (drug in c("AZD7762", "AZD6738", "Docetaxel", "Erlotinib")) {
  if (drug %in% as.character(dprs_df$DrugName)) {
    p <- single_scatter(
      dprs_df,
      drug,
      title = paste0("DPRS vs ", drug, " predicted IC50"),
      x_lab = "DPRS",
      color_title = "DPRS risk group"
    )
    save_plot(p, paste0("Drug_single_DPRS_", gsub("[^A-Za-z0-9]+", "_", drug)), 4.2, 3.8)
  }
}

summary_df <- data.frame(
  Item = c(
    "N_pdf",
    "N_png",
    "N_MO_selected_drugs",
    "N_DPRS_selected_drugs",
    "N_MO_top_drugs",
    "N_DPRS_top_drugs",
    "Figure_directory"
  ),
  Value = c(
    length(list.files(FIG_DIR, pattern = "\\.pdf$")),
    length(list.files(FIG_DIR, pattern = "\\.png$")),
    n_distinct(mo_df$DrugName),
    n_distinct(dprs_df$DrugName),
    n_distinct(mo_top_df$DrugName),
    n_distinct(dprs_top_df$DrugName),
    FIG_DIR
  )
)

save_csv(summary_df, file.path(FIG_DIR, "Drug_scatterfit_plot_summary.csv"))

cat("\nDone.\n")
cat("Figure directory:\n", FIG_DIR, "\n")
print(summary_df)
