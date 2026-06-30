############################################################
# Standardized drug sensitivity plots, v2
# Cleaner manuscript-style figures:
#   - no ns drugs in main effect plots
#   - duplicate drug names collapsed for display
#   - explicit significance stars
#   - comparable sensitivity index: log2(IC50_low / IC50_high)
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260531)

pkgs <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", "ComplexHeatmap",
  "circlize", "grid", "scales"
)

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
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(scales)
})

PROJECT_DIR <- file.path("D:/R_workspace", intToUtf8(c(0x8bc4, 0x5206)), "AD_DDR_project")
DRUG_DIR <- file.path(PROJECT_DIR, "drug")
TABLE_DIR <- file.path(DRUG_DIR, "03-res", "tables")
DATA_DIR <- file.path(DRUG_DIR, "02-data")
FIG_DIR <- file.path(DRUG_DIR, "03-res", "figures", "publication_ready_v2_standard")
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
    TRUE ~ ""
  )
}

theme_std <- function(base_size = 8) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.major = element_line(color = "#E8EDF3", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "#2C3440", fill = NA, linewidth = 0.45),
      axis.text = element_text(color = "#111827"),
      axis.title = element_text(color = "#111827", face = "bold"),
      plot.title = element_text(color = "#111827", face = "bold", hjust = 0.5, size = rel(1.05)),
      strip.background = element_rect(fill = "#F3F6FA", color = "#AEB9CA", linewidth = 0.35),
      strip.text = element_text(color = "#111827", face = "bold"),
      legend.title = element_text(color = "#111827", face = "bold"),
      legend.text = element_text(color = "#111827")
    )
}

save_plot <- function(p, name, width, height) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = width, height = height, dpi = 600, bg = "white")
}

clean_group <- function(x, analysis) {
  x %>%
    mutate(
      Analysis = analysis,
      Median_High = ifelse(Group1 == "High", Median_Group1, Median_Group2),
      Median_Low = ifelse(Group1 == "Low", Median_Group1, Median_Group2),
      Mean_High = ifelse(Group1 == "High", Mean_Group1, Mean_Group2),
      Mean_Low = ifelse(Group1 == "Low", Mean_Group1, Mean_Group2),
      Log2_IC50_Ratio_High_vs_Low = log2((Median_High + 1e-8) / (Median_Low + 1e-8)),
      Sensitivity_Index = log2((Median_Low + 1e-8) / (Median_High + 1e-8)),
      Direction = ifelse(Sensitivity_Index > 0, "High group more sensitive", "Low group more sensitive"),
      Signif = p_star(Wilcox_FDR)
    )
}

collapse_for_display <- function(tab) {
  tab %>%
    arrange(Wilcox_FDR, Wilcox_P) %>%
    group_by(DrugName) %>%
    slice(1) %>%
    ungroup()
}

choose_named_drugs <- function(tab, wanted) {
  out <- lapply(wanted, function(nm) {
    hit <- tab %>% filter(toupper(DrugName) == toupper(nm)) %>% arrange(Wilcox_FDR, Wilcox_P)
    if (nrow(hit) == 0) return(NULL)
    hit[1, , drop = FALSE]
  }) %>% bind_rows()
  out %>% distinct(DrugName, .keep_all = TRUE)
}

make_box_data <- function(pred, anno, selected, group_col) {
  pred %>%
    select(Sample, all_of(selected$DrugKey)) %>%
    pivot_longer(-Sample, names_to = "DrugKey", values_to = "Predicted_IC50") %>%
    left_join(anno, by = "Sample") %>%
    left_join(selected %>% select(DrugKey, DrugName, Wilcox_FDR, Signif, Direction), by = "DrugKey") %>%
    mutate(
      Group = factor(.data[[group_col]], levels = c("Low", "High")),
      Log10_IC50 = log10(Predicted_IC50 + 1),
      DrugName = factor(DrugName, levels = selected$DrugName)
    ) %>%
    filter(!is.na(Group), is.finite(Log10_IC50))
}

boxplot_panel <- function(plot_df, selected, title) {
  lab_df <- plot_df %>%
    group_by(DrugName) %>%
    summarise(
      ymin = min(Log10_IC50, na.rm = TRUE),
      ymax = max(Log10_IC50, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      yr = pmax(ymax - ymin, 0.05),
      y = ymax + 0.07 * yr,
      y_text = ymax + 0.11 * yr
    ) %>%
    left_join(selected %>% select(DrugName, Signif), by = "DrugName")

  ggplot(plot_df, aes(Group, Log10_IC50, fill = Group)) +
    geom_violin(width = 0.72, alpha = 0.34, color = NA, trim = TRUE) +
    geom_boxplot(width = 0.32, outlier.shape = NA, linewidth = 0.35, alpha = 0.88) +
    geom_point(
      aes(color = Group),
      position = position_jitter(width = 0.08, height = 0, seed = 20260531),
      size = 0.23,
      alpha = 0.18,
      show.legend = FALSE
    ) +
    geom_segment(
      data = lab_df,
      aes(x = 1, xend = 2, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.32,
      color = "#111827"
    ) +
    geom_text(
      data = lab_df,
      aes(x = 1.5, y = y_text, label = Signif),
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "bold",
      color = "#111827"
    ) +
    facet_wrap(~DrugName, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
    scale_color_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
    labs(
      x = NULL,
      y = "Predicted IC50 (log10)",
      title = title,
      fill = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_std(8) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(face = "bold"),
      panel.spacing = unit(1.0, "lines"),
      plot.margin = margin(7, 7, 7, 16)
    )
}

effect_dotplot <- function(tab, title, top_n = 28) {
  plot_tab <- tab %>%
    collapse_for_display() %>%
    filter(Wilcox_FDR < 0.05) %>%
    arrange(Wilcox_FDR, Wilcox_P) %>%
    slice_head(n = top_n) %>%
    mutate(
      DrugName = factor(DrugName, levels = rev(DrugName)),
      Direction = factor(Direction, levels = c("High group more sensitive", "Low group more sensitive")),
      NegLogFDR = pmin(-log10(Wilcox_FDR + 1e-300), 18)
    )

  ggplot(plot_tab, aes(Sensitivity_Index, DrugName)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, color = "#7C8798") +
    geom_segment(aes(x = 0, xend = Sensitivity_Index, yend = DrugName), linewidth = 0.5, color = "#8893A5") +
    geom_point(
      aes(size = NegLogFDR, fill = Direction),
      shape = 21,
      color = "#111827",
      stroke = 0.25,
      alpha = 0.95
    ) +
    geom_text(
      aes(label = Signif, x = Sensitivity_Index + ifelse(Sensitivity_Index >= 0, 0.035, -0.035)),
      size = 2.4,
      fontface = "bold",
      color = "#111827"
    ) +
    facet_grid(PathwayName ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(
      values = c("High group more sensitive" = "#C84630", "Low group more sensitive" = "#3B75AF"),
      labels = c("High group more sensitive", "Low group more sensitive")
    ) +
    scale_size_continuous(name = "-log10(FDR)", range = c(2.0, 6.2)) +
    labs(
      x = "Sensitivity index: log2(median IC50 Low / median IC50 High)",
      y = NULL,
      title = title,
      fill = NULL
    ) +
    theme_std(8) +
    theme(
      legend.position = "right",
      strip.text.y = element_text(angle = 0, size = 7),
      panel.spacing.y = unit(0.28, "lines")
    )
}

cor_lollipop <- function(cor_tab, title, top_n = 24) {
  plot_tab <- cor_tab %>%
    arrange(Spearman_FDR, Spearman_P) %>%
    group_by(DrugName) %>%
    slice(1) %>%
    ungroup() %>%
    filter(Spearman_FDR < 0.05) %>%
    slice_head(n = top_n) %>%
    mutate(
      Signif = p_star(Spearman_FDR),
      DrugName = factor(DrugName, levels = rev(DrugName)),
      NegLogFDR = pmin(-log10(Spearman_FDR + 1e-300), 18),
      Direction = ifelse(Spearman_Rho < 0, "Higher score: lower IC50", "Higher score: higher IC50")
    )

  ggplot(plot_tab, aes(Spearman_Rho, DrugName)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, color = "#7C8798") +
    geom_segment(aes(x = 0, xend = Spearman_Rho, yend = DrugName), linewidth = 0.55, color = "#8893A5") +
    geom_point(aes(size = NegLogFDR, fill = Direction), shape = 21, color = "#111827", stroke = 0.25) +
    geom_text(
      aes(label = Signif, x = Spearman_Rho + ifelse(Spearman_Rho >= 0, 0.018, -0.018)),
      size = 2.3,
      fontface = "bold",
      color = "#111827"
    ) +
    facet_grid(PathwayName ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = c("Higher score: lower IC50" = "#C84630", "Higher score: higher IC50" = "#3B75AF")) +
    scale_size_continuous(name = "-log10(FDR)", range = c(2.0, 5.8)) +
    labs(x = "Spearman rho with predicted IC50", y = NULL, title = title, fill = NULL) +
    theme_std(8) +
    theme(
      legend.position = "right",
      strip.text.y = element_text(angle = 0, size = 7),
      panel.spacing.y = unit(0.28, "lines")
    )
}

combined_heatmap <- function(mo, dprs) {
  mo2 <- mo %>%
    collapse_for_display() %>%
    select(DrugName, PathwayName, MO_Index = Sensitivity_Index, MO_FDR = Wilcox_FDR)
  dp2 <- dprs %>%
    collapse_for_display() %>%
    select(DrugName, DPRS_Index = Sensitivity_Index, DPRS_FDR = Wilcox_FDR)

  hm <- full_join(mo2, dp2, by = "DrugName") %>%
    mutate(
      MinFDR = pmin(MO_FDR, DPRS_FDR, na.rm = TRUE),
      PathwayName = ifelse(is.na(PathwayName), "Other", PathwayName)
    ) %>%
    filter(MinFDR < 0.05) %>%
    arrange(PathwayName, MinFDR)

  mat <- as.matrix(hm[, c("MO_Index", "DPRS_Index")])
  rownames(mat) <- hm$DrugName
  colnames(mat) <- c("MO-DDRscore", "DPRS")

  fdr <- as.matrix(hm[, c("MO_FDR", "DPRS_FDR")])
  rownames(fdr) <- hm$DrugName
  colnames(fdr) <- colnames(mat)

  max_abs <- quantile(abs(mat), 0.95, na.rm = TRUE)
  max_abs <- max(max_abs, 0.25, na.rm = TRUE)
  col_fun <- circlize::colorRamp2(c(-max_abs, 0, max_abs), c("#3B75AF", "white", "#C84630"))

  pathway_levels <- unique(hm$PathwayName)
  pathway_cols <- setNames(
    c("#355C7D", "#6C5B7B", "#C06C84", "#F67280", "#2A9D8F", "#E9C46A", "#8D99AE", "#457B9D")[seq_along(pathway_levels)],
    pathway_levels
  )

  ha <- rowAnnotation(
    Pathway = hm$PathwayName,
    col = list(Pathway = pathway_cols),
    gp = gpar(col = NA),
    annotation_name_gp = gpar(fontface = "bold", fontsize = 8)
  )

  ht <- Heatmap(
    mat,
    name = "Sensitivity\nindex",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = hm$PathwayName,
    left_annotation = ha,
    rect_gp = gpar(col = "white", lwd = 0.8),
    row_title_gp = gpar(fontsize = 8, fontface = "bold"),
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 9, fontface = "bold"),
    heatmap_legend_param = list(
      title_gp = gpar(fontface = "bold", fontsize = 8),
      labels_gp = gpar(fontsize = 7),
      legend_height = unit(32, "mm"),
      at = c(-round(max_abs, 2), 0, round(max_abs, 2)),
      labels = c("Low group\nmore sensitive", "0", "High group\nmore sensitive")
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      lab <- p_star(fdr[i, j])
      if (!is.na(fdr[i, j]) && fdr[i, j] < 0.05) {
        grid.text(lab, x, y, gp = gpar(fontsize = 7, fontface = "bold", col = "#111827"))
      }
    }
  )

  pdf(file.path(FIG_DIR, "Drug_E_standard_combined_sensitivity_heatmap.pdf"), width = 6.2, height = 8.4, useDingbats = FALSE)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  png(file.path(FIG_DIR, "Drug_E_standard_combined_sensitivity_heatmap.png"), width = 2800, height = 3800, res = 520, type = "cairo")
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  save_csv(hm, file.path(FIG_DIR, "Drug_E_standard_heatmap_data.csv"))
}

############################################################
# Load inputs
############################################################

pred <- fread(file.path(TABLE_DIR, "Predicted_IC50_TCGA_all_drugs.csv"), data.table = FALSE, check.names = FALSE)
anno <- fread(file.path(DATA_DIR, "Drug_analysis_sample_annotation.csv"), data.table = FALSE, check.names = FALSE)
mo <- fread(file.path(TABLE_DIR, "Drug_sensitivity_MO_DDRscore_group_comparison.csv"), data.table = FALSE, check.names = FALSE) %>%
  clean_group("MO-DDRscore")
dprs <- fread(file.path(TABLE_DIR, "Drug_sensitivity_DPRS_group_comparison.csv"), data.table = FALSE, check.names = FALSE) %>%
  clean_group("DPRS")
mo_cor <- fread(file.path(TABLE_DIR, "Drug_sensitivity_MO_DDRscore_spearman.csv"), data.table = FALSE, check.names = FALSE)
dprs_cor <- fread(file.path(TABLE_DIR, "Drug_sensitivity_DPRS_spearman.csv"), data.table = FALSE, check.names = FALSE)

save_csv(mo, file.path(FIG_DIR, "Drug_MO_DDRscore_plot_ready_v2.csv"))
save_csv(dprs, file.path(FIG_DIR, "Drug_DPRS_plot_ready_v2.csv"))

############################################################
# A/C. Boxplots for representative drugs
############################################################

mo_selected <- choose_named_drugs(
  mo,
  c("AZD7762", "MK-1775", "Wee1 Inhibitor", "AZD6738", "Olaparib", "Cisplatin", "Paclitaxel", "Pictilisib")
)
dprs_selected <- choose_named_drugs(
  dprs,
  c("AZD7762", "AZD6738", "MK-1775", "Wee1 Inhibitor", "Docetaxel", "Paclitaxel", "Gemcitabine", "Erlotinib")
)

mo_box <- make_box_data(pred, anno, mo_selected, "MO_DDRscore_group")
dprs_box <- make_box_data(pred, anno, dprs_selected, "DPRS_RiskGroup")

p_a <- boxplot_panel(mo_box, mo_selected, "MO-DDRscore groups")
p_c <- boxplot_panel(dprs_box, dprs_selected, "DPRS risk groups")
save_plot(p_a, "Drug_A_standard_MO_DDRscore_representative_boxplot", 7.6, 5.4)
save_plot(p_c, "Drug_C_standard_DPRS_representative_boxplot", 7.6, 5.4)
save_csv(mo_box, file.path(FIG_DIR, "Drug_A_standard_MO_DDRscore_boxplot_data.csv"))
save_csv(dprs_box, file.path(FIG_DIR, "Drug_C_standard_DPRS_boxplot_data.csv"))

############################################################
# B/D. Ranked sensitivity index plots
############################################################

p_b <- effect_dotplot(mo, "MO-DDRscore-associated drug sensitivity", top_n = 28)
p_d <- effect_dotplot(dprs, "DPRS-associated drug sensitivity", top_n = 28)
save_plot(p_b, "Drug_B_standard_MO_DDRscore_ranked_sensitivity_index", 7.8, 8.6)
save_plot(p_d, "Drug_D_standard_DPRS_ranked_sensitivity_index", 7.8, 8.6)

############################################################
# E. Combined heatmap
############################################################

combined_heatmap(mo, dprs)

############################################################
# F/G. Correlation plots
############################################################

p_f <- cor_lollipop(mo_cor, "MO-DDRscore correlation with predicted IC50", top_n = 24)
p_g <- cor_lollipop(dprs_cor, "DPRS correlation with predicted IC50", top_n = 24)
save_plot(p_f, "Drug_F_standard_MO_DDRscore_IC50_correlation", 7.4, 7.8)
save_plot(p_g, "Drug_G_standard_DPRS_IC50_correlation", 7.4, 7.8)

summary_df <- data.frame(
  Item = c(
    "N_pdf",
    "N_png",
    "N_MO_selected_box_drugs",
    "N_DPRS_selected_box_drugs",
    "N_MO_significant_unique_drugs",
    "N_DPRS_significant_unique_drugs",
    "Figure_directory"
  ),
  Value = c(
    length(list.files(FIG_DIR, pattern = "\\.pdf$")),
    length(list.files(FIG_DIR, pattern = "\\.png$")),
    nrow(mo_selected),
    nrow(dprs_selected),
    nrow(collapse_for_display(mo) %>% filter(Wilcox_FDR < 0.05)),
    nrow(collapse_for_display(dprs) %>% filter(Wilcox_FDR < 0.05)),
    FIG_DIR
  )
)

save_csv(summary_df, file.path(FIG_DIR, "Drug_standard_plot_summary.csv"))

cat("\nDone.\n")
cat("Figure directory:\n", FIG_DIR, "\n")
print(summary_df)
