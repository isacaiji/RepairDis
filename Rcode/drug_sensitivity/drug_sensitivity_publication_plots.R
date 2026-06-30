############################################################
# Publication-ready drug sensitivity plots
# Input: oncoPredict drug sensitivity tables
# Output: PDF + PNG figures
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260531)

pkgs <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", "ggpubr",
  "ComplexHeatmap", "circlize", "scales", "grid"
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
  library(ggpubr)
  library(ComplexHeatmap)
  library(circlize)
  library(scales)
  library(grid)
})

PROJECT_DIR <- file.path(
  "D:/R_workspace",
  intToUtf8(c(0x8bc4, 0x5206)),
  "AD_DDR_project"
)

DRUG_DIR <- file.path(PROJECT_DIR, "drug")
TABLE_DIR <- file.path(DRUG_DIR, "03-res", "tables")
DATA_DIR <- file.path(DRUG_DIR, "02-data")
FIG_DIR <- file.path(DRUG_DIR, "03-res", "figures", "publication_ready")
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

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.text = element_text(color = "#1F2933"),
      axis.title = element_text(color = "#0B1F3A", face = "bold"),
      plot.title = element_text(color = "#0B1F3A", face = "bold", hjust = 0.5),
      strip.background = element_rect(fill = "#F5F7FA", color = "#B8C2D6", linewidth = 0.4),
      strip.text = element_text(color = "#0B1F3A", face = "bold"),
      legend.title = element_text(face = "bold", color = "#0B1F3A"),
      legend.text = element_text(color = "#1F2933"),
      panel.grid.major.y = element_line(color = "#E8EDF3", linewidth = 0.25),
      panel.grid.minor = element_blank()
    )
}

save_plot <- function(plot, filename, width, height) {
  pdf_file <- file.path(FIG_DIR, paste0(filename, ".pdf"))
  png_file <- file.path(FIG_DIR, paste0(filename, ".png"))
  ggsave(pdf_file, plot, width = width, height = height, device = cairo_pdf)
  ggsave(png_file, plot, width = width, height = height, dpi = 450, bg = "white")
}

clean_group_result <- function(x, analysis_label) {
  x %>%
    mutate(
      Analysis = analysis_label,
      Median_High = ifelse(Group1 == "High", Median_Group1, Median_Group2),
      Median_Low = ifelse(Group1 == "Low", Median_Group1, Median_Group2),
      Mean_High = ifelse(Group1 == "High", Mean_Group1, Mean_Group2),
      Mean_Low = ifelse(Group1 == "Low", Mean_Group1, Mean_Group2),
      Delta_High_minus_Low = Median_High - Median_Low,
      High_sensitivity_effect = Median_Low - Median_High,
      Sensitivity_Direction = ifelse(
        High_sensitivity_effect > 0,
        "High group more sensitive",
        "Low group more sensitive"
      ),
      Signif = p_star(Wilcox_FDR)
    )
}

make_drug_label <- function(x) {
  dup <- duplicated(x$DrugName) | duplicated(x$DrugName, fromLast = TRUE)
  ifelse(dup, paste0(x$DrugName, " [", x$DrugID, "]"), x$DrugName)
}

pick_by_name <- function(tab, drug_names, n_fallback = 12) {
  picked <- lapply(drug_names, function(nm) {
    hit <- tab %>% filter(toupper(DrugName) == toupper(nm)) %>% arrange(Wilcox_FDR)
    if (nrow(hit) == 0) return(NULL)
    hit[1, , drop = FALSE]
  }) %>% bind_rows()

  if (nrow(picked) < n_fallback) {
    extra <- tab %>%
      filter(Wilcox_FDR < 0.05, !DrugKey %in% picked$DrugKey) %>%
      arrange(Wilcox_FDR) %>%
      slice_head(n = n_fallback - nrow(picked))
    picked <- bind_rows(picked, extra)
  }

  picked %>% distinct(DrugKey, .keep_all = TRUE) %>% slice_head(n = n_fallback)
}

make_violin_data <- function(pred, anno, selected, group_col, score_label) {
  pred %>%
    select(Sample, all_of(selected$DrugKey)) %>%
    pivot_longer(-Sample, names_to = "DrugKey", values_to = "Predicted_IC50") %>%
    left_join(anno, by = "Sample") %>%
    left_join(selected %>% select(DrugKey, DrugLabel, Wilcox_FDR, Signif), by = "DrugKey") %>%
    mutate(
      Group = factor(.data[[group_col]], levels = c("Low", "High")),
      LogIC50 = log10(Predicted_IC50 + 1),
      ScoreLabel = score_label
    ) %>%
    filter(!is.na(Group), is.finite(LogIC50))
}

make_violin_plot <- function(plot_df, selected, title) {
  ypos <- plot_df %>%
    group_by(DrugLabel) %>%
    summarise(
      y = max(LogIC50, na.rm = TRUE) + 0.08 * diff(range(LogIC50, na.rm = TRUE)) + 0.03,
      .groups = "drop"
    ) %>%
    left_join(selected %>% select(DrugLabel, Signif), by = "DrugLabel")

  ggplot(plot_df, aes(Group, LogIC50, fill = Group)) +
    geom_violin(width = 0.88, trim = TRUE, alpha = 0.78, color = NA) +
    geom_boxplot(width = 0.22, outlier.shape = NA, linewidth = 0.35, color = "#1F2933", alpha = 0.88) +
    geom_jitter(aes(color = Group), width = 0.12, size = 0.35, alpha = 0.28, show.legend = FALSE) +
    geom_segment(
      data = ypos,
      aes(x = 1, xend = 2, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.3,
      color = "#1F2933"
    ) +
    geom_text(
      data = ypos,
      aes(x = 1.5, y = y, label = Signif),
      inherit.aes = FALSE,
      size = 4.2,
      fontface = "bold",
      color = "#0B1F3A",
      vjust = -0.25
    ) +
    facet_wrap(~DrugLabel, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
    scale_color_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
    labs(
      x = NULL,
      y = "log10(predicted IC50 + 1)",
      title = title,
      fill = "Group"
    ) +
    coord_cartesian(clip = "off") +
    theme_pub(10) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(face = "bold"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

make_bubble_plot <- function(tab, title, top_n = 40) {
  plot_tab <- tab %>%
    arrange(Wilcox_FDR) %>%
    slice_head(n = top_n) %>%
    mutate(
      DrugLabel = make_drug_label(.),
      DrugLabel = factor(DrugLabel, levels = rev(DrugLabel)),
      PathwayName = factor(PathwayName, levels = unique(PathwayName[order(PathwayName)])),
      NegLogFDR = pmin(-log10(Wilcox_FDR + 1e-300), 20)
    )

  ggplot(plot_tab, aes(PathwayName, DrugLabel)) +
    geom_point(
      aes(size = NegLogFDR, fill = High_sensitivity_effect),
      shape = 21,
      color = "#243447",
      stroke = 0.25,
      alpha = 0.95
    ) +
    geom_text(aes(label = Signif), size = 2.6, fontface = "bold", color = "#111827") +
    scale_fill_gradient2(
      low = "#2F6DB3",
      mid = "white",
      high = "#C84630",
      midpoint = 0,
      name = "High-sensitivity\neffect"
    ) +
    scale_size_continuous(name = "-log10(FDR)", range = c(2.2, 8.2)) +
    labs(
      x = NULL,
      y = NULL,
      title = title
    ) +
    theme_pub(10) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      panel.grid.major.x = element_line(color = "#EEF2F7", linewidth = 0.25),
      legend.position = "right"
    )
}

make_correlation_plot <- function(cor_tab, title, top_n_each = 18) {
  plot_tab <- cor_tab %>%
    group_by(Analysis) %>%
    arrange(Spearman_FDR, .by_group = TRUE) %>%
    slice_head(n = top_n_each) %>%
    ungroup() %>%
    mutate(
      Signif = p_star(Spearman_FDR),
      DrugLabel = make_drug_label(.),
      DrugFacetLabel = paste(DrugLabel, Analysis, sep = "___"),
      DrugFacetLabel = factor(DrugFacetLabel, levels = rev(unique(DrugFacetLabel))),
      NegLogFDR = pmin(-log10(Spearman_FDR + 1e-300), 20)
    )

  label_map <- setNames(plot_tab$DrugLabel, plot_tab$DrugFacetLabel)

  ggplot(plot_tab, aes(Spearman_Rho, DrugFacetLabel)) +
    geom_vline(xintercept = 0, linetype = 2, color = "#8A94A6", linewidth = 0.35) +
    geom_segment(
      aes(x = 0, xend = Spearman_Rho, yend = DrugFacetLabel, color = PathwayName),
      linewidth = 0.75,
      alpha = 0.75
    ) +
    geom_point(aes(size = NegLogFDR, fill = Spearman_Rho), shape = 21, color = "#1F2933", stroke = 0.25) +
    geom_text(aes(label = Signif), nudge_x = 0.025, size = 2.7, fontface = "bold", color = "#111827") +
    facet_grid(Analysis ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = label_map) +
    scale_fill_gradient2(low = "#2F6DB3", mid = "white", high = "#C84630", midpoint = 0) +
    scale_size_continuous(range = c(2.2, 6.5), name = "-log10(FDR)") +
    labs(
      x = "Spearman correlation with predicted IC50",
      y = NULL,
      title = title,
      fill = "Rho",
      color = "Pathway"
    ) +
    theme_pub(9) +
    theme(
      legend.position = "right",
      strip.text.y = element_text(angle = 0)
    )
}

save_heatmap <- function(tab_mo, tab_dprs) {
  mo <- tab_mo %>%
    select(DrugKey, DrugName, DrugID, PutativeTarget, PathwayName, MO_effect = High_sensitivity_effect, MO_FDR = Wilcox_FDR)
  dp <- tab_dprs %>%
    select(DrugKey, DPRS_effect = High_sensitivity_effect, DPRS_FDR = Wilcox_FDR)

  hm <- full_join(mo, dp, by = "DrugKey") %>%
    mutate(
      MinFDR = pmin(MO_FDR, DPRS_FDR, na.rm = TRUE),
      DrugLabel = make_drug_label(.)
    ) %>%
    arrange(PathwayName, MinFDR) %>%
    filter(MinFDR < 0.05 | row_number() <= 45)

  mat <- as.matrix(hm[, c("MO_effect", "DPRS_effect")])
  rownames(mat) <- hm$DrugLabel
  colnames(mat) <- c("MO-DDRscore", "DPRS")
  mat[!is.finite(mat)] <- NA

  fdr_mat <- as.matrix(hm[, c("MO_FDR", "DPRS_FDR")])
  rownames(fdr_mat) <- hm$DrugLabel
  colnames(fdr_mat) <- colnames(mat)

  max_abs <- quantile(abs(mat), 0.96, na.rm = TRUE)
  max_abs <- max(max_abs, 0.1, na.rm = TRUE)
  col_fun <- circlize::colorRamp2(c(-max_abs, 0, max_abs), c("#2F6DB3", "white", "#C84630"))

  pathway_cols <- setNames(
    c("#1B4965", "#5FA8D3", "#7A4EAB", "#D95F02", "#2A9D8F", "#E9C46A", "#7F8C8D", "#A23E48")[seq_along(unique(hm$PathwayName))],
    unique(hm$PathwayName)
  )

  ha <- rowAnnotation(
    Pathway = hm$PathwayName,
    col = list(Pathway = pathway_cols),
    annotation_name_gp = gpar(fontface = "bold", fontsize = 9),
    gp = gpar(col = NA)
  )

  ht <- Heatmap(
    mat,
    name = "High-sensitivity\neffect",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = hm$PathwayName,
    left_annotation = ha,
    row_title_gp = gpar(fontsize = 8, fontface = "bold"),
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    rect_gp = gpar(col = "white", lwd = 0.8),
    heatmap_legend_param = list(
      title_gp = gpar(fontface = "bold", fontsize = 9),
      labels_gp = gpar(fontsize = 8),
      legend_height = unit(35, "mm")
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      fdr <- fdr_mat[i, j]
      lab <- p_star(fdr)
      if (!is.na(fdr) && fdr < 0.05) {
        grid.text(lab, x, y, gp = gpar(fontsize = 8, fontface = "bold", col = "#111827"))
      }
    }
  )

  pdf(file.path(FIG_DIR, "Drug_E_combined_sensitivity_effect_ComplexHeatmap.pdf"), width = 7.2, height = 9.2, useDingbats = FALSE)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  png(file.path(FIG_DIR, "Drug_E_combined_sensitivity_effect_ComplexHeatmap.png"), width = 3200, height = 4100, res = 450, type = "cairo")
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  save_csv(hm, file.path(FIG_DIR, "Drug_E_heatmap_plot_data.csv"))
}

make_pathway_summary_plot <- function(tab_all) {
  plot_tab <- tab_all %>%
    filter(Wilcox_FDR < 0.05) %>%
    count(Analysis, PathwayName, Sensitivity_Direction) %>%
    mutate(
      Direction = ifelse(Sensitivity_Direction == "High group more sensitive", "High more sensitive", "Low more sensitive"),
      n_signed = ifelse(Direction == "High more sensitive", n, -n)
    )

  ggplot(plot_tab, aes(n_signed, reorder(PathwayName, abs(n_signed), FUN = max), fill = Direction)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.25) +
    geom_vline(xintercept = 0, color = "#1F2933", linewidth = 0.35) +
    facet_wrap(~Analysis, ncol = 1) +
    scale_fill_manual(values = c("High more sensitive" = "#C84630", "Low more sensitive" = "#3B75AF")) +
    scale_x_continuous(labels = abs) +
    labs(
      x = "Number of significant drugs (FDR < 0.05)",
      y = NULL,
      title = "Pathway-level drug sensitivity patterns",
      fill = NULL
    ) +
    theme_pub(10) +
    theme(legend.position = "top")
}

############################################################
# Load data
############################################################

pred <- fread(file.path(TABLE_DIR, "Predicted_IC50_TCGA_all_drugs.csv"), data.table = FALSE, check.names = FALSE)
anno <- fread(file.path(DATA_DIR, "Drug_analysis_sample_annotation.csv"), data.table = FALSE, check.names = FALSE)

mo <- fread(file.path(TABLE_DIR, "Drug_sensitivity_MO_DDRscore_group_comparison.csv"), data.table = FALSE, check.names = FALSE) %>%
  clean_group_result("MO-DDRscore High vs Low")
dprs <- fread(file.path(TABLE_DIR, "Drug_sensitivity_DPRS_group_comparison.csv"), data.table = FALSE, check.names = FALSE) %>%
  clean_group_result("DPRS High vs Low")

mo_cor <- fread(file.path(TABLE_DIR, "Drug_sensitivity_MO_DDRscore_spearman.csv"), data.table = FALSE, check.names = FALSE) %>%
  mutate(Analysis = "MO-DDRscore")
dprs_cor <- fread(file.path(TABLE_DIR, "Drug_sensitivity_DPRS_spearman.csv"), data.table = FALSE, check.names = FALSE) %>%
  mutate(Analysis = "DPRS")

mo$DrugLabel <- make_drug_label(mo)
dprs$DrugLabel <- make_drug_label(dprs)

save_csv(mo, file.path(FIG_DIR, "Drug_MO_DDRscore_group_plot_ready_table.csv"))
save_csv(dprs, file.path(FIG_DIR, "Drug_DPRS_group_plot_ready_table.csv"))

############################################################
# A/B. Representative drug violin-box plots
############################################################

mo_preferred <- c(
  "Pictilisib", "AZD8055", "Temsirolimus", "Selumetinib",
  "AZD7762", "MK-1775", "Wee1 Inhibitor", "AZD6738",
  "Olaparib", "Talazoparib", "Cisplatin", "Paclitaxel"
)

dprs_preferred <- c(
  "AZD7762", "AZD6738", "MK-1775", "Wee1 Inhibitor",
  "Docetaxel", "Paclitaxel", "Cisplatin", "Gemcitabine",
  "Erlotinib", "Gefitinib", "Selumetinib", "VE-822"
)

mo_selected <- pick_by_name(mo, mo_preferred, n_fallback = 12)
dprs_selected <- pick_by_name(dprs, dprs_preferred, n_fallback = 12)

mo_violin_df <- make_violin_data(pred, anno, mo_selected, "MO_DDRscore_group", "MO-DDRscore")
dprs_violin_df <- make_violin_data(pred, anno, dprs_selected, "DPRS_RiskGroup", "DPRS")

p_mo_violin <- make_violin_plot(
  mo_violin_df,
  mo_selected,
  "Predicted drug sensitivity by MO-DDRscore group"
)
p_dprs_violin <- make_violin_plot(
  dprs_violin_df,
  dprs_selected,
  "Predicted drug sensitivity by DPRS risk group"
)

save_plot(p_mo_violin, "Drug_A_MO_DDRscore_representative_IC50_violin_box", 10.5, 8.2)
save_plot(p_dprs_violin, "Drug_B_DPRS_representative_IC50_violin_box", 10.5, 8.2)
save_csv(mo_violin_df, file.path(FIG_DIR, "Drug_A_MO_DDRscore_violin_plot_data.csv"))
save_csv(dprs_violin_df, file.path(FIG_DIR, "Drug_B_DPRS_violin_plot_data.csv"))

############################################################
# C/D. Global bubble plots
############################################################

p_mo_bubble <- make_bubble_plot(mo, "MO-DDRscore-associated predicted drug sensitivity", top_n = 40)
p_dprs_bubble <- make_bubble_plot(dprs, "DPRS-associated predicted drug sensitivity", top_n = 40)

save_plot(p_mo_bubble, "Drug_C_MO_DDRscore_global_bubble", 10.5, 9.0)
save_plot(p_dprs_bubble, "Drug_D_DPRS_global_bubble", 10.5, 9.0)

############################################################
# E. Combined ComplexHeatmap
############################################################

save_heatmap(mo, dprs)

############################################################
# F. Continuous correlation lollipop
############################################################

cor_all <- bind_rows(mo_cor, dprs_cor) %>%
  mutate(Signif = p_star(Spearman_FDR))
p_cor <- make_correlation_plot(
  cor_all,
  "Continuous-score associations with predicted IC50",
  top_n_each = 18
)
save_plot(p_cor, "Drug_F_continuous_score_IC50_correlation_lollipop", 9.2, 10.2)
save_csv(cor_all, file.path(FIG_DIR, "Drug_F_correlation_plot_ready_table.csv"))

############################################################
# G. Pathway summary
############################################################

group_all <- bind_rows(mo, dprs)
p_pathway <- make_pathway_summary_plot(group_all)
save_plot(p_pathway, "Drug_G_pathway_level_sensitivity_summary", 8.6, 6.8)
save_csv(group_all, file.path(FIG_DIR, "Drug_group_all_plot_ready_table.csv"))

summary_df <- data.frame(
  Item = c(
    "N_figures_pdf",
    "N_figures_png",
    "N_MO_DDRscore_selected_drugs",
    "N_DPRS_selected_drugs",
    "N_MO_DDRscore_FDR_lt_0.05",
    "N_DPRS_FDR_lt_0.05",
    "Figure_directory"
  ),
  Value = c(
    length(list.files(FIG_DIR, pattern = "\\.pdf$")),
    length(list.files(FIG_DIR, pattern = "\\.png$")),
    nrow(mo_selected),
    nrow(dprs_selected),
    sum(mo$Wilcox_FDR < 0.05, na.rm = TRUE),
    sum(dprs$Wilcox_FDR < 0.05, na.rm = TRUE),
    FIG_DIR
  )
)

save_csv(summary_df, file.path(FIG_DIR, "Drug_publication_plot_summary.csv"))

cat("\nDone.\n")
cat("Figure directory:\n", FIG_DIR, "\n")
print(summary_df)
