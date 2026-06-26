############################################################
# RepairDis Figure 2 final panels
# - all asterisks are generated from FDR, not raw P
# - Figure 2B keeps one representative algorithm per cell type
# - Figure 2E checkpoint heatmap is reordered by effect pattern
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
})

BASE_DIR <- "D:/R_workspace/summary/Repairium_pancancer_analysis_plan"
PLOT_DIR <- file.path(BASE_DIR, "03-res/plots/张岩-图")
FIG_DIR <- file.path(PLOT_DIR, "Figure2")
SEL_DIR <- file.path(BASE_DIR, "03-res/plots_data/main_figure_selected_data")
LOCK_DATA <- file.path(BASE_DIR, "03-res/plots/final_candidate_v1_locked/plot_data")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOCK_DATA, recursive = TRUE, showWarnings = FALSE)

pal_low <- "#20AEB3"
pal_high <- "#D45F5F"
pal_dark <- "#10243C"
pal_grid <- "#E7EDF3"
pal_grey <- "#C9D1DA"

cancer_order <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
  "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
  "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
)

sig_from_fdr <- function(fdr, ns_blank = FALSE) {
  out <- dplyr::case_when(
    is.na(fdr) ~ "",
    fdr < 1e-4 ~ "****",
    fdr < 1e-3 ~ "***",
    fdr < 1e-2 ~ "**",
    fdr < 5e-2 ~ "*",
    TRUE ~ "ns"
  )
  if (ns_blank) out[out == "ns"] <- ""
  out
}

cap <- function(x, lim = 2.6) pmax(pmin(x, lim), -lim)

order_heat_data <- function(df, row_col = "Feature_Display", column_col = "Cancer", value_col = "PlotValue") {
  df <- df %>%
    mutate(
      "{row_col}" := as.character(.data[[row_col]]),
      "{column_col}" := as.character(.data[[column_col]])
    )

  mat_df <- df %>%
    select(all_of(c(row_col, column_col, value_col))) %>%
    group_by(.data[[row_col]], .data[[column_col]]) %>%
    summarise(Value = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = all_of(column_col), values_from = Value, values_fill = 0)

  mat <- as.matrix(mat_df[, setdiff(colnames(mat_df), row_col), drop = FALSE])
  rownames(mat) <- mat_df[[row_col]]

  row_order <- rownames(mat)
  col_order <- colnames(mat)
  if (nrow(mat) > 1) row_order <- rownames(mat)[hclust(dist(mat))$order]
  if (ncol(mat) > 1) col_order <- colnames(mat)[hclust(dist(t(mat)))$order]

  df %>%
    mutate(
      "{row_col}" := factor(.data[[row_col]], levels = rev(row_order)),
      "{column_col}" := factor(.data[[column_col]], levels = col_order)
    )
}

save_plot <- function(p, stem, width, height) {
  pdf_file <- file.path(FIG_DIR, paste0(stem, ".pdf"))
  png_file <- file.path(FIG_DIR, paste0(stem, ".png"))
  tryCatch(
    ggsave(pdf_file, p, width = width, height = height, useDingbats = FALSE),
    error = function(e) {
      fallback <- file.path(FIG_DIR, paste0(stem, "_FDR_FIXED.pdf"))
      message("PDF locked, writing fallback: ", fallback)
      ggsave(fallback, p, width = width, height = height, useDingbats = FALSE)
    }
  )
  tryCatch(
    ggsave(png_file, p, width = width, height = height, dpi = 450, bg = "white"),
    error = function(e) {
      fallback <- file.path(FIG_DIR, paste0(stem, "_FDR_FIXED.png"))
      message("PNG locked, writing fallback: ", fallback)
      ggsave(fallback, p, width = width, height = height, dpi = 450, bg = "white")
    }
  )
}

theme_repair <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(color = pal_dark),
      plot.title = element_text(face = "bold", size = base_size + 5, hjust = 0),
      plot.subtitle = element_text(size = base_size + 0.5, color = "#617083", hjust = 0),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = pal_grid, linewidth = 0.28),
      panel.border = element_rect(color = "#CAD5E0", fill = NA, linewidth = 0.45),
      axis.text = element_text(color = "#263849"),
      axis.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "#EEF4F8", color = "#D7E1EA", linewidth = 0.45),
      strip.text = element_text(face = "bold", color = pal_dark),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(8, 12, 8, 8)
    )
}

heat_panel <- function(df, title, subtitle, stem, width, height, lim = 2.6) {
  p <- ggplot(df, aes(Cancer, Feature_Display, fill = PlotValue)) +
    geom_tile(color = "white", linewidth = 0.55, width = 0.94, height = 0.90) +
    geom_text(aes(label = Significance), size = 3.45, fontface = "bold", color = pal_dark) +
    scale_fill_gradient2(low = pal_low, mid = "white", high = pal_high,
                         midpoint = 0, limits = c(-lim, lim),
                         oob = scales::squish, name = "High-Low\nscaled effect") +
    labs(title = title, subtitle = subtitle, x = "Cancer type", y = NULL) +
    theme_repair(9.3) +
    theme(
      axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 8.0, face = "bold"),
      axis.text.y = element_text(size = 8.4, face = "bold"),
      legend.key.height = unit(18, "pt")
    )
  save_plot(p, stem, width, height)
  invisible(p)
}

############################################################
# Figure 2A. ESTIMATE TME profiles
############################################################

estimate <- fread(file.path(SEL_DIR, "Figure2A_ESTIMATE_selected_stats.csv"), data.table = FALSE) %>%
  mutate(
    Cancer = factor(Cancer, levels = cancer_order),
    Feature_Display = dplyr::case_when(
      grepl("ESTIMATE", Feature_Display, ignore.case = TRUE) ~ "ESTIMATEScore",
      grepl("Immune", Feature_Display, ignore.case = TRUE) ~ "ImmuneScore",
      grepl("Stromal", Feature_Display, ignore.case = TRUE) ~ "StromalScore",
      grepl("Purity", Feature_Display, ignore.case = TRUE) ~ "TumorPurity",
      TRUE ~ Feature_Display
    ),
    Feature_Display = factor(Feature_Display, levels = rev(c("ESTIMATEScore", "ImmuneScore", "StromalScore", "TumorPurity"))),
    FDR = as.numeric(FDR),
    PlotValue = cap(as.numeric(EffectScaledWithinFeature), 2.5),
    Significance = sig_from_fdr(FDR, ns_blank = TRUE)
  ) %>%
  filter(!is.na(Cancer), !is.na(Feature_Display), is.finite(PlotValue)) %>%
  order_heat_data()

fwrite(estimate, file.path(LOCK_DATA, "Figure2A_ESTIMATE_FDR_plot_data.csv"))
heat_panel(
  estimate,
  "A  ESTIMATE-based tumor microenvironment profiles",
  "MO-DDRscore-associated high-low differences in representative cancer types",
  "Figure2A_ESTIMATE",
  7.6, 3.8, 2.5
)

############################################################
# Figure 2B. TME cell type remodeling, one method per cell type
############################################################

canonical <- data.frame(
  Feature = c(
    "Neutrophils_MCPcounter",
    "NK_cells_MCPcounter",
    "Macrophages_M2_quantiseq",
    "CD4_Tcells_EPIC",
    "Dendritic_cells_quantiseq",
    "Plasma_cells_xCell"
  ),
  Feature_Label = c(
    "Neutrophils\n(MCP-counter)",
    "NK cells\n(MCP-counter)",
    "Macrophages M2\n(quanTIseq)",
    "CD4 T cells\n(EPIC)",
    "Dendritic cells\n(quanTIseq)",
    "Plasma cells\n(xCell)"
  ),
  stringsAsFactors = FALSE
)

cell <- fread(file.path(SEL_DIR, "Figure2B_immune_cell_selected_stats.csv"), data.table = FALSE) %>%
  inner_join(canonical, by = "Feature") %>%
  mutate(
    Cancer = factor(Cancer, levels = cancer_order),
    Feature_Label = factor(Feature_Label, levels = canonical$Feature_Label),
    LowValue = as.numeric(Median_Low),
    HighValue = as.numeric(Median_High),
    FDR = as.numeric(FDR),
    Significance = sig_from_fdr(FDR, ns_blank = TRUE),
    DirectionCol = ifelse(HighValue >= LowValue, "Higher in High", "Higher in Low")
  ) %>%
  filter(!is.na(Cancer), is.finite(LowValue), is.finite(HighValue))

cell_scale <- cell %>%
  group_by(Feature_Label) %>%
  summarise(mu = mean(c(LowValue, HighValue), na.rm = TRUE),
            sdv = stats::sd(c(LowValue, HighValue), na.rm = TRUE),
            .groups = "drop")

cell <- cell %>%
  left_join(cell_scale, by = "Feature_Label") %>%
  mutate(
    z_low = ifelse(is.finite(sdv) & sdv > 0, (LowValue - mu) / sdv, 0),
    z_high = ifelse(is.finite(sdv) & sdv > 0, (HighValue - mu) / sdv, 0)
  )

fwrite(cell, file.path(LOCK_DATA, "Figure2B_TME_cell_type_one_method_FDR_plot_data.csv"))

p2b <- ggplot(cell, aes(y = Cancer)) +
  geom_segment(aes(x = z_low, xend = z_high, yend = Cancer, color = DirectionCol),
               linewidth = 0.72, alpha = 0.78) +
  geom_point(aes(x = z_low), color = pal_low, fill = pal_low, shape = 21,
             size = 2.0, stroke = 0.25, alpha = 0.96) +
  geom_point(aes(x = z_high), color = pal_high, fill = pal_high, shape = 21,
             size = 2.0, stroke = 0.25, alpha = 0.96) +
  geom_vline(xintercept = 0, linetype = 2, color = "#AAB5C3", linewidth = 0.35) +
  geom_text(aes(x = 2.52, label = Significance), size = 3.35, fontface = "bold", color = pal_dark) +
  facet_wrap(~ Feature_Label, ncol = 3, scales = "free_x") +
  scale_color_manual(values = c("Higher in High" = pal_high, "Higher in Low" = pal_low), name = NULL) +
  coord_cartesian(xlim = c(-2.65, 2.82), clip = "off") +
  labs(
    title = "B  TME cell-type remodeling",
    subtitle = "One representative deconvolution output was retained for each cell type",
    x = "Feature-wise standardized median score",
    y = NULL
  ) +
  theme_repair(9.5) +
  theme(
    axis.text.y = element_text(size = 7.5),
    strip.text = element_text(size = 9.2),
    legend.position = "top"
  )

save_plot(p2b, "Figure2B_TME_cell_type", 10.8, 5.7)
save_plot(p2b, "Figure2B_TME_cell_type_one_method", 10.8, 5.7)

############################################################
# Figure 2C/D. Immune suppression and exclusion
############################################################

barrier <- fread(file.path(SEL_DIR, "Figure2C_immune_suppression_exclusion_selected_stats.csv"), data.table = FALSE) %>%
  mutate(
    Cancer = factor(Cancer, levels = cancer_order),
    FDR = as.numeric(FDR),
    Significance = sig_from_fdr(FDR, ns_blank = TRUE),
    PlotValue = cap(as.numeric(EffectScaledWithinFeature), 2.6)
  ) %>%
  filter(!is.na(Cancer), is.finite(PlotValue))

supp <- barrier %>%
  filter(Barrier_Axis == "Immune suppression") %>%
  order_heat_data()
fwrite(supp, file.path(LOCK_DATA, "Figure2C_immune_suppression_FDR_plot_data.csv"))
heat_panel(
  supp,
  "C  Immune-suppressive programs",
  "TIDE/IOBR-derived suppressive features associated with MO-DDRscore groups",
  "Figure2C_immune_suppression",
  7.2, 3.6, 2.6
)

excl <- barrier %>%
  filter(Barrier_Axis == "Immune exclusion") %>%
  order_heat_data()
fwrite(excl, file.path(LOCK_DATA, "Figure2D_immune_exclusion_FDR_plot_data.csv"))
heat_panel(
  excl,
  "D  Immune-exclusion programs",
  "TIDE/IOBR-derived exclusion features associated with MO-DDRscore groups",
  "Figure2D_immune_exclusion",
  7.2, 3.6, 2.6
)

############################################################
# Figure 2E. Checkpoint biomarker landscape, reordered
############################################################

checkpoint <- fread(file.path(SEL_DIR, "Figure2D_checkpoint_selected_stats.csv"), data.table = FALSE) %>%
  mutate(
    Cancer = as.character(Cancer),
    Feature_Display = dplyr::recode(as.character(Feature_Display), "PDCD1lg2" = "PDCD1LG2"),
    FDR = as.numeric(FDR),
    Significance = sig_from_fdr(FDR, ns_blank = TRUE),
    PlotValue = cap(as.numeric(EffectScaledWithinFeature), 2.6)
  ) %>%
  filter(is.finite(PlotValue), Cancer %in% cancer_order)

mat <- checkpoint %>%
  select(Feature_Display, Cancer, PlotValue) %>%
  tidyr::pivot_wider(names_from = Cancer, values_from = PlotValue, values_fill = 0)
mat_m <- as.matrix(mat[, -1, drop = FALSE])
rownames(mat_m) <- mat$Feature_Display

row_order <- rownames(mat_m)
col_order <- colnames(mat_m)
if (nrow(mat_m) > 1) row_order <- rownames(mat_m)[hclust(dist(mat_m))$order]
if (ncol(mat_m) > 1) col_order <- colnames(mat_m)[hclust(dist(t(mat_m)))$order]

checkpoint <- checkpoint %>%
  mutate(
    Feature_Display = factor(Feature_Display, levels = rev(row_order)),
    Cancer = factor(Cancer, levels = col_order)
  )
fwrite(checkpoint, file.path(LOCK_DATA, "Figure2E_checkpoint_reordered_FDR_plot_data.csv"))

heat_panel(
  checkpoint,
  "E  Checkpoint biomarker landscape",
  "Checkpoint expression differences reordered by high-low effect pattern",
  "Figure2E_checkpoint_biomarker",
  7.8, 4.4, 2.6
)

############################################################
# Figure 2F. TIDE predicted response, responder in red bottom
############################################################

tide_resp <- fread(file.path(SEL_DIR, "Figure2E_TIDE_selected_response_summary.csv"), data.table = FALSE) %>%
  mutate(
    Cancer = as.character(Cancer),
    Group = factor(Group, levels = c("Low", "High")),
    StatusLabel = factor(StatusLabel, levels = c("Non-responder", "Predicted responder")),
    Percent = as.numeric(Percent)
  )

tide_or <- fread(file.path(SEL_DIR, "Figure2E_TIDE_selected_OR_stats.csv"), data.table = FALSE) %>%
  mutate(
    Cancer = as.character(Cancer),
    FDR = as.numeric(FDR),
    Significance = sig_from_fdr(FDR, ns_blank = TRUE),
    Delta = as.numeric(High_Responder_Percent) - as.numeric(Low_Responder_Percent)
  )

tide_order <- tide_or %>% arrange(desc(Delta)) %>% pull(Cancer)
tide_resp <- tide_resp %>% mutate(Cancer = factor(Cancer, levels = tide_order))
tide_or <- tide_or %>% mutate(Cancer = factor(Cancer, levels = tide_order))

ann <- tide_or %>%
  transmute(Cancer, Group = factor("High", levels = c("Low", "High")), y = 103, Significance)

fwrite(tide_resp, file.path(LOCK_DATA, "Figure2F_TIDE_response_FDR_plot_data.csv"))
fwrite(tide_or, file.path(LOCK_DATA, "Figure2F_TIDE_response_OR_FDR_stats.csv"))

p2f <- ggplot(tide_resp, aes(x = Cancer, y = Percent, fill = StatusLabel)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.30) +
  geom_text(data = ann, aes(x = Cancer, y = y, label = Significance),
            inherit.aes = FALSE, color = pal_dark, fontface = "bold", size = 3.45) +
  facet_grid(. ~ Group) +
  scale_fill_manual(
    values = c("Predicted responder" = pal_high, "Non-responder" = pal_grey),
    breaks = c("Predicted responder", "Non-responder"),
    name = NULL
  ) +
  scale_y_continuous(limits = c(0, 108), breaks = c(0, 25, 50, 75, 100),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title = "F  TIDE-predicted immunotherapy response",
    subtitle = "Predicted responder fraction in MO-DDRscore-low and -high tumors",
    x = "Cancer type",
    y = "Predicted response proportion"
  ) +
  theme_repair(9.5) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 8.2, face = "bold"),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

save_plot(p2f, "Figure2F_TIDE_response", 7.6, 4.6)

############################################################
# Figure 2G. Drug/pathway sankey, keep current final panel
############################################################

src_candidates <- c(
  file.path(BASE_DIR, "03-res/plots/final_candidate_v1_locked/Figure2/Figure2F_DDR_gene_pathway_drug_sankey_FINAL_LOCKED.pdf"),
  file.path(BASE_DIR, "03-res/plots/final_candidate_v1_locked/Figure2/Figure2F_DDR_gene_pathway_drug_sankey_FINAL_LOCKED.tiff"),
  file.path(BASE_DIR, "03-res/drug_network/figures/DDR_gene_pathway_drug_sankey_CTM_style.pdf"),
  file.path(FIG_DIR, "Figure2G_DDR_gene_pathway_drug_sankey.pdf")
)
src <- src_candidates[file.exists(src_candidates)][1]
if (!is.na(src) && file.exists(src)) {
  ext <- tools::file_ext(src)
  file.copy(src, file.path(FIG_DIR, paste0("Figure2G_DDR_gene_pathway_drug_sankey.", ext)), overwrite = TRUE)
}

fwrite(
  data.frame(
    Panel = c("Figure2A", "Figure2B", "Figure2C", "Figure2D", "Figure2E", "Figure2F"),
    Significance_rule = "Asterisks were recalculated from FDR: * <0.05, ** <0.01, *** <0.001, **** <0.0001.",
    stringsAsFactors = FALSE
  ),
  file.path(FIG_DIR, "Figure2_significance_rule.csv")
)

cat("Figure 2 final panels regenerated in:\n", FIG_DIR, "\n")
