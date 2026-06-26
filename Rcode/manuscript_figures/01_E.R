############################################################
# RepairDis Figure 1E: aging hallmark landscape
# Data source:
#   03-res/Figure2_pan_cancer/D_Aging_hallmark/<Cancer>/
#   *_Aging_hallmark_group_comparison.csv
#
# Plot:
#   1) Main figure: same style as the accepted 2026-06-14 panel
#      31 cancers, excluding MESO and UVM to match the locked main figure.
#   2) Supplementary/full version: all 33 cancers.
#
# Significance labels are based on FDR only.
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

BASE_DIR <- "D:/R_workspace/summary/Repairium_pancancer_analysis_plan"
RAW_DIR <- file.path(BASE_DIR, "03-res", "Figure2_pan_cancer", "D_Aging_hallmark")
PLOT_DIR <- file.path(BASE_DIR, "03-res", "plots", "\u5f20\u5ca9-\u56fe")
FIG_DIR <- file.path(PLOT_DIR, "Figure1")
DATA_OUT_DIR <- file.path(FIG_DIR, "plot_data")
SUP_DIR <- file.path(PLOT_DIR, "sup", "Figure1")
ROOT_PLOT_DIR <- file.path(BASE_DIR, "03-res", "plots")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUP_DIR, recursive = TRUE, showWarnings = FALSE)

cancer_order_all <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
  "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
  "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
)

cancer_order_main <- cancer_order_all

feature_order <- c(
  "Genomic instability",
  "Telomere attrition",
  "Epigenetic alterations",
  "Loss of proteostasis",
  "Disabled macroautophagy",
  "Deregulated nutrient sensing",
  "Mitochondrial dysfunction",
  "Cellular senescence",
  "Stem cell exhaustion",
  "Altered intercellular communication",
  "Chronic inflammation",
  "Extracellular matrix changes"
)

category_map <- c(
  "Genomic instability" = "Genome maintenance",
  "Telomere attrition" = "Genome maintenance",
  "Epigenetic alterations" = "Genome maintenance",
  "Loss of proteostasis" = "Proteostasis / metabolism",
  "Disabled macroautophagy" = "Proteostasis / metabolism",
  "Deregulated nutrient sensing" = "Proteostasis / metabolism",
  "Mitochondrial dysfunction" = "Proteostasis / metabolism",
  "Cellular senescence" = "Cell fate",
  "Stem cell exhaustion" = "Cell fate",
  "Altered intercellular communication" = "Communication / inflammation",
  "Chronic inflammation" = "Communication / inflammation",
  "Extracellular matrix changes" = "Communication / inflammation"
)

category_colors <- c(
  "Genome maintenance" = "#C73A4A",
  "Proteostasis / metabolism" = "#1F9E91",
  "Cell fate" = "#E6BE4F",
  "Communication / inflammation" = "#7455E8"
)

category_short <- c(
  "Genome maintenance" = "Genome maintenance",
  "Proteostasis / metabolism" = "Proteostasis / metabolism",
  "Cell fate" = "Cell fate",
  "Communication / inflammation" = "communication / inflammation"
)

pal_low <- "#3262AA"
pal_mid <- "#F3EFE8"
pal_high <- "#CA4253"
pal_dark <- "#183047"
pal_grid <- "#E6ECF2"

sig_star <- function(fdr) {
  dplyr::case_when(
    is.na(fdr) ~ "",
    fdr < 1e-3 ~ "***",
    fdr < 1e-2 ~ "**",
    fdr < 5e-2 ~ "*",
    TRUE ~ ""
  )
}

cap <- function(x, lim = 2.6) {
  pmax(pmin(as.numeric(x), lim), -lim)
}

save_all <- function(plot, out_dir, stem, width, height) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_file <- file.path(out_dir, paste0(stem, ".pdf"))
  png_file <- file.path(out_dir, paste0(stem, ".png"))
  tiff_file <- file.path(out_dir, paste0(stem, ".tiff"))
  tryCatch(
    ggsave(pdf_file, plot, width = width, height = height, useDingbats = FALSE, bg = "white"),
    error = function(e) {
      ggsave(file.path(out_dir, paste0(stem, "_UPDATED.pdf")), plot,
             width = width, height = height, useDingbats = FALSE, bg = "white")
    }
  )
  ggsave(png_file, plot, width = width, height = height, dpi = 520, bg = "white")
  ggsave(tiff_file, plot, width = width, height = height, dpi = 600, bg = "white", compression = "lzw")
  invisible(pdf_file)
}

comparison_files <- Sys.glob(file.path(RAW_DIR, "*", "*_Aging_hallmark_group_comparison.csv"))
if (length(comparison_files) == 0) {
  stop("No aging hallmark group-comparison files found under: ", RAW_DIR)
}

aging_all <- rbindlist(lapply(comparison_files, fread), fill = TRUE) %>%
  as.data.frame() %>%
  mutate(
    Cancer = as.character(Cancer),
    Feature = as.character(Feature),
    Effect = as.numeric(Effect),
    P = as.numeric(P),
    FDR = as.numeric(FDR),
    Category = unname(category_map[Feature]),
    Row = length(feature_order) - match(Feature, feature_order) + 1,
    EffectPlot = cap(Effect, 2.6),
    SigLabel = sig_star(FDR)
  ) %>%
  filter(
    Cancer %in% cancer_order_all,
    Feature %in% feature_order,
    is.finite(EffectPlot),
    !is.na(Row)
  )

fwrite(
  aging_all %>%
    arrange(match(Feature, feature_order), match(Cancer, cancer_order_all)),
  file.path(DATA_OUT_DIR, "Figure1E_aging_hallmark_landscape_all33_plot_data.csv")
)

feature_rows <- data.frame(
  Feature = feature_order,
  Row = length(feature_order) - seq_along(feature_order) + 1,
  Category = unname(category_map[feature_order]),
  stringsAsFactors = FALSE
)

draw_aging_landscape <- function(plot_data, cancer_order, title_suffix = "", full33 = FALSE) {
  order_data <- plot_data %>%
    filter(Cancer %in% cancer_order, Feature %in% feature_order)

  order_mat <- matrix(
    NA_real_,
    nrow = length(feature_order),
    ncol = length(cancer_order),
    dimnames = list(feature_order, cancer_order)
  )
  order_i <- match(order_data$Feature, feature_order)
  order_j <- match(order_data$Cancer, cancer_order)
  order_mat[cbind(order_i, order_j)] <- order_data$EffectPlot
  order_mat[!is.finite(order_mat)] <- 0

  color_order_data <- plot_data %>%
    filter(Cancer %in% cancer_order, Feature %in% feature_order) %>%
    mutate(W = ifelse(FDR < 0.05, 1, 0.25))

  cancer_order_panel <- color_order_data %>%
    group_by(Cancer) %>%
    summarise(ColorScore = sum(EffectPlot * W, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(ColorScore)) %>%
    pull(Cancer)

  feature_order_panel <- color_order_data %>%
    group_by(Feature) %>%
    summarise(ColorScore = sum(EffectPlot * W, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(ColorScore)) %>%
    pull(Feature)

  feature_rows_panel <- data.frame(
    Feature = feature_order_panel,
    Row = length(feature_order_panel) - seq_along(feature_order_panel) + 1,
    stringsAsFactors = FALSE
  )

  plot_data <- plot_data %>%
    filter(Cancer %in% cancer_order) %>%
    mutate(
      Cancer = factor(Cancer, levels = cancer_order_panel),
      Feature = factor(Feature, levels = feature_order_panel),
      Row = length(feature_order_panel) - match(as.character(Feature), feature_order_panel) + 1
    )

  count_df <- plot_data %>%
    group_by(Feature, Row) %>%
    summarise(
      Low_n = sum(FDR < 0.05 & Effect < 0, na.rm = TRUE),
      High_n = sum(FDR < 0.05 & Effect > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Low_plot = -Low_n,
      Feature = factor(Feature, levels = feature_order_panel)
    )

  x_count_lim <- max(c(count_df$Low_n, count_df$High_n), na.rm = TRUE)
  if (!is.finite(x_count_lim) || x_count_lim < 1) x_count_lim <- 1
  x_count_lim <- max(8, ceiling(x_count_lim / 5) * 5)

  p_heat <- ggplot(plot_data, aes(x = Cancer, y = Row, fill = EffectPlot)) +
    geom_tile(color = "white", linewidth = 0.42, width = 0.93, height = 0.88) +
    geom_text(aes(label = SigLabel), color = pal_dark, fontface = "bold", size = ifelse(full33, 3.15, 3.45)) +
    scale_y_continuous(
      breaks = feature_rows_panel$Row,
      labels = feature_rows_panel$Feature,
      limits = c(0.5, length(feature_order_panel) + 0.5),
      expand = c(0, 0)
    ) +
    scale_fill_gradient2(
      low = pal_low,
      mid = pal_mid,
      high = pal_high,
      midpoint = 0,
      limits = c(-2.6, 2.6),
      breaks = c(-2, -1, 0, 1, 2),
      oob = scales::squish,
      name = "High - Low\nactivity"
    ) +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      text = element_text(color = pal_dark),
      axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = ifelse(full33, 6.4, 7.2), face = "bold"),
      axis.text.y = element_text(size = ifelse(full33, 7.4, 8.2), face = "bold", color = "#2A3E50"),
      axis.ticks = element_blank(),
      axis.line.x = element_line(color = "#2A3E50", linewidth = 0.45),
      axis.line.y = element_line(color = "#2A3E50", linewidth = 0.45),
      panel.grid = element_blank(),
      plot.margin = margin(40, 4, 28, 0),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 8.2),
      legend.text = element_text(size = 7.5),
      legend.key.height = unit(20, "pt")
    )

  p_count <- ggplot(count_df, aes(y = Row)) +
    geom_vline(xintercept = 0, color = "#68798A", linewidth = 0.42) +
    geom_segment(
      data = count_df %>% filter(Low_n > 0),
      aes(x = 0, xend = Low_plot, yend = Row),
      color = pal_low,
      linewidth = ifelse(full33, 5.8, 6.4),
      lineend = "butt",
      alpha = 0.98
    ) +
    geom_segment(
      data = count_df %>% filter(High_n > 0),
      aes(x = 0, xend = High_n, yend = Row),
      color = "#C93E50",
      linewidth = ifelse(full33, 5.8, 6.4),
      lineend = "butt",
      alpha = 0.98
    ) +
    geom_text(
      data = count_df %>% filter(Low_n > 0),
      aes(x = Low_plot - 1.60, label = Low_n),
      color = pal_low,
      fontface = "bold",
      size = 2.35
    ) +
    geom_text(
      data = count_df %>% filter(High_n > 0),
      aes(x = High_n + 1.60, label = High_n),
      color = "#C93E50",
      fontface = "bold",
      size = 2.35
    ) +
    annotate("text", x = -x_count_lim * 0.55, y = length(feature_order) + 0.65,
             label = "No. in low", color = pal_low, fontface = "bold", size = 3.0) +
    annotate("text", x = x_count_lim * 0.55, y = length(feature_order) + 0.65,
             label = "No. in high", color = "#C93E50", fontface = "bold", size = 3.0) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(
      xlim = c(-x_count_lim - 0.1, x_count_lim + 0.1),
      ylim = c(0.5, length(feature_order_panel) + 0.5),
      clip = "off"
    ) +
    theme_void() +
    theme(plot.margin = margin(40, 18, 28, 0))

  combined <- p_heat + p_count +
    plot_layout(widths = c(0.82, 0.18), guides = "collect")

  combined +
    plot_annotation(
      title = paste0("E  Aging hallmark landscape", title_suffix),
      theme = theme(
        plot.title = element_text(face = "bold", size = 21, color = pal_dark, hjust = 0.02),
        plot.subtitle = element_text(size = 10.2, color = "#718093", hjust = 0.02),
        plot.margin = margin(4, 8, 4, 4)
      )
    )
}

main_data <- aging_all %>%
  filter(Cancer %in% cancer_order_main)

fwrite(
  main_data %>% arrange(match(Feature, feature_order), match(Cancer, cancer_order_main)),
  file.path(DATA_OUT_DIR, "Figure1E_aging_hallmark_landscape_plot_data.csv")
)

main_counts <- main_data %>%
  group_by(Feature) %>%
  summarise(
    Low_cancers = sum(FDR < 0.05 & Effect < 0, na.rm = TRUE),
    High_cancers = sum(FDR < 0.05 & Effect > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(match(Feature, feature_order))

fwrite(main_counts, file.path(DATA_OUT_DIR, "Figure1E_aging_hallmark_direction_counts.csv"))

main_p <- draw_aging_landscape(main_data, cancer_order_main, title_suffix = "", full33 = FALSE)
save_all(main_p, FIG_DIR, "Figure1E_aging_hallmark_landscape", width = 17.2, height = 6.2)

full_p <- draw_aging_landscape(aging_all, cancer_order_all, title_suffix = " (all 33 cancer types)", full33 = TRUE)
save_all(full_p, SUP_DIR, "FigureS1E_aging_hallmark_landscape_all33", width = 17.2, height = 6.2)

# Keep a root-level copy because older figure-collection scripts expect this name.
for (ext in c("pdf", "png", "tiff")) {
  src <- file.path(FIG_DIR, paste0("Figure1E_aging_hallmark_landscape.", ext))
  if (file.exists(src)) {
    file.copy(src, file.path(ROOT_PLOT_DIR, paste0("Figure1E_aging_hallmark_landscape.", ext)), overwrite = TRUE)
  }
}

summary_df <- data.frame(
  Item = c(
    "Data_source",
    "N_cancers_all",
    "N_cancers_main",
    "Excluded_from_main",
    "N_aging_hallmarks",
    "N_tiles_main",
    "N_main_FDR_lt_0.05",
    "Significance_rule"
  ),
  Value = c(
    RAW_DIR,
    length(unique(aging_all$Cancer)),
    length(unique(main_data$Cancer)),
    "None",
    length(unique(aging_all$Feature)),
    nrow(main_data),
    sum(main_data$FDR < 0.05, na.rm = TRUE),
    "FDR based: * < 0.05, ** < 0.01, *** < 0.001"
  )
)
fwrite(summary_df, file.path(DATA_OUT_DIR, "Figure1E_aging_hallmark_landscape_run_summary.csv"))

cat("\nDone.\n")
cat("Main output:\n", file.path(FIG_DIR, "Figure1E_aging_hallmark_landscape.pdf"), "\n")
cat("Full 33-cancer supplementary output:\n", file.path(SUP_DIR, "FigureS1E_aging_hallmark_landscape_all33.pdf"), "\n")
print(summary_df)
