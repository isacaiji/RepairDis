############################################################
# Figure 3C-D candidate: immune program heatmaps
# - Main panels: compact heatmaps using High - Low standardized score
# - Significance stars: FDR-based only
# - Supplementary panels: paired group medians from the same result table
# NOTE: This script does not recompute statistics.
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

BASE_DIR <- "D:/R_workspace/summary/Repairium_pancancer_analysis_plan"
PLOT_DIR <- file.path(BASE_DIR, "03-res", "plots", "张岩-图")
DATA_DIR <- file.path(PLOT_DIR, "plot_data")
FIG_DIR  <- file.path(PLOT_DIR, "Figure3")
SUP_DIR  <- file.path(PLOT_DIR, "sup", "Figure3")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUP_DIR, recursive = TRUE, showWarnings = FALSE)

pal_low  <- "#24B7B6"
pal_high <- "#D85B68"
pal_mid  <- "#F7F3EE"
pal_dark <- "#10243C"
pal_grid <- "#E6EDF3"

fdr_star <- function(x) {
  dplyr::case_when(
    is.na(x) ~ "",
    x < 0.0001 ~ "****",
    x < 0.001  ~ "***",
    x < 0.01   ~ "**",
    x < 0.05   ~ "*",
    TRUE ~ ""
  )
}

squish_range <- function(x, lower = -2.5, upper = 2.5) {
  pmin(pmax(x, lower), upper)
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

read_result_table <- function(file) {
  x <- data.table::fread(file, data.table = FALSE, check.names = FALSE)
  x %>%
    mutate(
      Cancer = as.character(Cancer),
      Feature_Display = as.character(Feature_Display),
      FDR = safe_numeric(FDR),
      Diff_Z = safe_numeric(Diff_Z),
      Low_Z = safe_numeric(Low_Z),
      High_Z = safe_numeric(High_Z),
      Median_Low = safe_numeric(Median_Low),
      Median_High = safe_numeric(Median_High),
      FDR_star = fdr_star(FDR),
      Diff_Z_plot = squish_range(Diff_Z)
    ) %>%
    filter(
      !is.na(Cancer), Cancer != "",
      !is.na(Feature_Display), Feature_Display != "",
      is.finite(Diff_Z)
    )
}

order_for_heatmap <- function(df) {
  cancer_order <- df %>%
    group_by(Cancer) %>%
    summarise(
      MeanDiff = mean(Diff_Z, na.rm = TRUE),
      SigN = sum(FDR < 0.05, na.rm = TRUE),
      HighN = sum(Diff_Z > 0 & FDR < 0.05, na.rm = TRUE),
      LowN = sum(Diff_Z < 0 & FDR < 0.05, na.rm = TRUE),
      AbsSignal = mean(abs(Diff_Z), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(DirectionBalance = HighN - LowN) %>%
    arrange(desc(DirectionBalance), desc(MeanDiff), desc(SigN), desc(AbsSignal), Cancer) %>%
    pull(Cancer)

  feature_order <- df %>%
    group_by(Feature_Display) %>%
    summarise(
      MeanDiff = mean(Diff_Z, na.rm = TRUE),
      SigN = sum(FDR < 0.05, na.rm = TRUE),
      HighN = sum(Diff_Z > 0 & FDR < 0.05, na.rm = TRUE),
      LowN = sum(Diff_Z < 0 & FDR < 0.05, na.rm = TRUE),
      AbsSignal = mean(abs(Diff_Z), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(DirectionBalance = HighN - LowN) %>%
    arrange(desc(DirectionBalance), desc(MeanDiff), desc(SigN), desc(AbsSignal), Feature_Display) %>%
    pull(Feature_Display)

  list(cancer = cancer_order, feature = rev(feature_order))
}

make_heatmap_panel <- function(file, stem, panel, title, subtitle = NULL,
                               width = 12.4, height = 5.5) {
  df <- read_result_table(file)
  ord <- order_for_heatmap(df)

  plot_df <- df %>%
    mutate(
      Cancer = factor(Cancer, levels = ord$cancer),
      Feature_Display = factor(Feature_Display, levels = ord$feature)
    ) %>%
    arrange(Feature_Display, Cancer)

  data.table::fwrite(
    plot_df,
    file.path(DATA_DIR, paste0(stem, "_plot_data.csv"))
  )

  p <- ggplot(plot_df, aes(x = Cancer, y = Feature_Display, fill = Diff_Z_plot)) +
    geom_tile(color = "white", linewidth = 1.15) +
    geom_text(
      aes(label = FDR_star),
      color = pal_dark,
      fontface = "bold",
      size = 3.45,
      vjust = 0.55
    ) +
    scale_fill_gradient2(
      low = pal_low,
      mid = pal_mid,
      high = pal_high,
      midpoint = 0,
      limits = c(-2.5, 2.5),
      breaks = c(-2, -1, 0, 1, 2),
      name = "High - Low\nstandardized\nscore"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      title = paste0(panel, "  ", title),
      subtitle = subtitle,
      x = "Cancer type",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(color = pal_dark, face = "bold", size = 24, hjust = 0),
      plot.subtitle = element_text(color = "#718093", size = 11, hjust = 0),
      axis.text.x = element_text(color = pal_dark, face = "bold", size = 10.5, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(color = pal_dark, face = "bold", size = 11),
      axis.title.x = element_text(color = pal_dark, face = "bold", size = 12, margin = margin(t = 8)),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "#CBD8E3", fill = NA, linewidth = 0.7),
      legend.title = element_text(color = pal_dark, face = "bold", size = 11),
      legend.text = element_text(color = pal_dark, size = 10),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.45, "cm"),
      plot.margin = margin(12, 18, 12, 12)
    )

  ggsave(file.path(FIG_DIR, paste0(stem, ".png")), p, width = width, height = height, dpi = 600, bg = "white")
  ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), p, width = width, height = height, bg = "white", useDingbats = FALSE)
  ggsave(file.path(FIG_DIR, paste0(stem, ".tiff")), p, width = width, height = height, dpi = 600, bg = "white", compression = "lzw")

  invisible(plot_df)
}

make_median_supp_panel <- function(file, stem, panel, title, width = 14, height = 7.5) {
  df <- read_result_table(file)
  ord <- order_for_heatmap(df)

  long_df <- df %>%
    transmute(
      Cancer = factor(Cancer, levels = ord$cancer),
      Feature_Display = factor(Feature_Display, levels = ord$feature),
      FDR,
      FDR_star,
      Low = Median_Low,
      High = Median_High
    ) %>%
    pivot_longer(cols = c(Low, High), names_to = "Group", values_to = "Median_Value") %>%
    mutate(
      Group = factor(Group, levels = c("Low", "High")),
      Group_X = as.numeric(Group)
    )

  data.table::fwrite(
    long_df,
    file.path(DATA_DIR, paste0(stem, "_plot_data.csv"))
  )

  p <- ggplot(long_df, aes(x = Group, y = Median_Value, color = Group, group = Cancer)) +
    geom_line(color = "#8A96A3", linewidth = 0.45, alpha = 0.62) +
    geom_point(size = 1.65, alpha = 0.95) +
    facet_grid(Feature_Display ~ Cancer, scales = "free_y") +
    scale_color_manual(values = c(Low = pal_low, High = pal_high), name = "MO-DDRscore") +
    labs(
      title = paste0(panel, "  ", title),
      subtitle = "Paired group medians from the result table; full sample-level distributions require the original long source files",
      x = NULL,
      y = "Median feature value"
    ) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(color = pal_dark, face = "bold", size = 18, hjust = 0),
      plot.subtitle = element_text(color = "#718093", size = 9, hjust = 0),
      strip.background = element_rect(fill = "#EEF3F7", color = "#D6E0EA", linewidth = 0.45),
      strip.text = element_text(color = pal_dark, face = "bold", size = 7.4),
      axis.text.x = element_text(color = pal_dark, face = "bold", size = 6.7, angle = 45, hjust = 1),
      axis.text.y = element_text(color = pal_dark, size = 6.4),
      axis.title.y = element_text(color = pal_dark, face = "bold", size = 9),
      panel.grid.major = element_line(color = pal_grid, linewidth = 0.22),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "#D6E0EA", fill = NA, linewidth = 0.45),
      legend.position = "bottom",
      legend.title = element_text(color = pal_dark, face = "bold"),
      plot.margin = margin(10, 12, 10, 10)
    )

  ggsave(file.path(SUP_DIR, paste0(stem, ".png")), p, width = width, height = height, dpi = 600, bg = "white")
  ggsave(file.path(SUP_DIR, paste0(stem, ".pdf")), p, width = width, height = height, bg = "white", useDingbats = FALSE)
  ggsave(file.path(SUP_DIR, paste0(stem, ".tiff")), p, width = width, height = height, dpi = 600, bg = "white", compression = "lzw")

  invisible(long_df)
}

make_heatmap_panel(
  file = file.path(DATA_DIR, "Figure3C_immune_suppression_main_data.csv"),
  stem = "Figure3C_immune_suppression_heatmap_candidate",
  panel = "C",
  title = "Immune-suppressive programs",
  subtitle = "FDR-significant high-low differences in selected immune-suppressive programs"
)

make_heatmap_panel(
  file = file.path(DATA_DIR, "Figure3D_immune_exclusion_main_data.csv"),
  stem = "Figure3D_immune_exclusion_heatmap_candidate",
  panel = "D",
  title = "Immune-exclusion programs",
  subtitle = "FDR-significant high-low differences in selected immune-exclusion programs"
)

make_median_supp_panel(
  file = file.path(DATA_DIR, "Figure3C_immune_suppression_main_data.csv"),
  stem = "FigureS3C_immune_suppression_group_median_candidate",
  panel = "S3C",
  title = "Immune-suppressive programs"
)

make_median_supp_panel(
  file = file.path(DATA_DIR, "Figure3D_immune_exclusion_main_data.csv"),
  stem = "FigureS3D_immune_exclusion_group_median_candidate",
  panel = "S3D",
  title = "Immune-exclusion programs"
)

cat("Script finished.\n")
cat("Main figure output directory:\n", FIG_DIR, "\n")
cat("Supplementary figure output directory:\n", SUP_DIR, "\n")
