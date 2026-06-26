rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

BASE_DIR <- "D:/R_workspace/summary/Repairium_pancancer_analysis_plan"
PLOT_DIR <- file.path(BASE_DIR, "03-res", "plots", "\u5f20\u5ca9-\u56fe")
FIG_DIR <- file.path(PLOT_DIR, "Figure3")
DATA_DIR <- file.path(PLOT_DIR, "plot_data")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

pal_low <- "#24B7B6"
pal_high <- "#D85B68"
pal_dark <- "#10243C"
pal_grid <- "#E6EEF5"
pal_axis <- "#6D7F91"
pal_grey <- "#8D969F"

save_preview <- function(p, stem, width, height) {
  ggsave(file.path(FIG_DIR, paste0(stem, ".png")), p, width = width, height = height, dpi = 420, bg = "white")
  ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), p, width = width, height = height, bg = "white")
  ggsave(file.path(FIG_DIR, paste0(stem, ".tiff")), p, width = width, height = height, dpi = 420, bg = "white", compression = "lzw")
}

clean_main <- function(file) {
  fread(file, data.table = FALSE, check.names = FALSE) %>%
    mutate(
      Cancer = as.character(Cancer),
      Feature_Display = as.character(Feature_Display),
      Feature_Order = as.numeric(Feature_Order),
      Low_Z = as.numeric(Low_Z),
      High_Z = as.numeric(High_Z),
      Diff_Z = as.numeric(Diff_Z),
      FDR = as.numeric(FDR),
      FDR_star = ifelse(is.na(FDR_star), "", as.character(FDR_star)),
      Direction = case_when(
        Diff_Z > 0 ~ "Higher in High",
        Diff_Z < 0 ~ "Higher in Low",
        TRUE ~ "No difference"
      )
    ) %>%
    filter(
      is.finite(Low_Z),
      is.finite(High_Z),
      is.finite(Diff_Z),
      is.finite(FDR)
    )
}

make_order <- function(df) {
  df %>%
    group_by(Cancer) %>%
    summarise(
      sig_count = sum(FDR < 0.05, na.rm = TRUE),
      signed_effect = mean(Diff_Z, na.rm = TRUE),
      abs_effect = mean(abs(Diff_Z), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(sig_count), desc(signed_effect), desc(abs_effect)) %>%
    pull(Cancer)
}

make_combo_plot <- function(file, stem, panel, title, subtitle, width, height) {
  df <- clean_main(file)
  cancer_order <- make_order(df)
  feature_order <- df %>%
    distinct(Feature_Display, Feature_Order) %>%
    arrange(Feature_Order) %>%
    pull(Feature_Display)

  df <- df %>%
    mutate(
      Cancer = factor(Cancer, levels = rev(cancer_order)),
      Feature_Display = factor(Feature_Display, levels = feature_order)
    )

  score_long <- bind_rows(
    df %>% transmute(Cancer, Feature_Display, Group = "Low", Score = Low_Z),
    df %>% transmute(Cancer, Feature_Display, Group = "High", Score = High_Z)
  ) %>%
    mutate(Group = factor(Group, levels = c("Low", "High")))

  x_score <- max(abs(score_long$Score), na.rm = TRUE)
  x_score <- max(2.5, ceiling(x_score * 2) / 2)
  x_eff <- max(abs(df$Diff_Z), na.rm = TRUE)
  x_eff <- max(1.5, ceiling(x_eff * 2) / 2)

  star_df <- df %>%
    filter(FDR_star != "") %>%
    mutate(
      Star_X = ifelse(Diff_Z >= 0, x_eff * 0.92, -x_eff * 0.92),
      Hjust = ifelse(Diff_Z >= 0, 0.5, 0.5)
    )

  base_theme <- theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = pal_grid, linewidth = 0.42),
      panel.grid.major.x = element_line(color = pal_grid, linewidth = 0.42),
      panel.border = element_rect(color = "#CCD8E3", linewidth = 0.72),
      axis.text = element_text(color = pal_dark, face = "bold"),
      axis.title = element_text(color = pal_dark, face = "bold"),
      strip.background = element_rect(fill = "#EEF3F8", color = "#CCD8E3", linewidth = 0.72),
      strip.text.y.left = element_text(color = pal_dark, face = "bold", angle = 0, size = 10),
      strip.placement = "outside",
      legend.position = "bottom",
      legend.title = element_text(color = pal_dark, face = "bold"),
      legend.text = element_text(color = "#35475A"),
      plot.margin = margin(4, 4, 4, 4)
    )

  p_score <- ggplot(score_long, aes(x = Score, y = Cancer, fill = Group)) +
    geom_vline(xintercept = 0, color = "#AAB8C6", linewidth = 0.55, linetype = "22") +
    geom_col(
      position = position_dodge(width = 0.68),
      width = 0.28,
      color = "white",
      linewidth = 0.18,
      alpha = 0.94
    ) +
    facet_grid(Feature_Display ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = c(Low = pal_low, High = pal_high), name = "MO-DDRscore") +
    scale_x_continuous(limits = c(-x_score, x_score), breaks = pretty(c(-x_score, x_score), n = 5)) +
    labs(x = "Standardized median score", y = NULL) +
    base_theme +
    theme(
      strip.text.y.left = element_text(margin = margin(r = 7)),
      axis.text.y = element_text(size = 8.8),
      legend.position = "bottom"
    )

  p_effect <- ggplot(df, aes(y = Cancer)) +
    geom_vline(xintercept = 0, color = "#AAB8C6", linewidth = 0.55, linetype = "22") +
    geom_segment(
      aes(x = 0, xend = Diff_Z, yend = Cancer),
      color = pal_grey,
      linewidth = 0.82,
      lineend = "round"
    ) +
    geom_point(aes(x = Diff_Z, color = Direction), size = 3.0, alpha = 0.96) +
    geom_text(
      data = star_df,
      aes(x = Star_X, y = Cancer, label = FDR_star, hjust = Hjust),
      inherit.aes = FALSE,
      color = pal_dark,
      fontface = "bold",
      size = 3.35
    ) +
    facet_grid(Feature_Display ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_color_manual(
      values = c("Higher in Low" = pal_low, "Higher in High" = pal_high, "No difference" = "#C9D1D9"),
      name = "Effect direction"
    ) +
    scale_x_continuous(limits = c(-x_eff, x_eff), breaks = pretty(c(-x_eff, x_eff), n = 5)) +
    labs(x = "High - Low difference", y = NULL) +
    base_theme +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text.y.left = element_blank(),
      strip.background = element_blank(),
      legend.position = "bottom"
    )

  p <- p_score + p_effect + plot_layout(widths = c(1.32, 0.92), guides = "collect")

  p <- p +
    plot_annotation(
      title = paste0(panel, "  ", title),
      subtitle = subtitle,
      theme = theme(
        plot.title = element_text(face = "bold", size = 21, color = pal_dark, hjust = 0.02),
        plot.subtitle = element_text(size = 10.5, color = "#718093", hjust = 0.02),
        plot.margin = margin(6, 8, 6, 8)
      )
    )

  fwrite(df, file.path(DATA_DIR, paste0(stem, "_plot_data.csv")))
  save_preview(p, stem, width, height)
  invisible(p)
}

make_combo_plot(
  file = file.path(DATA_DIR, "Figure3C_immune_suppression_main_data.csv"),
  stem = "Figure3C_immune_suppression_forest_preview",
  panel = "C",
  title = "Immune-suppressive programs",
  subtitle = "Left: standardized low/high group medians; right: paired High-Low difference with FDR-based significance",
  width = 13.5,
  height = 7.1
)

make_combo_plot(
  file = file.path(DATA_DIR, "Figure3D_immune_exclusion_main_data.csv"),
  stem = "Figure3D_immune_exclusion_forest_preview",
  panel = "D",
  title = "Immune-exclusion programs",
  subtitle = "Left: standardized low/high group medians; right: paired High-Low difference with FDR-based significance",
  width = 13.5,
  height = 9.2
)

cat("Done.\n")
cat("Preview outputs:\n")
cat(file.path(FIG_DIR, "Figure3C_immune_suppression_forest_preview.png"), "\n")
cat(file.path(FIG_DIR, "Figure3D_immune_exclusion_forest_preview.png"), "\n")
