rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

BASE_DIR <- "D:/R_workspace/summary/Repairium_pancancer_analysis_plan"
SRC_FILE <- file.path(
  BASE_DIR, "03-res", "Figure2_pan_cancer", "H_Immune_suppression",
  "PanCancer_Immune_suppression_group_comparison_all.csv"
)
PLOT_DIR <- file.path(BASE_DIR, "03-res", "plots", "张岩-图")
FIG_DIR <- file.path(PLOT_DIR, "Figure3")
SUP_DIR <- file.path(PLOT_DIR, "sup", "Figure3")
DATA_DIR <- file.path(PLOT_DIR, "plot_data")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

pal_low <- "#24B7B6"
pal_high <- "#D85B68"
pal_dark <- "#10243C"
pal_grid <- "#E6EEF5"
pal_mid <- "#F7F4EF"

cancer_order_all <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
  "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
  "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
)

feature_order <- c("MDSC", "Tregs", "TAM M2", "T-cell dysfunction")

fdr_star <- function(x) {
  dplyr::case_when(
    is.na(x) ~ "",
    x < 1e-4 ~ "****",
    x < 1e-3 ~ "***",
    x < 1e-2 ~ "**",
    x < 0.05 ~ "*",
    TRUE ~ ""
  )
}

map_feature <- function(feature) {
  f <- toupper(as.character(feature))
  dplyr::case_when(
    grepl("DYSFUNCTION", f) ~ "T-cell dysfunction",
    grepl("MDSC", f) ~ "MDSC",
    grepl("TREG|REGULATORY", f) ~ "Tregs",
    grepl("TAM.*M2|M2", f) ~ "TAM M2",
    TRUE ~ NA_character_
  )
}

safe_spread <- function(x) {
  s <- mad(x, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(s) || s == 0) s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) s <- 1
  s
}

order_cancers_by_color <- function(df) {
  df %>%
    group_by(Cancer) %>%
    summarise(
      Mean_Diff = mean(Diff_Z, na.rm = TRUE),
      Sig_n = sum(FDR < 0.05, na.rm = TRUE),
      Abs_Diff = sum(abs(Diff_Z), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Mean_Diff), desc(Sig_n), desc(Abs_Diff)) %>%
    pull(Cancer)
}

make_heatmap <- function(plot_df, cancer_order, title_text, output_stem,
                         out_dir, width, height, star_size = 5.55) {
  plot_df <- plot_df %>%
    mutate(
      Cancer = factor(Cancer, levels = cancer_order),
      Feature_Display = factor(Feature_Display, levels = rev(feature_order)),
      Diff_Z_plot = pmax(pmin(Diff_Z, 2.5), -2.5)
    )

  p <- ggplot(plot_df, aes(x = Cancer, y = Feature_Display)) +
    geom_tile(aes(fill = Diff_Z_plot), color = "white", linewidth = 1.05) +
    geom_text(aes(label = FDR_star), color = pal_dark, fontface = "bold", size = star_size) +
    scale_fill_gradient2(
      low = pal_low, mid = pal_mid, high = pal_high, midpoint = 0,
      limits = c(-2.5, 2.5), oob = scales::squish,
      name = "High - Low\ndifference"
    ) +
    labs(title = title_text, x = "Cancer type", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 27, color = pal_dark, hjust = 0),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", color = pal_dark, size = 11),
      axis.text.y = element_text(face = "bold", color = pal_dark, size = 12),
      axis.title.x = element_text(face = "bold", color = pal_dark, size = 14),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "#CBD8E3", fill = NA, linewidth = 0.7),
      legend.title = element_text(face = "bold", color = pal_dark, size = 11),
      legend.text = element_text(color = pal_dark, size = 10),
      plot.margin = margin(10, 16, 8, 10)
    )

  ggsave(file.path(out_dir, paste0(output_stem, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(out_dir, paste0(output_stem, ".png")), p, width = width, height = height, dpi = 450)
  ggsave(file.path(out_dir, paste0(output_stem, ".tiff")), p, width = width, height = height, dpi = 450, compression = "lzw")
  p
}

raw <- fread(SRC_FILE, data.table = FALSE, check.names = FALSE)

dat <- raw %>%
  mutate(
    Cancer = as.character(Cancer),
    Feature = as.character(Feature),
    Method = as.character(Method),
    Feature_Display = map_feature(Feature),
    Median_Low = as.numeric(Median_Low),
    Median_High = as.numeric(Median_High),
    Effect = as.numeric(Effect),
    P = as.numeric(P),
    FDR = as.numeric(FDR)
  ) %>%
  filter(Cancer %in% cancer_order_all, Feature_Display %in% feature_order)

scale_tbl <- dat %>%
  group_by(Feature_Display) %>%
  summarise(
    Center = median(c(Median_Low, Median_High), na.rm = TRUE),
    Spread = safe_spread(c(Median_Low, Median_High)),
    .groups = "drop"
  )

plot_df <- dat %>%
  left_join(scale_tbl, by = "Feature_Display") %>%
  mutate(
    Low_Z = (Median_Low - Center) / Spread,
    High_Z = (Median_High - Center) / Spread,
    Diff_Z = High_Z - Low_Z,
    FDR_star = fdr_star(FDR)
  ) %>%
  arrange(Cancer, Feature_Display, FDR, desc(abs(Diff_Z))) %>%
  group_by(Cancer, Feature_Display) %>%
  slice(1) %>%
  ungroup()

selected_cancers <- c(
  "LUAD", "UCEC", "BRCA", "HNSC",
  "LUSC", "BLCA", "ESCA", "LGG"
)

main_df <- plot_df %>% filter(Cancer %in% selected_cancers)
main_cancer_order <- order_cancers_by_color(main_df)
full_cancer_order <- order_cancers_by_color(plot_df)

fwrite(plot_df, file.path(DATA_DIR, "Figure3C_immune_suppression_full_plot_data.csv"))
fwrite(main_df, file.path(DATA_DIR, "Figure3C_immune_suppression_main_plot_data.csv"))
fwrite(
  data.frame(Cancer = main_cancer_order),
  file.path(DATA_DIR, "Figure3C_immune_suppression_selected_cancers.csv")
)

make_heatmap(
  main_df, main_cancer_order,
  "C  Immune-suppressive programs",
  "Figure3C_immune_suppression",
  FIG_DIR, width = 12.0, height = 5.2
)

make_heatmap(
  plot_df, full_cancer_order,
  "Figure S3C  Immune-suppressive programs",
  "FigureS3C_immune_suppression_full_heatmap",
  SUP_DIR, width = 18.5, height = 5.3, star_size = 3.25
)

legend_txt <- c(
  "Figure 3C. Immune-suppressive programs across MO-DDRscore-defined tumor states.",
  "Representative cancer types were selected according to the number and magnitude of FDR-significant differences between MO-DDRscore-high and -low tumors.",
  "Tile colors indicate feature-wise standardized median-score differences between high- and low-score tumors (High - Low).",
  "Asterisks denote FDR-adjusted significance: * FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001, **** FDR < 0.0001.",
  "",
  "Figure S3C. Full pan-cancer immune-suppressive program landscape across all available cancer types."
)
writeLines(legend_txt, file.path(DATA_DIR, "Figure3C_legend.txt"), useBytes = TRUE)

cat("Done: Figure3C immune-suppressive programs\n")
cat("Main output:", file.path(FIG_DIR, "Figure3C_immune_suppression.png"), "\n")
cat("Supplementary output:", file.path(SUP_DIR, "FigureS3C_immune_suppression_full_heatmap.png"), "\n")
