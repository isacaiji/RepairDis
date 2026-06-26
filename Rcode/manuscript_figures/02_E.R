############################################################
# Figure 1F: TMB and MATH differences by MO-DDRscore group
# - Keeps cancer types with at least one FDR-significant TMB/MATH result.
# - Asterisks are generated from FDR, not raw P.
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

ROOT_DIR <- "D:/R_workspace/summary/Repairium_pancancer_analysis_plan"
DATA_DIR <- file.path(ROOT_DIR, "03-res", "Figure2_pan_cancer", "E_TMB_MATH")
PLOT_ROOT <- file.path(ROOT_DIR, "03-res", "plots")
OUT_ROOT <- file.path(PLOT_ROOT, "张岩-图")
FIG_DIR <- file.path(OUT_ROOT, "Figure1")
DATA_OUT_DIR <- file.path(OUT_ROOT, "plot_data")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cancer_order <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
  "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
  "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
)

feature_map <- c(
  "log1p_TMB_per_Mb" = "TMB",
  "MATH" = "MATH"
)

pal_low <- "#27B7B8"
pal_high <- "#D45A63"
pal_dark <- "#18324A"
pal_grid <- "#E9EEF3"

sig_from_fdr <- function(fdr) {
  dplyr::case_when(
    is.na(fdr) ~ "",
    fdr < 0.001 ~ "***",
    fdr < 0.01 ~ "**",
    fdr < 0.05 ~ "*",
    TRUE ~ ""
  )
}

theme_repair <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = pal_dark),
      axis.text = element_text(color = pal_dark),
      axis.title = element_text(color = pal_dark, face = "bold"),
      axis.line = element_line(color = "#637485", linewidth = 0.45),
      axis.ticks = element_line(color = "#637485", linewidth = 0.35),
      strip.background = element_rect(fill = "#EFF3F7", color = "#D5DEE8", linewidth = 0.45),
      strip.text = element_text(color = pal_dark, face = "bold", size = base_size),
      legend.title = element_blank(),
      legend.text = element_text(color = pal_dark),
      plot.title = element_text(face = "bold", color = pal_dark, size = base_size + 5, hjust = 0),
      plot.subtitle = element_text(color = "#708296", size = base_size + 0.5, hjust = 0),
      plot.margin = margin(8, 10, 8, 8)
    )
}

save_all <- function(plot, stem, width, height) {
  pdf_file <- file.path(FIG_DIR, paste0(stem, ".pdf"))
  png_file <- file.path(FIG_DIR, paste0(stem, ".png"))
  tiff_file <- file.path(FIG_DIR, paste0(stem, ".tiff"))

  ggsave(pdf_file, plot, width = width, height = height, device = cairo_pdf)
  ggsave(png_file, plot, width = width, height = height, dpi = 600, bg = "white")
  ggsave(tiff_file, plot, width = width, height = height, dpi = 600, compression = "lzw", bg = "white")

  file.copy(pdf_file, file.path(PLOT_ROOT, paste0(stem, ".pdf")), overwrite = TRUE)
  file.copy(png_file, file.path(PLOT_ROOT, paste0(stem, ".png")), overwrite = TRUE)
  file.copy(tiff_file, file.path(PLOT_ROOT, paste0(stem, ".tiff")), overwrite = TRUE)
}

stat_files <- file.path(DATA_DIR, cancer_order, paste0(cancer_order, "_TMB_MATH_group_comparison.csv"))
stat_files <- stat_files[file.exists(stat_files)]
if (length(stat_files) == 0) stop("No TMB/MATH group comparison files found.")

stat_df <- dplyr::bind_rows(lapply(stat_files, function(f) {
  data.table::fread(f, data.table = FALSE, check.names = FALSE)
})) %>%
  dplyr::filter(Feature %in% names(feature_map)) %>%
  dplyr::mutate(
    Cancer = factor(Cancer, levels = cancer_order),
    Feature_Display = factor(unname(feature_map[Feature]), levels = c("TMB", "MATH")),
    P = as.numeric(P),
    FDR = as.numeric(FDR),
    Effect = as.numeric(Effect),
    Significance = sig_from_fdr(FDR),
    Direction = dplyr::case_when(
      Effect > 0 ~ "Higher in High",
      Effect < 0 ~ "Higher in Low",
      TRUE ~ "No median difference"
    )
  ) %>%
  dplyr::filter(!is.na(Cancer), !is.na(Feature_Display), is.finite(FDR))

keep_cancers <- stat_df %>%
  dplyr::group_by(Cancer) %>%
  dplyr::summarise(Any_FDR_lt_0.05 = any(FDR < 0.05, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(Any_FDR_lt_0.05) %>%
  dplyr::arrange(Cancer) %>%
  dplyr::pull(Cancer) %>%
  as.character()

if (length(keep_cancers) == 0) stop("No cancer type has FDR-significant TMB or MATH result.")

long_files <- file.path(DATA_DIR, keep_cancers, paste0(keep_cancers, "_TMB_MATH_long.csv"))
long_files <- long_files[file.exists(long_files)]

long_df <- dplyr::bind_rows(lapply(long_files, function(f) {
  data.table::fread(f, data.table = FALSE, check.names = FALSE)
})) %>%
  dplyr::filter(Feature %in% names(feature_map), Cancer %in% keep_cancers) %>%
  dplyr::mutate(
    Cancer = factor(Cancer, levels = keep_cancers),
    Feature_Display = factor(unname(feature_map[Feature]), levels = c("TMB", "MATH")),
    MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High")),
    Value = as.numeric(Value)
  ) %>%
  dplyr::filter(!is.na(Feature_Display), !is.na(MO_DDRscore_group), is.finite(Value))

plot_df <- long_df %>%
  dplyr::group_by(Cancer, Feature, Feature_Display, MO_DDRscore_group) %>%
  dplyr::summarise(
    N = dplyr::n(),
    Median = median(Value, na.rm = TRUE),
    Q1 = quantile(Value, 0.25, na.rm = TRUE),
    Q3 = quantile(Value, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    stat_df %>% dplyr::select(Cancer, Feature, P, FDR, Effect, Significance, Direction),
    by = c("Cancer", "Feature")
  )

feature_y <- plot_df %>%
  dplyr::group_by(Feature_Display) %>%
  dplyr::summarise(
    y_min = min(Q1, Median, na.rm = TRUE),
    y_max = max(Q3, Median, na.rm = TRUE),
    y_star = y_min + 0.92 * (y_max - y_min),
    .groups = "drop"
  )

label_df <- plot_df %>%
  dplyr::group_by(Cancer, Feature_Display) %>%
  dplyr::summarise(
    Significance = dplyr::first(Significance),
    FDR = dplyr::first(FDR),
    .groups = "drop"
  ) %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::left_join(feature_y, by = "Feature_Display")

data.table::fwrite(stat_df, file.path(DATA_OUT_DIR, "Figure1F_TMB_MATH_all_stats_FDR.csv"))
data.table::fwrite(
  stat_df %>% dplyr::filter(Cancer %in% keep_cancers),
  file.path(DATA_OUT_DIR, "Figure1F_TMB_MATH_retained_cancer_stats_FDR.csv")
)
data.table::fwrite(plot_df, file.path(DATA_OUT_DIR, "Figure1F_TMB_MATH_plot_data.csv"))

p <- ggplot(plot_df, aes(x = Cancer, y = Median, color = MO_DDRscore_group)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.32, color = "#B8C3CC") +
  geom_linerange(
    aes(ymin = Q1, ymax = Q3),
    position = position_dodge(width = 0.48),
    linewidth = 0.55,
    alpha = 0.88
  ) +
  geom_point(
    aes(fill = MO_DDRscore_group),
    position = position_dodge(width = 0.48),
    shape = 21,
    size = 2.15,
    stroke = 0.35,
    alpha = 0.96
  ) +
  geom_text(
    data = label_df,
    aes(x = Cancer, y = y_star, label = Significance),
    inherit.aes = FALSE,
    color = pal_dark,
    fontface = "bold",
    size = 3.0,
    vjust = 0
  ) +
  facet_grid(Feature_Display ~ ., scales = "free_y", switch = "y") +
  scale_color_manual(values = c(Low = pal_low, High = pal_high), name = "MO-DDRscore") +
  scale_fill_manual(values = c(Low = pal_low, High = pal_high), name = "MO-DDRscore") +
  labs(
    title = "F  Tumor mutation burden and MATH",
    subtitle = "Cancer types with at least one FDR-significant TMB or MATH difference are shown",
    x = NULL,
    y = "Median value with interquartile range"
  ) +
  theme_repair(10) +
  theme(
    panel.grid.major.y = element_line(color = pal_grid, linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8.3, face = "bold"),
    axis.text.y = element_text(size = 8.5),
    legend.position = "top",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, size = 10.5, face = "bold"),
    strip.background = element_rect(fill = "#EFF3F7", color = "#D5DEE8", linewidth = 0.45)
  )

save_all(p, "Figure1F_TMB_MATH", width = max(7.8, length(keep_cancers) * 0.30 + 2.1), height = 4.8)

summary_df <- data.frame(
  Item = c(
    "N_cancers_total",
    "N_cancers_retained",
    "Retain_rule",
    "Significance_rule"
  ),
  Value = c(
    length(cancer_order),
    length(keep_cancers),
    "Retain cancer if TMB or MATH has FDR < 0.05",
    "Asterisks from FDR only: * < 0.05, ** < 0.01, *** < 0.001"
  )
)
data.table::fwrite(summary_df, file.path(DATA_OUT_DIR, "Figure1F_TMB_MATH_run_summary.csv"))

cat("Done.\n")
cat("Retained cancers:", paste(keep_cancers, collapse = ", "), "\n")
cat("Output:\n", file.path(FIG_DIR, "Figure1F_TMB_MATH.pdf"), "\n")
