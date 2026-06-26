rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

BASE_DIR <- "D:/R_workspace/summary/Repairium_pancancer_analysis_plan"
SRC_FILE <- file.path(
  BASE_DIR, "03-res", "Figure2_pan_cancer", "H_Immune_suppression",
  "PanCancer_Immune_suppression_group_comparison_all.csv"
)
PLOT_DIR <- file.path(BASE_DIR, "03-res", "plots", "\u5f20\u5ca9-\u56fe")
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
pal_grey <- "#8E969C"

cancer_order_all <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
  "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
  "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
)

fdr_star <- function(x) {
  dplyr::case_when(
    is.na(x) ~ "",
    x < 1e-4 ~ "****",
    x < 1e-3 ~ "***",
    x < 1e-2 ~ "**",
    x < 5e-2 ~ "*",
    TRUE ~ ""
  )
}

safe_z <- function(x, center, scale) {
  if (!is.finite(scale) || scale == 0) scale <- 1
  (x - center) / scale
}

save_plot <- function(p, stem, width, height) {
  ggsave(file.path(FIG_DIR, paste0(stem, ".png")), p, width = width, height = height, dpi = 420, bg = "white")
  ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), p, width = width, height = height, bg = "white")
  ggsave(file.path(FIG_DIR, paste0(stem, ".tiff")), p, width = width, height = height, dpi = 420, bg = "white", compression = "lzw")
}

save_sup_plot <- function(p, stem, width, height) {
  ggsave(file.path(SUP_DIR, paste0(stem, ".png")), p, width = width, height = height, dpi = 420, bg = "white")
  ggsave(file.path(SUP_DIR, paste0(stem, ".pdf")), p, width = width, height = height, bg = "white")
  ggsave(file.path(SUP_DIR, paste0(stem, ".tiff")), p, width = width, height = height, dpi = 420, bg = "white", compression = "lzw")
}

raw <- fread(SRC_FILE, data.table = FALSE, check.names = FALSE) %>%
  mutate(
    Feature_ID = paste(Feature, Method, sep = "__"),
    Cancer = factor(Cancer, levels = cancer_order_all),
    Median_Low = as.numeric(Median_Low),
    Median_High = as.numeric(Median_High),
    Effect = as.numeric(Effect),
    P = as.numeric(P),
    FDR = as.numeric(FDR),
    FDR_star = fdr_star(FDR)
  ) %>%
  filter(!is.na(Cancer))

feature_map <- data.frame(
  Feature_ID = c(
    "Dysfunction__TIDEpy",
    "MDSC__TIDEpy",
    "Tregs_quantiseq__quanTIseq",
    "TAM M2__TIDEpy"
  ),
  Feature_Display = c("T-cell dysfunction", "MDSC", "Tregs", "TAM M2"),
  Feature_Order = 1:4,
  stringsAsFactors = FALSE
)

standardized <- raw %>%
  inner_join(feature_map, by = "Feature_ID") %>%
  group_by(Feature_ID) %>%
  mutate(
    .center = mean(c(Median_Low, Median_High), na.rm = TRUE),
    .scale = sd(c(Median_Low, Median_High), na.rm = TRUE),
    Low_Z = safe_z(Median_Low, .center[1], .scale[1]),
    High_Z = safe_z(Median_High, .center[1], .scale[1]),
    Diff_Z = High_Z - Low_Z
  ) %>%
  ungroup()

selection_tbl <- standardized %>%
  group_by(Cancer) %>%
  summarise(
    FDR_sig_count = sum(FDR < 0.05, na.rm = TRUE),
    mean_abs_difference = mean(abs(Diff_Z), na.rm = TRUE),
    mean_signed_difference = mean(Diff_Z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(FDR_sig_count), desc(mean_abs_difference))

main_cancers <- selection_tbl %>%
  filter(as.character(Cancer) != "LUAD") %>%
  slice_head(n = 7) %>%
  pull(Cancer) %>%
  as.character()
main_cancers <- unique(c(main_cancers, "LUAD"))

cancer_order <- standardized %>%
  filter(as.character(Cancer) %in% main_cancers) %>%
  group_by(Cancer) %>%
  summarise(
    signed_score = mean(Diff_Z, na.rm = TRUE),
    sig_count = sum(FDR < 0.05, na.rm = TRUE),
    abs_score = mean(abs(Diff_Z), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(sig_count), desc(signed_score), desc(abs_score)) %>%
  pull(Cancer) %>%
  as.character()

main_df <- standardized %>%
  filter(as.character(Cancer) %in% cancer_order) %>%
  mutate(
    Cancer = factor(as.character(Cancer), levels = rev(cancer_order)),
    Feature_Display = factor(Feature_Display, levels = feature_map$Feature_Display)
  )

fwrite(selection_tbl, file.path(DATA_DIR, "Figure3C_immune_suppression_cancer_selection.csv"))
fwrite(main_df, file.path(DATA_DIR, "Figure3C_immune_suppression_main_data.csv"))

point_df <- bind_rows(
  main_df %>% transmute(Cancer, Feature_Display, Group = "Low", Score = Low_Z),
  main_df %>% transmute(Cancer, Feature_Display, Group = "High", Score = High_Z)
) %>%
  mutate(Group = factor(Group, levels = c("Low", "High")))

x_abs <- max(abs(c(main_df$Low_Z, main_df$High_Z)), na.rm = TRUE)
x_lim <- max(2.5, ceiling((x_abs + 0.15) * 2) / 2)
x_lim <- min(max(x_lim, 2.5), 4.5)
star_x <- x_lim - 0.12

star_df <- main_df %>%
  filter(FDR_star != "") %>%
  mutate(Star_X = star_x)

p_main <- ggplot(main_df, aes(y = Cancer)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.55, color = "#A9B7C6") +
  geom_segment(aes(x = Low_Z, xend = High_Z, yend = Cancer), color = pal_grey, linewidth = 0.72, alpha = 0.9) +
  geom_point(
    data = point_df,
    aes(x = Score, y = Cancer, color = Group),
    size = 2.65,
    alpha = 0.98
  ) +
  geom_text(
    data = star_df,
    aes(x = Star_X, y = Cancer, label = FDR_star),
    inherit.aes = FALSE,
    color = pal_dark,
    fontface = "bold",
    size = 3.35
  ) +
  facet_wrap(~Feature_Display, nrow = 1) +
  scale_color_manual(values = c(Low = pal_low, High = pal_high), name = "MO-DDRscore") +
  scale_x_continuous(limits = c(-x_lim, x_lim), breaks = pretty(c(-x_lim, x_lim), n = 5)) +
  labs(
    title = "C  Immune-suppressive programs",
    x = "Standardized median program score",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 23, color = pal_dark),
    strip.background = element_rect(fill = "#EDF3F7", color = "#D3DEE8", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 11.5, color = pal_dark),
    panel.border = element_rect(color = "#CAD7E3", linewidth = 0.8),
    panel.grid.major = element_line(color = pal_grid, linewidth = 0.7),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "#2C4056", face = "bold"),
    axis.title.x = element_text(color = pal_dark, face = "bold", size = 13),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", color = pal_dark),
    legend.text = element_text(color = "#3A4C62"),
    plot.margin = margin(10, 18, 8, 8)
  )

save_plot(p_main, "Figure3C_immune_suppression", width = 10.2, height = 4.9)

# Supplementary full landscape: all 33 cancers and all nonredundant suppression axes.
heat_df <- standardized %>%
  mutate(
    Feature_Display = factor(Feature_Display, levels = rev(feature_map$Feature_Display)),
    Cancer_chr = as.character(Cancer)
  )

mat <- heat_df %>%
  select(Feature_Display, Cancer_chr, Diff_Z) %>%
  tidyr::pivot_wider(names_from = Cancer_chr, values_from = Diff_Z) %>%
  as.data.frame()
rownames(mat) <- as.character(mat$Feature_Display)
mat$Feature_Display <- NULL
mat[is.na(mat)] <- 0
if (ncol(mat) > 2) {
  hc <- hclust(dist(t(as.matrix(mat))), method = "ward.D2")
  cancer_order_heat <- colnames(mat)[hc$order]
} else {
  cancer_order_heat <- colnames(mat)
}

heat_df <- heat_df %>%
  mutate(Cancer_chr = factor(Cancer_chr, levels = cancer_order_heat))

fwrite(heat_df, file.path(DATA_DIR, "FigureS3C_immune_suppression_full_heatmap_data.csv"))

p_sup <- ggplot(heat_df, aes(x = Cancer_chr, y = Feature_Display, fill = Diff_Z)) +
  geom_tile(color = "white", linewidth = 0.85) +
  geom_text(aes(label = FDR_star), color = pal_dark, fontface = "bold", size = 4.15) +
  scale_fill_gradient2(
    low = pal_low,
    mid = "#F7F3ED",
    high = pal_high,
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    name = "High - Low\ndifference"
  ) +
  labs(
    title = "Figure S3C  Full immune-suppressive program landscape",
    x = "Cancer type",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 18, color = pal_dark),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "#CAD7E3", linewidth = 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "#2C4056", face = "bold"),
    axis.text.y = element_text(color = "#2C4056", face = "bold"),
    axis.title.x = element_text(color = pal_dark, face = "bold"),
    legend.title = element_text(face = "bold", color = pal_dark),
    legend.text = element_text(color = "#2C4056"),
    plot.margin = margin(8, 10, 8, 8)
  )

save_sup_plot(p_sup, "FigureS3C_immune_suppression_full_heatmap", width = 11.5, height = 3.8)

cat("Done Figure3C.\n")
cat("Main:", file.path(FIG_DIR, "Figure3C_immune_suppression.png"), "\n")
cat("Supplement:", file.path(SUP_DIR, "FigureS3C_immune_suppression_full_heatmap.png"), "\n")
