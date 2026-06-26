############################################################
# RepairDis Figure 3E
# Checkpoint biomarker landscape
# - Main panel keeps all 33 cancer types
# - Cancer columns are reordered by the overall High-Low checkpoint pattern
# - Asterisks are generated from FDR, not raw P
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
SOURCE_FILE <- file.path(
  BASE_DIR,
  "03-res", "Figure2_pan_cancer", "J_Immune_checkpoint_biomarkers",
  "PanCancer_checkpoint_biomarkers_group_comparison_all.csv"
)
PLOT_DIR <- file.path(BASE_DIR, "03-res", "plots", "张岩-图")
FIG_DIR <- file.path(PLOT_DIR, "Figure3")
SUP_DIR <- file.path(PLOT_DIR, "sup", "Figure3")
DATA_DIR <- file.path(PLOT_DIR, "plot_data")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

pal_low <- "#20AEB3"
pal_high <- "#D65A67"
pal_dark <- "#10243C"
pal_grid <- "#E7EDF3"
pal_mid <- "#F8F3EE"

cancer_order <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
  "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
  "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
)

sig_from_fdr <- function(fdr, ns_blank = TRUE) {
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

cap <- function(x, lim = 2.5) {
  pmax(pmin(x, lim), -lim)
}

safe_scale <- function(x) {
  x <- as.numeric(x)
  sx <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / sx)
}

clean_checkpoint_label <- function(x) {
  x <- as.character(x)
  dplyr::recode(
    x,
    "CD40lg" = "CD40LG",
    "TNFrsf4" = "TNFRSF4",
    "IFNg" = "IFNG",
    "PDCD1lg2" = "PDCD1LG2",
    "Icos" = "ICOS",
    "Icoslg" = "ICOSLG",
    "Btla" = "BTLA",
    "Cd28" = "CD28",
    "Ctl" = "CTL",
    "Msi Score" = "MSI score",
    .default = x
  )
}

feature_meta <- data.frame(
  Feature_Display = c(
    "CTL", "CD8", "IFNG",
    "CD40LG", "CD40", "CD80", "CD86", "CD28", "ICOS", "ICOSLG", "TNFRSF4",
    "CD274", "PDCD1LG2", "PDCD1", "CTLA4", "TIGIT", "LAG3", "HAVCR2",
    "IDO1", "BTLA", "CD276", "CD47",
    "MSI score"
  ),
  Feature_Group = c(
    rep("Immune effector", 3),
    rep("Co-stimulatory", 8),
    rep("Checkpoint / inhibitory", 11),
    "Response biomarker"
  ),
  Group_Order = c(
    rep(1, 3),
    rep(2, 8),
    rep(3, 11),
    4
  ),
  stringsAsFactors = FALSE
)

theme_repair <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(color = pal_dark),
      plot.title = element_text(face = "bold", size = base_size + 7, hjust = 0),
      plot.subtitle = element_text(size = base_size + 1.0, color = "#718093", hjust = 0),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = pal_grid, linewidth = 0.24),
      panel.border = element_rect(color = "#CAD5E0", fill = NA, linewidth = 0.45),
      axis.text = element_text(color = "#263849"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(color = "#263849"),
      plot.margin = margin(8, 14, 8, 8)
    )
}

save_plot <- function(p, stem, width, height) {
  pdf_file <- file.path(FIG_DIR, paste0(stem, ".pdf"))
  png_file <- file.path(FIG_DIR, paste0(stem, ".png"))
  tiff_file <- file.path(FIG_DIR, paste0(stem, ".tiff"))

  ggsave(pdf_file, p, width = width, height = height, useDingbats = FALSE)
  ggsave(png_file, p, width = width, height = height, dpi = 450, bg = "white")
  ggsave(tiff_file, p, width = width, height = height, dpi = 600, bg = "white", compression = "lzw")
}

if (!file.exists(SOURCE_FILE)) {
  stop("Cannot find input file: ", SOURCE_FILE)
}

checkpoint <- fread(SOURCE_FILE, data.table = FALSE, colClasses = "character") %>%
  mutate(
    Cancer = as.character(Cancer),
    Feature = as.character(Feature),
    Display = as.character(Display),
    Feature_Display = clean_checkpoint_label(Display),
    Median_Low = as.numeric(Median_Low),
    Median_High = as.numeric(Median_High),
    Effect = as.numeric(Effect),
    P = as.numeric(P),
    FDR = as.numeric(FDR)
  ) %>%
  filter(
    Cancer %in% cancer_order,
    is.finite(Effect),
    Feature_Display %in% feature_meta$Feature_Display
  ) %>%
  left_join(feature_meta, by = "Feature_Display") %>%
  group_by(Feature_Display) %>%
  mutate(
    EffectScaledWithinFeature = safe_scale(Effect),
    PlotValue = cap(EffectScaledWithinFeature, 2.5)
  ) %>%
  ungroup() %>%
  mutate(
    Significance = sig_from_fdr(FDR, ns_blank = TRUE),
    Direction = dplyr::case_when(
      is.finite(FDR) & FDR < 0.05 & Effect > 0 ~ "Higher in High",
      is.finite(FDR) & FDR < 0.05 & Effect < 0 ~ "Higher in Low",
      TRUE ~ "NS"
    ),
    NegLogFDR = -log10(pmax(FDR, .Machine$double.xmin))
  )

if (nrow(checkpoint) == 0) {
  stop("No checkpoint biomarker statistics were loaded.")
}

# Column order: cancers with stronger checkpoint elevation in the
# MO-DDRscore-high group are placed to the left; low-group enriched cancers
# are placed to the right. This keeps all 33 cancer types while making the
# pan-cancer pattern readable.
cancer_rank <- checkpoint %>%
  group_by(Cancer) %>%
  summarise(
    Net_FDR_weighted_direction = sum(
      ifelse(is.finite(FDR) & FDR < 0.05, sign(Effect) * pmin(NegLogFDR, 6), 0),
      na.rm = TRUE
    ),
    Mean_scaled_effect = mean(PlotValue, na.rm = TRUE),
    N_high = sum(is.finite(FDR) & FDR < 0.05 & Effect > 0, na.rm = TRUE),
    N_low = sum(is.finite(FDR) & FDR < 0.05 & Effect < 0, na.rm = TRUE),
    N_sig = N_high + N_low,
    .groups = "drop"
  ) %>%
  arrange(
    desc(Net_FDR_weighted_direction),
    desc(Mean_scaled_effect),
    desc(N_sig),
    Cancer
  )

col_order <- cancer_rank$Cancer

# Row order: keep functional blocks, then order each block by its average
# High-Low pattern. This is more interpretable than unsupervised clustering for
# checkpoint genes.
feature_rank <- checkpoint %>%
  group_by(Feature_Display, Feature_Group, Group_Order) %>%
  summarise(
    Mean_scaled_effect = mean(PlotValue, na.rm = TRUE),
    N_high = sum(is.finite(FDR) & FDR < 0.05 & Effect > 0, na.rm = TRUE),
    N_low = sum(is.finite(FDR) & FDR < 0.05 & Effect < 0, na.rm = TRUE),
    N_sig = N_high + N_low,
    .groups = "drop"
  ) %>%
  arrange(Group_Order, desc(Mean_scaled_effect), desc(N_sig), Feature_Display)

row_order <- feature_rank$Feature_Display

plot_df <- checkpoint %>%
  mutate(
    Cancer = factor(Cancer, levels = col_order),
    Feature_Display = factor(Feature_Display, levels = rev(row_order))
  ) %>%
  arrange(Feature_Display, Cancer)

fwrite(
  plot_df,
  file.path(DATA_DIR, "Figure3E_checkpoint_biomarker_landscape_plot_data.csv")
)

fwrite(
  plot_df %>%
    transmute(
      Cancer,
      Feature = Feature_Display,
      Median_Low,
      Median_High,
      High_minus_Low_effect = Effect,
      FDR,
      Significance
    ),
  file.path(DATA_DIR, "Figure3E_checkpoint_biomarker_FDR_labels.txt"),
  sep = "\t"
)

fwrite(
  cancer_rank,
  file.path(DATA_DIR, "Figure3E_checkpoint_biomarker_cancer_order.csv")
)

fwrite(
  feature_rank,
  file.path(DATA_DIR, "Figure3E_checkpoint_biomarker_feature_order.csv")
)

p <- ggplot(plot_df, aes(x = Cancer, y = Feature_Display, fill = PlotValue)) +
  geom_tile(color = "white", linewidth = 0.50, width = 0.94, height = 0.90) +
  geom_text(
    aes(label = Significance),
    color = pal_dark,
    size = 2.45,
    fontface = "bold"
  ) +
  scale_fill_gradient2(
    low = pal_low,
    mid = pal_mid,
    high = pal_high,
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    name = "High - Low\nrelative\ndifference",
  ) +
  labs(
    title = "E  Checkpoint biomarker landscape",
    x = "Cancer type",
    y = NULL
  ) +
  theme_repair(9.5) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 7.2, face = "bold"),
    axis.text.y = element_text(size = 7.7, face = "bold"),
    panel.grid = element_blank(),
    legend.key.height = unit(25, "pt"),
    legend.key.width = unit(8, "pt"),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 19, hjust = 0),
    plot.subtitle = element_text(size = 9.5, color = "#718093", hjust = 0)
  )

save_plot(p, "Figure3E_checkpoint_biomarker_landscape", 9.6, 5.8)

legend_en <- c(
  "Figure 3E. Checkpoint biomarker landscape.",
  paste0(
    "Heatmap showing MO-DDRscore-high versus -low differences in immune checkpoint ",
    "and immunotherapy-related biomarkers across all 33 cancer types. The color scale ",
    "shows the direction and relative magnitude of marker-level differences between ",
    "MO-DDRscore-high and -low tumors; red denotes higher levels in the high-score group, ",
    "whereas teal denotes higher levels in the low-score group. Asterisks indicate FDR-adjusted ",
    "significance: * FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001, **** FDR < 0.0001."
  )
)

legend_cn <- c(
  "图3E. 免疫检查点标志物全景图。",
  paste0(
    "热图展示33种癌症中MO-DDRscore高低组在免疫检查点及免疫治疗相关标志物上的差异。",
    "颜色表示按每个标志物标准化后的High-Low表达差异；红色表示MO-DDRscore高分组更高，",
    "青绿色表示MO-DDRscore低分组更高。星号基于FDR校正：* FDR < 0.05，",
    "** FDR < 0.01，*** FDR < 0.001，**** FDR < 0.0001。"
  )
)

writeLines(c(legend_en, "", legend_cn), file.path(FIG_DIR, "Figure3E_legend.txt"))

cat("Done: Figure 3E checkpoint biomarker landscape script.\n")
cat("Output directory:\n", FIG_DIR, "\n")
