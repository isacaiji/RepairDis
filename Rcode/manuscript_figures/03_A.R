############################################################
# RepairDis Figure 3A
# ESTIMATE-based tumor microenvironment profiles
# - Main panel: 33 cancer type heatmap
# - Supplementary panel: full pan-cancer violin plot
# - All significance labels are recalculated from FDR
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

BASE_DIR <- "D:/R_workspace/summary/Repairium_pancancer_analysis_plan"
SOURCE_DIR <- file.path(BASE_DIR, "03-res", "Figure2_pan_cancer", "F_ESTIMATE_TME")
PLOT_DIR <- file.path(BASE_DIR, "03-res", "plots", "\u5f20\u5ca9-\u56fe")
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

feature_order <- c("ESTIMATEScore", "StromalScore", "ImmuneScore", "TumorPurity")

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

clean_estimate_feature <- function(x) {
  dplyr::case_when(
    grepl("^ESTIMATEScore", x, ignore.case = TRUE) ~ "ESTIMATEScore",
    grepl("^StromalScore", x, ignore.case = TRUE) ~ "StromalScore",
    grepl("^ImmuneScore", x, ignore.case = TRUE) ~ "ImmuneScore",
    grepl("^TumorPurity", x, ignore.case = TRUE) ~ "TumorPurity",
    TRUE ~ as.character(x)
  )
}

read_one_stat <- function(cancer) {
  file <- file.path(SOURCE_DIR, cancer, paste0(cancer, "_IOBR_estimate_group_comparison.csv"))
  if (!file.exists(file)) return(NULL)
  fread(file, data.table = FALSE, colClasses = "character") %>%
    mutate(Cancer = cancer)
}

stat_df <- bind_rows(lapply(cancer_order, read_one_stat)) %>%
  mutate(
    Feature_Display = clean_estimate_feature(Feature),
    Effect = as.numeric(Effect),
    FDR = as.numeric(FDR),
    P = as.numeric(P)
  ) %>%
  filter(Cancer %in% cancer_order, Feature_Display %in% feature_order)

if (nrow(stat_df) == 0) {
  stop("No ESTIMATE statistics were loaded. Please check SOURCE_DIR.")
}

heat_df <- stat_df %>%
  group_by(Feature_Display) %>%
  mutate(
    EffectScaledWithinFeature = as.numeric(scale(Effect)),
    PlotValue = cap(EffectScaledWithinFeature, 2.5)
  ) %>%
  ungroup() %>%
  mutate(
    Significance = sig_from_fdr(FDR),
    Feature_Display = factor(Feature_Display, levels = rev(feature_order))
  )

get_col_order <- function(df) {
  mat_df <- df %>%
    select(Cancer, Feature_Display, PlotValue) %>%
    mutate(Feature_Display = as.character(Feature_Display)) %>%
    pivot_wider(names_from = Cancer, values_from = PlotValue, values_fill = 0)
  mat <- as.matrix(mat_df[, setdiff(colnames(mat_df), "Feature_Display"), drop = FALSE])
  rownames(mat) <- mat_df$Feature_Display
  out <- colnames(mat)
  if (ncol(mat) > 1) {
    out <- colnames(mat)[hclust(dist(t(mat)), method = "average")$order]
  }
  out
}

get_tme_depletion_order <- function(df) {
  score_df <- df %>%
    select(Cancer, Feature_Display, PlotValue) %>%
    mutate(
      Feature_Display = as.character(Feature_Display),
      PlotValue = as.numeric(PlotValue)
    ) %>%
    pivot_wider(names_from = Feature_Display, values_from = PlotValue, values_fill = 0)

  for (nm in feature_order) {
    if (!nm %in% colnames(score_df)) score_df[[nm]] <- 0
  }

  score_df %>%
    mutate(
      TME_DepletionScore =
        -ESTIMATEScore - StromalScore - ImmuneScore + TumorPurity,
      ImmuneStromalDepletion =
        -ESTIMATEScore - StromalScore - ImmuneScore,
      PurityIncrease = TumorPurity
    ) %>%
    arrange(
      desc(TME_DepletionScore),
      desc(ImmuneStromalDepletion),
      desc(PurityIncrease),
      Cancer
    ) %>%
    pull(Cancer)
}

full_col_order <- get_tme_depletion_order(heat_df)

sig_cancers <- heat_df %>%
  group_by(Cancer) %>%
  summarise(Any_FDR_significant = any(is.finite(FDR) & FDR < 0.05), .groups = "drop") %>%
  filter(Any_FDR_significant) %>%
  pull(Cancer)

heat_main_df <- heat_df %>%
  filter(Cancer %in% sig_cancers)

# Reorder displayed cancer types by a biologically interpretable TME-depletion
# pattern: lower ESTIMATE/stromal/immune scores plus higher tumor purity in the
# MO-DDRscore-high group are placed toward the left.
col_order <- get_tme_depletion_order(heat_main_df)

heat_main_df <- heat_main_df %>%
  mutate(Cancer = factor(Cancer, levels = col_order))

fwrite(
  heat_df %>%
    mutate(Cancer = factor(Cancer, levels = full_col_order)) %>%
    arrange(Feature_Display, Cancer),
  file.path(DATA_DIR, "Figure3A_ESTIMATE_full33_heatmap_data.csv")
)

fwrite(
  heat_main_df %>%
    arrange(Feature_Display, Cancer),
  file.path(DATA_DIR, "Figure3A_ESTIMATE_FDR_significant_heatmap_data.csv")
)

fwrite(
  data.frame(
    Cancer = cancer_order,
    Kept_in_main = cancer_order %in% sig_cancers,
    stringsAsFactors = FALSE
  ),
  file.path(DATA_DIR, "Figure3A_ESTIMATE_main_kept_cancers.csv")
)

fwrite(
  data.frame(
    Panel = "Figure3A",
    Significance_rule = "Asterisks were recalculated from FDR: * <0.05, ** <0.01, *** <0.001, **** <0.0001.",
    Main_filter = "Cancer types with no FDR-significant ESTIMATE feature were omitted from the main heatmap.",
    stringsAsFactors = FALSE
  ),
  file.path(DATA_DIR, "Figure3A_ESTIMATE_significance_rule.csv")
)

# Legacy name kept for downstream PowerPoint links.
fwrite(
  heat_main_df %>%
    arrange(Feature_Display, Cancer),
  file.path(DATA_DIR, "Figure3A_ESTIMATE_full_heatmap_data.csv")
)

p_heat <- ggplot(heat_main_df, aes(x = Cancer, y = Feature_Display, fill = PlotValue)) +
  geom_tile(color = "white", linewidth = 0.82) +
  geom_text(aes(label = Significance), color = pal_dark, fontface = "bold", size = 3.45) +
  scale_fill_gradient2(
    low = pal_low,
    mid = pal_mid,
    high = pal_high,
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    name = "High - Low
 difference"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title = "A  ESTIMATE",
    x = "Cancer type",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 23, color = pal_dark, hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = pal_dark, face = "bold", size = 10.5),
    axis.text.y = element_text(color = pal_dark, face = "bold", size = 12),
    axis.title.x = element_text(color = pal_dark, face = "bold", size = 13, margin = margin(t = 8)),
    panel.grid = element_blank(),
    legend.title = element_text(face = "bold", color = pal_dark, size = 11),
    legend.text = element_text(color = pal_dark, size = 10),
    legend.key.height = unit(30, "pt"),
    legend.key.width = unit(10, "pt"),
    plot.margin = margin(8, 14, 8, 8)
  )

main_width <- max(9.2, length(col_order) * 0.42 + 2.8)
ggsave(file.path(FIG_DIR, "Figure3A_ESTIMAT.pdf"), p_heat, width = main_width, height = 4.4)
ggsave(file.path(FIG_DIR, "Figure3A_ESTIMATE.png"), p_heat, width = main_width, height = 4.4, dpi = 600)
ggsave(file.path(FIG_DIR, "Figure3A_ESTIMATE.tiff"), p_heat, width = main_width, height = 4.4, dpi = 600, compression = "lzw")

############################################################
# Supplementary full violin plot
############################################################

read_one_long <- function(cancer) {
  file <- file.path(SOURCE_DIR, cancer, paste0(cancer, "_F_ESTIMATE_TME_plot_df_used.csv"))
  if (!file.exists(file)) {
    file <- file.path(SOURCE_DIR, cancer, paste0(cancer, "_IOBR_estimate_long.csv"))
  }
  if (!file.exists(file)) return(NULL)
  fread(file, data.table = FALSE, colClasses = "character") %>%
    mutate(Cancer = cancer)
}

# Full violin keeps all 33 cancers but places cancers with available
# ESTIMATE statistics according to the same effect-pattern reorder.
violin_col_order <- c(full_col_order, setdiff(cancer_order, full_col_order))

long_df <- bind_rows(lapply(cancer_order, read_one_long)) %>%
  mutate(
    Feature_Display = clean_estimate_feature(Feature),
    Group = factor(Group, levels = c("Low", "High")),
    Cancer = factor(Cancer, levels = violin_col_order),
    Value = as.numeric(Value)
  ) %>%
  filter(Feature_Display %in% feature_order, Group %in% c("Low", "High"), is.finite(Value), !is.na(Cancer))

stat_for_star <- heat_df %>%
  mutate(
    Feature_Display = as.character(Feature_Display),
    Cancer = as.character(Cancer)
  ) %>%
  select(Cancer, Feature_Display, Significance)

star_df <- long_df %>%
  group_by(Cancer, Feature_Display) %>%
  summarise(
    y = max(Value, na.rm = TRUE),
    span = diff(range(Value, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(y = y + ifelse(is.finite(span) & span > 0, span * 0.08, abs(y) * 0.08 + 0.1)) %>%
  left_join(stat_for_star, by = c("Cancer", "Feature_Display")) %>%
  filter(Significance != "") %>%
  group_by(Feature_Display) %>%
  mutate(y = max(y, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    Cancer = factor(Cancer, levels = violin_col_order),
    Feature_Display = factor(Feature_Display, levels = rev(feature_order))
  )

long_df <- long_df %>%
  mutate(
    Cancer = factor(Cancer, levels = violin_col_order),
    Feature_Display = factor(Feature_Display, levels = rev(feature_order))
  )

fwrite(long_df, file.path(DATA_DIR, "FigureS3A_ESTIMATE_full_violin_data.csv"))
fwrite(star_df, file.path(DATA_DIR, "FigureS3A_ESTIMATE_full_violin_FDR_stars.csv"))

make_violin_plot <- function(plot_df, star_plot_df, title_text, x_levels, width_out) {
  plot_df <- plot_df %>%
    mutate(
      Cancer = factor(as.character(Cancer), levels = x_levels),
      Feature_Display = factor(as.character(Feature_Display), levels = rev(feature_order))
    ) %>%
    filter(!is.na(Cancer), !is.na(Feature_Display))

  star_plot_df <- star_plot_df %>%
    mutate(
      Cancer = factor(as.character(Cancer), levels = x_levels),
      Feature_Display = factor(as.character(Feature_Display), levels = rev(feature_order))
    ) %>%
    filter(!is.na(Cancer), !is.na(Feature_Display))

  ggplot(plot_df, aes(x = Cancer, y = Value, fill = Group)) +
  geom_violin(
    position = position_dodge(width = 0.78),
    width = 0.72,
    trim = TRUE,
    color = "white",
    linewidth = 0.18,
    alpha = 0.86
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.78),
    width = 0.13,
    outlier.shape = NA,
    color = pal_dark,
    linewidth = 0.25,
    alpha = 0.72
  ) +
  geom_text(
    data = star_plot_df,
    aes(x = Cancer, y = y, label = Significance),
    inherit.aes = FALSE,
    color = pal_dark,
    fontface = "bold",
    size = 3.0
  ) +
  facet_grid(Feature_Display ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = c(Low = pal_low, High = pal_high), name = "MO-DDRscore") +
  scale_x_discrete(expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.14))) +
  labs(
    title = title_text,
    x = "Cancer type",
    y = "ESTIMATE feature value"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 20, color = pal_dark, hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = pal_dark, face = "bold", size = 8.8),
    axis.text.y = element_text(color = pal_dark, size = 8.8),
    axis.title = element_text(color = pal_dark, face = "bold", size = 12),
    strip.background = element_rect(fill = "#EDF2F7", color = "#D7E0EA", linewidth = 0.45),
    strip.text.y.left = element_text(face = "bold", color = pal_dark, angle = 0, size = 10.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "#AEBBCC", fill = NA, linewidth = 0.55),
    legend.position = "top",
    legend.title = element_text(face = "bold", color = pal_dark),
    legend.text = element_text(color = pal_dark),
    plot.margin = margin(8, 12, 8, 8)
  )
}

p_violin <- make_violin_plot(
  long_df,
  star_df,
  "Figure S3A  Full ESTIMATE-based tumor microenvironment profiles",
  violin_col_order,
  15.8
)

ggsave(file.path(SUP_DIR, "FigureS3A_ESTIMATE_full_violin.pdf"), p_violin, width = 15.8, height = 8.8)
ggsave(file.path(SUP_DIR, "FigureS3A_ESTIMATE_full_violin.png"), p_violin, width = 15.8, height = 8.8, dpi = 600)

sig_violin_col_order <- col_order
sig_long_df <- long_df %>%
  mutate(Cancer = as.character(Cancer)) %>%
  filter(Cancer %in% sig_violin_col_order) %>%
  mutate(Cancer = factor(Cancer, levels = sig_violin_col_order))

sig_star_df <- star_df %>%
  mutate(Cancer = as.character(Cancer)) %>%
  filter(Cancer %in% sig_violin_col_order) %>%
  mutate(Cancer = factor(Cancer, levels = sig_violin_col_order))

fwrite(sig_long_df, file.path(DATA_DIR, "FigureS3B_ESTIMATE_FDR_significant_violin_data.csv"))
fwrite(sig_star_df, file.path(DATA_DIR, "FigureS3B_ESTIMATE_FDR_significant_violin_FDR_stars.csv"))

p_sig_violin <- make_violin_plot(
  sig_long_df,
  sig_star_df,
  "Figure S3B  FDR-significant ESTIMATE-based tumor microenvironment profiles",
  sig_violin_col_order,
  10.8
)

ggsave(file.path(SUP_DIR, "FigureS3B_ESTIMATE_FDR_significant_violin.pdf"), p_sig_violin, width = 10.8, height = 8.8)
ggsave(file.path(SUP_DIR, "FigureS3B_ESTIMATE_FDR_significant_violin.png"), p_sig_violin, width = 10.8, height = 8.8, dpi = 600)

cat("Done.\n")
cat("Main figure:\n", file.path(FIG_DIR, "Figure3A_ESTIMATE_full_heatmap.pdf"), "\n")
cat("Supplementary figure:\n", file.path(SUP_DIR, "FigureS3A_ESTIMATE_full_violin.pdf"), "\n")
cat("Supplementary significant-only figure:\n", file.path(SUP_DIR, "FigureS3B_ESTIMATE_FDR_significant_violin.pdf"), "\n")
