############################################################
# RepairDis Figure 3B
# TME cell-type remodeling
# - Main panel: FDR-driven representative dumbbell plot
# - Supplementary panel: full representative cell-type heatmap
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
SOURCE_DIR <- file.path(BASE_DIR, "03-res", "Figure2_pan_cancer", "G_TME_cell_types")
PLOT_DIR <- file.path(BASE_DIR, "03-res", "plots", "张岩-图")
FIG_DIR <- file.path(PLOT_DIR, "Figure3")
SUP_DIR <- file.path(PLOT_DIR, "sup", "Figure3")
DATA_DIR <- file.path(PLOT_DIR, "plot_data")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

unlink(file.path(SUP_DIR, c(
  "FigureS3B_TME_cell_type_full_representative_heatmap.pdf",
  "FigureS3B_TME_cell_type_full_representative_heatmap.png"
)))

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

N_MAIN_FEATURES <- 8
N_MAIN_CANCERS <- 8

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

clean_display <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    grepl("^b cells?$", x, ignore.case = TRUE) ~ "B cells",
    grepl("^t cells?$", x, ignore.case = TRUE) ~ "T cells",
    grepl("^cd4 tcells?$", x, ignore.case = TRUE) ~ "CD4 T cells",
    grepl("^cd8 tcells?$", x, ignore.case = TRUE) ~ "CD8 T cells",
    TRUE ~ x
  )
}

feature_label <- function(display, method) {
  paste0(clean_display(display), "\n(", method, ")")
}

main_feature_label <- function(display) {
  clean_display(display)
}

read_estimate_main_cancers <- function() {
  f1 <- file.path(DATA_DIR, "Figure3A_ESTIMATE_main_kept_cancers.csv")
  f2 <- file.path(DATA_DIR, "Figure3A_ESTIMATE_FDR_significant_heatmap_data.csv")

  if (file.exists(f1)) {
    x <- fread(f1, data.table = FALSE, colClasses = "character")
    if (all(c("Cancer", "Kept_in_main") %in% colnames(x))) {
      keep <- tolower(as.character(x$Kept_in_main)) %in% c("true", "t", "1", "yes")
      out <- x$Cancer[keep]
      out <- out[out %in% cancer_order]
      if (length(out) > 0) return(out)
    }
  }

  if (file.exists(f2)) {
    x <- fread(f2, data.table = FALSE, colClasses = "character")
    if ("Cancer" %in% colnames(x)) {
      out <- unique(x$Cancer)
      out <- out[out %in% cancer_order]
      if (length(out) > 0) return(out)
    }
  }

  cancer_order
}

source_file <- file.path(SOURCE_DIR, "PanCancer_TME_cell_types_group_comparison_all.csv")
if (!file.exists(source_file)) {
  stop("Cannot find input file: ", source_file)
}

stat_df <- fread(source_file, data.table = FALSE, colClasses = "character") %>%
  mutate(
    Cancer = as.character(Cancer),
    Feature = as.character(Feature),
    Display = clean_display(Display),
    Method = as.character(Method),
    Canonical = as.character(Canonical),
    Median_Low = as.numeric(Median_Low),
    Median_High = as.numeric(Median_High),
    Effect = as.numeric(Effect),
    P = as.numeric(P),
    FDR = as.numeric(FDR)
  ) %>%
  filter(
    Cancer %in% cancer_order,
    is.finite(Median_Low),
    is.finite(Median_High),
    !is.na(Feature),
    !is.na(Canonical)
  )

if (nrow(stat_df) == 0) {
  stop("No TME cell-type statistics were loaded.")
}

median_long <- bind_rows(
  stat_df %>%
    transmute(Cancer, Feature, Display, Method, Source, Canonical, Group = "Low", Median = Median_Low),
  stat_df %>%
    transmute(Cancer, Feature, Display, Method, Source, Canonical, Group = "High", Median = Median_High)
) %>%
  group_by(Feature) %>%
  mutate(StandardizedMedian = safe_scale(Median)) %>%
  ungroup()

plot_base <- median_long %>%
  select(Cancer, Feature, Display, Method, Source, Canonical, Group, StandardizedMedian) %>%
  pivot_wider(names_from = Group, values_from = StandardizedMedian) %>%
  left_join(
    stat_df %>%
      select(Cancer, Feature, Display, Method, Source, Canonical, P, FDR),
    by = c("Cancer", "Feature", "Display", "Method", "Source", "Canonical")
  ) %>%
  mutate(
    Low = as.numeric(Low),
    High = as.numeric(High),
    Z_Effect = High - Low,
    PlotEffect = cap(Z_Effect, 2.5),
    Significance = sig_from_fdr(FDR, ns_blank = TRUE),
    Direction = ifelse(Z_Effect >= 0, "Higher in High", "Higher in Low"),
    Feature_Label = feature_label(Display, Method),
    Feature_Main_Label = main_feature_label(Display)
  )

# One representative output is retained for each canonical cell type. Ranking is
# based on the number of FDR-significant cancer types, then standardized effect.
feature_rank_all <- plot_base %>%
  group_by(Canonical, Feature, Display, Method, Source, Feature_Label, Feature_Main_Label) %>%
  summarise(
    N_FDR_sig_cancers = sum(is.finite(FDR) & FDR < 0.05, na.rm = TRUE),
    Median_abs_standardized_effect = median(abs(Z_Effect), na.rm = TRUE),
    Mean_abs_standardized_effect = mean(abs(Z_Effect), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Canonical, desc(N_FDR_sig_cancers), desc(Median_abs_standardized_effect)) %>%
  group_by(Canonical) %>%
  mutate(Representative_for_canonical = row_number() == 1) %>%
  ungroup()

representative_rank <- feature_rank_all %>%
  filter(Representative_for_canonical) %>%
  arrange(desc(N_FDR_sig_cancers), desc(Median_abs_standardized_effect), Feature_Label)

selected_features <- representative_rank %>%
  slice_head(n = N_MAIN_FEATURES) %>%
  mutate(Selected_for_main = TRUE)

selected_feature_ids <- selected_features$Feature

estimate_main_cancers <- read_estimate_main_cancers()

selected_cancers <- plot_base %>%
  filter(Feature %in% selected_feature_ids, Cancer %in% estimate_main_cancers) %>%
  group_by(Cancer) %>%
  summarise(
    N_FDR_sig_selected_features = sum(is.finite(FDR) & FDR < 0.05, na.rm = TRUE),
    N_Higher_in_High = sum(is.finite(FDR) & FDR < 0.05 & Z_Effect > 0, na.rm = TRUE),
    N_Higher_in_Low = sum(is.finite(FDR) & FDR < 0.05 & Z_Effect < 0, na.rm = TRUE),
    Direction_balance = N_Higher_in_High - N_Higher_in_Low,
    Cell_remodeling_strength = sum(ifelse(is.finite(FDR) & FDR < 0.05, abs(Z_Effect), 0), na.rm = TRUE),
    Signed_remodeling_score = sum(ifelse(is.finite(FDR) & FDR < 0.05, Z_Effect, 0), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(N_FDR_sig_selected_features), desc(Cell_remodeling_strength), Cancer) %>%
  slice_head(n = N_MAIN_CANCERS)

selected_cancer_order <- selected_cancers %>%
  arrange(
    desc(Signed_remodeling_score),
    desc(Direction_balance),
    desc(N_FDR_sig_selected_features),
    desc(Cell_remodeling_strength),
    Cancer
  ) %>%
  pull(Cancer)

main_df <- plot_base %>%
  filter(Feature %in% selected_feature_ids, Cancer %in% selected_cancer_order) %>%
  left_join(
    selected_features %>%
      select(Feature, Feature_rank = N_FDR_sig_cancers, Feature_strength = Median_abs_standardized_effect),
    by = "Feature"
  ) %>%
  mutate(
    Cancer = factor(Cancer, levels = rev(selected_cancer_order)),
    Feature_Main_Label = factor(
      Feature_Main_Label,
      levels = selected_features$Feature_Main_Label
    ),
    Feature_Label = factor(
      Feature_Label,
      levels = selected_features$Feature_Label
    ),
    Significance = sig_from_fdr(FDR, ns_blank = TRUE),
    Direction = factor(Direction, levels = c("Higher in High", "Higher in Low"))
  )

fwrite(
  feature_rank_all %>%
    mutate(Selected_for_main = Feature %in% selected_feature_ids) %>%
    arrange(desc(Representative_for_canonical), desc(N_FDR_sig_cancers), desc(Median_abs_standardized_effect)),
  file.path(DATA_DIR, "Figure3B_TME_cell_type_feature_selection_ranking.csv")
)

fwrite(
  selected_features,
  file.path(DATA_DIR, "Figure3B_TME_cell_type_selected_features.csv")
)

fwrite(
  selected_cancers,
  file.path(DATA_DIR, "Figure3B_TME_cell_type_selected_cancers.csv")
)

fwrite(
  main_df,
  file.path(DATA_DIR, "Figure3B_TME_cell_type_main_plot_data.csv")
)

x_min <- -2.5
x_max <- 4.2
star_x <- 4.0

main_point_df <- bind_rows(
  main_df %>%
    transmute(
      Cancer, Feature_Main_Label,
      Score = Low,
      Group = "Low"
    ),
  main_df %>%
    transmute(
      Cancer, Feature_Main_Label,
      Score = High,
      Group = "High"
    )
) %>%
  mutate(
    Group = factor(Group, levels = c("Low", "High")),
    Cancer = factor(Cancer, levels = levels(main_df$Cancer)),
    Feature_Main_Label = factor(Feature_Main_Label, levels = levels(main_df$Feature_Main_Label))
  )

p_main <- ggplot(main_df, aes(y = Cancer)) +
  geom_vline(xintercept = 0, color = "#AEBAC8", linewidth = 0.55, linetype = 2) +
  geom_segment(
    aes(x = Low, xend = High, yend = Cancer),
    color = "#8D949C",
    linewidth = 0.68,
    alpha = 0.86,
    lineend = "round"
  ) +
  geom_point(
    data = main_point_df,
    aes(x = Score, y = Cancer, color = Group),
    inherit.aes = FALSE,
    size = 2.25,
    alpha = 0.96
  ) +
  geom_text(
    aes(x = star_x, label = Significance),
    color = pal_dark,
    fontface = "bold",
    size = 3.45,
    hjust = 0.5
  ) +
  facet_wrap(~ Feature_Main_Label, ncol = 4) +
  scale_color_manual(
    values = c("Low" = pal_low, "High" = pal_high),
    breaks = c("Low", "High"),
    name = "MO-DDRscore"
  ) +
  scale_x_continuous(breaks = c(-2.5, 0, 2.5, 4), limits = c(x_min, x_max)) +
  labs(
    title = "B  TME cell-type remodeling",
    x = "Standardized median cell score",
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 21, color = pal_dark, hjust = 0),
    plot.subtitle = element_text(size = 9.8, color = "#718093", hjust = 0),
    strip.background = element_rect(fill = "#EEF3F7", color = "#D5DEE8", linewidth = 0.5),
    strip.text = element_text(face = "bold", color = pal_dark, size = 10.5, lineheight = 0.95),
    panel.grid.major = element_line(color = pal_grid, linewidth = 0.45),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "#CAD6E2", fill = NA, linewidth = 0.65),
    axis.text = element_text(color = "#2D3E50", face = "bold"),
    axis.text.x = element_text(size = 9.2),
    axis.text.y = element_text(size = 9.2),
    axis.title.x = element_text(face = "bold", color = pal_dark, size = 12),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold", color = pal_dark),
    legend.text = element_text(color = "#34495E"),
    plot.margin = margin(8, 22, 8, 8)
  )

ggsave(file.path(FIG_DIR, "Figure3B_TME_cell_type.pdf"), p_main, width = 12.8, height = 7.2)
ggsave(file.path(FIG_DIR, "Figure3B_TME_cell_type.png"), p_main, width = 12.8, height = 7.2, dpi = 600)
ggsave(file.path(FIG_DIR, "Figure3B_TME_cell_type.tiff"), p_main, width = 12.8, height = 7.2, dpi = 600, compression = "lzw")

############################################################
# Supplementary heatmap: all representative cell-type features
############################################################

rep_feature_ids <- representative_rank$Feature

full_heat_df <- plot_base %>%
  filter(Feature %in% rep_feature_ids) %>%
  mutate(
    PlotEffect = cap(Z_Effect, 2.5),
    OrderEffect = ifelse(is.finite(FDR) & FDR < 0.05, PlotEffect, 0),
    Significance = sig_from_fdr(FDR, ns_blank = TRUE)
  )

make_order <- function(df, id_col, value_col = "OrderEffect") {
  id_col <- rlang::ensym(id_col)
  mat_df <- df %>%
    select(Cancer, !!id_col, all_of(value_col)) %>%
    mutate(ID = as.character(!!id_col)) %>%
    select(Cancer, ID, all_of(value_col)) %>%
    pivot_wider(names_from = Cancer, values_from = all_of(value_col), values_fill = 0)
  mat <- as.matrix(mat_df[, setdiff(colnames(mat_df), "ID"), drop = FALSE])
  rownames(mat) <- mat_df$ID
  if (nrow(mat) <= 1) return(rownames(mat))
  rownames(mat)[hclust(dist(mat), method = "average")$order]
}

heat_feature_order <- make_order(full_heat_df, Feature_Label)
heat_cancer_order <- {
  mat_df <- full_heat_df %>%
  select(Cancer, Feature_Label, PlotEffect) %>%
    left_join(
      full_heat_df %>% select(Cancer, Feature_Label, OrderEffect),
      by = c("Cancer", "Feature_Label")
    ) %>%
    select(Cancer, Feature_Label, OrderEffect) %>%
    pivot_wider(names_from = Cancer, values_from = OrderEffect, values_fill = 0)
  mat <- as.matrix(mat_df[, setdiff(colnames(mat_df), "Feature_Label"), drop = FALSE])
  if (ncol(mat) <= 1) {
    colnames(mat)
  } else {
    hc <- hclust(dist(t(mat)), method = "average")
    cancer_weight <- colMeans(mat, na.rm = TRUE)
    dend <- reorder(as.dendrogram(hc), wts = cancer_weight, agglo.FUN = mean)
    labels(dend)
  }
}

full_heat_df <- full_heat_df %>%
  mutate(
    Cancer = factor(Cancer, levels = heat_cancer_order),
    Feature_Label = factor(Feature_Label, levels = rev(heat_feature_order))
  )

fwrite(
  full_heat_df,
  file.path(DATA_DIR, "FigureS3C_TME_cell_type_full_representative_heatmap_data.csv")
)

p_full_heat <- ggplot(full_heat_df, aes(x = Cancer, y = Feature_Label, fill = PlotEffect)) +
  geom_tile(color = "white", linewidth = 0.48) +
  geom_text(aes(label = Significance), color = pal_dark, fontface = "bold", size = 2.55) +
  scale_fill_gradient2(
    low = pal_low,
    mid = pal_mid,
    high = pal_high,
    midpoint = 0,
    limits = c(-2.5, 2.5),
    breaks = c(-2, -1, 0, 1, 2),
    name = "High - Low\ndifference"
  ) +
  labs(
    title = "Figure S3C  Full TME cell-type remodeling landscape",
    x = "Cancer type",
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 17, color = pal_dark, hjust = 0),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, face = "bold", color = "#2D3E50"),
    axis.text.y = element_text(face = "bold", color = "#2D3E50", size = 8.2),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "#CAD6E2", fill = NA, linewidth = 0.65),
    legend.title = element_text(face = "bold", color = pal_dark),
    plot.margin = margin(8, 10, 8, 8)
  )

ggsave(file.path(SUP_DIR, "FigureS3C_TME_cell_type_full_representative_heatmap.pdf"), p_full_heat, width = 13.5, height = 8.8)
ggsave(file.path(SUP_DIR, "FigureS3C_TME_cell_type_full_representative_heatmap.png"), p_full_heat, width = 13.5, height = 8.8, dpi = 600)

rule_df <- data.frame(
  Item = c(
    "Input_file",
    "Main_feature_selection",
    "Main_cancer_pool",
    "Main_cancer_selection",
    "Significance_label",
    "Displayed_value"
  ),
  Value = c(
    source_file,
    "One representative output per canonical cell type; ranked by FDR-significant cancer count and median absolute standardized effect; top 8 shown.",
    "Cancer types retained in Figure 3A ESTIMATE heatmap when available; otherwise all cancer types.",
    "Top 8 cancer types ranked by FDR-significant selected features and total absolute standardized effect.",
    "Asterisks are generated from FDR: * < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001.",
    "Feature-wise standardized group median score for Low and High MO-DDRscore groups."
  ),
  stringsAsFactors = FALSE
)

fwrite(rule_df, file.path(DATA_DIR, "Figure3B_TME_cell_type_significance_and_selection_rule.csv"))

cat("\nDone.\n")
cat("Main figure:\n", file.path(FIG_DIR, "Figure3B_TME_cell_type.pdf"), "\n")
cat("Supplementary figures:\n", SUP_DIR, "\n")
cat("Plot data:\n", DATA_DIR, "\n")
