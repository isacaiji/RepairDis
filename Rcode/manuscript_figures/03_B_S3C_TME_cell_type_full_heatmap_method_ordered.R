############################################################
# Supplementary Figure S3C
# TME cell-type full heatmap, ordered by deconvolution method
# This standalone script does not overwrite the existing S3C file.
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

BASE_DIR <- "D:/R_workspace/summary/Repairium_pancancer_analysis_plan"
PLOT_DIR <- file.path(BASE_DIR, "03-res", "plots", "\u5f20\u5ca9-\u56fe")
SOURCE_DIR <- file.path(BASE_DIR, "03-res", "Figure2_pan_cancer", "G_TME_cell_types")
DATA_DIR <- file.path(PLOT_DIR, "plot_data")
OUT_DIR <- file.path(PLOT_DIR, "sup", "manuscript_ready", "FigureS3")

dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Change this line if a different method order is preferred.
METHOD_ORDER <- c("xCell", "MCP-counter", "EPIC", "quanTIseq")

cancer_order <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
  "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
  "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
)

pal_low <- "#20AEB3"
pal_high <- "#D65A67"
pal_dark <- "#10243C"

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
    Source = as.character(Source),
    Canonical = as.character(Canonical),
    Median_Low = as.numeric(Median_Low),
    Median_High = as.numeric(Median_High),
    P = as.numeric(P),
    FDR = as.numeric(FDR)
  ) %>%
  filter(
    Cancer %in% cancer_order,
    Method %in% METHOD_ORDER,
    is.finite(Median_Low),
    is.finite(Median_High),
    !is.na(Display),
    !is.na(Canonical)
  )

median_long <- bind_rows(
  stat_df %>%
    transmute(Cancer, Feature, Display, Method, Source, Canonical,
              Group = "Low", Median = Median_Low),
  stat_df %>%
    transmute(Cancer, Feature, Display, Method, Source, Canonical,
              Group = "High", Median = Median_High)
) %>%
  group_by(Feature) %>%
  mutate(StandardizedMedian = safe_scale(Median)) %>%
  ungroup()

plot_df <- median_long %>%
  select(Cancer, Feature, Display, Method, Source, Canonical, Group, StandardizedMedian) %>%
  pivot_wider(names_from = Group, values_from = StandardizedMedian) %>%
  left_join(
    stat_df %>% select(Cancer, Feature, Display, Method, Source, Canonical, P, FDR),
    by = c("Cancer", "Feature", "Display", "Method", "Source", "Canonical")
  ) %>%
  mutate(
    Method = factor(Method, levels = METHOD_ORDER),
    Cancer = factor(Cancer, levels = cancer_order),
    Z_Effect = High - Low,
    PlotEffect = cap(Z_Effect, 2.5),
    Significance = sig_from_fdr(FDR, ns_blank = TRUE),
    MethodDisplay = paste0(Display, "\n(", as.character(Method), ")")
  )

feature_order <- plot_df %>%
  distinct(Method, Display, MethodDisplay) %>%
  mutate(Method = factor(Method, levels = METHOD_ORDER)) %>%
  arrange(Method, Display) %>%
  pull(MethodDisplay) %>%
  unique()

plot_df <- plot_df %>%
  mutate(MethodDisplay = factor(MethodDisplay, levels = rev(feature_order)))

fwrite(
  plot_df %>%
    arrange(Method, Display, Cancer) %>%
    select(
      Cancer, Feature, Display, Method, Source, Canonical,
      Low, High, P, FDR, Z_Effect, PlotEffect, Significance, MethodDisplay
    ),
  file.path(DATA_DIR, "FigureS3C_TME_cell_type_full_heatmap_method_ordered_data.csv")
)

p <- ggplot(plot_df, aes(x = Cancer, y = MethodDisplay, fill = PlotEffect)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = Significance), color = pal_dark, fontface = "bold", size = 2.25) +
  facet_grid(Method ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(
    low = pal_low,
    mid = "#F8F3EE",
    high = pal_high,
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    name = "High - Low\ndifference"
  ) +
  labs(
    title = "TME cell-type remodeling by deconvolution method",
    x = "Cancer type",
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.y = unit(0.10, "lines"),
    strip.placement = "outside",
    strip.background.y = element_rect(fill = "#EDF3F8", color = "#D5DEE8", linewidth = 0.45),
    strip.text.y.left = element_text(face = "bold", color = pal_dark, angle = 0, size = 8.5),
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 7.0, face = "bold"),
    axis.text.y = element_text(size = 6.3, face = "bold", color = "#2D3E50", lineheight = 0.78),
    axis.title.x = element_text(face = "bold", color = pal_dark, size = 9.5),
    legend.title = element_text(face = "bold", color = pal_dark, size = 8.5),
    legend.text = element_text(color = pal_dark, size = 7.5),
    plot.title = element_text(face = "bold", color = pal_dark, size = 16, hjust = 0),
    plot.margin = margin(6, 8, 6, 6)
  )

out_stem <- "FigureS3C_TME_cell_type_full_heatmap_method_ordered"

ggsave(file.path(OUT_DIR, paste0(out_stem, ".pdf")), p, width = 13.8, height = 9.8, device = cairo_pdf)
ggsave(file.path(OUT_DIR, paste0(out_stem, ".png")), p, width = 13.8, height = 9.8, dpi = 600, bg = "white")
ggsave(
  file.path(OUT_DIR, paste0(out_stem, ".tiff")),
  p, width = 13.8, height = 9.8, dpi = 600, compression = "lzw", bg = "white"
)

cat("Done: method-ordered Supplementary Figure S3C\n")
cat("Output directory:", OUT_DIR, "\n")
cat("Method order:", paste(METHOD_ORDER, collapse = " -> "), "\n")
