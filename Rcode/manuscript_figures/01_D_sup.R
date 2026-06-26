############################################################
# RepairDis Figure S1D: full cancer hallmark landscape
# - Full 50 MSigDB Hallmark programs across 33 cancer types
# - Data source: Figure2_pan_cancer/C_Cancer_hallmark
# - Asterisks are generated from FDR, not raw P
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(grid)
})

BASE_DIR <- "D:/R_workspace/summary/Repairium_pancancer_analysis_plan"
RAW_DIR <- file.path(BASE_DIR, "03-res/Figure2_pan_cancer/C_Cancer_hallmark")
PLOT_DIR <- file.path(BASE_DIR, "03-res/plots/张岩-图")
SUP_DIR <- file.path(PLOT_DIR, "sup", "Figure1")
DATA_OUT_DIR <- file.path(SUP_DIR, "plot_data")

dir.create(SUP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

pal_low <- "#20AEB3"
pal_high <- "#D45F5F"
pal_high2 <- "#B23A48"
pal_dark <- "#10243C"
pal_grid <- "#E7EDF3"

cancer_order <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
  "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
  "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
)

nice_feature <- function(x) {
  x <- gsub("^HALLMARK_", "", as.character(x))
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  x <- gsub("Dna", "DNA", x)
  x <- gsub("E2f", "E2F", x)
  x <- gsub("G2m", "G2M", x)
  x <- gsub("Tgf", "TGF", x)
  x <- gsub("Mtorc1", "MTORC1", x)
  x <- gsub("Myc", "MYC", x)
  x <- gsub("Pi3k Akt Mtor", "PI3K-AKT-MTOR", x)
  x
}

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

theme_repair <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(color = pal_dark),
      plot.title = element_text(face = "bold", size = base_size + 5, hjust = 0),
      plot.subtitle = element_text(size = base_size + 0.5, color = "#617083", hjust = 0),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "#CAD5E0", fill = NA, linewidth = 0.45),
      axis.text = element_text(color = "#263849"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(color = "#263849"),
      plot.margin = margin(8, 14, 8, 8)
    )
}

save_plot <- function(p, stem, width, height) {
  pdf_file <- file.path(SUP_DIR, paste0(stem, ".pdf"))
  png_file <- file.path(SUP_DIR, paste0(stem, ".png"))
  tiff_file <- file.path(SUP_DIR, paste0(stem, ".tiff"))

  ggsave(pdf_file, p, width = width, height = height, useDingbats = FALSE,
         limitsize = FALSE)
  ggsave(png_file, p, width = width, height = height, dpi = 450, bg = "white",
         limitsize = FALSE)
  ggsave(tiff_file, p, width = width, height = height, dpi = 600, bg = "white",
         compression = "lzw", limitsize = FALSE)
}

read_one_gsea <- function(cancer) {
  f <- file.path(RAW_DIR, cancer, paste0(cancer, "_Hallmark_GSEA.csv"))
  if (!file.exists(f)) return(NULL)
  x <- data.table::fread(f, data.table = FALSE, check.names = FALSE)
  if (!all(c("Cancer", "Feature") %in% colnames(x))) return(NULL)
  x
}

gsea_all <- dplyr::bind_rows(lapply(cancer_order, read_one_gsea))
if (nrow(gsea_all) == 0) {
  stop("No per-cancer Hallmark GSEA csv files were found in: ", RAW_DIR)
}

p_file <- file.path(RAW_DIR, "C_Cancer_hallmark_GSEA_pvalues.txt")
if (!file.exists(p_file)) {
  p_file <- file.path(RAW_DIR, "C_Cancer_hallmark_all_pvalues.txt")
}
if (!file.exists(p_file)) {
  stop("Cannot find pan-cancer hallmark pvalue table in: ", RAW_DIR)
}

p_all <- data.table::fread(p_file, data.table = FALSE, check.names = FALSE)
required_p <- c("Cancer", "Feature", "P", "FDR")
if (!all(required_p %in% colnames(p_all))) {
  stop("Pvalue table must contain columns: ", paste(required_p, collapse = ", "))
}

merged <- gsea_all %>%
  mutate(
    Cancer = as.character(Cancer),
    Feature = as.character(Feature)
  ) %>%
  select(-any_of(c("P", "FDR", "Significance", "P_label", "FDR_label"))) %>%
  left_join(
    p_all %>%
      mutate(Cancer = as.character(Cancer), Feature = as.character(Feature)) %>%
      select(Cancer, Feature, P, FDR),
    by = c("Cancer", "Feature")
  ) %>%
  mutate(
    NES = if ("NES" %in% colnames(.)) {
      suppressWarnings(as.numeric(.data[["NES"]]))
    } else {
      suppressWarnings(as.numeric(.data[["Effect"]]))
    },
    P = as.numeric(P),
    FDR = as.numeric(FDR),
    Feature_Display = nice_feature(Feature),
    Significance = sig_from_fdr(FDR, ns_blank = TRUE),
    FDR_significance_rule = "Asterisks from FDR only: * <0.05, ** <0.01, *** <0.001, **** <0.0001"
  ) %>%
  filter(Cancer %in% cancer_order, is.finite(NES), is.finite(FDR))

feature_order <- merged %>%
  group_by(Feature, Feature_Display) %>%
  summarise(
    Median_NES = median(NES, na.rm = TRUE),
    Significant_Cancers = sum(FDR < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Significant_Cancers), desc(Median_NES), Feature_Display) %>%
  pull(Feature_Display)

plot_df <- merged %>%
  mutate(
    Cancer = factor(Cancer, levels = cancer_order),
    Feature_Display = factor(Feature_Display, levels = rev(feature_order))
  ) %>%
  arrange(Feature_Display, Cancer)

data.table::fwrite(
  plot_df,
  file.path(DATA_OUT_DIR, "FigureS1D_full_cancer_hallmark_landscape_plot_data.csv")
)

summary_df <- data.frame(
  Item = c(
    "Data_source",
    "N_cancers",
    "N_hallmarks_plotted",
    "N_tiles_plotted",
    "NES_min_plotted",
    "NES_max_plotted",
    "N_FDR_lt_0.05_plotted"
  ),
  Value = c(
    RAW_DIR,
    length(unique(plot_df$Cancer)),
    length(unique(plot_df$Feature)),
    nrow(plot_df),
    signif(min(plot_df$NES, na.rm = TRUE), 5),
    signif(max(plot_df$NES, na.rm = TRUE), 5),
    sum(plot_df$FDR < 0.05, na.rm = TRUE)
  )
)
data.table::fwrite(
  summary_df,
  file.path(DATA_OUT_DIR, "FigureS1D_full_cancer_hallmark_landscape_run_summary.csv")
)

legend_mode <- if (all(plot_df$NES >= 0, na.rm = TRUE)) {
  "positive_only"
} else if (all(plot_df$NES <= 0, na.rm = TRUE)) {
  "negative_only"
} else {
  "signed"
}

base_plot <- ggplot(plot_df, aes(Cancer, Feature_Display, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.20, width = 0.96, height = 0.92) +
  geom_text(aes(label = Significance), size = 1.25, fontface = "bold",
            color = pal_dark) +
  labs(
    title = "Figure S1D  Full cancer hallmark programs",
    subtitle = "All 50 MSigDB Hallmark GSEA results in MO-DDRscore-high versus -low tumors",
    x = "Cancer type",
    y = NULL
  ) +
  theme_repair(7.8) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 7.2,
                               face = "bold"),
    axis.text.y = element_text(size = 6.3, face = "bold"),
    legend.position = "right",
    legend.key.height = unit(30, "pt"),
    legend.key.width = unit(7, "pt")
  )

if (legend_mode == "positive_only") {
  p <- base_plot +
    scale_fill_gradient(
      low = "#FFF3EA",
      high = pal_high2,
      name = "NES",
      breaks = scales::pretty_breaks(n = 4)
    )
} else if (legend_mode == "negative_only") {
  p <- base_plot +
    scale_fill_gradient(
      low = "#E8FAF8",
      high = "#08777A",
      name = "NES",
      breaks = scales::pretty_breaks(n = 4)
    )
} else {
  lim <- max(abs(range(plot_df$NES, na.rm = TRUE)))
  p <- base_plot +
    scale_fill_gradient2(
      low = pal_low,
      mid = "white",
      high = pal_high2,
      midpoint = 0,
      limits = c(-lim, lim),
      name = "NES",
      breaks = scales::pretty_breaks(n = 5)
    )
}

save_plot(p, "FigureS1D_full_cancer_hallmark_landscape", 12.6, 13.2)

cat("\nDone.\n")
cat("Legend mode:", legend_mode, "\n")
cat("Output:\n", file.path(SUP_DIR, "FigureS1D_full_cancer_hallmark_landscape.pdf"), "\n")
print(summary_df)
