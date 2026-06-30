############################################################
# linkET-standard gene-immune cell relationship plot
# Similar to literature-style Fig. 6I:
# qcorrplot + geom_square + geom_mark + geom_couple
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260513)

############################
# 0. Packages
############################

pkgs <- c("data.table", "dplyr", "tidyr", "ggplot2", "linkET", "scales", "stringr")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p == "linkET") {
      remotes::install_github("Hy4m/linkET", upgrade = "never", dependencies = TRUE)
    } else {
      install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(linkET)
  library(scales)
  library(stringr)
})

############################
# 1. Paths and parameters
############################

PROJECT_DIR <- "D:/R_workspace/\u8bc4\u5206/AD_DDR_project"
FIG2_DIR <- file.path(PROJECT_DIR, "02_Figure2_MultiOmics_Immune")
BASIC_DIR <- file.path(FIG2_DIR, "Immune_basic_MO_DDRscore")
GENE_DIR <- file.path(FIG2_DIR, "Gene_Immune_Correlation_MO_DDRweight")
OUT_DIR <- file.path(FIG2_DIR, "Immune_publication_figures_MO_DDRscore")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DECONV_FILE <- file.path(BASIC_DIR, "Immune_deconvolution_long.csv")
GENE_EXPR_FILE <- file.path(GENE_DIR, "Top_MO_DDRweight_gene_expression_long.csv")
WEIGHT_FILE <- file.path(GENE_DIR, "Top_MO_DDRweight_genes_used.csv")

FOCUSED_GENE_SET <- c("UBE2T", "DNA2", "RECQL4", "CHEK1", "XRCC2", "RAD51", "BRIP1", "BRCA1")
EDGE_FDR_CUTOFF <- 0.05
EDGE_ABS_RHO_CUTOFF <- 0.15

############################
# 2. Helpers
############################

read_csv <- function(file) {
  data.table::fread(file, data.table = FALSE, check.names = FALSE)
}

save_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file)
}

num <- function(x) suppressWarnings(as.numeric(x))

safe_spearman <- function(x, y) {
  x <- num(x)
  y <- num(y)
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 10 || length(unique(x)) < 3 || length(unique(y)) < 3) {
    return(c(N = length(x), Rho = NA_real_, P = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
  c(N = length(x), Rho = unname(ct$estimate), P = ct$p.value)
}

save_plot <- function(plot, name, width, height) {
  pdf_file <- file.path(OUT_DIR, paste0(name, ".pdf"))
  png_file <- file.path(OUT_DIR, paste0(name, ".png"))
  tryCatch(
    ggsave(pdf_file, plot, width = width, height = height, device = cairo_pdf, bg = "white"),
    error = function(e) ggsave(pdf_file, plot, width = width, height = height, device = "pdf", bg = "white")
  )
  ggsave(png_file, plot, width = width, height = height, dpi = 500, bg = "white")
}

############################
# 3. Load and shape data
############################

cell_map <- c(
  "B_cells_quantiseq" = "B cells",
  "Dendritic_cells_quantiseq" = "Dendritic cells",
  "Macrophages_M1_quantiseq" = "M1 macrophages",
  "Macrophages_M2_quantiseq" = "M2 macrophages",
  "Monocytes_quantiseq" = "Monocytes",
  "Neutrophils_quantiseq" = "Neutrophils",
  "NK_cells_quantiseq" = "NK cells",
  "T_cells_CD4_quantiseq" = "CD4 T cells",
  "T_cells_CD8_quantiseq" = "CD8 T cells",
  "Tregs_quantiseq" = "Tregs"
)

cell_order <- unname(cell_map)

deconv <- read_csv(DECONV_FILE) %>%
  mutate(
    Sample = gsub("\\.", "-", Sample),
    Method = as.character(Method),
    Feature = as.character(Feature),
    Score = num(Score)
  ) %>%
  filter(Method == "quantiseq", Feature %in% names(cell_map), is.finite(Score)) %>%
  mutate(Cell = unname(cell_map[Feature])) %>%
  distinct(Sample, Cell, .keep_all = TRUE)

cell_wide <- deconv %>%
  select(Sample, Cell, Score) %>%
  pivot_wider(names_from = Cell, values_from = Score) %>%
  select(Sample, all_of(cell_order)) %>%
  filter(if_all(all_of(cell_order), is.finite))

gene_expr <- read_csv(GENE_EXPR_FILE) %>%
  mutate(
    Sample = gsub("\\.", "-", Sample),
    Gene = as.character(Gene),
    Expression = num(Expression)
  )

weight <- read_csv(WEIGHT_FILE) %>%
  mutate(Gene = as.character(Gene), WeightValue = num(WeightValue)) %>%
  arrange(desc(WeightValue))

gene_keep <- intersect(FOCUSED_GENE_SET, weight$Gene)
gene_expr <- gene_expr %>%
  filter(Gene %in% gene_keep, Sample %in% cell_wide$Sample)

cat("Samples:", length(intersect(unique(cell_wide$Sample), unique(gene_expr$Sample))), "\n")
cat("Genes:", paste(gene_keep, collapse = ", "), "\n")

############################
# 4. Gene-cell Spearman table
############################

gene_cell <- gene_expr %>%
  select(Sample, Gene, GeneExpression = Expression) %>%
  inner_join(deconv %>% select(Sample, Cell, CellScore = Score),
             by = "Sample", relationship = "many-to-many") %>%
  group_by(Gene, Cell) %>%
  summarise(
    N = safe_spearman(GeneExpression, CellScore)["N"],
    Rho = safe_spearman(GeneExpression, CellScore)["Rho"],
    P = safe_spearman(GeneExpression, CellScore)["P"],
    .groups = "drop"
  ) %>%
  group_by(Cell) %>%
  mutate(FDR_within_cell = p.adjust(P, method = "BH")) %>%
  ungroup() %>%
  left_join(weight %>% select(Gene, WeightValue), by = "Gene") %>%
  mutate(
    Direction = case_when(
      FDR_within_cell >= EDGE_FDR_CUTOFF ~ "Not significant",
      Rho >= 0 ~ "Positive",
      TRUE ~ "Negative"
    ),
    AbsRho = abs(Rho),
    Strength = cut(
      AbsRho,
      breaks = c(-Inf, 0.20, 0.25, 0.30, Inf),
      labels = c("< 0.20", "0.20-0.25", "0.25-0.30", ">= 0.30")
    )
  ) %>%
  arrange(desc(AbsRho), P)

edge_df <- gene_cell %>%
  filter(FDR_within_cell < EDGE_FDR_CUTOFF, AbsRho >= EDGE_ABS_RHO_CUTOFF) %>%
  mutate(
    Gene = factor(Gene, levels = gene_keep),
    Cell = factor(Cell, levels = cell_order),
    Direction = factor(Direction, levels = c("Positive", "Negative"))
  )

save_csv(gene_cell, file.path(OUT_DIR, "linkET_gene_quantiseq_cell_correlation_all.csv"))
save_csv(edge_df, file.path(OUT_DIR, "linkET_gene_quantiseq_cell_correlation_edges.csv"))

############################
# 5. linkET-standard plot
############################

cell_mat <- cell_wide %>%
  select(all_of(cell_order))

cor_obj <- linkET::correlate(cell_mat, method = "spearman")

ink <- "#17213A"
muted <- "#687386"
pos_col <- "#E88B2D"
neg_col <- "#078C80"

p <- linkET::qcorrplot(cor_obj, type = "lower", diag = TRUE, grid_col = "white", grid_size = 0.35) +
  linkET::geom_square(aes(fill = r), color = "white", linewidth = 0.3) +
  linkET::geom_mark(
    aes(pvalue = p),
    only_mark = TRUE,
    sig_level = c(0.05, 0.01, 0.001),
    mark = c("*", "**", "***"),
    size = 3.0,
    color = ink,
    fontface = "bold"
  ) +
  linkET::geom_couple(
    data = edge_df,
    aes(from = Gene, to = Cell, colour = Direction, size = AbsRho),
    curvature = linkET::nice_curvature(0.12),
    nudge_x = 0.65,
    label.size = 5.2,
    label.colour = ink,
    label.fontface = "bold",
    alpha = 0.86
  ) +
  scale_fill_gradient2(
    low = "#2D6FAB",
    mid = "#F8F9FB",
    high = "#D54D3D",
    midpoint = 0,
    limits = c(-1, 1),
    oob = squish,
    name = "Cell-cell\ncorrelation"
  ) +
  scale_color_manual(
    values = c(Positive = pos_col, Negative = neg_col),
    name = "Gene-cell\nassociation"
  ) +
  scale_size_continuous(
    range = c(0.35, 2.4),
    breaks = c(0.20, 0.25, 0.30, 0.35),
    name = "abs(rho)"
  ) +
  labs(
    title = "MO-DDRweight gene-immune cell association map",
    subtitle = "linkET qcorrplot + geom_couple; arcs indicate significant gene-cell Spearman associations",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 10, base_family = "Arial") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "#FBFCFE", color = NA),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", color = ink, size = 15),
    plot.subtitle = element_text(color = muted, size = 9.4, margin = margin(b = 8)),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, color = "#263242", size = 9),
    axis.text.y = element_text(color = "#263242", size = 9),
    legend.position = "right",
    legend.title = element_text(face = "bold", color = ink, size = 8.8),
    legend.text = element_text(color = "#263242", size = 8.3),
    plot.margin = margin(10, 25, 10, 10)
  )

save_plot(p, "Figure6I_linkET_MO_DDRweight_gene_immune_cell_network", 10.8, 6.8)

cat("\nDone.\n")
cat("Output directory:\n", OUT_DIR, "\n")
cat("linkET figure:\n",
    file.path(OUT_DIR, "Figure6I_linkET_MO_DDRweight_gene_immune_cell_network.pdf"), "\n")
