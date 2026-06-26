############################################################
# Basic immune analysis for MO-DDRscore high/low groups
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260513)

############################
# 0. Packages
############################

options(repos = c(
  IOBR = "https://iobr.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
))

pkgs <- c(
  "data.table", "dplyr", "tidyr", "ggplot2",
  "tibble", "stringr", "IOBR"
)

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(stringr)
  library(IOBR)
})

############################
# 1. Paths
############################

PROJECT_DIR <- "D:/R_workspace/评分/AD_DDR_project"
PROC_DIR <- file.path(PROJECT_DIR, "01_processed")

OUT_DIR <- file.path(
  PROJECT_DIR,
  "02_Figure2_MultiOmics_Immune",
  "Immune_basic_MO_DDRscore"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

EXPR_FILE  <- file.path(PROC_DIR, "LUAD_tcga_expr_tpm_matrix.rds")
SCORE_FILE <- file.path(PROC_DIR, "LUAD_MO_DDRscore.csv")

############################
# 2. Helper functions
############################

clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub(" ///.*$| //.*$|;.*$|,.*$|\\s*\\([^\\)]*\\)$", "", x)
  x[x %in% c("", "NA", "---", "NULL", "N/A")] <- NA
  x
}

patient_id <- function(x) substr(gsub("\\.", "-", x), 1, 12)
sample_type <- function(x) substr(gsub("\\.", "-", x), 14, 15)
is_tumor <- function(x) sample_type(x) %in% c("01", "02", "03", "05", "06", "07")

safe_wilcox <- function(x, g) {
  keep <- is.finite(x) & !is.na(g)
  x <- x[keep]
  g <- droplevels(factor(g[keep], levels = c("Low", "High")))
  if (length(unique(g)) < 2) return(NA_real_)
  if (length(unique(x)) < 2) return(NA_real_)
  tryCatch(wilcox.test(x ~ g)$p.value, error = function(e) NA_real_)
}

group_test <- function(df, value_col = "Score") {
  x <- df[[value_col]]
  g <- df$MO_DDRscore_group
  
  low <- x[g == "Low"]
  high <- x[g == "High"]
  
  tibble(
    N_Low = sum(is.finite(low)),
    N_High = sum(is.finite(high)),
    Median_Low = median(low, na.rm = TRUE),
    Median_High = median(high, na.rm = TRUE),
    Delta_High_minus_Low = Median_High - Median_Low,
    P = safe_wilcox(x, g)
  )
}

standardize_iobr_result <- function(res, method_name) {
  res <- as.data.frame(res, check.names = FALSE)
  id_col <- intersect(c("ID", "Sample", "sample", "samples"), colnames(res))[1]
  if (is.na(id_col)) {
    stop("Cannot identify sample ID column for method: ", method_name)
  }
  colnames(res)[colnames(res) == id_col] <- "Sample"
  
  res %>%
    tidyr::pivot_longer(
      cols = -Sample,
      names_to = "Feature",
      values_to = "Score"
    ) %>%
    mutate(
      Method = method_name,
      Score = suppressWarnings(as.numeric(Score))
    ) %>%
    filter(is.finite(Score))
}

plot_box <- function(plot_df, title, file, ncol = 4) {
  if (nrow(plot_df) == 0) return(NULL)
  
  plot_df <- plot_df %>%
    mutate(Feature_label = paste(Method, Feature, sep = ": "))
  
  p <- ggplot(plot_df, aes(x = MO_DDRscore_group, y = Score, fill = MO_DDRscore_group)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 0.45, alpha = 0.35) +
    scale_fill_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
    facet_wrap(~ Feature_label, scales = "free_y", ncol = ncol) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 8),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(x = NULL, y = "Score", title = title)
  
  ggsave(file, p, width = 10, height = 7, useDingbats = FALSE)
}

############################
# 3. Load expression and score
############################

expr <- readRDS(EXPR_FILE)
expr <- as.matrix(expr)
storage.mode(expr) <- "numeric"

rownames(expr) <- clean_gene(rownames(expr))
colnames(expr) <- gsub("\\.", "-", colnames(expr))

keep_gene <- !is.na(rownames(expr)) & rownames(expr) != ""
expr <- expr[keep_gene, , drop = FALSE]

# Average duplicated gene symbols
gene_count <- table(rownames(expr))
expr <- rowsum(expr, group = rownames(expr), reorder = FALSE)
expr <- sweep(expr, 1, as.numeric(gene_count[rownames(expr)]), "/")

expr[!is.finite(expr)] <- 0
expr[expr < 0] <- 0

score <- data.table::fread(SCORE_FILE, data.table = FALSE, check.names = FALSE)
score$Sample <- gsub("\\.", "-", score$Sample)

if (!"Patient" %in% colnames(score)) {
  score$Patient <- patient_id(score$Sample)
}

if (!"MO_DDRscore_raw" %in% colnames(score)) {
  if ("MO_DDRscore" %in% colnames(score)) {
    score$MO_DDRscore_raw <- score$MO_DDRscore
  } else {
    stop("Cannot find MO_DDRscore_raw or MO_DDRscore in score file.")
  }
}

score$MO_DDRscore_raw <- as.numeric(score$MO_DDRscore_raw)

if (!"MO_DDRscore_group" %in% colnames(score)) {
  cutoff <- median(score$MO_DDRscore_raw, na.rm = TRUE)
  score$MO_DDRscore_group <- ifelse(score$MO_DDRscore_raw >= cutoff, "High", "Low")
}

score$MO_DDRscore_group <- factor(score$MO_DDRscore_group, levels = c("Low", "High"))

score_df <- score %>%
  filter(
    Sample %in% colnames(expr),
    is_tumor(Sample),
    MO_DDRscore_group %in% c("Low", "High"),
    is.finite(MO_DDRscore_raw)
  ) %>%
  arrange(Patient, Sample) %>%
  distinct(Patient, .keep_all = TRUE) %>%
  select(Sample, Patient, MO_DDRscore_raw, MO_DDRscore_group)

expr_tpm <- expr[, score_df$Sample, drop = FALSE]

cat("Samples used:", ncol(expr_tpm), "\n")
print(table(score_df$MO_DDRscore_group))

data.table::fwrite(score_df, file.path(OUT_DIR, "Immune_sample_annotation.csv"))

############################
# 4. IOBR deconvolution
############################
# Input should be TPM-like expression, gene x sample, not log2 transformed.

run_iobr <- function(method_name, fun) {
  cat("Running:", method_name, "\n")
  out_file <- file.path(OUT_DIR, paste0("IOBR_", method_name, ".csv"))
  err_file <- file.path(OUT_DIR, paste0("IOBR_", method_name, "_ERROR.txt"))
  
  res <- tryCatch(
    {
      x <- fun()
      x <- as.data.frame(x, check.names = FALSE)
      data.table::fwrite(x, out_file)
      x
    },
    error = function(e) {
      writeLines(conditionMessage(e), err_file)
      NULL
    }
  )
  
  res
}

res_list <- list()

res_list$estimate <- run_iobr(
  "estimate",
  function() IOBR::deconvo_tme(eset = expr_tpm, method = "estimate")
)

res_list$mcpcounter <- run_iobr(
  "mcpcounter",
  function() IOBR::deconvo_tme(eset = expr_tpm, method = "mcpcounter")
)

res_list$xcell <- run_iobr(
  "xcell",
  function() IOBR::deconvo_tme(eset = expr_tpm, method = "xcell", arrays = FALSE)
)

res_list$epic <- run_iobr(
  "epic",
  function() IOBR::deconvo_tme(eset = expr_tpm, method = "epic", arrays = FALSE)
)

res_list$quantiseq <- run_iobr(
  "quantiseq",
  function() IOBR::deconvo_tme(
    eset = expr_tpm,
    method = "quantiseq",
    tumor = TRUE,
    arrays = FALSE,
    scale_mrna = TRUE
  )
)

res_list$cibersort <- run_iobr(
  "cibersort",
  function() IOBR::deconvo_tme(
    eset = expr_tpm,
    method = "cibersort",
    arrays = FALSE,
    perm = 100
  )
)

res_list$cibersort_abs <- run_iobr(
  "cibersort_abs",
  function() IOBR::deconvo_tme(
    eset = expr_tpm,
    method = "cibersort_abs",
    arrays = FALSE,
    perm = 100
  )
)

res_list$timer <- run_iobr(
  "timer",
  function() IOBR::deconvo_tme(
    eset = expr_tpm,
    method = "timer",
    group_list = rep("luad", ncol(expr_tpm))
  )
)

res_list$ips <- run_iobr(
  "ips",
  function() IOBR::deconvo_tme(
    eset = expr_tpm,
    method = "ips",
    plot = FALSE
  )
)

############################
# 5. Combine immune deconvolution results
############################

immune_long <- dplyr::bind_rows(lapply(names(res_list), function(nm) {
  if (is.null(res_list[[nm]])) return(NULL)
  standardize_iobr_result(res_list[[nm]], nm)
}))

immune_long <- immune_long %>%
  left_join(score_df, by = "Sample")

data.table::fwrite(
  immune_long,
  file.path(OUT_DIR, "Immune_deconvolution_long.csv")
)

immune_stat <- immune_long %>%
  group_by(Method, Feature) %>%
  group_modify(~ group_test(.x, "Score")) %>%
  ungroup() %>%
  group_by(Method) %>%
  mutate(FDR = p.adjust(P, method = "BH")) %>%
  ungroup() %>%
  arrange(Method, FDR, P)

data.table::fwrite(
  immune_stat,
  file.path(OUT_DIR, "Immune_deconvolution_group_comparison.csv")
)

############################
# 6. Basic plots
############################

# ESTIMATE plot
estimate_plot <- immune_long %>%
  filter(Method == "estimate")

plot_box(
  estimate_plot,
  "ESTIMATE features by MO-DDRscore group",
  file.path(OUT_DIR, "Fig_ESTIMATE_by_MO_DDRscore_group.pdf"),
  ncol = 2
)

# Top differential features from each method
top_features <- immune_stat %>%
  filter(!is.na(FDR)) %>%
  group_by(Method) %>%
  arrange(FDR, P) %>%
  slice_head(n = 6) %>%
  ungroup() %>%
  select(Method, Feature)

top_plot <- immune_long %>%
  inner_join(top_features, by = c("Method", "Feature"))

plot_box(
  top_plot,
  "Top immune deconvolution features by MO-DDRscore group",
  file.path(OUT_DIR, "Fig_Top_immune_features_by_MO_DDRscore_group.pdf"),
  ncol = 3
)

############################
# 7. Standard immune marker expression
############################
# These are pre-specified common immune marker panels.
# Here we compare single-gene expression, not self-defined immune scores.

marker_panels <- list(
  Checkpoint = c(
    "PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3",
    "TIGIT", "HAVCR2", "IDO1", "ICOS", "TNFRSF9"
  ),
  HLA_APM = c(
    "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1",
    "TAP2", "TAPBP", "NLRC5", "PSMB8", "PSMB9"
  ),
  Chemokine_Receptor = c(
    "CXCL9", "CXCL10", "CXCL11", "CCL2", "CCL5",
    "CCR5", "CXCR3", "CXCR4", "CXCL12"
  ),
  Cytolytic = c("GZMA", "PRF1")
)

marker_df <- dplyr::bind_rows(lapply(names(marker_panels), function(panel) {
  genes <- intersect(marker_panels[[panel]], rownames(expr_tpm))
  if (length(genes) == 0) return(NULL)
  
  mat <- log2(expr_tpm[genes, , drop = FALSE] + 1)
  
  as.data.frame(mat, check.names = FALSE) %>%
    rownames_to_column("Gene") %>%
    pivot_longer(
      cols = -Gene,
      names_to = "Sample",
      values_to = "Expression"
    ) %>%
    mutate(Panel = panel)
}))

marker_df <- marker_df %>%
  left_join(score_df, by = "Sample")

data.table::fwrite(
  marker_df,
  file.path(OUT_DIR, "Immune_marker_expression_long.csv")
)

marker_stat <- marker_df %>%
  group_by(Panel, Gene) %>%
  group_modify(~ group_test(.x, "Expression")) %>%
  ungroup() %>%
  group_by(Panel) %>%
  mutate(FDR = p.adjust(P, method = "BH")) %>%
  ungroup() %>%
  arrange(Panel, FDR, P)

data.table::fwrite(
  marker_stat,
  file.path(OUT_DIR, "Immune_marker_group_comparison.csv")
)

plot_marker_panel <- function(panel_name, width = 9, height = 5) {
  df <- marker_df %>% filter(Panel == panel_name)
  if (nrow(df) == 0) return(NULL)
  
  p <- ggplot(df, aes(x = MO_DDRscore_group, y = Expression, fill = MO_DDRscore_group)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 0.45, alpha = 0.35) +
    scale_fill_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 5) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(x = NULL, y = "log2(TPM + 1)", title = panel_name)
  
  ggsave(
    file.path(OUT_DIR, paste0("Fig_", panel_name, "_by_MO_DDRscore_group.pdf")),
    p,
    width = width,
    height = height,
    useDingbats = FALSE
  )
}

plot_marker_panel("Checkpoint")
plot_marker_panel("HLA_APM")
plot_marker_panel("Chemokine_Receptor")
plot_marker_panel("Cytolytic", width = 5, height = 4)

############################
# 8. Optional: cytolytic activity score
############################
# Literature-used CYT score: average expression of GZMA and PRF1.

if (all(c("GZMA", "PRF1") %in% rownames(expr_tpm))) {
  cyt_score <- colMeans(log2(expr_tpm[c("GZMA", "PRF1"), , drop = FALSE] + 1))
  
  cyt_df <- data.frame(
    Sample = names(cyt_score),
    CYT_score = as.numeric(cyt_score)
  ) %>%
    left_join(score_df, by = "Sample")
  
  data.table::fwrite(
    cyt_df,
    file.path(OUT_DIR, "Cytolytic_activity_GZMA_PRF1.csv")
  )
  
  cyt_p <- safe_wilcox(cyt_df$CYT_score, cyt_df$MO_DDRscore_group)
  
  p_cyt <- ggplot(cyt_df, aes(x = MO_DDRscore_group, y = CYT_score, fill = MO_DDRscore_group)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 0.6, alpha = 0.45) +
    scale_fill_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      x = NULL,
      y = "CYT score",
      title = paste0("Cytolytic activity, Wilcoxon P = ", signif(cyt_p, 3))
    )
  
  ggsave(
    file.path(OUT_DIR, "Fig_CYT_score_by_MO_DDRscore_group.pdf"),
    p_cyt,
    width = 4.5,
    height = 4,
    useDingbats = FALSE
  )
}

cat("\nDone.\n")
cat("Output directory:\n", OUT_DIR, "\n")

















library(data.table)

out_dir <- "D:/R_workspace/评分/AD_DDR_project/02_Figure2_MultiOmics_Immune/Immune_official_followup_MO_DDRscore"

x <- fread(file.path(out_dir, "TIDE_upload_gene_by_sample_TPM.csv"), data.table = FALSE)

fwrite(
  x,
  file.path(out_dir, "TIDE_upload_gene_by_sample_TPM.tsv"),
  sep = "\t"
)

