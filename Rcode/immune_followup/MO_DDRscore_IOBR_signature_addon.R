############################################################
# IOBR published-signature add-on analysis
# MO-DDRscore high/low groups
#
# Purpose:
# Add reference-paper-like immune suppression / immune exclusion /
# biomarker / immune response signature analyses using official IOBR
# built-in signature collections.
#
# This script does NOT define custom gene sets.
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
  "tibble", "stringr", "IOBR", "pheatmap"
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
  library(pheatmap)
})

############################
# 1. Paths
############################

PROJECT_DIR <- "D:/R_workspace/评分/AD_DDR_project"
PROC_DIR <- file.path(PROJECT_DIR, "01_processed")
FIG2_DIR <- file.path(PROJECT_DIR, "02_Figure2_MultiOmics_Immune")
BASIC_DIR <- file.path(FIG2_DIR, "Immune_basic_MO_DDRscore")
OUT_DIR <- file.path(FIG2_DIR, "Immune_IOBR_signature_MO_DDRscore")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

EXPR_FILE <- file.path(PROC_DIR, "LUAD_tcga_expr_tpm_matrix.rds")
ANNOT_FILE <- file.path(BASIC_DIR, "Immune_sample_annotation.csv")

############################
# 2. Helper functions
############################

save_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file)
}

clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub(" ///.*$| //.*$|;.*$|,.*$|\\s*\\([^\\)]*\\)$", "", x)
  x[x %in% c("", "NA", "---", "NULL", "N/A")] <- NA
  x
}

patient_id <- function(x) substr(gsub("\\.", "-", x), 1, 12)
safe_num <- function(x) suppressWarnings(as.numeric(x))

safe_wilcox <- function(x, g) {
  keep <- is.finite(x) & !is.na(g)
  x <- x[keep]
  g <- droplevels(factor(g[keep], levels = c("Low", "High")))
  if (length(unique(g)) < 2 || length(unique(x)) < 2) return(NA_real_)
  tryCatch(wilcox.test(x ~ g)$p.value, error = function(e) NA_real_)
}

group_compare <- function(df, value_col = "Score") {
  x <- safe_num(df[[value_col]])
  g <- factor(df$MO_DDRscore_group, levels = c("Low", "High"))
  low <- x[g == "Low"]
  high <- x[g == "High"]
  data.frame(
    N_Low = sum(is.finite(low)),
    N_High = sum(is.finite(high)),
    Median_Low = median(low, na.rm = TRUE),
    Median_High = median(high, na.rm = TRUE),
    Delta_High_minus_Low = median(high, na.rm = TRUE) - median(low, na.rm = TRUE),
    P = safe_wilcox(x, g)
  )
}

safe_spearman <- function(x, y) {
  x <- safe_num(x)
  y <- safe_num(y)
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 10 || length(unique(x)) < 3 || length(unique(y)) < 3) {
    return(c(Rho = NA_real_, P = NA_real_))
  }
  out <- suppressWarnings(cor.test(x, y, method = "spearman"))
  c(Rho = unname(out$estimate), P = out$p.value)
}

zscore <- function(x) {
  x <- safe_num(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

load_iobr_data <- function(obj_name) {
  # IOBR stores multiple built-in data objects. Different versions expose them
  # either from namespace or via data().
  out <- tryCatch(get(obj_name, envir = asNamespace("IOBR")), error = function(e) NULL)
  if (!is.null(out)) return(out)

  env <- new.env(parent = emptyenv())
  ok <- tryCatch({
    data(list = obj_name, package = "IOBR", envir = env)
    TRUE
  }, error = function(e) FALSE)
  if (ok && exists(obj_name, envir = env)) return(get(obj_name, envir = env))
  NULL
}

############################
# 3. Load expression and annotation
############################

anno <- data.table::fread(ANNOT_FILE, data.table = FALSE, check.names = FALSE)
anno$Sample <- gsub("\\.", "-", anno$Sample)
if (!"Patient" %in% colnames(anno)) {
  anno$Patient <- patient_id(anno$Sample)
}
anno <- anno %>%
  mutate(
    Patient = as.character(Patient),
    MO_DDRscore_raw = safe_num(MO_DDRscore_raw),
    MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High"))
  )

expr <- readRDS(EXPR_FILE)
expr <- as.matrix(expr)
storage.mode(expr) <- "numeric"
rownames(expr) <- clean_gene(rownames(expr))
colnames(expr) <- gsub("\\.", "-", colnames(expr))
expr <- expr[!is.na(rownames(expr)) & rownames(expr) != "", , drop = FALSE]

# Average duplicated gene symbols.
gene_count <- table(rownames(expr))
expr <- rowsum(expr, group = rownames(expr), reorder = FALSE)
expr <- sweep(expr, 1, as.numeric(gene_count[rownames(expr)]), "/")
expr[!is.finite(expr)] <- 0
expr[expr < 0] <- 0

samples <- intersect(anno$Sample, colnames(expr))
anno <- anno %>%
  filter(Sample %in% samples) %>%
  arrange(match(Sample, samples))
expr_tpm <- expr[, anno$Sample, drop = FALSE]

cat("Samples:", ncol(expr_tpm), "\n")
print(table(anno$MO_DDRscore_group))

############################
# 4. Load IOBR built-in signatures
############################

signature_collection <- load_iobr_data("signature_collection")
sig_group <- load_iobr_data("sig_group")

if (is.null(signature_collection)) {
  stop("Cannot load IOBR::signature_collection. Please update/reinstall IOBR.")
}

if (!is.list(signature_collection)) {
  stop("IOBR signature_collection is not a list in this version.")
}

all_sig_names <- names(signature_collection)

sig_info <- data.frame(
  Signature = all_sig_names,
  IOBR_group = NA_character_,
  stringsAsFactors = FALSE
)

if (!is.null(sig_group) && is.list(sig_group)) {
  group_map <- bind_rows(lapply(names(sig_group), function(g) {
    data.frame(IOBR_group = g, Signature = as.character(sig_group[[g]]),
               stringsAsFactors = FALSE)
  })) %>%
    filter(Signature %in% all_sig_names) %>%
    group_by(Signature) %>%
    summarise(IOBR_group = paste(unique(IOBR_group), collapse = ";"), .groups = "drop")

  sig_info <- sig_info %>%
    select(Signature) %>%
    left_join(group_map, by = "Signature")
}

save_csv(sig_info, file.path(OUT_DIR, "IOBR_available_signature_info.csv"))

############################
# 5. Select official signatures relevant to the reference paper
############################

# Selection rule:
# 1) Prefer IOBR groups whose group names match immune suppression,
#    exclusion, biomarkers, immune response, TME, or ICB.
# 2) Also include signatures whose published signature names match the same
#    immunotherapy-related keywords.
# This does not create new gene sets; it only filters IOBR's built-in list.

group_keywords <- paste(
  c(
    "suppress", "suppression",
    "exclusion", "exclude",
    "rejection", "reject",
    "response", "responder",
    "biomarker", "checkpoint",
    "immunotherapy", "ICB", "ICI",
    "TME", "CAF", "TGF", "MDSC", "TAM",
    "dysfunction", "exhaust"
  ),
  collapse = "|"
)

name_keywords <- paste(
  c(
    "suppress", "exclusion", "rejection", "response",
    "biomarker", "checkpoint", "PD1", "PD_1", "CTLA4",
    "ICB", "ICI", "TME", "CAF", "TGF", "MDSC", "TAM",
    "dysfunction", "exhaust", "inflamed", "APM",
    "IFNG", "IFN", "cytolytic"
  ),
  collapse = "|"
)

selected_from_group <- sig_info %>%
  filter(!is.na(IOBR_group), grepl(group_keywords, IOBR_group, ignore.case = TRUE)) %>%
  pull(Signature)

selected_from_name <- sig_info %>%
  filter(grepl(name_keywords, Signature, ignore.case = TRUE)) %>%
  pull(Signature)

selected_signatures <- unique(c(selected_from_group, selected_from_name))
selected_signatures <- intersect(selected_signatures, all_sig_names)

if (length(selected_signatures) < 5) {
  warning("Few signatures were selected by keyword. Falling back to all IOBR signatures.")
  selected_signatures <- all_sig_names
}

selected_collection <- signature_collection[selected_signatures]

selection_table <- sig_info %>%
  filter(Signature %in% selected_signatures) %>%
  mutate(
    N_genes_total = lengths(signature_collection[Signature]),
    N_genes_matched = vapply(signature_collection[Signature], function(g) {
      sum(clean_gene(g) %in% rownames(expr_tpm))
    }, numeric(1))
  ) %>%
  arrange(IOBR_group, Signature)

save_csv(selection_table, file.path(OUT_DIR, "IOBR_selected_official_signatures.csv"))

cat("Selected IOBR official signatures:", length(selected_collection), "\n")

############################
# 6. Calculate IOBR signature scores
############################

# IOBR official scoring. PCA is commonly used in IOBR examples for published
# signatures. If PCA fails for a signature because of too few matched genes,
# the script falls back to z-score scoring for the full selected collection.

score_raw <- tryCatch(
  {
    IOBR::calculate_sig_score(
      eset = expr_tpm,
      signature = selected_collection,
      method = "pca",
      mini_gene_count = 3
    )
  },
  error = function(e) {
    message("PCA scoring failed: ", conditionMessage(e))
    message("Falling back to z-score scoring.")
    IOBR::calculate_sig_score(
      eset = expr_tpm,
      signature = selected_collection,
      method = "zscore",
      mini_gene_count = 3
    )
  }
)

score_df <- as.data.frame(score_raw, check.names = FALSE)
id_col <- intersect(c("ID", "Sample", "sample", "samples"), colnames(score_df))[1]
if (is.na(id_col)) {
  stop("Cannot identify sample ID column in IOBR signature score output.")
}
colnames(score_df)[colnames(score_df) == id_col] <- "Sample"
score_df$Sample <- gsub("\\.", "-", score_df$Sample)

save_csv(score_df, file.path(OUT_DIR, "IOBR_official_signature_scores_wide.csv"))

score_long <- score_df %>%
  pivot_longer(
    cols = -Sample,
    names_to = "Signature",
    values_to = "Score"
  ) %>%
  mutate(Score = safe_num(Score)) %>%
  filter(is.finite(Score)) %>%
  left_join(selection_table %>% select(Signature, IOBR_group, N_genes_total, N_genes_matched),
            by = "Signature") %>%
  left_join(anno %>% select(Sample, Patient, MO_DDRscore_raw, MO_DDRscore_group),
            by = "Sample")

save_csv(score_long, file.path(OUT_DIR, "IOBR_official_signature_scores_long.csv"))

############################
# 7. Group comparison and continuous correlation
############################

sig_stat <- score_long %>%
  group_by(Signature, IOBR_group) %>%
  group_modify(~ group_compare(.x, "Score")) %>%
  ungroup() %>%
  mutate(
    FDR = p.adjust(P, method = "BH"),
    Direction = ifelse(Delta_High_minus_Low > 0, "High_up", "High_down")
  ) %>%
  arrange(FDR, P)

sig_cor <- score_long %>%
  group_by(Signature, IOBR_group) %>%
  summarise(
    N = n(),
    Rho = safe_spearman(Score, MO_DDRscore_raw)["Rho"],
    P = safe_spearman(Score, MO_DDRscore_raw)["P"],
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(P, method = "BH"),
    Direction = ifelse(Rho > 0, "Positive", "Negative")
  ) %>%
  arrange(FDR, P)

save_csv(sig_stat, file.path(OUT_DIR, "IOBR_official_signature_group_comparison.csv"))
save_csv(sig_cor, file.path(OUT_DIR, "IOBR_official_signature_vs_MO_DDRscore_spearman.csv"))

############################
# 8. Reference-style plots
############################

top_n <- 40
top_sig <- sig_stat %>%
  filter(is.finite(FDR)) %>%
  arrange(FDR, P) %>%
  slice_head(n = top_n) %>%
  mutate(
    Label = ifelse(is.na(IOBR_group) | IOBR_group == "",
                   Signature,
                   paste0(Signature, " [", IOBR_group, "]")),
    Label = factor(Label, levels = rev(Label)),
    LogFDR = -log10(pmax(FDR, 1e-300))
  )

if (nrow(top_sig) > 0) {
  p <- ggplot(top_sig, aes(x = Delta_High_minus_Low, y = Label)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
    geom_point(aes(size = LogFDR, color = Direction), alpha = 0.9) +
    scale_color_manual(values = c(High_down = "#3B75AF", High_up = "#C84630")) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      x = "Median difference: High - Low",
      y = NULL,
      color = "Direction",
      size = "-log10(FDR)",
      title = "IOBR published immune signatures"
    )
  ggsave(
    file.path(OUT_DIR, "Fig_IOBR_signature_delta_dotplot.pdf"),
    p,
    width = 9,
    height = max(6, 0.22 * nrow(top_sig)),
    useDingbats = FALSE
  )
}

# Make grouped boxplot panels similar to reference-paper Figure 4B-D.
panel_keywords <- list(
  Immune_Suppression = "suppress|TGF|CAF|MDSC|TAM|exhaust|dysfunction",
  Immune_Exclusion = "exclusion|exclude|rejection|reject|CAF|TGF|EMT",
  Biomarkers_Response = "biomarker|response|responder|checkpoint|PD1|PD_1|CTLA4|ICB|ICI|inflamed|APM|IFNG|cytolytic"
)

for (panel in names(panel_keywords)) {
  pat <- panel_keywords[[panel]]
  panel_sigs <- sig_stat %>%
    filter(
      grepl(pat, Signature, ignore.case = TRUE) |
        (!is.na(IOBR_group) & grepl(pat, IOBR_group, ignore.case = TRUE))
    ) %>%
    arrange(FDR, P) %>%
    slice_head(n = 18) %>%
    pull(Signature)

  panel_df <- score_long %>%
    filter(Signature %in% panel_sigs) %>%
    mutate(Signature = factor(Signature, levels = panel_sigs))

  if (nrow(panel_df) == 0) next

  p <- ggplot(panel_df, aes(x = MO_DDRscore_group, y = Score, fill = MO_DDRscore_group)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.9) +
    geom_jitter(width = 0.15, size = 0.35, alpha = 0.25) +
    scale_fill_manual(values = c(Low = "#E3B23C", High = "#1E5A86")) +
    facet_wrap(~ Signature, scales = "free_y", ncol = 6) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      legend.position = "top",
      axis.text.x = element_text(angle = 0, color = "black"),
      axis.text.y = element_text(color = "black"),
      strip.text = element_text(size = 7, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(x = NULL, y = "Signature score", fill = "MO-DDRscore", title = panel)

  ggsave(
    file.path(OUT_DIR, paste0("Fig_IOBR_", panel, "_boxplot_panel.pdf")),
    p,
    width = 13,
    height = 7,
    useDingbats = FALSE
  )
}

# Heatmap of top changed signatures.
top_heat_sigs <- sig_stat %>%
  filter(is.finite(FDR)) %>%
  arrange(FDR, P) %>%
  slice_head(n = 35) %>%
  pull(Signature)

heat_df <- score_long %>%
  filter(Signature %in% top_heat_sigs) %>%
  group_by(Signature) %>%
  mutate(Score_z = zscore(Score)) %>%
  ungroup()

if (nrow(heat_df) > 0) {
  heat_mat <- heat_df %>%
    select(Signature, Sample, Score_z) %>%
    pivot_wider(names_from = Sample, values_from = Score_z) %>%
    column_to_rownames("Signature") %>%
    as.matrix()

  sample_order <- anno %>%
    arrange(MO_DDRscore_group, MO_DDRscore_raw) %>%
    pull(Sample)
  sample_order <- intersect(sample_order, colnames(heat_mat))
  heat_mat <- heat_mat[, sample_order, drop = FALSE]

  ann_col <- anno %>%
    select(Sample, MO_DDRscore_group, MO_DDRscore_raw) %>%
    filter(Sample %in% sample_order) %>%
    arrange(match(Sample, sample_order)) %>%
    column_to_rownames("Sample")

  pdf(file.path(OUT_DIR, "Fig_IOBR_top_signature_heatmap.pdf"), width = 10, height = 9)
  pheatmap::pheatmap(
    heat_mat,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_colnames = FALSE,
    annotation_col = ann_col,
    color = colorRampPalette(c("#3B75AF", "white", "#C84630"))(100),
    breaks = seq(-2.5, 2.5, length.out = 101),
    main = "Top IOBR published immune signatures"
  )
  dev.off()
}

############################
# 9. Summary
############################

summary_df <- data.frame(
  Item = c(
    "N_samples",
    "N_available_IOBR_signatures",
    "N_selected_IOBR_signatures",
    "N_group_FDR_lt_0.05",
    "N_correlation_FDR_lt_0.05"
  ),
  Value = c(
    nrow(anno),
    length(all_sig_names),
    length(selected_collection),
    sum(sig_stat$FDR < 0.05, na.rm = TRUE),
    sum(sig_cor$FDR < 0.05, na.rm = TRUE)
  )
)

save_csv(summary_df, file.path(OUT_DIR, "IOBR_signature_addon_summary.csv"))

cat("\nDone.\n")
cat("Output directory:\n", OUT_DIR, "\n")
print(summary_df)
