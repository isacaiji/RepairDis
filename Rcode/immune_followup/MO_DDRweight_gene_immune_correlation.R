############################################################
# Top MO-DDRweight genes vs immune features correlation
# Data-table only; plotting will be handled later.
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260513)

############################
# 0. Packages
############################

pkgs <- c("data.table", "dplyr", "tidyr", "tibble", "stringr")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

############################
# 1. Paths and parameters
############################

PROJECT_DIR <- "D:/R_workspace/评分/AD_DDR_project"
PROC_DIR <- file.path(PROJECT_DIR, "01_processed")
FIG2_DIR <- file.path(PROJECT_DIR, "02_Figure2_MultiOmics_Immune")
BASIC_DIR <- file.path(FIG2_DIR, "Immune_basic_MO_DDRscore")
FOLLOW_DIR <- file.path(FIG2_DIR, "Immune_official_followup_MO_DDRscore")
IOBR_SIG_DIR <- file.path(FIG2_DIR, "Immune_IOBR_signature_MO_DDRscore")
OUT_DIR <- file.path(FIG2_DIR, "Gene_Immune_Correlation_MO_DDRweight")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

EXPR_FILE <- file.path(PROC_DIR, "LUAD_tcga_expr_tpm_matrix.rds")
WEIGHT_FILE <- file.path(PROC_DIR, "MO_DDRweight_gene_table.csv")
ANNOT_FILE <- file.path(BASIC_DIR, "Immune_sample_annotation.csv")

TOP_N <- 30

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

find_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

normalize_sample_col <- function(df) {
  id_col <- find_col(df, c(
    "Sample", "sample", "ID", "id", "Tumor_Sample_Barcode",
    "Tumor.Sample.Barcode", "Patient", "patient", "submitter_id",
    "V1", "X", "...1"
  ))
  if (!is.na(id_col)) {
    colnames(df)[colnames(df) == id_col] <- "Sample"
    df$Sample <- gsub("\\.", "-", df$Sample)
    df$Patient <- patient_id(df$Sample)
  }
  df
}

safe_spearman <- function(x, y) {
  x <- safe_num(x)
  y <- safe_num(y)
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 10 || length(unique(x)) < 3 || length(unique(y)) < 3) {
    return(c(N = length(x), Rho = NA_real_, P = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
  c(N = length(x), Rho = unname(ct$estimate), P = ct$p.value)
}

add_feature <- function(feature_df, source_name, feature_name, sample_col = "Sample",
                        value_col = "Value") {
  if (is.null(feature_df) || nrow(feature_df) == 0) return(NULL)
  feature_df %>%
    transmute(
      Sample = .data[[sample_col]],
      Source = source_name,
      Feature = feature_name,
      Value = safe_num(.data[[value_col]])
    ) %>%
    filter(is.finite(Value))
}

############################
# 3. Load annotation, expression, and MO-DDRweight
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
expr_log <- log2(expr_tpm + 1)

weight <- data.table::fread(WEIGHT_FILE, data.table = FALSE, check.names = FALSE)
weight$Gene <- clean_gene(weight$Gene)
weight <- weight %>% filter(!is.na(Gene), Gene %in% rownames(expr_log))

weight_col <- find_col(weight, c("MO_DDRweight", "Raw_MO_DDRweight", "MO_DDR_weight"))
if (is.na(weight_col)) stop("Cannot identify MO-DDRweight column.")

weight <- weight %>%
  mutate(WeightValue = safe_num(.data[[weight_col]])) %>%
  filter(is.finite(WeightValue)) %>%
  arrange(desc(WeightValue))

top_genes <- weight %>%
  slice_head(n = TOP_N) %>%
  pull(Gene)

save_csv(weight %>% slice_head(n = TOP_N),
         file.path(OUT_DIR, "Top_MO_DDRweight_genes_used.csv"))

cat("Samples:", ncol(expr_log), "\n")
cat("Top genes:", length(top_genes), "\n")

############################
# 4. Build immune feature matrix
############################

feature_long_list <- list()

# 4.1 Result-driven axes from official follow-up analysis.
axis_file <- file.path(FOLLOW_DIR, "Official_deconvolution_result_driven_axes_long.csv")
if (file.exists(axis_file)) {
  axis_long <- data.table::fread(axis_file, data.table = FALSE, check.names = FALSE) %>%
    mutate(Sample = gsub("\\.", "-", Sample), Score = safe_num(Score))
  feature_long_list$axes <- axis_long %>%
    transmute(Sample, Source = "ImmuneAxis", Feature = Axis, Value = Score)
}

# 4.2 TIDEpy output merged by previous script.
tide_file <- file.path(FOLLOW_DIR, "Official_TIDE_merged.csv")
if (file.exists(tide_file)) {
  tide <- data.table::fread(tide_file, data.table = FALSE, check.names = FALSE) %>%
    mutate(Sample = gsub("\\.", "-", Sample))
  tide_features <- intersect(
    c("TIDE", "Dysfunction", "Exclusion", "MDSC", "CAF", "TAM M2",
      "IFNG", "CD274", "CD8", "CTL", "MSI Score"),
    colnames(tide)
  )
  feature_long_list$tide <- bind_rows(lapply(tide_features, function(ff) {
    add_feature(tide, "TIDEpy", ff, value_col = ff)
  }))
}

# 4.3 TMB and MATH.
tmb_file <- file.path(FOLLOW_DIR, "Official_TMB_merged.csv")
if (file.exists(tmb_file)) {
  tmb <- data.table::fread(tmb_file, data.table = FALSE, check.names = FALSE) %>%
    normalize_sample_col()
  if ("TMB_value" %in% colnames(tmb)) {
    feature_long_list$tmb <- add_feature(tmb, "MutationBiomarker", "TMB", value_col = "TMB_value")
  }
}

math_file <- file.path(FOLLOW_DIR, "Official_MATH_score.csv")
if (file.exists(math_file)) {
  math <- data.table::fread(math_file, data.table = FALSE, check.names = FALSE)
  if (!"Sample" %in% colnames(math)) {
    math <- math %>% left_join(anno %>% select(Patient, Sample), by = "Patient")
  }
  if ("MATH" %in% colnames(math)) {
    feature_long_list$math <- add_feature(math, "MutationBiomarker", "MATH", value_col = "MATH")
  }
}

# 4.4 Immune marker expression: checkpoint, HLA/APM, CXCL axis, CYT genes.
marker_file <- file.path(BASIC_DIR, "Immune_marker_expression_long.csv")
if (file.exists(marker_file)) {
  marker <- data.table::fread(marker_file, data.table = FALSE, check.names = FALSE) %>%
    mutate(Sample = gsub("\\.", "-", Sample), Expression = safe_num(Expression))
  marker_keep <- c(
    "CD274", "PDCD1LG2", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "TNFRSF9",
    "TAP1", "TAP2", "TAPBP", "NLRC5", "PSMB8", "PSMB9", "HLA-A", "HLA-B", "HLA-C", "B2M",
    "CXCL9", "CXCL10", "CXCL11", "CCL2", "CCL5", "CXCR3", "CXCR4",
    "GZMA", "PRF1"
  )
  feature_long_list$markers <- marker %>%
    filter(Gene %in% marker_keep) %>%
    transmute(Sample, Source = paste0("Marker_", Panel), Feature = Gene, Value = Expression)
}

# 4.5 Selected official MSigDB modules.
msig_file <- file.path(FOLLOW_DIR, "Official_MSigDB_ssGSEA_long.csv")
if (file.exists(msig_file)) {
  msig <- data.table::fread(msig_file, data.table = FALSE, check.names = FALSE) %>%
    mutate(Sample = gsub("\\.", "-", Sample), Score = safe_num(Score))
  msig_keep <- c(
    "GO_B_CELL_ACTIVATION", "GO_BCR_SIGNALING", "GO_PLASMA_CELL_DIFFERENTIATION",
    "GO_ANTIGEN_PRESENTATION", "REACTOME_ANTIGEN_PRESENTATION",
    "HALLMARK_IFNG", "HALLMARK_IFNA", "REACTOME_PD1_SIGNALING",
    "HALLMARK_TGF_BETA", "HALLMARK_EMT", "GO_NEUTROPHIL_ACTIVATION",
    "GO_MYELOID_LEUKOCYTE_ACTIVATION", "HALLMARK_INFLAMMATION"
  )
  feature_long_list$msig <- msig %>%
    filter(Module %in% msig_keep) %>%
    group_by(Sample, Module) %>%
    summarise(Value = mean(Score, na.rm = TRUE), .groups = "drop") %>%
    transmute(Sample, Source = "MSigDB_ssGSEA", Feature = Module, Value = Value)
}

# 4.6 IOBR published signatures, selected high-value immune signatures.
iobr_file <- file.path(IOBR_SIG_DIR, "IOBR_official_signature_scores_long.csv")
if (file.exists(iobr_file)) {
  iobr <- data.table::fread(iobr_file, data.table = FALSE, check.names = FALSE) %>%
    mutate(Sample = gsub("\\.", "-", Sample), Score = safe_num(Score))
  iobr_keep_patterns <- paste(
    c(
      "B_cells", "TLS", "MHC_Class_II", "MDSC", "WNT_target",
      "TMEscoreA", "TMEscoreA_CIR", "TMEscore_plus", "TMEscore_CIR",
      "TIP_Release", "TIP_Infiltration", "Antigen_Processing",
      "TCR_signaling", "Natural_Killer_Cell_Cytotoxicity",
      "CD8_T_cells", "DDR", "Mismatch_Repair", "Homologous_recombination",
      "Cell_cycle", "DNA_replication"
    ),
    collapse = "|"
  )
  feature_long_list$iobr <- iobr %>%
    filter(grepl(iobr_keep_patterns, Signature, ignore.case = TRUE)) %>%
    transmute(Sample, Source = "IOBR_signature", Feature = Signature, Value = Score)
}

feature_long <- bind_rows(feature_long_list) %>%
  filter(Sample %in% anno$Sample, is.finite(Value)) %>%
  distinct(Source, Feature, Sample, .keep_all = TRUE) %>%
  left_join(anno %>% select(Sample, Patient, MO_DDRscore_raw, MO_DDRscore_group), by = "Sample")

save_csv(feature_long, file.path(OUT_DIR, "Immune_features_used_long.csv"))

feature_summary <- feature_long %>%
  distinct(Source, Feature) %>%
  count(Source, name = "N_features") %>%
  arrange(Source)
save_csv(feature_summary, file.path(OUT_DIR, "Immune_features_used_summary.csv"))

############################
# 5. Gene expression long table
############################

gene_expr_long <- as.data.frame(expr_log[top_genes, , drop = FALSE], check.names = FALSE) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "Expression"
  ) %>%
  mutate(Expression = safe_num(Expression)) %>%
  left_join(weight %>% select(Gene, WeightValue), by = "Gene") %>%
  left_join(anno %>% select(Sample, Patient, MO_DDRscore_raw, MO_DDRscore_group), by = "Sample")

save_csv(gene_expr_long, file.path(OUT_DIR, "Top_MO_DDRweight_gene_expression_long.csv"))

############################
# 6. Spearman correlations
############################

cor_df <- gene_expr_long %>%
  select(Gene, Sample, GeneExpression = Expression, WeightValue) %>%
  inner_join(feature_long %>% select(Sample, Source, Feature, ImmuneValue = Value),
             by = "Sample", relationship = "many-to-many") %>%
  group_by(Gene, WeightValue, Source, Feature) %>%
  summarise(
    N = safe_spearman(GeneExpression, ImmuneValue)["N"],
    Rho = safe_spearman(GeneExpression, ImmuneValue)["Rho"],
    P = safe_spearman(GeneExpression, ImmuneValue)["P"],
    .groups = "drop"
  ) %>%
  group_by(Source, Feature) %>%
  mutate(FDR_within_feature = p.adjust(P, method = "BH")) %>%
  ungroup() %>%
  group_by(Gene) %>%
  mutate(FDR_within_gene = p.adjust(P, method = "BH")) %>%
  ungroup() %>%
  mutate(
    Direction = ifelse(Rho > 0, "Positive", "Negative"),
    AbsRho = abs(Rho)
  ) %>%
  arrange(desc(AbsRho), P)

save_csv(cor_df, file.path(OUT_DIR, "Top_MO_DDRweight_gene_immune_feature_spearman.csv"))

############################
# 7. Focused summary tables
############################

focus_features <- c(
  "Global_TME", "B_Plasma_axis", "Myeloid_Neutrophil_axis",
  "quanTIseq_M1_M2_ratio", "IPS_axis",
  "TIDE", "Dysfunction", "Exclusion", "MDSC", "CAF", "TAM M2",
  "TMB", "MATH",
  "CD274", "PDCD1LG2", "LAG3", "TNFRSF9",
  "TAP1", "TAP2", "NLRC5", "PSMB8", "PSMB9",
  "CXCL9", "CXCL10", "CXCL11"
)

focused_cor <- cor_df %>%
  filter(Feature %in% focus_features) %>%
  mutate(FeatureID = paste(Source, Feature, sep = "__")) %>%
  arrange(Feature, desc(AbsRho), P)

save_csv(focused_cor, file.path(OUT_DIR, "Focused_gene_immune_feature_spearman.csv"))

top_per_feature <- focused_cor %>%
  group_by(Source, Feature) %>%
  arrange(desc(AbsRho), P) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  arrange(Source, Feature, desc(AbsRho))

save_csv(top_per_feature, file.path(OUT_DIR, "Top5_genes_per_key_immune_feature.csv"))

top_per_gene <- focused_cor %>%
  group_by(Gene) %>%
  arrange(desc(AbsRho), P) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  arrange(Gene, desc(AbsRho))

save_csv(top_per_gene, file.path(OUT_DIR, "Top10_immune_features_per_MO_DDRweight_gene.csv"))

# Gene-level compact matrix for later heatmap plotting.
rho_matrix <- focused_cor %>%
  select(Gene, FeatureID, Rho) %>%
  group_by(Gene, FeatureID) %>%
  summarise(Rho = mean(Rho, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = FeatureID, values_from = Rho) %>%
  arrange(match(Gene, top_genes))

save_csv(rho_matrix, file.path(OUT_DIR, "Focused_gene_immune_rho_matrix.csv"))

############################
# 8. Score-level correlations for reference
############################

score_cor <- feature_long %>%
  group_by(Source, Feature) %>%
  summarise(
    N = safe_spearman(Value, MO_DDRscore_raw)["N"],
    Rho = safe_spearman(Value, MO_DDRscore_raw)["Rho"],
    P = safe_spearman(Value, MO_DDRscore_raw)["P"],
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(P, method = "BH"),
         Direction = ifelse(Rho > 0, "Positive", "Negative"),
         AbsRho = abs(Rho)) %>%
  arrange(desc(AbsRho), P)

save_csv(score_cor, file.path(OUT_DIR, "MO_DDRscore_immune_feature_spearman.csv"))

############################
# 9. Summary
############################

summary_df <- data.frame(
  Item = c(
    "N_samples",
    "N_top_MO_DDRweight_genes",
    "N_immune_features",
    "N_gene_feature_tests",
    "N_focused_gene_feature_tests",
    "N_gene_feature_FDR_within_feature_lt_0.05",
    "N_score_feature_FDR_lt_0.05"
  ),
  Value = c(
    nrow(anno),
    length(top_genes),
    nrow(feature_long %>% distinct(Source, Feature)),
    nrow(cor_df),
    nrow(focused_cor),
    sum(cor_df$FDR_within_feature < 0.05, na.rm = TRUE),
    sum(score_cor$FDR < 0.05, na.rm = TRUE)
  )
)

save_csv(summary_df, file.path(OUT_DIR, "Gene_immune_correlation_summary.csv"))

cat("\nDone.\n")
cat("Output directory:\n", OUT_DIR, "\n")
print(summary_df)
