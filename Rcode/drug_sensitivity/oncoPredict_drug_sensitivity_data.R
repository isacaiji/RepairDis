############################################################
# Drug sensitivity prediction for MO-DDRscore and DPRS
# Data-only script; plotting is handled separately.
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260531)

required_pkgs <- c("data.table", "dplyr", "tidyr", "readxl", "oncoPredict")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(oncoPredict)
})

############################################################
# 1. Paths and parameters
############################################################

# Avoid non-ASCII path literals in Rscript on Windows.
PROJECT_DIR <- file.path(
  "D:/R_workspace",
  intToUtf8(c(0x8bc4, 0x5206)),
  "AD_DDR_project"
)

DATA_DIR <- file.path(PROJECT_DIR, "00_data")
PROC_DIR <- file.path(PROJECT_DIR, "01_processed")
DPRS_DIR <- file.path(PROJECT_DIR, "04-5_Mime1Matched_StepCoxFixed")

DRUG_DIR <- file.path(PROJECT_DIR, "drug")
SCRIPT_DIR <- file.path(DRUG_DIR, "01-script")
DATA_OUT_DIR <- file.path(DRUG_DIR, "02-data")
RES_DIR <- file.path(DRUG_DIR, "03-res")
TABLE_DIR <- file.path(RES_DIR, "tables")
ONCOPREDICT_DIR <- file.path(RES_DIR, "oncoPredict_run")

dir.create(SCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ONCOPREDICT_DIR, recursive = TRUE, showWarnings = FALSE)

CELL_EXPR_FILE <- file.path(DATA_DIR, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
MODEL_FILE <- file.path(DATA_DIR, "Model.csv")
GDSC2_FILE <- file.path(DATA_DIR, "GDSC2_fitted_dose_response_27Oct23.xlsx")
TCGA_EXPR_FILE <- file.path(PROC_DIR, "LUAD_tcga_expr_tpm_matrix.rds")
MO_SCORE_FILE <- file.path(PROC_DIR, "LUAD_MO_DDRscore.csv")
DPRS_FILE <- file.path(DPRS_DIR, "Fig4C_DPRS_all_sets.csv")

MIN_CELL_LINES_PER_DRUG <- 30
RUN_PRIORITY_DRUGS_ONLY <- TRUE

PRIORITY_DRUG_PATTERNS <- c(
  "PARP", "OLAPARIB", "TALAZOPARIB", "NIRAPARIB", "RUCAPARIB",
  "ATR", "AZD6738", "VE-822", "BERZOSERTIB",
  "WEE1", "AZD1775", "ADAVOSERTIB",
  "CHK", "CHEK", "PREXASERTIB",
  "CISPLATIN", "CARBOPLATIN", "OXALIPLATIN", "PACLITAXEL",
  "DOCETAXEL", "GEMCITABINE", "PEMETREXED", "ETOPOSIDE",
  "ERLOTINIB", "GEFITINIB", "AFATINIB", "OSIMERTINIB",
  "MEK", "TRAMETINIB", "SELUMETINIB", "PI3K", "AKT", "MTOR"
)

############################################################
# 2. Helper functions
############################################################

save_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file)
}

patient_id <- function(x) substr(gsub("\\.", "-", as.character(x)), 1, 12)

sample_type <- function(x) substr(gsub("\\.", "-", as.character(x)), 14, 15)

is_tumor <- function(x) sample_type(x) %in% c("01", "02", "03", "05", "06", "07")

clean_gene <- function(x) {
  x <- as.character(x)
  x <- gsub("\\s*\\([^\\)]*\\)\\s*$", "", x)
  x <- gsub(" ///.*$| //.*$|;.*$|,.*$", "", x)
  x <- toupper(trimws(x))
  x[x %in% c("", "NA", "---", "NULL", "N/A")] <- NA_character_
  x
}

safe_name <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[x == "" | is.na(x)] <- "UNKNOWN"
  x
}

collapse_rows_by_mean <- function(mat) {
  rn <- rownames(mat)
  keep <- !is.na(rn) & rn != ""
  mat <- mat[keep, , drop = FALSE]
  rn <- rn[keep]
  rownames(mat) <- rn
  if (!anyDuplicated(rn)) return(mat)

  sums <- rowsum(mat, group = rn, reorder = FALSE, na.rm = TRUE)
  counts <- rowsum(matrix(1, nrow = length(rn), ncol = 1), group = rn, reorder = FALSE)
  sums / as.numeric(counts[, 1])
}

wilcox_compare_drugs <- function(pred_mat, anno, group_col, score_col, label) {
  stopifnot(all(rownames(pred_mat) %in% anno$Sample))
  anno <- anno[match(rownames(pred_mat), anno$Sample), , drop = FALSE]

  res <- lapply(colnames(pred_mat), function(drug) {
    value <- as.numeric(pred_mat[, drug])
    df <- data.frame(
      Sample = rownames(pred_mat),
      Group = anno[[group_col]],
      Score = anno[[score_col]],
      Predicted_IC50 = value,
      stringsAsFactors = FALSE
    )
    df <- df[is.finite(df$Predicted_IC50) & !is.na(df$Group), , drop = FALSE]
    df$Group <- as.character(df$Group)

    groups <- sort(unique(df$Group))
    if (length(groups) != 2) {
      return(data.frame(
        Analysis = label, DrugKey = drug, Group1 = NA, Group2 = NA,
        N_Group1 = NA, N_Group2 = NA, Median_Group1 = NA, Median_Group2 = NA,
        Mean_Group1 = NA, Mean_Group2 = NA, Delta_Median_Group2_minus_Group1 = NA,
        Wilcox_P = NA, Sensitive_Group = NA, stringsAsFactors = FALSE
      ))
    }

    g1 <- groups[1]
    g2 <- groups[2]
    x1 <- df$Predicted_IC50[df$Group == g1]
    x2 <- df$Predicted_IC50[df$Group == g2]
    p <- tryCatch(wilcox.test(x2, x1)$p.value, error = function(e) NA_real_)
    med1 <- median(x1, na.rm = TRUE)
    med2 <- median(x2, na.rm = TRUE)

    data.frame(
      Analysis = label,
      DrugKey = drug,
      Group1 = g1,
      Group2 = g2,
      N_Group1 = sum(df$Group == g1),
      N_Group2 = sum(df$Group == g2),
      Median_Group1 = med1,
      Median_Group2 = med2,
      Mean_Group1 = mean(x1, na.rm = TRUE),
      Mean_Group2 = mean(x2, na.rm = TRUE),
      Delta_Median_Group2_minus_Group1 = med2 - med1,
      Wilcox_P = p,
      Sensitive_Group = ifelse(med2 < med1, g2, g1),
      stringsAsFactors = FALSE
    )
  })

  out <- bind_rows(res)
  out$Wilcox_FDR <- p.adjust(out$Wilcox_P, method = "BH")
  out
}

spearman_drugs <- function(pred_mat, anno, score_col, label) {
  anno <- anno[match(rownames(pred_mat), anno$Sample), , drop = FALSE]

  res <- lapply(colnames(pred_mat), function(drug) {
    value <- as.numeric(pred_mat[, drug])
    score <- as.numeric(anno[[score_col]])
    keep <- is.finite(value) & is.finite(score)
    if (sum(keep) < 10) {
      return(data.frame(
        Analysis = label, DrugKey = drug, N = sum(keep),
        Spearman_Rho = NA, Spearman_P = NA, stringsAsFactors = FALSE
      ))
    }
    ct <- suppressWarnings(cor.test(score[keep], value[keep], method = "spearman"))
    data.frame(
      Analysis = label,
      DrugKey = drug,
      N = sum(keep),
      Spearman_Rho = unname(ct$estimate),
      Spearman_P = ct$p.value,
      stringsAsFactors = FALSE
    )
  })

  out <- bind_rows(res)
  out$Spearman_FDR <- p.adjust(out$Spearman_P, method = "BH")
  out
}

add_drug_meta <- function(x, meta) {
  x %>%
    left_join(meta, by = "DrugKey") %>%
    relocate(DrugName, DrugID, PutativeTarget, PathwayName, .after = DrugKey)
}

############################################################
# 3. Sample annotation
############################################################

cat("Project:", PROJECT_DIR, "\n")
cat("Loading sample annotation...\n")

mo_score <- fread(MO_SCORE_FILE, data.table = FALSE, check.names = FALSE)
mo_score$Sample <- gsub("\\.", "-", mo_score$Sample)
if (!"Patient" %in% colnames(mo_score)) {
  mo_score$Patient <- patient_id(mo_score$Sample)
} else {
  mo_score$Patient <- as.character(mo_score$Patient)
}
if (!"MO_DDRscore_raw" %in% colnames(mo_score) && "MO_DDRscore" %in% colnames(mo_score)) {
  mo_score$MO_DDRscore_raw <- mo_score$MO_DDRscore
}
mo_score <- mo_score %>%
  mutate(
    MO_DDRscore_raw = as.numeric(MO_DDRscore_raw),
    MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High"))
  ) %>%
  filter(MO_DDRscore_group %in% c("Low", "High"), is.finite(MO_DDRscore_raw))

dprs <- fread(DPRS_FILE, data.table = FALSE, check.names = FALSE) %>%
  mutate(
    Patient = as.character(ID),
    DPRS = as.numeric(DPRS),
    DPRS_RiskGroup = factor(RiskGroup, levels = c("Low", "High"))
  ) %>%
  filter(Dataset %in% c("Training", "Testing"), DPRS_RiskGroup %in% c("Low", "High")) %>%
  arrange(Patient, Dataset) %>%
  distinct(Patient, .keep_all = TRUE) %>%
  select(Patient, DPRS, DPRS_RiskGroup, DPRS_Dataset = Dataset)

sample_anno <- mo_score %>%
  filter(is_tumor(Sample)) %>%
  arrange(Patient, Sample) %>%
  distinct(Patient, .keep_all = TRUE) %>%
  left_join(dprs, by = "Patient") %>%
  select(
    Sample, Patient, MO_DDRscore_raw, MO_DDRscore_group,
    DPRS, DPRS_RiskGroup, DPRS_Dataset
  )

save_csv(sample_anno, file.path(DATA_OUT_DIR, "Drug_analysis_sample_annotation.csv"))

cat("Samples for MO-DDRscore:", sum(!is.na(sample_anno$MO_DDRscore_group)), "\n")
print(table(sample_anno$MO_DDRscore_group, useNA = "ifany"))
cat("Samples for DPRS:", sum(!is.na(sample_anno$DPRS_RiskGroup)), "\n")
print(table(sample_anno$DPRS_RiskGroup, useNA = "ifany"))

############################################################
# 4. TCGA test expression matrix
############################################################

cat("Loading TCGA expression...\n")

tcga_expr <- readRDS(TCGA_EXPR_FILE)
tcga_expr <- as.matrix(tcga_expr)
storage.mode(tcga_expr) <- "numeric"
colnames(tcga_expr) <- gsub("\\.", "-", colnames(tcga_expr))
rownames(tcga_expr) <- clean_gene(rownames(tcga_expr))
tcga_expr <- collapse_rows_by_mean(tcga_expr)

sample_anno <- sample_anno %>% filter(Sample %in% colnames(tcga_expr))
tcga_expr <- tcga_expr[, sample_anno$Sample, drop = FALSE]

if (max(tcga_expr, na.rm = TRUE) > 100) {
  tcga_expr <- log2(tcga_expr + 1)
}
tcga_expr[!is.finite(tcga_expr)] <- 0

cat("TCGA expression dim:", paste(dim(tcga_expr), collapse = " x "), "\n")

############################################################
# 5. GDSC2 / DepMap training matrices
############################################################

cat("Loading DepMap cell-line expression...\n")

cell_expr_dt <- fread(CELL_EXPR_FILE, data.table = TRUE, check.names = FALSE)
model_col <- "ModelID"
default_col <- "IsDefaultEntryForModel"

if (default_col %in% colnames(cell_expr_dt)) {
  cell_expr_dt <- cell_expr_dt[get(default_col) %in% c(TRUE, "TRUE", "True", "Yes", "YES", "1")]
}

cell_model_ids <- as.character(cell_expr_dt[[model_col]])
meta_cols <- intersect(
  c("", "SequencingID", "ModelID", "IsDefaultEntryForModel", "ModelConditionID", "IsDefaultEntryForMC"),
  colnames(cell_expr_dt)
)
gene_cols <- setdiff(colnames(cell_expr_dt), meta_cols)
gene_names <- clean_gene(gene_cols)
valid_gene_cols <- !is.na(gene_names) & gene_names != ""
gene_cols <- gene_cols[valid_gene_cols]
gene_names <- gene_names[valid_gene_cols]

cell_expr <- as.matrix(cell_expr_dt[, ..gene_cols])
storage.mode(cell_expr) <- "numeric"
rownames(cell_expr) <- cell_model_ids
colnames(cell_expr) <- gene_names
rm(cell_expr_dt)
gc()

training_expr <- t(cell_expr)
rm(cell_expr)
gc()
training_expr <- collapse_rows_by_mean(training_expr)
training_expr[!is.finite(training_expr)] <- 0

cat("DepMap expression dim:", paste(dim(training_expr), collapse = " x "), "\n")

cat("Loading GDSC2 response...\n")

model_map <- fread(MODEL_FILE, data.table = FALSE, check.names = FALSE) %>%
  mutate(
    ModelID = as.character(ModelID),
    SangerModelID = as.character(SangerModelID),
    COSMICID = as.character(COSMICID)
  ) %>%
  select(ModelID, SangerModelID, COSMICID, CellLineName, OncotreeLineage, OncotreePrimaryDisease) %>%
  distinct(ModelID, .keep_all = TRUE)

gdsc2 <- readxl::read_excel(GDSC2_FILE)
gdsc2 <- as.data.frame(gdsc2, stringsAsFactors = FALSE)
gdsc2 <- gdsc2 %>%
  mutate(
    SANGER_MODEL_ID = as.character(SANGER_MODEL_ID),
    COSMIC_ID = as.character(COSMIC_ID),
    DRUG_ID = as.character(DRUG_ID),
    LN_IC50 = as.numeric(LN_IC50),
    DRUG_NAME = as.character(DRUG_NAME),
    PUTATIVE_TARGET = as.character(PUTATIVE_TARGET),
    PATHWAY_NAME = as.character(PATHWAY_NAME)
  ) %>%
  filter(DATASET == "GDSC2", is.finite(LN_IC50), !is.na(DRUG_ID), !is.na(DRUG_NAME))

gdsc2_sid <- gdsc2 %>%
  left_join(model_map, by = c("SANGER_MODEL_ID" = "SangerModelID"))

gdsc2_unmapped <- gdsc2_sid %>% filter(is.na(ModelID))
if (nrow(gdsc2_unmapped) > 0) {
  gdsc2_cosmic <- gdsc2_unmapped %>%
    select(-ModelID, -COSMICID, -CellLineName, -OncotreeLineage, -OncotreePrimaryDisease) %>%
    left_join(model_map, by = c("COSMIC_ID" = "COSMICID"))
  gdsc2_mapped <- bind_rows(
    gdsc2_sid %>% filter(!is.na(ModelID)),
    gdsc2_cosmic %>% filter(!is.na(ModelID))
  )
} else {
  gdsc2_mapped <- gdsc2_sid %>% filter(!is.na(ModelID))
}

drug_meta <- gdsc2_mapped %>%
  mutate(
    DrugKey = paste0("Drug_", DRUG_ID, "__", safe_name(DRUG_NAME)),
    DrugName = DRUG_NAME,
    DrugID = DRUG_ID,
    PutativeTarget = PUTATIVE_TARGET,
    PathwayName = PATHWAY_NAME
  ) %>%
  group_by(DrugKey, DrugName, DrugID, PutativeTarget, PathwayName) %>%
  summarise(
    N_CellLines_GDSC2 = n_distinct(ModelID),
    Median_LN_IC50 = median(LN_IC50, na.rm = TRUE),
    .groups = "drop"
  )

priority_pattern <- paste(PRIORITY_DRUG_PATTERNS, collapse = "|")
drug_meta <- drug_meta %>%
  mutate(
    IsPriorityDrug = grepl(
      priority_pattern,
      toupper(paste(DrugName, PutativeTarget, PathwayName)),
      perl = TRUE
    )
  )

resp_long <- gdsc2_mapped %>%
  mutate(
    DrugKey = paste0("Drug_", DRUG_ID, "__", safe_name(DRUG_NAME)),
    IC50 = exp(LN_IC50)
  ) %>%
  filter(ModelID %in% colnames(training_expr), is.finite(IC50)) %>%
  group_by(ModelID, DrugKey) %>%
  summarise(IC50 = mean(IC50, na.rm = TRUE), .groups = "drop")

training_ptype_df <- resp_long %>%
  pivot_wider(names_from = DrugKey, values_from = IC50)

training_ptype <- as.data.frame(training_ptype_df)
rownames(training_ptype) <- training_ptype$ModelID
training_ptype$ModelID <- NULL
training_ptype <- as.matrix(training_ptype)
storage.mode(training_ptype) <- "numeric"

drug_n <- colSums(is.finite(training_ptype))
keep_drugs <- names(drug_n)[drug_n >= MIN_CELL_LINES_PER_DRUG]
if (RUN_PRIORITY_DRUGS_ONLY) {
  keep_drugs <- intersect(
    keep_drugs,
    drug_meta$DrugKey[drug_meta$IsPriorityDrug]
  )
}
training_ptype <- training_ptype[, keep_drugs, drop = FALSE]
drug_meta <- drug_meta %>% filter(DrugKey %in% colnames(training_ptype))

common_cell_lines <- intersect(colnames(training_expr), rownames(training_ptype))
training_expr <- training_expr[, common_cell_lines, drop = FALSE]
training_ptype <- training_ptype[common_cell_lines, , drop = FALSE]

common_genes <- intersect(rownames(training_expr), rownames(tcga_expr))
training_expr <- training_expr[common_genes, , drop = FALSE]
tcga_expr <- tcga_expr[common_genes, , drop = FALSE]

cat("Common genes:", length(common_genes), "\n")
cat("Common cell lines:", length(common_cell_lines), "\n")
cat("Drugs retained:", ncol(training_ptype), "\n")

save_csv(drug_meta, file.path(DATA_OUT_DIR, "GDSC2_drug_metadata_retained.csv"))
saveRDS(training_expr, file.path(DATA_OUT_DIR, "GDSC2_DepMap_training_expression_gene_by_cellline.rds"))
saveRDS(training_ptype, file.path(DATA_OUT_DIR, "GDSC2_training_IC50_cellline_by_drug.rds"))
saveRDS(tcga_expr, file.path(DATA_OUT_DIR, "TCGA_LUAD_test_expression_log2TPM_gene_by_sample.rds"))

training_summary <- data.frame(
  Item = c(
    "N_TCGA_samples",
    "N_training_cell_lines",
    "N_common_genes",
    "N_GDSC2_drugs_retained",
    "MIN_CELL_LINES_PER_DRUG",
    "RUN_PRIORITY_DRUGS_ONLY"
  ),
  Value = c(
    ncol(tcga_expr),
    ncol(training_expr),
    nrow(training_expr),
    ncol(training_ptype),
    MIN_CELL_LINES_PER_DRUG,
    RUN_PRIORITY_DRUGS_ONLY
  )
)
save_csv(training_summary, file.path(DATA_OUT_DIR, "Drug_training_matrix_summary.csv"))

############################################################
# 6. oncoPredict
############################################################

cat("Running oncoPredict::calcPhenotype...\n")

old_wd <- getwd()
setwd(ONCOPREDICT_DIR)
pred_ic50 <- oncoPredict::calcPhenotype(
  trainingExprData = training_expr,
  trainingPtype = training_ptype,
  testExprData = tcga_expr,
  batchCorrect = "standardize",
  powerTransformPhenotype = TRUE,
  removeLowVaryingGenes = 0.2,
  minNumSamples = 10,
  selection = 1,
  printOutput = TRUE,
  pcr = FALSE,
  removeLowVaringGenesFrom = "homogenizeData",
  report_pc = FALSE,
  cc = FALSE,
  percent = 80,
  rsq = FALSE,
  folder = FALSE
)
setwd(old_wd)

pred_ic50 <- as.matrix(pred_ic50)
storage.mode(pred_ic50) <- "numeric"
pred_ic50 <- pred_ic50[sample_anno$Sample, , drop = FALSE]

saveRDS(pred_ic50, file.path(RES_DIR, "Predicted_IC50_TCGA_all_drugs.rds"))
save_csv(
  data.frame(Sample = rownames(pred_ic50), pred_ic50, check.names = FALSE),
  file.path(TABLE_DIR, "Predicted_IC50_TCGA_all_drugs.csv")
)

############################################################
# 7. Group comparison and correlation
############################################################

cat("Running group comparisons...\n")

mo_anno <- sample_anno %>%
  filter(!is.na(MO_DDRscore_group), is.finite(MO_DDRscore_raw), Sample %in% rownames(pred_ic50))
mo_pred <- pred_ic50[mo_anno$Sample, , drop = FALSE]

dprs_anno <- sample_anno %>%
  filter(!is.na(DPRS_RiskGroup), is.finite(DPRS), Sample %in% rownames(pred_ic50))
dprs_pred <- pred_ic50[dprs_anno$Sample, , drop = FALSE]

mo_group <- wilcox_compare_drugs(
  pred_mat = mo_pred,
  anno = mo_anno,
  group_col = "MO_DDRscore_group",
  score_col = "MO_DDRscore_raw",
  label = "MO-DDRscore High vs Low"
) %>%
  add_drug_meta(drug_meta)

dprs_group <- wilcox_compare_drugs(
  pred_mat = dprs_pred,
  anno = dprs_anno,
  group_col = "DPRS_RiskGroup",
  score_col = "DPRS",
  label = "DPRS High vs Low"
) %>%
  add_drug_meta(drug_meta)

mo_cor <- spearman_drugs(
  pred_mat = mo_pred,
  anno = mo_anno,
  score_col = "MO_DDRscore_raw",
  label = "MO-DDRscore continuous"
) %>%
  add_drug_meta(drug_meta)

dprs_cor <- spearman_drugs(
  pred_mat = dprs_pred,
  anno = dprs_anno,
  score_col = "DPRS",
  label = "DPRS continuous"
) %>%
  add_drug_meta(drug_meta)

save_csv(mo_group, file.path(TABLE_DIR, "Drug_sensitivity_MO_DDRscore_group_comparison.csv"))
save_csv(dprs_group, file.path(TABLE_DIR, "Drug_sensitivity_DPRS_group_comparison.csv"))
save_csv(mo_cor, file.path(TABLE_DIR, "Drug_sensitivity_MO_DDRscore_spearman.csv"))
save_csv(dprs_cor, file.path(TABLE_DIR, "Drug_sensitivity_DPRS_spearman.csv"))

all_group <- bind_rows(mo_group, dprs_group)
all_cor <- bind_rows(mo_cor, dprs_cor)

save_csv(all_group, file.path(TABLE_DIR, "Drug_sensitivity_all_group_comparisons.csv"))
save_csv(all_cor, file.path(TABLE_DIR, "Drug_sensitivity_all_spearman_correlations.csv"))

priority_group <- all_group %>%
  filter(IsPriorityDrug | Wilcox_FDR < 0.05) %>%
  arrange(Analysis, Wilcox_FDR, Wilcox_P)
priority_cor <- all_cor %>%
  filter(IsPriorityDrug | Spearman_FDR < 0.05) %>%
  arrange(Analysis, Spearman_FDR, Spearman_P)

save_csv(priority_group, file.path(TABLE_DIR, "Drug_sensitivity_priority_and_significant_group_results.csv"))
save_csv(priority_cor, file.path(TABLE_DIR, "Drug_sensitivity_priority_and_significant_correlations.csv"))

summary_df <- data.frame(
  Item = c(
    "N_TCGA_samples_predicted",
    "N_MO_DDRscore_samples",
    "N_DPRS_samples",
    "N_drugs_predicted",
    "N_MO_DDRscore_FDR_lt_0.05",
    "N_DPRS_FDR_lt_0.05",
    "N_MO_DDRscore_priority_drugs",
    "N_DPRS_priority_drugs",
    "Output_directory"
  ),
  Value = c(
    nrow(pred_ic50),
    nrow(mo_anno),
    nrow(dprs_anno),
    ncol(pred_ic50),
    sum(mo_group$Wilcox_FDR < 0.05, na.rm = TRUE),
    sum(dprs_group$Wilcox_FDR < 0.05, na.rm = TRUE),
    sum(mo_group$IsPriorityDrug, na.rm = TRUE),
    sum(dprs_group$IsPriorityDrug, na.rm = TRUE),
    TABLE_DIR
  )
)

save_csv(summary_df, file.path(TABLE_DIR, "Drug_sensitivity_run_summary.csv"))

cat("\nDone.\n")
cat("Output directory:\n", TABLE_DIR, "\n")
print(summary_df)
