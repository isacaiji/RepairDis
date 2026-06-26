############################################################
# Manual DPRS machine-learning pipeline
# Only StepCox logic fixed; other model logic/parameters unchanged
#
# Core logic:
#   1. Train/Test split first
#   2. Training-only DEG as raw candidate genes
#   3. Training-only uni-Cox filtering, matching Mime1 unicox.filter.for.candi = TRUE
#   4. RSF selector uses minimal depth var.select(), matching Mime1-style RSF selector
#   5. StepCox is corrected:
#        forward: null model -> forward selection
#        both:    null model -> both direction selection
#        backward: full model -> backward selection
#   6. No manual deletion of models; final model ranked by Average C-index
############################################################

options(stringsAsFactors = FALSE)
set.seed(20260513)

############################################################
# 0. Parameters
############################################################

PROJECT_DIR <- "D:/R_workspace/评分/AD_DDR_project"
DATA_DIR <- file.path(PROJECT_DIR, "00_data")
PROC_DIR <- file.path(PROJECT_DIR, "01_processed")
DB_DIR <- file.path(PROJECT_DIR, "05_database_tables")
FIG4_DIR <- file.path(PROJECT_DIR, "04-5_Mime1Matched_StepCoxFixed")

dir.create(FIG4_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DB_DIR, recursive = TRUE, showWarnings = FALSE)

FOCUS_CANCER <- "LUAD"
TCGA_DIR <- "D:/R/R_workspace/梁老师文件/TCGA"

TCGA_EXPR_FILE <- file.path(
  TCGA_DIR,
  "mRNA_exp_TPM_only_TCGA/mRNA_exp_TPM_only_TCGA",
  paste0("TCGA-", FOCUS_CANCER, ".gene_expression_TPM.tsv")
)

TCGA_CLIN_PROCESSED_FILE <- file.path(
  PROC_DIR,
  paste0(FOCUS_CANCER, "_clinical_processed.csv")
)

MO_SCORE_FILE <- file.path(
  PROC_DIR,
  paste0(FOCUS_CANCER, "_MO_DDRscore.csv")
)

GSE72094_EXPR_FILE <- file.path(DATA_DIR, "GSE72094_expression.csv")
GSE72094_CLIN_FILE <- file.path(DATA_DIR, "GSE72094_clinical.csv")
GSE68465_EXPR_FILE <- file.path(DATA_DIR, "GSE68465_expression.csv")
GSE68465_CLIN_FILE <- file.path(DATA_DIR, "GSE68465_clinical.csv")

TRAIN_RATIO <- 0.80
DEG_ADJ_P <- 0.05
DEG_LOGFC <- 1.00
UNICOX_P <- 0.02
MAX_ML_GENES <- 120
NODE_SIZE <- 15
SEED <- 20260513

COX_P_CUTOFF <- 0.02

RF_NTREE <- 3000
RF_NODESIZE <- NODE_SIZE

RSF_IMP_CUTOFF <- 0.20

STEP_MAX_GENES <- 60
USE_STEP_MAX_GENES <- FALSE

GBM_NTREES <- 10000
GBM_DEPTH <- 3
GBM_MINOBSINNODE <- 10
GBM_SHRINKAGE <- 0.001
CV_FOLDS <- 10

STEPCOX_K_TYPE <- "AIC"
STEPCOX_K_NUMERIC <- 2

set.seed(SEED)

############################################################
# 1. Packages
############################################################

pkgs <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "ggplot2",
  "survival",
  "limma",
  "randomForestSRC",
  "glmnet",
  "plsRcox",
  "gbm",
  "CoxBoost",
  "survivalsvm",
  "superpc",
  "survminer",
  "timeROC"
)

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(survival)
  library(limma)
  library(randomForestSRC)
  library(glmnet)
  library(plsRcox)
  library(gbm)
  library(CoxBoost)
  library(survivalsvm)
  library(superpc)
  library(survminer)
  library(timeROC)
})

############################################################
# 2. Helper functions
############################################################

save_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file)
}

safe_ggsave <- function(file, p, w = 7, h = 5) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = file,
    plot = p,
    width = w,
    height = h,
    device = "pdf",
    useDingbats = FALSE
  )
}

clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub(" ///.*$| //.*$|;.*$|,.*$|\\s*\\([^\\)]*\\)$", "", x)
  x[x %in% c("", "NA", "---", "NULL", "N/A", "NAN")] <- NA
  x
}

patient_id <- function(x) {
  substr(gsub("\\.", "-", as.character(x)), 1, 12)
}

sample_type <- function(x) {
  substr(gsub("\\.", "-", as.character(x)), 14, 15)
}

sample_class <- function(x) {
  code <- sample_type(x)
  ifelse(
    code %in% c("01", "02", "03", "05", "06", "07"),
    "Tumor",
    ifelse(code %in% c("10", "11", "12", "13", "14"), "Normal", "Other")
  )
}

zscore_vector <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

calc_cindex_risk <- function(time, status, score) {
  df <- data.frame(
    time = as.numeric(time),
    status = as.numeric(status),
    score = as.numeric(score)
  )
  df <- df[is.finite(df$time) & is.finite(df$status) & is.finite(df$score), ]
  if (nrow(df) < 30 || sum(df$status == 1) < 5) return(NA_real_)
  
  as.numeric(
    survival::concordance(
      survival::Surv(time, status) ~ score,
      data = df,
      reverse = TRUE
    )$concordance
  )
}

get_time_points <- function(time) {
  if (median(time, na.rm = TRUE) > 100) {
    c(365, 1095, 1825)
  } else {
    c(1, 3, 5)
  }
}

model_complexity <- function(x) {
  sapply(strsplit(x, "\\+"), length)
}

############################################################
# 3. Load required objects if missing
############################################################

TCGA_EXPR_RDS <- file.path(
  PROC_DIR,
  paste0(FOCUS_CANCER, "_tcga_expr_tpm_matrix.rds")
)

if (exists("tcga_expr")) {
  cat("tcga_expr already exists in environment. Skip reading expression file.\n")
} else if (file.exists(TCGA_EXPR_RDS)) {
  cat("tcga_expr not found in environment. Reading cached RDS...\n")
  tcga_expr <- readRDS(TCGA_EXPR_RDS)
} else {
  cat("tcga_expr not found and cached RDS not found. Reading raw expression file...\n")
  
  expr_raw <- data.table::fread(
    TCGA_EXPR_FILE,
    data.table = FALSE,
    check.names = FALSE
  )
  
  gene_col <- colnames(expr_raw)[1]
  expr_raw[[gene_col]] <- clean_gene(expr_raw[[gene_col]])
  
  expr_raw <- expr_raw %>%
    dplyr::filter(!is.na(.data[[gene_col]]), .data[[gene_col]] != "") %>%
    dplyr::group_by(.data[[gene_col]]) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::everything(),
        ~ suppressWarnings(mean(as.numeric(.x), na.rm = TRUE))
      ),
      .groups = "drop"
    )
  
  gene_vec <- expr_raw[[gene_col]]
  expr_raw[[gene_col]] <- NULL
  
  tcga_expr <- as.matrix(expr_raw)
  rownames(tcga_expr) <- gene_vec
  storage.mode(tcga_expr) <- "numeric"
  
  saveRDS(tcga_expr, TCGA_EXPR_RDS)
  cat("tcga_expr cached to:", TCGA_EXPR_RDS, "\n")
}

colnames(tcga_expr) <- gsub("\\.", "-", colnames(tcga_expr))

if (!exists("tumor_samples")) {
  cat("tumor_samples not found. Generating from tcga_expr columns...\n")
  tumor_samples <- colnames(tcga_expr)[sample_class(colnames(tcga_expr)) == "Tumor"]
} else {
  cat("tumor_samples already exists. Skip generating tumor_samples.\n")
  tumor_samples <- gsub("\\.", "-", tumor_samples)
}

if (!exists("tcga_clin")) {
  cat("tcga_clin not found. Reading processed clinical file...\n")
  tcga_clin <- data.table::fread(
    TCGA_CLIN_PROCESSED_FILE,
    data.table = FALSE,
    check.names = FALSE
  )
} else {
  cat("tcga_clin already exists. Skip reading clinical file.\n")
}

tcga_clin <- tcga_clin %>%
  dplyr::mutate(
    Patient = as.character(Patient),
    time = as.numeric(time),
    status = as.numeric(status)
  ) %>%
  dplyr::select(Patient, time, status, age, gender, stage)

if (!exists("mo_score")) {
  cat("mo_score not found. Reading MO-DDRscore file...\n")
  mo_score <- data.table::fread(
    MO_SCORE_FILE,
    data.table = FALSE,
    check.names = FALSE
  )
} else {
  cat("mo_score already exists. Skip reading MO-DDRscore file.\n")
}

mo_score <- mo_score %>%
  dplyr::mutate(
    Sample = gsub("\\.", "-", as.character(Sample)),
    Patient = as.character(Patient),
    SampleType = sample_type(Sample),
    MO_DDRscore_raw = as.numeric(MO_DDRscore_raw),
    MO_DDRscore_group = as.character(MO_DDRscore_group)
  ) %>%
  dplyr::select(Sample, Patient, SampleType, MO_DDRscore_group, MO_DDRscore_raw)

cat("\nLoaded objects summary:\n")
cat("tcga_expr:", nrow(tcga_expr), "genes x", ncol(tcga_expr), "samples\n")
cat("tumor_samples:", length(tumor_samples), "\n")
cat("tcga_clin:", nrow(tcga_clin), "patients\n")
cat("mo_score:", nrow(mo_score), "samples\n")

############################################################
# 4. Remove old outputs
############################################################

old_files <- c(
  "ManualML_all_model_Cindex_long.csv",
  "Fig4B_all_model_Cindex_summary.csv",
  "Fig4B_final_selected_model_info.csv",
  "Fig4B_Cindex_heatmap.pdf",
  "Fig4C_DPRS_all_sets.csv",
  "Fig4C_final_signature_genes.csv",
  "Fig4C_candidate_gene_lollipop.pdf",
  "Fig4D_risk_score_distribution.pdf",
  "Fig4D_survival_status_distribution.pdf",
  "Fig4D_selected_model_Cindex_by_dataset.csv",
  "ManualML_combo_vs_base_risk_equal_check.csv",
  "ManualML_model_feature_counts.csv"
)

unlink(file.path(FIG4_DIR, old_files), force = TRUE)

unlink(
  list.files(
    FIG4_DIR,
    pattern = "Fig4E_.*_KM\\.pdf|Fig4F_.*_timeROC\\.pdf|Fig4F_.*_timeROC_AUC\\.csv",
    full.names = TRUE
  ),
  force = TRUE
)

############################################################
# 5. Build score survival table
############################################################

score_df <- mo_score %>%
  dplyr::filter(
    Sample %in% tumor_samples,
    !is.na(MO_DDRscore_group),
    is.finite(MO_DDRscore_raw)
  ) %>%
  dplyr::arrange(Patient, SampleType) %>%
  dplyr::distinct(Patient, .keep_all = TRUE) %>%
  dplyr::inner_join(tcga_clin, by = "Patient") %>%
  dplyr::filter(
    is.finite(time),
    time > 0,
    status %in% c(0, 1)
  )

save_csv(
  score_df,
  file.path(FIG4_DIR, "Fig4A_TCGA_score_survival_table.csv")
)

cat("\nTCGA score survival table:\n")
cat("Patients:", nrow(score_df), "\n")
cat("Events:", sum(score_df$status == 1), "\n")

############################################################
# 6. Train/Test split first
############################################################

set.seed(SEED)

event_ids <- score_df$Patient[score_df$status == 1]
cens_ids <- score_df$Patient[score_df$status == 0]

train_ids <- c(
  sample(event_ids, floor(TRAIN_RATIO * length(event_ids))),
  sample(cens_ids, floor(TRAIN_RATIO * length(cens_ids)))
)

test_ids <- setdiff(score_df$Patient, train_ids)

train_meta <- score_df %>%
  dplyr::filter(Patient %in% train_ids) %>%
  dplyr::arrange(Patient)

test_meta <- score_df %>%
  dplyr::filter(Patient %in% test_ids) %>%
  dplyr::arrange(Patient)

train_samples <- train_meta$Sample
test_samples <- test_meta$Sample

cat("\nTraining:", nrow(train_meta), "patients /", sum(train_meta$status == 1), "events\n")
cat("Testing:", nrow(test_meta), "patients /", sum(test_meta$status == 1), "events\n")

############################################################
# 7. Training-only DEG
############################################################

train_group <- factor(train_meta$MO_DDRscore_group, levels = c("Low", "High"))
expr_train <- log2(tcga_expr[, train_samples, drop = FALSE] + 1)

design <- model.matrix(~0 + train_group)
colnames(design) <- c("Low", "High")

fit <- limma::lmFit(expr_train, design)
cont <- limma::makeContrasts(High - Low, levels = design)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont))

train_deg <- limma::topTable(fit2, number = Inf, adjust.method = "BH")
train_deg$Gene <- rownames(train_deg)
train_deg$FDR <- train_deg$adj.P.Val

save_csv(
  train_deg,
  file.path(FIG4_DIR, "Fig4A_training_only_DEG.csv")
)

deg_genes <- train_deg %>%
  dplyr::filter(
    is.finite(FDR),
    FDR < DEG_ADJ_P,
    abs(logFC) > DEG_LOGFC
  ) %>%
  dplyr::arrange(FDR, dplyr::desc(abs(logFC))) %>%
  dplyr::pull(Gene) %>%
  unique()

cat("\nTraining-only DEG genes:", length(deg_genes), "\n")

############################################################
# 8. GEO common genes
############################################################

read_geo_genes <- function(expr_file) {
  x <- data.table::fread(expr_file, data.table = FALSE, check.names = FALSE)
  unique(na.omit(clean_gene(x[[1]])))
}

geo_common_genes <- Reduce(
  intersect,
  list(
    rownames(tcga_expr),
    read_geo_genes(GSE72094_EXPR_FILE),
    read_geo_genes(GSE68465_EXPR_FILE)
  )
)

############################################################
# 9. Candidate genes = DEG raw candidate + uni-Cox filter
# Matching Mime1:
#   raw candidate genes -> univariable Cox filter P < 0.05
############################################################

cat("\nSelecting candidate genes by training-only DEG + training-only uni-Cox...\n")

deg_pool <- train_deg %>%
  dplyr::filter(
    Gene %in% geo_common_genes,
    is.finite(FDR),
    FDR < DEG_ADJ_P,
    abs(logFC) > DEG_LOGFC
  ) %>%
  dplyr::arrange(FDR, dplyr::desc(abs(logFC))) %>%
  dplyr::pull(Gene) %>%
  unique()

cat("Raw candidate genes from training-only DEG:", length(deg_pool), "\n")

if (length(deg_pool) < 10) {
  stop("Too few DEG genes after filtering. Please check DEG_ADJ_P / DEG_LOGFC.")
}

train_expr_cox <- as.data.frame(
  t(log2(tcga_expr[deg_pool, train_samples, drop = FALSE] + 1)),
  check.names = FALSE
)

train_expr_cox$Sample <- rownames(train_expr_cox)
train_expr_cox$Patient <- patient_id(train_expr_cox$Sample)

train_surv <- train_expr_cox %>%
  dplyr::inner_join(
    tcga_clin[, c("Patient", "time", "status")],
    by = "Patient"
  ) %>%
  dplyr::filter(
    is.finite(time),
    time > 0,
    status %in% c(0, 1)
  )

cox_table <- dplyr::bind_rows(lapply(deg_pool, function(g) {
  
  df <- train_surv[, c("time", "status", g), drop = FALSE]
  colnames(df) <- c("time", "status", "expr")
  
  df <- df %>%
    dplyr::filter(
      is.finite(time),
      time > 0,
      status %in% c(0, 1),
      is.finite(expr)
    )
  
  if (nrow(df) < 30 || sum(df$status == 1) < 5 || sd(df$expr, na.rm = TRUE) == 0) {
    return(NULL)
  }
  
  fit <- tryCatch(
    survival::coxph(survival::Surv(time, status) ~ expr, data = df),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NULL)
  
  sm <- summary(fit)
  
  data.frame(
    Gene = g,
    CoxCoef = sm$coefficients[1, "coef"],
    HR = sm$coefficients[1, "exp(coef)"],
    CoxP = sm$coefficients[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}))

cox_table$CoxFDR <- p.adjust(cox_table$CoxP, method = "BH")
cox_table <- cox_table %>% dplyr::arrange(CoxP)

save_csv(
  cox_table,
  file.path(FIG4_DIR, "Fig4A_training_only_unicox_for_candidate_selection.csv")
)

candidate_genes <- cox_table %>%
  dplyr::filter(
    is.finite(CoxP),
    CoxP < COX_P_CUTOFF
  ) %>%
  dplyr::arrange(CoxP) %>%
  dplyr::pull(Gene) %>%
  unique()

cat("Uni-Cox filtered candidate genes:", length(candidate_genes), "\n")

if (length(candidate_genes) < 10) {
  stop("Too few genes after uni-Cox filtering. Please check COX_P_CUTOFF.")
}

if (length(candidate_genes) > MAX_ML_GENES) {
  candidate_genes <- candidate_genes[1:MAX_ML_GENES]
}

save_csv(
  data.frame(Gene = candidate_genes),
  file.path(FIG4_DIR, "Fig4A_candidate_genes_DEG_unicox_training_only.csv")
)

cat("Candidate genes for manual ML:", length(candidate_genes), "\n")
print(candidate_genes)

############################################################
# 10. Build TCGA and GEO ML input
############################################################

make_tcga_ml_input <- function(meta_df, gene_use) {
  
  samples <- meta_df$Sample
  
  x <- as.data.frame(
    t(log2(tcga_expr[gene_use, samples, drop = FALSE] + 1)),
    check.names = FALSE
  )
  
  x$Sample <- rownames(x)
  x$Patient <- patient_id(x$Sample)
  
  x <- x %>%
    dplyr::arrange(Patient, Sample) %>%
    dplyr::distinct(Patient, .keep_all = TRUE) %>%
    dplyr::inner_join(tcga_clin[, c("Patient", "time", "status")], by = "Patient") %>%
    dplyr::filter(is.finite(time), time > 0, status %in% c(0, 1))
  
  out <- x[, c("Patient", "time", "status", gene_use)]
  colnames(out)[1:3] <- c("ID", "OS.time", "OS")
  
  for (g in gene_use) {
    out[[g]] <- zscore_vector(out[[g]])
  }
  
  out
}

prepare_geo_ml_input <- function(expr_file, clin_file, gene_use, cohort_name) {
  
  expr <- data.table::fread(expr_file, data.table = FALSE, check.names = FALSE)
  clin <- data.table::fread(clin_file, data.table = FALSE, check.names = FALSE)
  
  gene_col <- colnames(expr)[1]
  expr[[gene_col]] <- clean_gene(expr[[gene_col]])
  
  expr <- expr[!is.na(expr[[gene_col]]) & expr[[gene_col]] != "", ]
  expr <- expr[!duplicated(expr[[gene_col]]), ]
  
  rownames(expr) <- expr[[gene_col]]
  expr[[gene_col]] <- NULL
  
  expr_mat <- as.matrix(expr)
  storage.mode(expr_mat) <- "numeric"
  
  common_genes <- intersect(gene_use, rownames(expr_mat))
  
  x <- as.data.frame(
    t(expr_mat[common_genes, , drop = FALSE]),
    check.names = FALSE
  )
  
  x$ID <- rownames(x)
  clin$Sample <- as.character(clin$Sample)
  
  out <- x %>%
    dplyr::inner_join(clin[, c("Sample", "time", "status")], by = c("ID" = "Sample")) %>%
    dplyr::filter(is.finite(time), time > 0, status %in% c(0, 1))
  
  missing_genes <- setdiff(gene_use, colnames(out))
  for (g in missing_genes) {
    out[[g]] <- 0
  }
  
  out <- out[, c("ID", "time", "status", gene_use)]
  colnames(out)[1:3] <- c("ID", "OS.time", "OS")
  
  for (g in gene_use) {
    out[[g]] <- zscore_vector(out[[g]])
  }
  
  cat(
    cohort_name, ":",
    nrow(out), "samples /",
    sum(out$OS == 1), "events /",
    length(common_genes), "matched genes\n"
  )
  
  out
}

train_df <- make_tcga_ml_input(train_meta, candidate_genes)
test_df <- make_tcga_ml_input(test_meta, candidate_genes)

gse72094_df <- prepare_geo_ml_input(
  GSE72094_EXPR_FILE,
  GSE72094_CLIN_FILE,
  candidate_genes,
  "GSE72094"
)

gse68465_df <- prepare_geo_ml_input(
  GSE68465_EXPR_FILE,
  GSE68465_CLIN_FILE,
  candidate_genes,
  "GSE68465"
)

save_csv(train_df, file.path(FIG4_DIR, "ManualML_input_Training.csv"))
save_csv(test_df, file.path(FIG4_DIR, "ManualML_input_Testing.csv"))
save_csv(gse72094_df, file.path(FIG4_DIR, "ManualML_input_GSE72094.csv"))
save_csv(gse68465_df, file.path(FIG4_DIR, "ManualML_input_GSE68465.csv"))

############################################################
# 11. Manual ML preparation
############################################################

gene_map <- data.frame(
  Gene = candidate_genes,
  SafeGene = make.names(candidate_genes, unique = TRUE),
  stringsAsFactors = FALSE
)

safe_genes <- gene_map$SafeGene
names(safe_genes) <- gene_map$Gene

rename_ml_genes <- function(df, gene_map) {
  df <- as.data.frame(df, check.names = FALSE)
  idx <- match(gene_map$Gene, colnames(df))
  keep <- !is.na(idx)
  colnames(df)[idx[keep]] <- gene_map$SafeGene[keep]
  df
}

train_ml <- rename_ml_genes(train_df, gene_map)
test_ml <- rename_ml_genes(test_df, gene_map)
gse72094_ml <- rename_ml_genes(gse72094_df, gene_map)
gse68465_ml <- rename_ml_genes(gse68465_df, gene_map)

val_ml_list <- list(
  Training = train_ml,
  Testing = test_ml,
  GSE72094 = gse72094_ml,
  GSE68465 = gse68465_ml
)

cox_anno <- cox_table %>%
  dplyr::left_join(gene_map, by = "Gene")

cox_rank_safe <- cox_anno %>%
  dplyr::filter(SafeGene %in% unname(safe_genes), is.finite(CoxP)) %>%
  dplyr::arrange(CoxP) %>%
  dplyr::pull(SafeGene) %>%
  unique()

safe_gene_back <- function(safe_vec) {
  gene_map$Gene[match(safe_vec, gene_map$SafeGene)]
}

make_model_df <- function(df, genes) {
  out <- df[, c("OS.time", "OS", genes), drop = FALSE]
  out$OS.time <- as.numeric(out$OS.time)
  out$OS <- as.numeric(out$OS)
  out
}

make_val_df <- function(df, genes) {
  out <- df[, c("ID", "OS.time", "OS", genes), drop = FALSE]
  out$OS.time <- as.numeric(out$OS.time)
  out$OS <- as.numeric(out$OS)
  out
}

cindex_from_rs <- function(dat) {
  
  dat <- dat %>%
    dplyr::filter(
      is.finite(OS.time),
      OS.time > 0,
      OS %in% c(0, 1),
      is.finite(RS)
    )
  
  if (nrow(dat) < 30 || sum(dat$OS == 1) < 5) return(NA_real_)
  
  fit <- tryCatch(
    survival::coxph(survival::Surv(OS.time, OS) ~ RS, data = dat),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NA_real_)
  
  as.numeric(summary(fit)$concordance[1])
}

eval_rs_list <- function(rs_list, model_name) {
  dplyr::bind_rows(lapply(names(rs_list), function(ds) {
    data.frame(
      Dataset = ds,
      Model = model_name,
      Cindex = cindex_from_rs(rs_list[[ds]]),
      stringsAsFactors = FALSE
    )
  }))
}

predict_to_rs <- function(val_df, rs) {
  data.frame(
    ID = as.character(val_df$ID),
    OS.time = as.numeric(val_df$OS.time),
    OS = as.numeric(val_df$OS),
    RS = as.numeric(rs),
    stringsAsFactors = FALSE
  )
}

result_long <- data.frame()
risk_store <- list()
feature_store <- list()
selector_store <- list()

add_result <- function(model_name, rs_list, used_genes) {
  
  tmp <- eval_rs_list(rs_list, model_name)
  result_long <<- dplyr::bind_rows(result_long, tmp)
  
  risk_store[[model_name]] <<- dplyr::bind_rows(lapply(names(rs_list), function(ds) {
    x <- rs_list[[ds]]
    data.frame(
      ID = as.character(x$ID),
      time = as.numeric(x$OS.time),
      status = as.numeric(x$OS),
      DPRS = as.numeric(x$RS),
      Dataset = ds,
      Model = model_name,
      stringsAsFactors = FALSE
    )
  }))
  
  feature_store[[model_name]] <<- unique(used_genes)
}

############################################################
# 12. Learner functions
############################################################

fit_predict_rsf <- function(genes) {
  
  est <- make_model_df(train_ml, genes)
  
  fit <- randomForestSRC::rfsrc(
    Surv(OS.time, OS) ~ .,
    data = est,
    ntree = RF_NTREE,
    nodesize = RF_NODESIZE,
    splitrule = "logrank",
    importance = TRUE,
    proximity = FALSE,
    forest = TRUE,
    seed = SEED
  )
  
  rs_list <- lapply(names(val_ml_list), function(ds) {
    v <- make_val_df(val_ml_list[[ds]], genes)
    
    pred <- predict(
      fit,
      newdata = v[, c("OS.time", "OS", genes), drop = FALSE]
    )$predicted
    
    predict_to_rs(v, pred)
  })
  
  names(rs_list) <- names(val_ml_list)
  rs_list
}

fit_predict_glmnet <- function(genes, alpha_value) {
  
  est <- make_model_df(train_ml, genes)
  x <- as.matrix(est[, genes, drop = FALSE])
  y <- survival::Surv(est$OS.time, est$OS)
  
  set.seed(SEED)
  
  fit <- glmnet::cv.glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = alpha_value,
    nfolds = CV_FOLDS
  )
  
  rs_list <- lapply(names(val_ml_list), function(ds) {
    v <- make_val_df(val_ml_list[[ds]], genes)
    rs <- as.numeric(
      predict(
        fit,
        type = "link",
        newx = as.matrix(v[, genes, drop = FALSE]),
        s = fit$lambda.min
      )
    )
    predict_to_rs(v, rs)
  })
  
  names(rs_list) <- names(val_ml_list)
  
  if (alpha_value > 0) {
    cf <- glmnet::coef.glmnet(fit, s = fit$lambda.min)
    selected_genes <- rownames(cf)[as.numeric(cf) != 0]
    selected_genes <- unique(selected_genes)
    
    if (length(selected_genes) >= 2) {
      attr(rs_list, "selected_genes") <- selected_genes
    }
  }
  
  rs_list
}

fit_predict_coxboost <- function(genes) {
  
  est <- make_model_df(train_ml, genes)
  x <- as.matrix(est[, genes, drop = FALSE])
  
  pen <- CoxBoost::optimCoxBoostPenalty(
    time = est$OS.time,
    status = est$OS,
    x = x,
    trace = FALSE,
    start.penalty = 500,
    parallel = FALSE
  )
  
  cv.res <- CoxBoost::cv.CoxBoost(
    time = est$OS.time,
    status = est$OS,
    x = x,
    maxstepno = 500,
    K = CV_FOLDS,
    type = "verweij",
    penalty = pen$penalty
  )
  
  fit <- CoxBoost::CoxBoost(
    time = est$OS.time,
    status = est$OS,
    x = x,
    stepno = cv.res$optimal.step,
    penalty = pen$penalty
  )
  
  rs_list <- lapply(names(val_ml_list), function(ds) {
    
    v <- make_val_df(val_ml_list[[ds]], genes)
    
    rs <- as.numeric(
      predict(
        fit,
        newdata = as.matrix(v[, genes, drop = FALSE]),
        newtime = v$OS.time,
        newstatus = v$OS,
        type = "lp"
      )
    )
    
    predict_to_rs(v, rs)
  })
  
  names(rs_list) <- names(val_ml_list)
  
  cf <- coef(fit)
  
  if (is.null(names(cf))) {
    names(cf) <- genes[seq_along(cf)]
  }
  
  selected_genes <- names(cf)[as.numeric(cf) != 0]
  selected_genes <- unique(selected_genes)
  
  if (length(selected_genes) >= 2) {
    attr(rs_list, "selected_genes") <- selected_genes
  }
  
  rs_list
}

fit_predict_plsrcox <- function(genes) {
  
  est <- make_model_df(train_ml, genes)
  nt_use <- min(3, length(genes))
  
  fit <- plsRcox::plsRcox(
    Xplan = est[, genes, drop = FALSE],
    time = est$OS.time,
    event = est$OS,
    nt = nt_use
  )
  
  rs_list <- lapply(names(val_ml_list), function(ds) {
    v <- make_val_df(val_ml_list[[ds]], genes)
    rs <- as.numeric(
      predict(
        fit,
        type = "lp",
        newdata = v[, genes, drop = FALSE]
      )
    )
    predict_to_rs(v, rs)
  })
  
  names(rs_list) <- names(val_ml_list)
  rs_list
}

fit_predict_superpc <- function(genes) {
  
  est <- make_model_df(train_ml, genes)
  
  data_spc <- list(
    x = t(est[, genes, drop = FALSE]),
    y = est$OS.time,
    censoring.status = est$OS,
    featurenames = genes
  )
  
  fit <- superpc::superpc.train(
    data = data_spc,
    type = "survival",
    s0.perc = 0.5
  )
  
  cv.fit <- superpc::superpc.cv(
    fit,
    data_spc,
    n.threshold = 20,
    n.fold = CV_FOLDS,
    n.components = 3,
    min.features = min(5, length(genes)),
    max.features = length(genes),
    compute.fullcv = TRUE,
    compute.preval = TRUE
  )
  
  thr <- cv.fit$thresholds[which.max(cv.fit$scor[1, ])]
  
  rs_list <- lapply(names(val_ml_list), function(ds) {
    v <- make_val_df(val_ml_list[[ds]], genes)
    
    test_spc <- list(
      x = t(v[, genes, drop = FALSE]),
      y = v$OS.time,
      censoring.status = v$OS,
      featurenames = genes
    )
    
    pred <- superpc::superpc.predict(
      fit,
      data_spc,
      test_spc,
      threshold = thr,
      n.components = 1
    )
    
    predict_to_rs(v, as.numeric(pred$v.pred))
  })
  
  names(rs_list) <- names(val_ml_list)
  rs_list
}

fit_predict_gbm <- function(genes) {
  
  est <- make_model_df(train_ml, genes)
  
  fit0 <- gbm::gbm(
    formula = Surv(OS.time, OS) ~ .,
    data = est,
    distribution = "coxph",
    n.trees = GBM_NTREES,
    interaction.depth = GBM_DEPTH,
    n.minobsinnode = GBM_MINOBSINNODE,
    shrinkage = GBM_SHRINKAGE,
    cv.folds = CV_FOLDS,
    n.cores = 4,
    verbose = FALSE
  )
  
  best <- which.min(fit0$cv.error)
  if (!is.finite(best) || best < 1) best <- GBM_NTREES
  
  fit <- gbm::gbm(
    formula = Surv(OS.time, OS) ~ .,
    data = est,
    distribution = "coxph",
    n.trees = best,
    interaction.depth = GBM_DEPTH,
    n.minobsinnode = GBM_MINOBSINNODE,
    shrinkage = GBM_SHRINKAGE,
    cv.folds = CV_FOLDS,
    n.cores = 4,
    verbose = FALSE
  )
  
  rs_list <- lapply(names(val_ml_list), function(ds) {
    v <- make_val_df(val_ml_list[[ds]], genes)
    rs <- as.numeric(
      predict(
        fit,
        newdata = v[, c("OS.time", "OS", genes), drop = FALSE],
        n.trees = best,
        type = "link"
      )
    )
    predict_to_rs(v, rs)
  })
  
  names(rs_list) <- names(val_ml_list)
  rs_list
}

fit_predict_svm <- function(genes) {
  
  est <- make_model_df(train_ml, genes)
  
  fit <- survivalsvm::survivalsvm(
    Surv(OS.time, OS) ~ .,
    data = est,
    gamma.mu = 1
  )
  
  rs_list <- lapply(names(val_ml_list), function(ds) {
    v <- make_val_df(val_ml_list[[ds]], genes)
    rs <- as.numeric(
      predict(
        fit,
        v[, c("OS.time", "OS", genes), drop = FALSE]
      )$predicted
    )
    predict_to_rs(v, rs)
  })
  
  names(rs_list) <- names(val_ml_list)
  rs_list
}

get_stepcox_genes <- function(genes) {
  
  genes_use <- intersect(cox_rank_safe, genes)
  
  if (USE_STEP_MAX_GENES) {
    genes_use <- head(genes_use, min(STEP_MAX_GENES, length(genes_use)))
  }
  
  genes_use
}

get_stepcox_k <- function(est) {
  
  if (STEPCOX_K_TYPE == "BIC") {
    return(log(nrow(est)))
  }
  
  if (STEPCOX_K_TYPE == "AIC") {
    return(STEPCOX_K_NUMERIC)
  }
  
  stop("Unknown STEPCOX_K_TYPE. Please use 'AIC' or 'BIC'.")
}

fit_predict_stepcox <- function(genes, direction) {
  
  genes_use <- get_stepcox_genes(genes)
  
  if (length(genes_use) < 2) {
    stop("Too few genes for StepCox.")
  }
  
  est <- make_model_df(train_ml, genes_use)
  step_k <- get_stepcox_k(est)
  
  if (direction == "forward") {
    
    null_fit <- survival::coxph(
      Surv(OS.time, OS) ~ 1,
      data = est
    )
    
    upper_formula <- as.formula(
      paste0("~ ", paste(genes_use, collapse = " + "))
    )
    
    fit <- step(
      null_fit,
      scope = list(lower = ~1, upper = upper_formula),
      direction = "forward",
      trace = 0,
      k = step_k
    )
    
  } else if (direction == "backward") {
    
    full_fit <- survival::coxph(
      Surv(OS.time, OS) ~ .,
      data = est
    )
    
    fit <- step(
      full_fit,
      direction = "backward",
      trace = 0,
      k = step_k
    )
    
  } else if (direction == "both") {
    
    null_fit <- survival::coxph(
      Surv(OS.time, OS) ~ 1,
      data = est
    )
    
    upper_formula <- as.formula(
      paste0("~ ", paste(genes_use, collapse = " + "))
    )
    
    fit <- step(
      null_fit,
      scope = list(lower = ~1, upper = upper_formula),
      direction = "both",
      trace = 0,
      k = step_k
    )
    
  } else {
    stop("Unknown StepCox direction: ", direction)
  }
  coef_vec <- stats::coef(fit)
  selected_genes <- names(coef_vec)[is.finite(coef_vec)]
  selected_genes <- unique(selected_genes)
  
  if (length(selected_genes) < 2) {
    stop(paste0("StepCox[", direction, "] selected fewer than 2 genes."))
  }
  
  rs_list <- lapply(names(val_ml_list), function(ds) {
    
    v <- make_val_df(val_ml_list[[ds]], selected_genes)
    
    rs <- as.numeric(
      predict(
        fit,
        type = "risk",
        newdata = v[, c("OS.time", "OS", selected_genes), drop = FALSE]
      )
    )
    
    predict_to_rs(v, rs)
  })
  
  names(rs_list) <- names(val_ml_list)
  attr(rs_list, "selected_genes") <- selected_genes
  
  rs_list
}

############################################################
# 13. Selector functions
############################################################

select_by_rsf <- function(genes) {
  
  est <- make_model_df(train_ml, genes)
  
  fit <- randomForestSRC::rfsrc(
    Surv(OS.time, OS) ~ .,
    data = est,
    ntree = RF_NTREE,
    nodesize = RF_NODESIZE,
    splitrule = "logrank",
    importance = TRUE,
    proximity = FALSE,
    forest = TRUE,
    seed = SEED
  )
  
  vs <- tryCatch(
    randomForestSRC::var.select(
      object = fit,
      conservative = "high",
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  
  selected <- NULL
  
  if (!is.null(vs) && "topvars" %in% names(vs)) {
    selected <- vs$topvars
  }
  
  selected <- unique(selected)
  selected <- selected[selected %in% genes]
  
  selector_method <- "RSF_minimal_depth_var.select"
  
  if (length(selected) < 2) {
    
    selector_method <- "RSF_VIMP_fallback"
    
    imp_obj <- randomForestSRC::vimp.rfsrc(fit)$importance
    
    imp <- as.numeric(imp_obj)
    names(imp) <- names(imp_obj)
    
    if (is.null(names(imp)) || any(is.na(names(imp)))) {
      names(imp) <- genes[seq_along(imp)]
    }
    
    imp_df <- data.frame(
      SafeGene = names(imp),
      VIMP = imp,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::filter(
        SafeGene %in% genes,
        is.finite(VIMP)
      ) %>%
      dplyr::mutate(
        Gene = safe_gene_back(SafeGene),
        VIMP_Positive = ifelse(VIMP > 0, VIMP, 0)
      )
    
    max_pos <- max(imp_df$VIMP_Positive, na.rm = TRUE)
    
    if (!is.finite(max_pos) || max_pos <= 0) {
      stop("RSF minimal depth and VIMP both failed to select genes.")
    }
    
    imp_df <- imp_df %>%
      dplyr::mutate(
        RSF_RelativeImportance = VIMP_Positive / max_pos,
        RSF_Selected = RSF_RelativeImportance > RSF_IMP_CUTOFF
      ) %>%
      dplyr::arrange(dplyr::desc(RSF_RelativeImportance))
    
    selected <- imp_df %>%
      dplyr::filter(RSF_Selected) %>%
      dplyr::pull(SafeGene) %>%
      unique()
  }
  
  if (length(selected) < 2) {
    stop("RSF selected fewer than 2 genes.")
  }
  
  out_df <- data.frame(
    SafeGene = selected,
    Gene = safe_gene_back(selected),
    Selector = selector_method,
    stringsAsFactors = FALSE
  )
  
  save_csv(
    out_df,
    file.path(FIG4_DIR, "Fig4A_RSF_selected_genes.csv")
  )
  
  cat("RSF selector selected genes:", length(selected), "|", selector_method, "\n")
  cat(paste(safe_gene_back(selected), collapse = " + "), "\n")
  
  selected
}

select_by_lasso <- function(genes) {
  
  est <- make_model_df(train_ml, genes)
  x <- as.matrix(est[, genes, drop = FALSE])
  y <- survival::Surv(est$OS.time, est$OS)
  
  set.seed(SEED)
  
  fit <- glmnet::cv.glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = 1,
    nfolds = CV_FOLDS
  )
  
  cf <- glmnet::coef.glmnet(fit, s = fit$lambda.min)
  selected <- rownames(cf)[as.numeric(cf) != 0]
  selected <- unique(selected)
  
  if (length(selected) < 2) {
    stop("Lasso selected fewer than 2 genes.")
  }
  
  cat("Lasso selector selected genes:", length(selected), "\n")
  cat(paste(safe_gene_back(selected), collapse = " + "), "\n")
  
  selected
}

select_by_coxboost <- function(genes) {
  
  est <- make_model_df(train_ml, genes)
  x <- as.matrix(est[, genes, drop = FALSE])
  
  pen <- CoxBoost::optimCoxBoostPenalty(
    time = est$OS.time,
    status = est$OS,
    x = x,
    trace = FALSE,
    start.penalty = 500,
    parallel = FALSE
  )
  
  cv.res <- CoxBoost::cv.CoxBoost(
    time = est$OS.time,
    status = est$OS,
    x = x,
    maxstepno = 500,
    K = CV_FOLDS,
    type = "verweij",
    penalty = pen$penalty
  )
  
  fit <- CoxBoost::CoxBoost(
    time = est$OS.time,
    status = est$OS,
    x = x,
    stepno = cv.res$optimal.step,
    penalty = pen$penalty
  )
  
  cf <- coef(fit)
  
  if (is.null(names(cf))) {
    names(cf) <- genes[seq_along(cf)]
  }
  
  selected <- names(cf)[as.numeric(cf) != 0]
  selected <- unique(selected)
  
  if (length(selected) < 2) {
    stop("CoxBoost selected fewer than 2 genes.")
  }
  
  cat("CoxBoost selector selected genes:", length(selected), "\n")
  cat(paste(safe_gene_back(selected), collapse = " + "), "\n")
  
  selected
}

select_by_stepcox <- function(genes, direction) {
  
  genes_use <- get_stepcox_genes(genes)
  
  if (length(genes_use) < 2) {
    stop("Too few genes for StepCox.")
  }
  
  est <- make_model_df(train_ml, genes_use)
  step_k <- get_stepcox_k(est)
  
  if (direction == "forward") {
    
    null_fit <- survival::coxph(
      Surv(OS.time, OS) ~ 1,
      data = est
    )
    
    upper_formula <- as.formula(
      paste0("~ ", paste(genes_use, collapse = " + "))
    )
    
    fit <- step(
      null_fit,
      scope = list(lower = ~1, upper = upper_formula),
      direction = "forward",
      trace = 0,
      k = step_k
    )
    
  } else if (direction == "backward") {
    
    full_fit <- survival::coxph(
      Surv(OS.time, OS) ~ .,
      data = est
    )
    
    fit <- step(
      full_fit,
      direction = "backward",
      trace = 0,
      k = step_k
    )
    
  } else if (direction == "both") {
    
    null_fit <- survival::coxph(
      Surv(OS.time, OS) ~ 1,
      data = est
    )
    
    upper_formula <- as.formula(
      paste0("~ ", paste(genes_use, collapse = " + "))
    )
    
    fit <- step(
      null_fit,
      scope = list(lower = ~1, upper = upper_formula),
      direction = "both",
      trace = 0,
      k = step_k
    )
    
  } else {
    stop("Unknown StepCox direction: ", direction)
  }
  
  coef_vec <- stats::coef(fit)
  selected <- names(coef_vec)[is.finite(coef_vec)]
  selected <- unique(selected)
  
  if (length(selected) < 2) {
    stop(paste0("StepCox[", direction, "] selected fewer than 2 genes."))
  }
  
  cat("StepCox[", direction, "] selector selected genes:", length(selected),
      " | criterion: ", STEPCOX_K_TYPE, "\n", sep = "")
  cat(paste(safe_gene_back(selected), collapse = " + "), "\n")
  
  selected
}

############################################################
# 14. Universal learner wrapper
############################################################

run_learner <- function(learner_name, genes) {
  
  if (length(genes) < 2) return(NULL)
  
  out <- tryCatch({
    
    if (learner_name == "RSF") {
      
      fit_predict_rsf(genes)
      
    } else if (learner_name == "Ridge") {
      
      fit_predict_glmnet(genes, 0)
      
    } else if (learner_name == "Lasso") {
      
      fit_predict_glmnet(genes, 1)
      
    } else if (grepl("^Enet", learner_name)) {
      
      alpha_value <- as.numeric(gsub("Enet\\[α=|\\]", "", learner_name))
      fit_predict_glmnet(genes, alpha_value)
      
    } else if (learner_name == "CoxBoost") {
      
      fit_predict_coxboost(genes)
      
    } else if (learner_name == "plsRcox") {
      
      fit_predict_plsrcox(genes)
      
    } else if (learner_name == "SuperPC") {
      
      fit_predict_superpc(genes)
      
    } else if (learner_name == "GBM") {
      
      fit_predict_gbm(genes)
      
    } else if (learner_name == "survival-SVM") {
      
      fit_predict_svm(genes)
      
    } else if (grepl("^StepCox", learner_name)) {
      
      direction <- gsub("StepCox\\[|\\]", "", learner_name)
      fit_predict_stepcox(genes, direction)
      
    } else {
      
      NULL
    }
    
  }, error = function(e) {
    message("Skip model: ", learner_name, " | ", e$message)
    NULL
  })
  
  out
}

############################################################
# 15. Run base models
############################################################

base_learners <- c(
  "RSF",
  "CoxBoost",
  paste0("Enet[α=", seq(0.1, 0.9, 0.1), "]"),
  "Ridge",
  "Lasso",
  "plsRcox",
  "SuperPC",
  "GBM",
  "survival-SVM",
  "StepCox[both]",
  "StepCox[backward]",
  "StepCox[forward]"
)

for (learner in base_learners) {
  
  cat("Running base model:", learner, "\n")
  
  rs_list <- run_learner(learner, safe_genes)
  
  if (!is.null(rs_list)) {
    
    used <- attr(rs_list, "selected_genes")
    if (is.null(used)) used <- safe_genes
    
    add_result(learner, rs_list, used)
  }
}

############################################################
# 16. Run combination models
############################################################

selector_list <- list()

cat("Selecting genes by RSF...\n")
selector_list[["RSF"]] <- tryCatch(select_by_rsf(safe_genes), error = function(e) {
  message("RSF selector failed: ", e$message)
  NULL
})

cat("Selecting genes by Lasso...\n")
selector_list[["Lasso"]] <- tryCatch(select_by_lasso(safe_genes), error = function(e) {
  message("Lasso selector failed: ", e$message)
  NULL
})

cat("Selecting genes by CoxBoost...\n")
selector_list[["CoxBoost"]] <- tryCatch(select_by_coxboost(safe_genes), error = function(e) {
  message("CoxBoost selector failed: ", e$message)
  NULL
})

for (direction in c("both", "backward", "forward")) {
  
  nm <- paste0("StepCox[", direction, "]")
  
  cat("Selecting genes by ", nm, "...\n")
  
  selector_list[[nm]] <- tryCatch(
    select_by_stepcox(safe_genes, direction),
    error = function(e) {
      message(nm, " selector failed: ", e$message)
      NULL
    }
  )
}

selector_count_df <- data.frame(
  Selector = names(selector_list),
  GeneN = sapply(selector_list, function(x) if (is.null(x)) NA_integer_ else length(unique(x))),
  stringsAsFactors = FALSE
)

save_csv(
  selector_count_df,
  file.path(FIG4_DIR, "ManualML_selector_gene_counts.csv")
)

for (nm in names(selector_list)) {
  if (!is.null(selector_list[[nm]])) {
    selector_store[[nm]] <- selector_list[[nm]]
  }
}

combo_learners <- c(
  "RSF",
  "CoxBoost",
  paste0("Enet[α=", seq(0.1, 0.9, 0.1), "]"),
  "Ridge",
  "Lasso",
  "plsRcox",
  "SuperPC",
  "GBM",
  "survival-SVM",
  "StepCox[both]",
  "StepCox[backward]",
  "StepCox[forward]"
)

for (selector_name in names(selector_list)) {
  
  genes_sel <- selector_list[[selector_name]]
  
  if (is.null(genes_sel) || length(genes_sel) < 2) next
  
  for (learner in combo_learners) {
    
    if (selector_name == learner) next
    
    model_name <- paste0(selector_name, " + ", learner)
    
    cat("Running combination model:", model_name, "| genes:", length(genes_sel), "\n")
    
    rs_list <- run_learner(learner, genes_sel)
    
    if (!is.null(rs_list)) {
      
      used <- attr(rs_list, "selected_genes")
      if (is.null(used)) used <- genes_sel
      
      selector_store[[model_name]] <- genes_sel
      
      add_result(model_name, rs_list, used)
    }
  }
}

############################################################
# 17. C-index summary and final model selection
############################################################

save_csv(
  result_long,
  file.path(FIG4_DIR, "ManualML_all_model_Cindex_long.csv")
)

cindex_summary <- result_long %>%
  dplyr::mutate(Cindex = as.numeric(Cindex)) %>%
  tidyr::pivot_wider(
    id_cols = Model,
    names_from = Dataset,
    values_from = Cindex
  )

dataset_cols <- intersect(
  c("Training", "Testing", "GSE72094", "GSE68465"),
  colnames(cindex_summary)
)

cindex_summary <- cindex_summary %>%
  dplyr::mutate(
    dplyr::across(dplyr::all_of(dataset_cols), as.numeric),
    Average = rowMeans(
      dplyr::across(dplyr::all_of(dataset_cols)),
      na.rm = TRUE
    ),
    ExternalMean = rowMeans(
      dplyr::across(
        dplyr::all_of(intersect(c("GSE72094", "GSE68465"), dataset_cols))
      ),
      na.rm = TRUE
    ),
    Train_Test_Gap = if (all(c("Training", "Testing") %in% dataset_cols)) {
      abs(Training - Testing)
    } else {
      NA_real_
    },
    Complexity = model_complexity(Model)
  ) %>%
  dplyr::arrange(
    dplyr::desc(Average),
    dplyr::desc(Testing),
    dplyr::desc(Training),
    Train_Test_Gap,
    Complexity
  )

feature_count_df <- data.frame(
  Model = names(feature_store),
  FeatureN = sapply(feature_store, function(x) length(unique(x))),
  stringsAsFactors = FALSE
)

save_csv(
  feature_count_df,
  file.path(FIG4_DIR, "ManualML_model_feature_counts.csv")
)

cindex_summary <- cindex_summary %>%
  dplyr::left_join(feature_count_df, by = "Model")

save_csv(
  cindex_summary,
  file.path(FIG4_DIR, "Fig4B_all_model_Cindex_summary.csv")
)

final_model <- cindex_summary$Model[1]

selected_model_info <- cindex_summary %>%
  dplyr::filter(Model == final_model)

save_csv(
  selected_model_info,
  file.path(FIG4_DIR, "Fig4B_final_selected_model_info.csv")
)

cat("\nFinal DPRS model selected by average C-index:", final_model, "\n")
print(selected_model_info)

############################################################
# 17.1 Diagnostic: combo vs base learner equality
############################################################

check_model_equal <- function(m1, m2) {
  
  if (!m1 %in% names(risk_store)) return(NULL)
  if (!m2 %in% names(risk_store)) return(NULL)
  
  x1 <- risk_store[[m1]] %>%
    dplyr::arrange(Dataset, ID) %>%
    dplyr::select(Dataset, ID, DPRS)
  
  x2 <- risk_store[[m2]] %>%
    dplyr::arrange(Dataset, ID) %>%
    dplyr::select(Dataset, ID, DPRS)
  
  merged <- dplyr::inner_join(
    x1,
    x2,
    by = c("Dataset", "ID"),
    suffix = c("_combo", "_base")
  )
  
  if (nrow(merged) == 0) return(NULL)
  
  data.frame(
    ComboModel = m1,
    BaseLearner = m2,
    N = nrow(merged),
    Cor = suppressWarnings(cor(merged$DPRS_combo, merged$DPRS_base, use = "complete.obs")),
    MaxAbsDiff = max(abs(merged$DPRS_combo - merged$DPRS_base), na.rm = TRUE),
    IdenticalRisk = isTRUE(max(abs(merged$DPRS_combo - merged$DPRS_base), na.rm = TRUE) == 0),
    stringsAsFactors = FALSE
  )
}

combo_models <- cindex_summary$Model[grepl("\\+", cindex_summary$Model)]

risk_equal_table <- dplyr::bind_rows(lapply(combo_models, function(m) {
  
  parts <- trimws(strsplit(m, "\\+")[[1]])
  base_learner <- parts[length(parts)]
  
  check_model_equal(m, base_learner)
}))

if (nrow(risk_equal_table) > 0) {
  risk_equal_table <- risk_equal_table %>%
    dplyr::arrange(dplyr::desc(IdenticalRisk), dplyr::desc(Cor))
  
  save_csv(
    risk_equal_table,
    file.path(FIG4_DIR, "ManualML_combo_vs_base_risk_equal_check.csv")
  )
  
  cat("\nIdentical combo-base risk models:\n")
  print(risk_equal_table %>% dplyr::filter(IdenticalRisk))
}

############################################################
# 18. C-index heatmap
############################################################

ci_show <- cindex_summary %>%
  dplyr::slice_head(n = 35)

ci_long <- ci_show %>%
  dplyr::select(Model, dplyr::all_of(dataset_cols), Average) %>%
  tidyr::pivot_longer(
    cols = -Model,
    names_to = "Dataset",
    values_to = "Cindex"
  )

ci_long$Model <- factor(ci_long$Model, levels = rev(ci_show$Model))
ci_long$Dataset <- factor(
  ci_long$Dataset,
  levels = c(dataset_cols, "Average")
)

p_ci <- ggplot(ci_long, aes(Dataset, Model, fill = Cindex)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", Cindex)), size = 2.4) +
  scale_fill_gradient2(
    low = "#D7E4EF",
    mid = "white",
    high = "#B2182B",
    midpoint = 0.65,
    name = "C-index"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 6, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(title = "Manual ML model performance ranked by average C-index")

safe_ggsave(
  file.path(FIG4_DIR, "Fig4B_Cindex_heatmap.pdf"),
  p_ci,
  7.5,
  8.5
)

############################################################
# 19. Extract final DPRS risk scores
############################################################

final_score <- risk_store[[final_model]]

final_score <- final_score %>%
  dplyr::filter(
    is.finite(time),
    time > 0,
    status %in% c(0, 1),
    is.finite(DPRS)
  )

train_tmp <- final_score %>%
  dplyr::filter(Dataset == "Training")

cox_dir <- survival::coxph(
  survival::Surv(time, status) ~ DPRS,
  data = train_tmp
)

if (as.numeric(coef(cox_dir)[1]) < 0) {
  final_score$DPRS <- -final_score$DPRS
}

final_score <- final_score %>%
  dplyr::group_by(Dataset) %>%
  dplyr::mutate(
    Cutoff = median(DPRS, na.rm = TRUE),
    RiskGroup = ifelse(DPRS >= Cutoff, "High", "Low")
  ) %>%
  dplyr::ungroup()

final_score$RiskGroup <- factor(final_score$RiskGroup, levels = c("Low", "High"))
final_score$Dataset <- factor(
  final_score$Dataset,
  levels = c("Training", "Testing", "GSE72094", "GSE68465")
)

save_csv(
  final_score,
  file.path(FIG4_DIR, "Fig4C_DPRS_all_sets.csv")
)

save_csv(
  final_score,
  file.path(DB_DIR, "DPRS_score_table.csv")
)

############################################################
# 20. Final signature genes used by selected model
############################################################

final_safe_genes <- unique(feature_store[[final_model]])
final_safe_genes <- final_safe_genes[!is.na(final_safe_genes)]

final_genes <- safe_gene_back(final_safe_genes)

signature_genes <- data.frame(
  Gene = final_genes,
  SafeGene = final_safe_genes,
  Model = final_model,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(!is.na(Gene), Gene != "") %>%
  dplyr::distinct(Gene, .keep_all = TRUE) %>%
  dplyr::left_join(
    cox_anno %>%
      dplyr::select(Gene, CoxCoef, HR, CoxP, CoxFDR),
    by = "Gene"
  ) %>%
  dplyr::mutate(
    GeneRole = dplyr::case_when(
      is.na(CoxCoef) ~ "Unknown",
      CoxCoef >= 0 ~ "Risk",
      CoxCoef < 0 ~ "Protective"
    ),
    SignatureType = paste0(
      "Final features actually used by selected model: ", final_model
    )
  ) %>%
  dplyr::arrange(CoxP)

save_csv(
  signature_genes,
  file.path(FIG4_DIR, "Fig4C_final_signature_genes.csv")
)

save_csv(
  signature_genes,
  file.path(DB_DIR, "DPRS_signature_gene_table.csv")
)

cat("\nFinal model:", final_model, "\n")
cat("Final signature genes actually used by model:", nrow(signature_genes), "\n")
cat(paste(signature_genes$Gene, collapse = " + "), "\n")

plot_genes <- signature_genes %>%
  dplyr::filter(is.finite(CoxCoef)) %>%
  dplyr::mutate(
    PlotCoef = CoxCoef,
    AbsCoef = abs(CoxCoef)
  ) %>%
  dplyr::arrange(dplyr::desc(AbsCoef)) %>%
  dplyr::slice_head(n = 50)

if (nrow(plot_genes) > 0) {
  
  p_gene <- ggplot(
    plot_genes,
    aes(
      x = reorder(Gene, PlotCoef),
      y = PlotCoef,
      color = GeneRole
    )
  ) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_segment(
      aes(xend = Gene, y = 0, yend = PlotCoef),
      color = "grey60",
      linewidth = 0.85
    ) +
    geom_point(size = 2.8) +
    scale_color_manual(
      values = c(
        "Risk" = "#C53030",
        "Protective" = "#2B6CB0",
        "Unknown" = "grey60"
      )
    ) +
    coord_flip() +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      x = NULL,
      y = "Training-only univariate Cox coefficient",
      title = paste0("Final signature genes: ", final_model)
    )
  
  safe_ggsave(
    file.path(FIG4_DIR, "Fig4C_candidate_gene_lollipop_top50.pdf"),
    p_gene,
    6.5,
    max(4.5, 0.25 * nrow(plot_genes) + 2)
  )
}

############################################################
# 21. Risk score and survival status plots
############################################################

risk_rank <- final_score %>%
  dplyr::filter(is.finite(DPRS)) %>%
  dplyr::arrange(Dataset, DPRS) %>%
  dplyr::group_by(Dataset) %>%
  dplyr::mutate(Rank = dplyr::row_number()) %>%
  dplyr::ungroup()

p_risk <- ggplot(risk_rank, aes(Rank, DPRS, color = RiskGroup)) +
  geom_point(size = 0.8, alpha = 0.9) +
  facet_wrap(~Dataset, scales = "free_x", nrow = 1) +
  scale_color_manual(values = c(Low = "#2B6CB0", High = "#C53030")) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = "Patient rank",
    y = "DPRS",
    title = paste0("Risk score distribution: ", final_model)
  )

safe_ggsave(
  file.path(FIG4_DIR, "Fig4D_risk_score_distribution.pdf"),
  p_risk,
  12,
  4
)

p_status <- ggplot(risk_rank, aes(Rank, time, color = factor(status))) +
  geom_point(size = 0.8, alpha = 0.9) +
  facet_wrap(~Dataset, scales = "free_x", nrow = 1) +
  scale_color_manual(
    values = c("0" = "#2B6CB0", "1" = "#C53030"),
    labels = c("Alive/Censored", "Dead")
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = "Patient rank",
    y = "Survival time",
    color = "Status",
    title = "Survival status distribution"
  )

safe_ggsave(
  file.path(FIG4_DIR, "Fig4D_survival_status_distribution.pdf"),
  p_status,
  12,
  4
)

score_cindex <- final_score %>%
  dplyr::group_by(Dataset) %>%
  dplyr::summarise(
    N = dplyr::n(),
    Events = sum(status == 1),
    Cindex = calc_cindex_risk(time, status, DPRS),
    .groups = "drop"
  )

save_csv(
  score_cindex,
  file.path(FIG4_DIR, "Fig4D_selected_model_Cindex_by_dataset.csv")
)

print(score_cindex)

############################################################
# 22. KM and timeROC
############################################################

for (ds in unique(as.character(final_score$Dataset))) {
  
  df <- final_score %>%
    dplyr::filter(Dataset == ds)
  
  fit <- survival::survfit(
    survival::Surv(time, status) ~ RiskGroup,
    data = df
  )
  
  p_km <- survminer::ggsurvplot(
    fit,
    data = df,
    pval = TRUE,
    risk.table = FALSE,
    palette = c("#2B6CB0", "#C53030"),
    title = paste0(ds, " DPRS"),
    legend.title = "",
    legend.labs = c("Low", "High"),
    xlab = "Time",
    ylab = "Overall survival probability",
    ggtheme = theme_bw()
  )$plot
  
  safe_ggsave(
    file.path(FIG4_DIR, paste0("Fig4E_", ds, "_KM.pdf")),
    p_km,
    6,
    5
  )
  
  times_use <- get_time_points(df$time)
  
  roc <- timeROC::timeROC(
    T = df$time,
    delta = df$status,
    marker = df$DPRS,
    cause = 1,
    weighting = "marginal",
    times = times_use,
    ROC = TRUE
  )
  
  auc_df <- data.frame(
    Dataset = ds,
    Model = final_model,
    Time = c("1-year", "3-year", "5-year"),
    TimeValue = times_use,
    AUC = as.numeric(roc$AUC)
  )
  
  save_csv(
    auc_df,
    file.path(FIG4_DIR, paste0("Fig4F_", ds, "_timeROC_AUC.csv"))
  )
  
  pdf(
    file.path(FIG4_DIR, paste0("Fig4F_", ds, "_timeROC.pdf")),
    width = 6,
    height = 5,
    useDingbats = FALSE
  )
  
  plot(roc, time = times_use[1], col = "#C53030", title = FALSE)
  plot(roc, time = times_use[2], add = TRUE, col = "#2B6CB0")
  plot(roc, time = times_use[3], add = TRUE, col = "#2F855A")
  
  legend(
    "bottomright",
    legend = paste0(auc_df$Time, " AUC=", sprintf("%.3f", auc_df$AUC)),
    col = c("#C53030", "#2B6CB0", "#2F855A"),
    lwd = 2,
    bty = "n"
  )
  
  title(main = paste0(ds, " DPRS time-dependent ROC"))
  dev.off()
}

############################################################
# 23. Done
############################################################

cat("\nManual DPRS machine-learning pipeline finished.\n")
cat("Final model:", final_model, "\n")
cat("Candidate genes:", length(candidate_genes), "\n")
cat("Final feature genes:", length(final_genes), "\n")
cat("Output directory:", FIG4_DIR, "\n")

cat("\nMain outputs:\n")
cat(" - Fig4A_TCGA_score_survival_table.csv\n")
cat(" - Fig4A_training_only_DEG.csv\n")
cat(" - Fig4A_candidate_genes_DEG_unicox_training_only.csv\n")
cat(" - Fig4A_training_only_unicox_for_candidate_selection.csv\n")
cat(" - Fig4A_RSF_selected_genes.csv\n")
cat(" - ManualML_selector_gene_counts.csv\n")
cat(" - ManualML_all_model_Cindex_long.csv\n")
cat(" - Fig4B_all_model_Cindex_summary.csv\n")
cat(" - Fig4B_final_selected_model_info.csv\n")
cat(" - ManualML_combo_vs_base_risk_equal_check.csv\n")
cat(" - Fig4B_Cindex_heatmap.pdf\n")
cat(" - Fig4C_DPRS_all_sets.csv\n")
cat(" - Fig4C_final_signature_genes.csv\n")
cat(" - Fig4D_selected_model_Cindex_by_dataset.csv\n")
cat(" - Fig4E_*_KM.pdf\n")
cat(" - Fig4F_*_timeROC.pdf\n")

