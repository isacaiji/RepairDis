############################################################
# Published DDR-related signatures comparison
#
# Purpose:
#   Compare ML-DDRscore with published DDR-related prognostic signatures.
#
# Required before running:
#   1) Run your MO-DDRscore / MLsig pipeline first.
#   2) Current R session should contain:
#        tcga_expr: gene x sample expression matrix
#        tcga_clin: Patient, time, status
#   3) Put published signature table here:
#        D:/R_workspace/评分/AD_DDR_project/00_data/Published_signatures.csv
#
# Published_signatures.csv columns:
#   Signature,Gene,Coefficient,Source,Type,Year,Endpoint,Note,PMID
#
# Output:
#   D:/R_workspace/评分/AD_DDR_project/MLsig/Published_signature_benchmark
############################################################

options(stringsAsFactors = FALSE)
set.seed(20260513)

############################
# 0. Paths
############################

PROJECT_DIR <- "D:/R_workspace/评分/AD_DDR_project"
DATA_DIR <- file.path(PROJECT_DIR, "00_data")
MLSIG_DIR <- file.path(PROJECT_DIR, "MLsig")
OUT_DIR <- file.path(MLSIG_DIR, "Published_signature_benchmark")
FIG_DIR <- file.path(OUT_DIR, "fig")
TAB_DIR <- file.path(OUT_DIR, "table")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

SIGNATURE_FILE <- file.path(DATA_DIR, "Published_signatures.csv")

############################
# 1. Packages and helpers
############################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(survival)
})

save_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file)
}

save_pdf <- function(file, p, w = 7, h = 5) {
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

sample_type_code <- function(x) {
  x <- gsub("\\.", "-", as.character(x))
  ifelse(nchar(x) >= 15, substr(x, 14, 15), NA_character_)
}

zscore_vector <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  sdx <- sd(x, na.rm = TRUE)
  if (all(is.na(x)) || !is.finite(sdx) || sdx == 0) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

calc_cindex <- function(time, status, score) {
  df <- data.frame(
    time = as.numeric(time),
    status = as.numeric(status),
    score = as.numeric(score)
  )
  
  df <- df[
    is.finite(df$time) &
      is.finite(df$status) &
      is.finite(df$score),
  ]
  
  if (nrow(df) < 30 || sum(df$status == 1) < 5) {
    return(NA_real_)
  }
  
  # score 越高 = 风险越高 = 生存越差
  cidx <- survival::concordance(
    survival::Surv(time, status) ~ score,
    data = df,
    reverse = TRUE
  )$concordance
  
  as.numeric(cidx)
}

get_time_points <- function(time) {
  # TCGA/GEO 通常是 days；如果是 years，则用 1/3/5
  if (median(time, na.rm = TRUE) > 100) {
    c(365, 1095, 1825)
  } else {
    c(1, 3, 5)
  }
}

make_display_name <- function(sig, pmid_vec = NULL) {
  if (sig == "ML_DDRscore") {
    return("DPRS")
  }
  if (!is.null(pmid_vec) && sig %in% names(pmid_vec)) {
    pmid <- pmid_vec[[sig]]
    if (!is.na(pmid) && pmid != "") return(as.character(pmid))
  }
  sig
}

# 少量基因别名。主要防止 HIST3H2A 在表达矩阵里写成 H3C15。
gene_alias_map <- c(
  "HIST3H2A" = "H3C15"
)

############################
# 2. Required objects
############################

stopifnot(exists("tcga_expr"))
stopifnot(exists("tcga_clin"))

tcga_expr <- as.matrix(tcga_expr)
rownames(tcga_expr) <- clean_gene(rownames(tcga_expr))
tcga_expr <- tcga_expr[!is.na(rownames(tcga_expr)) & rownames(tcga_expr) != "", , drop = FALSE]
tcga_expr <- tcga_expr[!duplicated(rownames(tcga_expr)), , drop = FALSE]

tcga_clin$Patient <- as.character(tcga_clin$Patient)

cat("TCGA expr genes:", nrow(tcga_expr), "\n")
cat("TCGA expr samples:", ncol(tcga_expr), "\n")
cat("TCGA clinical patients:", nrow(tcga_clin), "\n")

############################
# 3. Read published signatures
############################

sig_tbl <- data.table::fread(
  SIGNATURE_FILE,
  data.table = FALSE,
  check.names = FALSE,
  fill = TRUE
)

colnames(sig_tbl) <- trimws(colnames(sig_tbl))

if (!"PMID" %in% colnames(sig_tbl)) {
  sig_tbl$PMID <- NA_character_
}

required_sig_cols <- c(
  "Signature", "Gene", "Coefficient",
  "Source", "Type", "Year", "Endpoint", "Note", "PMID"
)

missing_sig_cols <- setdiff(required_sig_cols, colnames(sig_tbl))
if (length(missing_sig_cols) > 0) {
  stop("Published signature table missing columns: ",
       paste(missing_sig_cols, collapse = ", "))
}

sig_tbl <- sig_tbl %>%
  dplyr::select(all_of(required_sig_cols)) %>%
  mutate(across(everything(), ~ trimws(as.character(.x)))) %>%
  filter(
    !(is.na(Signature) | Signature == "") |
      !(is.na(Gene) | Gene == "")
  ) %>%
  mutate(
    Signature = trimws(as.character(Signature)),
    Gene = clean_gene(Gene),
    Coefficient = suppressWarnings(as.numeric(Coefficient)),
    Source = trimws(as.character(Source)),
    Type = trimws(as.character(Type)),
    Year = suppressWarnings(as.integer(Year)),
    Endpoint = trimws(as.character(Endpoint)),
    Note = trimws(as.character(Note)),
    PMID = gsub("[^0-9]", "", as.character(PMID)),
    PMID = ifelse(PMID == "", NA, PMID)
  ) %>%
  filter(
    !is.na(Signature),
    Signature != "",
    !is.na(Gene),
    Gene != "",
    is.finite(Coefficient)
  ) %>%
  distinct(Signature, Gene, .keep_all = TRUE)

# 只保留目前决定进入主比较的 6 个模型
keep_sigs <- c(
  "Li_2022_NSCLC_DDR6",
  "Zhao_2022_EarlyLUAD_DRG16",
  "Chen_2021_LUAD_DRG6",
  "Li_2022_RespRes_DCG4",
  "Hu_2020_LUAD_DNARepair13",
  "Yang_2020_SurgLUAD_DNARepair6"
)

sig_tbl <- sig_tbl %>%
  filter(Signature %in% keep_sigs)

sig_meta <- sig_tbl %>%
  group_by(Signature) %>%
  summarise(
    Source = first(Source),
    Type = first(Type),
    Year = first(Year),
    Endpoint = first(Endpoint),
    Note = first(Note),
    PMID = first(PMID),
    SignatureGeneN = n_distinct(Gene),
    Genes = paste(unique(Gene), collapse = ";"),
    .groups = "drop"
  )

save_csv(sig_tbl, file.path(TAB_DIR, "Published_signatures_cleaned.csv"))
save_csv(sig_meta, file.path(TAB_DIR, "Published_signatures_metadata.csv"))

cat("Published signatures loaded:\n")
print(table(sig_tbl$Signature))

cat("\nPublished signature metadata:\n")
print(sig_meta)

############################
# 4. Read train/test split and best MLsig run
############################

FIG4_SIMPLE_DIR <- file.path(PROJECT_DIR, "04")

final_model_info_file <- file.path(FIG4_SIMPLE_DIR, "Fig4B_final_selected_model_info.csv")
if (!file.exists(final_model_info_file)) {
  stop("Cannot find: ", final_model_info_file)
}

global_best <- data.table::fread(
  final_model_info_file,
  data.table = FALSE,
  check.names = FALSE
)

best_model <- as.character(global_best$Model[1])
best_dir <- FIG4_SIMPLE_DIR

cat("\nBest ML model from Fig4 module:", best_model, "\n")

train_file <- file.path(best_dir, "Mime1_input_Training.csv")
test_file  <- file.path(best_dir, "Mime1_input_Testing.csv")

if (!file.exists(train_file)) stop("Cannot find: ", train_file)
if (!file.exists(test_file)) stop("Cannot find: ", test_file)

train_input <- data.table::fread(
  train_file,
  data.table = FALSE,
  check.names = FALSE
)

test_input <- data.table::fread(
  test_file,
  data.table = FALSE,
  check.names = FALSE
)

train_patients <- as.character(train_input$ID)
test_patients <- as.character(test_input$ID)

cat("Training patients:", length(train_patients), "\n")
cat("Testing patients:", length(test_patients), "\n")

############################
# 5. Build TCGA expression object by patient
############################

make_tcga_expr_by_patient <- function(patient_ids) {
  all_samples <- colnames(tcga_expr)
  all_patients <- patient_id(all_samples)
  
  use_samples <- sapply(patient_ids, function(pid) {
    smp <- all_samples[all_patients == pid]
    if (length(smp) == 0) return(NA_character_)
    
    # 优先选 tumor sample 01
    st <- sample_type_code(smp)
    smp_tumor <- smp[st == "01"]
    if (length(smp_tumor) > 0) return(smp_tumor[1])
    
    # 如果列名本身就是 patient ID 或没有 01 信息，则直接选第一个
    smp[1]
  })
  
  valid <- !is.na(use_samples)
  patient_ids <- patient_ids[valid]
  use_samples <- use_samples[valid]
  
  expr <- tcga_expr[, use_samples, drop = FALSE]
  colnames(expr) <- patient_ids
  
  # 如果你的 tcga_expr 已经是 log2(TPM+1)，这里可能会再次 log。
  # 为了兼容 TPM 原始矩阵，这里按数值范围判断。
  expr_num <- expr
  storage.mode(expr_num) <- "numeric"
  if (quantile(expr_num, 0.95, na.rm = TRUE) > 50) {
    expr_log <- log2(expr_num + 1)
  } else {
    expr_log <- expr_num
  }
  
  clin <- tcga_clin %>%
    filter(Patient %in% patient_ids) %>%
    distinct(Patient, .keep_all = TRUE)
  
  list(expr = expr_log, clin = clin)
}

tcga_train <- make_tcga_expr_by_patient(train_patients)
tcga_test  <- make_tcga_expr_by_patient(test_patients)

cat("TCGA Training expr:", nrow(tcga_train$expr), "genes x", ncol(tcga_train$expr), "samples\n")
cat("TCGA Testing expr:", nrow(tcga_test$expr), "genes x", ncol(tcga_test$expr), "samples\n")

############################
# 6. Build GEO expression object
############################

read_geo_expr_clin <- function(expr_file, clin_file) {
  if (!file.exists(expr_file)) stop("Cannot find: ", expr_file)
  if (!file.exists(clin_file)) stop("Cannot find: ", clin_file)
  
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
  
  if (quantile(expr_mat, 0.95, na.rm = TRUE) > 50) {
    expr_mat <- log2(expr_mat + 1)
  }
  
  clin$Sample <- as.character(clin$Sample)
  
  list(expr = expr_mat, clin = clin)
}

gse72094 <- read_geo_expr_clin(
  file.path(DATA_DIR, "GSE72094_expression.csv"),
  file.path(DATA_DIR, "GSE72094_clinical.csv")
)

gse68465 <- read_geo_expr_clin(
  file.path(DATA_DIR, "GSE68465_expression.csv"),
  file.path(DATA_DIR, "GSE68465_clinical.csv")
)

cohort_list <- list(
  Training = tcga_train,
  Testing = tcga_test,
  GSE72094 = gse72094,
  GSE68465 = gse68465
)

############################
# 7. Calculate published signature scores
############################

calc_signature_scores_one_cohort <- function(expr, clin, dataset_name, sig_tbl) {
  
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  
  expr_z <- t(apply(expr, 1, zscore_vector))
  rownames(expr_z) <- rownames(expr)
  colnames(expr_z) <- colnames(expr)
  
  score_list <- list()
  overlap_list <- list()
  
  for (sig in unique(sig_tbl$Signature)) {
    ss <- sig_tbl %>% filter(Signature == sig)
    
    # 基因别名映射
    ss$Gene_used <- ss$Gene
    for (old in names(gene_alias_map)) {
      new <- gene_alias_map[[old]]
      if (old %in% ss$Gene_used &&
          !(old %in% rownames(expr_z)) &&
          new %in% rownames(expr_z)) {
        ss$Gene_used[ss$Gene_used == old] <- new
      }
    }
    
    genes <- intersect(ss$Gene_used, rownames(expr_z))
    
    # 至少匹配 2 个基因，否则跳过
    if (length(genes) < 2) next
    
    coef_vec <- ss$Coefficient[match(genes, ss$Gene_used)]
    names(coef_vec) <- genes
    
    score <- as.numeric(crossprod(coef_vec, expr_z[genes, , drop = FALSE]))
    
    score_list[[sig]] <- data.frame(
      ID = colnames(expr_z),
      Signature = sig,
      Score = score,
      Dataset = dataset_name,
      stringsAsFactors = FALSE
    )
    
    overlap_list[[sig]] <- data.frame(
      Signature = sig,
      PMID = first(ss$PMID),
      Source = first(ss$Source),
      Type = first(ss$Type),
      Year = first(ss$Year),
      Dataset = dataset_name,
      SignatureGeneN = nrow(ss),
      MatchedGeneN = length(genes),
      MatchRate = length(genes) / nrow(ss),
      OriginalGenes = paste(ss$Gene, collapse = ";"),
      UsedGenes = paste(ss$Gene_used, collapse = ";"),
      MatchedGenes = paste(genes, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  
  score_df <- bind_rows(score_list)
  overlap_df <- bind_rows(overlap_list)
  
  if (nrow(score_df) == 0) {
    return(list(score = data.frame(), overlap = overlap_df))
  }
  
  if (dataset_name %in% c("Training", "Testing")) {
    clin2 <- clin %>%
      transmute(
        ID = as.character(Patient),
        time = as.numeric(time),
        status = as.numeric(status)
      )
  } else {
    clin2 <- clin %>%
      transmute(
        ID = as.character(Sample),
        time = as.numeric(time),
        status = as.numeric(status)
      )
  }
  
  score_df <- score_df %>%
    left_join(clin2, by = "ID") %>%
    filter(is.finite(time), time > 0, status %in% c(0, 1), is.finite(Score))
  
  list(score = score_df, overlap = overlap_df)
}

all_score <- list()
all_overlap <- list()

for (ds in names(cohort_list)) {
  cat("Calculating published scores:", ds, "\n")
  res <- calc_signature_scores_one_cohort(
    expr = cohort_list[[ds]]$expr,
    clin = cohort_list[[ds]]$clin,
    dataset_name = ds,
    sig_tbl = sig_tbl
  )
  all_score[[ds]] <- res$score
  all_overlap[[ds]] <- res$overlap
}

published_score <- bind_rows(all_score)
published_overlap <- bind_rows(all_overlap)

save_csv(published_score, file.path(TAB_DIR, "Published_signature_scores_long_raw.csv"))
save_csv(published_overlap, file.path(TAB_DIR, "Published_signature_gene_overlap.csv"))

cat("\nPublished signature gene overlap:\n")
print(published_overlap %>% dplyr::select(Signature, Dataset, SignatureGeneN, MatchedGeneN, MatchRate))

############################
# 8. Add your ML-DDRscore
############################

our_score_files <- c(
  file.path(best_dir, "Fig4C_ML_DDRscore_all_sets.csv"),
  file.path(MLSIG_DIR, "Final_Fig4_reference_style", "Fig4_score_survival_groups_reference_style.csv"),
  file.path(MLSIG_DIR, "Final_Fig4_plots", "Fig4C_ML_DDRscore_all_sets_for_plot.csv"),
  file.path(MLSIG_DIR, "Final_Fig4_plots", "Fig4C_ML_DDRscore_all_sets.csv")
)

our_score_file <- our_score_files[file.exists(our_score_files)][1]

if (is.na(our_score_file)) {
  stop("Cannot find ML-DDRscore output file. Please run Fig4 plotting module first.")
}

cat("\nUsing ML-DDRscore file:\n", our_score_file, "\n")

our_score <- data.table::fread(
  our_score_file,
  data.table = FALSE,
  check.names = FALSE
)

score_col <- intersect(c("riskscore", "ML_DDRscore", "RiskScore", "risk_score"), colnames(our_score))[1]
id_col <- intersect(c("ID", "Sample", "Patient"), colnames(our_score))[1]

if (is.na(score_col)) {
  stop("Cannot identify score column in ML-DDRscore file.")
}
if (is.na(id_col)) {
  stop("Cannot identify ID column in ML-DDRscore file.")
}

our_score2 <- our_score %>%
  transmute(
    ID = as.character(.data[[id_col]]),
    time = as.numeric(time),
    status = as.numeric(status),
    Dataset = as.character(Dataset),
    Signature = "ML_DDRscore",
    Score = as.numeric(.data[[score_col]])
  ) %>%
  filter(is.finite(time), time > 0, status %in% c(0, 1), is.finite(Score))

all_scores <- bind_rows(
  our_score2,
  published_score
)

############################
# 9. Direction correction and median grouping
############################

# 方向校正：以 Training 中 Cox 系数为准。
# 如果某 signature 在 Training 中 coef < 0，则把 score 乘以 -1，
# 保证 score 越高，风险越高。
all_scores <- all_scores %>%
  group_by(Signature) %>%
  group_modify(~{
    x <- .x
    
    tr <- x %>% filter(Dataset == "Training")
    fit <- tryCatch(
      survival::coxph(Surv(time, status) ~ Score, data = tr),
      error = function(e) NULL
    )
    
    if (!is.null(fit) && length(coef(fit)) > 0 && is.finite(as.numeric(coef(fit)[1]))) {
      if (as.numeric(coef(fit)[1]) < 0) {
        x$Score <- -x$Score
      }
    }
    
    x
  }) %>%
  ungroup()

all_scores <- all_scores %>%
  group_by(Signature, Dataset) %>%
  mutate(
    Cutoff = median(Score, na.rm = TRUE),
    RiskGroup = ifelse(Score >= Cutoff, "High", "Low")
  ) %>%
  ungroup()

save_csv(all_scores, file.path(TAB_DIR, "All_signature_scores_with_ML_DDRscore.csv"))

############################
# 10. C-index, Cox HR, timeROC AUC
############################

metric_table <- all_scores %>%
  group_by(Signature, Dataset) %>%
  summarise(
    N = n(),
    Events = sum(status == 1),
    Cindex = calc_cindex(time, status, Score),
    HR = {
      fit <- tryCatch(coxph(Surv(time, status) ~ RiskGroup, data = cur_data()), error = function(e) NULL)
      if (is.null(fit)) NA_real_ else as.numeric(summary(fit)$coefficients[1, "exp(coef)"])
    },
    CoxP = {
      fit <- tryCatch(coxph(Surv(time, status) ~ RiskGroup, data = cur_data()), error = function(e) NULL)
      if (is.null(fit)) NA_real_ else as.numeric(summary(fit)$coefficients[1, "Pr(>|z|)"])
    },
    .groups = "drop"
  )

metric_table <- metric_table %>%
  left_join(sig_meta, by = "Signature") %>%
  mutate(
    Source = ifelse(Signature == "ML_DDRscore", "This_study", Source),
    Type = ifelse(Signature == "ML_DDRscore", "MO_DDRscore_derived_ML_model", Type),
    Year = ifelse(Signature == "ML_DDRscore", 2026, Year),
    Endpoint = ifelse(Signature == "ML_DDRscore", "OS", Endpoint),
    Note = ifelse(Signature == "ML_DDRscore", "this_study", Note),
    PMID = ifelse(Signature == "ML_DDRscore", NA, PMID)
  )

save_csv(metric_table, file.path(TAB_DIR, "Benchmark_Cindex_HR_by_dataset.csv"))

cindex_wide <- metric_table %>%
  dplyr::select(Signature, Dataset, Cindex) %>%
  pivot_wider(names_from = Dataset, values_from = Cindex)

dataset_cols <- intersect(c("Training", "Testing", "GSE72094", "GSE68465"), colnames(cindex_wide))

cindex_wide <- cindex_wide %>%
  mutate(
    Average = rowMeans(across(all_of(dataset_cols)), na.rm = TRUE),
    ExternalMean = rowMeans(
      across(all_of(intersect(c("GSE72094", "GSE68465"), dataset_cols))),
      na.rm = TRUE
    ),
    MinValidation = pmin(Testing, GSE72094, GSE68465, na.rm = TRUE)
  ) %>%
  left_join(sig_meta, by = "Signature") %>%
  mutate(
    Source = ifelse(Signature == "ML_DDRscore", "This_study", Source),
    Type = ifelse(Signature == "ML_DDRscore", "MO_DDRscore_derived_ML_model", Type),
    Year = ifelse(Signature == "ML_DDRscore", 2026, Year),
    Endpoint = ifelse(Signature == "ML_DDRscore", "OS", Endpoint),
    Note = ifelse(Signature == "ML_DDRscore", "this_study", Note),
    PMID = ifelse(Signature == "ML_DDRscore", NA, PMID)
  ) %>%
  arrange(desc(Average), desc(ExternalMean))

save_csv(cindex_wide, file.path(TAB_DIR, "Benchmark_Cindex_wide_summary.csv"))

cat("\nC-index summary:\n")
print(cindex_wide)

############################
# 10B. timeROC AUC and ROC curve coordinates
############################

auc_list <- list()
roc_curve_list <- list()

if (requireNamespace("timeROC", quietly = TRUE)) {
  
  for (sig in unique(all_scores$Signature)) {
    for (ds in unique(all_scores$Dataset)) {
      
      dd <- all_scores %>%
        filter(Signature == sig, Dataset == ds) %>%
        filter(is.finite(time), is.finite(status), is.finite(Score))
      
      if (nrow(dd) < 40 || sum(dd$status == 1) < 10) next
      
      times <- get_time_points(dd$time)
      time_labels <- c("1 year", "3 years", "5 years")
      
      roc <- tryCatch(
        timeROC::timeROC(
          T = dd$time,
          delta = dd$status,
          marker = dd$Score,
          cause = 1,
          weighting = "marginal",
          times = times,
          ROC = TRUE
        ),
        error = function(e) {
          message("timeROC failed: ", sig, " - ", ds, " : ", e$message)
          NULL
        }
      )
      
      if (is.null(roc)) next
      
      auc_now <- data.frame(
        Signature = sig,
        Dataset = ds,
        Time = time_labels,
        TimeValue = times,
        AUC = as.numeric(roc$AUC),
        stringsAsFactors = FALSE
      )
      
      auc_list[[paste(sig, ds, sep = "__")]] <- auc_now
      
      if (!is.null(roc$FP) && !is.null(roc$TP)) {
        for (i in seq_along(times)) {
          curve_now <- data.frame(
            Signature = sig,
            Dataset = ds,
            Time = time_labels[i],
            TimeValue = times[i],
            FP = as.numeric(roc$FP[, i]),
            TP = as.numeric(roc$TP[, i]),
            AUC = as.numeric(roc$AUC[i]),
            stringsAsFactors = FALSE
          )
          
          curve_now <- curve_now %>%
            filter(is.finite(FP), is.finite(TP)) %>%
            arrange(FP, TP)
          
          roc_curve_list[[paste(sig, ds, time_labels[i], sep = "__")]] <- curve_now
        }
      }
    }
  }
  
} else {
  message("Package timeROC is not installed. Skip time-dependent ROC.")
}

auc_table <- bind_rows(auc_list)
roc_curve_long <- bind_rows(roc_curve_list)

save_csv(auc_table, file.path(TAB_DIR, "Benchmark_timeROC_AUC_long.csv"))
save_csv(roc_curve_long, file.path(TAB_DIR, "Benchmark_timeROC_curve_long.csv"))

############################
# 11. Plots
############################

plot_df <- cindex_wide %>%
  arrange(desc(Average), desc(ExternalMean))

############################
# 11A. ROC comparison plots
############################

plot_roc_compare <- function(dataset_to_plot,
                             times_to_plot = c("1 year", "3 years"),
                             out_file = NULL) {
  
  if (!exists("roc_curve_long") || nrow(roc_curve_long) == 0) {
    message("roc_curve_long is empty. Skip ROC comparison plot.")
    return(NULL)
  }
  
  dd_curve <- roc_curve_long %>%
    dplyr::filter(
      Dataset == dataset_to_plot,
      Time %in% times_to_plot
    )
  
  if (nrow(dd_curve) == 0) {
    message("No ROC curve data for dataset: ", dataset_to_plot)
    return(NULL)
  }
  
  dd_auc <- auc_table %>%
    dplyr::filter(
      Dataset == dataset_to_plot,
      Time %in% times_to_plot
    )
  
  pmid_map <- sig_meta$PMID
  names(pmid_map) <- sig_meta$Signature
  
  dd_auc <- dd_auc %>%
    dplyr::mutate(
      DisplayName = sapply(Signature, make_display_name, pmid_vec = pmid_map),
      DisplayName = ifelse(Signature == "ML_DDRscore", "DPRS", DisplayName),
      Label = paste0(DisplayName, ": ", sprintf("%.2f", AUC))
    )
  
  dd_curve <- dd_curve %>%
    dplyr::left_join(
      dd_auc %>% dplyr::select(Signature, Time, DisplayName, Label),
      by = c("Signature", "Time")
    )
  
  sig_order <- dd_auc %>%
    dplyr::group_by(Signature) %>%
    dplyr::summarise(MeanAUC = mean(AUC, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(IsOur = Signature == "ML_DDRscore") %>%
    dplyr::arrange(desc(IsOur), desc(MeanAUC)) %>%
    dplyr::pull(Signature)
  
  dd_curve$Signature <- factor(dd_curve$Signature, levels = sig_order)
  dd_auc$Signature <- factor(dd_auc$Signature, levels = sig_order)
  
  base_cols <- c(
    "#E41A1C", "#A65628", "#4DAF4A", "#377EB8",
    "#984EA3", "#FF00AA", "#0000CC", "#999999",
    "#FF7F00", "#66C2A5", "#FC8D62"
  )
  
  use_cols <- base_cols[seq_along(sig_order)]
  names(use_cols) <- sig_order
  
  # 标签放到左上角，避免被右边面板截断
  label_df <- dd_auc %>%
    dplyr::arrange(Time, desc(Signature == "ML_DDRscore"), desc(AUC)) %>%
    dplyr::group_by(Time) %>%
    dplyr::mutate(
      x = 0.06,
      y = 0.94 - (dplyr::row_number() - 1) * 0.070
    ) %>%
    dplyr::ungroup()
  
  p <- ggplot(dd_curve, aes(x = FP, y = TP, color = Signature, group = Signature)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      linewidth = 0.8,
      color = "grey75"
    ) +
    geom_step(linewidth = 1.05, alpha = 0.95) +
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = Label, color = Signature),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0.5,
      size = 3.3,
      fontface = "bold"
    ) +
    facet_wrap(~Time, nrow = 1) +
    scale_color_manual(values = use_cols, guide = "none") +
    coord_fixed(
      ratio = 1,
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE,
      clip = "on"
    ) +
    theme_bw(base_size = 13) +
    theme(
      strip.background = element_rect(fill = "grey88", color = "grey88"),
      strip.text = element_text(face = "bold", size = 15, color = "black"),
      axis.text = element_text(color = "black", size = 11),
      axis.title = element_text(color = "black", face = "bold", size = 15),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.8),
      panel.spacing = grid::unit(0.9, "lines"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.margin = margin(8, 12, 8, 8)
    ) +
    labs(
      x = "1 - Specificity",
      y = "Sensitivity",
      title = paste0("Time-dependent ROC comparison - ", dataset_to_plot)
    )
  
  if (!is.null(out_file)) {
    # 根据 panel 数自动调大画布
    n_time <- length(unique(dd_curve$Time))
    if (n_time == 2) {
      save_pdf(out_file, p, w = 12.5, h = 5.2)
    } else if (n_time == 3) {
      save_pdf(out_file, p, w = 16.5, h = 5.3)
    } else {
      save_pdf(out_file, p, w = 6.2 * n_time, h = 5.2)
    }
  }
  
  p
}

if (nrow(roc_curve_long) > 0) {
  for (ds in c("Training", "Testing", "GSE72094", "GSE68465")) {
    plot_roc_compare(
      dataset_to_plot = ds,
      times_to_plot = c("1 year", "3 years"),
      out_file = file.path(FIG_DIR, paste0("Benchmark_ROC_compare_", ds, "_1y_3y.pdf"))
    )
  }
  
  for (ds in c("Training", "Testing", "GSE72094", "GSE68465")) {
    plot_roc_compare(
      dataset_to_plot = ds,
      times_to_plot = c("1 year", "3 years", "5 years"),
      out_file = file.path(FIG_DIR, paste0("Benchmark_ROC_compare_", ds, "_1y_3y_5y.pdf"))
    )
  }
}

############################
# 11B. C-index heatmap
############################

heat_cols <- intersect(c("Training", "Testing", "GSE72094", "GSE68465", "Average", "ExternalMean"), colnames(plot_df))

heat_long <- plot_df %>%
  dplyr::select(Signature, all_of(heat_cols)) %>%
  pivot_longer(cols = -Signature, names_to = "Dataset", values_to = "Cindex")

heat_long$Signature <- factor(heat_long$Signature, levels = rev(plot_df$Signature))
heat_long$Dataset <- factor(heat_long$Dataset, levels = heat_cols)

p_heat <- ggplot(heat_long, aes(Dataset, Signature, fill = Cindex)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.3f", Cindex)), size = 3) +
  scale_fill_gradientn(
    colours = c("#D7E4EF", "#FFFFFF", "#F4A582", "#B2182B"),
    limits = c(0.45, 0.90),
    oob = scales::squish,
    name = "C-index"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid = element_blank()
  ) +
  labs(title = "DPRS versus published DDR-related signatures")

save_pdf(file.path(FIG_DIR, "Benchmark_Cindex_heatmap.pdf"), p_heat, 8, 5.5)

############################
# 11C. Average C-index barplot
############################

bar_df <- plot_df %>%
  mutate(
    Signature = factor(Signature, levels = rev(Signature)),
    Highlight = ifelse(Signature == "ML_DDRscore", "DPRS", "Published")
  )

p_bar <- ggplot(bar_df, aes(Average, Signature, fill = Highlight)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = sprintf("%.3f", Average)),
    hjust = -0.08,
    size = 3.3,
    fontface = "bold"
  ) +
  geom_vline(xintercept = 0.65, linetype = 2, color = "grey45") +
  scale_fill_manual(values = c("DPRS" = "#C53030", "Published" = "grey75")) +
  coord_cartesian(xlim = c(0.45, max(bar_df$Average, na.rm = TRUE) + 0.06)) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold")
  ) +
  labs(
    x = "Average C-index across cohorts",
    y = NULL,
    title = "Average C-index comparison"
  )

save_pdf(file.path(FIG_DIR, "Benchmark_average_Cindex_barplot.pdf"), p_bar, 7.8, 4.8)

############################
# 11D. External C-index barplot
############################

p_ext <- ggplot(bar_df, aes(ExternalMean, Signature, fill = Highlight)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = sprintf("%.3f", ExternalMean)),
    hjust = -0.08,
    size = 3.3,
    fontface = "bold"
  ) +
  geom_vline(xintercept = 0.60, linetype = 2, color = "grey45") +
  scale_fill_manual(values = c("DPRS" = "#C53030", "Published" = "grey75")) +
  coord_cartesian(xlim = c(0.45, max(bar_df$ExternalMean, na.rm = TRUE) + 0.06)) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold")
  ) +
  labs(
    x = "External mean C-index",
    y = NULL,
    title = "External validation C-index comparison"
  )

save_pdf(file.path(FIG_DIR, "Benchmark_external_Cindex_barplot.pdf"), p_ext, 7.8, 4.8)

############################
# 11E. C-index dotplot by dataset
############################

cindex_dot <- cindex_wide %>%
  dplyr::select(Signature, all_of(dataset_cols), Average, ExternalMean) %>%
  pivot_longer(
    cols = all_of(dataset_cols),
    names_to = "Dataset",
    values_to = "Cindex"
  ) %>%
  filter(is.finite(Cindex)) %>%
  mutate(
    Signature = factor(
      Signature,
      levels = cindex_wide %>% arrange(Average) %>% pull(Signature)
    ),
    Highlight = ifelse(Signature == "ML_DDRscore", "DPRS", "Published")
  )

p_cindex_dot <- ggplot(
  cindex_dot,
  aes(x = Cindex, y = Signature, color = Dataset, shape = Highlight)
) +
  geom_vline(xintercept = 0.60, linetype = 2, color = "grey65", linewidth = 0.6) +
  geom_vline(xintercept = 0.65, linetype = 2, color = "grey45", linewidth = 0.6) +
  geom_point(size = 3.2, alpha = 0.95) +
  scale_shape_manual(values = c("DPRS" = 17, "Published" = 16)) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "C-index",
    y = NULL,
    title = "C-index comparison across cohorts"
  )

save_pdf(
  file.path(FIG_DIR, "Benchmark_Cindex_dotplot_by_dataset.pdf"),
  p_cindex_dot,
  w = 8,
  h = 5
)

############################
# 11F. Average C-index barplot with PMID labels
############################

bar_df_pmid <- cindex_wide %>%
  arrange(Average) %>%
  mutate(
    DisplayName = ifelse(
      Signature == "ML_DDRscore",
      "DPRS",
      paste0(Signature, "\nPMID: ", PMID)
    ),
    DisplayName = factor(DisplayName, levels = DisplayName),
    Highlight = ifelse(Signature == "ML_DDRscore", "DPRS", "Published")
  )

p_bar_pmid <- ggplot(bar_df_pmid, aes(x = Average, y = DisplayName, fill = Highlight)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = sprintf("%.3f", Average)),
    hjust = -0.08,
    size = 3.4,
    fontface = "bold"
  ) +
  geom_vline(xintercept = 0.65, linetype = 2, color = "grey45") +
  scale_fill_manual(values = c("DPRS" = "#C53030", "Published" = "grey75")) +
  coord_cartesian(xlim = c(0.45, max(bar_df_pmid$Average, na.rm = TRUE) + 0.06)) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold")
  ) +
  labs(
    x = "Average C-index across cohorts",
    y = NULL,
    title = "DPRS versus published DDR-related signatures"
  )

save_pdf(
  file.path(FIG_DIR, "Benchmark_average_Cindex_barplot_with_PMID.pdf"),
  p_bar_pmid,
  w = 8.5,
  h = 5.5
)

############################
# 11G. Mean timeROC AUC barplot
############################

if (nrow(auc_table) > 0) {
  auc_plot <- auc_table %>%
    group_by(Signature, Time) %>%
    summarise(MeanAUC = mean(AUC, na.rm = TRUE), .groups = "drop") %>%
    left_join(plot_df[, c("Signature", "Average")], by = "Signature") %>%
    arrange(desc(Average))
  
  auc_plot$Signature <- factor(
    auc_plot$Signature,
    levels = rev(unique(plot_df$Signature))
  )
  
  p_auc <- ggplot(auc_plot, aes(MeanAUC, Signature, fill = Time)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "black", linewidth = 0.2) +
    theme_classic(base_size = 11) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold")
    ) +
    labs(
      x = "Mean time-dependent AUC across cohorts",
      y = NULL,
      fill = NULL,
      title = "Mean time-dependent AUC comparison"
    )
  
  save_pdf(file.path(FIG_DIR, "Benchmark_mean_timeROC_AUC_barplot.pdf"), p_auc, 8, 5)
}

############################
# 12. KM plots for ML-DDRscore and top published signatures
############################

top_publish <- plot_df %>%
  filter(Signature != "ML_DDRscore") %>%
  slice_head(n = 3) %>%
  pull(Signature)

km_sigs <- c("ML_DDRscore", top_publish)

if (requireNamespace("survminer", quietly = TRUE)) {
  for (sig in km_sigs) {
    for (ds in c("Training", "Testing", "GSE72094", "GSE68465")) {
      dd <- all_scores %>% filter(Signature == sig, Dataset == ds)
      if (nrow(dd) < 40 || sum(dd$status == 1) < 8) next
      
      fit <- survival::survfit(Surv(time, status) ~ RiskGroup, data = dd)
      
      p <- survminer::ggsurvplot(
        fit,
        data = dd,
        pval = TRUE,
        risk.table = TRUE,
        risk.table.height = 0.23,
        palette = c("#2B6CB0", "#C53030"),
        legend.title = "",
        legend.labs = c("Low", "High"),
        title = paste0(ifelse(sig == "ML_DDRscore", "DPRS", sig), " - ", ds),
        xlab = "Time",
        ylab = "Overall survival probability",
        ggtheme = theme_bw(base_size = 12)
      )
      
      pdf(
        file.path(FIG_DIR, paste0("KM_", sig, "_", ds, ".pdf")),
        width = 5.6,
        height = 6.2,
        useDingbats = FALSE
      )
      print(p)
      dev.off()
    }
  }
} else {
  message("Package survminer is not installed. Skip KM plots.")
}

############################
# 13. Final message
############################

cat("\nPublished signature benchmark finished.\n")
cat("Output:\n", OUT_DIR, "\n")
cat("Main tables:\n")
cat(" - table/Published_signatures_cleaned.csv\n")
cat(" - table/Published_signatures_metadata.csv\n")
cat(" - table/Published_signature_gene_overlap.csv\n")
cat(" - table/All_signature_scores_with_ML_DDRscore.csv\n")
cat(" - table/Benchmark_Cindex_wide_summary.csv\n")
cat(" - table/Benchmark_Cindex_HR_by_dataset.csv\n")
cat(" - table/Benchmark_timeROC_AUC_long.csv\n")
cat(" - table/Benchmark_timeROC_curve_long.csv\n")
cat("Main figures:\n")
cat(" - fig/Benchmark_Cindex_heatmap.pdf\n")
cat(" - fig/Benchmark_average_Cindex_barplot.pdf\n")
cat(" - fig/Benchmark_external_Cindex_barplot.pdf\n")
cat(" - fig/Benchmark_Cindex_dotplot_by_dataset.pdf\n")
cat(" - fig/Benchmark_average_Cindex_barplot_with_PMID.pdf\n")
cat(" - fig/Benchmark_ROC_compare_*_1y_3y.pdf\n")
cat(" - fig/Benchmark_ROC_compare_*_1y_3y_5y.pdf\n")