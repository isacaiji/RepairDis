############################################################
# 03_ML_DPRS_published_signature_benchmark.R
#
# Purpose:
#   Compare DPRS with published DDR-related prognostic signatures.
#
# Required before running:
#   1. Run 01_ML_DPRS_modeling_main.R first.
#   2. Current R session should contain:
#        tcga_expr: gene x sample expression matrix
#        tcga_clin: Patient, time, status
#   3. Put published signature table here:
#        D:/R_workspace/评分/AD_DDR_project/00_data/Published_signatures.csv
#
# Published_signatures.csv columns:
#   Signature,Gene,Coefficient,Source,Type,Year,Endpoint,Note,PMID
#
# Input:
#   D:/R_workspace/评分/AD_DDR_project/04-ML
#
# Output:
#   D:/R_workspace/评分/AD_DDR_project/04-ML/03_Published_signature_benchmark
############################################################

options(stringsAsFactors = FALSE)
set.seed(20260513)

############################
# 0. Paths
############################

PROJECT_DIR <- "D:/R_workspace/评分/AD_DDR_project"
DATA_DIR <- file.path(PROJECT_DIR, "00_data")
ML_DIR <- file.path(PROJECT_DIR, "04-ML")

OUT_DIR <- file.path(ML_DIR, "03_Published_signature_benchmark")
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

SIGNATURE_FILE <- file.path(DATA_DIR, "Published_signatures.csv")

############################
# 1. Packages
############################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(survival)
})

############################
# 2. Helper functions
############################

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

make_display_name <- function(sig, pmid_map = NULL) {
  if (sig == "ML_DDRscore") return("DPRS")
  
  if (!is.null(pmid_map) && sig %in% names(pmid_map)) {
    pmid <- pmid_map[[sig]]
    if (!is.na(pmid) && pmid != "") {
      return(paste0(sig, "\nPMID: ", pmid))
    }
  }
  
  sig
}

gene_alias_map <- c(
  "HIST3H2A" = "H3C15"
)

############################
# 3. Required objects
############################

stopifnot(exists("tcga_expr"))
stopifnot(exists("tcga_clin"))

tcga_expr <- as.matrix(tcga_expr)
rownames(tcga_expr) <- clean_gene(rownames(tcga_expr))
tcga_expr <- tcga_expr[!is.na(rownames(tcga_expr)) & rownames(tcga_expr) != "", , drop = FALSE]
tcga_expr <- tcga_expr[!duplicated(rownames(tcga_expr)), , drop = FALSE]

tcga_clin$Patient <- as.character(tcga_clin$Patient)

cat("TCGA expression:", nrow(tcga_expr), "genes x", ncol(tcga_expr), "samples\n")
cat("TCGA clinical:", nrow(tcga_clin), "patients\n")

############################
# 4. Read published signatures
############################

if (!file.exists(SIGNATURE_FILE)) {
  stop("Cannot find: ", SIGNATURE_FILE)
}

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
  dplyr::select(dplyr::all_of(required_sig_cols)) %>%
  dplyr::mutate(dplyr::across(everything(), ~ trimws(as.character(.x)))) %>%
  dplyr::filter(
    !(is.na(Signature) | Signature == "") |
      !(is.na(Gene) | Gene == "")
  ) %>%
  dplyr::mutate(
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
  dplyr::filter(
    !is.na(Signature),
    Signature != "",
    !is.na(Gene),
    Gene != "",
    is.finite(Coefficient)
  ) %>%
  dplyr::distinct(Signature, Gene, .keep_all = TRUE)

# 当前主比较保留的 published signatures
keep_sigs <- c(
  "Li_2022_NSCLC_DDR6",
  "Zhao_2022_EarlyLUAD_DRG16",
  "Chen_2021_LUAD_DRG6",
  "Li_2022_RespRes_DCG4",
  "Hu_2020_LUAD_DNARepair13",
  "Yang_2020_SurgLUAD_DNARepair6"
)

sig_tbl <- sig_tbl %>%
  dplyr::filter(Signature %in% keep_sigs)

sig_meta <- sig_tbl %>%
  dplyr::group_by(Signature) %>%
  dplyr::summarise(
    Source = dplyr::first(Source),
    Type = dplyr::first(Type),
    Year = dplyr::first(Year),
    Endpoint = dplyr::first(Endpoint),
    Note = dplyr::first(Note),
    PMID = dplyr::first(PMID),
    SignatureGeneN = dplyr::n_distinct(Gene),
    Genes = paste(unique(Gene), collapse = ";"),
    .groups = "drop"
  )

save_csv(sig_tbl, file.path(TAB_DIR, "Published_signatures_cleaned.csv"))
save_csv(sig_meta, file.path(TAB_DIR, "Published_signatures_metadata.csv"))

cat("\nPublished signatures loaded:\n")
print(table(sig_tbl$Signature))

############################
# 5. Read DPRS / split files from ML_DIR
############################

final_model_info_file <- file.path(ML_DIR, "Fig4B_final_selected_model_info.csv")
train_file <- file.path(ML_DIR, "Mime1_input_Training.csv")
test_file <- file.path(ML_DIR, "Mime1_input_Testing.csv")
dprs_file <- file.path(ML_DIR, "Fig4C_ML_DDRscore_all_sets.csv")

if (!file.exists(final_model_info_file)) stop("Cannot find: ", final_model_info_file)
if (!file.exists(train_file)) stop("Cannot find: ", train_file)
if (!file.exists(test_file)) stop("Cannot find: ", test_file)
if (!file.exists(dprs_file)) stop("Cannot find: ", dprs_file)

global_best <- data.table::fread(
  final_model_info_file,
  data.table = FALSE,
  check.names = FALSE
)

best_model <- as.character(global_best$Model[1])

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

cat("\nBest ML model:", best_model, "\n")
cat("Training patients:", length(train_patients), "\n")
cat("Testing patients:", length(test_patients), "\n")

############################
# 6. Build TCGA expression by patient
############################

make_tcga_expr_by_patient <- function(patient_ids) {
  
  all_samples <- colnames(tcga_expr)
  all_patients <- patient_id(all_samples)
  
  use_samples <- sapply(patient_ids, function(pid) {
    smp <- all_samples[all_patients == pid]
    if (length(smp) == 0) return(NA_character_)
    
    st <- sample_type_code(smp)
    smp_tumor <- smp[st == "01"]
    if (length(smp_tumor) > 0) return(smp_tumor[1])
    
    smp[1]
  })
  
  valid <- !is.na(use_samples)
  patient_ids <- patient_ids[valid]
  use_samples <- use_samples[valid]
  
  expr <- tcga_expr[, use_samples, drop = FALSE]
  colnames(expr) <- patient_ids
  
  expr_num <- expr
  storage.mode(expr_num) <- "numeric"
  
  if (quantile(expr_num, 0.95, na.rm = TRUE) > 50) {
    expr_log <- log2(expr_num + 1)
  } else {
    expr_log <- expr_num
  }
  
  clin <- tcga_clin %>%
    dplyr::filter(Patient %in% patient_ids) %>%
    dplyr::distinct(Patient, .keep_all = TRUE)
  
  list(expr = expr_log, clin = clin)
}

tcga_train <- make_tcga_expr_by_patient(train_patients)
tcga_test <- make_tcga_expr_by_patient(test_patients)

cat("TCGA Training:", nrow(tcga_train$expr), "genes x", ncol(tcga_train$expr), "samples\n")
cat("TCGA Testing:", nrow(tcga_test$expr), "genes x", ncol(tcga_test$expr), "samples\n")

############################
# 7. Build GEO expression objects
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
# 8. Calculate published signature scores
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
    
    ss <- sig_tbl %>% dplyr::filter(Signature == sig)
    
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
      PMID = dplyr::first(ss$PMID),
      Source = dplyr::first(ss$Source),
      Type = dplyr::first(ss$Type),
      Year = dplyr::first(ss$Year),
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
  
  score_df <- dplyr::bind_rows(score_list)
  overlap_df <- dplyr::bind_rows(overlap_list)
  
  if (nrow(score_df) == 0) {
    return(list(score = data.frame(), overlap = overlap_df))
  }
  
  if (dataset_name %in% c("Training", "Testing")) {
    clin2 <- clin %>%
      dplyr::transmute(
        ID = as.character(Patient),
        time = as.numeric(time),
        status = as.numeric(status)
      )
  } else {
    clin2 <- clin %>%
      dplyr::transmute(
        ID = as.character(Sample),
        time = as.numeric(time),
        status = as.numeric(status)
      )
  }
  
  score_df <- score_df %>%
    dplyr::left_join(clin2, by = "ID") %>%
    dplyr::filter(
      is.finite(time),
      time > 0,
      status %in% c(0, 1),
      is.finite(Score)
    )
  
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

published_score <- dplyr::bind_rows(all_score)
published_overlap <- dplyr::bind_rows(all_overlap)

save_csv(
  published_score,
  file.path(TAB_DIR, "Published_signature_scores_long_raw.csv")
)

save_csv(
  published_overlap,
  file.path(TAB_DIR, "Published_signature_gene_overlap.csv")
)

############################
# 9. Add DPRS score
############################

our_score <- data.table::fread(
  dprs_file,
  data.table = FALSE,
  check.names = FALSE
)

score_col <- intersect(
  c("ML_DDRscore", "riskscore", "RiskScore", "risk_score"),
  colnames(our_score)
)[1]

id_col <- intersect(
  c("ID", "Sample", "Patient"),
  colnames(our_score)
)[1]

if (is.na(score_col)) stop("Cannot identify DPRS score column.")
if (is.na(id_col)) stop("Cannot identify DPRS ID column.")

our_score2 <- our_score %>%
  dplyr::transmute(
    ID = as.character(.data[[id_col]]),
    time = as.numeric(time),
    status = as.numeric(status),
    Dataset = as.character(Dataset),
    Signature = "ML_DDRscore",
    Score = as.numeric(.data[[score_col]]),
    Model = if ("Model" %in% colnames(our_score)) {
      as.character(.data[["Model"]])
    } else {
      best_model
    }
  ) %>%
  dplyr::filter(
    is.finite(time),
    time > 0,
    status %in% c(0, 1),
    is.finite(Score)
  )

published_score$Model <- published_score$Signature

all_scores <- dplyr::bind_rows(
  our_score2,
  published_score
)

############################
# 10. Direction correction and median grouping
############################

all_scores <- all_scores %>%
  dplyr::group_by(Signature) %>%
  dplyr::group_modify(~{
    
    x <- .x
    tr <- x %>% dplyr::filter(Dataset == "Training")
    
    fit <- tryCatch(
      survival::coxph(survival::Surv(time, status) ~ Score, data = tr),
      error = function(e) NULL
    )
    
    if (!is.null(fit) &&
        length(coef(fit)) > 0 &&
        is.finite(as.numeric(coef(fit)[1]))) {
      if (as.numeric(coef(fit)[1]) < 0) {
        x$Score <- -x$Score
      }
    }
    
    x
  }) %>%
  dplyr::ungroup()

all_scores <- all_scores %>%
  dplyr::group_by(Signature, Dataset) %>%
  dplyr::mutate(
    Cutoff = median(Score, na.rm = TRUE),
    RiskGroup = ifelse(Score >= Cutoff, "High", "Low")
  ) %>%
  dplyr::ungroup()

all_scores$RiskGroup <- factor(all_scores$RiskGroup, levels = c("Low", "High"))

save_csv(
  all_scores,
  file.path(TAB_DIR, "All_signature_scores_with_DPRS.csv")
)

############################
# 11. C-index and Cox HR
############################

metric_table <- all_scores %>%
  dplyr::group_by(Signature, Dataset) %>%
  dplyr::summarise(
    N = dplyr::n(),
    Events = sum(status == 1),
    Cindex = calc_cindex(time, status, Score),
    HR = {
      fit <- tryCatch(
        survival::coxph(survival::Surv(time, status) ~ RiskGroup, data = dplyr::cur_data()),
        error = function(e) NULL
      )
      if (is.null(fit)) NA_real_ else as.numeric(summary(fit)$coefficients[1, "exp(coef)"])
    },
    CoxP = {
      fit <- tryCatch(
        survival::coxph(survival::Surv(time, status) ~ RiskGroup, data = dplyr::cur_data()),
        error = function(e) NULL
      )
      if (is.null(fit)) NA_real_ else as.numeric(summary(fit)$coefficients[1, "Pr(>|z|)"])
    },
    .groups = "drop"
  )

metric_table <- metric_table %>%
  dplyr::left_join(sig_meta, by = "Signature") %>%
  dplyr::mutate(
    Source = ifelse(Signature == "ML_DDRscore", "This_study", Source),
    Type = ifelse(Signature == "ML_DDRscore", "MO_DDRscore_derived_ML_model", Type),
    Year = ifelse(Signature == "ML_DDRscore", 2026, Year),
    Endpoint = ifelse(Signature == "ML_DDRscore", "OS", Endpoint),
    Note = ifelse(Signature == "ML_DDRscore", "this_study", Note),
    PMID = ifelse(Signature == "ML_DDRscore", NA, PMID)
  )

save_csv(
  metric_table,
  file.path(TAB_DIR, "Benchmark_Cindex_HR_by_dataset.csv")
)

cindex_wide <- metric_table %>%
  dplyr::select(Signature, Dataset, Cindex) %>%
  tidyr::pivot_wider(names_from = Dataset, values_from = Cindex)

dataset_cols <- intersect(
  c("Training", "Testing", "GSE72094", "GSE68465"),
  colnames(cindex_wide)
)

cindex_wide <- cindex_wide %>%
  dplyr::mutate(
    Average = rowMeans(dplyr::across(dplyr::all_of(dataset_cols)), na.rm = TRUE),
    ExternalMean = rowMeans(
      dplyr::across(
        dplyr::all_of(intersect(c("GSE72094", "GSE68465"), dataset_cols))
      ),
      na.rm = TRUE
    ),
    MinValidation = pmin(Testing, GSE72094, GSE68465, na.rm = TRUE)
  ) %>%
  dplyr::left_join(sig_meta, by = "Signature") %>%
  dplyr::mutate(
    Source = ifelse(Signature == "ML_DDRscore", "This_study", Source),
    Type = ifelse(Signature == "ML_DDRscore", "MO_DDRscore_derived_ML_model", Type),
    Year = ifelse(Signature == "ML_DDRscore", 2026, Year),
    Endpoint = ifelse(Signature == "ML_DDRscore", "OS", Endpoint),
    Note = ifelse(Signature == "ML_DDRscore", "this_study", Note),
    PMID = ifelse(Signature == "ML_DDRscore", NA, PMID)
  ) %>%
  dplyr::arrange(dplyr::desc(Average), dplyr::desc(ExternalMean))

save_csv(
  cindex_wide,
  file.path(TAB_DIR, "Benchmark_Cindex_wide_summary.csv")
)

cat("\nBenchmark C-index summary:\n")
print(cindex_wide)

############################
# 12. timeROC AUC and ROC coordinates
############################

auc_list <- list()
roc_curve_list <- list()

if (requireNamespace("timeROC", quietly = TRUE)) {
  
  for (sig in unique(all_scores$Signature)) {
    for (ds in unique(all_scores$Dataset)) {
      
      dd <- all_scores %>%
        dplyr::filter(Signature == sig, Dataset == ds) %>%
        dplyr::filter(is.finite(time), is.finite(status), is.finite(Score))
      
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
          ) %>%
            dplyr::filter(is.finite(FP), is.finite(TP)) %>%
            dplyr::arrange(FP, TP)
          
          roc_curve_list[[paste(sig, ds, time_labels[i], sep = "__")]] <- curve_now
        }
      }
    }
  }
  
} else {
  message("Package timeROC is not installed. Skip time-dependent ROC.")
}

auc_table <- dplyr::bind_rows(auc_list)
roc_curve_long <- dplyr::bind_rows(roc_curve_list)

save_csv(
  auc_table,
  file.path(TAB_DIR, "Benchmark_timeROC_AUC_long.csv")
)

save_csv(
  roc_curve_long,
  file.path(TAB_DIR, "Benchmark_timeROC_curve_long.csv")
)

############################
# 13. Plot: C-index comparison heatmap
############################

pmid_map <- sig_meta$PMID
names(pmid_map) <- sig_meta$Signature

plot_cindex <- cindex_wide %>%
  dplyr::mutate(
    DisplayName = vapply(Signature, make_display_name, character(1), pmid_map = pmid_map),
    DisplayName = ifelse(Signature == "ML_DDRscore", "DPRS", DisplayName),
    IsDPRS = Signature == "ML_DDRscore"
  ) %>%
  dplyr::arrange(dplyr::desc(IsDPRS), dplyr::desc(Average))

plot_cindex_long <- plot_cindex %>%
  dplyr::select(Signature, DisplayName, IsDPRS, dplyr::all_of(dataset_cols), Average) %>%
  tidyr::pivot_longer(
    cols = c(dplyr::all_of(dataset_cols), Average),
    names_to = "Dataset",
    values_to = "Cindex"
  )

plot_cindex_long$DisplayName <- factor(
  plot_cindex_long$DisplayName,
  levels = rev(plot_cindex$DisplayName)
)

plot_cindex_long$Dataset <- factor(
  plot_cindex_long$Dataset,
  levels = c(dataset_cols, "Average")
)

p_cindex_heatmap <- ggplot(
  plot_cindex_long,
  aes(x = Dataset, y = DisplayName, fill = Cindex)
) +
  geom_tile(color = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.3f", Cindex)), size = 3.2, color = "black") +
  scale_fill_gradient2(
    low = "#4195C1",
    mid = "white",
    high = "#CB5746",
    midpoint = 0.65,
    limits = c(0.45, max(plot_cindex_long$Cindex, na.rm = TRUE) + 0.02),
    name = "C-index"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5, face = "bold"),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "DPRS versus published DDR-related signatures"
  )

save_pdf(
  file.path(FIG_DIR, "Benchmark_Cindex_heatmap.pdf"),
  p_cindex_heatmap,
  8.5,
  max(4.5, 0.42 * length(unique(plot_cindex_long$DisplayName)) + 2)
)

############################
# 14. Plot: Average C-index barplot
############################

p_avg_cindex <- plot_cindex %>%
  dplyr::mutate(
    DisplayName = factor(DisplayName, levels = rev(DisplayName)),
    Group = ifelse(Signature == "ML_DDRscore", "DPRS", "Published signatures")
  ) %>%
  ggplot(aes(x = Average, y = DisplayName, fill = Group)) +
  geom_col(width = 0.70, color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = sprintf("%.3f", Average)),
    hjust = -0.08,
    size = 3.3,
    fontface = "bold"
  ) +
  scale_fill_manual(values = c("DPRS" = "#C53030", "Published signatures" = "grey75")) +
  coord_cartesian(xlim = c(0.45, max(plot_cindex$Average, na.rm = TRUE) + 0.08)) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "Average C-index",
    y = NULL,
    title = "Average C-index comparison"
  )

save_pdf(
  file.path(FIG_DIR, "Benchmark_average_Cindex_barplot.pdf"),
  p_avg_cindex,
  8,
  max(4.5, 0.42 * nrow(plot_cindex) + 2)
)

############################
# 15. Plot: timeROC AUC summary
############################

if (nrow(auc_table) > 0) {
  
  auc_plot <- auc_table %>%
    dplyr::mutate(
      DisplayName = ifelse(Signature == "ML_DDRscore", "DPRS", Signature),
      Group = ifelse(Signature == "ML_DDRscore", "DPRS", "Published signatures")
    )
  
  save_csv(
    auc_plot,
    file.path(TAB_DIR, "Benchmark_timeROC_AUC_for_plot.csv")
  )
  
  p_auc <- ggplot(
    auc_plot,
    aes(x = Time, y = AUC, fill = DisplayName)
  ) +
    geom_col(
      position = position_dodge(width = 0.78),
      width = 0.70,
      color = "black",
      linewidth = 0.18
    ) +
    facet_wrap(~Dataset, nrow = 1) +
    coord_cartesian(ylim = c(0.45, max(auc_plot$AUC, na.rm = TRUE) + 0.06)) +
    theme_classic(base_size = 11) +
    theme(
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(color = "black", face = "bold"),
      legend.title = element_blank(),
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      x = NULL,
      y = "timeROC AUC",
      title = "Time-dependent AUC comparison"
    )
  
  save_pdf(
    file.path(FIG_DIR, "Benchmark_timeROC_AUC_barplot.pdf"),
    p_auc,
    12,
    5.5
  )
}

############################
# 16. Plot: ROC comparison curves
############################

plot_roc_compare <- function(dataset_to_plot,
                             times_to_plot = c("1 year", "3 years", "5 years"),
                             out_file) {
  
  if (nrow(roc_curve_long) == 0) return(NULL)
  
  dd <- roc_curve_long %>%
    dplyr::filter(
      Dataset == dataset_to_plot,
      Time %in% times_to_plot
    ) %>%
    dplyr::mutate(
      DisplayName = ifelse(Signature == "ML_DDRscore", "DPRS", Signature),
      Label = paste0(DisplayName, " | ", Time, " AUC=", sprintf("%.3f", AUC))
    )
  
  if (nrow(dd) == 0) return(NULL)
  
  p <- ggplot(dd, aes(x = FP, y = TP, color = Label)) +
    geom_line(linewidth = 0.85, alpha = 0.95) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold"),
      legend.title = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 6.8),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      x = "False positive rate",
      y = "True positive rate",
      title = paste0(dataset_to_plot, " ROC comparison")
    )
  
  save_pdf(out_file, p, 8.8, 5.8)
}

if (nrow(roc_curve_long) > 0) {
  
  for (ds in c("Training", "Testing", "GSE72094", "GSE68465")) {
    
    plot_roc_compare(
      dataset_to_plot = ds,
      times_to_plot = c("1 year", "3 years"),
      out_file = file.path(FIG_DIR, paste0("Benchmark_ROC_compare_", ds, "_1y_3y.pdf"))
    )
    
    plot_roc_compare(
      dataset_to_plot = ds,
      times_to_plot = c("1 year", "3 years", "5 years"),
      out_file = file.path(FIG_DIR, paste0("Benchmark_ROC_compare_", ds, "_1y_3y_5y.pdf"))
    )
  }
}

############################
# 17. Optional KM plots for DPRS and top published signatures
############################

top_publish <- cindex_wide %>%
  dplyr::filter(Signature != "ML_DDRscore") %>%
  dplyr::arrange(dplyr::desc(Average)) %>%
  dplyr::slice_head(n = 3) %>%
  dplyr::pull(Signature)

km_sigs <- c("ML_DDRscore", top_publish)

if (requireNamespace("survminer", quietly = TRUE)) {
  
  for (sig in km_sigs) {
    for (ds in c("Training", "Testing", "GSE72094", "GSE68465")) {
      
      dd <- all_scores %>%
        dplyr::filter(Signature == sig, Dataset == ds)
      
      if (nrow(dd) < 40 || sum(dd$status == 1) < 8) next
      
      fit <- survival::survfit(
        survival::Surv(time, status) ~ RiskGroup,
        data = dd
      )
      
      title_sig <- ifelse(sig == "ML_DDRscore", "DPRS", sig)
      
      p_km <- survminer::ggsurvplot(
        fit,
        data = dd,
        pval = TRUE,
        risk.table = FALSE,
        palette = c("#2B6CB0", "#C53030"),
        title = paste0(title_sig, " - ", ds),
        legend.title = "",
        legend.labs = c("Low", "High"),
        xlab = "Time",
        ylab = "Overall survival probability",
        ggtheme = theme_bw(base_size = 12)
      )$plot
      
      safe_name <- gsub("[^A-Za-z0-9]+", "_", paste(sig, ds, sep = "_"))
      
      save_pdf(
        file.path(FIG_DIR, paste0("Benchmark_KM_", safe_name, ".pdf")),
        p_km,
        5.8,
        5
      )
    }
  }
}

############################
# 18. Output manifest
############################

benchmark_manifest <- data.frame(
  File = c(
    "Published_signatures_cleaned.csv",
    "Published_signatures_metadata.csv",
    "Published_signature_gene_overlap.csv",
    "All_signature_scores_with_DPRS.csv",
    "Benchmark_Cindex_HR_by_dataset.csv",
    "Benchmark_Cindex_wide_summary.csv",
    "Benchmark_timeROC_AUC_long.csv",
    "Benchmark_timeROC_curve_long.csv",
    "Benchmark_Cindex_heatmap.pdf",
    "Benchmark_average_Cindex_barplot.pdf",
    "Benchmark_timeROC_AUC_barplot.pdf",
    "Benchmark_ROC_compare_*",
    "Benchmark_KM_*"
  ),
  Description = c(
    "Cleaned published signature table used for comparison",
    "Metadata of published signatures",
    "Gene overlap between signatures and each cohort expression matrix",
    "Combined DPRS and published signature scores",
    "C-index, HR and Cox P value by signature and dataset",
    "Wide C-index summary table",
    "Long table of time-dependent AUC values",
    "ROC curve coordinates for plotting",
    "Heatmap comparing C-index values",
    "Average C-index barplot",
    "Time-dependent AUC barplot",
    "ROC comparison curves by cohort and time point",
    "KM curves for DPRS and top published signatures"
  ),
  Path = c(
    file.path(TAB_DIR, "Published_signatures_cleaned.csv"),
    file.path(TAB_DIR, "Published_signatures_metadata.csv"),
    file.path(TAB_DIR, "Published_signature_gene_overlap.csv"),
    file.path(TAB_DIR, "All_signature_scores_with_DPRS.csv"),
    file.path(TAB_DIR, "Benchmark_Cindex_HR_by_dataset.csv"),
    file.path(TAB_DIR, "Benchmark_Cindex_wide_summary.csv"),
    file.path(TAB_DIR, "Benchmark_timeROC_AUC_long.csv"),
    file.path(TAB_DIR, "Benchmark_timeROC_curve_long.csv"),
    file.path(FIG_DIR, "Benchmark_Cindex_heatmap.pdf"),
    file.path(FIG_DIR, "Benchmark_average_Cindex_barplot.pdf"),
    file.path(FIG_DIR, "Benchmark_timeROC_AUC_barplot.pdf"),
    file.path(FIG_DIR, "Benchmark_ROC_compare_*"),
    file.path(FIG_DIR, "Benchmark_KM_*")
  ),
  stringsAsFactors = FALSE
)

save_csv(
  benchmark_manifest,
  file.path(OUT_DIR, "ML_DPRS_published_signature_benchmark_manifest.csv")
)

sink(file.path(OUT_DIR, "ML_DPRS_published_signature_benchmark_session_info.txt"))
cat("DPRS published signature benchmark finished.\n")
cat("Input ML directory:", ML_DIR, "\n")
cat("Output directory:", OUT_DIR, "\n")
cat("Best model:", best_model, "\n\n")
cat("Published signatures:\n")
print(sig_meta)
cat("\nC-index summary:\n")
print(cindex_wide)
cat("\nOutput manifest:\n")
print(benchmark_manifest)
cat("\nSession info:\n")
print(sessionInfo())
sink()

cat("\nDPRS published signature benchmark finished.\n")
cat("Output directory:", OUT_DIR, "\n")
cat("\nMain outputs:\n")
cat(" - tables/Benchmark_Cindex_wide_summary.csv\n")
cat(" - tables/Benchmark_timeROC_AUC_long.csv\n")
cat(" - figures/Benchmark_Cindex_heatmap.pdf\n")
cat(" - figures/Benchmark_average_Cindex_barplot.pdf\n")
cat(" - figures/Benchmark_timeROC_AUC_barplot.pdf\n")
cat(" - figures/Benchmark_ROC_compare_*.pdf\n")
cat(" - figures/Benchmark_KM_*.pdf\n")