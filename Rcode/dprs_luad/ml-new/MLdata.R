############################################################
# 01_ML_DPRS_modeling_main.R
#
# Purpose:
#   Build ML-based DPRS prognostic model and generate data tables only.
#
# Output:
#   1. Training/Test split files
#   2. Training-only DEG table
#   3. Training-only uni-Cox table
#   4. Candidate gene table
#   5. Mime1 input tables
#   6. Mime1 result RDS
#   7. C-index summary table
#   8. Final selected model info
#   9. DPRS score table
#   10. Final signature genes table
#
# No plotting in this script.
############################################################

options(stringsAsFactors = FALSE)

############################
# 0. Parameters
############################

PROJECT_DIR <- "D:/R_workspace/评分/AD_DDR_project"
DATA_DIR    <- file.path(PROJECT_DIR, "00_data")
DB_DIR      <- file.path(PROJECT_DIR, "05_database_tables")
FIG4_DIR    <- file.path(PROJECT_DIR, "04-ML")

dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG4_DIR, recursive = TRUE, showWarnings = FALSE)

TRAIN_RATIO <- 0.80
DEG_ADJ_P   <- 0.05
DEG_LOGFC   <- 1.20
UNICOX_P    <- 0.03
MAX_ML_GENES <- 60
NODE_SIZE   <- 20
SEED        <- 20260513

TARGET_MODEL <- "StepCox[forward] + RSF"

set.seed(SEED)

############################
# 1. Packages
############################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(survival)
  library(limma)
  library(Mime1)
})

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
  x[x %in% c("", "NA", "---", "NULL", "N/A", "NAN")] <- NA
  x
}

patient_id <- function(x) {
  substr(gsub("\\.", "-", as.character(x)), 1, 12)
}

zscore_vector <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  sdx <- sd(x, na.rm = TRUE)
  if (all(is.na(x)) || !is.finite(sdx) || sdx == 0) {
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

read_geo_genes <- function(expr_file) {
  x <- data.table::fread(expr_file, data.table = FALSE, check.names = FALSE)
  unique(na.omit(clean_gene(x[[1]])))
}

############################
# 3. Required objects
############################
# Must exist in current R session:
#   tcga_expr: gene x sample expression matrix
#   tumor_samples: TCGA tumor sample IDs
#   tcga_clin: Patient, time, status, age, gender, stage
#   mo_score: Sample, Patient, SampleType, MO_DDRscore_group, MO_DDRscore_raw

stopifnot(exists("tcga_expr"))
stopifnot(exists("tumor_samples"))
stopifnot(exists("tcga_clin"))
stopifnot(exists("mo_score"))

tcga_expr <- as.matrix(tcga_expr)
rownames(tcga_expr) <- clean_gene(rownames(tcga_expr))
tcga_expr <- tcga_expr[!is.na(rownames(tcga_expr)) & rownames(tcga_expr) != "", , drop = FALSE]
tcga_expr <- tcga_expr[!duplicated(rownames(tcga_expr)), , drop = FALSE]

tcga_clin$Patient <- as.character(tcga_clin$Patient)
mo_score$Patient <- as.character(mo_score$Patient)

score_df <- mo_score %>%
  dplyr::filter(
    Sample %in% tumor_samples,
    !is.na(MO_DDRscore_group),
    is.finite(MO_DDRscore_raw)
  ) %>%
  dplyr::arrange(Patient, SampleType) %>%
  dplyr::distinct(Patient, .keep_all = TRUE) %>%
  dplyr::inner_join(
    tcga_clin[, c("Patient", "time", "status", "age", "gender", "stage")],
    by = "Patient"
  ) %>%
  dplyr::filter(
    is.finite(time),
    time > 0,
    status %in% c(0, 1)
  )

save_csv(
  score_df,
  file.path(FIG4_DIR, "Fig4A_TCGA_score_survival_table.csv")
)

cat("TCGA patients used:", nrow(score_df), "\n")
cat("Events:", sum(score_df$status == 1), "\n")

############################
# 4. Train/Test split
############################

set.seed(SEED)

event_ids <- score_df$Patient[score_df$status == 1]
cens_ids  <- score_df$Patient[score_df$status == 0]

train_ids <- c(
  sample(event_ids, floor(TRAIN_RATIO * length(event_ids))),
  sample(cens_ids,  floor(TRAIN_RATIO * length(cens_ids)))
)

test_ids <- setdiff(score_df$Patient, train_ids)

train_samples <- score_df$Sample[score_df$Patient %in% train_ids]
test_samples  <- score_df$Sample[score_df$Patient %in% test_ids]

split_info <- data.frame(
  Patient = score_df$Patient,
  Sample = score_df$Sample,
  Dataset = ifelse(score_df$Patient %in% train_ids, "Training", "Testing"),
  time = score_df$time,
  status = score_df$status,
  MO_DDRscore_raw = score_df$MO_DDRscore_raw,
  MO_DDRscore_group = score_df$MO_DDRscore_group,
  stringsAsFactors = FALSE
)

save_csv(
  split_info,
  file.path(FIG4_DIR, "Fig4A_train_test_split_info.csv")
)

cat("Training:", length(train_ids), "patients /",
    sum(score_df$status[score_df$Patient %in% train_ids] == 1), "events\n")

cat("Testing:", length(test_ids), "patients /",
    sum(score_df$status[score_df$Patient %in% test_ids] == 1), "events\n")

############################
# 5. Training-only DEG
############################

train_group <- factor(
  score_df$MO_DDRscore_group[match(train_samples, score_df$Sample)],
  levels = c("Low", "High")
)

expr_train <- log2(tcga_expr[, train_samples, drop = FALSE] + 1)

design <- model.matrix(~0 + train_group)
colnames(design) <- c("Low", "High")

fit <- limma::lmFit(expr_train, design)
cont <- limma::makeContrasts(High - Low, levels = design)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont))

train_deg <- limma::topTable(fit2, number = Inf, adjust.method = "BH")
train_deg$Gene <- rownames(train_deg)
train_deg$FDR <- train_deg$adj.P.Val

train_deg <- train_deg %>%
  dplyr::select(Gene, everything()) %>%
  dplyr::arrange(FDR, dplyr::desc(abs(logFC)))

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

if (length(deg_genes) < 30) {
  deg_genes <- train_deg %>%
    dplyr::arrange(FDR, dplyr::desc(abs(logFC))) %>%
    dplyr::slice_head(n = 500) %>%
    dplyr::pull(Gene) %>%
    unique()
}

cat("DEG genes after fallback:", length(deg_genes), "\n")

############################
# 6. GEO common genes
############################

gse72094_expr_file <- file.path(DATA_DIR, "GSE72094_expression.csv")
gse72094_clin_file <- file.path(DATA_DIR, "GSE72094_clinical.csv")
gse68465_expr_file <- file.path(DATA_DIR, "GSE68465_expression.csv")
gse68465_clin_file <- file.path(DATA_DIR, "GSE68465_clinical.csv")

geo_common_genes <- Reduce(
  intersect,
  list(
    rownames(tcga_expr),
    read_geo_genes(gse72094_expr_file),
    read_geo_genes(gse68465_expr_file)
  )
)

deg_genes <- intersect(deg_genes, geo_common_genes)

if (length(deg_genes) < 10) {
  stop("Too few DEG genes after intersecting with GEO common genes: ", length(deg_genes))
}

save_csv(
  data.frame(Gene = deg_genes),
  file.path(FIG4_DIR, "Fig4A_DEG_genes_common_with_GEO.csv")
)

cat("DEG genes common with GEO:", length(deg_genes), "\n")

############################
# 7. Training-only uni-Cox filtering
############################

train_expr_cox <- as.data.frame(
  t(log2(tcga_expr[deg_genes, train_samples, drop = FALSE] + 1)),
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

cox_table <- dplyr::bind_rows(lapply(deg_genes, function(g) {
  
  df <- train_surv[, c("time", "status", g)]
  colnames(df) <- c("time", "status", "expr")
  
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

if (nrow(cox_table) == 0) {
  stop("No valid uni-Cox result was generated.")
}

cox_table$CoxFDR <- p.adjust(cox_table$CoxP, method = "BH")

cox_table <- cox_table %>%
  dplyr::arrange(CoxP)

save_csv(
  cox_table,
  file.path(FIG4_DIR, "Fig4A_training_only_unicox.csv")
)

candidate_genes <- cox_table %>%
  dplyr::filter(CoxP < UNICOX_P) %>%
  dplyr::pull(Gene) %>%
  unique()

if (length(candidate_genes) < 10) {
  candidate_genes <- head(cox_table$Gene, 50)
}

candidate_genes <- head(unique(candidate_genes), MAX_ML_GENES)

save_csv(
  data.frame(Gene = candidate_genes),
  file.path(FIG4_DIR, "Fig4A_candidate_genes_train_only.csv")
)

cat("Candidate genes:", length(candidate_genes), "\n")

############################
# 8. Build TCGA Mime input
############################

make_tcga_mime <- function(samples, gene_use) {
  
  x <- as.data.frame(
    t(log2(tcga_expr[gene_use, samples, drop = FALSE] + 1)),
    check.names = FALSE
  )
  
  x$Sample <- rownames(x)
  x$Patient <- patient_id(x$Sample)
  
  x <- x %>%
    dplyr::arrange(Patient, Sample) %>%
    dplyr::distinct(Patient, .keep_all = TRUE) %>%
    dplyr::inner_join(
      tcga_clin[, c("Patient", "time", "status")],
      by = "Patient"
    ) %>%
    dplyr::filter(
      is.finite(time),
      time > 0,
      status %in% c(0, 1)
    )
  
  out <- x[, c("Patient", "time", "status", gene_use)]
  colnames(out)[1:3] <- c("ID", "OS.time", "OS")
  
  for (g in gene_use) {
    out[[g]] <- zscore_vector(out[[g]])
  }
  
  out
}

mime_all <- make_tcga_mime(score_df$Sample, candidate_genes)

train_df_mime <- mime_all[mime_all$ID %in% train_ids, ]
test_df_mime  <- mime_all[mime_all$ID %in% test_ids, ]

save_csv(
  train_df_mime,
  file.path(FIG4_DIR, "Mime1_input_Training.csv")
)

save_csv(
  test_df_mime,
  file.path(FIG4_DIR, "Mime1_input_Testing.csv")
)

cat("Mime Training:", nrow(train_df_mime), "samples\n")
cat("Mime Testing:", nrow(test_df_mime), "samples\n")

############################
# 9. Build GEO Mime input
############################

prepare_geo_mime <- function(expr_file, clin_file, gene_use, cohort_name) {
  
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
  
  if (length(common_genes) < 2) {
    stop(cohort_name, " has too few matched candidate genes: ", length(common_genes))
  }
  
  x <- as.data.frame(
    t(expr_mat[common_genes, , drop = FALSE]),
    check.names = FALSE
  )
  
  x$ID <- rownames(x)
  
  clin$Sample <- as.character(clin$Sample)
  
  out <- x %>%
    dplyr::inner_join(
      clin[, c("Sample", "time", "status")],
      by = c("ID" = "Sample")
    ) %>%
    dplyr::filter(
      is.finite(time),
      time > 0,
      status %in% c(0, 1)
    )
  
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
    cohort_name, ":", nrow(out), "samples /",
    sum(out$OS == 1), "events /",
    length(common_genes), "matched genes\n"
  )
  
  out
}

gse72094_mime <- prepare_geo_mime(
  gse72094_expr_file,
  gse72094_clin_file,
  candidate_genes,
  "GSE72094"
)

gse68465_mime <- prepare_geo_mime(
  gse68465_expr_file,
  gse68465_clin_file,
  candidate_genes,
  "GSE68465"
)

save_csv(
  gse72094_mime,
  file.path(FIG4_DIR, "Mime1_input_GSE72094.csv")
)

save_csv(
  gse68465_mime,
  file.path(FIG4_DIR, "Mime1_input_GSE68465.csv")
)

list_train_vali_Data <- list(
  Training = train_df_mime,
  Testing = test_df_mime,
  GSE72094 = gse72094_mime,
  GSE68465 = gse68465_mime
)

############################
# 10. Run Mime1
############################

set.seed(SEED)

res.mime <- Mime1::ML.Dev.Prog.Sig(
  train_data = train_df_mime,
  list_train_vali_Data = list_train_vali_Data,
  candidate_genes = candidate_genes,
  mode = "all",
  unicox.filter.for.candi = TRUE,
  unicox_p_cutoff = UNICOX_P,
  nodesize = NODE_SIZE,
  seed = SEED
)

saveRDS(
  res.mime,
  file.path(FIG4_DIR, "Mime1_ML_DDR_res.rds")
)

save_csv(
  data.frame(ObjectNames = names(res.mime)),
  file.path(FIG4_DIR, "Mime1_result_object_names.csv")
)

############################
# 11. C-index summary
############################

cindex_raw <- as.data.frame(res.mime$Cindex.res, check.names = FALSE)

if (all(c("ID", "Model", "Cindex") %in% colnames(cindex_raw))) {
  
  cindex_summary <- cindex_raw %>%
    dplyr::mutate(
      ID = as.character(ID),
      Model = as.character(Model),
      Cindex = as.numeric(Cindex)
    ) %>%
    dplyr::group_by(Model, ID) %>%
    dplyr::summarise(
      Cindex = mean(Cindex, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      id_cols = Model,
      names_from = ID,
      values_from = Cindex
    )
  
} else {
  
  cindex_summary <- cindex_raw
  
  if (!"Model" %in% colnames(cindex_summary)) {
    cindex_summary$Model <- rownames(cindex_summary)
  }
}

dataset_cols <- intersect(
  c("Training", "Testing", "GSE72094", "GSE68465"),
  colnames(cindex_summary)
)

if (length(dataset_cols) < 2) {
  stop("Cannot find dataset C-index columns in cindex_summary.")
}

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
    }
  ) %>%
  dplyr::arrange(
    dplyr::desc(Average),
    dplyr::desc(ExternalMean),
    Train_Test_Gap
  )

save_csv(
  cindex_summary,
  file.path(FIG4_DIR, "Fig4B_all_model_Cindex_summary.csv")
)

############################
# 12. Force final model
############################

if (!TARGET_MODEL %in% cindex_summary$Model) {
  cat("Target model not found in cindex_summary.\n")
  cat("Similar models:\n")
  print(grep("StepCox.*RSF|RSF.*StepCox", cindex_summary$Model, value = TRUE))
  stop("Cannot find target model: ", TARGET_MODEL)
}

best_model <- TARGET_MODEL

selected_model_info <- cindex_summary %>%
  dplyr::filter(Model == best_model)

stable_models <- cindex_summary %>%
  dplyr::filter(
    Training >= 0.60,
    Testing >= 0.58,
    ExternalMean >= 0.58,
    Training <= 0.90,
    Train_Test_Gap <= 0.30
  ) %>%
  dplyr::arrange(
    dplyr::desc(Average),
    dplyr::desc(ExternalMean),
    dplyr::desc(Testing),
    Train_Test_Gap
  )

if (nrow(stable_models) == 0) {
  stable_models <- cindex_summary %>%
    dplyr::filter(
      Testing >= 0.55,
      ExternalMean >= 0.55,
      Train_Test_Gap <= 0.20
    ) %>%
    dplyr::arrange(
      dplyr::desc(Average),
      dplyr::desc(ExternalMean),
      Train_Test_Gap
    )
}

if (nrow(stable_models) == 0) {
  stable_models <- cindex_summary %>%
    dplyr::arrange(
      dplyr::desc(Average),
      Train_Test_Gap
    )
}

stable_models <- dplyr::bind_rows(
  selected_model_info,
  stable_models %>% dplyr::filter(Model != best_model)
) %>%
  dplyr::distinct(Model, .keep_all = TRUE)

save_csv(
  stable_models,
  file.path(FIG4_DIR, "Fig4B_stable_model_candidates.csv")
)

save_csv(
  selected_model_info,
  file.path(FIG4_DIR, "Fig4B_final_selected_model_info.csv")
)

cat("Final model:", best_model, "\n")
print(selected_model_info)

############################
# 13. Extract DPRS risk score
############################

risk_obj <- res.mime$riskscore

if (!best_model %in% names(risk_obj)) {
  cat("Model not found in res.mime$riskscore.\n")
  cat("Similar models:\n")
  print(grep("StepCox.*RSF|RSF.*StepCox", names(risk_obj), value = TRUE))
  stop("Cannot find model in res.mime$riskscore: ", best_model)
}

model_risk <- risk_obj[[best_model]]

extract_mime_score <- function(df, dataset_name) {
  
  df <- as.data.frame(df, check.names = FALSE)
  
  id_col <- grep("^ID$|sample|patient", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(id_col)) id_col <- colnames(df)[1]
  
  time_col <- grep("^OS.time$|^time$|survival", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(time_col)) time_col <- colnames(df)[2]
  
  status_col <- grep("^OS$|^status$|event", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(status_col)) status_col <- colnames(df)[3]
  
  score_col <- grep("^RS$|risk|score", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  
  if (is.na(score_col)) {
    candidate_cols <- setdiff(colnames(df), c(id_col, time_col, status_col))
    numeric_cols <- candidate_cols[sapply(df[, candidate_cols, drop = FALSE], is.numeric)]
    score_col <- numeric_cols[1]
  }
  
  if (is.na(score_col)) {
    stop("Cannot identify risk score column for dataset: ", dataset_name)
  }
  
  data.frame(
    ID = as.character(df[[id_col]]),
    time = as.numeric(df[[time_col]]),
    status = as.numeric(df[[status_col]]),
    ML_DDRscore = as.numeric(df[[score_col]]),
    Dataset = dataset_name,
    Model = best_model,
    stringsAsFactors = FALSE
  )
}

final_score <- dplyr::bind_rows(lapply(names(model_risk), function(ds) {
  extract_mime_score(model_risk[[ds]], ds)
}))

final_score <- final_score %>%
  dplyr::filter(
    is.finite(time),
    time > 0,
    status %in% c(0, 1),
    is.finite(ML_DDRscore)
  )

train_tmp <- final_score %>%
  dplyr::filter(Dataset == "Training")

cox_dir <- survival::coxph(
  survival::Surv(time, status) ~ ML_DDRscore,
  data = train_tmp
)

if (as.numeric(coef(cox_dir)[1]) < 0) {
  final_score$ML_DDRscore <- -final_score$ML_DDRscore
}

final_score <- final_score %>%
  dplyr::group_by(Dataset) %>%
  dplyr::mutate(
    Cutoff = median(ML_DDRscore, na.rm = TRUE),
    RiskGroup = ifelse(ML_DDRscore >= Cutoff, "High", "Low")
  ) %>%
  dplyr::ungroup()

final_score$RiskGroup <- factor(final_score$RiskGroup, levels = c("Low", "High"))
final_score$Dataset <- factor(
  final_score$Dataset,
  levels = c("Training", "Testing", "GSE72094", "GSE68465")
)

save_csv(
  final_score,
  file.path(FIG4_DIR, "Fig4C_ML_DDRscore_all_sets.csv")
)

save_csv(
  final_score,
  file.path(DB_DIR, "ML_DDRscore_table.csv")
)

score_cindex <- final_score %>%
  dplyr::group_by(Dataset) %>%
  dplyr::summarise(
    N = dplyr::n(),
    Events = sum(status == 1),
    Cindex = calc_cindex_risk(time, status, ML_DDRscore),
    .groups = "drop"
  )

save_csv(
  score_cindex,
  file.path(FIG4_DIR, "Fig4D_selected_model_Cindex_by_dataset.csv")
)

############################
# 14. Final signature genes
############################

gene_cols <- setdiff(colnames(train_df_mime), c("ID", "OS.time", "OS"))

gene_map <- data.frame(
  Gene = gene_cols,
  SafeGene = make.names(gene_cols, unique = TRUE),
  stringsAsFactors = FALSE
)

train_safe <- train_df_mime
colnames(train_safe)[match(gene_map$Gene, colnames(train_safe))] <- gene_map$SafeGene

null_fit <- survival::coxph(
  survival::Surv(OS.time, OS) ~ 1,
  data = train_safe
)

upper_formula <- as.formula(
  paste0("~ ", paste(gene_map$SafeGene, collapse = " + "))
)

step_fit <- step(
  null_fit,
  scope = list(lower = ~1, upper = upper_formula),
  direction = "forward",
  trace = 0
)

selected_safe_genes <- names(coef(step_fit))

if (length(selected_safe_genes) == 0) {
  stop("StepCox[forward] selected zero genes.")
}

signature_genes <- gene_map %>%
  dplyr::filter(SafeGene %in% selected_safe_genes) %>%
  dplyr::mutate(
    StepCoxCoef = as.numeric(coef(step_fit)[SafeGene])
  ) %>%
  dplyr::select(Gene, SafeGene, StepCoxCoef) %>%
  dplyr::left_join(
    cox_table %>%
      dplyr::select(Gene, CoxCoef, HR, CoxP, CoxFDR),
    by = "Gene"
  ) %>%
  dplyr::mutate(
    Model = best_model,
    GeneRole = ifelse(StepCoxCoef >= 0, "Risk", "Protective"),
    SignatureType = "StepCox_forward_selected_RSF_input",
    Note = "Final DPRS feature genes selected by StepCox[forward] and used as RSF input"
  ) %>%
  dplyr::arrange(dplyr::desc(abs(StepCoxCoef)))

save_csv(
  signature_genes,
  file.path(FIG4_DIR, "Fig4C_final_signature_genes.csv")
)

save_csv(
  signature_genes,
  file.path(DB_DIR, "ML_DDR_signature_table.csv")
)

############################
# 15. Output manifest
############################

output_manifest <- data.frame(
  File = c(
    "Fig4A_TCGA_score_survival_table.csv",
    "Fig4A_train_test_split_info.csv",
    "Fig4A_training_only_DEG.csv",
    "Fig4A_DEG_genes_common_with_GEO.csv",
    "Fig4A_training_only_unicox.csv",
    "Fig4A_candidate_genes_train_only.csv",
    "Mime1_input_Training.csv",
    "Mime1_input_Testing.csv",
    "Mime1_input_GSE72094.csv",
    "Mime1_input_GSE68465.csv",
    "Mime1_ML_DDR_res.rds",
    "Mime1_result_object_names.csv",
    "Fig4B_all_model_Cindex_summary.csv",
    "Fig4B_stable_model_candidates.csv",
    "Fig4B_final_selected_model_info.csv",
    "Fig4C_ML_DDRscore_all_sets.csv",
    "Fig4C_final_signature_genes.csv",
    "Fig4D_selected_model_Cindex_by_dataset.csv"
  ),
  Description = c(
    "TCGA samples with MO-DDRscore group and survival information",
    "Training/testing split information",
    "Training-only DEG result between MO-DDRscore high and low groups",
    "DEG genes retained in both GEO validation cohorts",
    "Training-only univariate Cox result",
    "Candidate genes used for ML modeling",
    "Mime1 input table for training cohort",
    "Mime1 input table for testing cohort",
    "Mime1 input table for GSE72094",
    "Mime1 input table for GSE68465",
    "Raw Mime1 result object",
    "Names of objects stored in Mime1 result",
    "C-index summary for all ML models",
    "Stable model candidates after filtering",
    "Manually selected final DPRS model information",
    "Final DPRS score table across all cohorts",
    "Final DPRS feature genes selected by StepCox[forward]",
    "C-index of selected final DPRS model by cohort"
  ),
  Path = c(
    file.path(FIG4_DIR, "Fig4A_TCGA_score_survival_table.csv"),
    file.path(FIG4_DIR, "Fig4A_train_test_split_info.csv"),
    file.path(FIG4_DIR, "Fig4A_training_only_DEG.csv"),
    file.path(FIG4_DIR, "Fig4A_DEG_genes_common_with_GEO.csv"),
    file.path(FIG4_DIR, "Fig4A_training_only_unicox.csv"),
    file.path(FIG4_DIR, "Fig4A_candidate_genes_train_only.csv"),
    file.path(FIG4_DIR, "Mime1_input_Training.csv"),
    file.path(FIG4_DIR, "Mime1_input_Testing.csv"),
    file.path(FIG4_DIR, "Mime1_input_GSE72094.csv"),
    file.path(FIG4_DIR, "Mime1_input_GSE68465.csv"),
    file.path(FIG4_DIR, "Mime1_ML_DDR_res.rds"),
    file.path(FIG4_DIR, "Mime1_result_object_names.csv"),
    file.path(FIG4_DIR, "Fig4B_all_model_Cindex_summary.csv"),
    file.path(FIG4_DIR, "Fig4B_stable_model_candidates.csv"),
    file.path(FIG4_DIR, "Fig4B_final_selected_model_info.csv"),
    file.path(FIG4_DIR, "Fig4C_ML_DDRscore_all_sets.csv"),
    file.path(FIG4_DIR, "Fig4C_final_signature_genes.csv"),
    file.path(FIG4_DIR, "Fig4D_selected_model_Cindex_by_dataset.csv")
  ),
  stringsAsFactors = FALSE
)

save_csv(
  output_manifest,
  file.path(FIG4_DIR, "ML_DPRS_modeling_output_manifest.csv")
)

sink(file.path(FIG4_DIR, "ML_DPRS_modeling_session_info.txt"))
cat("ML-DPRS modeling finished.\n")
cat("Project directory:", PROJECT_DIR, "\n")
cat("Output directory:", FIG4_DIR, "\n")
cat("Final model:", best_model, "\n\n")
cat("Parameters:\n")
cat("TRAIN_RATIO =", TRAIN_RATIO, "\n")
cat("DEG_ADJ_P =", DEG_ADJ_P, "\n")
cat("DEG_LOGFC =", DEG_LOGFC, "\n")
cat("UNICOX_P =", UNICOX_P, "\n")
cat("MAX_ML_GENES =", MAX_ML_GENES, "\n")
cat("NODE_SIZE =", NODE_SIZE, "\n")
cat("SEED =", SEED, "\n\n")
cat("Selected model info:\n")
print(selected_model_info)
cat("\nSelected model C-index by dataset:\n")
print(score_cindex)
cat("\nFinal signature genes:\n")
print(signature_genes)
cat("\nSession info:\n")
print(sessionInfo())
sink()

cat("\nML-DPRS modeling finished.\n")
cat("Final model:", best_model, "\n")
cat("Output directory:", FIG4_DIR, "\n")
cat("\nMain outputs:\n")
cat(" - Fig4B_all_model_Cindex_summary.csv\n")
cat(" - Fig4B_final_selected_model_info.csv\n")
cat(" - Fig4C_ML_DDRscore_all_sets.csv\n")
cat(" - Fig4C_final_signature_genes.csv\n")
cat(" - Fig4D_selected_model_Cindex_by_dataset.csv\n")
cat(" - ML_DPRS_modeling_output_manifest.csv\n")