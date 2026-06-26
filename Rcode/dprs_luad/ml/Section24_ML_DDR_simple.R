############################################################
# Section 24 simple version
# ML-DDR prognostic signature by Mime1
# Key rule:
#   Train/Test split first
#   Training-only DEG
#   Training-only Cox filtering
#   Testing/GEO only validation
############################################################

set.seed(20260513)

############################
# 0. Parameters
############################

PROJECT_DIR <- "D:/R_workspace/评分/AD_DDR_project"
DATA_DIR <- file.path(PROJECT_DIR, "00_data")
DB_DIR <- file.path(PROJECT_DIR, "05_database_tables")
FIG4_DIR <- file.path(PROJECT_DIR, "04")
dir.create(FIG4_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DB_DIR, recursive = TRUE, showWarnings = FALSE)

TRAIN_RATIO <- 0.80
DEG_ADJ_P <- 0.05
DEG_LOGFC <- 1.20
UNICOX_P <- 0.03
MAX_ML_GENES <- 60
NODE_SIZE <- 20
SEED <- 20260513

############################
# 1. Packages
############################

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(limma)
library(Mime1)

save_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file)
}

safe_ggsave <- function(file, p, w = 7, h = 5) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file, p, width = w, height = h, device = "pdf", useDingbats = FALSE)
}

clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub(" ///.*$| //.*$|;.*$|,.*$|\\s*\\([^\\)]*\\)$", "", x)
  x[x %in% c("", "NA", "---", "NULL", "N/A", "NAN")] <- NA
  x
}

patient_id <- function(x) substr(gsub("\\.", "-", as.character(x)), 1, 12)
sample_type <- function(x) substr(gsub("\\.", "-", as.character(x)), 14, 15)

zscore_vector <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

calc_cindex <- function(time, status, score) {
  df <- data.frame(time = as.numeric(time), status = as.numeric(status), score = as.numeric(score))
  df <- df[is.finite(df$time) & is.finite(df$status) & is.finite(df$score), ]
  if (nrow(df) < 30 || sum(df$status == 1) < 5) return(NA_real_)
  as.numeric(survival::concordance(Surv(time, status) ~ score, data = df)$concordance)
}

############################
# 2. Required objects from previous pipeline
############################
# Must exist:
#   tcga_expr: gene x sample TPM matrix
#   tumor_samples: TCGA tumor sample IDs
#   tcga_clin: Patient, time, status, age, gender, stage
#   mo_score: Sample, Patient, MO_DDRscore_group, MO_DDRscore_raw

stopifnot(exists("tcga_expr"))
stopifnot(exists("tumor_samples"))
stopifnot(exists("tcga_clin"))
stopifnot(exists("mo_score"))

score_df <- mo_score %>%
  filter(Sample %in% tumor_samples,
         !is.na(MO_DDRscore_group),
         is.finite(MO_DDRscore_raw)) %>%
  arrange(Patient, SampleType) %>%
  distinct(Patient, .keep_all = TRUE) %>%
  inner_join(tcga_clin[, c("Patient", "time", "status", "age", "gender", "stage")],
             by = "Patient") %>%
  filter(is.finite(time), time > 0, status %in% c(0, 1))

############################
# 3. Split TCGA first
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

cat("Training:", length(train_ids), "patients /",
    sum(score_df$status[score_df$Patient %in% train_ids] == 1), "events\n")
cat("Testing:", length(test_ids), "patients /",
    sum(score_df$status[score_df$Patient %in% test_ids] == 1), "events\n")

save_csv(score_df, file.path(FIG4_DIR, "Fig4A_TCGA_score_survival_table.csv"))

############################
# 4. Training-only DEG
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

save_csv(train_deg, file.path(FIG4_DIR, "Fig4A_training_only_DEG.csv"))

deg_genes <- train_deg %>%
  filter(is.finite(FDR), FDR < DEG_ADJ_P, abs(logFC) > DEG_LOGFC) %>%
  arrange(FDR, desc(abs(logFC))) %>%
  pull(Gene) %>%
  unique()

if (length(deg_genes) < 30) {
  deg_genes <- train_deg %>%
    arrange(FDR, desc(abs(logFC))) %>%
    slice_head(n = 500) %>%
    pull(Gene) %>%
    unique()
}

############################
# 5. GEO files and common genes
############################

gse72094_expr_file <- file.path(DATA_DIR, "GSE72094_expression.csv")
gse72094_clin_file <- file.path(DATA_DIR, "GSE72094_clinical.csv")
gse68465_expr_file <- file.path(DATA_DIR, "GSE68465_expression.csv")
gse68465_clin_file <- file.path(DATA_DIR, "GSE68465_clinical.csv")

read_geo_genes <- function(expr_file) {
  x <- data.table::fread(expr_file, data.table = FALSE, check.names = FALSE)
  unique(na.omit(clean_gene(x[[1]])))
}

geo_common_genes <- Reduce(
  intersect,
  list(
    rownames(tcga_expr),
    read_geo_genes(gse72094_expr_file),
    read_geo_genes(gse68465_expr_file)
  )
)

deg_genes <- intersect(deg_genes, geo_common_genes)

############################
# 6. Training-only Cox filtering
############################

train_expr_cox <- as.data.frame(t(log2(tcga_expr[deg_genes, train_samples, drop = FALSE] + 1)),
                                check.names = FALSE)
train_expr_cox$Sample <- rownames(train_expr_cox)
train_expr_cox$Patient <- patient_id(train_expr_cox$Sample)

train_surv <- train_expr_cox %>%
  inner_join(tcga_clin[, c("Patient", "time", "status")], by = "Patient")

cox_table <- bind_rows(lapply(deg_genes, function(g) {
  df <- train_surv[, c("time", "status", g)]
  colnames(df) <- c("time", "status", "expr")
  fit <- tryCatch(coxph(Surv(time, status) ~ expr, data = df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  sm <- summary(fit)
  data.frame(
    Gene = g,
    CoxCoef = sm$coefficients[1, "coef"],
    HR = sm$coefficients[1, "exp(coef)"],
    CoxP = sm$coefficients[1, "Pr(>|z|)"]
  )
}))

cox_table$CoxFDR <- p.adjust(cox_table$CoxP, method = "BH")
cox_table <- cox_table %>% arrange(CoxP)

save_csv(cox_table, file.path(FIG4_DIR, "Fig4A_training_only_unicox.csv"))

candidate_genes <- cox_table %>%
  filter(CoxP < UNICOX_P) %>%
  pull(Gene) %>%
  unique()

if (length(candidate_genes) < 10) {
  candidate_genes <- head(cox_table$Gene, 50)
}

candidate_genes <- head(unique(candidate_genes), MAX_ML_GENES)

save_csv(data.frame(Gene = candidate_genes),
         file.path(FIG4_DIR, "Fig4A_candidate_genes_train_only.csv"))

cat("Candidate genes:", length(candidate_genes), "\n")

############################
# 7. Build TCGA Mime input
############################

make_tcga_mime <- function(samples, gene_use) {
  x <- as.data.frame(t(log2(tcga_expr[gene_use, samples, drop = FALSE] + 1)),
                     check.names = FALSE)
  x$Sample <- rownames(x)
  x$Patient <- patient_id(x$Sample)
  x <- x %>%
    arrange(Patient, Sample) %>%
    distinct(Patient, .keep_all = TRUE) %>%
    inner_join(tcga_clin[, c("Patient", "time", "status")], by = "Patient") %>%
    filter(is.finite(time), time > 0, status %in% c(0, 1))

  out <- x[, c("Patient", "time", "status", gene_use)]
  colnames(out)[1:3] <- c("ID", "OS.time", "OS")
  for (g in gene_use) out[[g]] <- zscore_vector(out[[g]])
  out
}

mime_all <- make_tcga_mime(score_df$Sample, candidate_genes)
train_df_mime <- mime_all[mime_all$ID %in% train_ids, ]
test_df_mime  <- mime_all[mime_all$ID %in% test_ids, ]

save_csv(train_df_mime, file.path(FIG4_DIR, "Mime1_input_Training.csv"))
save_csv(test_df_mime, file.path(FIG4_DIR, "Mime1_input_Testing.csv"))

############################
# 8. Build GEO Mime input
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

  x <- as.data.frame(t(expr_mat[common_genes, , drop = FALSE]), check.names = FALSE)
  x$ID <- rownames(x)

  clin$Sample <- as.character(clin$Sample)

  out <- x %>%
    inner_join(clin[, c("Sample", "time", "status")], by = c("ID" = "Sample")) %>%
    filter(is.finite(time), time > 0, status %in% c(0, 1))

  missing_genes <- setdiff(gene_use, colnames(out))
  for (g in missing_genes) out[[g]] <- 0

  out <- out[, c("ID", "time", "status", gene_use)]
  colnames(out)[1:3] <- c("ID", "OS.time", "OS")

  for (g in gene_use) out[[g]] <- zscore_vector(out[[g]])

  cat(cohort_name, ":", nrow(out), "samples /", sum(out$OS == 1), "events /",
      length(common_genes), "genes\n")

  out
}

gse72094_mime <- prepare_geo_mime(
  gse72094_expr_file, gse72094_clin_file, candidate_genes, "GSE72094"
)

gse68465_mime <- prepare_geo_mime(
  gse68465_expr_file, gse68465_clin_file, candidate_genes, "GSE68465"
)

list_train_vali_Data <- list(
  Training = train_df_mime,
  Testing = test_df_mime,
  GSE72094 = gse72094_mime,
  GSE68465 = gse68465_mime
)

############################
# 9. Run Mime1
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

saveRDS(res.mime, file.path(FIG4_DIR, "Mime1_ML_DDR_res.rds"))
save_csv(data.frame(ObjectNames = names(res.mime)),
         file.path(FIG4_DIR, "Mime1_result_object_names.csv"))

############################
# 10. C-index table and force final model
############################

# ------------------------------------------------------------
# 10.0 Safety: if running from Section 10 directly, reload needed files
# ------------------------------------------------------------

if (!exists("FIG4_DIR")) {
  PROJECT_DIR <- "D:/R_workspace/评分/AD_DDR_project"
  FIG4_DIR <- file.path(PROJECT_DIR, "04")
  DATA_DIR <- file.path(PROJECT_DIR, "00_data")
  DB_DIR <- file.path(PROJECT_DIR, "05_database_tables")
}

dir.create(FIG4_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DB_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("save_csv")) {
  save_csv <- function(x, file) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(x, file)
  }
}

if (!exists("safe_ggsave")) {
  safe_ggsave <- function(file, p, w = 7, h = 5) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(file, p, width = w, height = h, device = "pdf", useDingbats = FALSE)
  }
}

if (!exists("calc_cindex_risk")) {
  calc_cindex_risk <- function(time, status, score) {
    df <- data.frame(
      time = as.numeric(time),
      status = as.numeric(status),
      score = as.numeric(score)
    )
    df <- df[is.finite(df$time) & is.finite(df$status) & is.finite(df$score), ]
    if (nrow(df) < 30 || sum(df$status == 1) < 5) return(NA_real_)
    
    # score 越高 = 风险越高，所以 reverse = TRUE
    as.numeric(
      survival::concordance(
        survival::Surv(time, status) ~ score,
        data = df,
        reverse = TRUE
      )$concordance
    )
  }
}

if (!exists("get_time_points")) {
  get_time_points <- function(time) {
    if (median(time, na.rm = TRUE) > 100) {
      c(365, 1095, 1825)
    } else {
      c(1, 3, 5)
    }
  }
}

# 如果不是从前面连续运行，而是只从第10部分开始跑，则读取 rds
if (!exists("res.mime")) {
  rds_file <- file.path(FIG4_DIR, "Mime1_ML_DDR_res.rds")
  if (!file.exists(rds_file)) {
    stop("Cannot find Mime1 result rds: ", rds_file)
  }
  res.mime <- readRDS(rds_file)
}

# 如果当前环境没有 train_df_mime / cox_table / candidate_genes，就从文件读
if (!exists("train_df_mime")) {
  train_file <- file.path(FIG4_DIR, "Mime1_input_Training.csv")
  if (!file.exists(train_file)) stop("Cannot find: ", train_file)
  train_df_mime <- data.table::fread(train_file, data.table = FALSE, check.names = FALSE)
}

if (!exists("cox_table")) {
  unicox_file <- file.path(FIG4_DIR, "Fig4A_training_only_unicox.csv")
  if (!file.exists(unicox_file)) stop("Cannot find: ", unicox_file)
  cox_table <- data.table::fread(unicox_file, data.table = FALSE, check.names = FALSE)
}

if (!exists("candidate_genes")) {
  candidate_file <- file.path(FIG4_DIR, "Fig4A_candidate_genes_train_only.csv")
  if (!file.exists(candidate_file)) stop("Cannot find: ", candidate_file)
  candidate_genes <- data.table::fread(candidate_file, data.table = FALSE)$Gene
}

# ------------------------------------------------------------
# 10.1 Build C-index summary from Mime1 result
# ------------------------------------------------------------

cindex_raw <- as.data.frame(res.mime$Cindex.res, check.names = FALSE)

if (all(c("ID", "Model", "Cindex") %in% colnames(cindex_raw))) {
  cindex_summary <- cindex_raw %>%
    dplyr::mutate(
      ID = as.character(ID),
      Model = as.character(Model),
      Cindex = as.numeric(Cindex)
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
    Average = rowMeans(dplyr::across(dplyr::all_of(dataset_cols)), na.rm = TRUE),
    ExternalMean = rowMeans(
      dplyr::across(dplyr::all_of(intersect(c("GSE72094", "GSE68465"), dataset_cols))),
      na.rm = TRUE
    ),
    Train_Test_Gap = if (all(c("Training", "Testing") %in% dataset_cols)) {
      abs(Training - Testing)
    } else {
      NA_real_
    }
  ) %>%
  dplyr::arrange(dplyr::desc(Average), dplyr::desc(ExternalMean))

save_csv(
  cindex_summary,
  file.path(FIG4_DIR, "Fig4B_all_model_Cindex_summary.csv")
)

# ------------------------------------------------------------
# 10.2 Stable model candidates
# ------------------------------------------------------------

stable_models <- cindex_summary %>%
  dplyr::filter(
    Training >= 0.60,
    Testing >= 0.58,
    ExternalMean >= 0.58,
    Training <= 0.90,
    Train_Test_Gap <= 0.3
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
    dplyr::arrange(dplyr::desc(Average), Train_Test_Gap)
}

# ------------------------------------------------------------
# 10.3 Force DPRS final model
# ------------------------------------------------------------

target_model <- "StepCox[forward] + RSF"

if (!target_model %in% cindex_summary$Model) {
  cat("Target model not found in cindex_summary.\n")
  cat("Similar models are:\n")
  print(grep("StepCox.*RSF|RSF.*StepCox", cindex_summary$Model, value = TRUE))
  stop("Cannot find target model: ", target_model)
}

best_model <- target_model

selected_model_info <- cindex_summary %>%
  dplyr::filter(Model == best_model)

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

cat("Best model manually set to:", best_model, "\n")
print(selected_model_info)


############################
# 11. C-index heatmap
############################

ci_show <- stable_models %>%
  dplyr::slice_head(n = 40)

ci_long <- ci_show %>%
  dplyr::select(Model, dplyr::all_of(dataset_cols), Average, ExternalMean) %>%
  tidyr::pivot_longer(
    cols = -Model,
    names_to = "Dataset",
    values_to = "Cindex"
  )

ci_long$Model <- factor(ci_long$Model, levels = rev(ci_show$Model))
ci_long$Dataset <- factor(
  ci_long$Dataset,
  levels = c(dataset_cols, "Average", "ExternalMean")
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
  labs(title = "Model performance across cohorts")

safe_ggsave(
  file.path(FIG4_DIR, "Fig4B_Cindex_heatmap.pdf"),
  p_ci,
  7.2,
  9
)


############################
# 12. Extract risk scores from Mime1
############################

risk_obj <- res.mime$riskscore

if (!best_model %in% names(risk_obj)) {
  cat("best_model not found in res.mime$riskscore.\n")
  cat("Available StepCox + RSF-like models:\n")
  print(grep("StepCox.*RSF|RSF.*StepCox", names(risk_obj), value = TRUE))
  stop("best_model not found in res.mime$riskscore: ", best_model)
}

cat("Extracting DPRS risk scores from model:", best_model, "\n")

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

# 统一方向：Training 中分数越高 = 风险越高
train_tmp <- final_score %>%
  dplyr::filter(Dataset == "Training")

cox_dir <- survival::coxph(
  survival::Surv(time, status) ~ ML_DDRscore,
  data = train_tmp
)

if (as.numeric(coef(cox_dir)[1]) < 0) {
  final_score$ML_DDRscore <- -final_score$ML_DDRscore
}

# 每个队列内部按中位数分组
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

cat("DPRS score table saved.\n")
print(table(final_score$Dataset, final_score$RiskGroup))


############################
# 13. Final DPRS signature gene table
# For StepCox[forward] + RSF:
# final feature genes = genes selected by StepCox[forward],
# then used as input features for RSF.
############################

cat("Generating final signature genes for:", best_model, "\n")

gene_cols <- setdiff(colnames(train_df_mime), c("ID", "OS.time", "OS"))

gene_map <- data.frame(
  Gene = gene_cols,
  SafeGene = make.names(gene_cols, unique = TRUE),
  stringsAsFactors = FALSE
)

train_safe <- train_df_mime
colnames(train_safe)[match(gene_map$Gene, colnames(train_safe))] <- gene_map$SafeGene

# StepCox[forward]：从空模型开始，向前选择
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
  stop("StepCox[forward] selected zero genes. Please check candidate_genes or training data.")
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

cat("Final DPRS signature genes saved:\n")
print(signature_genes)
cat("Number of final signature genes:", nrow(signature_genes), "\n")

# Final signature gene lollipop
p_gene <- ggplot(
  signature_genes,
  aes(
    x = reorder(Gene, StepCoxCoef),
    y = StepCoxCoef,
    color = GeneRole
  )
) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_segment(
    aes(xend = Gene, y = 0, yend = StepCoxCoef),
    color = "grey60",
    linewidth = 0.85
  ) +
  geom_point(size = 2.8) +
  scale_color_manual(values = c("Risk" = "#C53030", "Protective" = "#2B6CB0")) +
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = NULL,
    y = "StepCox coefficient",
    title = paste0("Final DPRS signature genes: ", best_model)
  )

safe_ggsave(
  file.path(FIG4_DIR, "Fig4C_signature_gene_lollipop.pdf"),
  p_gene,
  6.5,
  max(4.5, 0.28 * nrow(signature_genes) + 2)
)


############################
# 14. Risk plot, KM, timeROC
############################

risk_rank <- final_score %>%
  dplyr::filter(is.finite(ML_DDRscore)) %>%
  dplyr::arrange(Dataset, ML_DDRscore) %>%
  dplyr::group_by(Dataset) %>%
  dplyr::mutate(Rank = dplyr::row_number()) %>%
  dplyr::ungroup()

p_risk <- ggplot(risk_rank, aes(Rank, ML_DDRscore, color = RiskGroup)) +
  geom_point(size = 0.8, alpha = 0.9) +
  facet_wrap(~Dataset, scales = "free_x", nrow = 1) +
  scale_color_manual(values = c(Low = "#2B6CB0", High = "#C53030")) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    legend.title = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = "Patient rank",
    y = "DPRS",
    title = paste0("Risk score distribution: ", best_model)
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
    axis.title = element_text(color = "black", face = "bold"),
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
    Cindex = calc_cindex_risk(time, status, ML_DDRscore),
    .groups = "drop"
  )

save_csv(
  score_cindex,
  file.path(FIG4_DIR, "Fig4D_selected_model_Cindex_by_dataset.csv")
)

cat("Selected model C-index by dataset:\n")
print(score_cindex)

# KM and timeROC
for (ds in unique(as.character(final_score$Dataset))) {
  
  df <- final_score %>%
    dplyr::filter(Dataset == ds)
  
  if (nrow(df) < 40 || sum(df$status == 1) < 8) {
    message("Skip KM/timeROC for ", ds, ": too few samples/events.")
    next
  }
  
  if (requireNamespace("survminer", quietly = TRUE)) {
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
  }
  
  if (requireNamespace("timeROC", quietly = TRUE)) {
    
    times_use <- get_time_points(df$time)
    
    roc <- timeROC::timeROC(
      T = df$time,
      delta = df$status,
      marker = df$ML_DDRscore,
      cause = 1,
      weighting = "marginal",
      times = times_use,
      ROC = TRUE
    )
    
    auc_df <- data.frame(
      Dataset = ds,
      Model = best_model,
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
}

cat("\nSection 24 replacement from Section 10 finished.\n")
cat("Final DPRS model:", best_model, "\n")
cat("Output directory:", FIG4_DIR, "\n")
cat("Main outputs:\n")
cat(" - Fig4B_final_selected_model_info.csv\n")
cat(" - Fig4C_ML_DDRscore_all_sets.csv\n")
cat(" - Fig4C_final_signature_genes.csv\n")
cat(" - Fig4C_signature_gene_lollipop.pdf\n")
cat(" - Fig4D_selected_model_Cindex_by_dataset.csv\n")
cat(" - Fig4E_*_KM.pdf\n")
cat(" - Fig4F_*_timeROC.pdf\n")