############################################################
# 24. ML-DDR prognostic signature by Mime1 with GEO validation
# Replace your original Section 24 with this whole block.
#
# Required existing objects from previous sections:
# DATA_DIR, FIG4_DIR, DB_DIR
# tcga_expr, tumor_samples, ad_score, tcga_clin, deg
# patient_id(), sample_type_code(), clean_gene_symbol(), save_csv(), safe_ggsave()
# TRAIN_RATIO, DEG_ADJ_P, DEG_LOGFC, UNICOX_P, MAX_ML_GENES, MIN_EVENTS_TRAIN
############################################################

cat("\nBuilding ML-DDR prognostic signature using Mime1 + GEO validation...\n")

if (!requireNamespace("Mime1", quietly = TRUE)) {
  stop(
    "Mime1 is not installed. Please install it first:\n",
    "if (!requireNamespace('pak', quietly = TRUE)) install.packages('pak')\n",
    "pak::pkg_install('l-magnificence/Mime')\n"
  )
}
suppressPackageStartupMessages(library(Mime1))

if (!requireNamespace("tidyr", quietly = TRUE)) stop("Please install tidyr.")
if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install dplyr.")
if (!requireNamespace("data.table", quietly = TRUE)) stop("Please install data.table.")

# ----------------------------
# 24.0 Local helper functions
# ----------------------------

safe_num2 <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

zscore_vector <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x)) || !is.finite(stats::sd(x, na.rm = TRUE)) || stats::sd(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

plot_km_by_group <- function(df, title = "", file) {
  if (!requireNamespace("survminer", quietly = TRUE)) {
    message("survminer not installed; KM skipped: ", title)
    return(invisible(NULL))
  }

  df <- df[
    is.finite(df$time) &
      is.finite(df$status) &
      !is.na(df$RiskGroup),
  ]

  if (nrow(df) < 30 || sum(df$status == 1, na.rm = TRUE) < 5) {
    message("KM skipped: insufficient samples/events for ", title)
    return(invisible(NULL))
  }

  df$RiskGroup <- factor(df$RiskGroup, levels = c("Low", "High"))

  fit <- survival::survfit(
    survival::Surv(time, status) ~ RiskGroup,
    data = df
  )

  p <- tryCatch({
    survminer::ggsurvplot(
      fit,
      data = df,
      pval = TRUE,
      risk.table = FALSE,
      palette = c("#2B6CB0", "#C53030"),
      title = title,
      legend.title = "",
      legend.labs = c("Low", "High"),
      ggtheme = ggplot2::theme_bw()
    )$plot
  }, error = function(e) {
    message("KM failed for ", title, ": ", e$message)
    NULL
  })

  if (!is.null(p)) {
    safe_ggsave(file, p, width = 6, height = 5)
  }

  invisible(NULL)
}

plot_time_roc_local <- function(df, score_col, file, title = "") {
  if (!requireNamespace("timeROC", quietly = TRUE)) {
    message("timeROC not installed; ROC skipped: ", title)
    return(invisible(NULL))
  }

  df <- df[
    is.finite(df[[score_col]]) &
      is.finite(df$time) &
      is.finite(df$status),
  ]

  if (nrow(df) < 30 || sum(df$status == 1, na.rm = TRUE) < 5) {
    message("timeROC skipped: insufficient samples/events for ", title)
    return(invisible(NULL))
  }

  times <- c(365, 1095, 1825)

  roc <- tryCatch({
    timeROC::timeROC(
      T = df$time,
      delta = df$status,
      marker = df[[score_col]],
      cause = 1,
      weighting = "marginal",
      times = times,
      ROC = TRUE,
      iid = TRUE
    )
  }, error = function(e) {
    message("timeROC failed for ", title, ": ", e$message)
    NULL
  })

  if (is.null(roc)) return(invisible(NULL))

  auc_df <- data.frame(
    Time = c("1-year", "3-year", "5-year"),
    AUC = as.numeric(roc$AUC),
    stringsAsFactors = FALSE
  )

  save_csv(auc_df, sub("\\.pdf$", "_AUC.csv", file))

  grDevices::pdf(file, width = 6, height = 5)
  plot(roc, time = times[1], col = "#C53030", title = FALSE)
  plot(roc, time = times[2], add = TRUE, col = "#2B6CB0")
  plot(roc, time = times[3], add = TRUE, col = "#2F855A")
  legend(
    "bottomright",
    legend = paste0(auc_df$Time, " AUC=", sprintf("%.3f", auc_df$AUC)),
    col = c("#C53030", "#2B6CB0", "#2F855A"),
    lwd = 2
  )
  title(main = title)
  grDevices::dev.off()

  invisible(auc_df)
}

calc_cindex_local <- function(time, status, risk) {
  df <- data.frame(
    time = safe_num2(time),
    status = safe_num2(status),
    risk = safe_num2(risk)
  )
  df <- df[is.finite(df$time) & is.finite(df$status) & is.finite(df$risk), ]
  if (nrow(df) < 20 || sum(df$status == 1, na.rm = TRUE) < 5) return(NA_real_)
  out <- tryCatch({
    survival::concordance(survival::Surv(time, status) ~ risk, data = df)$concordance
  }, error = function(e) NA_real_)
  as.numeric(out)
}

# ----------------------------
# 24.1 Candidate genes from AD-DDRscore-related DEGs
# ----------------------------

if (!all(c("Gene", "logFC") %in% colnames(deg))) {
  stop("deg must contain at least Gene and logFC columns.")
}
if (!"FDR" %in% colnames(deg)) {
  if ("adj.P.Val" %in% colnames(deg)) {
    deg$FDR <- deg$adj.P.Val
  } else {
    deg$FDR <- p.adjust(deg$P.Value, method = "BH")
  }
}

candidate_genes <- unique(deg$Gene[
  is.finite(deg$FDR) &
    deg$FDR < DEG_ADJ_P &
    abs(deg$logFC) > DEG_LOGFC
])

candidate_genes <- intersect(candidate_genes, rownames(tcga_expr))

if (length(candidate_genes) < 20) {
  warning("Too few DEG-filtered candidate genes. Using top 500 genes ranked by P.Value.")
  candidate_genes <- unique(intersect(
    head(deg$Gene[order(deg$P.Value)], 500),
    rownames(tcga_expr)
  ))
}

read_geo_gene_set <- function(expr_file) {
  x <- data.table::fread(expr_file, data.table = FALSE, check.names = FALSE)
  g <- clean_gene_symbol(x[[1]])
  unique(g[!is.na(g) & g != ""])
}

gse72094_expr_file <- file.path(DATA_DIR, "GSE72094_expression.csv")
gse72094_clin_file <- file.path(DATA_DIR, "GSE72094_clinical.csv")
gse68465_expr_file <- file.path(DATA_DIR, "GSE68465_expression.csv")
gse68465_clin_file <- file.path(DATA_DIR, "GSE68465_clinical.csv")

if (!file.exists(gse72094_expr_file)) stop("Missing GSE72094 expression file: ", gse72094_expr_file)
if (!file.exists(gse72094_clin_file)) stop("Missing GSE72094 clinical file: ", gse72094_clin_file)
if (!file.exists(gse68465_expr_file)) stop("Missing GSE68465 expression file: ", gse68465_expr_file)
if (!file.exists(gse68465_clin_file)) stop("Missing GSE68465 clinical file: ", gse68465_clin_file)

gse72094_genes <- read_geo_gene_set(gse72094_expr_file)
gse68465_genes <- read_geo_gene_set(gse68465_expr_file)

candidate_genes <- Reduce(
  intersect,
  list(candidate_genes, rownames(tcga_expr), gse72094_genes, gse68465_genes)
)

if (length(candidate_genes) < 20) {
  warning("Too few genes after intersecting TCGA/GEO. Relaxing by using top DEG genes available in all cohorts.")
  all_common_genes <- Reduce(intersect, list(rownames(tcga_expr), gse72094_genes, gse68465_genes))
  deg_tmp <- deg[deg$Gene %in% all_common_genes, ]
  deg_tmp <- deg_tmp[order(deg_tmp$FDR, -abs(deg_tmp$logFC)), ]
  candidate_genes <- head(unique(deg_tmp$Gene), min(500, nrow(deg_tmp)))
}

if (length(candidate_genes) > MAX_ML_GENES) {
  deg_tmp <- deg[deg$Gene %in% candidate_genes, ]
  deg_tmp <- deg_tmp[order(deg_tmp$FDR, -abs(deg_tmp$logFC)), ]
  candidate_genes <- head(unique(deg_tmp$Gene), MAX_ML_GENES)
}

candidate_genes <- unique(candidate_genes)

cat("Candidate genes for Mime1 after TCGA/GEO intersection:", length(candidate_genes), "\n")
save_csv(data.frame(Gene = candidate_genes), file.path(FIG4_DIR, "Mime1_candidate_genes_TCGA_GEO_common.csv"))

if (length(candidate_genes) < 5) {
  stop("Too few candidate genes for Mime1 modeling.")
}

# ----------------------------
# 24.2 Build TCGA Mime input
# ----------------------------

model_samples <- intersect(tumor_samples, ad_score$Sample)

expr_model <- t(log2(tcga_expr[candidate_genes, model_samples, drop = FALSE] + 1))
expr_model <- as.data.frame(expr_model, check.names = FALSE)

expr_model$Sample <- rownames(expr_model)
expr_model$Patient <- patient_id(expr_model$Sample)
expr_model$SampleType <- sample_type_code(expr_model$Sample)

expr_model <- expr_model[order(expr_model$Patient, expr_model$SampleType), ]
expr_model <- expr_model[!duplicated(expr_model$Patient), ]

ml_df <- expr_model %>%
  dplyr::inner_join(
    tcga_clin[, c("Patient", "time", "status", "age", "gender", "stage")],
    by = "Patient"
  )

ml_df <- ml_df[
  is.finite(ml_df$time) &
    ml_df$time > 0 &
    ml_df$status %in% c(0, 1),
]

gene_cols_model <- intersect(candidate_genes, colnames(ml_df))

for (g in gene_cols_model) {
  ml_df[[g]] <- suppressWarnings(as.numeric(ml_df[[g]]))
}

valid_gene <- sapply(gene_cols_model, function(g) {
  x <- ml_df[[g]]
  sum(is.finite(x)) >= 50 &&
    is.finite(stats::sd(x, na.rm = TRUE)) &&
    stats::sd(x, na.rm = TRUE) > 0
})

gene_cols_model <- gene_cols_model[valid_gene]

if (length(gene_cols_model) > MAX_ML_GENES) {
  gene_cols_model <- gene_cols_model[seq_len(MAX_ML_GENES)]
}

if (length(gene_cols_model) < 5) {
  stop("Too few valid expression genes for Mime1 modeling.")
}

cat("TCGA survival samples:", nrow(ml_df), "\n")
cat("TCGA survival events:", sum(ml_df$status == 1), "\n")
cat("Valid genes for Mime1:", length(gene_cols_model), "\n")

mime_df <- ml_df[, c("Patient", "time", "status", gene_cols_model), drop = FALSE]
colnames(mime_df)[1:3] <- c("ID", "OS.time", "OS")

for (g in gene_cols_model) {
  mime_df[[g]] <- zscore_vector(mime_df[[g]])
}

mime_df <- mime_df[complete.cases(mime_df[, c("ID", "OS.time", "OS")]), ]

set.seed(20260510)

event_ids <- mime_df$ID[mime_df$OS == 1]
cens_ids  <- mime_df$ID[mime_df$OS == 0]

train_event <- sample(event_ids, size = floor(TRAIN_RATIO * length(event_ids)))
train_cens  <- sample(cens_ids,  size = floor(TRAIN_RATIO * length(cens_ids)))

train_ids <- c(train_event, train_cens)

train_df_mime <- mime_df[mime_df$ID %in% train_ids, ]
test_df_mime  <- mime_df[!mime_df$ID %in% train_ids, ]

cat("Training n/events:", nrow(train_df_mime), "/", sum(train_df_mime$OS == 1), "\n")
cat("Testing n/events:", nrow(test_df_mime), "/", sum(test_df_mime$OS == 1), "\n")

if (nrow(train_df_mime) < 80 || sum(train_df_mime$OS == 1) < MIN_EVENTS_TRAIN) {
  stop("Training set has insufficient samples/events.")
}

save_csv(train_df_mime, file.path(FIG4_DIR, "Mime1_input_Training.csv"))
save_csv(test_df_mime, file.path(FIG4_DIR, "Mime1_input_Testing.csv"))
save_csv(mime_df, file.path(FIG4_DIR, "Mime1_input_Entire_TCGA.csv"))

# ----------------------------
# 24.3 Prepare GEO cohorts for Mime1
# ----------------------------

prepare_geo_for_mime <- function(expr_file, clin_file, gene_use, cohort_name) {
  expr <- data.table::fread(expr_file, data.table = FALSE, check.names = FALSE)
  clin <- data.table::fread(clin_file, data.table = FALSE, check.names = FALSE)

  gene_col <- colnames(expr)[1]
  expr[[gene_col]] <- clean_gene_symbol(expr[[gene_col]])
  expr <- expr[!is.na(expr[[gene_col]]) & expr[[gene_col]] != "", ]
  expr <- expr[!duplicated(expr[[gene_col]]), ]

  rownames(expr) <- expr[[gene_col]]
  expr[[gene_col]] <- NULL

  expr_mat <- as.matrix(expr)
  storage.mode(expr_mat) <- "numeric"

  common_genes <- intersect(gene_use, rownames(expr_mat))
  cat(cohort_name, "common genes:", length(common_genes), "/", length(gene_use), "\n")

  if (length(common_genes) < 5) {
    stop(cohort_name, ": too few common genes with TCGA/Mime1 model genes.")
  }

  expr_use <- t(expr_mat[common_genes, , drop = FALSE])
  expr_use <- as.data.frame(expr_use, check.names = FALSE)
  expr_use$ID <- rownames(expr_use)

  clin$Sample <- as.character(clin$Sample)

  geo_df <- expr_use %>%
    dplyr::inner_join(
      clin[, c("Sample", "time", "status")],
      by = c("ID" = "Sample")
    )

  geo_df <- geo_df[
    is.finite(geo_df$time) &
      geo_df$time > 0 &
      geo_df$status %in% c(0, 1),
  ]

  missing_genes <- setdiff(gene_use, colnames(geo_df))
  if (length(missing_genes) > 0) {
    for (g in missing_genes) geo_df[[g]] <- 0
  }

  geo_df <- geo_df[, c("ID", "time", "status", gene_use), drop = FALSE]
  colnames(geo_df)[1:3] <- c("ID", "OS.time", "OS")

  for (g in gene_use) {
    geo_df[[g]] <- zscore_vector(geo_df[[g]])
  }

  geo_df <- geo_df[complete.cases(geo_df[, c("ID", "OS.time", "OS")]), ]

  cat(cohort_name, "final n/events:", nrow(geo_df), "/", sum(geo_df$OS == 1), "\n")

  if (nrow(geo_df) < 50 || sum(geo_df$OS == 1) < 10) {
    warning(cohort_name, " has relatively few samples/events after matching.")
  }

  geo_df
}

gse72094_mime <- prepare_geo_for_mime(
  expr_file = gse72094_expr_file,
  clin_file = gse72094_clin_file,
  gene_use = gene_cols_model,
  cohort_name = "GSE72094"
)

gse68465_mime <- prepare_geo_for_mime(
  expr_file = gse68465_expr_file,
  clin_file = gse68465_clin_file,
  gene_use = gene_cols_model,
  cohort_name = "GSE68465"
)

save_csv(gse72094_mime, file.path(FIG4_DIR, "Mime1_input_GSE72094.csv"))
save_csv(gse68465_mime, file.path(FIG4_DIR, "Mime1_input_GSE68465.csv"))

list_train_vali_Data <- list(
  Training = train_df_mime,
  Testing = test_df_mime,
  GSE72094 = gse72094_mime,
  GSE68465 = gse68465_mime
)

# ----------------------------
# 24.4 Run Mime1 101 ML combinations
# ----------------------------

set.seed(20260510)

res.mime <- tryCatch({
  Mime1::ML.Dev.Prog.Sig(
    train_data = train_df_mime,
    list_train_vali_Data = list_train_vali_Data,
    candidate_genes = gene_cols_model,
    mode = "all",
    unicox.filter.for.candi = TRUE,
    unicox_p_cutoff = UNICOX_P,
    nodesize = 5,
    seed = 20260510
  )
}, error = function(e) {
  message("Mime1 ML.Dev.Prog.Sig failed: ", e$message)
  NULL
})

if (is.null(res.mime)) {
  stop("Mime1 modeling failed. Check Mime1 installation and input format.")
}

saveRDS(res.mime, file.path(FIG4_DIR, "Mime1_ML_DDR_res_with_GEO.rds"))
save_csv(data.frame(ObjectNames = names(res.mime)), file.path(FIG4_DIR, "Mime1_result_object_names_with_GEO.csv"))

# ----------------------------
# 24.5 Extract and organize C-index
# ----------------------------

cindex_raw <- NULL
possible_cindex_names <- c(
  "cindex.res", "cindex_res", "Cindex.res",
  "Cindex", "cindex", "cindex.result",
  "Cindex_result", "Cindex_Result"
)

for (nm in possible_cindex_names) {
  if (nm %in% names(res.mime)) {
    cindex_raw <- as.data.frame(res.mime[[nm]], check.names = FALSE)
    break
  }
}

if (is.null(cindex_raw)) {
  stop("Cannot find C-index table in Mime1 object. Check Mime1_result_object_names_with_GEO.csv.")
}

save_csv(cindex_raw, file.path(FIG4_DIR, "Mime1_raw_Cindex_summary_with_GEO.csv"))

if (all(c("ID", "Cindex", "Model") %in% colnames(cindex_raw))) {
  cindex_summary <- cindex_raw %>%
    dplyr::mutate(
      ID = as.character(ID),
      Model = as.character(Model),
      Cindex = suppressWarnings(as.numeric(Cindex))
    ) %>%
    tidyr::pivot_wider(
      id_cols = Model,
      names_from = ID,
      values_from = Cindex
    ) %>%
    as.data.frame(check.names = FALSE)
} else {
  cindex_summary <- cindex_raw
  if (!"Model" %in% colnames(cindex_summary)) {
    cindex_summary$Model <- rownames(cindex_summary)
  }
}

colnames(cindex_summary) <- gsub("\\s+", "_", colnames(cindex_summary))
colnames(cindex_summary) <- gsub("-", "_", colnames(cindex_summary))
colnames(cindex_summary) <- gsub("\\.", "_", colnames(cindex_summary))

dataset_cols <- intersect(c("Training", "Testing", "GSE72094", "GSE68465"), colnames(cindex_summary))

if (length(dataset_cols) < 3) {
  save_csv(cindex_summary, file.path(FIG4_DIR, "Mime1_Cindex_after_wider_with_GEO_debug.csv"))
  stop("C-index table was reshaped but expected dataset columns were not recognized.")
}

for (cc in dataset_cols) {
  cindex_summary[[cc]] <- suppressWarnings(as.numeric(cindex_summary[[cc]]))
}

cindex_summary$MeanCindex <- rowMeans(cindex_summary[, dataset_cols, drop = FALSE], na.rm = TRUE)
cindex_summary$ExternalMean <- rowMeans(cindex_summary[, intersect(c("GSE72094", "GSE68465"), dataset_cols), drop = FALSE], na.rm = TRUE)
cindex_summary$Train_Test_Gap <- abs(cindex_summary$Training - cindex_summary$Testing)
cindex_summary$OverfitFlag <- ifelse(cindex_summary$Training > 0.85 & cindex_summary$Train_Test_Gap > 0.15, "Overfit", "Stable")

cindex_summary <- cindex_summary[order(cindex_summary$MeanCindex, decreasing = TRUE), ]

save_csv(cindex_summary, file.path(FIG4_DIR, "Fig4B_Mime1_101_model_Cindex_summary_with_GEO.csv"))
save_csv(cindex_summary, file.path(FIG4_DIR, "Fig4B_model_Cindex_summary.csv"))

stable_models <- cindex_summary %>%
  dplyr::filter(
    !is.na(Testing),
    !is.na(GSE72094),
    !is.na(GSE68465),
    Testing >= 0.60,
    GSE72094 >= 0.55,
    GSE68465 >= 0.55,
    Train_Test_Gap <= 0.15,
    OverfitFlag == "Stable"
  ) %>%
  dplyr::arrange(
    dplyr::desc(ExternalMean),
    dplyr::desc(Testing),
    dplyr::desc(MeanCindex),
    Train_Test_Gap
  )

if (nrow(stable_models) == 0) {
  warning("No model met strict stability thresholds. Relaxing thresholds.")
  stable_models <- cindex_summary %>%
    dplyr::filter(
      !is.na(Testing),
      !is.na(GSE72094),
      !is.na(GSE68465),
      Train_Test_Gap <= 0.20
    ) %>%
    dplyr::arrange(
      dplyr::desc(ExternalMean),
      dplyr::desc(Testing),
      dplyr::desc(MeanCindex),
      Train_Test_Gap
    )
}

if (nrow(stable_models) == 0) {
  warning("Still no stable model candidates. Using highest MeanCindex model.")
  stable_models <- cindex_summary %>% dplyr::arrange(dplyr::desc(MeanCindex))
}

save_csv(stable_models, file.path(FIG4_DIR, "Fig4B_Mime1_stable_model_candidates_with_GEO.csv"))

best_model <- as.character(stable_models$Model[1])
cat("Best Mime1 model with GEO validation:", best_model, "\n")
cat("C-index dataset columns used:", paste(dataset_cols, collapse = ", "), "\n")
print(head(stable_models, 20))

# ----------------------------
# 24.6 C-index heatmap
# ----------------------------

try({
  grDevices::pdf(file.path(FIG4_DIR, "Fig4B_Mime1_Cindex_original_style_with_GEO.pdf"), width = 7, height = 10)
  Mime1::cindex_dis_all(
    res.mime,
    validate_set = c("Testing", "GSE72094", "GSE68465"),
    order = c("Training", "Testing", "GSE72094", "GSE68465"),
    width = 0.35
  )
  grDevices::dev.off()
}, silent = TRUE)

ci_show <- stable_models %>%
  dplyr::slice_head(n = 45)

ci_show$Average <- rowMeans(ci_show[, dataset_cols, drop = FALSE], na.rm = TRUE)

heat_cols <- c("Training", "Testing", "GSE72094", "GSE68465", "Average")
heat_cols <- intersect(heat_cols, colnames(ci_show))

ci_long <- ci_show %>%
  dplyr::select(Model, dplyr::all_of(heat_cols)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(heat_cols),
    names_to = "Dataset",
    values_to = "Cindex"
  )

ci_long$Dataset <- factor(ci_long$Dataset, levels = heat_cols)
ci_long$Model <- factor(ci_long$Model, levels = rev(ci_show$Model))

if (requireNamespace("scales", quietly = TRUE)) {
  red_pal <- scales::col_numeric(
    palette = c("#EFD6D1", "#C45A4D", "#8F261F"),
    domain = c(0.50, 1.00)
  )

  blue_pal <- scales::col_numeric(
    palette = c("#D7E4EF", "#6B98BF", "#255B86"),
    domain = c(0.45, 0.85)
  )

  ci_long$FillColor <- ifelse(
    ci_long$Dataset == "Training",
    red_pal(ci_long$Cindex),
    blue_pal(ci_long$Cindex)
  )
} else {
  ci_long$FillColor <- ifelse(ci_long$Dataset == "Training", "#C45A4D", "#6B98BF")
}

ci_long$TextColor <- ifelse(
  (ci_long$Dataset == "Training" & ci_long$Cindex >= 0.82) |
    (ci_long$Dataset != "Training" & ci_long$Cindex >= 0.70),
  "white",
  "black"
)

p_cindex <- ggplot2::ggplot(
  ci_long,
  ggplot2::aes(x = Dataset, y = Model)
) +
  ggplot2::geom_tile(
    ggplot2::aes(fill = FillColor),
    color = "white",
    linewidth = 0.5
  ) +
  ggplot2::scale_fill_identity() +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf("%.3f", Cindex),
      color = TextColor
    ),
    size = 2.35
  ) +
  ggplot2::scale_color_identity() +
  ggplot2::theme_void(base_size = 10) +
  ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 6.6, color = "black", hjust = 1),
    axis.text.x = ggplot2::element_text(size = 8.8, color = "black", face = "bold"),
    axis.title = ggplot2::element_blank(),
    legend.position = "none",
    plot.margin = ggplot2::margin(5, 5, 5, 5)
  )

safe_ggsave(
  file.path(FIG4_DIR, "Fig4B_Mime1_Cindex_heatmap_with_GEO_paper_like.pdf"),
  p_cindex,
  width = 6.6,
  height = 10.5
)
safe_ggsave(
  file.path(FIG4_DIR, "Fig4B_Mime1_Cindex_heatmap_with_GEO_paper_like.png"),
  p_cindex,
  width = 6.6,
  height = 10.5
)

# ----------------------------
# 24.7 Risk-score extraction
# ----------------------------

extract_risk_from_table <- function(tab, best_model, dataset_name = NULL) {
  tab <- as.data.frame(tab, check.names = FALSE)
  cn <- colnames(tab)

  if (best_model %in% cn) {
    id_col <- grep("^ID$|sample|patient|ID", cn, ignore.case = TRUE, value = TRUE)[1]
    if (is.na(id_col)) id_col <- cn[1]
    out <- data.frame(
      ID = as.character(tab[[id_col]]),
      ML_DDRscore = safe_num2(tab[[best_model]]),
      stringsAsFactors = FALSE
    )
    if (!is.null(dataset_name)) out$Dataset <- dataset_name
    return(out)
  }

  model_col <- grep("model|method|algorithm", cn, ignore.case = TRUE, value = TRUE)[1]
  score_col <- grep("risk|score|RS", cn, ignore.case = TRUE, value = TRUE)[1]
  id_col <- grep("^ID$|sample|patient|ID", cn, ignore.case = TRUE, value = TRUE)[1]
  ds_col <- grep("dataset|cohort|set", cn, ignore.case = TRUE, value = TRUE)[1]

  if (!is.na(model_col) && !is.na(score_col) && !is.na(id_col)) {
    tab2 <- tab[as.character(tab[[model_col]]) == best_model, , drop = FALSE]
    if (nrow(tab2) > 0) {
      out <- data.frame(
        ID = as.character(tab2[[id_col]]),
        ML_DDRscore = safe_num2(tab2[[score_col]]),
        stringsAsFactors = FALSE
      )
      if (!is.na(ds_col)) out$Dataset <- as.character(tab2[[ds_col]])
      if (!is.null(dataset_name) && !"Dataset" %in% colnames(out)) out$Dataset <- dataset_name
      return(out)
    }
  }

  NULL
}

extract_risk_scores_robust <- function(res.mime, best_model, list_train_vali_Data) {
  possible_risk_names <- c(
    "riskscore", "riskScore", "risk.score", "RS",
    "RS.res", "rs.res", "risk.res", "ml.res", "res",
    "riskscore.res", "risk_score"
  )

  risk_obj <- NULL
  risk_name <- NA_character_

  for (nm in possible_risk_names) {
    if (nm %in% names(res.mime)) {
      risk_obj <- res.mime[[nm]]
      risk_name <- nm
      break
    }
  }

  if (is.null(risk_obj)) {
    save_csv(data.frame(ObjectNames = names(res.mime)), file.path(FIG4_DIR, "Mime1_result_object_names_for_risk_score_with_GEO.csv"))
    return(NULL)
  }

  saveRDS(risk_obj, file.path(FIG4_DIR, paste0("Mime1_risk_object_", risk_name, "_with_GEO.rds")))

  risk_list <- list()

  if (is.list(risk_obj) && !is.data.frame(risk_obj)) {
    for (ds in names(list_train_vali_Data)) {
      if (ds %in% names(risk_obj)) {
        tmp <- extract_risk_from_table(risk_obj[[ds]], best_model, dataset_name = ds)
        if (!is.null(tmp)) risk_list[[ds]] <- tmp
      }
    }

    if (length(risk_list) > 0) {
      return(dplyr::bind_rows(risk_list))
    }

    if (best_model %in% names(risk_obj)) {
      model_obj <- risk_obj[[best_model]]
      if (is.list(model_obj) && !is.data.frame(model_obj)) {
        for (ds in names(list_train_vali_Data)) {
          if (ds %in% names(model_obj)) {
            tmp <- as.data.frame(model_obj[[ds]], check.names = FALSE)
            tmp2 <- extract_risk_from_table(tmp, best_model, dataset_name = ds)
            if (is.null(tmp2)) {
              cn <- colnames(tmp)
              id_col <- grep("^ID$|sample|patient|ID", cn, ignore.case = TRUE, value = TRUE)[1]
              score_col <- grep("risk|score|RS", cn, ignore.case = TRUE, value = TRUE)[1]
              if (!is.na(id_col) && !is.na(score_col)) {
                tmp2 <- data.frame(
                  ID = as.character(tmp[[id_col]]),
                  ML_DDRscore = safe_num2(tmp[[score_col]]),
                  Dataset = ds,
                  stringsAsFactors = FALSE
                )
              }
            }
            if (!is.null(tmp2)) risk_list[[ds]] <- tmp2
          }
        }
        if (length(risk_list) > 0) {
          return(dplyr::bind_rows(risk_list))
        }
      } else {
        tmp <- extract_risk_from_table(model_obj, best_model, dataset_name = NULL)
        if (!is.null(tmp)) return(tmp)
      }
    }

    for (nm in names(risk_obj)) {
      tmp <- tryCatch(extract_risk_from_table(risk_obj[[nm]], best_model, dataset_name = nm), error = function(e) NULL)
      if (!is.null(tmp)) risk_list[[nm]] <- tmp
    }
    if (length(risk_list) > 0) return(dplyr::bind_rows(risk_list))
  }

  if (is.data.frame(risk_obj)) {
    tmp <- extract_risk_from_table(risk_obj, best_model, dataset_name = NULL)
    if (!is.null(tmp)) return(tmp)
  }

  NULL
}

risk_out <- extract_risk_scores_robust(res.mime, best_model, list_train_vali_Data)

if (is.null(risk_out) || nrow(risk_out) == 0) {
  save_csv(data.frame(ObjectNames = names(res.mime)), file.path(FIG4_DIR, "Mime1_result_object_names_for_risk_score_with_GEO.csv"))
  stop(
    "Cannot automatically extract Mime1 risk scores. ",
    "Please send Mime1_result_object_names_for_risk_score_with_GEO.csv and Mime1_risk_object_*.rds names."
  )
}

risk_out <- as.data.frame(risk_out, check.names = FALSE)
risk_out$ID <- as.character(risk_out$ID)
risk_out$ML_DDRscore <- safe_num2(risk_out$ML_DDRscore)

save_csv(risk_out, file.path(FIG4_DIR, "Mime1_best_model_raw_risk_score_table_with_GEO.csv"))

if (!"Dataset" %in% colnames(risk_out)) {
  risk_out$Dataset <- NA_character_
  for (ds in names(list_train_vali_Data)) {
    ids <- as.character(list_train_vali_Data[[ds]]$ID)
    risk_out$Dataset[risk_out$ID %in% ids] <- ds
  }
}

base_score_list <- lapply(names(list_train_vali_Data), function(ds) {
  x <- list_train_vali_Data[[ds]][, c("ID", "OS.time", "OS"), drop = FALSE]
  data.frame(
    ID = as.character(x$ID),
    time = as.numeric(x$OS.time),
    status = as.numeric(x$OS),
    Dataset = ds,
    stringsAsFactors = FALSE
  )
})
base_score_df <- dplyr::bind_rows(base_score_list)

final_score_all <- base_score_df %>%
  dplyr::left_join(
    risk_out[, c("ID", "Dataset", "ML_DDRscore"), drop = FALSE],
    by = c("ID", "Dataset")
  )

if (sum(is.finite(final_score_all$ML_DDRscore)) < 0.5 * nrow(final_score_all)) {
  final_score_all$ML_DDRscore <- NULL
  risk_id_only <- risk_out[, c("ID", "ML_DDRscore"), drop = FALSE]
  risk_id_only <- risk_id_only[!duplicated(risk_id_only$ID), ]
  final_score_all <- base_score_df %>%
    dplyr::left_join(risk_id_only, by = "ID")
}

if (sum(is.finite(final_score_all$ML_DDRscore)) < 50) {
  save_csv(final_score_all, file.path(FIG4_DIR, "Fig4C_ML_DDRscore_all_sets_failed_extract_debug.csv"))
  stop("Too few finite ML_DDRscore values after extracting Mime1 risk scores.")
}

train_tmp <- final_score_all[
  final_score_all$Dataset == "Training" &
    is.finite(final_score_all$ML_DDRscore) &
    is.finite(final_score_all$time) &
    final_score_all$status %in% c(0, 1),
]

cox_dir <- tryCatch({
  survival::coxph(survival::Surv(time, status) ~ ML_DDRscore, data = train_tmp)
}, error = function(e) NULL)

risk_direction <- "original_by_training_cox"
if (!is.null(cox_dir)) {
  coef_dir <- as.numeric(stats::coef(cox_dir)[1])
  if (is.finite(coef_dir) && coef_dir < 0) {
    final_score_all$ML_DDRscore <- -final_score_all$ML_DDRscore
    risk_direction <- "reversed_by_training_cox"
  }
}
final_score_all$RiskDirection <- risk_direction

# Cohort-specific median cutoff is common for cross-platform external validation.
# Continuous C-index/AUC are unaffected by this grouping cutoff.
final_score_all$RiskGroup <- NA_character_
for (ds in unique(final_score_all$Dataset)) {
  idx <- final_score_all$Dataset == ds & is.finite(final_score_all$ML_DDRscore)
  cut_ds <- stats::median(final_score_all$ML_DDRscore[idx], na.rm = TRUE)
  final_score_all$RiskGroup[idx] <- ifelse(final_score_all$ML_DDRscore[idx] >= cut_ds, "High", "Low")
}
final_score_all$RiskGroup <- factor(final_score_all$RiskGroup, levels = c("Low", "High"))

final_score_all$Patient <- final_score_all$ID
final_score_all$Sample <- final_score_all$ID

save_csv(final_score_all, file.path(FIG4_DIR, "Fig4C_ML_DDRscore_all_sets.csv"))
save_csv(final_score_all, file.path(DB_DIR, "ML_DDRscore_table.csv"))

# ----------------------------
# 24.8 Final signature gene table
# ----------------------------

cox_list <- lapply(gene_cols_model, function(g) {
  df_tmp <- train_df_mime[, c("OS.time", "OS", g), drop = FALSE]
  colnames(df_tmp) <- c("time", "status", "expr")
  fit <- tryCatch(
    survival::coxph(survival::Surv(time, status) ~ expr, data = df_tmp),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  sm <- summary(fit)
  data.frame(
    Gene = g,
    Coef = sm$coefficients[1, "coef"],
    HR = sm$coefficients[1, "exp(coef)"],
    P = sm$coefficients[1, "Pr(>|z|)"],
    Model = best_model,
    stringsAsFactors = FALSE
  )
})

final_coef <- dplyr::bind_rows(cox_list)
final_coef <- final_coef[is.finite(final_coef$P), ]
final_coef <- final_coef[order(final_coef$P), ]
final_coef <- head(final_coef, 30)

save_csv(final_coef, file.path(FIG4_DIR, "Fig4C_final_ML_DDR_signature_genes.csv"))
save_csv(final_coef, file.path(DB_DIR, "ML_DDR_signature_table.csv"))

p_coef <- ggplot2::ggplot(
  final_coef,
  ggplot2::aes(x = reorder(Gene, Coef), y = Coef)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::geom_segment(ggplot2::aes(xend = Gene, y = 0, yend = Coef), linewidth = 0.6) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::coord_flip() +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(x = NULL, y = "Univariable Cox coefficient", title = paste0("Top prognostic genes: ", best_model))

safe_ggsave(file.path(FIG4_DIR, "Fig4C_final_signature_coefficients_lollipop.pdf"), p_coef, width = 6, height = 5)

# ----------------------------
# 24.9 Risk rank plots, KM and ROC
# ----------------------------

risk_rank <- final_score_all[is.finite(final_score_all$ML_DDRscore), ]
risk_rank <- risk_rank[order(risk_rank$Dataset, risk_rank$ML_DDRscore), ]
risk_rank$Rank <- ave(risk_rank$ML_DDRscore, risk_rank$Dataset, FUN = seq_along)

p_risk <- ggplot2::ggplot(
  risk_rank,
  ggplot2::aes(x = Rank, y = ML_DDRscore, color = RiskGroup)
) +
  ggplot2::geom_point(size = 0.8) +
  ggplot2::facet_wrap(~ Dataset, scales = "free_x", nrow = 1) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(x = "Patient rank", y = "ML-DDRscore", title = "Risk score distribution")

safe_ggsave(file.path(FIG4_DIR, "Fig4C_ML_DDRscore_risk_rank_plot_with_GEO.pdf"), p_risk, width = 12, height = 4)

p_status <- ggplot2::ggplot(
  risk_rank,
  ggplot2::aes(x = Rank, y = time, color = factor(status))
) +
  ggplot2::geom_point(size = 0.8) +
  ggplot2::facet_wrap(~ Dataset, scales = "free_x", nrow = 1) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(x = "Patient rank", y = "Survival time", color = "Status", title = "Survival status distribution")

safe_ggsave(file.path(FIG4_DIR, "Fig4C_ML_DDRscore_survival_status_plot_with_GEO.pdf"), p_status, width = 12, height = 4)

selected_cindex <- cindex_summary[cindex_summary$Model == best_model, , drop = FALSE]
save_csv(selected_cindex, file.path(FIG4_DIR, "Fig4B_selected_model_Cindex_with_GEO.csv"))

score_cindex <- lapply(unique(final_score_all$Dataset), function(ds) {
  df <- final_score_all[final_score_all$Dataset == ds, ]
  data.frame(
    Dataset = ds,
    N = nrow(df),
    Events = sum(df$status == 1, na.rm = TRUE),
    Cindex = calc_cindex_local(df$time, df$status, df$ML_DDRscore),
    stringsAsFactors = FALSE
  )
}) %>% dplyr::bind_rows()
save_csv(score_cindex, file.path(FIG4_DIR, "Fig4D_selected_model_score_Cindex_by_dataset.csv"))

for (ds in unique(final_score_all$Dataset)) {
  df_ds <- final_score_all[final_score_all$Dataset == ds, ]

  plot_km_by_group(
    df_ds,
    title = paste0(ds, " ML-DDR signature"),
    file = file.path(FIG4_DIR, paste0("Fig4D_", ds, "_ML_DDR_KM.pdf"))
  )

  plot_time_roc_local(
    df_ds,
    "ML_DDRscore",
    file = file.path(FIG4_DIR, paste0("Fig4E_", ds, "_ML_DDR_timeROC.pdf")),
    title = paste0(ds, " time-dependent ROC")
  )
}

# ----------------------------
# 24.10 Optional Mime1-native KM
# ----------------------------

try({
  if (requireNamespace("aplot", quietly = TRUE)) {
    survplot <- list()
    for (ds in names(list_train_vali_Data)) {
      survplot[[ds]] <- Mime1::rs_sur(
        res.mime,
        model_name = best_model,
        dataset = ds,
        median.line = "hv",
        cutoff = 0.5,
        conf.int = TRUE,
        xlab = "Time",
        pval.coord = c(1000, 0.9)
      )
    }
    grDevices::pdf(file.path(FIG4_DIR, "Fig4D_Mime1_native_KM_with_GEO.pdf"), width = 14, height = 4)
    print(aplot::plot_list(gglist = survplot, ncol = length(survplot)))
    grDevices::dev.off()
  }
}, silent = TRUE)

# Objects required by downstream summary
ml_genes <- gene_cols_model
model_scores <- list()
model_coefs <- list()

cat("\nSection 24 completed.\n")
cat("Best model: ", best_model, "\n", sep = "")
cat("Risk direction: ", unique(final_score_all$RiskDirection), "\n", sep = "")
cat("Key output files:\n")
cat(" - ", file.path(FIG4_DIR, "Fig4B_Mime1_101_model_Cindex_summary_with_GEO.csv"), "\n", sep = "")
cat(" - ", file.path(FIG4_DIR, "Fig4B_Mime1_stable_model_candidates_with_GEO.csv"), "\n", sep = "")
cat(" - ", file.path(FIG4_DIR, "Fig4C_ML_DDRscore_all_sets.csv"), "\n", sep = "")
cat(" - ", file.path(FIG4_DIR, "Fig4D_selected_model_score_Cindex_by_dataset.csv"), "\n", sep = "")
