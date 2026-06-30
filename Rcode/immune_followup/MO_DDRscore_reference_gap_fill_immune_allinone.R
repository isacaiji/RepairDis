############################################################
# Reference-gap immune analyses for MO-DDRscore / NMF subtype
# Purpose:
#   Fill the immune-analysis modules that were present in the
#   reference paper but not yet fully covered in the current
#   MO-DDRscore immune workflow.
#
# This script does NOT fabricate unavailable external resources.
# Analyses that require official external files are run only when
# those files are present in 00_data/immune_external.
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260513)

############################
# 0. Packages
############################

pkgs <- c("data.table", "dplyr", "tidyr", "stringr", "survival")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(survival)
})

############################
# 1. Paths
############################

PROJECT_DIR <- "D:/R_workspace/评分/AD_DDR_project"
PROC_DIR <- file.path(PROJECT_DIR, "01_processed")
FIG2_DIR <- file.path(PROJECT_DIR, "02_Figure2_MultiOmics_Immune")

BASIC_DIR <- file.path(FIG2_DIR, "Immune_basic_MO_DDRscore")
FOLLOW_DIR <- file.path(FIG2_DIR, "Immune_official_followup_MO_DDRscore")
IOBR_DIR <- file.path(FIG2_DIR, "Immune_IOBR_signature_MO_DDRscore")
NMF_DIR <- file.path(FIG2_DIR, "NMF_DDR_subtype_quick")

EXT_DIR <- Sys.getenv(
  "IMMUNE_EXTERNAL_DIR",
  unset = file.path(PROJECT_DIR, "00_data", "immune_external")
)
OUT_DIR <- Sys.getenv(
  "IMMUNE_GAP_OUT_DIR",
  unset = file.path(FIG2_DIR, "Immune_reference_gap_fill_allinone")
)

dir.create(EXT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SCORE_FILE <- file.path(PROC_DIR, "LUAD_MO_DDRscore.csv")
EXPR_FILE <- file.path(PROC_DIR, "LUAD_tcga_expr_tpm_matrix.rds")
CLIN_FILE <- file.path(PROC_DIR, "LUAD_clinical_processed.csv")
NMF_FILE <- file.path(NMF_DIR, "DDR_NMF_subtype_with_MO_DDRscore.csv")

DECONV_LONG_FILE <- file.path(BASIC_DIR, "Immune_deconvolution_long.csv")
MARKER_LONG_FILE <- file.path(BASIC_DIR, "Immune_marker_expression_long.csv")
IOBR_LONG_FILE <- file.path(IOBR_DIR, "IOBR_official_signature_scores_long.csv")
TIDE_FILE <- file.path(FOLLOW_DIR, "Official_TIDE_merged.csv")
TMB_FILE <- file.path(FOLLOW_DIR, "Official_TMB_merged.csv")
MATH_FILE <- file.path(FOLLOW_DIR, "Official_MATH_score.csv")

EXT_IMMUNE_SUBTYPE_FILE <- file.path(EXT_DIR, "TCGA_pan_cancer_immune_subtype_annotation.csv")
EXT_TCIA_IPS_FILE <- file.path(EXT_DIR, "TCIA_IPS_LUAD.csv")
EXT_TNB_FILE <- file.path(EXT_DIR, "TNB_neoantigen_burden.csv")
EXT_MSI_STEMNESS_FILE <- file.path(EXT_DIR, "MSI_stemness_annotation.csv")
EXT_CHECKPOINT_FILE <- file.path(EXT_DIR, "official_checkpoint_gene_panel.csv")
EXT_IMMUNOTHERAPY_FILE <- file.path(EXT_DIR, "immunotherapy_response_cohort.csv")

############################
# 2. Helpers
############################

read_csv <- function(file) {
  data.table::fread(file, data.table = FALSE, check.names = FALSE)
}

save_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file)
}

num <- function(x) suppressWarnings(as.numeric(x))

norm_sample <- function(x) gsub("\\.", "-", as.character(x))

patient_id <- function(x) substr(norm_sample(x), 1, 12)

clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub(" ///.*$| //.*$|;.*$|,.*$|\\s*\\([^\\)]*\\)$", "", x)
  x[x %in% c("", "NA", "---", "NULL", "N/A")] <- NA
  x
}

sig_star <- function(p) {
  p <- num(p)
  ifelse(
    is.na(p), "",
    ifelse(p < 0.001, "***",
           ifelse(p < 0.01, "**",
                  ifelse(p < 0.05, "*", "ns")))
  )
}

calc_cramers_v <- function(tab) {
  if (min(dim(tab)) < 2 || sum(tab) == 0) return(NA_real_)
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
  n <- sum(tab)
  denom <- n * min(nrow(tab) - 1, ncol(tab) - 1)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  sqrt(as.numeric(chi$statistic) / denom)
}

standardize_patient_col <- function(df) {
  cn <- colnames(df)
  hit <- intersect(c("Patient", "patient", "PATIENT", "bcr_patient_barcode", "case_submitter_id"), cn)
  if (length(hit) == 0) {
    hit <- cn[grepl("patient|barcode|case", cn, ignore.case = TRUE)][1]
  }
  if (is.na(hit) || length(hit) == 0) stop("Cannot identify patient column.")
  df$Patient <- patient_id(df[[hit]])
  df
}

compare_long <- function(df, group_col, value_col, feature_cols, level_order = NULL) {
  stopifnot(group_col %in% colnames(df), value_col %in% colnames(df))
  d <- df %>%
    filter(!is.na(.data[[group_col]]), is.finite(num(.data[[value_col]]))) %>%
    mutate(
      .group = as.character(.data[[group_col]]),
      .value = num(.data[[value_col]])
    )
  if (!is.null(level_order)) {
    d$.group <- factor(d$.group, levels = level_order)
  }
  d <- d %>% filter(!is.na(.group))
  if (nrow(d) == 0) return(data.frame())

  missing_cols <- setdiff(feature_cols, colnames(d))
  if (length(missing_cols) > 0) stop("Missing feature columns: ", paste(missing_cols, collapse = ", "))

  d$.key <- do.call(paste, c(d[feature_cols], sep = "\r"))
  pieces <- split(d, d$.key)

  out <- lapply(pieces, function(x) {
    g <- unique(as.character(x$.group))
    g <- g[!is.na(g)]
    if (!is.null(level_order)) g <- intersect(level_order, g) else g <- sort(g)
    if (length(g) != 2) return(NULL)

    x1 <- x$.value[as.character(x$.group) == g[1]]
    x2 <- x$.value[as.character(x$.group) == g[2]]
    if (sum(is.finite(x1)) < 3 || sum(is.finite(x2)) < 3) return(NULL)

    p <- tryCatch(stats::wilcox.test(x2, x1)$p.value, error = function(e) NA_real_)
    feature_values <- x[1, feature_cols, drop = FALSE]
    data.frame(
      feature_values,
      Group1 = g[1],
      Group2 = g[2],
      N_Group1 = sum(is.finite(x1)),
      N_Group2 = sum(is.finite(x2)),
      Median_Group1 = median(x1, na.rm = TRUE),
      Median_Group2 = median(x2, na.rm = TRUE),
      Delta_Group2_minus_Group1 = median(x2, na.rm = TRUE) - median(x1, na.rm = TRUE),
      P = p,
      stringsAsFactors = FALSE
    )
  })

  out <- bind_rows(out)
  if (nrow(out) == 0) return(out)
  out$FDR <- p.adjust(out$P, method = "BH")
  out$Significance <- sig_star(out$FDR)
  out
}

categorical_test <- function(df, group_col, cat_col, level_order = NULL) {
  stopifnot(group_col %in% colnames(df), cat_col %in% colnames(df))
  d <- df %>%
    filter(!is.na(.data[[group_col]]), !is.na(.data[[cat_col]])) %>%
    mutate(.group = as.character(.data[[group_col]]), .cat = as.character(.data[[cat_col]]))
  if (!is.null(level_order)) d$.group <- factor(d$.group, levels = level_order)
  d <- d %>% filter(!is.na(.group))
  if (length(unique(d$.group)) < 2 || length(unique(d$.cat)) < 2) return(data.frame())

  tab <- table(d$.group, d$.cat)
  fisher_p <- tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)
  chisq_p <- tryCatch(suppressWarnings(chisq.test(tab, correct = FALSE)$p.value), error = function(e) NA_real_)

  counts <- as.data.frame(tab)
  colnames(counts) <- c("Group", "Category", "N")
  counts <- counts %>%
    group_by(Group) %>%
    mutate(Proportion = N / sum(N)) %>%
    ungroup()

  list(
    counts = counts,
    test = data.frame(
      GroupVariable = group_col,
      CategoryVariable = cat_col,
      Fisher_P = fisher_p,
      ChiSquare_P = chisq_p,
      Cramers_V = calc_cramers_v(tab),
      stringsAsFactors = FALSE
    )
  )
}

survival_test <- function(df, group_col, time_col = "time", status_col = "status") {
  stopifnot(group_col %in% colnames(df), time_col %in% colnames(df), status_col %in% colnames(df))
  d <- df %>%
    mutate(
      .time = num(.data[[time_col]]),
      .status = num(.data[[status_col]]),
      .group = as.factor(.data[[group_col]])
    ) %>%
    filter(is.finite(.time), .time > 0, .status %in% c(0, 1), !is.na(.group))
  if (nrow(d) < 20 || length(unique(d$.group)) < 2) return(data.frame())

  sd <- tryCatch(survdiff(Surv(.time, .status) ~ .group, data = d), error = function(e) NULL)
  logrank_p <- if (is.null(sd)) NA_real_ else pchisq(sd$chisq, df = length(sd$n) - 1, lower.tail = FALSE)

  cox <- tryCatch(coxph(Surv(.time, .status) ~ .group, data = d), error = function(e) NULL)
  cox_global_p <- NA_real_
  hr <- NA_real_
  hr_low <- NA_real_
  hr_high <- NA_real_
  if (!is.null(cox)) {
    ss <- summary(cox)
    cox_global_p <- tryCatch(as.numeric(ss$logtest["pvalue"]), error = function(e) NA_real_)
    if (nlevels(d$.group) == 2 && nrow(ss$conf.int) >= 1) {
      hr <- as.numeric(ss$conf.int[1, "exp(coef)"])
      hr_low <- as.numeric(ss$conf.int[1, "lower .95"])
      hr_high <- as.numeric(ss$conf.int[1, "upper .95"])
    }
  }

  data.frame(
    GroupVariable = group_col,
    N = nrow(d),
    Events = sum(d$.status == 1),
    N_Groups = nlevels(d$.group),
    GroupCounts = paste(names(table(d$.group)), as.integer(table(d$.group)), sep = "=", collapse = "; "),
    Logrank_P = logrank_p,
    Cox_Global_P = cox_global_p,
    HR_Second_vs_First = hr,
    HR95_Low = hr_low,
    HR95_High = hr_high,
    stringsAsFactors = FALSE
  )
}

find_numeric_feature_cols <- function(df, exclude = character()) {
  setdiff(names(df)[sapply(df, function(x) is.numeric(x) || all(is.na(num(x)) | is.finite(num(x))))], exclude)
}

attach_anno_patient <- function(df, keep = c("Patient", "MO_DDRscore_group", "MO_DDRscore_raw", "NMF_subtype")) {
  anno_cols <- c("MO_DDRscore_group", "MO_DDRscore_raw", "MO_DDRscore", "NMF_subtype", "NMF_cluster_raw")
  df %>%
    select(-any_of(anno_cols)) %>%
    left_join(anno %>% select(all_of(keep)), by = "Patient")
}

attach_anno_sample <- function(df) {
  anno_cols <- c("Patient", "MO_DDRscore_group", "MO_DDRscore_raw", "MO_DDRscore", "NMF_subtype", "NMF_cluster_raw")
  df %>%
    select(-any_of(anno_cols)) %>%
    left_join(anno %>% select(Sample, Patient, MO_DDRscore_group, MO_DDRscore_raw, NMF_subtype), by = "Sample")
}

############################
# 3. Load core annotations
############################

score <- read_csv(SCORE_FILE) %>%
  mutate(
    Sample = norm_sample(Sample),
    Patient = patient_id(Patient),
    MO_DDRscore_raw = num(MO_DDRscore_raw),
    MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High"))
  ) %>%
  filter(SampleClass == "Tumor", !is.na(MO_DDRscore_group), is.finite(MO_DDRscore_raw)) %>%
  arrange(Patient, Sample) %>%
  distinct(Patient, .keep_all = TRUE)

nmf <- read_csv(NMF_FILE) %>%
  mutate(
    Sample = norm_sample(Sample),
    Patient = patient_id(Patient),
    MO_DDRscore_raw = num(MO_DDRscore_raw),
    MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High")),
    NMF_subtype = factor(NMF_subtype, levels = c("C1", "C2"))
  ) %>%
  filter(!is.na(NMF_subtype)) %>%
  distinct(Patient, .keep_all = TRUE)

anno <- score %>%
  select(Sample, Patient, MO_DDRscore_raw, MO_DDRscore_group) %>%
  left_join(nmf %>% select(Patient, NMF_cluster_raw, NMF_subtype), by = "Patient")

save_csv(anno, file.path(OUT_DIR, "Master_MO_DDRscore_NMF_annotation.csv"))

cat("Samples in MO-DDRscore annotation:", nrow(anno), "\n")
print(table(anno$MO_DDRscore_group, useNA = "ifany"))
print(table(anno$NMF_subtype, useNA = "ifany"))

############################
# 4. NMF vs MO-DDRscore concordance
############################

tab_nmf_score <- table(nmf$NMF_subtype, nmf$MO_DDRscore_group)
concordance <- data.frame(
  Item = c("N_samples", "Fisher_P", "ChiSquare_P", "Cramers_V"),
  Value = c(
    nrow(nmf),
    tryCatch(fisher.test(tab_nmf_score)$p.value, error = function(e) NA_real_),
    tryCatch(suppressWarnings(chisq.test(tab_nmf_score, correct = FALSE)$p.value), error = function(e) NA_real_),
    calc_cramers_v(tab_nmf_score)
  )
)

save_csv(as.data.frame.matrix(tab_nmf_score), file.path(OUT_DIR, "NMF_subtype_vs_MO_DDRscore_group_count_matrix.csv"))
save_csv(concordance, file.path(OUT_DIR, "NMF_subtype_vs_MO_DDRscore_group_test.csv"))

############################
# 5. NMF subtype immune deconvolution
############################

if (file.exists(DECONV_LONG_FILE)) {
  deconv_long <- read_csv(DECONV_LONG_FILE) %>%
    mutate(
      Sample = norm_sample(Sample),
      Patient = patient_id(Patient),
      Score = num(Score)
    ) %>%
    select(Sample, Patient, Method, Feature, Score) %>%
    attach_anno_patient() %>%
    filter(!is.na(NMF_subtype))

  deconv_nmf <- compare_long(
    deconv_long,
    group_col = "NMF_subtype",
    value_col = "Score",
    feature_cols = c("Method", "Feature"),
    level_order = c("C1", "C2")
  ) %>%
    arrange(FDR, P)

  save_csv(deconv_long, file.path(OUT_DIR, "NMF_immune_deconvolution_long.csv"))
  save_csv(deconv_nmf, file.path(OUT_DIR, "NMF_immune_deconvolution_group_comparison.csv"))
}

############################
# 6. NMF subtype official IOBR signature programs
############################

if (file.exists(IOBR_LONG_FILE)) {
  iobr_long <- read_csv(IOBR_LONG_FILE) %>%
    mutate(
      Sample = norm_sample(Sample),
      Patient = patient_id(Patient),
      Score = num(Score)
    ) %>%
    select(Sample, Patient, Signature, IOBR_group, N_genes_total, N_genes_matched, Score) %>%
    attach_anno_patient() %>%
    filter(!is.na(NMF_subtype))

  iobr_nmf <- compare_long(
    iobr_long,
    group_col = "NMF_subtype",
    value_col = "Score",
    feature_cols = c("Signature", "IOBR_group"),
    level_order = c("C1", "C2")
  ) %>%
    arrange(FDR, P)

  save_csv(iobr_long, file.path(OUT_DIR, "NMF_IOBR_signature_scores_long.csv"))
  save_csv(iobr_nmf, file.path(OUT_DIR, "NMF_IOBR_signature_group_comparison.csv"))
}

############################
# 7. NMF subtype TIDE / TMB / MATH
############################

if (file.exists(TIDE_FILE)) {
  tide <- read_csv(TIDE_FILE) %>%
    mutate(
      Sample = norm_sample(Sample),
      Patient = patient_id(Patient)
    ) %>%
    attach_anno_patient()

  tide_features <- intersect(
    c("TIDE", "IFNG", "MSI Score", "CD274", "CD8", "Dysfunction", "Exclusion", "MDSC", "CAF", "TAM M2", "CTL"),
    colnames(tide)
  )

  tide_long <- tide %>%
    select(Sample, Patient, MO_DDRscore_group, MO_DDRscore_raw, NMF_subtype, all_of(tide_features)) %>%
    pivot_longer(cols = all_of(tide_features), names_to = "Feature", values_to = "Score") %>%
    mutate(Score = num(Score)) %>%
    filter(!is.na(NMF_subtype), is.finite(Score))

  tide_nmf <- compare_long(
    tide_long,
    group_col = "NMF_subtype",
    value_col = "Score",
    feature_cols = c("Feature"),
    level_order = c("C1", "C2")
  ) %>%
    arrange(FDR, P)

  save_csv(tide_long, file.path(OUT_DIR, "NMF_TIDE_features_long.csv"))
  save_csv(tide_nmf, file.path(OUT_DIR, "NMF_TIDE_features_group_comparison.csv"))

  cat_results <- list()
  for (cat_col in intersect(c("Responder", "No benefits", "CTL.flag"), colnames(tide))) {
    res <- categorical_test(tide, "NMF_subtype", cat_col, level_order = c("C1", "C2"))
    if (length(res) > 0) {
      save_csv(res$counts, file.path(OUT_DIR, paste0("NMF_TIDE_", make.names(cat_col), "_counts.csv")))
      cat_results[[cat_col]] <- res$test
    }
  }
  if (length(cat_results) > 0) {
    save_csv(bind_rows(cat_results), file.path(OUT_DIR, "NMF_TIDE_categorical_tests.csv"))
  }
}

if (file.exists(TMB_FILE)) {
  tmb <- read_csv(TMB_FILE) %>%
    mutate(
      Sample = norm_sample(Sample),
      Patient = patient_id(Patient),
      TMB_value = num(TMB_value),
      TMB_log1p = log1p(TMB_value)
    ) %>%
    attach_anno_patient()

  tmb_nmf <- compare_long(
    tmb %>% transmute(Patient, NMF_subtype, Feature = "TMB_log1p", Score = TMB_log1p),
    group_col = "NMF_subtype",
    value_col = "Score",
    feature_cols = c("Feature"),
    level_order = c("C1", "C2")
  )

  save_csv(tmb, file.path(OUT_DIR, "NMF_TMB_merged.csv"))
  save_csv(tmb_nmf, file.path(OUT_DIR, "NMF_TMB_group_comparison.csv"))
}

if (file.exists(MATH_FILE)) {
  math <- read_csv(MATH_FILE) %>%
    mutate(
      Patient = patient_id(Patient),
      MATH = num(MATH)
    ) %>%
    attach_anno_patient()

  math_nmf <- compare_long(
    math %>% transmute(Patient, NMF_subtype, Feature = "MATH", Score = MATH),
    group_col = "NMF_subtype",
    value_col = "Score",
    feature_cols = c("Feature"),
    level_order = c("C1", "C2")
  )

  save_csv(math, file.path(OUT_DIR, "NMF_MATH_merged.csv"))
  save_csv(math_nmf, file.path(OUT_DIR, "NMF_MATH_group_comparison.csv"))
}

############################
# 8. Official checkpoint panel
############################

checkpoint_source <- "none"
checkpoint_genes <- character()

if (file.exists(EXT_CHECKPOINT_FILE)) {
  cp <- read_csv(EXT_CHECKPOINT_FILE)
  checkpoint_genes <- clean_gene(cp[[1]])
  checkpoint_genes <- unique(checkpoint_genes[!is.na(checkpoint_genes)])
  checkpoint_source <- "external official_checkpoint_gene_panel.csv"
} else if (requireNamespace("IOBR", quietly = TRUE)) {
  suppressPackageStartupMessages(library(IOBR))
  data(signature_collection, package = "IOBR")
  if (exists("signature_collection") && "Immune_Checkpoint" %in% names(signature_collection)) {
    checkpoint_genes <- clean_gene(signature_collection$Immune_Checkpoint)
    checkpoint_genes <- unique(checkpoint_genes[!is.na(checkpoint_genes)])
    checkpoint_source <- "IOBR::signature_collection$Immune_Checkpoint"
  }
}

checkpoint_meta <- data.frame(
  Source = checkpoint_source,
  N_genes = length(checkpoint_genes),
  Gene = paste(checkpoint_genes, collapse = ";"),
  stringsAsFactors = FALSE
)
save_csv(checkpoint_meta, file.path(OUT_DIR, "Official_checkpoint_panel_source.csv"))

if (length(checkpoint_genes) > 0 && file.exists(EXPR_FILE)) {
  expr <- readRDS(EXPR_FILE)
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  rownames(expr) <- clean_gene(rownames(expr))
  colnames(expr) <- norm_sample(colnames(expr))
  expr <- expr[!is.na(rownames(expr)) & rownames(expr) != "", , drop = FALSE]
  expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]

  matched_cp <- intersect(checkpoint_genes, rownames(expr))
  unmatched_cp <- setdiff(checkpoint_genes, matched_cp)

  cp_expr <- log2(expr[matched_cp, intersect(anno$Sample, colnames(expr)), drop = FALSE] + 1)
  cp_long <- as.data.frame(cp_expr) %>%
    tibble::rownames_to_column("Gene") %>%
    pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Expression") %>%
    mutate(Sample = norm_sample(Sample), Expression = num(Expression)) %>%
    attach_anno_sample()

  cp_mo <- compare_long(
    cp_long,
    group_col = "MO_DDRscore_group",
    value_col = "Expression",
    feature_cols = c("Gene"),
    level_order = c("Low", "High")
  ) %>%
    arrange(FDR, P)

  cp_nmf <- compare_long(
    cp_long,
    group_col = "NMF_subtype",
    value_col = "Expression",
    feature_cols = c("Gene"),
    level_order = c("C1", "C2")
  ) %>%
    arrange(FDR, P)

  save_csv(data.frame(Gene = matched_cp), file.path(OUT_DIR, "Official_checkpoint_genes_matched.csv"))
  save_csv(data.frame(Gene = unmatched_cp), file.path(OUT_DIR, "Official_checkpoint_genes_unmatched.csv"))
  save_csv(cp_long, file.path(OUT_DIR, "Official_checkpoint_expression_long.csv"))
  save_csv(cp_mo, file.path(OUT_DIR, "Official_checkpoint_MO_DDRscore_group_comparison.csv"))
  save_csv(cp_nmf, file.path(OUT_DIR, "Official_checkpoint_NMF_subtype_group_comparison.csv"))
}

############################
# 9. TMB / MATH survival linkage
############################

if (file.exists(CLIN_FILE)) {
  clin <- read_csv(CLIN_FILE) %>%
    mutate(
      Patient = patient_id(Patient),
      time = num(time),
      status = num(status)
    ) %>%
    filter(is.finite(time), time > 0, status %in% c(0, 1))

  surv_df <- clin %>%
    left_join(anno %>% select(Patient, MO_DDRscore_group, MO_DDRscore_raw, NMF_subtype), by = "Patient")

  if (exists("tmb")) {
    tmb_surv <- tmb %>%
      select(Patient, TMB_value, TMB_log1p) %>%
      distinct(Patient, .keep_all = TRUE) %>%
      mutate(
        TMB_group = ifelse(TMB_value >= median(TMB_value, na.rm = TRUE), "TMB-high", "TMB-low")
      )
    surv_df <- surv_df %>% left_join(tmb_surv, by = "Patient")
  }

  if (exists("math")) {
    math_surv <- math %>%
      select(Patient, MATH) %>%
      distinct(Patient, .keep_all = TRUE) %>%
      mutate(
        MATH_group = ifelse(MATH >= median(MATH, na.rm = TRUE), "MATH-high", "MATH-low")
      )
    surv_df <- surv_df %>% left_join(math_surv, by = "Patient")
  }

  surv_df <- surv_df %>%
    mutate(
      MO_TMB_group = ifelse(!is.na(MO_DDRscore_group) & !is.na(TMB_group),
                            paste(MO_DDRscore_group, TMB_group, sep = "_"), NA),
      NMF_TMB_group = ifelse(!is.na(NMF_subtype) & !is.na(TMB_group),
                             paste(NMF_subtype, TMB_group, sep = "_"), NA),
      MO_MATH_group = ifelse(!is.na(MO_DDRscore_group) & !is.na(MATH_group),
                             paste(MO_DDRscore_group, MATH_group, sep = "_"), NA),
      NMF_MATH_group = ifelse(!is.na(NMF_subtype) & !is.na(MATH_group),
                              paste(NMF_subtype, MATH_group, sep = "_"), NA)
    )

  surv_vars <- intersect(
    c("MO_DDRscore_group", "NMF_subtype", "TMB_group", "MATH_group",
      "MO_TMB_group", "NMF_TMB_group", "MO_MATH_group", "NMF_MATH_group"),
    colnames(surv_df)
  )

  surv_results <- bind_rows(lapply(surv_vars, function(v) survival_test(surv_df, v)))

  save_csv(surv_df, file.path(OUT_DIR, "Survival_MO_NMF_TMB_MATH_merged.csv"))
  save_csv(surv_results, file.path(OUT_DIR, "Survival_MO_NMF_TMB_MATH_tests.csv"))
}

############################
# 10. Optional official external modules
############################

external_status <- data.frame(
  Module = c(
    "TCGA pan-cancer immune subtype C1-C6",
    "TCIA IPS PD1/CTLA4 modes",
    "Tumor neoantigen burden TNB",
    "MSI / stemness biomarkers",
    "External immunotherapy response cohort",
    "Expanded official checkpoint panel"
  ),
  Expected_File = c(
    EXT_IMMUNE_SUBTYPE_FILE,
    EXT_TCIA_IPS_FILE,
    EXT_TNB_FILE,
    EXT_MSI_STEMNESS_FILE,
    EXT_IMMUNOTHERAPY_FILE,
    EXT_CHECKPOINT_FILE
  ),
  Required_Minimum_Columns = c(
    "Patient, Immune_Subtype",
    "Patient, IPS-like numeric columns",
    "Patient, TNB",
    "Patient, MSI/stemness numeric or categorical columns",
    "Patient or Sample, response column, optional time/status",
    "Gene"
  ),
  Status = c(
    ifelse(file.exists(EXT_IMMUNE_SUBTYPE_FILE), "FOUND", "MISSING"),
    ifelse(file.exists(EXT_TCIA_IPS_FILE), "FOUND", "MISSING"),
    ifelse(file.exists(EXT_TNB_FILE), "FOUND", "MISSING"),
    ifelse(file.exists(EXT_MSI_STEMNESS_FILE), "FOUND", "MISSING"),
    ifelse(file.exists(EXT_IMMUNOTHERAPY_FILE), "FOUND", "MISSING"),
    ifelse(file.exists(EXT_CHECKPOINT_FILE), "FOUND", "MISSING; using IOBR checkpoint if available")
  ),
  stringsAsFactors = FALSE
)

save_csv(external_status, file.path(OUT_DIR, "README_required_external_files.csv"))

run_external_compare <- function(file, module_name, group_vars = c("MO_DDRscore_group", "NMF_subtype")) {
  x <- read_csv(file)
  x <- standardize_patient_col(x)
  x <- attach_anno_patient(x)

  num_cols <- names(x)[sapply(x, function(z) {
    zz <- num(z)
    sum(is.finite(zz)) >= 20
  })]
  num_cols <- setdiff(num_cols, c("MO_DDRscore_raw"))

  if (length(num_cols) > 0) {
    long <- x %>%
      select(Patient, all_of(group_vars[group_vars %in% colnames(x)]), all_of(num_cols)) %>%
      pivot_longer(cols = all_of(num_cols), names_to = "Feature", values_to = "Score") %>%
      mutate(Score = num(Score))

    for (gv in group_vars[group_vars %in% colnames(long)]) {
      lv <- if (gv == "MO_DDRscore_group") c("Low", "High") else if (gv == "NMF_subtype") c("C1", "C2") else NULL
      cmp <- compare_long(long, gv, "Score", c("Feature"), level_order = lv)
      save_csv(cmp, file.path(OUT_DIR, paste0(module_name, "_", gv, "_numeric_comparison.csv")))
    }
  }

  cat_cols <- setdiff(names(x), c("Patient", "Sample", "MO_DDRscore_raw", num_cols))
  cat_cols <- cat_cols[sapply(x[cat_cols], function(z) length(unique(na.omit(z))) >= 2 && length(unique(na.omit(z))) <= 12)]
  for (cc in cat_cols) {
    for (gv in group_vars[group_vars %in% colnames(x)]) {
      lv <- if (gv == "MO_DDRscore_group") c("Low", "High") else if (gv == "NMF_subtype") c("C1", "C2") else NULL
      res <- categorical_test(x, gv, cc, level_order = lv)
      if (length(res) > 0) {
        base <- paste0(module_name, "_", gv, "_", make.names(cc))
        save_csv(res$counts, file.path(OUT_DIR, paste0(base, "_counts.csv")))
        save_csv(res$test, file.path(OUT_DIR, paste0(base, "_test.csv")))
      }
    }
  }

  save_csv(x, file.path(OUT_DIR, paste0(module_name, "_merged.csv")))
}

if (file.exists(EXT_IMMUNE_SUBTYPE_FILE)) {
  run_external_compare(EXT_IMMUNE_SUBTYPE_FILE, "External_immune_subtype")
}

if (file.exists(EXT_TCIA_IPS_FILE)) {
  run_external_compare(EXT_TCIA_IPS_FILE, "External_TCIA_IPS")
}

if (file.exists(EXT_TNB_FILE)) {
  run_external_compare(EXT_TNB_FILE, "External_TNB")

  if (file.exists(CLIN_FILE)) {
    tnb <- standardize_patient_col(read_csv(EXT_TNB_FILE))
    tnb_col <- names(tnb)[grepl("^TNB$|neoantigen|neo_antigen", names(tnb), ignore.case = TRUE)][1]
    if (!is.na(tnb_col)) {
      tnb_surv <- tnb %>%
        transmute(Patient, TNB = num(.data[[tnb_col]])) %>%
        filter(is.finite(TNB)) %>%
        mutate(TNB_group = ifelse(TNB >= median(TNB, na.rm = TRUE), "TNB-high", "TNB-low"))
      clin_tnb <- read_csv(CLIN_FILE) %>%
        mutate(Patient = patient_id(Patient), time = num(time), status = num(status)) %>%
        attach_anno_patient(keep = c("Patient", "MO_DDRscore_group", "NMF_subtype")) %>%
        left_join(tnb_surv, by = "Patient") %>%
        mutate(
          MO_TNB_group = ifelse(!is.na(MO_DDRscore_group) & !is.na(TNB_group),
                                paste(MO_DDRscore_group, TNB_group, sep = "_"), NA),
          NMF_TNB_group = ifelse(!is.na(NMF_subtype) & !is.na(TNB_group),
                                 paste(NMF_subtype, TNB_group, sep = "_"), NA)
        )
      tnb_surv_res <- bind_rows(lapply(
        intersect(c("TNB_group", "MO_TNB_group", "NMF_TNB_group"), colnames(clin_tnb)),
        function(v) survival_test(clin_tnb, v)
      ))
      save_csv(clin_tnb, file.path(OUT_DIR, "External_TNB_survival_merged.csv"))
      save_csv(tnb_surv_res, file.path(OUT_DIR, "External_TNB_survival_tests.csv"))
    }
  }
}

if (file.exists(EXT_MSI_STEMNESS_FILE)) {
  run_external_compare(EXT_MSI_STEMNESS_FILE, "External_MSI_stemness")
}

if (file.exists(EXT_IMMUNOTHERAPY_FILE)) {
  run_external_compare(EXT_IMMUNOTHERAPY_FILE, "External_immunotherapy_response")
}

############################
# 11. Summary
############################

files_now <- list.files(OUT_DIR, full.names = FALSE)
summary_df <- data.frame(
  Item = c(
    "N_MO_DDRscore_samples",
    "N_NMF_samples",
    "N_output_files",
    "External_files_found",
    "External_files_missing"
  ),
  Value = c(
    nrow(anno),
    sum(!is.na(anno$NMF_subtype)),
    length(files_now),
    sum(external_status$Status == "FOUND"),
    sum(grepl("MISSING", external_status$Status))
  ),
  stringsAsFactors = FALSE
)

save_csv(summary_df, file.path(OUT_DIR, "Allinone_reference_gap_summary.csv"))

cat("\nDone.\n")
cat("Output directory:\n", OUT_DIR, "\n")
print(summary_df)
cat("\nExternal file status:\n")
print(external_status)
