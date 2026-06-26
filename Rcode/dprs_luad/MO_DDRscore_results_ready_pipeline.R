
############################################################
# MO-DDRscore clean pipeline
# Multi-omics Context-Adaptive DNA Damage Repair Score
#
# Stable version:
# 1) Build gene-level DDR weights from TCGA-LUAD multi-omics context
# 2) Entropy-based adaptive evidence fusion
# 3) Calculate sample-level MO-DDRscore
# 4) Fig2: score construction / baseline comparison / SPIDR external validation
# 5) Fig3: DEG, pathway, mutation, CNV, immune, checkpoint, survival
#
# SPIDR is NOT used to construct the score.
# SPIDR is only used as an independent external functional-dependency validation.
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260510)

############################
# 0. Paths and parameters
############################

PROJECT_DIR <- "D:/R_workspace/评分/AD_DDR_project"
DATA_DIR <- file.path(PROJECT_DIR, "00_data")
PROC_DIR <- file.path(PROJECT_DIR, "01_processed")
FIG2_DIR <- file.path(PROJECT_DIR, "02_Figure2_score_validation")
FIG3_DIR <- file.path(PROJECT_DIR, "03_Figure3_multiomics")
DB_DIR <- file.path(PROJECT_DIR, "05_database_tables")

dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PROC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG2_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG3_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DB_DIR, recursive = TRUE, showWarnings = FALSE)

FOCUS_CANCER <- "LUAD"
TCGA_DIR <- "D:/R/R_workspace/梁老师文件/TCGA"

TCGA_EXPR_FILE <- file.path(TCGA_DIR, "mRNA_exp_TPM_only_TCGA/mRNA_exp_TPM_only_TCGA", paste0("TCGA-", FOCUS_CANCER, ".gene_expression_TPM.tsv"))
TCGA_CLIN_FILE <- file.path(TCGA_DIR, "clinical", paste0("TCGA.", FOCUS_CANCER, ".sampleMap"), paste0(FOCUS_CANCER, "_clinicalMatrix"))
TCGA_MUT_FILE  <- file.path(TCGA_DIR, "mutation", paste0(FOCUS_CANCER, ".txt"))
TCGA_CNV_FILE  <- file.path(TCGA_DIR, "cnv", paste0(FOCUS_CANCER, ".txt"))

DDR_GENE_FILE <- file.path(DATA_DIR, "DDR_236_genes.csv")
PPI_FILE <- file.path(DATA_DIR, "string_interactions_short_0.7.tsv")
SPIDR_FILE <- file.path(DATA_DIR, "spidr_supp_table3_clean.csv")

DEPMAP_EXPR_FILE <- file.path(DATA_DIR, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
DEPMAP_MODEL_FILE <- file.path(DATA_DIR, "Model.csv")
GDSC_FILE <- file.path(DATA_DIR, "GDSC2_fitted_dose_response_27Oct23.xlsx")

DEG_LOGFC <- 1.20  # 标准差异阈值：|log2FC| > 1.2
DEG_FDR <- 0.05
MAX_HEATMAP_GENES <- 50

# 分组策略：
# median = 推荐主分析，避免用生存结局选择cutoff导致偏倚
# surv_cutpoint = 仅建议作为敏感性分析/补充图，不建议用于多组学主分组
SCORE_GROUP_METHOD <- "median"

############################
# 1. Packages and helpers
############################

for (p in c("data.table", "dplyr", "tidyr", "ggplot2", "survival")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(survival)
})

has_pkg <- function(p) requireNamespace(p, quietly = TRUE)

safe_fread <- function(file, ...) {
  if (!file.exists(file)) stop("File not found: ", file)
  data.table::fread(file, data.table = FALSE, check.names = FALSE, ...)
}

save_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file)
}

safe_ggsave <- function(file, p, w = 7, h = 5) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file, p, width = w, height = h, useDingbats = FALSE)
}

clean_gene <- function(x) {
  x <- as.character(x)
  x <- gsub(" ///.*$", "", x)
  x <- gsub(" //.*$", "", x)
  x <- gsub(";.*$", "", x)
  x <- gsub(",.*$", "", x)
  x <- gsub("\\s*\\([^\\)]*\\)$", "", x)
  x <- toupper(trimws(x))
  x[x %in% c("", "NA", "---", "NULL", "N/A", "NAN")] <- NA
  x
}

clean_na <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "[Not Available]", "[Not Applicable]", "not reported", "Not Reported", "null", "NULL", "--")] <- NA
  x
}

to_num <- function(x) suppressWarnings(as.numeric(clean_na(x)))
patient_id <- function(x) substr(gsub("\\.", "-", as.character(x)), 1, 12)
sample_type <- function(x) substr(gsub("\\.", "-", as.character(x)), 14, 15)
sample_class <- function(x) {
  code <- sample_type(x)
  ifelse(code %in% c("01","02","03","05","06","07"), "Tumor",
         ifelse(code %in% c("10","11","12","13","14"), "Normal", "Other"))
}

find_col <- function(df, patterns) {
  cn <- colnames(df)
  for (pa in patterns) {
    hit <- grep(pa, cn, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

minmax01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) return(rep(0, length(x)))
  rg <- range(x, na.rm = TRUE)
  if (!all(is.finite(rg)) || rg[1] == rg[2]) return(rep(0, length(x)))
  (x - rg[1]) / (rg[2] - rg[1])
}

z_rows <- function(m) {
  m <- as.matrix(m); storage.mode(m) <- "numeric"
  mu <- rowMeans(m, na.rm = TRUE)
  sdv <- apply(m, 1, sd, na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- 1
  sweep(sweep(m, 1, mu, "-"), 1, sdv, "/")
}

z_cols <- function(m) {
  m <- as.matrix(m); storage.mode(m) <- "numeric"
  mu <- colMeans(m, na.rm = TRUE)
  sdv <- apply(m, 2, sd, na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- 1
  sweep(sweep(m, 2, mu, "-"), 2, sdv, "/")
}

safe_cor <- function(x, y, method = "spearman") {
  x <- suppressWarnings(as.numeric(x)); y <- suppressWarnings(as.numeric(y))
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10) return(c(cor = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method))
  c(cor = unname(ct$estimate), p = ct$p.value)
}

calc_cindex <- function(time, status, score) {
  df <- data.frame(time = as.numeric(time), status = as.numeric(status), score = as.numeric(score))
  df <- df[is.finite(df$time) & is.finite(df$status) & is.finite(df$score), ]
  if (nrow(df) < 20 || sum(df$status == 1, na.rm = TRUE) < 5) return(NA_real_)
  suppressWarnings(tryCatch(survival::concordance(Surv(time, status) ~ score, data = df)$concordance, error = function(e) NA_real_))
}


make_score_group <- function(score_df, clin_df = NULL,
                             score_col = "AD_DDRscore",
                             method = "median",
                             out_file = NULL) {
  df <- score_df
  df$AD_DDRscore_group <- NA_character_

  tumor_idx <- which(df$SampleClass == "Tumor" & is.finite(df[[score_col]]))

  if (method == "surv_cutpoint") {
    if (!has_pkg("survminer")) {
      warning("survminer is not installed; fallback to median cutoff.")
      method <- "median"
    } else if (is.null(clin_df)) {
      warning("clin_df is NULL; fallback to median cutoff.")
      method <- "median"
    } else {
      tmp <- df[tumor_idx, c("Sample", "Patient", score_col), drop = FALSE] %>%
        dplyr::inner_join(clin_df[, c("Patient", "time", "status")], by = "Patient")
      tmp <- tmp[is.finite(tmp[[score_col]]) & is.finite(tmp$time) & tmp$status %in% c(0, 1), ]

      if (nrow(tmp) < 50 || sum(tmp$status == 1, na.rm = TRUE) < 10) {
        warning("Too few survival samples/events for surv_cutpoint; fallback to median cutoff.")
        method <- "median"
      } else {
        colnames(tmp)[colnames(tmp) == score_col] <- "score_tmp"
        cut_obj <- survminer::surv_cutpoint(
          tmp,
          time = "time",
          event = "status",
          variables = "score_tmp",
          minprop = 0.30
        )
        cut_value <- cut_obj$cutpoint$cutpoint[1]
        df$AD_DDRscore_group[tumor_idx] <- ifelse(
          df[[score_col]][tumor_idx] >= cut_value,
          "High", "Low"
        )
        if (!is.null(out_file)) {
          save_csv(
            data.frame(
              Method = "surv_cutpoint",
              Cutoff = cut_value,
              N = nrow(tmp),
              Events = sum(tmp$status == 1, na.rm = TRUE)
            ),
            out_file
          )
        }
      }
    }
  }

  if (method == "median") {
    cut_value <- stats::median(df[[score_col]][tumor_idx], na.rm = TRUE)
    df$AD_DDRscore_group[tumor_idx] <- ifelse(
      df[[score_col]][tumor_idx] >= cut_value,
      "High", "Low"
    )
    if (!is.null(out_file)) {
      save_csv(
        data.frame(
          Method = "median",
          Cutoff = cut_value,
          N = length(tumor_idx),
          Events = NA_integer_
        ),
        out_file
      )
    }
  }

  df$AD_DDRscore_group <- factor(df$AD_DDRscore_group, levels = c("Low", "High"))
  df
}

run_standard_survival <- function(df,
                                  score_col = "AD_DDRscore",
                                  group_col = "AD_DDRscore_group",
                                  out_prefix,
                                  title = "MO-DDRscore survival") {
  df <- df[
    is.finite(df[[score_col]]) &
      is.finite(df$time) &
      df$status %in% c(0, 1) &
      !is.na(df[[group_col]]),
  ]

  if (nrow(df) < 30 || sum(df$status == 1, na.rm = TRUE) < 5) {
    message("Survival skipped: insufficient samples/events.")
    return(invisible(NULL))
  }

  df[[group_col]] <- factor(df[[group_col]], levels = c("Low", "High"))

  # Continuous Cox
  cox_cont <- tryCatch(
    summary(survival::coxph(survival::Surv(time, status) ~ df[[score_col]], data = df)),
    error = function(e) NULL
  )

  # Group Cox
  cox_group <- tryCatch(
    summary(survival::coxph(survival::Surv(time, status) ~ df[[group_col]], data = df)),
    error = function(e) NULL
  )

  surv_table <- list()

  if (!is.null(cox_cont)) {
    surv_table[[length(surv_table) + 1]] <- data.frame(
      Model = "Continuous_score",
      HR = cox_cont$coefficients[1, "exp(coef)"],
      P = cox_cont$coefficients[1, "Pr(>|z|)"],
      CI_low = cox_cont$conf.int[1, "lower .95"],
      CI_high = cox_cont$conf.int[1, "upper .95"]
    )
  }

  if (!is.null(cox_group)) {
    surv_table[[length(surv_table) + 1]] <- data.frame(
      Model = "High_vs_Low",
      HR = cox_group$coefficients[1, "exp(coef)"],
      P = cox_group$coefficients[1, "Pr(>|z|)"],
      CI_low = cox_group$conf.int[1, "lower .95"],
      CI_high = cox_group$conf.int[1, "upper .95"]
    )
  }

  if (length(surv_table) > 0) {
    save_csv(dplyr::bind_rows(surv_table), paste0(out_prefix, "_cox.csv"))
  }

  # KM with survminer
  if (has_pkg("survminer")) {
    fit <- survival::survfit(
      stats::as.formula(paste0("survival::Surv(time, status) ~ ", group_col)),
      data = df
    )
    p <- survminer::ggsurvplot(
      fit,
      data = df,
      pval = TRUE,
      conf.int = FALSE,
      risk.table = TRUE,
      risk.table.height = 0.22,
      palette = c("#2B6CB0", "#C53030"),
      legend.title = "",
      legend.labs = c("Low", "High"),
      title = title,
      ggtheme = ggplot2::theme_bw()
    )
    ggplot2::ggsave(
      paste0(out_prefix, "_KM.pdf"),
      p$plot,
      width = 6,
      height = 5,
      useDingbats = FALSE
    )
    ggplot2::ggsave(
      paste0(out_prefix, "_KM_with_risktable.pdf"),
      print(p),
      width = 6,
      height = 6.2,
      useDingbats = FALSE
    )
  }

  # Time-dependent ROC if available
  if (has_pkg("timeROC")) {
    roc_obj <- tryCatch(
      timeROC::timeROC(
        T = df$time,
        delta = df$status,
        marker = df[[score_col]],
        cause = 1,
        times = c(365, 1095, 1825),
        iid = TRUE
      ),
      error = function(e) NULL
    )

    if (!is.null(roc_obj)) {
      roc_df <- data.frame(
        Time = c("1-year", "3-year", "5-year"),
        AUC = as.numeric(roc_obj$AUC)
      )
      save_csv(roc_df, paste0(out_prefix, "_timeROC_AUC.csv"))

      grDevices::pdf(paste0(out_prefix, "_timeROC.pdf"), width = 5.5, height = 5)
      plot(roc_obj, time = 365, col = "#2B6CB0", title = FALSE)
      plot(roc_obj, time = 1095, add = TRUE, col = "#D69E2E")
      plot(roc_obj, time = 1825, add = TRUE, col = "#C53030")
      legend(
        "bottomright",
        legend = paste0(roc_df$Time, " AUC=", sprintf("%.3f", roc_df$AUC)),
        col = c("#2B6CB0", "#D69E2E", "#C53030"),
        lwd = 2,
        bty = "n"
      )
      title(title)
      grDevices::dev.off()
    }
  }

  invisible(NULL)
}


plot_km <- function(df, score_col, title, file) {
  if (!has_pkg("survminer")) return(invisible(NULL))
  df <- df[is.finite(df[[score_col]]) & is.finite(df$time) & is.finite(df$status), ]
  if (nrow(df) < 30 || sum(df$status == 1) < 5) return(invisible(NULL))
  df$RiskGroup <- ifelse(df[[score_col]] >= median(df[[score_col]], na.rm = TRUE), "High", "Low")
  df$RiskGroup <- factor(df$RiskGroup, levels = c("Low", "High"))
  fit <- survival::survfit(Surv(time, status) ~ RiskGroup, data = df)
  p <- survminer::ggsurvplot(
    fit, data = df, pval = TRUE, risk.table = FALSE,
    palette = c("#2B6CB0", "#C53030"),
    legend.title = "", legend.labs = c("Low", "High"),
    title = title, ggtheme = theme_bw()
  )$plot
  safe_ggsave(file, p, 6, 5)
}

entropy_fusion <- function(df, id_col = "Gene") {
  feats <- setdiff(colnames(df), id_col)
  mat <- as.matrix(df[, feats, drop = FALSE])
  storage.mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0
  mat_norm <- apply(mat, 2, minmax01)
  if (is.null(dim(mat_norm))) mat_norm <- matrix(mat_norm, ncol = 1)
  colnames(mat_norm) <- feats
  rownames(mat_norm) <- df[[id_col]]
  eps <- 1e-12
  p <- sweep(mat_norm + eps, 2, colSums(mat_norm + eps), "/")
  p[!is.finite(p)] <- 0
  k <- 1 / log(nrow(p))
  ent <- -k * colSums(p * log(p + eps), na.rm = TRUE)
  div <- 1 - ent
  div[!is.finite(div) | div < 0] <- 0
  if (sum(div) == 0) {
    w <- rep(1 / length(feats), length(feats))
  } else {
    w <- div / sum(div)
  }
  names(w) <- feats
  list(weights = w, normalized = mat_norm, entropy = ent, divergence = div)
}

############################
# 2. DDR genes
############################

fallback_genes <- c(
  "ABRAXAS1","ALKBH1","ALKBH2","ALKBH3","APEX1","APEX2","APLF","APTX","ATAD5",
  "ATM","ATR","ATRIP","ATRX","AUNIP","BABAM1","BABAM2","BARD1","BLM","BOD1L1",
  "BRAT1","BRCA1","BRCA2","BRCC3","BRIP1","CCNH","CDK7","CETN2","CGAS","CHAF1A",
  "CHEK1","CHEK2","CLSPN","DCLRE1A","DCLRE1B","DCLRE1C","DDB1","DDB2","DMC1",
  "DNA2","DNPH1","DNTT","DUT","EME1","EME2","ENDOV","ERCC1","ERCC2","ERCC3",
  "ERCC4","ERCC5","ERCC6","ERCC8","EXO1","EXO5","FAAP100","FAAP20","FAAP24",
  "FAN1","FANCA","FANCB","FANCC","FANCD2","FANCE","FANCF","FANCG","FANCI",
  "FANCL","FANCM","GEN1","GTF2H1","GTF2H2","GTF2H3","GTF2H4","GTF2H5","H2AX",
  "HMGB1","HUS1","LIG1","LIG3","LIG4","MAD2L2","MBD4","MGMT","MLH1","MLH3",
  "MMS22L","MRE11","MSH2","MSH3","MSH4","MSH5","MSH6","MUS81","NABP1","NABP2",
  "NBN","NEIL1","NEIL2","NEIL3","NHEJ1","NTHL1","OGG1","PALB2","PARP1","PARP2",
  "PARP3","PARP4","PCNA","PMS1","PMS2","POLD1","POLD2","POLD3","POLD4","POLE",
  "POLE2","POLE3","POLE4","POLB","POLDIP2","POLG","POLH","POLI","POLK","POLL",
  "POLM","POLN","POLQ","PRKDC","RAD1","RAD9A","RAD17","RAD18","RAD23A","RAD23B",
  "RAD50","RAD51","RAD51B","RAD51C","RAD51D","RAD52","RAD54B","RAD54L","RBBP8",
  "RECQL","RECQL4","RECQL5","REV1","REV3L","RFC1","RFC2","RFC3","RFC4","RFC5",
  "RIF1","RNF4","RNF8","RNF168","RPA1","RPA2","RPA3","RPA4","RRM2B","SETMAR",
  "SHLD1","SHLD2","SHLD3","SLX1A","SLX1B","SLX4","SMC5","SMC6","SMUG1","SPRTN",
  "SWI5","SWSAP1","TDG","TDP1","TDP2","TELO2","TOP3A","TOPBP1","TP53","TP53BP1",
  "TREX1","TREX2","UBE2A","UBE2B","UBE2N","UBE2T","UBE2V2","UNG","USP1","UVSSA",
  "WDR48","WRN","XAB2","XPA","XPC","XRCC1","XRCC2","XRCC3","XRCC4","XRCC5",
  "XRCC6","ZSWIM7","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MCM8","MCM9",
  "CDC45","CDC6","CDK1","CDK2","CCNA2","CCNB1","CCNE1","TYMS","DTL","WDHD1",
  "HELLS","CHAF1B","SMC2","NCAPG2","BUB1B","KIF11","KIF18A","KIF20B"
)

if (file.exists(DDR_GENE_FILE)) {
  tmp <- safe_fread(DDR_GENE_FILE)
  gc <- find_col(tmp, c("^Gene$", "gene", "symbol"))
  if (is.na(gc)) gc <- colnames(tmp)[1]
  DDR_GENES <- unique(na.omit(clean_gene(tmp[[gc]])))
} else {
  DDR_GENES <- unique(na.omit(clean_gene(fallback_genes)))
}
save_csv(data.frame(Gene = DDR_GENES), file.path(PROC_DIR, "DDR_genes_used.csv"))

############################
# 3. Read TCGA data
############################

read_expr <- function(file) {
  x <- safe_fread(file)
  gc <- find_col(x, c("^Gene$", "gene", "symbol", "Gene Symbol", "gene_name"))
  if (is.na(gc)) gc <- colnames(x)[1]
  genes <- clean_gene(x[[gc]])
  matdf <- x[, setdiff(colnames(x), gc), drop = FALSE]
  matdf$Gene <- genes
  matdf <- matdf[!is.na(matdf$Gene) & matdf$Gene != "", ]
  matdf <- matdf %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ suppressWarnings(mean(as.numeric(.x), na.rm = TRUE))), .groups = "drop")
  rn <- matdf$Gene
  matdf$Gene <- NULL
  mat <- as.matrix(as.data.frame(matdf, check.names = FALSE))
  rownames(mat) <- rn
  storage.mode(mat) <- "numeric"
  mat
}

read_clin <- function(file) {
  x <- safe_fread(file)
  pc <- find_col(x, c("^_PATIENT$", "bcr_patient_barcode", "patient_id", "Patient", "sampleID"))
  if (is.na(pc)) pc <- colnames(x)[1]
  vc <- find_col(x, c("^vital_status$", "vital.status", "death", "dead"))
  dc <- find_col(x, c("^days_to_death$", "days.to.death"))
  fc <- find_col(x, c("^days_to_last_followup$", "days_to_last_follow_up", "days.to.last.followup", "days_to_followup"))
  if (is.na(vc) || (is.na(dc) && is.na(fc))) stop("Cannot infer clinical survival columns.")
  vital <- toupper(clean_na(x[[vc]]))
  death <- if (!is.na(dc)) to_num(x[[dc]]) else rep(NA_real_, nrow(x))
  follow <- if (!is.na(fc)) to_num(x[[fc]]) else rep(NA_real_, nrow(x))
  status <- ifelse(vital %in% c("DECEASED", "DEAD", "DIED"), 1, ifelse(vital %in% c("LIVING", "ALIVE"), 0, NA_real_))
  time <- ifelse(status == 1 & is.finite(death) & death > 0, death, follow)
  agec <- find_col(x, c("age_at_initial_pathologic_diagnosis", "^age$", "age_at_diagnosis"))
  genc <- find_col(x, c("^gender$", "sex"))
  stgc <- find_col(x, c("pathologic_stage", "tumor_stage", "^stage$"))
  out <- data.frame(
    Patient = patient_id(x[[pc]]),
    time = as.numeric(time),
    status = as.numeric(status),
    age = if (!is.na(agec)) to_num(x[[agec]]) else NA_real_,
    gender = if (!is.na(genc)) as.character(x[[genc]]) else NA_character_,
    stage = if (!is.na(stgc)) as.character(x[[stgc]]) else NA_character_
  )
  out <- out[!is.na(out$Patient) & out$Patient != "" & is.finite(out$time) & out$time > 0 & out$status %in% c(0, 1), ]
  out %>% dplyr::arrange(Patient, dplyr::desc(status), dplyr::desc(time)) %>% dplyr::distinct(Patient, .keep_all = TRUE)
}

read_mut <- function(file) {
  if (!file.exists(file)) return(NULL)
  x <- tryCatch(safe_fread(file), error = function(e) NULL)
  if (is.null(x) || nrow(x) == 0) return(NULL)
  gc <- find_col(x, c("^Hugo_Symbol$", "^Gene$", "gene", "symbol"))
  sc <- find_col(x, c("^Tumor_Sample_Barcode$", "sample", "barcode", "Tumor"))
  if (is.na(gc) || is.na(sc)) return(NULL)
  x$Hugo_Symbol <- clean_gene(x[[gc]])
  x$Tumor_Sample_Barcode <- as.character(x[[sc]])
  x$Patient <- patient_id(x$Tumor_Sample_Barcode)
  if (!"Variant_Classification" %in% colnames(x)) x$Variant_Classification <- "Mutation"
  if (!"Variant_Type" %in% colnames(x)) x$Variant_Type <- "SNP"
  if (!"Start_Position" %in% colnames(x)) x$Start_Position <- 1
  if (!"End_Position" %in% colnames(x)) x$End_Position <- x$Start_Position
  if (!"Reference_Allele" %in% colnames(x)) x$Reference_Allele <- "N"
  if (!"Tumor_Seq_Allele2" %in% colnames(x)) x$Tumor_Seq_Allele2 <- "N"
  x
}

read_cnv <- function(file) {
  if (!file.exists(file)) return(NULL)
  x <- tryCatch(safe_fread(file), error = function(e) NULL)
  if (is.null(x) || nrow(x) == 0) return(NULL)
  gc <- find_col(x, c("Gene Symbol", "^Gene$", "gene", "symbol"))
  if (is.na(gc)) gc <- colnames(x)[1]
  genes <- clean_gene(x[[gc]])
  matdf <- x[, setdiff(colnames(x), gc), drop = FALSE]
  matdf$Gene <- genes
  matdf <- matdf[!is.na(matdf$Gene) & matdf$Gene != "", ]
  matdf <- matdf %>% dplyr::group_by(Gene) %>% dplyr::summarise(dplyr::across(dplyr::everything(), ~ suppressWarnings(mean(as.numeric(.x), na.rm = TRUE))), .groups = "drop")
  rn <- matdf$Gene
  matdf$Gene <- NULL
  mat <- as.matrix(as.data.frame(matdf, check.names = FALSE))
  rownames(mat) <- rn
  storage.mode(mat) <- "numeric"
  mat
}

cat("\nLoading data...\n")
tcga_expr <- read_expr(TCGA_EXPR_FILE)
tumor_samples <- colnames(tcga_expr)[sample_class(colnames(tcga_expr)) == "Tumor"]
normal_samples <- colnames(tcga_expr)[sample_class(colnames(tcga_expr)) == "Normal"]
tcga_clin <- read_clin(TCGA_CLIN_FILE)
mut <- read_mut(TCGA_MUT_FILE)
cnv_mat <- read_cnv(TCGA_CNV_FILE)

cat("Expression:", nrow(tcga_expr), "genes x", ncol(tcga_expr), "samples\n")
cat("Tumor:", length(tumor_samples), " Normal:", length(normal_samples), "\n")
cat("Clinical:", nrow(tcga_clin), " Events:", sum(tcga_clin$status == 1), "\n")
save_csv(tcga_clin, file.path(PROC_DIR, paste0(FOCUS_CANCER, "_clinical_processed.csv")))

############################
# 4. MO-DDRweight construction
############################

cat("\nConstructing MO-DDRweight...\n")
ddr_genes <- intersect(DDR_GENES, rownames(tcga_expr))
if (length(ddr_genes) < 20) stop("Too few DDR genes in expression matrix.")

# E1: expression dysregulation
expr_log <- log2(tcga_expr[ddr_genes, , drop = FALSE] + 1)
if (length(normal_samples) >= 5) {
  tumor_mean <- rowMeans(expr_log[, tumor_samples, drop = FALSE], na.rm = TRUE)
  normal_mean <- rowMeans(expr_log[, normal_samples, drop = FALSE], na.rm = TRUE)
  logFC <- tumor_mean - normal_mean
  pval <- sapply(ddr_genes, function(g) {
    tryCatch(wilcox.test(expr_log[g, tumor_samples], expr_log[g, normal_samples])$p.value, error = function(e) NA_real_)
  })
  fdr <- p.adjust(pval, "BH")
  ExprDysregulation <- abs(logFC) * (-log10(fdr + 1e-300))
} else {
  logFC <- rep(NA_real_, length(ddr_genes)); pval <- fdr <- rep(NA_real_, length(ddr_genes))
  tumor_expr <- expr_log[, tumor_samples, drop = FALSE]
  ExprDysregulation <- minmax01(rowMeans(tumor_expr, na.rm = TRUE)) * minmax01(apply(tumor_expr, 1, sd, na.rm = TRUE))
}
expr_ev <- data.frame(Gene = ddr_genes, TumorNormal_logFC = logFC, TumorNormal_P = pval, TumorNormal_FDR = fdr, ExprDysregulation = ExprDysregulation)

# E2: genomic alteration
gen_ev <- data.frame(Gene = ddr_genes, MutationFreq = 0, AmpFreq = 0, DelFreq = 0, CNVAlterFreq = 0)
if (!is.null(mut)) {
  pats <- unique(patient_id(tumor_samples))
  mf <- mut %>% dplyr::filter(Hugo_Symbol %in% ddr_genes, Patient %in% pats) %>%
    dplyr::distinct(Hugo_Symbol, Patient) %>%
    dplyr::group_by(Hugo_Symbol) %>%
    dplyr::summarise(MutationFreq = dplyr::n_distinct(Patient) / length(pats), .groups = "drop")
  gen_ev$MutationFreq <- mf$MutationFreq[match(gen_ev$Gene, mf$Hugo_Symbol)]
  gen_ev$MutationFreq[is.na(gen_ev$MutationFreq)] <- 0
}
if (!is.null(cnv_mat)) {
  cg <- intersect(ddr_genes, rownames(cnv_mat))
  cs <- intersect(colnames(cnv_mat), tumor_samples)
  if (length(cs) == 0) cs <- colnames(cnv_mat)
  if (length(cg) > 0 && length(cs) > 5) {
    gen_ev$AmpFreq[match(cg, gen_ev$Gene)] <- rowMeans(cnv_mat[cg, cs, drop = FALSE] > 0.2, na.rm = TRUE)
    gen_ev$DelFreq[match(cg, gen_ev$Gene)] <- rowMeans(cnv_mat[cg, cs, drop = FALSE] < -0.2, na.rm = TRUE)
    gen_ev$CNVAlterFreq[match(cg, gen_ev$Gene)] <- rowMeans(abs(cnv_mat[cg, cs, drop = FALSE]) > 0.2, na.rm = TRUE)
  }
}
gen_ev$GenomicAlteration <- minmax01(gen_ev$MutationFreq) + minmax01(gen_ev$CNVAlterFreq)
gen_ev$GenomicAlteration <- minmax01(gen_ev$GenomicAlteration)

# E3: PPI network context
read_ppi <- function(file) {
  if (!file.exists(file)) return(NULL)
  p <- tryCatch(safe_fread(file), error = function(e) NULL)
  if (is.null(p) || nrow(p) == 0) return(NULL)
  g1 <- find_col(p, c("protein1", "gene1", "Gene1", "preferredName_A", "node1", "from"))
  g2 <- find_col(p, c("protein2", "gene2", "Gene2", "preferredName_B", "node2", "to"))
  sc <- find_col(p, c("combined_score", "score", "confidence", "weight"))
  if (is.na(g1) || is.na(g2)) return(NULL)
  out <- data.frame(GeneA = clean_gene(p[[g1]]), GeneB = clean_gene(p[[g2]]), PPI_score = if (!is.na(sc)) as.numeric(p[[sc]]) else 1)
  out <- out[!is.na(out$GeneA) & !is.na(out$GeneB) & out$GeneA != out$GeneB, ]
  out$PPI_score[!is.finite(out$PPI_score)] <- 1
  out$PPI_score <- minmax01(out$PPI_score)
  out
}
ppi <- read_ppi(PPI_FILE)
net_ev <- data.frame(Gene = ddr_genes, PPI_Degree = 0, PPI_WeightedDegree = 0, NetworkContext = 0)
if (!is.null(ppi)) {
  ppi_ddr <- ppi %>% dplyr::filter(GeneA %in% ddr_genes, GeneB %in% ddr_genes)
  if (nrow(ppi_ddr) > 0) {
    dl <- dplyr::bind_rows(ppi_ddr %>% dplyr::select(Gene = GeneA, PPI_score), ppi_ddr %>% dplyr::select(Gene = GeneB, PPI_score))
    ds <- dl %>% dplyr::group_by(Gene) %>% dplyr::summarise(PPI_Degree = dplyr::n(), PPI_WeightedDegree = sum(PPI_score, na.rm = TRUE), .groups = "drop")
    net_ev$PPI_Degree <- ds$PPI_Degree[match(net_ev$Gene, ds$Gene)]
    net_ev$PPI_WeightedDegree <- ds$PPI_WeightedDegree[match(net_ev$Gene, ds$Gene)]
    net_ev$PPI_Degree[is.na(net_ev$PPI_Degree)] <- 0
    net_ev$PPI_WeightedDegree[is.na(net_ev$PPI_WeightedDegree)] <- 0
    net_ev$NetworkContext <- minmax01(log1p(net_ev$PPI_Degree)) + minmax01(net_ev$PPI_WeightedDegree)
    net_ev$NetworkContext <- minmax01(net_ev$NetworkContext)
  }
}

# E4: co-expression coherence
tumor_ddr <- log2(tcga_expr[ddr_genes, tumor_samples, drop = FALSE] + 1)
tumor_ddr_z <- t(z_rows(tumor_ddr))
cc <- suppressWarnings(cor(tumor_ddr_z, method = "spearman", use = "pairwise.complete.obs"))
diag(cc) <- NA
mean_abs <- rowMeans(abs(cc), na.rm = TRUE)
top_abs <- apply(abs(cc), 1, function(x) {
  x <- sort(x[is.finite(x)], decreasing = TRUE)
  if (length(x) == 0) return(0)
  mean(head(x, max(1, ceiling(length(x) * 0.10))), na.rm = TRUE)
})
co_ev <- data.frame(Gene = ddr_genes, MeanAbsDDRCoexpression = mean_abs, Top10AbsDDRCoexpression = top_abs)
co_ev$CoexpressionCoherence <- minmax01(co_ev$MeanAbsDDRCoexpression) + minmax01(co_ev$Top10AbsDDRCoexpression)
co_ev$CoexpressionCoherence <- minmax01(co_ev$CoexpressionCoherence)

evidence <- expr_ev %>% dplyr::left_join(gen_ev, by = "Gene") %>% dplyr::left_join(net_ev, by = "Gene") %>% dplyr::left_join(co_ev, by = "Gene")
evidence_cols <- c("ExprDysregulation", "GenomicAlteration", "NetworkContext", "CoexpressionCoherence")
for (v in evidence_cols) evidence[[v]][!is.finite(evidence[[v]])] <- 0

ef <- entropy_fusion(evidence[, c("Gene", evidence_cols)], "Gene")
layer_weights <- ef$weights
norm_ev <- as.data.frame(ef$normalized, check.names = FALSE)
norm_ev$Gene <- rownames(ef$normalized)
AD_DDRweight <- evidence %>% dplyr::left_join(norm_ev, by = "Gene", suffix = c("", "_Norm"))

AD_DDRweight$Raw_AD_DDRweight <- 0
for (v in evidence_cols) {
  AD_DDRweight$Raw_AD_DDRweight <- AD_DDRweight$Raw_AD_DDRweight + AD_DDRweight[[paste0(v, "_Norm")]] * layer_weights[v]
}
AD_DDRweight$AD_DDRweight <- minmax01(AD_DDRweight$Raw_AD_DDRweight)
AD_DDRweight <- AD_DDRweight[order(AD_DDRweight$AD_DDRweight, decreasing = TRUE), ]

layer_df <- data.frame(EvidenceLayer = names(layer_weights), AdaptiveWeight = as.numeric(layer_weights), Entropy = as.numeric(ef$entropy[names(layer_weights)]), Divergence = as.numeric(ef$divergence[names(layer_weights)]))
save_csv(evidence, file.path(PROC_DIR, "MO_DDR_gene_level_evidence_table.csv"))
save_csv(layer_df, file.path(PROC_DIR, "MO_DDR_entropy_adaptive_layer_weights.csv"))
save_csv(AD_DDRweight, file.path(PROC_DIR, "AD_DDRweight_gene_table.csv"))
save_csv(AD_DDRweight, file.path(DB_DIR, "gene_weight_table.csv"))
print(layer_df)

############################
# 5. Fig2
############################

# evidence heatmap
top_evi <- head(AD_DDRweight$Gene, min(40, nrow(AD_DDRweight)))
ev_mat <- as.matrix(AD_DDRweight[match(top_evi, AD_DDRweight$Gene), evidence_cols])
ev_mat <- apply(ev_mat, 2, minmax01)
rownames(ev_mat) <- top_evi
if (has_pkg("pheatmap")) {
  pdf(file.path(FIG2_DIR, "Fig2B_gene_level_multiomics_evidence_heatmap.pdf"), 6.5, 8)
  pheatmap::pheatmap(ev_mat, cluster_cols = FALSE, color = colorRampPalette(c("#F7FBFF", "#6BAED6", "#08306B"))(100), fontsize_row = 7, main = "Gene-level multi-omics evidence")
  dev.off()
}
p_layer <- ggplot(layer_df, aes(x = reorder(EvidenceLayer, AdaptiveWeight), y = AdaptiveWeight)) +
  geom_col(fill = "#4A5568") + coord_flip() + theme_bw(base_size = 12) +
  labs(x = NULL, y = "Entropy-derived adaptive weight", title = "Adaptive evidence-layer contribution")
safe_ggsave(file.path(FIG2_DIR, "Fig2C_entropy_adaptive_evidence_weights.pdf"), p_layer, 5.5, 4)

p_top <- ggplot(head(AD_DDRweight, 30), aes(x = reorder(Gene, AD_DDRweight), y = AD_DDRweight)) +
  geom_segment(aes(xend = Gene, y = 0, yend = AD_DDRweight), color = "grey60") +
  geom_point(size = 2.5, color = "#C53030") + coord_flip() + theme_bw(base_size = 12) +
  labs(x = NULL, y = "MO-DDRweight", title = "Top MO-DDRweighted genes")
safe_ggsave(file.path(FIG2_DIR, "Fig2D_top_MO_DDRweight_genes.pdf"), p_top, 6, 6)

p_dist <- ggplot(AD_DDRweight, aes(x = AD_DDRweight)) + geom_histogram(bins = 30, fill = "#4575B4", color = "white") +
  theme_bw(base_size = 12) + labs(x = "MO-DDRweight", y = "Gene count", title = "MO-DDRweight distribution")
safe_ggsave(file.path(FIG2_DIR, "Fig2D_MO_DDRweight_distribution.pdf"), p_dist, 5, 4)

############################
# 6. Sample score and baselines
############################

compute_score <- function(expr, samples, wt, score_name, ref_samples = tumor_samples) {
  genes <- intersect(wt$Gene, rownames(expr))
  expr_log <- log2(expr[genes, samples, drop = FALSE] + 1)
  ref_log <- log2(expr[genes, ref_samples, drop = FALSE] + 1)
  mu <- rowMeans(ref_log, na.rm = TRUE)
  sdv <- apply(ref_log, 1, sd, na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- 1
  z <- sweep(sweep(expr_log, 1, mu, "-"), 1, sdv, "/")
  w <- wt$AD_DDRweight[match(genes, wt$Gene)]
  w[is.na(w)] <- 0
  sc <- colSums(z * w, na.rm = TRUE) / sum(abs(w), na.rm = TRUE)
  out <- data.frame(Sample = names(sc), Patient = patient_id(names(sc)), SampleType = sample_type(names(sc)), SampleClass = sample_class(names(sc)))
  out[[score_name]] <- as.numeric(sc)
  out
}

score_samples <- unique(c(tumor_samples, normal_samples))
ad_score <- compute_score(tcga_expr, score_samples, AD_DDRweight[, c("Gene", "AD_DDRweight")], "AD_DDRscore")

eq <- compute_score(tcga_expr, score_samples, data.frame(Gene = ddr_genes, AD_DDRweight = 1), "Equal_DDRscore")
exprwt <- AD_DDRweight[, c("Gene", "ExprDysregulation_Norm")]; colnames(exprwt) <- c("Gene", "AD_DDRweight")
expronly <- compute_score(tcga_expr, score_samples, exprwt, "ExpressionOnly_DDRscore")
nonet <- data.frame(Gene = AD_DDRweight$Gene, AD_DDRweight = rowMeans(AD_DDRweight[, c("ExprDysregulation_Norm", "GenomicAlteration_Norm")], na.rm = TRUE))
nonets <- compute_score(tcga_expr, score_samples, nonet, "NoNetwork_DDRscore")

ad_score <- ad_score %>% left_join(eq[, c("Sample", "Equal_DDRscore")], by = "Sample") %>% left_join(expronly[, c("Sample", "ExpressionOnly_DDRscore")], by = "Sample") %>% left_join(nonets[, c("Sample", "NoNetwork_DDRscore")], by = "Sample")
ad_score <- make_score_group(
  score_df = ad_score,
  clin_df = tcga_clin,
  score_col = "AD_DDRscore",
  method = SCORE_GROUP_METHOD,
  out_file = file.path(PROC_DIR, "MO_DDRscore_group_cutoff.csv")
)
save_csv(ad_score, file.path(PROC_DIR, paste0(FOCUS_CANCER, "_AD_DDRscore.csv")))
save_csv(ad_score, file.path(DB_DIR, "AD_DDRscore_table.csv"))

p_tn <- ggplot(ad_score[ad_score$SampleClass %in% c("Tumor", "Normal"), ], aes(SampleClass, AD_DDRscore, fill = SampleClass)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) + geom_jitter(width = 0.15, size = 0.6, alpha = 0.35) +
  theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = NULL, y = "MO-DDRscore", title = "MO-DDRscore in tumor and normal")
safe_ggsave(file.path(FIG2_DIR, "Fig2E_MO_DDRscore_tumor_normal_boxplot.pdf"), p_tn, 4.5, 4)

p_hist <- ggplot(ad_score[ad_score$SampleClass == "Tumor", ], aes(AD_DDRscore, fill = AD_DDRscore_group)) +
  geom_histogram(bins = 35, color = "white", alpha = 0.85) + theme_bw(base_size = 12) +
  labs(x = "MO-DDRscore", y = "Tumor sample count", fill = "Group")
safe_ggsave(file.path(FIG2_DIR, "Fig2E_MO_DDRscore_distribution_groups.pdf"), p_hist, 5, 4)

############################
# 7. Multi-omics table and baseline comparison
############################

tumor_ad <- ad_score[ad_score$SampleClass == "Tumor", ] %>% arrange(Patient, SampleType) %>% distinct(Patient, .keep_all = TRUE)
ad_score_clin <- tumor_ad %>% inner_join(tcga_clin, by = "Patient")

tmb_df <- NULL
if (!is.null(mut)) {
  tmb_df <- mut %>% filter(!is.na(Patient), !is.na(Hugo_Symbol)) %>% group_by(Patient) %>%
    summarise(MutationCount = n(), MutatedGeneCount = n_distinct(Hugo_Symbol), DDR_MutatedGeneCount = n_distinct(Hugo_Symbol[Hugo_Symbol %in% DDR_GENES]), .groups = "drop")
}
cnv_df <- NULL
if (!is.null(cnv_mat)) {
  cnv_df <- data.frame(Sample = colnames(cnv_mat), Patient = patient_id(colnames(cnv_mat)),
                       CNVBurden = colMeans(abs(cnv_mat) > 0.2, na.rm = TRUE),
                       AmpBurden = colMeans(cnv_mat > 0.2, na.rm = TRUE),
                       DelBurden = colMeans(cnv_mat < -0.2, na.rm = TRUE)) %>%
    group_by(Patient) %>% summarise(CNVBurden = mean(CNVBurden, na.rm = TRUE), AmpBurden = mean(AmpBurden, na.rm = TRUE), DelBurden = mean(DelBurden, na.rm = TRUE), .groups = "drop")
}
multiomics_score <- ad_score_clin
if (!is.null(tmb_df)) multiomics_score <- multiomics_score %>% left_join(tmb_df, by = "Patient")
if (!is.null(cnv_df)) multiomics_score <- multiomics_score %>% left_join(cnv_df, by = "Patient")
for (v in c("MutationCount", "MutatedGeneCount", "DDR_MutatedGeneCount", "CNVBurden", "AmpBurden", "DelBurden")) if (!v %in% colnames(multiomics_score)) multiomics_score[[v]] <- NA_real_
save_csv(multiomics_score, file.path(PROC_DIR, paste0(FOCUS_CANCER, "_multiomics_score_table.csv")))
save_csv(multiomics_score, file.path(DB_DIR, "multiomics_score_table.csv"))

score_sets <- list(
  DNA_Repair = c("BRCA1","BRCA2","RAD51","RAD50","MRE11","NBN","ATM","ATR","CHEK1","CHEK2","ERCC1","ERCC2","MSH2","MSH6","MLH1","PMS2","XRCC1","XRCC5","XRCC6","PRKDC","PARP1","LIG1","LIG3","LIG4","FANCA","FANCD2","FANCI"),
  Cell_Cycle_G2M = c("AURKA","AURKB","BUB1","BUB1B","CCNB1","CCNB2","CDC20","CDC25A","CDC25B","CDC25C","CDK1","MKI67","PLK1","TOP2A"),
  DNA_Replication = c("MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","PCNA","POLE","POLD1","RFC1","RFC2","RPA1","RPA2","CDC6","CDC45"),
  Checkpoint = c("ATM","ATR","CHEK1","CHEK2","TP53","TP53BP1","CLSPN","TOPBP1","HUS1","RAD9A","RAD17")
)
expr_tumor_log <- log2(tcga_expr[, tumor_samples, drop = FALSE] + 1)
expr_tumor_z <- z_rows(expr_tumor_log)
ps <- list()
for (nm in names(score_sets)) {
  gs <- intersect(clean_gene(score_sets[[nm]]), rownames(expr_tumor_z))
  if (length(gs) >= 3) ps[[nm]] <- colMeans(expr_tumor_z[gs, , drop = FALSE], na.rm = TRUE)
}
path_score <- as.data.frame(ps, check.names = FALSE)
path_score$Sample <- rownames(path_score)
path_score <- path_score %>% left_join(ad_score[, c("Sample", "Patient", "AD_DDRscore", "Equal_DDRscore", "ExpressionOnly_DDRscore", "NoNetwork_DDRscore", "AD_DDRscore_group")], by = "Sample") %>%
  left_join(multiomics_score[, c("Patient", "MutationCount", "CNVBurden", "AmpBurden", "DelBurden")], by = "Patient")
save_csv(path_score, file.path(FIG2_DIR, "Fig2_pathway_genomic_score_table.csv"))

bio_feats <- setdiff(colnames(path_score), c("Sample","Patient","AD_DDRscore","Equal_DDRscore","ExpressionOnly_DDRscore","NoNetwork_DDRscore","AD_DDRscore_group"))
bio_feats <- bio_feats[sapply(bio_feats, function(v) is.numeric(path_score[[v]]) || is.integer(path_score[[v]]))]
score_cols <- c("AD_DDRscore", "Equal_DDRscore", "ExpressionOnly_DDRscore", "NoNetwork_DDRscore")
base_stats <- bind_rows(lapply(score_cols, function(sc) {
  bind_rows(lapply(bio_feats, function(ft) {
    cr <- safe_cor(path_score[[sc]], path_score[[ft]])
    data.frame(ScoreType = sc, Feature = ft, Correlation = cr["cor"], P = cr["p"], AbsCorrelation = abs(cr["cor"]))
  }))
}))
base_stats$FDR <- p.adjust(base_stats$P, "BH")
save_csv(base_stats, file.path(FIG2_DIR, "Fig2_MO_DDRscore_baseline_comparison.csv"))

p_base <- ggplot(base_stats, aes(x = reorder(Feature, AbsCorrelation), y = Correlation, fill = ScoreType)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_col(position = "dodge", width = 0.75) +
  coord_flip() + theme_bw(base_size = 11) + labs(x = NULL, y = "Spearman correlation", fill = "Score", title = "MO-DDRscore vs baseline scores")
safe_ggsave(file.path(FIG2_DIR, "Fig2E_MO_DDRscore_vs_baselines.pdf"), p_base, 8, 5.8)

cindex_df <- bind_rows(lapply(score_cols, function(sc) data.frame(ScoreType = sc, Cindex = calc_cindex(ad_score_clin$time, ad_score_clin$status, ad_score_clin[[sc]]))))
save_csv(cindex_df, file.path(FIG2_DIR, "Fig2E_score_level_Cindex_comparison.csv"))

############################
# 7.5 Formal ablation analysis for MO-DDRscore
############################

cat("\nRunning formal ablation analysis for MO-DDRscore...\n")

ABLATION_DIR <- file.path(FIG2_DIR, "Fig2G_ablation")
dir.create(ABLATION_DIR, recursive = TRUE, showWarnings = FALSE)

build_ablation_weight <- function(weight_table, layer_cols, score_name) {
  tmp <- weight_table[, c("Gene", layer_cols), drop = FALSE]
  for (cc in layer_cols) tmp[[cc]][!is.finite(tmp[[cc]])] <- 0

  if (length(layer_cols) == 1) {
    w <- minmax01(tmp[[layer_cols]])
    layer_weight_df <- data.frame(ScoreType = score_name, EvidenceLayer = layer_cols, LayerWeight = 1)
  } else {
    ef2 <- entropy_fusion(tmp, id_col = "Gene")
    w <- as.numeric(ef2$normalized %*% ef2$weights[colnames(ef2$normalized)])
    w <- minmax01(w)
    layer_weight_df <- data.frame(
      ScoreType = score_name,
      EvidenceLayer = names(ef2$weights),
      LayerWeight = as.numeric(ef2$weights)
    )
  }

  list(
    weight = data.frame(Gene = tmp$Gene, AD_DDRweight = w, ScoreType = score_name),
    layer_weight = layer_weight_df
  )
}

ablation_weight_list <- list()
ablation_layer_weight_list <- list()

ablation_weight_list[["Equal_weight"]] <- data.frame(
  Gene = AD_DDRweight$Gene,
  AD_DDRweight = 1,
  ScoreType = "Equal_weight"
)
ablation_layer_weight_list[["Equal_weight"]] <- data.frame(
  ScoreType = "Equal_weight",
  EvidenceLayer = "Equal",
  LayerWeight = 1
)

abl_defs <- list(
  Expression_only = c("ExprDysregulation_Norm"),
  Expr_Genomic = c("ExprDysregulation_Norm", "GenomicAlteration_Norm"),
  Expr_Genomic_Network = c("ExprDysregulation_Norm", "GenomicAlteration_Norm", "NetworkContext_Norm"),
  Full_MO_DDRscore = c("ExprDysregulation_Norm", "GenomicAlteration_Norm", "NetworkContext_Norm", "CoexpressionCoherence_Norm")
)

for (nm in names(abl_defs)) {
  bw <- build_ablation_weight(AD_DDRweight, abl_defs[[nm]], nm)
  ablation_weight_list[[nm]] <- bw$weight
  ablation_layer_weight_list[[nm]] <- bw$layer_weight
}

ablation_weight_df <- dplyr::bind_rows(ablation_weight_list)
ablation_layer_weight_df <- dplyr::bind_rows(ablation_layer_weight_list)

save_csv(ablation_weight_df, file.path(ABLATION_DIR, "Fig2G_ablation_gene_weights.csv"))
save_csv(ablation_layer_weight_df, file.path(ABLATION_DIR, "Fig2G_ablation_layer_weights.csv"))

ablation_score_list <- list()
for (nm in names(ablation_weight_list)) {
  score_col <- paste0(nm, "_score")
  tmp_score <- compute_score(
    expr = tcga_expr,
    samples = score_samples,
    wt = ablation_weight_list[[nm]][, c("Gene", "AD_DDRweight")],
    score_name = score_col
  )
  ablation_score_list[[nm]] <- tmp_score[, c("Sample", "Patient", "SampleType", "SampleClass", score_col)]
}

ablation_scores <- ablation_score_list[[1]]
for (i in 2:length(ablation_score_list)) {
  ablation_scores <- ablation_scores %>%
    dplyr::left_join(
      ablation_score_list[[i]][, c("Sample", paste0(names(ablation_score_list)[i], "_score"))],
      by = "Sample"
    )
}
save_csv(ablation_scores, file.path(ABLATION_DIR, "Fig2G_ablation_sample_scores.csv"))

ablation_score_cols <- grep("_score$", colnames(ablation_scores), value = TRUE)

ablation_eval_df <- path_score[, c("Sample", bio_feats), drop = FALSE] %>%
  dplyr::left_join(ablation_scores[, c("Sample", ablation_score_cols), drop = FALSE], by = "Sample")

ablation_cor_stats <- list()
for (sc in ablation_score_cols) {
  for (ft in bio_feats) {
    cr <- safe_cor(ablation_eval_df[[sc]], ablation_eval_df[[ft]])
    ablation_cor_stats[[paste(sc, ft, sep = "__")]] <- data.frame(
      ScoreType = gsub("_score$", "", sc),
      Feature = ft,
      Correlation = cr["cor"],
      P = cr["p"],
      AbsCorrelation = abs(cr["cor"])
    )
  }
}
ablation_cor_stats <- dplyr::bind_rows(ablation_cor_stats)
ablation_cor_stats$FDR <- p.adjust(ablation_cor_stats$P, "BH")
save_csv(ablation_cor_stats, file.path(ABLATION_DIR, "Fig2G_ablation_biological_correlation_stats.csv"))

ablation_surv_df <- ablation_scores %>%
  dplyr::filter(SampleClass == "Tumor") %>%
  dplyr::arrange(Patient, SampleType) %>%
  dplyr::distinct(Patient, .keep_all = TRUE) %>%
  dplyr::inner_join(tcga_clin, by = "Patient")

ablation_surv_stats <- list()
for (sc in ablation_score_cols) {
  cidx <- calc_cindex(ablation_surv_df$time, ablation_surv_df$status, ablation_surv_df[[sc]])
  cidx_adj <- ifelse(is.finite(cidx), max(cidx, 1 - cidx), NA_real_)
  direction <- ifelse(is.finite(cidx) && cidx < 0.5, "reversed", "original")

  cox_fit <- tryCatch(
    summary(coxph(Surv(time, status) ~ ablation_surv_df[[sc]], data = ablation_surv_df)),
    error = function(e) NULL
  )

  if (!is.null(cox_fit)) {
    ablation_surv_stats[[sc]] <- data.frame(
      ScoreType = gsub("_score$", "", sc),
      Cindex = cidx,
      AdjustedCindex = cidx_adj,
      Direction = direction,
      HR = cox_fit$coefficients[1, "exp(coef)"],
      CoxP = cox_fit$coefficients[1, "Pr(>|z|)"],
      CI_low = cox_fit$conf.int[1, "lower .95"],
      CI_high = cox_fit$conf.int[1, "upper .95"]
    )
  } else {
    ablation_surv_stats[[sc]] <- data.frame(
      ScoreType = gsub("_score$", "", sc),
      Cindex = cidx,
      AdjustedCindex = cidx_adj,
      Direction = direction,
      HR = NA_real_,
      CoxP = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_
    )
  }
}
ablation_surv_stats <- dplyr::bind_rows(ablation_surv_stats)
save_csv(ablation_surv_stats, file.path(ABLATION_DIR, "Fig2G_ablation_survival_stats.csv"))

score_order <- c("Equal_weight", "Expression_only", "Expr_Genomic", "Expr_Genomic_Network", "Full_MO_DDRscore")

ablation_summary <- ablation_cor_stats %>%
  dplyr::group_by(ScoreType) %>%
  dplyr::summarise(
    MeanAbsCorrelation = mean(AbsCorrelation, na.rm = TRUE),
    MedianAbsCorrelation = median(AbsCorrelation, na.rm = TRUE),
    SignificantFeatureCount = sum(FDR < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    ablation_surv_stats[, c("ScoreType", "Cindex", "AdjustedCindex", "HR", "CoxP")],
    by = "ScoreType"
  )
ablation_summary$ScoreType <- factor(ablation_summary$ScoreType, levels = score_order)
ablation_summary <- ablation_summary[order(ablation_summary$ScoreType), ]
save_csv(ablation_summary, file.path(ABLATION_DIR, "Fig2G_ablation_summary.csv"))

ablation_cor_mat <- ablation_cor_stats %>%
  dplyr::select(ScoreType, Feature, Correlation) %>%
  tidyr::pivot_wider(names_from = ScoreType, values_from = Correlation)
ablation_cor_mat_df <- as.data.frame(ablation_cor_mat)
rownames(ablation_cor_mat_df) <- ablation_cor_mat_df$Feature
ablation_cor_mat_df$Feature <- NULL
score_order2 <- intersect(score_order, colnames(ablation_cor_mat_df))
ablation_cor_mat_df <- ablation_cor_mat_df[, score_order2, drop = FALSE]

if (has_pkg("pheatmap")) {
  pdf(file.path(ABLATION_DIR, "Fig2G_ablation_biological_correlation_heatmap.pdf"), width = 7.5, height = 6.5)
  pheatmap::pheatmap(
    as.matrix(ablation_cor_mat_df),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    number_format = "%.2f",
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    main = "Ablation analysis: biological consistency"
  )
  dev.off()
}

ablation_summary_long <- ablation_summary %>%
  dplyr::select(ScoreType, MeanAbsCorrelation, AdjustedCindex) %>%
  tidyr::pivot_longer(cols = c("MeanAbsCorrelation", "AdjustedCindex"), names_to = "Metric", values_to = "Value")

p_ablation_summary <- ggplot(ablation_summary_long, aes(x = ScoreType, y = Value, fill = Metric)) +
  geom_col(position = "dodge", width = 0.72) +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = "Metric value", fill = NULL, title = "Ablation analysis of MO-DDRscore")
safe_ggsave(file.path(ABLATION_DIR, "Fig2G_ablation_summary_barplot.pdf"), p_ablation_summary, 7, 4.8)

p_ablation_cindex <- ggplot(ablation_surv_stats, aes(x = factor(ScoreType, levels = score_order), y = AdjustedCindex)) +
  geom_col(fill = "#2B6CB0", width = 0.72) +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = "Adjusted C-index", title = "Ablation analysis: survival association")
safe_ggsave(file.path(ABLATION_DIR, "Fig2G_ablation_survival_Cindex.pdf"), p_ablation_cindex, 6, 4.5)

cat("Formal ablation analysis finished.\n")

############################
# 8. SPIDR external validation
############################

read_spidr <- function(file) {
  if (!file.exists(file)) return(NULL)
  sp <- tryCatch(safe_fread(file), error = function(e) NULL)
  if (is.null(sp)) return(NULL)
  g1 <- if ("gene_a_base" %in% colnames(sp)) "gene_a_base" else if ("gene_a_raw" %in% colnames(sp)) "gene_a_raw" else find_col(sp, c("^GeneA$", "^gene_a$", "^Gene1$", "^gene1$"))
  g2 <- if ("gene_b_base" %in% colnames(sp)) "gene_b_base" else if ("gene_b_raw" %in% colnames(sp)) "gene_b_raw" else find_col(sp, c("^GeneB$", "^gene_b$", "^Gene2$", "^gene2$"))
  if (is.na(g1) || is.na(g2) || !"gemini_sensitive" %in% colnames(sp)) return(NULL)
  lab <- if ("sl_label_le_1_0" %in% colnames(sp)) "sl_label_le_1_0" else if ("strong_sl_label_le_1_5" %in% colnames(sp)) "strong_sl_label_le_1_5" else NA
  out <- data.frame(GeneA = clean_gene(sp[[g1]]), GeneB = clean_gene(sp[[g2]]), GeminiScore = as.numeric(sp$gemini_sensitive))
  out$SL_Label <- if (!is.na(lab)) as.integer(as.character(sp[[lab]])) else ifelse(out$GeminiScore <= -1, 1, 0)
  out <- out[!is.na(out$GeneA) & !is.na(out$GeneB) & out$GeneA != out$GeneB & is.finite(out$GeminiScore), ]
  out$PairKey <- ifelse(out$GeneA < out$GeneB, paste(out$GeneA, out$GeneB, sep = "__"), paste(out$GeneB, out$GeneA, sep = "__"))
  out <- out %>% group_by(PairKey) %>% summarise(GeneA = first(GeneA), GeneB = first(GeneB), GeminiScore = min(GeminiScore, na.rm = TRUE), SL_Label = max(SL_Label, na.rm = TRUE), .groups = "drop")
  out$SLStrength <- minmax01(pmax(0, -out$GeminiScore))
  out
}
spidr <- read_spidr(SPIDR_FILE)
if (!is.null(spidr)) {
  spidr_ddr <- spidr %>% filter(GeneA %in% ddr_genes, GeneB %in% ddr_genes)
  if (nrow(spidr_ddr) > 0) {
    sp_gene <- bind_rows(spidr_ddr %>% select(Gene = GeneA, SLStrength, SL_Label), spidr_ddr %>% select(Gene = GeneB, SLStrength, SL_Label)) %>%
      group_by(Gene) %>% summarise(SPIDR_PairCount = n(), SPIDR_PositivePairCount = sum(SL_Label == 1, na.rm = TRUE), SPIDR_PositivePairFraction = mean(SL_Label == 1, na.rm = TRUE), SPIDR_MeanSLStrength = mean(SLStrength, na.rm = TRUE), SPIDR_Top10SLStrength = mean(head(sort(SLStrength, decreasing = TRUE), max(1, ceiling(length(SLStrength) * 0.10))), na.rm = TRUE), .groups = "drop")
    spv <- AD_DDRweight %>% select(Gene, AD_DDRweight) %>% left_join(sp_gene, by = "Gene")
    for (v in setdiff(colnames(spv), c("Gene", "AD_DDRweight"))) spv[[v]][is.na(spv[[v]])] <- 0
    save_csv(spv, file.path(FIG2_DIR, "Fig2F_SPIDR_external_dependency_validation.csv"))
    sp_cor <- bind_rows(lapply(c("SPIDR_PositivePairCount","SPIDR_PositivePairFraction","SPIDR_MeanSLStrength","SPIDR_Top10SLStrength"), function(v) {
      cr <- safe_cor(spv$AD_DDRweight, spv[[v]])
      data.frame(Feature = v, Cor = cr["cor"], P = cr["p"])
    }))
    sp_cor$FDR <- p.adjust(sp_cor$P, "BH")
    save_csv(sp_cor, file.path(FIG2_DIR, "Fig2F_SPIDR_dependency_correlation_stats.csv"))
    p_sp <- ggplot(sp_cor, aes(x = reorder(Feature, Cor), y = Cor)) + geom_hline(yintercept = 0, linetype = 2) + geom_col(fill = "#805AD5") + coord_flip() + theme_bw(base_size = 12) + labs(x = NULL, y = "Spearman correlation with MO-DDRweight", title = "SPIDR external dependency validation")
    safe_ggsave(file.path(FIG2_DIR, "Fig2F_SPIDR_dependency_correlation.pdf"), p_sp, 6, 4.5)
  }
}

############################
# 9. DEG
############################

cat("\nRunning DEG analysis...\n")
de_samples <- intersect(ad_score$Sample[ad_score$SampleClass == "Tumor" & !is.na(ad_score$AD_DDRscore_group)], colnames(tcga_expr))
group <- factor(ad_score$AD_DDRscore_group[match(de_samples, ad_score$Sample)], levels = c("Low", "High"))
expr_de <- log2(tcga_expr[, de_samples, drop = FALSE] + 1)

if (has_pkg("limma")) {
  design <- model.matrix(~0 + group); colnames(design) <- levels(group)
  fit <- limma::lmFit(expr_de, design)
  cont <- limma::makeContrasts(High - Low, levels = design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont))
  deg <- limma::topTable(fit2, number = Inf, adjust.method = "BH")
  deg$Gene <- rownames(deg)
  deg <- deg[, c("Gene", setdiff(colnames(deg), "Gene"))]
  deg$FDR <- deg$adj.P.Val
} else {
  deg <- bind_rows(lapply(rownames(expr_de), function(g) {
    x <- expr_de[g, group == "High"]; y <- expr_de[g, group == "Low"]
    data.frame(Gene = g, logFC = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE), P.Value = tryCatch(t.test(x, y)$p.value, error = function(e) NA_real_))
  }))
  deg$FDR <- p.adjust(deg$P.Value, "BH")
  deg$adj.P.Val <- deg$FDR
}
deg$Significance <- "NS"
deg$Significance[deg$FDR < DEG_FDR & deg$logFC > DEG_LOGFC] <- "Up"
deg$Significance[deg$FDR < DEG_FDR & deg$logFC < -DEG_LOGFC] <- "Down"
deg <- deg[order(deg$FDR, -abs(deg$logFC)), ]
save_csv(deg, file.path(FIG3_DIR, "Fig3A_MO_DDRscore_high_vs_low_DEGs.csv"))
save_csv(deg, file.path(DB_DIR, "DEG_table.csv"))

p_vol <- ggplot(deg, aes(logFC, -log10(FDR + 1e-300), color = Significance)) +
  geom_point(size = 0.8, alpha = 0.8) +
  scale_color_manual(values = c(Down = "#2B6CB0", NS = "grey70", Up = "#C53030")) +
  geom_vline(xintercept = c(-DEG_LOGFC, DEG_LOGFC), linetype = 2, color = "grey50") +
  geom_hline(yintercept = -log10(DEG_FDR), linetype = 2, color = "grey50") +
  theme_bw(base_size = 12) + labs(x = "log2 fold change (High vs Low)", y = "-log10(FDR)", title = "DEGs by MO-DDRscore group")
safe_ggsave(file.path(FIG3_DIR, "Fig3A_DEG_volcano.pdf"), p_vol, 6, 5)

topg <- unique(c(head(deg$Gene[deg$Significance == "Up"], MAX_HEATMAP_GENES/2), head(deg$Gene[deg$Significance == "Down"], MAX_HEATMAP_GENES/2)))
topg <- intersect(topg, rownames(expr_de))
if (length(topg) >= 5 && has_pkg("pheatmap")) {
  hm <- z_rows(expr_de[topg, , drop = FALSE])
  hm[hm > 2.5] <- 2.5
  hm[hm < -2.5] <- -2.5

  sample_order_df <- data.frame(
    Sample = de_samples,
    MO_DDRscore_group = group,
    MO_DDRscore = ad_score$AD_DDRscore[match(de_samples, ad_score$Sample)],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(MO_DDRscore_group, MO_DDRscore)

  hm <- hm[, sample_order_df$Sample, drop = FALSE]

  anno <- data.frame(
    MO_DDRscore_group = sample_order_df$MO_DDRscore_group,
    MO_DDRscore = sample_order_df$MO_DDRscore
  )
  rownames(anno) <- sample_order_df$Sample

  gap_pos <- sum(sample_order_df$MO_DDRscore_group == "Low", na.rm = TRUE)
  if (gap_pos <= 0 || gap_pos >= nrow(sample_order_df)) gap_pos <- NULL

  save_csv(data.frame(Gene = rownames(hm), hm, check.names = FALSE), file.path(FIG3_DIR, "Fig3A_top_DEG_heatmap_matrix.csv"))

  pdf(file.path(FIG3_DIR, "Fig3A_top_DEG_heatmap_ordered_by_group.pdf"), 8, 7)
  pheatmap::pheatmap(
    hm,
    annotation_col = anno,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    gaps_col = gap_pos,
    show_colnames = FALSE,
    fontsize_row = 7,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    main = "Top DEGs ordered by MO-DDRscore group"
  )
  dev.off()
}

############################
# 10. Pathway scores
############################

cat("\nRunning pathway analysis...\n")
get_sets <- function() {
  if (has_pkg("msigdbr")) {
    ms <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", category = "H"), error = function(e) NULL)
    if (!is.null(ms) && nrow(ms) > 0) return(lapply(split(ms$gene_symbol, ms$gs_name), clean_gene))
  }
  list(
    HALLMARK_DNA_REPAIR = score_sets$DNA_Repair,
    HALLMARK_G2M_CHECKPOINT = score_sets$Cell_Cycle_G2M,
    HALLMARK_E2F_TARGETS = c("E2F1","E2F2","E2F3","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","PCNA","TYMS","TK1","CDC6","CDK2","CCNE1"),
    HALLMARK_MYC_TARGETS_V1 = c("MYC","NPM1","NCL","RPLP0","RPS3","LDHA","ODC1","CAD","HSPD1","HSPE1","MCM4","MCM5","MCM6"),
    HALLMARK_P53_PATHWAY = c("TP53","CDKN1A","MDM2","GADD45A","BAX","BBC3","PMAIP1","DDB2","RRM2B","SESN1","SESN2"),
    HALLMARK_HYPOXIA = c("HIF1A","VEGFA","CA9","LDHA","SLC2A1","ENO1","PGK1","BNIP3","NDRG1","EGLN3"),
    HALLMARK_INTERFERON_GAMMA_RESPONSE = c("STAT1","IRF1","CXCL9","CXCL10","CXCL11","GBP1","GBP5","IDO1","HLA-DRA","HLA-DRB1"),
    HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = c("VIM","CDH2","SNAI1","SNAI2","TWIST1","ZEB1","ZEB2","COL1A1","COL1A2","FN1","MMP2","MMP9")
  )
}
sets <- get_sets()
expr_z <- z_rows(expr_tumor_log)
hs <- list()
for (nm in names(sets)) {
  gs <- intersect(clean_gene(sets[[nm]]), rownames(expr_z))
  if (length(gs) >= 3) hs[[nm]] <- colMeans(expr_z[gs, , drop = FALSE], na.rm = TRUE)
}
hall_scores <- as.data.frame(hs, check.names = FALSE)
if (ncol(hall_scores) > 0) {
  hall_scores$Sample <- rownames(hall_scores)
  hall_scores <- hall_scores %>% left_join(ad_score[, c("Sample", "AD_DDRscore", "AD_DDRscore_group")], by = "Sample")
  save_csv(hall_scores, file.path(FIG3_DIR, "Fig3B_Hallmark_pathway_scores.csv"))
  save_csv(hall_scores, file.path(DB_DIR, "hallmark_pathway_score_table.csv"))
  hcols <- setdiff(colnames(hall_scores), c("Sample","AD_DDRscore","AD_DDRscore_group"))
  hstat <- bind_rows(lapply(hcols, function(pw) {
    dd <- hall_scores[is.finite(hall_scores[[pw]]) & !is.na(hall_scores$AD_DDRscore_group), ]
    cr <- safe_cor(dd$AD_DDRscore, dd[[pw]])
    data.frame(Pathway = pw, Cor = cr["cor"], CorP = cr["p"], Diff = median(dd[[pw]][dd$AD_DDRscore_group == "High"], na.rm = TRUE) - median(dd[[pw]][dd$AD_DDRscore_group == "Low"], na.rm = TRUE), WilcoxP = tryCatch(wilcox.test(dd[[pw]] ~ dd$AD_DDRscore_group)$p.value, error = function(e) NA_real_))
  }))
  hstat$FDR <- p.adjust(hstat$WilcoxP, "BH")
  hstat <- hstat[order(hstat$FDR, -abs(hstat$Diff)), ]
  save_csv(hstat, file.path(FIG3_DIR, "Fig3B_Hallmark_pathway_difference.csv"))
  p_h <- ggplot(head(hstat, 25), aes(x = reorder(Pathway, Diff), y = Diff)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(size = -log10(WilcoxP + 1e-300), color = Cor), alpha = 0.9) +
    coord_flip() + theme_bw(base_size = 11) + labs(x = NULL, y = "Median difference (High - Low)", size = "-log10(P)", color = "Spearman r", title = "Hallmark pathway alteration")
  safe_ggsave(file.path(FIG3_DIR, "Fig3B_Hallmark_pathway_bubble.pdf"), p_h, 7.2, 6)
  if (has_pkg("pheatmap")) {
    top_pw <- head(hstat$Pathway, min(20, nrow(hstat)))
    # 标准画法：按 MO-DDRscore 高低组分开，并在组内按 score 排序，不再让列聚类打乱分组
    hall_order_df <- hall_scores %>%
      dplyr::filter(!is.na(AD_DDRscore_group), is.finite(AD_DDRscore)) %>%
      dplyr::arrange(AD_DDRscore_group, AD_DDRscore)

    hm <- as.matrix(t(hall_order_df[, top_pw, drop = FALSE]))
    colnames(hm) <- hall_order_df$Sample
    hm <- z_rows(hm)
    hm[hm > 2.5] <- 2.5
    hm[hm < -2.5] <- -2.5

    anno <- data.frame(
      MO_DDRscore_group = hall_order_df$AD_DDRscore_group,
      MO_DDRscore = hall_order_df$AD_DDRscore
    )
    rownames(anno) <- hall_order_df$Sample

    gap_pos <- sum(hall_order_df$AD_DDRscore_group == "Low", na.rm = TRUE)
    if (gap_pos <= 0 || gap_pos >= nrow(hall_order_df)) gap_pos <- NULL

    pdf(file.path(FIG3_DIR, "Fig3B_Hallmark_pathway_heatmap.pdf"), 9, 6.5)
    pheatmap::pheatmap(
      hm,
      annotation_col = anno,
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      gaps_col = gap_pos,
      show_colnames = FALSE,
      color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      fontsize_row = 7,
      main = "Hallmark pathway scores ordered by MO-DDRscore group"
    )
    dev.off()
  }
}

############################
# 11. Mutation / CNV / immune
############################

cat("\nRunning mutation/CNV/immune analyses...\n")
driver <- c("TP53","KRAS","EGFR","STK11","KEAP1","BRAF","MET","ERBB2","ALK","ROS1","RET","NF1","PIK3CA","RBM10","SMARCA4")
if (!is.null(mut)) {
  pg <- tumor_ad %>% select(Patient, AD_DDRscore_group) %>% distinct()
  mut2 <- mut %>% left_join(pg, by = "Patient") %>% filter(!is.na(AD_DDRscore_group))
  ngrp <- table(pg$AD_DDRscore_group)
  dstat <- bind_rows(lapply(driver, function(g) {
    pats <- unique(mut2$Patient[mut2$Hugo_Symbol == g])
    ht <- as.numeric(ngrp["High"]); lt <- as.numeric(ngrp["Low"])
    hm <- sum(pg$Patient[pg$AD_DDRscore_group == "High"] %in% pats)
    lm <- sum(pg$Patient[pg$AD_DDRscore_group == "Low"] %in% pats)
    p <- tryCatch(fisher.test(matrix(c(hm, ht-hm, lm, lt-lm), nrow = 2))$p.value, error = function(e) NA_real_)
    data.frame(Gene = g, High_Freq = hm/ht, Low_Freq = lm/lt, Diff = hm/ht - lm/lt, P = p)
  }))
  dstat$FDR <- p.adjust(dstat$P, "BH")
  save_csv(dstat, file.path(FIG3_DIR, "Fig3C_LUAD_driver_mutation_frequency_stats.csv"))
  dlong <- dstat %>% select(Gene, High_Freq, Low_Freq) %>% pivot_longer(cols = c(High_Freq, Low_Freq), names_to = "Group", values_to = "Frequency")
  dlong$Group <- ifelse(dlong$Group == "High_Freq", "High", "Low")
  p_drv <- ggplot(dlong, aes(x = reorder(Gene, Frequency), y = Frequency, fill = Group)) + geom_col(position = "dodge") + coord_flip() + theme_bw(base_size = 12) + labs(x = NULL, y = "Mutation frequency", fill = "Group", title = "LUAD driver mutation frequency")
  safe_ggsave(file.path(FIG3_DIR, "Fig3C_driver_mutation_frequency_barplot.pdf"), p_drv, 6, 5)
  if (has_pkg("maftools")) {
    clin_anno <- pg %>% select(Tumor_Sample_Barcode = Patient, AD_DDRscore_group)
    maf_obj <- tryCatch(maftools::read.maf(maf = mut2, clinicalData = clin_anno), error = function(e) NULL)
    if (!is.null(maf_obj)) {
      pdf(file.path(FIG3_DIR, "Fig3C_mutation_waterfall_MO_DDRscore_group.pdf"), 9, 7)
      try(maftools::oncoplot(maf_obj, genes = intersect(driver, unique(mut2$Hugo_Symbol)), clinicalFeatures = "AD_DDRscore_group", sortByAnnotation = TRUE), silent = TRUE)
      dev.off()
    }
  }
}

if (!is.null(cnv_df)) {
  cnv_plot <- cnv_df %>% left_join(tumor_ad[, c("Patient", "AD_DDRscore_group")], by = "Patient") %>% filter(!is.na(AD_DDRscore_group))
  cnv_long <- cnv_plot %>% pivot_longer(cols = any_of(c("CNVBurden","AmpBurden","DelBurden")), names_to = "Feature", values_to = "Value")
  cnv_stat <- cnv_long %>% group_by(Feature) %>% summarise(Diff = median(Value[AD_DDRscore_group == "High"], na.rm = TRUE) - median(Value[AD_DDRscore_group == "Low"], na.rm = TRUE), P = tryCatch(wilcox.test(Value ~ AD_DDRscore_group)$p.value, error = function(e) NA_real_), .groups = "drop")
  cnv_stat$FDR <- p.adjust(cnv_stat$P, "BH")
  save_csv(cnv_stat, file.path(FIG3_DIR, "Fig3D_CNV_burden_group_stats.csv"))
  p_cnv <- ggplot(cnv_long, aes(AD_DDRscore_group, Value, fill = AD_DDRscore_group)) + geom_boxplot(outlier.shape = NA, alpha = 0.8) + geom_jitter(width = 0.15, size = 0.55, alpha = 0.35) + facet_wrap(~Feature, scales = "free_y") + theme_bw(base_size = 12) + theme(legend.position = "none") + labs(x = NULL, y = "Burden", title = "CNV burden by MO-DDRscore group")
  safe_ggsave(file.path(FIG3_DIR, "Fig3D_CNV_burden_boxplot.pdf"), p_cnv, 7, 4)
}

immune_sets <- list(
  Cytolytic_Activity = c("GZMA", "PRF1"),
  HLA_Class_I = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2"),
  HLA_Class_II = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1"),
  T_cell_CD8 = c("CD8A","CD8B","GZMB","NKG7","CXCL9","CXCL10"),
  B_cell = c("MS4A1","CD79A","CD79B","CD19"),
  NK_cell = c("NKG7","GNLY","KLRD1","KLRK1"),
  Macrophage = c("CD68","CD163","MSR1","CSF1R"),
  Treg = c("FOXP3","IL2RA","CTLA4","IKZF2")
)
expr_imm_z <- z_rows(expr_tumor_log)
im <- list()
for (nm in names(immune_sets)) {
  gs <- intersect(clean_gene(immune_sets[[nm]]), rownames(expr_imm_z))
  if (length(gs) >= 2) im[[nm]] <- colMeans(expr_imm_z[gs, , drop = FALSE], na.rm = TRUE)
}
immune_score_df <- as.data.frame(im, check.names = FALSE)
if (ncol(immune_score_df) > 0) {
  immune_score_df$Sample <- rownames(immune_score_df)
  immune_score_df <- immune_score_df %>% left_join(ad_score[, c("Sample","AD_DDRscore","AD_DDRscore_group")], by = "Sample")
  save_csv(immune_score_df, file.path(FIG3_DIR, "Fig3E_immune_HLA_CYT_scores.csv"))
  save_csv(immune_score_df, file.path(DB_DIR, "immune_score_table.csv"))
  imcols <- setdiff(colnames(immune_score_df), c("Sample","AD_DDRscore","AD_DDRscore_group"))
  istat <- bind_rows(lapply(imcols, function(v) {
    dd <- immune_score_df[is.finite(immune_score_df[[v]]) & !is.na(immune_score_df$AD_DDRscore_group), ]
    cr <- safe_cor(dd$AD_DDRscore, dd[[v]])
    data.frame(Feature = v, Cor = cr["cor"], CorP = cr["p"], Diff = median(dd[[v]][dd$AD_DDRscore_group == "High"], na.rm = TRUE) - median(dd[[v]][dd$AD_DDRscore_group == "Low"], na.rm = TRUE), WilcoxP = tryCatch(wilcox.test(dd[[v]] ~ dd$AD_DDRscore_group)$p.value, error = function(e) NA_real_))
  }))
  istat$FDR <- p.adjust(istat$WilcoxP, "BH")
  save_csv(istat, file.path(FIG3_DIR, "Fig3E_immune_HLA_CYT_stats.csv"))
  p_im <- ggplot(istat, aes(x = reorder(Feature, Diff), y = Diff)) + geom_hline(yintercept = 0, linetype = 2) + geom_point(aes(size = -log10(WilcoxP + 1e-300), color = Cor), alpha = 0.9) + coord_flip() + theme_bw(base_size = 12) + labs(x = NULL, y = "Median difference (High - Low)", size = "-log10(P)", color = "Correlation", title = "Immune/HLA/CYT score difference")
  safe_ggsave(file.path(FIG3_DIR, "Fig3E_immune_HLA_CYT_bubble.pdf"), p_im, 6.5, 5)
}

chk <- intersect(c("CD274","PDCD1","CTLA4","LAG3","TIGIT","HAVCR2","PDCD1LG2","ICOS","IDO1","CD80","CD86","TNFRSF9","TNFRSF4","CD40","CD40LG","VSIR","SIGLEC15"), rownames(tcga_expr))
if (length(chk) >= 3) {
  chk_df <- as.data.frame(t(log2(tcga_expr[chk, tumor_samples, drop = FALSE] + 1)), check.names = FALSE)
  chk_df$Sample <- rownames(chk_df)
  chk_df <- chk_df %>% left_join(ad_score[, c("Sample","AD_DDRscore_group")], by = "Sample")
  chk_long <- chk_df %>% pivot_longer(cols = all_of(chk), names_to = "Gene", values_to = "Expression")
  chk_stat <- chk_long %>% group_by(Gene) %>% summarise(Diff = median(Expression[AD_DDRscore_group == "High"], na.rm = TRUE) - median(Expression[AD_DDRscore_group == "Low"], na.rm = TRUE), P = tryCatch(wilcox.test(Expression ~ AD_DDRscore_group)$p.value, error = function(e) NA_real_), .groups = "drop")
  chk_stat$FDR <- p.adjust(chk_stat$P, "BH")
  save_csv(chk_stat, file.path(FIG3_DIR, "Fig3F_checkpoint_expression_stats.csv"))
  showg <- head(chk_stat$Gene[order(chk_stat$FDR)], min(12, nrow(chk_stat)))
  p_chk <- ggplot(chk_long[chk_long$Gene %in% showg, ], aes(AD_DDRscore_group, Expression, fill = AD_DDRscore_group)) + geom_boxplot(outlier.shape = NA, alpha = 0.85) + geom_jitter(width = 0.12, size = 0.45, alpha = 0.25) + facet_wrap(~Gene, scales = "free_y", ncol = 4) + theme_bw(base_size = 10) + theme(legend.position = "none") + labs(x = NULL, y = "log2(TPM+1)", title = "Immune checkpoint expression")
  safe_ggsave(file.path(FIG3_DIR, "Fig3F_checkpoint_expression_boxplot.pdf"), p_chk, 8, 6)
}

############################
# 12. Optional DepMap/GDSC
############################

cat("\nOptional DepMap/GDSC analysis...\n")
score_depmap_expression <- function(file, wt) {
  if (!file.exists(file)) return(NULL)
  dep <- tryCatch(fread(file, data.table = FALSE, check.names = FALSE), error = function(e) NULL)
  if (is.null(dep) || ncol(dep) < 5) return(NULL)
  idc <- find_col(dep, c("^ModelID$", "^DepMap_ID$", "^Model$", "^CellLine$", "^Unnamed: 0$", "^X$"))
  if (is.na(idc)) idc <- colnames(dep)[1]
  ids <- as.character(dep[[idc]])
  edf <- dep[, setdiff(colnames(dep), idc), drop = FALSE]
  genes <- clean_gene(gsub("\\s*\\(.*\\)$", "", colnames(edf)))
  keep <- !is.na(genes) & genes != ""
  edf <- edf[, keep, drop = FALSE]; genes <- genes[keep]
  mat <- sapply(edf, function(x) suppressWarnings(as.numeric(x)))
  if (is.null(dim(mat))) return(NULL)
  mat <- as.matrix(mat); rownames(mat) <- ids; colnames(mat) <- genes
  mat <- mat[, colSums(is.finite(mat)) > 0, drop = FALSE]
  ug <- unique(colnames(mat))
  mg <- sapply(ug, function(g) rowMeans(mat[, colnames(mat) == g, drop = FALSE], na.rm = TRUE))
  mg <- as.matrix(mg); rownames(mg) <- rownames(mat)
  common <- intersect(wt$Gene, colnames(mg))
  message("DepMap common weighted genes: ", length(common), "/", nrow(wt))
  if (length(common) < 5) return(NULL)
  mz <- z_cols(mg[, common, drop = FALSE]); mz[!is.finite(mz)] <- 0
  w <- wt$AD_DDRweight[match(colnames(mz), wt$Gene)]; w[is.na(w)] <- 0
  if (sum(abs(w)) == 0) return(NULL)
  data.frame(ModelID = rownames(mz), DepMap_MO_DDRscore = as.numeric(mz %*% w / sum(abs(w))))
}
dep_score <- score_depmap_expression(DEPMAP_EXPR_FILE, AD_DDRweight)
if (!is.null(dep_score)) {
  save_csv(dep_score, file.path(FIG3_DIR, "Fig3F_DepMap_MO_DDRscore.csv"))
  if (file.exists(GDSC_FILE) && has_pkg("readxl")) {
    gdsc <- tryCatch(as.data.frame(readxl::read_excel(GDSC_FILE), check.names = FALSE), error = function(e) NULL)
    if (!is.null(gdsc)) {
      mc <- find_col(gdsc, c("ModelID", "SangerModelID", "COSMIC_ID", "CELL_LINE_NAME", "Cell line"))
      dc <- find_col(gdsc, c("^DRUG_NAME$", "Drug", "drug_name"))
      ic <- find_col(gdsc, c("LN_IC50", "IC50", "AUC"))
      if (!is.na(mc) && !is.na(dc) && !is.na(ic)) {
        gu <- data.frame(ModelID = as.character(gdsc[[mc]]), Drug = as.character(gdsc[[dc]]), Response = as.numeric(gdsc[[ic]]))
        if (file.exists(DEPMAP_MODEL_FILE)) {
          mm <- tryCatch(fread(DEPMAP_MODEL_FILE, data.table = FALSE, check.names = FALSE), error = function(e) NULL)
          if (!is.null(mm)) {
            sc <- find_col(mm, c("SangerModelID", "Sanger"))
            dpc <- find_col(mm, c("^ModelID$", "DepMap_ID", "ModelID"))
            if (!is.na(sc) && !is.na(dpc)) {
              map <- data.frame(GDSC_ModelID = as.character(mm[[sc]]), ModelID2 = as.character(mm[[dpc]]))
              gu <- gu %>% left_join(map, by = c("ModelID" = "GDSC_ModelID")) %>% mutate(ModelID = ifelse(!is.na(ModelID2), ModelID2, ModelID)) %>% select(-ModelID2)
            }
          }
        }
        dd <- gu %>% group_by(ModelID, Drug) %>% summarise(Response = mean(Response, na.rm = TRUE), .groups = "drop") %>% inner_join(dep_score, by = "ModelID")
        dcor <- dd %>% group_by(Drug) %>% summarise(N = n(), Cor = ifelse(N >= 10, safe_cor(DepMap_MO_DDRscore, Response)["cor"], NA_real_), P = ifelse(N >= 10, safe_cor(DepMap_MO_DDRscore, Response)["p"], NA_real_), .groups = "drop") %>% filter(is.finite(Cor), is.finite(P))
        dcor$FDR <- p.adjust(dcor$P, "BH")
        dcor <- dcor[order(dcor$P), ]
        save_csv(dcor, file.path(FIG3_DIR, "Fig3F_GDSC_drug_correlation_all.csv"))
        pat <- paste(c("cisplatin","carboplatin","oxaliplatin","paclitaxel","docetaxel","gemcitabine","etoposide","olaparib","talazoparib","niraparib","veliparib","rucaparib","azd6738","berzosertib","atr","wee1","mk-1775","adavosertib","topotecan","irinotecan","doxorubicin"), collapse = "|")
        focus <- dcor[grepl(pat, dcor$Drug, ignore.case = TRUE), ]
        if (nrow(focus) >= 3) {
          save_csv(focus, file.path(FIG3_DIR, "Fig3F_GDSC_DDR_chemotherapy_focused_drugs.csv"))
          p_drug <- ggplot(head(focus, 30), aes(x = reorder(Drug, Cor), y = Cor)) + geom_hline(yintercept = 0, linetype = 2) + geom_point(aes(size = -log10(P + 1e-300), color = Cor), alpha = 0.9) + coord_flip() + theme_bw(base_size = 12) + labs(x = NULL, y = "Spearman correlation with drug response", size = "-log10(P)", color = "Correlation", title = "Focused DDR/chemotherapy drug association")
          safe_ggsave(file.path(FIG3_DIR, "Fig3F_GDSC_DDR_chemotherapy_focused_bubble.pdf"), p_drug, 7, 6)
        }
      }
    }
  }
}

############################
# 13. Integrated panel and survival
############################

base_int <- multiomics_score
if (exists("hall_scores") && is.data.frame(hall_scores)) {
  keep <- intersect(setdiff(colnames(hall_scores), c("Sample","AD_DDRscore","AD_DDRscore_group")),
                    c("HALLMARK_DNA_REPAIR","HALLMARK_G2M_CHECKPOINT","HALLMARK_E2F_TARGETS","HALLMARK_MYC_TARGETS_V1","HALLMARK_P53_PATHWAY","HALLMARK_HYPOXIA","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"))
  if (length(keep) > 0) base_int <- base_int %>% left_join(hall_scores[, c("Sample", keep), drop = FALSE], by = "Sample")
}
if (exists("immune_score_df") && is.data.frame(immune_score_df)) {
  imc <- setdiff(colnames(immune_score_df), c("Sample","AD_DDRscore","AD_DDRscore_group"))
  base_int <- base_int %>% left_join(immune_score_df[, c("Sample", imc), drop = FALSE], by = "Sample")
}
fint <- setdiff(colnames(base_int), c("Sample","Patient","AD_DDRscore","Equal_DDRscore","ExpressionOnly_DDRscore","NoNetwork_DDRscore","AD_DDRscore_group","SampleType","SampleClass","time","status","age","gender","stage"))
fint <- fint[sapply(fint, function(v) is.numeric(base_int[[v]]) || is.integer(base_int[[v]]))]
corr_panel <- bind_rows(lapply(fint, function(v) {
  cr <- safe_cor(base_int$AD_DDRscore, base_int[[v]])
  data.frame(Feature = v, Correlation = cr["cor"], P = cr["p"])
}))
if (nrow(corr_panel) > 0) {
  corr_panel$FDR <- p.adjust(corr_panel$P, "BH")
  corr_panel <- corr_panel[order(corr_panel$FDR, -abs(corr_panel$Correlation)), ]
  save_csv(corr_panel, file.path(FIG3_DIR, "Fig3_multiomics_MO_DDRscore_correlation_panel.csv"))
  p_corr <- ggplot(head(corr_panel, 30), aes(x = reorder(Feature, Correlation), y = Correlation)) +
    geom_hline(yintercept = 0, linetype = 2) + geom_point(aes(size = -log10(P + 1e-300), color = Correlation), alpha = 0.9) +
    coord_flip() + theme_bw(base_size = 11) + labs(x = NULL, y = "Spearman correlation with MO-DDRscore", size = "-log10(P)", color = "Correlation", title = "Integrated multi-omics correlation panel")
  safe_ggsave(file.path(FIG3_DIR, "Fig3_multiomics_MO_DDRscore_correlation_panel.pdf"), p_corr, 7, 7)
}

run_standard_survival(
  df = ad_score_clin,
  score_col = "AD_DDRscore",
  group_col = "AD_DDRscore_group",
  out_prefix = file.path(FIG3_DIR, "Fig3G_MO_DDRscore"),
  title = paste0(FOCUS_CANCER, " MO-DDRscore survival")
)

############################
# 14. Summary
############################

summary_df <- data.frame(
  Item = c("Cancer","Score_method","DDR_gene_count_used","TCGA_total_expression_samples","TCGA_tumor_samples","TCGA_normal_samples","Clinical_patients","Clinical_events","Weighted_DDR_genes","Tumor_MO_DDRscore_samples","Entropy_weight_ExprDysregulation","Entropy_weight_GenomicAlteration","Entropy_weight_NetworkContext","Entropy_weight_CoexpressionCoherence","DEG_up","DEG_down"),
  Value = c(FOCUS_CANCER,"MO-DDRscore_entropy_adaptive_multiomics_context_logFC1.2_standard_grouping",length(DDR_GENES),ncol(tcga_expr),length(tumor_samples),length(normal_samples),nrow(tcga_clin),sum(tcga_clin$status == 1),nrow(AD_DDRweight),sum(ad_score$SampleClass == "Tumor"),round(layer_weights["ExprDysregulation"],4),round(layer_weights["GenomicAlteration"],4),round(layer_weights["NetworkContext"],4),round(layer_weights["CoexpressionCoherence"],4),sum(deg$Significance == "Up"),sum(deg$Significance == "Down"))
)
save_csv(summary_df, file.path(PROC_DIR, "MO_DDR_clean_pipeline_summary.csv"))

############################
# 15. Results-ready key tables
############################

RESULT_DIR <- file.path(PROJECT_DIR, "06_results_ready_tables")
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

group_counts <- ad_score %>%
  dplyr::filter(SampleClass == "Tumor", !is.na(AD_DDRscore_group)) %>%
  dplyr::count(AD_DDRscore_group, name = "N")
save_csv(group_counts, file.path(RESULT_DIR, "Result_group_counts.csv"))

top_weight_genes <- AD_DDRweight %>%
  dplyr::arrange(dplyr::desc(AD_DDRweight)) %>%
  dplyr::select(Gene, AD_DDRweight, dplyr::contains("_Norm")) %>%
  head(20)
save_csv(top_weight_genes, file.path(RESULT_DIR, "Result_top20_MO_DDRweight_genes.csv"))

if (exists("hstat")) save_csv(head(hstat, 20), file.path(RESULT_DIR, "Result_top20_Hallmark_pathways.csv"))
if (exists("dstat")) save_csv(dstat[order(dstat$FDR), ], file.path(RESULT_DIR, "Result_driver_mutation_stats.csv"))
if (exists("cnv_stat")) save_csv(cnv_stat, file.path(RESULT_DIR, "Result_CNV_burden_stats.csv"))
if (exists("istat")) save_csv(istat[order(istat$FDR), ], file.path(RESULT_DIR, "Result_immune_stats.csv"))
if (exists("ablation_summary")) save_csv(ablation_summary, file.path(RESULT_DIR, "Result_ablation_summary.csv"))

result_key_numbers <- c(
  "MO-DDRscore results-ready summary",
  "=================================",
  paste0("Cancer: ", FOCUS_CANCER),
  paste0("DDR genes used: ", length(ddr_genes)),
  paste0("Tumor samples: ", length(tumor_samples)),
  paste0("Normal samples: ", length(normal_samples)),
  paste0("Clinical patients: ", nrow(tcga_clin), "; events: ", sum(tcga_clin$status == 1)),
  paste0("Low group n: ", group_counts$N[match("Low", group_counts$AD_DDRscore_group)]),
  paste0("High group n: ", group_counts$N[match("High", group_counts$AD_DDRscore_group)]),
  paste0("DEGs: Up=", sum(deg$Significance == "Up", na.rm = TRUE),
         ", Down=", sum(deg$Significance == "Down", na.rm = TRUE),
         ", Total=", sum(deg$Significance != "NS", na.rm = TRUE)),
  "",
  "Evidence layer weights:",
  paste0(layer_df$EvidenceLayer, ": ", round(layer_df$AdaptiveWeight, 4)),
  "",
  "Key output folders:",
  paste0("Fig2: ", FIG2_DIR),
  paste0("Fig3: ", FIG3_DIR),
  paste0("Result-ready tables: ", RESULT_DIR),
  paste0("Database tables: ", DB_DIR)
)
writeLines(result_key_numbers, file.path(RESULT_DIR, "Result_key_numbers_for_writing.txt"))

cat("\nMO-DDRscore pipeline finished successfully.\n")
cat("Fig2:", FIG2_DIR, "\n")
cat("Fig3:", FIG3_DIR, "\n")
cat("Processed:", PROC_DIR, "\n")
print(summary_df)
