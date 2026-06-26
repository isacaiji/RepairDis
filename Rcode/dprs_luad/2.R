############################################################
# Figure 2C: DEG heatmap with clinical annotations
# High vs Low MO-DDRscore
############################################################

#rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260513)

############################
# 0. Parameters
############################

PROJECT_DIR <- "D:/R_workspace/评分/AD_DDR_project"
PROC_DIR <- file.path(PROJECT_DIR, "01_processed")
OUT_DIR <- file.path(PROJECT_DIR, "02_Figure2_MultiOmics_Immune")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

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

TOP_N_GENES <- 60
DEG_ADJ_P <- 0.05
DEG_LOGFC <- 1.20

############################
# 1. Packages
############################

pkgs <- c("data.table", "dplyr", "limma", "ComplexHeatmap", "circlize", "grid")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p %in% c("ComplexHeatmap", "circlize")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(p, ask = FALSE, update = FALSE)
    } else {
      install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(limma)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

############################
# 2. Helper functions
############################

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

save_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file)
}

simple_stage <- function(x) {
  x <- toupper(as.character(x))
  dplyr::case_when(
    grepl("STAGE I[^I]|^I$|^IA|^IB", x) ~ "I",
    grepl("STAGE II[^I]|^II$|^IIA|^IIB", x) ~ "II",
    grepl("STAGE III|^III|^IIIA|^IIIB|^IIIC", x) ~ "III",
    grepl("STAGE IV|^IV|^IVA|^IVB", x) ~ "IV",
    TRUE ~ "Unknown"
  )
}

simple_gender <- function(x) {
  x <- toupper(as.character(x))
  dplyr::case_when(
    x %in% c("MALE", "M", "1") ~ "Male",
    x %in% c("FEMALE", "F", "0") ~ "Female",
    TRUE ~ "Unknown"
  )
}

simple_status <- function(x) {
  x <- as.numeric(x)
  dplyr::case_when(
    x == 1 ~ "Dead",
    x == 0 ~ "Alive",
    TRUE ~ "Unknown"
  )
}

############################
# 3. Load expression matrix
############################

TCGA_EXPR_RDS <- file.path(PROC_DIR, paste0(FOCUS_CANCER, "_tcga_expr_tpm_matrix.rds"))

if (file.exists(TCGA_EXPR_RDS)) {
  tcga_expr <- readRDS(TCGA_EXPR_RDS)
} else {
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
}

rownames(tcga_expr) <- clean_gene(rownames(tcga_expr))
colnames(tcga_expr) <- gsub("\\.", "-", colnames(tcga_expr))
tcga_expr <- tcga_expr[!is.na(rownames(tcga_expr)) & rownames(tcga_expr) != "", , drop = FALSE]
tcga_expr <- tcga_expr[!duplicated(rownames(tcga_expr)), , drop = FALSE]

############################
# 4. Load clinical and MO-DDRscore
############################

tcga_clin <- data.table::fread(
  TCGA_CLIN_PROCESSED_FILE,
  data.table = FALSE,
  check.names = FALSE
)

tcga_clin <- tcga_clin %>%
  dplyr::mutate(
    Patient = as.character(Patient),
    time = suppressWarnings(as.numeric(time)),
    status = suppressWarnings(as.numeric(status))
  )

mo_score <- data.table::fread(
  MO_SCORE_FILE,
  data.table = FALSE,
  check.names = FALSE
)

mo_score <- mo_score %>%
  dplyr::mutate(
    Sample = gsub("\\.", "-", as.character(Sample)),
    Patient = as.character(Patient),
    SampleType = sample_type(Sample),
    MO_DDRscore_raw = suppressWarnings(as.numeric(MO_DDRscore_raw)),
    MO_DDRscore_group = as.character(MO_DDRscore_group)
  )

tumor_samples <- colnames(tcga_expr)[sample_class(colnames(tcga_expr)) == "Tumor"]

plot_meta <- mo_score %>%
  dplyr::filter(
    Sample %in% tumor_samples,
    Sample %in% colnames(tcga_expr),
    MO_DDRscore_group %in% c("Low", "High"),
    is.finite(MO_DDRscore_raw)
  ) %>%
  dplyr::arrange(Patient, SampleType) %>%
  dplyr::distinct(Patient, .keep_all = TRUE) %>%
  dplyr::left_join(tcga_clin, by = "Patient")

# 自动识别临床列
pick_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

age_col <- pick_col(tcga_clin, c("age", "age_at_initial_pathologic_diagnosis", "age_at_diagnosis"))
gender_col <- pick_col(tcga_clin, c("gender", "sex"))
stage_col <- pick_col(tcga_clin, c("stage", "pathologic_stage", "ajcc_pathologic_stage", "clinical_stage"))

plot_meta <- plot_meta %>%
  dplyr::mutate(
    MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High")),
    Age = if (!is.na(age_col)) suppressWarnings(as.numeric(.data[[age_col]])) else NA_real_,
    Gender = if (!is.na(gender_col)) simple_gender(.data[[gender_col]]) else "Unknown",
    Stage = if (!is.na(stage_col)) simple_stage(.data[[stage_col]]) else "Unknown",
    Status = simple_status(status)
  )

cat("\nBefore removing Unknown:\n")
print(table(plot_meta$Gender, useNA = "ifany"))
print(table(plot_meta$Stage, useNA = "ifany"))
print(table(plot_meta$Status, useNA = "ifany"))

# 直接去掉临床注释不完整的样本
plot_meta <- plot_meta %>%
  dplyr::filter(
    is.finite(Age),
    Gender %in% c("Female", "Male"),
    Stage %in% c("I", "II", "III", "IV"),
    Status %in% c("Alive", "Dead")
  ) %>%
  dplyr::arrange(MO_DDRscore_group, MO_DDRscore_raw)

cat("\nAfter removing Unknown:\n")
print(table(plot_meta$Gender, useNA = "ifany"))
print(table(plot_meta$Stage, useNA = "ifany"))
print(table(plot_meta$Status, useNA = "ifany"))
cat("Final heatmap samples:", nrow(plot_meta), "\n")
cat("Samples used:", nrow(plot_meta), "\n")
print(table(plot_meta$MO_DDRscore_group))

############################
# 5. DEG analysis: High vs Low
############################

expr_log <- log2(tcga_expr[, plot_meta$Sample, drop = FALSE] + 1)

group_factor <- factor(plot_meta$MO_DDRscore_group, levels = c("Low", "High"))

design <- model.matrix(~0 + group_factor)
colnames(design) <- c("Low", "High")

fit <- limma::lmFit(expr_log, design)
contrast <- limma::makeContrasts(High - Low, levels = design)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrast))

deg <- limma::topTable(fit2, number = Inf, adjust.method = "BH")
deg$Gene <- rownames(deg)

deg <- deg %>%
  dplyr::mutate(
    Regulation = dplyr::case_when(
      adj.P.Val < DEG_ADJ_P & logFC >= DEG_LOGFC ~ "Up in High",
      adj.P.Val < DEG_ADJ_P & logFC <= -DEG_LOGFC ~ "Down in High",
      TRUE ~ "NS"
    )
  ) %>%
  dplyr::arrange(adj.P.Val)

save_csv(
  deg,
  file.path(OUT_DIR, "Figure2C_DEG_for_heatmap.csv")
)

UP_N <- 30
DOWN_N <- 30

deg_sig <- deg %>%
  dplyr::filter(
    adj.P.Val < DEG_ADJ_P,
    abs(logFC) >= DEG_LOGFC
  )

cat("\nDEG direction counts:\n")
print(table(deg_sig$Regulation))

deg_up_high <- deg_sig %>%
  dplyr::filter(logFC > 0) %>%
  dplyr::arrange(adj.P.Val, dplyr::desc(abs(logFC))) %>%
  dplyr::slice_head(n = UP_N)

deg_down_high <- deg_sig %>%
  dplyr::filter(logFC < 0) %>%
  dplyr::arrange(adj.P.Val, dplyr::desc(abs(logFC))) %>%
  dplyr::slice_head(n = DOWN_N)

heat_deg <- dplyr::bind_rows(
  deg_up_high,
  deg_down_high
) %>%
  dplyr::distinct(Gene, .keep_all = TRUE)

# 如果某一方向基因太少，就保留实际数量，不硬凑
heat_genes <- heat_deg$Gene
heat_genes <- intersect(heat_genes, rownames(expr_log))

save_csv(
  heat_deg,
  file.path(OUT_DIR, "Figure2C_DEG_heatmap_genes_balanced_up_down.csv")
)

cat("Heatmap genes:", length(heat_genes), "\n")
cat("Up in High:", sum(heat_deg$logFC > 0), "\n")
cat("Down in High:", sum(heat_deg$logFC < 0), "\n")

if (length(heat_genes) < 10) {
  heat_genes <- deg %>%
    dplyr::arrange(adj.P.Val) %>%
    dplyr::slice_head(n = TOP_N_GENES) %>%
    dplyr::pull(Gene)
}

heat_genes <- intersect(heat_genes, rownames(expr_log))

save_csv(
  deg %>% dplyr::filter(Gene %in% heat_genes),
  file.path(OUT_DIR, "Figure2C_DEG_heatmap_genes.csv")
)

cat("Heatmap genes:", length(heat_genes), "\n")

############################
# 6. Build heatmap matrix
############################

hm_mat <- expr_log[heat_genes, plot_meta$Sample, drop = FALSE]

hm_mat <- t(scale(t(hm_mat)))
hm_mat[hm_mat > 2] <- 2
hm_mat[hm_mat < -2] <- -2
hm_mat[is.na(hm_mat)] <- 0

############################
# 7. Column annotations - remove unused Unknown
############################
anno_df <- data.frame(
  Group = as.character(plot_meta$MO_DDRscore_group),
  Score = plot_meta$MO_DDRscore_raw,
  Age = plot_meta$Age,
  Gender = plot_meta$Gender,
  Stage = plot_meta$Stage,
  Status = plot_meta$Status,
  row.names = plot_meta$Sample,
  stringsAsFactors = FALSE
)

anno_df$Group <- factor(anno_df$Group, levels = c("Low", "High"))
anno_df$Gender <- factor(anno_df$Gender, levels = c("Female", "Male"))
anno_df$Stage <- factor(anno_df$Stage, levels = intersect(c("I", "II", "III", "IV"), unique(anno_df$Stage)))
anno_df$Status <- factor(anno_df$Status, levels = c("Alive", "Dead"))

top_anno <- HeatmapAnnotation(
  Group = anno_df$Group,
  Score = anno_df$Score,
  Age = anno_df$Age,
  Gender = anno_df$Gender,
  Stage = anno_df$Stage,
  Status = anno_df$Status,
  col = list(
    Group = c(Low = "#3B75AF", High = "#C84630"),
    Score = circlize::colorRamp2(
      quantile(anno_df$Score, c(0.05, 0.5, 0.95), na.rm = TRUE),
      c("#3B75AF", "white", "#C84630")
    ),
    Age = circlize::colorRamp2(
      quantile(anno_df$Age, c(0.05, 0.5, 0.95), na.rm = TRUE),
      c("#E8F1FA", "white", "#D6604D")
    ),
    Gender = c(Female = "#DDA0DD", Male = "#7FB3D5"),
    Stage = c(
      I = "#D9F0D3",
      II = "#ADDD8E",
      III = "#41AB5D",
      IV = "#005A32"
    )[levels(anno_df$Stage)],
    Status = c(Alive = "#4DBBD5", Dead = "#E64B35")
  ),
  annotation_name_side = "left",
  show_annotation_name = TRUE
)
############################
# 8. Draw heatmap
############################

col_fun <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#2B6CB0", "white", "#C53030")
)

ht <- Heatmap(
  hm_mat,
  name = "Expression\nZ-score",
  col = col_fun,
  top_annotation = top_anno,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 6),
  column_split = anno_df$Group,
  heatmap_legend_param = list(
    title = "Expression\nZ-score",
    legend_height = unit(4, "cm")
  ),
  column_title = "Differentially expressed genes between High and Low MO-DDRscore groups",
  column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
  row_title_gp = grid::gpar(fontsize = 9, fontface = "bold")
)
pdf(
  file.path(OUT_DIR, "Figure2C_DEG_heatmap_with_clinical_annotation.pdf"),
  width = 10,
  height = 9,
  useDingbats = FALSE
)

draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legend = TRUE
)

dev.off()

cat("\nDone.\n")
cat("Output:\n")
cat(file.path(OUT_DIR, "Figure2C_DEG_heatmap_with_clinical_annotation.pdf"), "\n")
cat(file.path(OUT_DIR, "Figure2C_DEG_heatmap_genes.csv"), "\n")