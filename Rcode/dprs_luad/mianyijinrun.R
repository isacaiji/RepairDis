############################################################
# Figure 2F: CIBERSORT immune infiltration
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

MO_SCORE_FILE <- file.path(
  PROC_DIR,
  paste0(FOCUS_CANCER, "_MO_DDRscore.csv")
)

# CIBERSORT 可信样本过滤
# 如果过滤后样本太少，可以改成 FALSE
USE_CIBERSORT_P_FILTER <- TRUE
CIBERSORT_P_CUTOFF <- 0.05

############################
# 1. Packages
############################

pkgs <- c("data.table", "dplyr", "tidyr", "ggplot2", "ggpubr", "pheatmap", "RColorBrewer")

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
  library(ggpubr)
  library(pheatmap)
  library(RColorBrewer)
})

# 需要你本地已经安装
suppressPackageStartupMessages({
  library(bseqsc)
  library(CIBERSORT)
})

data(LM22)

############################
# 2. Helper functions
############################

clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub(" ///.*$| //.*$|;.*$|,.*$|\\s*\\([^\\)]*\\)$", "", x)
  x[x %in% c("", "NA", "---", "NULL", "N/A", "NAN")] <- NA
  x
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
# 4. Load MO-DDRscore group
############################

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

score_df <- mo_score %>%
  dplyr::filter(
    Sample %in% tumor_samples,
    Sample %in% colnames(tcga_expr),
    MO_DDRscore_group %in% c("Low", "High"),
    is.finite(MO_DDRscore_raw)
  ) %>%
  dplyr::arrange(Patient, SampleType) %>%
  dplyr::distinct(Patient, .keep_all = TRUE)

score_df$MO_DDRscore_group <- factor(
  score_df$MO_DDRscore_group,
  levels = c("Low", "High")
)

cat("Samples for CIBERSORT:", nrow(score_df), "\n")
print(table(score_df$MO_DDRscore_group))

############################
# 5. Prepare CIBERSORT input
############################

# CIBERSORT 输入：gene x sample 表达矩阵
# 用 TPM 原始值，不用 log2(TPM+1)
ciber_expr <- tcga_expr[, score_df$Sample, drop = FALSE]

# 去掉低表达和常数基因，减少 CIBERSORT 报错
ciber_expr <- ciber_expr[rowMeans(ciber_expr, na.rm = TRUE) > 0, , drop = FALSE]
ciber_expr <- ciber_expr[apply(ciber_expr, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]

# 保留 LM22 中出现过的基因，提高匹配效率
lm22_genes <- clean_gene(rownames(LM22))
rownames(LM22) <- lm22_genes

common_genes <- intersect(rownames(ciber_expr), rownames(LM22))
ciber_expr <- ciber_expr[common_genes, , drop = FALSE]

cat("Genes matched with LM22:", nrow(ciber_expr), "\n")

# CIBERSORT 通常要求 gene symbol 行名
ciber_input <- as.data.frame(ciber_expr, check.names = FALSE)
ciber_input$GeneSymbol <- rownames(ciber_input)
ciber_input <- ciber_input[, c("GeneSymbol", setdiff(colnames(ciber_input), "GeneSymbol"))]

input_file <- file.path(OUT_DIR, "Fig2F_CIBERSORT_input_TPM.txt")

data.table::fwrite(
  ciber_input,
  input_file,
  sep = "\t",
  quote = FALSE
)

############################
# 6. Run CIBERSORT
############################

# 参考教程中使用 cibersort(sig_matrix = LM22, mixture_file = FPKM)
# 这里 mixture_file 使用 TPM 矩阵
source("D:/R_workspace/评分/AD_DDR_project/01-script/Cibersort.R")



DATA_DIR <- "D:/R_workspace/评分/AD_DDR_project/00_data"
LM22_FILE <- file.path(DATA_DIR, "LM22.txt")

lm22_df <- data.frame(
  GeneSymbol = rownames(LM22),
  LM22,
  check.names = FALSE
)

data.table::fwrite(
  lm22_df,
  LM22_FILE,
  sep = "\t",
  quote = FALSE
)

file.exists(LM22_FILE)


ciber_res <- CIBERSORT(
  sig_matrix = LM22_FILE,
  mixture_file = input_file,
  perm = 1000,
  QN = FALSE
)

# 某些版本返回 matrix，某些返回 data.frame
ciber_res <- as.data.frame(ciber_res, check.names = FALSE)
ciber_res$Sample <- rownames(ciber_res)

save_csv(
  ciber_res,
  file.path(OUT_DIR, "Fig2F_CIBERSORT_raw_result.csv")
)

############################
# 7. Merge group
############################

# CIBERSORT 前 22 列是免疫细胞比例，后面通常是 P-value / Correlation / RMSE
cell_cols <- setdiff(
  colnames(ciber_res),
  c("Sample", "P-value", "P.value", "Correlation", "RMSE")
)

# 只保留数值型细胞列
cell_cols <- cell_cols[sapply(ciber_res[, cell_cols, drop = FALSE], is.numeric)]

ciber_merged <- ciber_res %>%
  dplyr::left_join(
    score_df %>%
      dplyr::select(Sample, Patient, MO_DDRscore_group, MO_DDRscore_raw),
    by = "Sample"
  ) %>%
  dplyr::filter(!is.na(MO_DDRscore_group))

# 标准化 P-value 列名
if ("P-value" %in% colnames(ciber_merged)) {
  colnames(ciber_merged)[colnames(ciber_merged) == "P-value"] <- "CIBERSORT_P"
}
if ("P.value" %in% colnames(ciber_merged)) {
  colnames(ciber_merged)[colnames(ciber_merged) == "P.value"] <- "CIBERSORT_P"
}

save_csv(
  ciber_merged,
  file.path(OUT_DIR, "Fig2F_CIBERSORT_merged_with_group.csv")
)

# 是否按 CIBERSORT P-value 过滤可信样本
if (USE_CIBERSORT_P_FILTER && "CIBERSORT_P" %in% colnames(ciber_merged)) {
  
  ciber_plot_df <- ciber_merged %>%
    dplyr::filter(CIBERSORT_P < CIBERSORT_P_CUTOFF)
  
  cat("Samples after CIBERSORT P filter:", nrow(ciber_plot_df), "\n")
  print(table(ciber_plot_df$MO_DDRscore_group))
  
  if (nrow(ciber_plot_df) < 50) {
    message("Too few samples after CIBERSORT P filter. Use all samples instead.")
    ciber_plot_df <- ciber_merged
  }
  
} else {
  
  ciber_plot_df <- ciber_merged
}

############################
# 8. Long table and statistics
############################

ciber_long <- ciber_plot_df %>%
  dplyr::select(
    Sample,
    Patient,
    MO_DDRscore_group,
    MO_DDRscore_raw,
    dplyr::all_of(cell_cols)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(cell_cols),
    names_to = "Celltype",
    values_to = "Composition"
  )

save_csv(
  ciber_long,
  file.path(OUT_DIR, "Fig2F_CIBERSORT_long_table.csv")
)

ciber_stat <- ciber_long %>%
  dplyr::group_by(Celltype) %>%
  dplyr::summarise(
    P = tryCatch(
      wilcox.test(Composition ~ MO_DDRscore_group)$p.value,
      error = function(e) NA_real_
    ),
    Median_High = median(Composition[MO_DDRscore_group == "High"], na.rm = TRUE),
    Median_Low = median(Composition[MO_DDRscore_group == "Low"], na.rm = TRUE),
    Diff = Median_High - Median_Low,
    .groups = "drop"
  ) %>%
  dplyr::mutate(FDR = p.adjust(P, method = "BH")) %>%
  dplyr::arrange(P)

save_csv(
  ciber_stat,
  file.path(OUT_DIR, "Fig2F_CIBERSORT_statistics.csv")
)

############################
# 9. Boxplot: 22 immune cells
############################

ciber_long$Celltype <- factor(
  ciber_long$Celltype,
  levels = ciber_stat$Celltype
)

p_box <- ggplot(
  ciber_long,
  aes(x = Celltype, y = Composition, fill = MO_DDRscore_group)
) +
  geom_boxplot(
    position = position_dodge(0.75),
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    aes(color = MO_DDRscore_group),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    size = 0.35,
    alpha = 0.35,
    show.legend = FALSE
  ) +
  ggpubr::stat_compare_means(
    aes(group = MO_DDRscore_group),
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE,
    size = 3
  ) +
  scale_fill_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
  scale_color_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = NULL,
    y = "Estimated proportion",
    fill = "Group",
    title = "CIBERSORT immune cell infiltration"
  )

safe_ggsave(
  file.path(OUT_DIR, "Fig2F_CIBERSORT_22cell_boxplot.pdf"),
  p_box,
  11,
  5.5
)

############################
# 10. Stacked barplot
############################

ciber_bar <- ciber_long %>%
  dplyr::arrange(MO_DDRscore_group, MO_DDRscore_raw)

ciber_bar$Sample <- factor(
  ciber_bar$Sample,
  levels = unique(ciber_bar$Sample)
)

my_palette <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(length(cell_cols))

p_bar <- ggplot(
  ciber_bar,
  aes(x = Sample, y = Composition, fill = Celltype)
) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = my_palette) +
  facet_grid(~MO_DDRscore_group, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = NULL,
    y = "Estimated proportion",
    fill = "Celltype",
    title = "Immune cell composition by CIBERSORT"
  )

safe_ggsave(
  file.path(OUT_DIR, "Fig2F_CIBERSORT_stacked_barplot.pdf"),
  p_bar,
  12,
  6
)

############################
# 11. Heatmap
############################

heat_mat <- as.matrix(ciber_plot_df[, cell_cols, drop = FALSE])
rownames(heat_mat) <- ciber_plot_df$Sample

heat_mat_t <- t(heat_mat)
heat_mat_t <- t(scale(t(heat_mat_t)))
heat_mat_t[is.na(heat_mat_t)] <- 0
heat_mat_t[heat_mat_t > 2] <- 2
heat_mat_t[heat_mat_t < -2] <- -2

anno_col <- data.frame(
  Group = ciber_plot_df$MO_DDRscore_group,
  row.names = ciber_plot_df$Sample
)

pdf(
  file.path(OUT_DIR, "Fig2F_CIBERSORT_22cell_heatmap.pdf"),
  width = 10,
  height = 6,
  useDingbats = FALSE
)

pheatmap::pheatmap(
  heat_mat_t,
  annotation_col = anno_col,
  show_colnames = FALSE,
  color = colorRampPalette(c("#2B6CB0", "white", "#C53030"))(100),
  main = "CIBERSORT immune cell infiltration"
)

dev.off()

############################
# 12. Done
############################

cat("\nDone.\n")
cat("Main outputs:\n")
cat(file.path(OUT_DIR, "Fig2F_CIBERSORT_22cell_boxplot.pdf"), "\n")
cat(file.path(OUT_DIR, "Fig2F_CIBERSORT_stacked_barplot.pdf"), "\n")
cat(file.path(OUT_DIR, "Fig2F_CIBERSORT_22cell_heatmap.pdf"), "\n")
cat(file.path(OUT_DIR, "Fig2F_CIBERSORT_statistics.csv"), "\n")

