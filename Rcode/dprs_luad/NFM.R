############################################################
# Quick NMF DDR subtype vs MO-DDRscore
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260513)

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(NMF)

############################################################
# 1. 路径
############################################################

PROJECT_DIR <- "D:/R_workspace/评分/AD_DDR_project"
PROC_DIR <- file.path(PROJECT_DIR, "01_processed")
DATA_DIR <- file.path(PROJECT_DIR, "00_data")
OUT_DIR <- file.path(PROJECT_DIR, "02_Figure2_MultiOmics_Immune", "NMF_DDR_subtype_quick")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

EXPR_FILE  <- file.path(PROC_DIR, "LUAD_tcga_expr_tpm_matrix.rds")
SCORE_FILE <- file.path(PROC_DIR, "LUAD_MO_DDRscore.csv")

# 改成你的 236 DDR 基因文件
DDR_FILE <- file.path(DATA_DIR, "DDR_236_genes.csv")

############################################################
# 2. 小函数
############################################################

clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub(" ///.*$| //.*$|;.*$|,.*$|\\s*\\([^\\)]*\\)$", "", x)
  x[x %in% c("", "NA", "---", "NULL", "N/A")] <- NA
  x
}

patient_id <- function(x) substr(gsub("\\.", "-", x), 1, 12)

sample_type <- function(x) substr(gsub("\\.", "-", x), 14, 15)

is_tumor <- function(x) sample_type(x) %in% c("01", "02", "03", "05", "06", "07")

calc_cramers_v <- function(tab) {
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
  n <- sum(tab)
  sqrt(as.numeric(chi$statistic) / (n * min(nrow(tab) - 1, ncol(tab) - 1)))
}

############################################################
# 3. 读取表达矩阵和 MO-DDRscore
############################################################

expr <- readRDS(EXPR_FILE)
expr <- as.matrix(expr)
storage.mode(expr) <- "numeric"

rownames(expr) <- clean_gene(rownames(expr))
colnames(expr) <- gsub("\\.", "-", colnames(expr))

expr <- expr[!is.na(rownames(expr)) & rownames(expr) != "", , drop = FALSE]
expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]

score <- fread(SCORE_FILE, data.table = FALSE, check.names = FALSE)
score$Sample <- gsub("\\.", "-", score$Sample)

if (!"Patient" %in% colnames(score)) {
  score$Patient <- patient_id(score$Sample)
}

if (!"MO_DDRscore_raw" %in% colnames(score)) {
  if ("MO_DDRscore" %in% colnames(score)) {
    score$MO_DDRscore_raw <- score$MO_DDRscore
  } else {
    stop("score 文件里找不到 MO_DDRscore_raw 或 MO_DDRscore")
  }
}

score$MO_DDRscore_raw <- as.numeric(score$MO_DDRscore_raw)
score$MO_DDRscore_group <- factor(score$MO_DDRscore_group, levels = c("Low", "High"))

score_df <- score %>%
  filter(
    Sample %in% colnames(expr),
    is_tumor(Sample),
    MO_DDRscore_group %in% c("Low", "High"),
    is.finite(MO_DDRscore_raw)
  ) %>%
  arrange(Patient, Sample) %>%
  distinct(Patient, .keep_all = TRUE)

cat("Samples used:", nrow(score_df), "\n")
print(table(score_df$MO_DDRscore_group))

expr_log <- log2(expr[, score_df$Sample, drop = FALSE] + 1)

############################################################
# 4. 读取 DDR genes，并准备 NMF 矩阵
############################################################

ddr <- fread(DDR_FILE, data.table = FALSE, check.names = FALSE)
ddr_genes <- clean_gene(ddr[[1]])
ddr_genes <- unique(ddr_genes[!is.na(ddr_genes)])

matched_genes <- intersect(ddr_genes, rownames(expr_log))

cat("DDR genes loaded:", length(ddr_genes), "\n")
cat("DDR genes matched:", length(matched_genes), "\n")

if (length(matched_genes) < 30) {
  stop("匹配到的 DDR genes 太少，请检查 DDR_FILE")
}

ddr_expr <- expr_log[matched_genes, , drop = FALSE]

# 过滤低表达
ddr_expr <- ddr_expr[rowMeans(ddr_expr, na.rm = TRUE) > 0.1, , drop = FALSE]

# 选 MAD 最高的 150 个 DDR genes
gene_mad <- apply(ddr_expr, 1, mad, na.rm = TRUE)
top_genes <- names(sort(gene_mad, decreasing = TRUE))[1:min(236, length(gene_mad))]

nmf_mat <- ddr_expr[top_genes, , drop = FALSE]
nmf_mat[nmf_mat < 0] <- 0
nmf_mat[!is.finite(nmf_mat)] <- 0

cat("Genes used for NMF:", nrow(nmf_mat), "\n")
cat("Samples used for NMF:", ncol(nmf_mat), "\n")

fwrite(
  data.frame(Gene = top_genes, MAD = gene_mad[top_genes]),
  file.path(OUT_DIR, "NMF_used_DDR_genes.csv")
)

############################################################
# 5. NMF k = 2
############################################################

nmf_fit <- nmf(
  nmf_mat,
  rank = 2,
  method = "brunet",
  nrun = 50,
  seed = 20260513
)

# 稳定提取样本分型：用 H 矩阵最大成分归类
best_fit <- NMF::fit(nmf_fit)
H <- as.matrix(NMF::coef(best_fit))

if (is.null(colnames(H))) {
  colnames(H) <- colnames(nmf_mat)
}

raw_cluster <- apply(H, 2, which.max)
raw_cluster <- as.character(raw_cluster[colnames(nmf_mat)])

nmf_df <- data.frame(
  Sample = colnames(nmf_mat),
  NMF_cluster_raw = raw_cluster,
  stringsAsFactors = FALSE
) %>%
  left_join(
    score_df %>% select(Sample, Patient, MO_DDRscore_raw, MO_DDRscore_group),
    by = "Sample"
  )

# 按 MO-DDRscore 中位数重命名：低的是 C1，高的是 C2
cluster_order <- nmf_df %>%
  group_by(NMF_cluster_raw) %>%
  summarise(MedianScore = median(MO_DDRscore_raw), .groups = "drop") %>%
  arrange(MedianScore)

map <- setNames(c("C1", "C2"), cluster_order$NMF_cluster_raw)

nmf_df$NMF_subtype <- factor(map[nmf_df$NMF_cluster_raw], levels = c("C1", "C2"))

fwrite(nmf_df, file.path(OUT_DIR, "DDR_NMF_subtype_with_MO_DDRscore.csv"))

cat("\nNMF subtype counts:\n")
print(table(nmf_df$NMF_subtype))

cat("\nNMF subtype vs MO-DDRscore group:\n")
print(table(nmf_df$NMF_subtype, nmf_df$MO_DDRscore_group))

############################################################
# 6. 图 1：NMF consensus heatmap
# 保持原来 pheatmap 图形样式，只修改上方横条
############################################################

cons <- consensus(nmf_fit)

if (is.null(rownames(cons))) rownames(cons) <- colnames(nmf_mat)
if (is.null(colnames(cons))) colnames(cons) <- colnames(nmf_mat)

# 样本排序：先按 C1/C2，再按 MO-DDRscore 从低到高
sample_order <- nmf_df %>%
  arrange(NMF_subtype, MO_DDRscore_raw) %>%
  pull(Sample)

cons <- cons[sample_order, sample_order, drop = FALSE]

# 这里只保留两个横条：
# 1. NMF_subtype
# 2. MO_DDRscore_group
anno <- nmf_df %>%
  select(Sample, NMF_subtype, MO_DDRscore_group) %>%
  distinct(Sample, .keep_all = TRUE)

rownames(anno) <- anno$Sample
anno$Sample <- NULL

anno <- anno[sample_order, , drop = FALSE]

anno$NMF_subtype <- factor(anno$NMF_subtype, levels = c("C1", "C2"))
anno$MO_DDRscore_group <- factor(anno$MO_DDRscore_group, levels = c("Low", "High"))

anno_colors <- list(
  NMF_subtype = c(
    C1 = "#2874C5",
    C2 = "#C6524A"
  ),
  MO_DDRscore_group = c(
    Low = "#2874C5",
    High = "#C6524A"
  )
)

pdf(
  file.path(OUT_DIR, "Fig2A_NMF_consensus_heatmap_clean_same_style.pdf"),
  width = 6.5,
  height = 6,
  useDingbats = FALSE
)

pheatmap(
  cons,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = anno,
  annotation_row = anno,
  annotation_colors = anno_colors,
  color = colorRampPalette(c("white", "#08306B"))(100),
  border_color = NA,
  main = "NMF consensus matrix, k = 2"
)

dev.off()
############################################################
# 7. 图 2：C1/C2 的 MO-DDRscore 差异
############################################################

p1 <- ggplot(nmf_df, aes(x = NMF_subtype, y = MO_DDRscore_raw, fill = NMF_subtype)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.45) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = c(C1 = "#3B75AF", C2 = "#C84630")) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "NMF-derived DDR subtype",
    y = "MO-DDRscore",
    title = "MO-DDRscore across NMF-derived DDR subtypes"
  )

ggsave(
  file.path(OUT_DIR, "Fig2B_MO_DDRscore_by_NMF_subtype.pdf"),
  p1,
  width = 4.8,
  height = 4.3,
  useDingbats = FALSE
)

############################################################
# 8. 图 3：NMF subtype vs MO-DDRscore High/Low
############################################################

tab <- table(nmf_df$NMF_subtype, nmf_df$MO_DDRscore_group)

fisher_p <- fisher.test(tab)$p.value
chisq_p <- suppressWarnings(chisq.test(tab, correct = FALSE)$p.value)
cramer_v <- calc_cramers_v(tab)

prop_df <- nmf_df %>%
  count(NMF_subtype, MO_DDRscore_group) %>%
  group_by(NMF_subtype) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup()

p2 <- ggplot(prop_df, aes(x = NMF_subtype, y = Proportion, fill = MO_DDRscore_group)) +
  geom_col(width = 0.65, color = "white") +
  geom_text(
    aes(label = paste0(n, "\n", scales::percent(Proportion, accuracy = 0.1))),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 3.2
  ) +
  scale_fill_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "NMF-derived DDR subtype",
    y = "Proportion",
    fill = "MO-DDRscore group",
    title = paste0("Fisher P = ", signif(fisher_p, 3),
                   "; Cramer's V = ", signif(cramer_v, 3))
  )

ggsave(
  file.path(OUT_DIR, "Fig2C_NMF_subtype_vs_MO_DDRscore_group.pdf"),
  p2,
  width = 5.2,
  height = 4.5,
  useDingbats = FALSE
)

############################################################
# 9. 输出 summary
############################################################

wilcox_p <- wilcox.test(MO_DDRscore_raw ~ NMF_subtype, data = nmf_df)$p.value

summary_df <- data.frame(
  Item = c(
    "N_samples",
    "N_DDR_genes_loaded",
    "N_DDR_genes_matched",
    "N_genes_used_for_NMF",
    "C1_n",
    "C2_n",
    "Wilcox_P_MO_DDRscore_C1_vs_C2",
    "Fisher_P_NMF_vs_MO_DDRscore_group",
    "ChiSquare_P_NMF_vs_MO_DDRscore_group",
    "Cramers_V"
  ),
  Value = c(
    nrow(nmf_df),
    length(ddr_genes),
    length(matched_genes),
    nrow(nmf_mat),
    sum(nmf_df$NMF_subtype == "C1"),
    sum(nmf_df$NMF_subtype == "C2"),
    signif(wilcox_p, 5),
    signif(fisher_p, 5),
    signif(chisq_p, 5),
    signif(cramer_v, 5)
  )
)

fwrite(summary_df, file.path(OUT_DIR, "NMF_MO_DDRscore_summary.csv"))










############################################################
# 10. 生存分析：NMF subtype C1/C2
############################################################

library(survival)
library(survminer)

# 临床文件路径：优先用你项目里处理好的临床文件
CLINICAL_CANDIDATES <- c(
  file.path(PROC_DIR, "LUAD_clinical_processed.csv"),
  file.path(PROC_DIR, "LUAD_clinical.csv"),
  "D:/R/R_workspace/梁老师文件/TCGA/clinical/TCGA.LUAD.sampleMap/LUAD_clinicalMatrix",
  "D:/R/R_workspace/梁老师文件/TCGA/clinical/TCGA.LUAD.sampleMap/LUAD_clinicalMatrix.txt",
  "D:/R/R_workspace/梁老师文件/TCGA/clinical/TCGA.LUAD.sampleMap/LUAD_clinicalMatrix.tsv"
)

CLIN_FILE <- CLINICAL_CANDIDATES[file.exists(CLINICAL_CANDIDATES)][1]

if (is.na(CLIN_FILE)) {
  stop("没有找到临床文件，请手动设置 CLIN_FILE 路径")
}

cat("Clinical file:\n", CLIN_FILE, "\n")

clinical <- fread(CLIN_FILE, data.table = FALSE, check.names = FALSE)

# 自动识别病人 ID 列
id_col <- intersect(
  c("Patient", "patient", "sample", "Sample", "bcr_patient_barcode", "submitter_id"),
  colnames(clinical)
)[1]

if (is.na(id_col)) {
  id_col <- colnames(clinical)[1]
}

clinical$Patient <- substr(gsub("\\.", "-", as.character(clinical[[id_col]])), 1, 12)

# 自动处理生存状态
# TCGA clinicalMatrix 常见 vital_status: LIVING / DECEASED
if ("vital_status" %in% colnames(clinical)) {
  clinical$status <- ifelse(
    toupper(as.character(clinical$vital_status)) %in% c("DECEASED", "DEAD", "1"),
    1,
    0
  )
} else if ("status" %in% colnames(clinical)) {
  clinical$status <- as.numeric(clinical$status)
} else {
  stop("临床文件里找不到 vital_status 或 status")
}

# 自动处理生存时间
# OS time = 死亡者 days_to_death，否则 days_to_last_followup
if (all(c("days_to_death", "days_to_last_followup") %in% colnames(clinical))) {
  
  clinical$days_to_death <- suppressWarnings(as.numeric(clinical$days_to_death))
  clinical$days_to_last_followup <- suppressWarnings(as.numeric(clinical$days_to_last_followup))
  
  clinical$time <- ifelse(
    clinical$status == 1,
    clinical$days_to_death,
    clinical$days_to_last_followup
  )
  
} else if ("time" %in% colnames(clinical)) {
  
  clinical$time <- suppressWarnings(as.numeric(clinical$time))
  
} else if ("OS.time" %in% colnames(clinical)) {
  
  clinical$time <- suppressWarnings(as.numeric(clinical$OS.time))
  
} else {
  
  stop("临床文件里找不到 days_to_death / days_to_last_followup / time / OS.time")
}

clinical_use <- clinical %>%
  dplyr::select(Patient, time, status) %>%
  dplyr::filter(is.finite(time), time > 0, status %in% c(0, 1)) %>%
  dplyr::distinct(Patient, .keep_all = TRUE)

# 合并 NMF subtype
surv_df <- nmf_df %>%
  dplyr::select(Patient, Sample, NMF_subtype, MO_DDRscore_raw, MO_DDRscore_group) %>%
  dplyr::left_join(clinical_use, by = "Patient") %>%
  dplyr::filter(is.finite(time), time > 0, status %in% c(0, 1))

surv_df$NMF_subtype <- factor(surv_df$NMF_subtype, levels = c("C1", "C2"))

cat("Survival samples:", nrow(surv_df), "\n")
print(table(surv_df$NMF_subtype))
print(table(surv_df$status))

fwrite(
  surv_df,
  file.path(OUT_DIR, "NMF_subtype_survival_input.csv")
)

############################################################
# 10.1 KM 曲线
############################################################

fit <- survfit(Surv(time, status) ~ NMF_subtype, data = surv_df)

cox_fit <- coxph(Surv(time, status) ~ NMF_subtype, data = surv_df)
cox_sum <- summary(cox_fit)

hr <- cox_sum$coefficients[1, "exp(coef)"]
hr_low <- cox_sum$conf.int[1, "lower .95"]
hr_high <- cox_sum$conf.int[1, "upper .95"]
cox_p <- cox_sum$coefficients[1, "Pr(>|z|)"]

cox_res <- data.frame(
  Comparison = "C2 vs C1",
  HR = hr,
  CI95_low = hr_low,
  CI95_high = hr_high,
  Cox_P = cox_p,
  stringsAsFactors = FALSE
)

fwrite(
  cox_res,
  file.path(OUT_DIR, "NMF_subtype_survival_cox_result.csv")
)

p_km <- ggsurvplot(
  fit,
  data = surv_df,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = c("#3B75AF", "#C84630"),
  legend.title = "DDR subtype",
  legend.labs = c("C1", "C2"),
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  risk.table.height = 0.25,
  ggtheme = theme_bw(base_size = 12)
)

p_km$plot <- p_km$plot +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  ggtitle(
    paste0(
      "Overall survival by NMF-derived DDR subtype\n",
      "HR(C2 vs C1) = ", signif(hr, 3),
      " [", signif(hr_low, 3), "-", signif(hr_high, 3), "]"
    )
  )

pdf(
  file.path(OUT_DIR, "Fig2D_NMF_subtype_overall_survival_KM.pdf"),
  width = 6.2,
  height = 6.0,
  useDingbats = FALSE
)
print(p_km)
dev.off()

############################################################
# 10.2 可选：MO-DDRscore High/Low 的 KM，放补图
############################################################

fit_score <- survfit(Surv(time, status) ~ MO_DDRscore_group, data = surv_df)

cox_score <- coxph(Surv(time, status) ~ MO_DDRscore_group, data = surv_df)
cox_score_sum <- summary(cox_score)

score_hr <- cox_score_sum$coefficients[1, "exp(coef)"]
score_hr_low <- cox_score_sum$conf.int[1, "lower .95"]
score_hr_high <- cox_score_sum$conf.int[1, "upper .95"]
score_p <- cox_score_sum$coefficients[1, "Pr(>|z|)"]

score_cox_res <- data.frame(
  Comparison = "MO-DDRscore High vs Low",
  HR = score_hr,
  CI95_low = score_hr_low,
  CI95_high = score_hr_high,
  Cox_P = score_p,
  stringsAsFactors = FALSE
)

fwrite(
  score_cox_res,
  file.path(OUT_DIR, "MO_DDRscore_group_survival_cox_result.csv")
)

p_km_score <- ggsurvplot(
  fit_score,
  data = surv_df,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = c("#3B75AF", "#C84630"),
  legend.title = "MO-DDRscore group",
  legend.labs = c("Low", "High"),
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  risk.table.height = 0.25,
  ggtheme = theme_bw(base_size = 12)
)

p_km_score$plot <- p_km_score$plot +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  ggtitle(
    paste0(
      "Overall survival by MO-DDRscore group\n",
      "HR(High vs Low) = ", signif(score_hr, 3),
      " [", signif(score_hr_low, 3), "-", signif(score_hr_high, 3), "]"
    )
  )

pdf(
  file.path(OUT_DIR, "Supplementary_MO_DDRscore_group_overall_survival_KM.pdf"),
  width = 6.2,
  height = 6.0,
  useDingbats = FALSE
)
print(p_km_score)
dev.off()

############################################################
# 10.3 输出生存分析 summary
############################################################

surv_summary <- rbind(
  cox_res %>% dplyr::mutate(Module = "NMF_subtype"),
  score_cox_res %>% dplyr::mutate(Module = "MO_DDRscore_group")
)

fwrite(
  surv_summary,
  file.path(OUT_DIR, "Survival_analysis_summary.csv")
)

cat("\nSurvival analysis finished.\n")
cat("Main survival plot:\n")
cat(file.path(OUT_DIR, "Fig2D_NMF_subtype_overall_survival_KM.pdf"), "\n")
cat("Survival summary:\n")
cat(file.path(OUT_DIR, "Survival_analysis_summary.csv"), "\n")




cat("\nDone.\n")
cat("Output:", OUT_DIR, "\n")
print(summary_df)