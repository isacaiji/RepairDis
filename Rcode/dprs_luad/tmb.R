############################################################
# Add two panels for Figure 2:
#   Fig2A: DDR pathway activity by GSVA/ssGSEA
#   Fig2E: TMB by maftools
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260513)

############################
# 0. Paths
############################

PROJECT_DIR <- "D:/R_workspace/评分/AD_DDR_project"
PROC_DIR <- file.path(PROJECT_DIR, "01_processed")
DATA_DIR <- file.path(PROJECT_DIR, "00_data")
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

DDR_GENE_FILE <- file.path(DATA_DIR, "DDR_236_genes.csv")

MAF_FILE <- file.path(
  TCGA_DIR,
  "mutation",
  paste0(FOCUS_CANCER, ".txt")
)

############################
# 1. Packages
############################

cran_pkgs <- c("data.table", "dplyr", "tidyr", "ggplot2", "ggpubr", "pheatmap")
bioc_pkgs <- c("GSVA", "maftools")

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(pheatmap)
  library(GSVA)
  library(maftools)
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

patient_id <- function(x) {
  substr(gsub("\\.", "-", as.character(x)), 1, 12)
}

save_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file)
}

save_pdf <- function(file, plot, w = 7, h = 5) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = file,
    plot = plot,
    width = w,
    height = h,
    device = "pdf",
    useDingbats = FALSE
  )
}

safe_ssgsea <- function(expr_mat, gene_sets) {
  
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "numeric"
  
  gene_sets <- lapply(gene_sets, function(x) {
    intersect(unique(clean_gene(x)), rownames(expr_mat))
  })
  
  gene_sets <- gene_sets[sapply(gene_sets, length) >= 2]
  
  if (length(gene_sets) < 2) {
    stop("Too few DDR gene sets matched to expression matrix.")
  }
  
  res <- tryCatch(
    {
      GSVA::gsva(
        expr_mat,
        gene_sets,
        method = "ssgsea",
        kcdf = "Gaussian",
        abs.ranking = TRUE,
        verbose = FALSE
      )
    },
    error = function(e) {
      if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
        param <- GSVA::ssgseaParam(
          exprData = expr_mat,
          geneSets = gene_sets,
          normalize = TRUE
        )
        GSVA::gsva(param, verbose = FALSE)
      } else {
        stop(e)
      }
    }
  )
  
  as.matrix(res)
}

############################
# 3. Load expression
############################

TCGA_EXPR_RDS <- file.path(PROC_DIR, paste0(FOCUS_CANCER, "_tcga_expr_tpm_matrix.rds"))

if (file.exists(TCGA_EXPR_RDS)) {
  
  tcga_expr <- readRDS(TCGA_EXPR_RDS)
  
} else {
  
  expr_raw <- fread(TCGA_EXPR_FILE, data.table = FALSE, check.names = FALSE)
  gene_col <- colnames(expr_raw)[1]
  expr_raw[[gene_col]] <- clean_gene(expr_raw[[gene_col]])
  
  expr_raw <- expr_raw %>%
    filter(!is.na(.data[[gene_col]]), .data[[gene_col]] != "") %>%
    group_by(.data[[gene_col]]) %>%
    summarise(
      across(
        everything(),
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

mo_score <- fread(MO_SCORE_FILE, data.table = FALSE, check.names = FALSE)

mo_score <- mo_score %>%
  mutate(
    Sample = gsub("\\.", "-", as.character(Sample)),
    Patient = as.character(Patient),
    SampleType = sample_type(Sample),
    MO_DDRscore_raw = suppressWarnings(as.numeric(MO_DDRscore_raw)),
    MO_DDRscore_group = as.character(MO_DDRscore_group)
  )

tumor_samples <- colnames(tcga_expr)[sample_class(colnames(tcga_expr)) == "Tumor"]

score_df <- mo_score %>%
  filter(
    Sample %in% tumor_samples,
    Sample %in% colnames(tcga_expr),
    MO_DDRscore_group %in% c("Low", "High"),
    is.finite(MO_DDRscore_raw)
  ) %>%
  arrange(Patient, SampleType) %>%
  distinct(Patient, .keep_all = TRUE)

score_df$MO_DDRscore_group <- factor(
  score_df$MO_DDRscore_group,
  levels = c("Low", "High")
)

expr_log <- log2(tcga_expr[, score_df$Sample, drop = FALSE] + 1)

cat("Samples used:\n")
print(table(score_df$MO_DDRscore_group))

############################################################
# 5. Fig2A: DDR pathway activity by ssGSEA
############################################################

ddr_anno <- fread(DDR_GENE_FILE, data.table = FALSE, check.names = FALSE)

# 兼容 Gene / Pathway 列名
if (!"Gene" %in% colnames(ddr_anno)) {
  colnames(ddr_anno)[1] <- "Gene"
}

if (!"Primary_Pathway" %in% colnames(ddr_anno)) {
  pathway_col <- intersect(
    c("Pathway", "pathway", "Primary_pathway", "primary_pathway",
      "PrimaryPathway", "Category", "category", "Type", "type"),
    colnames(ddr_anno)
  )[1]
  
  if (is.na(pathway_col)) {
    stop("DDR gene file must contain Pathway or Primary_Pathway column.")
  } else {
    colnames(ddr_anno)[colnames(ddr_anno) == pathway_col] <- "Primary_Pathway"
  }
}

ddr_anno <- ddr_anno %>%
  mutate(
    Gene = clean_gene(Gene),
    Primary_Pathway = trimws(as.character(Primary_Pathway)),
    Primary_Pathway = case_when(
      Primary_Pathway %in% c("DNA damage checkpoint", "DNA damage checkpoint ") ~ "Checkpoint",
      Primary_Pathway %in% c("Replication stress", "Replication Stress") ~ "Replication_stress",
      TRUE ~ Primary_Pathway
    )
  ) %>%
  filter(
    !is.na(Gene),
    Gene != "",
    !is.na(Primary_Pathway),
    Primary_Pathway != ""
  ) %>%
  distinct(Gene, Primary_Pathway)

ddr_sets <- split(ddr_anno$Gene, ddr_anno$Primary_Pathway)

ddr_sets <- lapply(ddr_sets, function(x) {
  intersect(unique(clean_gene(x)), rownames(expr_log))
})

ddr_sets <- ddr_sets[sapply(ddr_sets, length) >= 2]

cat("\nDDR pathway sets used:\n")
print(sapply(ddr_sets, length))

ddr_ssgsea <- safe_ssgsea(expr_log, ddr_sets)

ddr_long <- as.data.frame(t(ddr_ssgsea), check.names = FALSE)
ddr_long$Sample <- rownames(ddr_long)

ddr_long <- ddr_long %>%
  pivot_longer(
    cols = -Sample,
    names_to = "DDR_Pathway",
    values_to = "Activity"
  ) %>%
  left_join(
    score_df %>% select(Sample, Patient, MO_DDRscore_group, MO_DDRscore_raw),
    by = "Sample"
  )

save_csv(
  ddr_long,
  file.path(OUT_DIR, "Fig2A_DDR_pathway_activity_ssGSEA_long.csv")
)

ddr_stat <- ddr_long %>%
  group_by(DDR_Pathway) %>%
  summarise(
    P = tryCatch(
      wilcox.test(Activity ~ MO_DDRscore_group)$p.value,
      error = function(e) NA_real_
    ),
    Median_High = median(Activity[MO_DDRscore_group == "High"], na.rm = TRUE),
    Median_Low = median(Activity[MO_DDRscore_group == "Low"], na.rm = TRUE),
    Diff = Median_High - Median_Low,
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(P, method = "BH")) %>%
  arrange(P)

save_csv(
  ddr_stat,
  file.path(OUT_DIR, "Fig2A_DDR_pathway_activity_statistics.csv")
)

ddr_long$DDR_Pathway <- factor(
  ddr_long$DDR_Pathway,
  levels = ddr_stat$DDR_Pathway
)

p_ddr_box <- ggplot(
  ddr_long,
  aes(x = MO_DDRscore_group, y = Activity, fill = MO_DDRscore_group)
) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.16, size = 0.45, alpha = 0.45) +
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    size = 3
  ) +
  facet_wrap(~DDR_Pathway, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = NULL,
    y = "ssGSEA score",
    title = "DDR pathway activity"
  )

save_pdf(
  file.path(OUT_DIR, "Fig2A_DDR_pathway_activity_boxplot.pdf"),
  p_ddr_box,
  9,
  7
)

# Heatmap version
ddr_mat <- ddr_ssgsea[, score_df$Sample, drop = FALSE]
ddr_mat <- t(scale(t(ddr_mat)))
ddr_mat[ddr_mat > 2] <- 2
ddr_mat[ddr_mat < -2] <- -2
ddr_mat[is.na(ddr_mat)] <- 0

anno_col <- data.frame(
  Group = score_df$MO_DDRscore_group,
  row.names = score_df$Sample
)

pdf(
  file.path(OUT_DIR, "Fig2A_DDR_pathway_activity_heatmap.pdf"),
  width = 8,
  height = 5,
  useDingbats = FALSE
)

pheatmap::pheatmap(
  ddr_mat,
  annotation_col = anno_col,
  show_colnames = FALSE,
  color = colorRampPalette(c("#2B6CB0", "white", "#C53030"))(100),
  main = "DDR pathway activity by ssGSEA"
)

dev.off()

############################################################
# 6. Fig2E: TMB by maftools official workflow
############################################################

if (!file.exists(MAF_FILE)) {
  stop("MAF_FILE not found: ", MAF_FILE)
}

maf_raw <- fread(MAF_FILE, data.table = FALSE, check.names = FALSE)

if (!all(c("Hugo_Symbol", "Tumor_Sample_Barcode") %in% colnames(maf_raw))) {
  stop(
    "For maftools official TMB workflow, mutation file must contain at least:\n",
    "Hugo_Symbol and Tumor_Sample_Barcode."
  )
}

maf_raw$Hugo_Symbol <- clean_gene(maf_raw$Hugo_Symbol)
maf_raw$Tumor_Sample_Barcode <- patient_id(maf_raw$Tumor_Sample_Barcode)

maf_raw <- maf_raw %>%
  filter(
    !is.na(Hugo_Symbol),
    Hugo_Symbol != "",
    Tumor_Sample_Barcode %in% score_df$Patient
  )

clinical_maf <- score_df %>%
  transmute(
    Tumor_Sample_Barcode = Patient,
    MO_DDRscore_group = as.character(MO_DDRscore_group)
  )

maf_obj <- maftools::read.maf(
  maf = maf_raw,
  clinicalData = clinical_maf,
  verbose = FALSE
)

# captureSize 默认 50Mb，这是 maftools::tmb 默认值。
# 如果你有 panel / WES 的真实可捕获区域大小，可以改 captureSize。
tmb_df <- maftools::tmb(
  maf = maf_obj,
  captureSize = 50,
  logScale = FALSE
)

tmb_df <- as.data.frame(tmb_df, check.names = FALSE)

# 兼容不同版本列名
if (!"Tumor_Sample_Barcode" %in% colnames(tmb_df)) {
  sample_col <- intersect(c("Tumor_Sample_Barcode", "sample", "Sample", "Tumor_Sample"), colnames(tmb_df))[1]
  if (is.na(sample_col)) {
    colnames(tmb_df)[1] <- "Tumor_Sample_Barcode"
  } else {
    colnames(tmb_df)[colnames(tmb_df) == sample_col] <- "Tumor_Sample_Barcode"
  }
}

tmb_df$Tumor_Sample_Barcode <- patient_id(tmb_df$Tumor_Sample_Barcode)

tmb_col <- intersect(c("total_perMB", "total", "TMB"), colnames(tmb_df))[1]

if (is.na(tmb_col)) {
  stop("Cannot identify TMB column from maftools::tmb output.")
}

tmb_plot_df <- tmb_df %>%
  left_join(
    score_df %>%
      transmute(
        Tumor_Sample_Barcode = Patient,
        MO_DDRscore_group,
        MO_DDRscore_raw
      ),
    by = "Tumor_Sample_Barcode"
  ) %>%
  filter(
    !is.na(MO_DDRscore_group),
    is.finite(.data[[tmb_col]])
  )

save_csv(
  tmb_plot_df,
  file.path(OUT_DIR, "Fig2E_TMB_maftools_table.csv")
)

tmb_stat <- tmb_plot_df %>%
  summarise(
    P = wilcox.test(.data[[tmb_col]] ~ MO_DDRscore_group)$p.value,
    Median_High = median(.data[[tmb_col]][MO_DDRscore_group == "High"], na.rm = TRUE),
    Median_Low = median(.data[[tmb_col]][MO_DDRscore_group == "Low"], na.rm = TRUE),
    Diff = Median_High - Median_Low
  )

save_csv(
  tmb_stat,
  file.path(OUT_DIR, "Fig2E_TMB_maftools_statistics.csv")
)

p_tmb <- ggplot(
  tmb_plot_df,
  aes(x = MO_DDRscore_group, y = .data[[tmb_col]], fill = MO_DDRscore_group)
) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.16, size = 0.75, alpha = 0.45) +
  ggpubr::stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    size = 4
  ) +
  scale_fill_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = NULL,
    y = "Tumor mutation burden",
    title = "Tumor mutation burden"
  )

save_pdf(
  file.path(OUT_DIR, "Fig2E_TMB_maftools_boxplot.pdf"),
  p_tmb,
  4.8,
  4.8
)

cat("\nDone.\n")
cat("Outputs:\n")
cat(file.path(OUT_DIR, "Fig2A_DDR_pathway_activity_boxplot.pdf"), "\n")
cat(file.path(OUT_DIR, "Fig2A_DDR_pathway_activity_heatmap.pdf"), "\n")
cat(file.path(OUT_DIR, "Fig2E_TMB_maftools_boxplot.pdf"), "\n")