############################################################
# Figure 2A: MO-DDRscore landscape
# Samples ranked by MO-DDRscore
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260513)

############################
# 0. Parameters
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

MUT_FILE <- file.path(TCGA_DIR, "mutation", paste0(FOCUS_CANCER, ".txt"))
CNV_FILE <- file.path(TCGA_DIR, "cnv", paste0(FOCUS_CANCER, ".txt"))

DDR_GENE_FILE <- file.path(DATA_DIR, "DDR_gene_pathway.csv")

CNV_AMP_TH <- 0.20
CNV_DEL_TH <- -0.20

############################
# 1. Packages
############################

pkgs <- c("data.table", "dplyr", "tidyr", "ggplot2", "patchwork", "GSVA")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p == "GSVA") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install("GSVA", ask = FALSE, update = FALSE)
    } else {
      install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(GSVA)
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

zscore_vec <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x)) || is.na(sd(x, na.rm = TRUE)) || sd(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

save_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file)
}

safe_ssgsea <- function(expr_mat, gene_sets) {
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "numeric"
  
  gene_sets <- lapply(gene_sets, function(x) {
    intersect(unique(clean_gene(x)), rownames(expr_mat))
  })
  gene_sets <- gene_sets[sapply(gene_sets, length) >= 2]
  
  if (length(gene_sets) == 0) {
    stop("No valid gene sets after matching expression matrix.")
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
        param <- GSVA::ssgseaParam(expr_mat, gene_sets, normalize = TRUE)
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
# 4. Load MO-DDRscore and define samples
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
    is.finite(MO_DDRscore_raw)
  ) %>%
  dplyr::arrange(Patient, SampleType) %>%
  dplyr::distinct(Patient, .keep_all = TRUE)

# 如果文件里没有 High/Low 分组，自动按中位数重新分
if (!"MO_DDRscore_group" %in% colnames(score_df) ||
    any(!score_df$MO_DDRscore_group %in% c("Low", "High"))) {
  
  cutoff <- median(score_df$MO_DDRscore_raw, na.rm = TRUE)
  
  score_df$MO_DDRscore_group <- ifelse(
    score_df$MO_DDRscore_raw >= cutoff,
    "High",
    "Low"
  )
}

score_df$MO_DDRscore_group <- factor(
  score_df$MO_DDRscore_group,
  levels = c("Low", "High")
)

expr_log <- log2(tcga_expr[, score_df$Sample, drop = FALSE] + 1)

cat("Samples used:", nrow(score_df), "\n")
print(table(score_df$MO_DDRscore_group))

############################
# 5. DDR activity
############################

default_ddr_sets <- list(
  HR = c("BRCA1","BRCA2","RAD51","RAD51B","RAD51C","RAD51D","RAD52","RAD54L","BARD1","PALB2","BRIP1","XRCC2","XRCC3","ATM","ATR","CHEK1","CHEK2","MRE11","RAD50","NBN","BLM"),
  NER = c("XPA","ERCC3","XPC","ERCC2","DDB1","DDB2","ERCC4","ERCC5","ERCC1","LIG1","PCNA","CETN2","RAD23B"),
  BER = c("OGG1","MUTYH","NTHL1","NEIL1","NEIL2","NEIL3","UNG","APEX1","APEX2","POLB","LIG3","XRCC1","PARP1","PARP2"),
  MMR = c("MLH1","MLH3","MSH2","MSH3","MSH6","PMS1","PMS2","EXO1","PCNA","RPA1","RPA2","RPA3"),
  NHEJ = c("XRCC5","XRCC6","PRKDC","DCLRE1C","LIG4","XRCC4","NHEJ1","MRE11","RAD50","NBN"),
  FA = c("FANCA","FANCB","FANCC","FANCD2","FANCE","FANCF","FANCG","FANCI","FANCL","FANCM","BRCA2","PALB2","BRIP1","RAD51C"),
  Checkpoint = c("ATM","ATR","ATRIP","CHEK1","CHEK2","TP53","MDM2","WEE1","CLSPN","HUS1","RAD1","RAD9A","TOPBP1","MDC1","H2AFX"),
  Replication_Stress = c("ATR","ATRIP","CHEK1","CLSPN","TIMELESS","TIPIN","RPA1","RPA2","RPA3","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","CDC45")
)

if (file.exists(DDR_GENE_FILE)) {
  ddr_anno <- data.table::fread(DDR_GENE_FILE, data.table = FALSE, check.names = FALSE)
  colnames(ddr_anno)[1] <- "Gene"
  if (!"Primary_Pathway" %in% colnames(ddr_anno)) {
    ddr_anno$Primary_Pathway <- "DDR"
  }
  
  ddr_anno <- ddr_anno %>%
    dplyr::mutate(
      Gene = clean_gene(Gene),
      Primary_Pathway = as.character(Primary_Pathway)
    ) %>%
    dplyr::filter(!is.na(Gene), Gene != "")
  
  ddr_sets <- split(ddr_anno$Gene, ddr_anno$Primary_Pathway)
} else {
  ddr_sets <- default_ddr_sets
}

ddr_sets <- lapply(ddr_sets, function(x) {
  intersect(unique(clean_gene(x)), rownames(expr_log))
})
ddr_sets <- ddr_sets[sapply(ddr_sets, length) >= 2]

ddr_ssgsea <- safe_ssgsea(expr_log, ddr_sets)
DDR_Activity <- colMeans(ddr_ssgsea, na.rm = TRUE)

############################
# 6. Immune activity
############################

immune_sets <- list(
  CD8_T_cell = c("CD8A","CD8B","GZMB","GZMH","PRF1","IFNG","NKG7"),
  Cytolytic_activity = c("GZMA","GZMB","PRF1","GNLY","NKG7"),
  IFN_gamma_response = c("IFNG","CXCL9","CXCL10","STAT1","IRF1","GBP1"),
  Macrophage = c("CD68","CD163","MRC1","MSR1","CSF1R"),
  Treg = c("FOXP3","IL2RA","CTLA4","IKZF2"),
  Checkpoint_activity = c("PDCD1","CD274","CTLA4","LAG3","TIGIT","HAVCR2")
)

immune_sets <- lapply(immune_sets, function(x) {
  intersect(unique(clean_gene(x)), rownames(expr_log))
})
immune_sets <- immune_sets[sapply(immune_sets, length) >= 2]

immune_ssgsea <- safe_ssgsea(expr_log, immune_sets)
Immune_Activity <- colMeans(immune_ssgsea, na.rm = TRUE)

############################
# 7. Mutation burden
############################

read_mutation_flexible <- function(file) {
  if (!file.exists(file)) return(NULL)
  
  x <- data.table::fread(file, data.table = FALSE, check.names = FALSE)
  nms <- colnames(x)
  
  gene_col <- intersect(c("Hugo_Symbol", "Gene", "gene", "Gene_Symbol", "symbol"), nms)[1]
  sample_col <- intersect(c("Tumor_Sample_Barcode", "Sample", "sample", "Tumor_Sample", "patient"), nms)[1]
  
  if (!is.na(gene_col) && !is.na(sample_col)) {
    out <- x %>%
      dplyr::transmute(
        Gene = clean_gene(.data[[gene_col]]),
        Sample = gsub("\\.", "-", as.character(.data[[sample_col]]))
      ) %>%
      dplyr::filter(!is.na(Gene), Gene != "", !is.na(Sample), Sample != "")
    
    out$Patient <- patient_id(out$Sample)
    return(out)
  }
  
  gene_col <- colnames(x)[1]
  x[[gene_col]] <- clean_gene(x[[gene_col]])
  
  long <- x %>%
    tidyr::pivot_longer(
      cols = -all_of(gene_col),
      names_to = "Sample",
      values_to = "Mut"
    ) %>%
    dplyr::transmute(
      Gene = .data[[gene_col]],
      Sample = gsub("\\.", "-", Sample),
      Mut = suppressWarnings(as.numeric(Mut))
    ) %>%
    dplyr::filter(!is.na(Gene), Gene != "", is.finite(Mut), Mut != 0)
  
  long$Patient <- patient_id(long$Sample)
  long
}

mut_long <- read_mutation_flexible(MUT_FILE)

if (!is.null(mut_long) && nrow(mut_long) > 0) {
  mut_burden <- mut_long %>%
    dplyr::distinct(Patient, Gene) %>%
    dplyr::count(Patient, name = "Mutation_Burden")
  
  score_df <- score_df %>%
    dplyr::left_join(mut_burden, by = "Patient") %>%
    dplyr::mutate(
      Mutation_Burden = ifelse(is.na(Mutation_Burden), 0, Mutation_Burden)
    )
} else {
  score_df$Mutation_Burden <- NA_real_
}

############################
# 8. CNV burden
############################

read_cnv_wide <- function(file) {
  if (!file.exists(file)) return(NULL)
  
  x <- data.table::fread(file, data.table = FALSE, check.names = FALSE)
  gene_col <- colnames(x)[1]
  x[[gene_col]] <- clean_gene(x[[gene_col]])
  x <- x[!is.na(x[[gene_col]]) & x[[gene_col]] != "", ]
  x <- x[!duplicated(x[[gene_col]]), ]
  
  rownames(x) <- x[[gene_col]]
  x[[gene_col]] <- NULL
  
  mat <- as.matrix(x)
  storage.mode(mat) <- "numeric"
  rownames(mat) <- clean_gene(rownames(mat))
  colnames(mat) <- gsub("\\.", "-", colnames(mat))
  mat
}

cnv_mat <- read_cnv_wide(CNV_FILE)

if (!is.null(cnv_mat)) {
  ddr_genes <- unique(unlist(ddr_sets))
  common_ddr_cnv <- intersect(ddr_genes, rownames(cnv_mat))
  common_samples_cnv <- intersect(score_df$Sample, colnames(cnv_mat))
  
  if (length(common_ddr_cnv) >= 2 && length(common_samples_cnv) >= 10) {
    cnv_use <- cnv_mat[common_ddr_cnv, common_samples_cnv, drop = FALSE]
    
    cnv_burden <- data.frame(
      Sample = common_samples_cnv,
      CNV_Burden = colMeans(abs(cnv_use) > CNV_AMP_TH, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    score_df <- score_df %>%
      dplyr::left_join(cnv_burden, by = "Sample")
  } else {
    score_df$CNV_Burden <- NA_real_
  }
} else {
  score_df$CNV_Burden <- NA_real_
}

############################
# 9. Build landscape data
############################

land_df <- score_df %>%
  dplyr::mutate(
    DDR_Activity = as.numeric(DDR_Activity[Sample]),
    Immune_Activity = as.numeric(Immune_Activity[Sample])
  ) %>%
  dplyr::arrange(MO_DDRscore_raw) %>%
  dplyr::mutate(Rank = dplyr::row_number())

land_df$MO_DDRscore_group <- factor(land_df$MO_DDRscore_group, levels = c("Low", "High"))

save_csv(
  land_df,
  file.path(OUT_DIR, "Figure2A_MO_DDRscore_landscape_data.csv")
)

land_long <- land_df %>%
  dplyr::select(
    Sample,
    Rank,
    MO_DDRscore_group,
    MO_DDRscore_raw,
    DDR_Activity,
    Mutation_Burden,
    CNV_Burden,
    Immune_Activity
  ) %>%
  tidyr::pivot_longer(
    cols = c(MO_DDRscore_raw, DDR_Activity, Mutation_Burden, CNV_Burden, Immune_Activity),
    names_to = "Feature",
    values_to = "Value"
  ) %>%
  dplyr::group_by(Feature) %>%
  dplyr::mutate(Value_z = zscore_vec(Value)) %>%
  dplyr::ungroup()

land_long$Feature <- factor(
  land_long$Feature,
  levels = c(
    "MO_DDRscore_raw",
    "DDR_Activity",
    "Mutation_Burden",
    "CNV_Burden",
    "Immune_Activity"
  ),
  labels = c(
    "MO-DDRscore",
    "DDR activity",
    "Mutation burden",
    "DDR-CNV burden",
    "Immune activity"
  )
)

############################
# 10. Plot landscape
############################

p_score <- ggplot(
  land_df,
  aes(x = Rank, y = MO_DDRscore_raw, color = MO_DDRscore_group)
) +
  geom_point(size = 0.9, alpha = 0.85) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.7) +
  scale_color_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = NULL,
    y = "MO-DDRscore",
    color = "Group",
    title = "MO-DDRscore-defined DDR state landscape"
  )

p_group <- ggplot(
  land_df,
  aes(x = Rank, y = 1, fill = MO_DDRscore_group)
) +
  geom_tile() +
  scale_fill_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
  theme_void(base_size = 10) +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 5, 0, 5)
  ) +
  labs(x = NULL, y = NULL)

p_heat <- ggplot(
  land_long,
  aes(x = Rank, y = Feature, fill = Value_z)
) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#3B75AF",
    mid = "white",
    high = "#C84630",
    midpoint = 0,
    na.value = "grey90",
    name = "Z-score"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    legend.position = "right",
    plot.margin = margin(0, 5, 5, 5)
  ) +
  labs(x = "Patients ranked by MO-DDRscore", y = NULL)

p_final <- p_score / p_group / p_heat +
  patchwork::plot_layout(heights = c(2.2, 0.25, 1.6))

ggsave(
  filename = file.path(OUT_DIR, "Figure2A_MO_DDRscore_landscape.pdf"),
  plot = p_final,
  width = 10,
  height = 5.8,
  device = "pdf",
  useDingbats = FALSE
)

ggsave(
  filename = file.path(OUT_DIR, "Figure2A_MO_DDRscore_landscape.png"),
  plot = p_final,
  width = 10,
  height = 5.8,
  dpi = 600
)

cat("\nDone.\n")
cat("Output files:\n")
cat(file.path(OUT_DIR, "Figure2A_MO_DDRscore_landscape.pdf"), "\n")
cat(file.path(OUT_DIR, "Figure2A_MO_DDRscore_landscape.png"), "\n")
cat(file.path(OUT_DIR, "Figure2A_MO_DDRscore_landscape_data.csv"), "\n")