############################################################
## Fig3 multi-omics + GDSC/DepMap drug sensitivity pipeline
##
## Main Figure:
##   Fig3A1 DEG volcano
##   Fig3A2 Top DEG heatmap
##   Fig3B  Focused Hallmark/Reactome GSEA
##   Fig3C  Focused pathway ssGSEA activity
##   Fig3D1 DDR-only mutation landscape
##   Fig3D2 DDR mutation burden
##   Fig3E1 ESTIMATE tumor microenvironment
##   Fig3E2 CIBERSORT LM22 immune infiltration, if official files exist
##   Fig3F  GDSC/DepMap drug sensitivity
##   Fig3G  KM survival
##
## Supplementary:
##   checkpoint expression
##   clinical distribution
##   gene-level CNA altered fraction
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260513)

############################################################
## 0. Paths / parameters
############################################################

PROJECT_DIR <- "D:/R_workspace/评分/AD_DDR_project"
DATA_DIR    <- file.path(PROJECT_DIR, "00_data")
PROC_DIR    <- file.path(PROJECT_DIR, "01_processed")
DB_DIR      <- file.path(PROJECT_DIR, "05_database_tables")

FIG3_DIR     <- file.path(PROJECT_DIR, "fig", "Fig3_multiomics_GDSC")
TAB3_DIR     <- file.path(PROJECT_DIR, "table", "Fig3_multiomics_GDSC")
SUPP_DIR     <- file.path(PROJECT_DIR, "fig", "Supplementary_Fig3_multiomics_GDSC")
SUPP_TAB_DIR <- file.path(PROJECT_DIR, "table", "Supplementary_Fig3_multiomics_GDSC")

dir.create(FIG3_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB3_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DB_DIR, recursive = TRUE, showWarnings = FALSE)

FOCUS_CANCER <- "LUAD"
TCGA_DIR <- "D:/R/R_workspace/梁老师文件/TCGA"

TCGA_EXPR_FILE <- file.path(
  TCGA_DIR,
  "mRNA_exp_TPM_only_TCGA/mRNA_exp_TPM_only_TCGA",
  paste0("TCGA-", FOCUS_CANCER, ".gene_expression_TPM.tsv")
)

TCGA_MUT_FILE <- file.path(TCGA_DIR, "mutation", paste0(FOCUS_CANCER, ".txt"))
TCGA_CNV_FILE <- file.path(TCGA_DIR, "cnv", paste0(FOCUS_CANCER, ".txt"))

DDR_GENE_FILE <- file.path(DATA_DIR, "DDR_236_genes.csv")
SCORE_FILE    <- file.path(PROC_DIR, paste0(FOCUS_CANCER, "_MO_DDRscore.csv"))
CLIN_FILE     <- file.path(PROC_DIR, paste0(FOCUS_CANCER, "_clinical_processed.csv"))

## CIBERSORT official files
CIBERSORT_R_FILE <- file.path(DATA_DIR, "CIBERSORT.R")
LM22_FILE        <- file.path(DATA_DIR, "LM22.txt")

## DepMap / GDSC files
DEPMAP_EXPR_FILE <- file.path(DATA_DIR, "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
MODEL_MAP_FILE   <- file.path(DATA_DIR, "Model.csv")
GDSC_FILE        <- file.path(DATA_DIR, "GDSC2_fitted_dose_response_27Oct23.xlsx")

LOW_COL  <- "#4DBBD5"
HIGH_COL <- "#E64B35"
BLUE_COL <- "#2166AC"
RED_COL  <- "#B2182B"

DEG_LOGFC <- 1.20
DEG_FDR   <- 0.05

TOP_VOLCANO_LABEL <- 12
TOP_UP_N <- 25
TOP_DOWN_N <- 25
TOP_MUT_GENES <- 25
TOP_CHECKPOINT_GENES <- 12
TOP_CIBERSORT_CELLS <- 12

CNV_THRESHOLD <- 0.2

############################################################
## 1. Packages
############################################################

install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

cran_pkgs <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", "ggrepel",
  "survival", "survminer", "msigdbr", "pheatmap",
  "patchwork", "readxl"
)

bioc_pkgs <- c(
  "limma", "clusterProfiler", "maftools", "GSVA"
)

for (p in cran_pkgs) install_if_missing(p, bioc = FALSE)
for (p in bioc_pkgs) install_if_missing(p, bioc = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(survival)
  library(survminer)
  library(msigdbr)
  library(pheatmap)
  library(patchwork)
  library(readxl)
  library(limma)
  library(clusterProfiler)
  library(maftools)
  library(GSVA)
})

has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

############################################################
## 2. Helper functions
############################################################

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
  ggplot2::ggsave(
    filename = file,
    plot = p,
    width = w,
    height = h,
    device = "pdf",
    useDingbats = FALSE
  )
}

open_pdf <- function(file, w = 7, h = 5) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(file, width = w, height = h, useDingbats = FALSE)
}

find_col <- function(df, patterns) {
  cn <- colnames(df)
  for (pa in patterns) {
    hit <- grep(pa, cn, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

pick_col <- function(df, candidates) {
  cn <- colnames(df)
  for (x in candidates) {
    hit <- cn[tolower(cn) == tolower(x)]
    if (length(hit) > 0) return(hit[1])
  }
  for (x in candidates) {
    hit <- grep(x, cn, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub(" ///.*$| //.*$|;.*$|,.*$|\\s*\\([^\\)]*\\)$", "", x)
  x[x %in% c("", "NA", "---", "NULL", "N/A", "NAN")] <- NA
  x
}

clean_depmap_gene <- function(x) {
  x <- as.character(x)
  x <- gsub("\\s*\\(.*?\\)$", "", x)
  clean_gene(x)
}

clean_na <- function(x) {
  x <- as.character(x)
  x[x %in% c(
    "", "NA", "NaN", "[Not Available]", "[Not Applicable]",
    "not reported", "Not Reported", "null", "NULL", "--"
  )] <- NA
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

z_vec <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

z_rows <- function(m) {
  m <- as.matrix(m)
  storage.mode(m) <- "numeric"
  mu <- rowMeans(m, na.rm = TRUE)
  sdv <- apply(m, 1, sd, na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- 1
  out <- sweep(sweep(m, 1, mu, "-"), 1, sdv, "/")
  out[!is.finite(out)] <- 0
  out
}

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_rect(fill = "grey95", color = "grey60"),
      strip.text = element_text(face = "bold", color = "black"),
      plot.title = element_text(face = "bold", hjust = 0)
    )
}

format_p <- function(p) {
  ifelse(
    is.na(p),
    "P = NA",
    ifelse(p < 0.001, "P < 0.001", paste0("P = ", signif(p, 3)))
  )
}

sig_label <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "NS"
  )
}

pretty_term <- function(x) {
  x <- gsub("^HALLMARK_", "Hallmark: ", x)
  x <- gsub("^REACTOME_", "Reactome: ", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  x <- gsub("Dna", "DNA", x)
  x <- gsub("Rna", "RNA", x)
  x <- gsub("G2 M", "G2/M", x)
  x <- gsub("G2m", "G2/M", x)
  x <- gsub("E2f", "E2F", x)
  x <- gsub("Atr", "ATR", x)
  x <- gsub("P53", "p53", x)
  x <- gsub("Nhej", "NHEJ", x)
  x
}

safe_cor <- function(x, y) {
  out <- tryCatch(
    suppressWarnings(cor.test(x, y, method = "spearman")),
    error = function(e) NULL
  )
  if (is.null(out)) {
    return(data.frame(Rho = NA_real_, CorP = NA_real_))
  }
  data.frame(
    Rho = unname(out$estimate),
    CorP = out$p.value
  )
}

############################################################
## 3. Read input data
############################################################

read_expr <- function(file) {
  x <- safe_fread(file)
  
  gc <- find_col(
    x,
    c("^Gene$", "gene", "symbol", "Gene Symbol", "gene_name", "Hugo_Symbol")
  )
  if (is.na(gc)) gc <- colnames(x)[1]
  
  genes <- clean_gene(x[[gc]])
  
  matdf <- x[, setdiff(colnames(x), gc), drop = FALSE]
  matdf$Gene <- genes
  matdf <- matdf[!is.na(matdf$Gene) & matdf$Gene != "", ]
  
  matdf <- matdf %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::everything(),
        ~ suppressWarnings(mean(as.numeric(.x), na.rm = TRUE))
      ),
      .groups = "drop"
    )
  
  rn <- matdf$Gene
  matdf$Gene <- NULL
  
  mat <- as.matrix(as.data.frame(matdf, check.names = FALSE))
  rownames(mat) <- rn
  storage.mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0
  
  message("Expression matrix: ", nrow(mat), " genes x ", ncol(mat), " samples")
  mat
}

read_ddr_genes <- function(file) {
  x <- safe_fread(file)
  gc <- find_col(x, c("^Gene$", "gene", "symbol", "Gene Symbol", "Hugo_Symbol"))
  if (is.na(gc)) gc <- colnames(x)[1]
  genes <- unique(na.omit(clean_gene(x[[gc]])))
  genes <- genes[genes != ""]
  message("DDR genes loaded: ", length(genes))
  genes
}

read_score <- function(file) {
  x <- safe_fread(file)
  
  if (!all(c("Patient", "MO_DDRscore_raw", "MO_DDRscore_group") %in% colnames(x))) {
    stop("Score file must contain Patient, MO_DDRscore_raw and MO_DDRscore_group.")
  }
  
  if ("SampleClass" %in% colnames(x)) {
    x <- x %>% dplyr::filter(SampleClass == "Tumor")
  }
  
  x$Patient <- patient_id(x$Patient)
  x$ScoreSample <- if ("Sample" %in% colnames(x)) as.character(x$Sample) else x$Patient
  
  x %>%
    dplyr::filter(is.finite(MO_DDRscore_raw), !is.na(MO_DDRscore_group)) %>%
    dplyr::mutate(
      MO_DDRscore_group = factor(as.character(MO_DDRscore_group), levels = c("Low", "High"))
    ) %>%
    dplyr::select(
      Patient,
      ScoreSample,
      MO_DDRscore_raw,
      MO_DDRscore_group,
      dplyr::any_of(c("MO_DDRscore_0_100"))
    ) %>%
    dplyr::arrange(Patient) %>%
    dplyr::distinct(Patient, .keep_all = TRUE)
}

read_clinical_processed <- function(file) {
  x <- safe_fread(file)
  
  req <- c("Patient", "time", "status")
  if (!all(req %in% colnames(x))) {
    stop("Clinical processed file must contain Patient, time and status.")
  }
  
  if (!"age" %in% colnames(x)) x$age <- NA_real_
  if (!"gender" %in% colnames(x)) x$gender <- NA_character_
  if (!"stage" %in% colnames(x)) x$stage <- NA_character_
  
  out <- x %>%
    dplyr::mutate(
      Patient = patient_id(Patient),
      time = suppressWarnings(as.numeric(time)),
      status = suppressWarnings(as.numeric(status)),
      age = suppressWarnings(as.numeric(age)),
      gender = clean_na(gender),
      stage = clean_na(stage)
    ) %>%
    dplyr::filter(is.finite(time), time > 0, status %in% c(0, 1)) %>%
    dplyr::arrange(Patient) %>%
    dplyr::distinct(Patient, .keep_all = TRUE)
  
  if (median(out$age, na.rm = TRUE) > 150) {
    out$age <- out$age / 365.25
  }
  
  out
}

read_mut <- function(file) {
  if (!file.exists(file)) {
    message("Mutation file not found: ", file)
    return(NULL)
  }
  
  x <- safe_fread(file)
  
  gc <- find_col(x, c("^Hugo_Symbol$", "^Gene$", "gene", "symbol"))
  sc <- find_col(x, c("^Tumor_Sample_Barcode$", "sample", "barcode", "Tumor"))
  
  if (is.na(gc) || is.na(sc)) {
    message("Mutation file does not contain gene/sample columns.")
    return(NULL)
  }
  
  x$Hugo_Symbol <- clean_gene(x[[gc]])
  x$Tumor_Sample_Barcode <- as.character(x[[sc]])
  x$Patient <- patient_id(x$Tumor_Sample_Barcode)
  
  if (!"Variant_Classification" %in% colnames(x)) x$Variant_Classification <- "Mutation"
  if (!"Variant_Type" %in% colnames(x)) x$Variant_Type <- "SNP"
  if (!"Chromosome" %in% colnames(x)) x$Chromosome <- "1"
  if (!"Start_Position" %in% colnames(x)) x$Start_Position <- 1
  if (!"End_Position" %in% colnames(x)) x$End_Position <- x$Start_Position
  if (!"Reference_Allele" %in% colnames(x)) x$Reference_Allele <- "N"
  if (!"Tumor_Seq_Allele2" %in% colnames(x)) x$Tumor_Seq_Allele2 <- "N"
  
  message("Mutation records loaded: ", nrow(x))
  x
}

read_cnv <- function(file) {
  if (!file.exists(file)) {
    message("CNV file not found: ", file)
    return(NULL)
  }
  
  x <- safe_fread(file)
  
  gc <- find_col(x, c("Gene Symbol", "^Gene$", "gene", "symbol", "Hugo_Symbol"))
  if (is.na(gc)) gc <- colnames(x)[1]
  
  genes <- clean_gene(x[[gc]])
  
  matdf <- x[, setdiff(colnames(x), gc), drop = FALSE]
  matdf$Gene <- genes
  matdf <- matdf[!is.na(matdf$Gene) & matdf$Gene != "", ]
  
  matdf <- matdf %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::everything(),
        ~ suppressWarnings(mean(as.numeric(.x), na.rm = TRUE))
      ),
      .groups = "drop"
    )
  
  rn <- matdf$Gene
  matdf$Gene <- NULL
  
  mat <- as.matrix(as.data.frame(matdf, check.names = FALSE))
  rownames(mat) <- rn
  storage.mode(mat) <- "numeric"
  
  message("CNV matrix loaded: ", nrow(mat), " genes x ", ncol(mat), " samples")
  mat
}

tcga_expr <- read_expr(TCGA_EXPR_FILE)
DDR_GENES <- read_ddr_genes(DDR_GENE_FILE)
score_df <- read_score(SCORE_FILE)
clin_df <- read_clinical_processed(CLIN_FILE)
mut <- read_mut(TCGA_MUT_FILE)
cnv_mat <- read_cnv(TCGA_CNV_FILE)

tumor_samples <- colnames(tcga_expr)[sample_class(colnames(tcga_expr)) == "Tumor"]

sample_patient <- data.frame(
  Sample = tumor_samples,
  Patient = patient_id(tumor_samples),
  stringsAsFactors = FALSE
)

sample_map <- sample_patient %>%
  dplyr::inner_join(score_df, by = "Patient") %>%
  dplyr::filter(Sample %in% colnames(tcga_expr)) %>%
  dplyr::arrange(MO_DDRscore_group, MO_DDRscore_raw) %>%
  dplyr::distinct(Patient, .keep_all = TRUE)

if (nrow(sample_map) < 50) {
  stop("Too few matched tumor samples: ", nrow(sample_map))
}

expr_use <- tcga_expr[, sample_map$Sample, drop = FALSE]
log_expr <- log2(expr_use + 1)

score_clin <- sample_map %>%
  dplyr::inner_join(clin_df, by = "Patient")

save_csv(sample_map, file.path(TAB3_DIR, "Fig3_sample_map_used.csv"))
save_csv(score_clin, file.path(TAB3_DIR, "Fig3_score_clinical_used.csv"))

cat("Matched tumor samples:", nrow(sample_map), "\n")
cat("Clinical matched patients:", nrow(score_clin), "\n")
print(table(sample_map$MO_DDRscore_group, useNA = "ifany"))

sample_map <- sample_map %>%
  dplyr::filter(!is.na(MO_DDRscore_group)) %>%
  dplyr::mutate(
    MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High"))
  ) %>%
  dplyr::arrange(MO_DDRscore_group, MO_DDRscore_raw)

fig_samples <- sample_map$Sample
expr_fig <- log_expr[, fig_samples, drop = FALSE]
expr_tpm_fig <- expr_use[, fig_samples, drop = FALSE]

group <- factor(sample_map$MO_DDRscore_group, levels = c("Low", "High"))

if (ncol(expr_fig) != length(group)) {
  stop("Sample mismatch in Fig3.")
}

############################################################
## 4. Fig3A1: DEG volcano
############################################################

cat("\nRunning Fig3A1 DEG volcano...\n")

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- limma::lmFit(expr_fig, design)
cont <- limma::makeContrasts(High - Low, levels = design)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont))

deg <- limma::topTable(fit2, number = Inf, adjust.method = "BH")
deg$Gene <- rownames(deg)
deg$FDR <- deg$adj.P.Val
deg$tstat <- deg$t

deg <- deg %>%
  dplyr::mutate(
    Direction = dplyr::case_when(
      FDR < DEG_FDR & logFC > DEG_LOGFC ~ "Up in High",
      FDR < DEG_FDR & logFC < -DEG_LOGFC ~ "Down in High",
      TRUE ~ "NS"
    ),
    negLogFDR = -log10(FDR + 1e-300)
  ) %>%
  dplyr::arrange(FDR, dplyr::desc(abs(logFC)))

save_csv(deg, file.path(TAB3_DIR, "Fig3A_DEG_high_vs_low.csv"))
save_csv(deg, file.path(DB_DIR, "DEG_table.csv"))

label_genes <- deg %>%
  dplyr::filter(Direction != "NS") %>%
  dplyr::arrange(FDR, dplyr::desc(abs(logFC))) %>%
  dplyr::slice_head(n = TOP_VOLCANO_LABEL)

p_volcano <- ggplot(deg, aes(x = logFC, y = negLogFDR)) +
  geom_point(aes(color = Direction), size = 0.8, alpha = 0.75) +
  geom_vline(xintercept = c(-DEG_LOGFC, DEG_LOGFC), linetype = 2, color = "grey50") +
  geom_hline(yintercept = -log10(DEG_FDR), linetype = 2, color = "grey50") +
  ggrepel::geom_text_repel(
    data = label_genes,
    aes(label = Gene),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.35,
    min.segment.length = 0
  ) +
  scale_color_manual(
    values = c(
      "Up in High" = HIGH_COL,
      "Down in High" = LOW_COL,
      "NS" = "grey75"
    )
  ) +
  theme_pub(11.5) +
  labs(
    x = "log2 fold change: High vs Low MO-DDRscore",
    y = "-log10(FDR)",
    color = NULL,
    title = "Differential expression associated with MO-DDRscore"
  )

safe_ggsave(
  file.path(FIG3_DIR, "Fig3A1_DEG_volcano.pdf"),
  p_volcano,
  6.2,
  5.2
)

############################################################
## 5. Fig3A2: top DEG heatmap
############################################################

cat("\nRunning Fig3A2 top DEG heatmap...\n")

top_up <- deg %>%
  dplyr::filter(Direction == "Up in High") %>%
  dplyr::arrange(FDR, dplyr::desc(logFC)) %>%
  dplyr::slice_head(n = TOP_UP_N) %>%
  dplyr::pull(Gene)

top_down <- deg %>%
  dplyr::filter(Direction == "Down in High") %>%
  dplyr::arrange(FDR, logFC) %>%
  dplyr::slice_head(n = TOP_DOWN_N) %>%
  dplyr::pull(Gene)

top_genes <- unique(c(top_up, top_down))
top_genes <- intersect(top_genes, rownames(expr_fig))

if (length(top_genes) < 10) {
  top_genes <- deg %>%
    dplyr::arrange(FDR, dplyr::desc(abs(logFC))) %>%
    dplyr::slice_head(n = 50) %>%
    dplyr::pull(Gene) %>%
    intersect(rownames(expr_fig))
}

save_csv(
  data.frame(Gene = top_genes),
  file.path(TAB3_DIR, "Fig3A_heatmap_genes_used.csv")
)

hm_mat <- z_rows(expr_fig[top_genes, fig_samples, drop = FALSE])
hm_mat[hm_mat > 2.5] <- 2.5
hm_mat[hm_mat < -2.5] <- -2.5

anno_col <- data.frame(
  MO_DDRscore_group = sample_map$MO_DDRscore_group
)

rownames(anno_col) <- sample_map$Sample

anno_colors <- list(
  MO_DDRscore_group = c(
    Low = LOW_COL,
    High = HIGH_COL
  )
)

gap_pos <- sum(sample_map$MO_DDRscore_group == "Low")
if (gap_pos <= 0 || gap_pos >= nrow(sample_map)) gap_pos <- NULL

open_pdf(
  file.path(FIG3_DIR, "Fig3A2_top_DEG_heatmap.pdf"),
  8.8,
  7.6
)

pheatmap::pheatmap(
  hm_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = gap_pos,
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_col = anno_col,
  annotation_colors = anno_colors,
  color = colorRampPalette(c(BLUE_COL, "white", RED_COL))(100),
  fontsize_row = 7,
  border_color = NA,
  main = "Top differentially expressed genes"
)

dev.off()

############################################################
## 6. Official focused MSigDB / Reactome gene sets
############################################################

cat("\nPreparing official focused MSigDB/Reactome gene sets...\n")

focused_sets <- c(
  "HALLMARK_DNA_REPAIR",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "REACTOME_DNA_REPAIR",
  "REACTOME_DNA_REPLICATION",
  "REACTOME_CELL_CYCLE_CHECKPOINTS",
  "REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT",
  "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS",
  "REACTOME_NUCLEOTIDE_EXCISION_REPAIR",
  "REACTOME_MISMATCH_REPAIR",
  "REACTOME_HOMOLOGOUS_RECOMBINATION_REPAIR_OF_REPLICATION_INDEPENDENT_DOUBLE_STRAND_BREAKS",
  "REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ",
  "REACTOME_FANCONI_ANEMIA_PATHWAY",
  "REACTOME_BASE_EXCISION_REPAIR",
  "REACTOME_TRANSLESION_SYNTHESIS_BY_Y_FAMILY_DNA_POLYMERASES_BYPASSES_LESIONS_ON_DNA_TEMPLATE"
)

msig_all <- msigdbr::msigdbr(species = "Homo sapiens")

if (!"gs_name" %in% colnames(msig_all)) {
  stop("msigdbr output does not contain column 'gs_name'. Please update msigdbr.")
}

term2gene <- msig_all %>%
  dplyr::filter(gs_name %in% focused_sets) %>%
  dplyr::transmute(
    gs_name = gs_name,
    gene_symbol = clean_gene(gene_symbol)
  ) %>%
  dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
  dplyr::distinct()

if (nrow(term2gene) == 0) {
  stop("No focused MSigDB/Reactome gene sets matched. Please check msigdbr version.")
}

save_csv(
  term2gene,
  file.path(TAB3_DIR, "Fig3_official_MSigDB_Reactome_gene_sets_used.csv")
)

############################################################
## 7. Fig3B: focused GSEA
############################################################

cat("\nRunning Fig3B focused GSEA...\n")

rank_df <- deg %>%
  dplyr::mutate(
    Gene = clean_gene(Gene),
    rank_stat = as.numeric(tstat)
  ) %>%
  dplyr::filter(!is.na(Gene), Gene != "", is.finite(rank_stat)) %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(rank_stat = mean(rank_stat, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(rank_stat))

gene_list <- rank_df$rank_stat
names(gene_list) <- rank_df$Gene
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_res <- clusterProfiler::GSEA(
  geneList = gene_list,
  TERM2GENE = term2gene,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  seed = TRUE
)

gsea_df <- as.data.frame(gsea_res@result)

save_csv(
  gsea_df,
  file.path(TAB3_DIR, "Fig3B_focused_GSEA_high_vs_low.csv")
)

if (nrow(gsea_df) > 0) {
  
  gsea_plot <- gsea_df %>%
    dplyr::mutate(
      Term = if ("Description" %in% colnames(.)) Description else ID,
      Label = pretty_term(Term),
      Direction = ifelse(NES >= 0, "High MO-DDRscore", "Low MO-DDRscore"),
      negLogFDR = -log10(p.adjust + 1e-300)
    ) %>%
    dplyr::filter(is.finite(NES), !is.na(p.adjust)) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = min(16, nrow(.))) %>%
    dplyr::mutate(Label = factor(Label, levels = rev(Label)))
  
  p_gsea <- ggplot(gsea_plot, aes(x = NES, y = Label)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
    geom_point(aes(size = negLogFDR, color = Direction), alpha = 0.92) +
    scale_color_manual(
      values = c(
        "High MO-DDRscore" = HIGH_COL,
        "Low MO-DDRscore" = LOW_COL
      )
    ) +
    theme_pub(11.5) +
    theme(panel.grid.major.y = element_blank()) +
    labs(
      x = "Normalized enrichment score",
      y = NULL,
      size = "-log10(FDR)",
      color = NULL,
      title = "Focused MSigDB/Reactome GSEA"
    )
  
  safe_ggsave(
    file.path(FIG3_DIR, "Fig3B_focused_MSigDB_Reactome_GSEA.pdf"),
    p_gsea,
    8.8,
    6.0
  )
}

############################################################
## 8. Fig3C: focused ssGSEA pathway activity
############################################################

cat("\nRunning Fig3C focused ssGSEA pathway activity...\n")

gene_sets <- split(term2gene$gene_symbol, term2gene$gs_name)
gene_sets <- lapply(gene_sets, function(x) intersect(unique(x), rownames(expr_fig)))
gene_sets <- gene_sets[sapply(gene_sets, length) >= 10]

run_ssgsea <- function(expr_mat, gene_sets) {
  expr_mat <- as.matrix(expr_mat)
  
  if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    param <- GSVA::ssgseaParam(
      exprData = expr_mat,
      geneSets = gene_sets,
      normalize = TRUE
    )
    return(GSVA::gsva(param))
  }
  
  GSVA::gsva(
    expr = expr_mat,
    gset.idx.list = gene_sets,
    method = "ssgsea",
    kcdf = "Gaussian",
    abs.ranking = FALSE,
    ssgsea.norm = TRUE,
    verbose = FALSE
  )
}

ssgsea_mat <- tryCatch(
  run_ssgsea(expr_fig, gene_sets),
  error = function(e) {
    message("ssGSEA failed: ", e$message)
    NULL
  }
)

if (!is.null(ssgsea_mat)) {
  
  ssgsea_df <- as.data.frame(t(ssgsea_mat), check.names = FALSE)
  ssgsea_df$Sample <- rownames(ssgsea_df)
  
  ssgsea_df <- ssgsea_df %>%
    dplyr::left_join(
      sample_map %>%
        dplyr::select(Sample, Patient, MO_DDRscore_group, MO_DDRscore_raw),
      by = "Sample"
    )
  
  save_csv(
    ssgsea_df,
    file.path(TAB3_DIR, "Fig3C_focused_pathway_ssGSEA_scores.csv")
  )
  
  pathway_long <- ssgsea_df %>%
    tidyr::pivot_longer(
      cols = all_of(names(gene_sets)),
      names_to = "Pathway",
      values_to = "Score"
    ) %>%
    dplyr::filter(is.finite(Score), !is.na(MO_DDRscore_group))
  
  pathway_stat <- pathway_long %>%
    dplyr::group_by(Pathway) %>%
    dplyr::summarise(
      Median_Low = median(Score[MO_DDRscore_group == "Low"], na.rm = TRUE),
      Median_High = median(Score[MO_DDRscore_group == "High"], na.rm = TRUE),
      Diff = Median_High - Median_Low,
      P = tryCatch(wilcox.test(Score ~ MO_DDRscore_group)$p.value, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      FDR = p.adjust(P, "BH"),
      Label = sig_label(FDR),
      Pretty = pretty_term(Pathway)
    ) %>%
    dplyr::arrange(FDR, dplyr::desc(abs(Diff)))
  
  save_csv(
    pathway_stat,
    file.path(TAB3_DIR, "Fig3C_focused_pathway_ssGSEA_stats.csv")
  )
  
  show_pathways <- pathway_stat %>%
    dplyr::slice_head(n = min(14, nrow(.))) %>%
    dplyr::pull(Pathway)
  
  order_pretty <- pathway_stat$Pretty[pathway_stat$Pathway %in% show_pathways]
  
  label_df <- pathway_long %>%
    dplyr::filter(Pathway %in% show_pathways) %>%
    dplyr::group_by(Pathway) %>%
    dplyr::summarise(y_pos = max(Score, na.rm = TRUE) * 1.05, .groups = "drop") %>%
    dplyr::left_join(pathway_stat[, c("Pathway", "Label", "Pretty")], by = "Pathway") %>%
    dplyr::mutate(
      Pretty = factor(Pretty, levels = rev(order_pretty)),
      x_pos = 1.5
    )
  
  p_pathway <- pathway_long %>%
    dplyr::filter(Pathway %in% show_pathways) %>%
    dplyr::left_join(pathway_stat[, c("Pathway", "Pretty")], by = "Pathway") %>%
    dplyr::mutate(
      Pretty = factor(Pretty, levels = rev(order_pretty))
    ) %>%
    ggplot(aes(x = MO_DDRscore_group, y = Score, fill = MO_DDRscore_group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.62, linewidth = 0.32) +
    geom_jitter(width = 0.12, size = 0.42, alpha = 0.25) +
    geom_text(
      data = label_df,
      aes(x = x_pos, y = y_pos, label = Label),
      inherit.aes = FALSE,
      size = 3,
      fontface = "bold"
    ) +
    facet_wrap(~Pretty, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c(Low = LOW_COL, High = HIGH_COL)) +
    theme_pub(10.5) +
    theme(legend.position = "none", strip.text = element_text(face = "bold")) +
    labs(
      x = NULL,
      y = "ssGSEA score",
      title = "Official MSigDB/Reactome pathway activity"
    )
  
  safe_ggsave(
    file.path(FIG3_DIR, "Fig3C_focused_pathway_ssGSEA_activity.pdf"),
    p_pathway,
    9.5,
    8
  )
}

############################################################
## 9. Fig3D1: DDR-only mutation landscape
############################################################

cat("\nRunning Fig3D1 DDR-only mutation landscape...\n")

if (!is.null(mut)) {
  
  pg <- sample_map %>%
    dplyr::select(Patient, Sample, MO_DDRscore_raw, MO_DDRscore_group) %>%
    dplyr::distinct(Patient, .keep_all = TRUE) %>%
    dplyr::arrange(MO_DDRscore_group, MO_DDRscore_raw)
  
  mut_ddr <- mut %>%
    dplyr::filter(
      Patient %in% pg$Patient,
      !is.na(Hugo_Symbol),
      Hugo_Symbol %in% DDR_GENES
    )
  
  save_csv(
    mut_ddr,
    file.path(TAB3_DIR, "Fig3D_DDR_only_mutation_maf_input.csv")
  )
  
  if (nrow(mut_ddr) > 0 && has_pkg("maftools")) {
    
    maf_tmp <- mut_ddr
    maf_tmp$Tumor_Sample_Barcode <- maf_tmp$Patient
    
    clinical_anno <- pg %>%
      dplyr::select(Tumor_Sample_Barcode = Patient, MO_DDRscore_group) %>%
      dplyr::mutate(MO_DDRscore_group = as.character(MO_DDRscore_group))
    
    maf_obj <- tryCatch(
      maftools::read.maf(
        maf = maf_tmp,
        clinicalData = clinical_anno,
        verbose = FALSE
      ),
      error = function(e) {
        message("maftools::read.maf failed: ", e$message)
        NULL
      }
    )
    
    if (!is.null(maf_obj)) {
      open_pdf(file.path(FIG3_DIR, "Fig3D1_DDR_only_mutation_oncoplot.pdf"), 10.8, 7.2)
      
      try(
        maftools::oncoplot(
          maf = maf_obj,
          top = TOP_MUT_GENES,
          clinicalFeatures = "MO_DDRscore_group",
          sortByAnnotation = TRUE,
          annotationColor = list(
            MO_DDRscore_group = c(Low = LOW_COL, High = HIGH_COL)
          ),
          removeNonMutated = TRUE,
          fontSize = 0.8
        ),
        silent = TRUE
      )
      
      dev.off()
    }
  }
  
  ############################################################
  ## Fig3D2: DDR mutation burden
  ############################################################
  
  nonsyn_class <- c(
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
    "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
    "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins"
  )
  
  mut_ddr_nonsyn <- mut_ddr %>%
    dplyr::filter(
      Variant_Classification %in% nonsyn_class
    )
  
  ddr_mut_burden <- pg %>%
    dplyr::select(Patient, MO_DDRscore_group, MO_DDRscore_raw) %>%
    dplyr::left_join(
      mut_ddr_nonsyn %>%
        dplyr::group_by(Patient) %>%
        dplyr::summarise(
          DDR_Mutated_Gene_N = length(unique(Hugo_Symbol)),
          DDR_Mutation_N = dplyr::n(),
          .groups = "drop"
        ),
      by = "Patient"
    ) %>%
    dplyr::mutate(
      DDR_Mutated_Gene_N = ifelse(is.na(DDR_Mutated_Gene_N), 0, DDR_Mutated_Gene_N),
      DDR_Mutation_N = ifelse(is.na(DDR_Mutation_N), 0, DDR_Mutation_N),
      MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High"))
    )
  
  save_csv(
    ddr_mut_burden,
    file.path(TAB3_DIR, "Fig3D2_DDR_mutation_burden_table.csv")
  )
  
  burden_long <- ddr_mut_burden %>%
    tidyr::pivot_longer(
      cols = c(DDR_Mutated_Gene_N, DDR_Mutation_N),
      names_to = "BurdenType",
      values_to = "Burden"
    )
  
  burden_stat <- burden_long %>%
    dplyr::group_by(BurdenType) %>%
    dplyr::summarise(
      Median_Low = median(Burden[MO_DDRscore_group == "Low"], na.rm = TRUE),
      Median_High = median(Burden[MO_DDRscore_group == "High"], na.rm = TRUE),
      Diff = Median_High - Median_Low,
      P = tryCatch(wilcox.test(Burden ~ MO_DDRscore_group)$p.value, error = function(e) NA_real_),
      y_pos = max(Burden, na.rm = TRUE) * 1.08 + 0.2,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      FDR = p.adjust(P, "BH"),
      Label = sig_label(FDR),
      x_pos = 1.5
    )
  
  save_csv(
    burden_stat,
    file.path(TAB3_DIR, "Fig3D2_DDR_mutation_burden_stats.csv")
  )
  
  burden_long$BurdenType <- factor(
    burden_long$BurdenType,
    levels = c("DDR_Mutated_Gene_N", "DDR_Mutation_N"),
    labels = c("Number of mutated DDR genes", "Number of DDR mutations")
  )
  
  burden_stat$BurdenType <- factor(
    burden_stat$BurdenType,
    levels = c("DDR_Mutated_Gene_N", "DDR_Mutation_N"),
    labels = c("Number of mutated DDR genes", "Number of DDR mutations")
  )
  
  p_burden <- ggplot(
    burden_long,
    aes(x = MO_DDRscore_group, y = Burden, fill = MO_DDRscore_group)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.62, linewidth = 0.35) +
    geom_jitter(width = 0.14, size = 0.55, alpha = 0.35) +
    geom_text(
      data = burden_stat,
      aes(x = x_pos, y = y_pos, label = Label),
      inherit.aes = FALSE,
      size = 4,
      fontface = "bold"
    ) +
    facet_wrap(~BurdenType, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c(Low = LOW_COL, High = HIGH_COL)) +
    theme_pub(11.5) +
    theme(legend.position = "none", strip.text = element_text(face = "bold")) +
    labs(
      x = NULL,
      y = "Count per patient",
      title = "DDR mutation burden associated with MO-DDRscore"
    )
  
  safe_ggsave(
    file.path(FIG3_DIR, "Fig3D2_DDR_mutation_burden.pdf"),
    p_burden,
    7.5,
    4.2
  )
}

############################################################
## 10. Fig3E1: ESTIMATE tumor microenvironment
############################################################

cat("\nRunning Fig3E1 ESTIMATE tumor microenvironment...\n")

get_existing_estimate <- function() {
  candidates <- c(
    file.path(PROJECT_DIR, "table", "Fig2_standard_multiomics", "Fig2_ESTIMATE_scores.csv"),
    file.path(PROJECT_DIR, "table", "Fig3_standard_multiomics", "Fig3D_ESTIMATE_scores.csv"),
    file.path(TAB3_DIR, "Fig3E_ESTIMATE_scores.csv")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

run_estimate_simple <- function(expr_log, sample_map) {
  
  if (!has_pkg("estimate")) {
    stop("Package estimate is not installed and no existing ESTIMATE score file was found.")
  }
  
  ESTIMATE_DIR <- file.path(TAB3_DIR, "ESTIMATE_work")
  dir.create(ESTIMATE_DIR, recursive = TRUE, showWarnings = FALSE)
  
  input_txt <- file.path(ESTIMATE_DIR, "estimate_input.txt")
  common_gct <- file.path(ESTIMATE_DIR, "estimate_common.gct")
  score_gct <- file.path(ESTIMATE_DIR, "estimate_score.gct")
  
  expr_est <- as.matrix(expr_log)
  storage.mode(expr_est) <- "numeric"
  rownames(expr_est) <- clean_gene(rownames(expr_est))
  expr_est <- expr_est[!is.na(rownames(expr_est)) & rownames(expr_est) != "", , drop = FALSE]
  expr_est <- expr_est[!duplicated(rownames(expr_est)), , drop = FALSE]
  expr_est[!is.finite(expr_est)] <- 0
  
  input_df <- data.frame(GeneSymbol = rownames(expr_est), expr_est, check.names = FALSE)
  data.table::fwrite(input_df, input_txt, sep = "\t", quote = FALSE, na = "0")
  
  estimate::filterCommonGenes(input.f = input_txt, output.f = common_gct, id = "GeneSymbol")
  estimate::estimateScore(input.ds = common_gct, output.ds = score_gct, platform = "illumina")
  
  est_raw <- data.table::fread(score_gct, skip = 2, data.table = FALSE, check.names = FALSE)
  rownames(est_raw) <- est_raw$NAME
  
  est_mat <- est_raw[, -c(1, 2), drop = FALSE]
  est_df <- as.data.frame(t(est_mat), check.names = FALSE)
  est_df$Sample <- rownames(est_df)
  est_df$Sample <- gsub("\\.", "-", est_df$Sample)
  
  est_df$ImmuneScore <- suppressWarnings(as.numeric(est_df$ImmuneScore))
  est_df$StromalScore <- suppressWarnings(as.numeric(est_df$StromalScore))
  est_df$ESTIMATEScore <- suppressWarnings(as.numeric(est_df$ESTIMATEScore))
  est_df$TumorPurity <- cos(0.6049872018 + 0.0001467884 * est_df$ESTIMATEScore)
  
  est_df %>%
    dplyr::select(Sample, ImmuneScore, StromalScore, ESTIMATEScore, TumorPurity) %>%
    dplyr::left_join(
      sample_map %>%
        dplyr::select(Sample, Patient, MO_DDRscore_raw, MO_DDRscore_group),
      by = "Sample"
    )
}

estimate_df <- NULL
existing_estimate <- get_existing_estimate()

############################################################
## Robustly read existing ESTIMATE file and merge MO-DDRscore group
############################################################

standardize_existing_estimate <- function(est_file, sample_map) {
  
  x <- data.table::fread(
    est_file,
    data.table = FALSE,
    check.names = FALSE
  )
  
  colnames(x) <- gsub("\\s+", "", colnames(x))
  
  score_need <- c("ImmuneScore", "StromalScore", "ESTIMATEScore")
  score_hit <- intersect(score_need, colnames(x))
  
  if (length(score_hit) < 3) {
    stop(
      "Existing ESTIMATE file does not contain ImmuneScore/StromalScore/ESTIMATEScore: ",
      est_file
    )
  }
  
  sample_col <- pick_col(
    x,
    c(
      "^Sample$",
      "^sample$",
      "sampleID",
      "SampleID",
      "Tumor_Sample_Barcode",
      "ID",
      "^V1$"
    )
  )
  
  if (is.na(sample_col)) {
    non_score_cols <- setdiff(
      colnames(x),
      c("ImmuneScore", "StromalScore", "ESTIMATEScore", "TumorPurity")
    )
    
    tcga_like <- non_score_cols[
      sapply(non_score_cols, function(cc) {
        any(grepl("^TCGA", as.character(x[[cc]])))
      })
    ]
    
    if (length(tcga_like) > 0) {
      sample_col <- tcga_like[1]
    }
  }
  
  if (!is.na(sample_col)) {
    x$Sample <- gsub("\\.", "-", as.character(x[[sample_col]]))
    x$Patient <- patient_id(x$Sample)
  } else if ("Patient" %in% colnames(x)) {
    x$Patient <- patient_id(x$Patient)
    x$Sample <- x$Patient
  } else {
    stop("Cannot identify Sample/Patient column in existing ESTIMATE file: ", est_file)
  }
  
  x$ImmuneScore <- suppressWarnings(as.numeric(x$ImmuneScore))
  x$StromalScore <- suppressWarnings(as.numeric(x$StromalScore))
  x$ESTIMATEScore <- suppressWarnings(as.numeric(x$ESTIMATEScore))
  
  if (!"TumorPurity" %in% colnames(x)) {
    x$TumorPurity <- cos(0.6049872018 + 0.0001467884 * x$ESTIMATEScore)
  } else {
    x$TumorPurity <- suppressWarnings(as.numeric(x$TumorPurity))
  }
  
  x_use <- x %>%
    dplyr::select(
      Sample,
      Patient,
      ImmuneScore,
      StromalScore,
      ESTIMATEScore,
      TumorPurity
    ) %>%
    dplyr::distinct(Sample, .keep_all = TRUE)
  
  sm <- sample_map %>%
    dplyr::select(
      Sample,
      Patient,
      MO_DDRscore_raw,
      MO_DDRscore_group
    ) %>%
    dplyr::distinct(Patient, .keep_all = TRUE)
  
  n_match_sample <- sum(x_use$Sample %in% sm$Sample)
  n_match_patient <- sum(x_use$Patient %in% sm$Patient)
  
  cat("\nESTIMATE merge diagnostic:\n")
  cat("Existing ESTIMATE rows:", nrow(x_use), "\n")
  cat("Matched by Sample:", n_match_sample, "\n")
  cat("Matched by Patient:", n_match_patient, "\n")
  
  if (n_match_sample >= n_match_patient && n_match_sample > 0) {
    
    estimate_df <- x_use %>%
      dplyr::select(
        Sample,
        ImmuneScore,
        StromalScore,
        ESTIMATEScore,
        TumorPurity
      ) %>%
      dplyr::left_join(
        sm %>%
          dplyr::select(
            Sample,
            Patient,
            MO_DDRscore_raw,
            MO_DDRscore_group
          ),
        by = "Sample"
      )
    
  } else {
    
    estimate_df <- x_use %>%
      dplyr::select(
        Patient,
        ImmuneScore,
        StromalScore,
        ESTIMATEScore,
        TumorPurity
      ) %>%
      dplyr::left_join(
        sm %>%
          dplyr::select(
            Patient,
            Sample,
            MO_DDRscore_raw,
            MO_DDRscore_group
          ),
        by = "Patient"
      )
  }
  
  estimate_df <- estimate_df %>%
    dplyr::filter(!is.na(MO_DDRscore_group)) %>%
    dplyr::mutate(
      MO_DDRscore_group = factor(
        MO_DDRscore_group,
        levels = c("Low", "High")
      )
    )
  
  cat("Final matched ESTIMATE rows:", nrow(estimate_df), "\n")
  print(table(estimate_df$MO_DDRscore_group, useNA = "ifany"))
  
  if (nrow(estimate_df) == 0) {
    stop("ESTIMATE file was loaded but no samples matched current sample_map.")
  }
  
  estimate_df
}

estimate_df <- NULL
existing_estimate <- get_existing_estimate()

if (!is.null(existing_estimate)) {
  
  message("Using existing ESTIMATE file: ", existing_estimate)
  
  estimate_df <- tryCatch(
    standardize_existing_estimate(existing_estimate, sample_map),
    error = function(e) {
      message("Existing ESTIMATE file failed: ", e$message)
      NULL
    }
  )
}

if (is.null(estimate_df)) {
  
  estimate_df <- tryCatch(
    run_estimate_simple(expr_fig, sample_map),
    error = function(e) {
      message("ESTIMATE skipped: ", e$message)
      NULL
    }
  )
}


if (!is.null(existing_estimate)) {
  message("Using existing ESTIMATE file: ", existing_estimate)
  
  estimate_df <- data.table::fread(existing_estimate, data.table = FALSE, check.names = FALSE)
  estimate_df$Sample <- gsub("\\.", "-", estimate_df$Sample)
  
  if (!"TumorPurity" %in% colnames(estimate_df)) {
    estimate_df$TumorPurity <- cos(0.6049872018 + 0.0001467884 * as.numeric(estimate_df$ESTIMATEScore))
  }
  
  estimate_df <- estimate_df %>%
    dplyr::left_join(
      sample_map %>%
        dplyr::select(Sample, Patient, MO_DDRscore_raw, MO_DDRscore_group),
      by = "Sample"
    )
} else {
  estimate_df <- tryCatch(
    run_estimate_simple(expr_fig, sample_map),
    error = function(e) {
      message("ESTIMATE skipped: ", e$message)
      NULL
    }
  )
}

if (!is.null(estimate_df)) {
  
  save_csv(estimate_df, file.path(TAB3_DIR, "Fig3E_ESTIMATE_scores.csv"))
  
  
  estimate_df <- estimate_df %>%
    dplyr::mutate(
      Patient = dplyr::coalesce(
        as.character(Patient.x),
        as.character(Patient.y),
        patient_id(Sample)
      ),
      MO_DDRscore_raw = dplyr::coalesce(
        suppressWarnings(as.numeric(MO_DDRscore_raw.x)),
        suppressWarnings(as.numeric(MO_DDRscore_raw.y))
      ),
      MO_DDRscore_group = dplyr::coalesce(
        as.character(MO_DDRscore_group.x),
        as.character(MO_DDRscore_group.y)
      )
    ) %>%
    dplyr::select(
      Sample,
      Patient,
      StromalScore,
      ImmuneScore,
      ESTIMATEScore,
      TumorPurity,
      MO_DDRscore_raw,
      MO_DDRscore_group
    ) %>%
    dplyr::filter(!is.na(MO_DDRscore_group)) %>%
    dplyr::mutate(
      MO_DDRscore_group = factor(
        MO_DDRscore_group,
        levels = c("Low", "High")
      )
    ) %>%
    dplyr::distinct(Patient, .keep_all = TRUE)

  estimate_long <- estimate_df %>%
    dplyr::filter(!is.na(MO_DDRscore_group)) %>%
    tidyr::pivot_longer(
      cols = c(ImmuneScore, StromalScore, ESTIMATEScore, TumorPurity),
      names_to = "Feature",
      values_to = "Score"
    ) %>%
    dplyr::filter(is.finite(Score))
  
  estimate_stat <- estimate_long %>%
    dplyr::group_by(Feature) %>%
    dplyr::summarise(
      Median_Low = median(Score[MO_DDRscore_group == "Low"], na.rm = TRUE),
      Median_High = median(Score[MO_DDRscore_group == "High"], na.rm = TRUE),
      Diff = Median_High - Median_Low,
      P = tryCatch(wilcox.test(Score ~ MO_DDRscore_group)$p.value, error = function(e) NA_real_),
      y_pos = max(Score, na.rm = TRUE) * 1.08,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      FDR = p.adjust(P, "BH"),
      Label = sig_label(FDR),
      x_pos = 1.5
    )
  
  save_csv(estimate_stat, file.path(TAB3_DIR, "Fig3E_ESTIMATE_stats.csv"))
  
  p_est <- ggplot(
    estimate_long,
    aes(x = MO_DDRscore_group, y = Score, fill = MO_DDRscore_group)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.62, linewidth = 0.35) +
    geom_jitter(width = 0.14, size = 0.5, alpha = 0.30) +
    geom_text(
      data = estimate_stat,
      aes(x = x_pos, y = y_pos, label = Label),
      inherit.aes = FALSE,
      size = 4,
      fontface = "bold"
    ) +
    facet_wrap(~Feature, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c(Low = LOW_COL, High = HIGH_COL)) +
    theme_pub(11.5) +
    theme(legend.position = "none", strip.text = element_text(face = "bold")) +
    labs(
      x = NULL,
      y = "ESTIMATE-derived value",
      title = "ESTIMATE tumor microenvironment features"
    )
  
  safe_ggsave(
    file.path(FIG3_DIR, "Fig3E1_ESTIMATE_tumor_microenvironment.pdf"),
    p_est,
    9.2,
    4.2
  )
}

############################################################
## 11. Fig3E2: CIBERSORT LM22 immune infiltration
############################################################

cat("\nRunning Fig3E2 CIBERSORT LM22 immune infiltration...\n")

run_cibersort_available <- file.exists(CIBERSORT_R_FILE) && file.exists(LM22_FILE)

if (!run_cibersort_available) {
  
  message("CIBERSORT.R or LM22.txt not found. Fig3E2 CIBERSORT skipped.")
  message("Expected files:")
  message("  ", CIBERSORT_R_FILE)
  message("  ", LM22_FILE)
  
} else {
  
  CIBER_DIR <- file.path(TAB3_DIR, "CIBERSORT_work")
  dir.create(CIBER_DIR, recursive = TRUE, showWarnings = FALSE)
  
  ciber_input <- file.path(CIBER_DIR, "LUAD_CIBERSORT_input_TPM.txt")
  
  expr_ciber <- expr_tpm_fig[, sample_map$Sample, drop = FALSE]
  rownames(expr_ciber) <- clean_gene(rownames(expr_ciber))
  expr_ciber <- expr_ciber[!is.na(rownames(expr_ciber)) & rownames(expr_ciber) != "", , drop = FALSE]
  expr_ciber <- expr_ciber[!duplicated(rownames(expr_ciber)), , drop = FALSE]
  expr_ciber[!is.finite(expr_ciber)] <- 0
  
  ciber_df <- data.frame(
    GeneSymbol = rownames(expr_ciber),
    expr_ciber,
    check.names = FALSE
  )
  
  data.table::fwrite(ciber_df, ciber_input, sep = "\t", quote = FALSE, na = "0")
  
  source(CIBERSORT_R_FILE)
  
  if (!exists("CIBERSORT")) {
    stop("CIBERSORT function was not found after sourcing CIBERSORT.R.")
  }
  
  ciber_res <- tryCatch(
    CIBERSORT(
      sig_matrix = LM22_FILE,
      mixture_file = ciber_input,
      perm = 100,
      QN = FALSE
    ),
    error = function(e) {
      message("CIBERSORT failed: ", e$message)
      NULL
    }
  )
  
  if (!is.null(ciber_res)) {
    
    ciber_res <- as.data.frame(ciber_res, check.names = FALSE)
    ciber_res$Sample <- rownames(ciber_res)
    ciber_res$Sample <- gsub("\\.", "-", ciber_res$Sample)
    
    cell_cols <- setdiff(
      colnames(ciber_res),
      c("P-value", "P.value", "Pvalue", "Correlation", "RMSE", "Sample")
    )
    
    ciber_df2 <- ciber_res %>%
      dplyr::left_join(
        sample_map %>%
          dplyr::select(Sample, Patient, MO_DDRscore_raw, MO_DDRscore_group),
        by = "Sample"
      )
    
    save_csv(ciber_df2, file.path(TAB3_DIR, "Fig3E_CIBERSORT_LM22_scores.csv"))
    
    ciber_long <- ciber_df2 %>%
      dplyr::filter(!is.na(MO_DDRscore_group)) %>%
      tidyr::pivot_longer(
        cols = all_of(cell_cols),
        names_to = "CellType",
        values_to = "Fraction"
      ) %>%
      dplyr::filter(is.finite(Fraction))
    
    ciber_stat <- ciber_long %>%
      dplyr::group_by(CellType) %>%
      dplyr::summarise(
        Median_Low = median(Fraction[MO_DDRscore_group == "Low"], na.rm = TRUE),
        Median_High = median(Fraction[MO_DDRscore_group == "High"], na.rm = TRUE),
        Diff = Median_High - Median_Low,
        P = tryCatch(wilcox.test(Fraction ~ MO_DDRscore_group)$p.value, error = function(e) NA_real_),
        y_pos = max(Fraction, na.rm = TRUE) * 1.08,
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        FDR = p.adjust(P, "BH"),
        Label = sig_label(FDR)
      ) %>%
      dplyr::arrange(FDR, dplyr::desc(abs(Diff)))
    
    save_csv(
      ciber_stat,
      file.path(TAB3_DIR, "Fig3E_CIBERSORT_immune_cell_stats.csv")
    )
    
    show_cells <- ciber_stat %>%
      dplyr::slice_head(n = min(TOP_CIBERSORT_CELLS, nrow(.))) %>%
      dplyr::pull(CellType)
    
    label_df <- ciber_stat %>%
      dplyr::filter(CellType %in% show_cells) %>%
      dplyr::mutate(
        CellType = factor(CellType, levels = show_cells),
        x_pos = 1.5
      )
    
    p_ciber <- ciber_long %>%
      dplyr::filter(CellType %in% show_cells) %>%
      dplyr::mutate(CellType = factor(CellType, levels = show_cells)) %>%
      ggplot(aes(x = MO_DDRscore_group, y = Fraction, fill = MO_DDRscore_group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.62, linewidth = 0.32) +
      geom_jitter(width = 0.12, size = 0.45, alpha = 0.25) +
      geom_text(
        data = label_df,
        aes(x = x_pos, y = y_pos, label = Label),
        inherit.aes = FALSE,
        size = 3.2,
        fontface = "bold"
      ) +
      facet_wrap(~CellType, scales = "free_y", ncol = 4) +
      scale_fill_manual(values = c(Low = LOW_COL, High = HIGH_COL)) +
      theme_pub(10) +
      theme(legend.position = "none", strip.text = element_text(face = "bold")) +
      labs(
        x = NULL,
        y = "Estimated cell fraction",
        title = "CIBERSORT LM22 immune cell infiltration"
      )
    
    safe_ggsave(
      file.path(FIG3_DIR, "Fig3E2_CIBERSORT_LM22_immune_cell_infiltration.pdf"),
      p_ciber,
      8.5,
      6.2
    )
  }
}

############################################################
## 12. Fig3F: DepMap/GDSC drug sensitivity
############################################################

cat("\nRunning Fig3F DepMap/GDSC drug sensitivity analysis...\n")

GDSC_DIR <- file.path(TAB3_DIR, "GDSC_DepMap_work")
dir.create(GDSC_DIR, recursive = TRUE, showWarnings = FALSE)

if (file.exists(DEPMAP_EXPR_FILE) && file.exists(MODEL_MAP_FILE) && file.exists(GDSC_FILE)) {
  
  dep_raw <- data.table::fread(
    DEPMAP_EXPR_FILE,
    data.table = FALSE,
    check.names = FALSE
  )
  
  model_col <- pick_col(
    dep_raw,
    c("ModelID", "DepMap_ID", "DepMapID", "model_id", "ID")
  )
  
  if (is.na(model_col)) {
    model_col <- colnames(dep_raw)[1]
  }
  
  dep_model_id <- as.character(dep_raw[[model_col]])
  
  gene_cols_raw <- setdiff(colnames(dep_raw), model_col)
  gene_symbols <- clean_depmap_gene(gene_cols_raw)
  
  keep_gene_cols <- !is.na(gene_symbols) & gene_symbols != ""
  gene_cols_raw <- gene_cols_raw[keep_gene_cols]
  gene_symbols <- gene_symbols[keep_gene_cols]
  
  dep_expr_mat <- as.matrix(dep_raw[, gene_cols_raw, drop = FALSE])
  storage.mode(dep_expr_mat) <- "numeric"
  dep_expr_mat[!is.finite(dep_expr_mat)] <- NA
  
  rownames(dep_expr_mat) <- dep_model_id
  
  dep_gene_by_model_sum <- rowsum(
    t(dep_expr_mat),
    group = gene_symbols,
    reorder = FALSE,
    na.rm = TRUE
  )
  
  gene_count <- as.numeric(table(gene_symbols)[rownames(dep_gene_by_model_sum)])
  dep_gene_by_model <- dep_gene_by_model_sum / gene_count
  
  colnames(dep_gene_by_model) <- dep_model_id
  
  cat("DepMap expression loaded:", nrow(dep_gene_by_model), "genes x", ncol(dep_gene_by_model), "models\n")
  
  deg_for_sig <- deg %>%
    dplyr::mutate(
      Gene = clean_gene(Gene),
      FDR = as.numeric(FDR),
      logFC = as.numeric(logFC)
    ) %>%
    dplyr::filter(
      !is.na(Gene),
      Gene != "",
      is.finite(FDR),
      is.finite(logFC)
    ) %>%
    dplyr::arrange(FDR, dplyr::desc(abs(logFC)))
  
  up_sig <- deg_for_sig %>%
    dplyr::filter(
      FDR < DEG_FDR,
      logFC > DEG_LOGFC,
      Gene %in% rownames(dep_gene_by_model)
    ) %>%
    dplyr::arrange(FDR, dplyr::desc(logFC)) %>%
    dplyr::slice_head(n = 100) %>%
    dplyr::pull(Gene) %>%
    unique()
  
  down_sig <- deg_for_sig %>%
    dplyr::filter(
      FDR < DEG_FDR,
      logFC < -DEG_LOGFC,
      Gene %in% rownames(dep_gene_by_model)
    ) %>%
    dplyr::arrange(FDR, logFC) %>%
    dplyr::slice_head(n = 100) %>%
    dplyr::pull(Gene) %>%
    unique()
  
  if (length(up_sig) < 10) {
    up_sig <- deg_for_sig %>%
      dplyr::filter(logFC > 0, Gene %in% rownames(dep_gene_by_model)) %>%
      dplyr::arrange(FDR, dplyr::desc(logFC)) %>%
      dplyr::slice_head(n = 50) %>%
      dplyr::pull(Gene) %>%
      unique()
  }
  
  if (length(down_sig) < 10) {
    down_sig <- deg_for_sig %>%
      dplyr::filter(logFC < 0, Gene %in% rownames(dep_gene_by_model)) %>%
      dplyr::arrange(FDR, logFC) %>%
      dplyr::slice_head(n = 50) %>%
      dplyr::pull(Gene) %>%
      unique()
  }
  
  if (length(up_sig) >= 5) {
    
    save_csv(
      data.frame(
        Gene = c(up_sig, down_sig),
        Direction = c(
          rep("High_up_signature", length(up_sig)),
          rep("High_down_signature", length(down_sig))
        )
      ),
      file.path(GDSC_DIR, "CellLine_MO_DDRscore_signature_genes.csv")
    )
    
    score_genes <- unique(c(up_sig, down_sig))
    dep_score_expr <- dep_gene_by_model[score_genes, , drop = FALSE]
    
    dep_score_expr_z <- t(apply(dep_score_expr, 1, z_vec))
    rownames(dep_score_expr_z) <- rownames(dep_score_expr)
    colnames(dep_score_expr_z) <- colnames(dep_score_expr)
    
    up_score <- colMeans(
      dep_score_expr_z[intersect(up_sig, rownames(dep_score_expr_z)), , drop = FALSE],
      na.rm = TRUE
    )
    
    if (length(down_sig) >= 5) {
      down_score <- colMeans(
        dep_score_expr_z[intersect(down_sig, rownames(dep_score_expr_z)), , drop = FALSE],
        na.rm = TRUE
      )
    } else {
      down_score <- rep(0, length(up_score))
      names(down_score) <- names(up_score)
    }
    
    cell_score <- data.frame(
      ModelID = names(up_score),
      UpSignatureScore = as.numeric(up_score),
      DownSignatureScore = as.numeric(down_score[names(up_score)]),
      CellLine_MO_DDRscore = as.numeric(up_score - down_score[names(up_score)]),
      stringsAsFactors = FALSE
    )
    
    cell_score$CellLine_MO_Group <- ifelse(
      cell_score$CellLine_MO_DDRscore >= median(cell_score$CellLine_MO_DDRscore, na.rm = TRUE),
      "High",
      "Low"
    )
    
    cell_score$CellLine_MO_Group <- factor(
      cell_score$CellLine_MO_Group,
      levels = c("Low", "High")
    )
    
    save_csv(
      cell_score,
      file.path(GDSC_DIR, "DepMap_CellLine_MO_DDRscore.csv")
    )
    
    model_map <- data.table::fread(
      MODEL_MAP_FILE,
      data.table = FALSE,
      check.names = FALSE
    )
    
    modelid_col <- pick_col(
      model_map,
      c("ModelID", "DepMap_ID", "DepMapID", "model_id")
    )
    
    sanger_col <- pick_col(
      model_map,
      c("SangerModelID", "SANGER_MODEL_ID", "Sanger_Model_ID", "sanger_model_id")
    )
    
    if (is.na(modelid_col) || is.na(sanger_col)) {
      stop("Model.csv must contain ModelID and SangerModelID columns.")
    }
    
    lineage_cols <- grep(
      "lineage|tissue|cancer|oncotree|disease|TCGA|primary",
      colnames(model_map),
      ignore.case = TRUE,
      value = TRUE
    )
    
    model_map_use <- model_map %>%
      dplyr::mutate(
        ModelID = as.character(.data[[modelid_col]]),
        SangerModelID = as.character(.data[[sanger_col]])
      ) %>%
      dplyr::select(
        ModelID,
        SangerModelID,
        dplyr::all_of(lineage_cols)
      ) %>%
      dplyr::distinct(ModelID, SangerModelID, .keep_all = TRUE)
    
    save_csv(
      model_map_use,
      file.path(GDSC_DIR, "ModelID_SangerModelID_mapping_used.csv")
    )
    
    gdsc_raw <- as.data.frame(
      readxl::read_excel(GDSC_FILE),
      check.names = FALSE
    )
    
    gdsc_sanger_col <- pick_col(
      gdsc_raw,
      c("SANGER_MODEL_ID", "SangerModelID", "Sanger_Model_ID", "sanger_model_id")
    )
    
    drug_col <- pick_col(
      gdsc_raw,
      c("DRUG_NAME", "DrugName", "Drug_Name", "drug_name", "Drug")
    )
    
    lnic50_col <- pick_col(
      gdsc_raw,
      c("LN_IC50", "ln_ic50", "IC50", "LNIC50")
    )
    
    auc_col <- pick_col(
      gdsc_raw,
      c("AUC", "auc")
    )
    
    if (is.na(gdsc_sanger_col) || is.na(drug_col) || is.na(lnic50_col)) {
      stop("GDSC file must contain SANGER_MODEL_ID, DRUG_NAME and LN_IC50 columns.")
    }
    
    gdsc_use <- gdsc_raw %>%
      dplyr::mutate(
        SangerModelID = as.character(.data[[gdsc_sanger_col]]),
        Drug = as.character(.data[[drug_col]]),
        LN_IC50 = suppressWarnings(as.numeric(.data[[lnic50_col]])),
        AUC = if (!is.na(auc_col)) suppressWarnings(as.numeric(.data[[auc_col]])) else NA_real_
      ) %>%
      dplyr::filter(
        !is.na(SangerModelID),
        !is.na(Drug),
        is.finite(LN_IC50)
      ) %>%
      dplyr::select(SangerModelID, Drug, LN_IC50, AUC, dplyr::everything())
    
    save_csv(
      gdsc_use,
      file.path(GDSC_DIR, "GDSC_response_cleaned.csv")
    )
    
    gdsc_merged <- gdsc_use %>%
      dplyr::left_join(model_map_use, by = "SangerModelID") %>%
      dplyr::left_join(cell_score, by = "ModelID") %>%
      dplyr::filter(
        !is.na(ModelID),
        is.finite(CellLine_MO_DDRscore),
        !is.na(CellLine_MO_Group),
        is.finite(LN_IC50)
      )
    
    save_csv(
      gdsc_merged,
      file.path(GDSC_DIR, "GDSC_DepMap_CellLine_MO_DDRscore_merged_all.csv")
    )
    
    analysis_df_all <- gdsc_merged %>%
      dplyr::mutate(AnalysisScope = "All_mapped_cell_lines")
    
    analysis_df_lung <- NULL
    
    if (length(lineage_cols) > 0) {
      
      lung_flag <- apply(
        gdsc_merged[, lineage_cols, drop = FALSE],
        1,
        function(z) {
          any(grepl(
            "LUAD|LUSC|NSCLC|LUNG|lung|Lung|lung adenocarcinoma|lung carcinoma",
            paste(z, collapse = " "),
            ignore.case = TRUE
          ))
        }
      )
      
      tmp_lung <- gdsc_merged[lung_flag, , drop = FALSE]
      
      if (nrow(tmp_lung) > 0) {
        analysis_df_lung <- tmp_lung %>%
          dplyr::mutate(AnalysisScope = "Lung_lineage_cell_lines")
      }
    }
    
    if (!is.null(analysis_df_lung) && length(unique(analysis_df_lung$ModelID)) >= 20) {
      analysis_df <- analysis_df_lung
      main_scope <- "Lung lineage cell lines"
    } else {
      analysis_df <- analysis_df_all
      main_scope <- "All mapped cell lines"
    }
    
    save_csv(
      analysis_df,
      file.path(GDSC_DIR, "GDSC_DepMap_analysis_dataset_used.csv")
    )
    
    target_drug_keywords <- c(
      "Cisplatin",
      "Carboplatin",
      "Oxaliplatin",
      "Paclitaxel",
      "Docetaxel",
      "Gemcitabine",
      "Etoposide",
      "Topotecan",
      "Irinotecan",
      "Camptothecin",
      "Olaparib",
      "Talazoparib",
      "Niraparib",
      "Rucaparib",
      "AZD2281",
      "Veliparib"
    )
    
    drug_pattern <- paste(target_drug_keywords, collapse = "|")
    
    analysis_focus <- analysis_df %>%
      dplyr::filter(
        grepl(drug_pattern, Drug, ignore.case = TRUE)
      )
    
    if (nrow(analysis_focus) == 0) {
      message("No focused DNA-damage/LUAD drugs matched. Use top tested drugs by sample size.")
      
      top_drugs <- analysis_df %>%
        dplyr::count(Drug, sort = TRUE) %>%
        dplyr::slice_head(n = 20) %>%
        dplyr::pull(Drug)
      
      analysis_focus <- analysis_df %>%
        dplyr::filter(Drug %in% top_drugs)
    }
    
    save_csv(
      analysis_focus,
      file.path(GDSC_DIR, "GDSC_focused_drug_records_used.csv")
    )
    
    drug_stats <- dplyr::bind_rows(
      lapply(sort(unique(analysis_focus$Drug)), function(dg) {
        
        df <- analysis_focus %>%
          dplyr::filter(Drug == dg) %>%
          dplyr::filter(
            is.finite(LN_IC50),
            is.finite(CellLine_MO_DDRscore),
            !is.na(CellLine_MO_Group)
          )
        
        n_total <- nrow(df)
        n_model <- length(unique(df$ModelID))
        n_low <- sum(df$CellLine_MO_Group == "Low")
        n_high <- sum(df$CellLine_MO_Group == "High")
        
        if (n_total < 10 || n_low < 3 || n_high < 3) {
          return(data.frame(
            Drug = dg,
            N = n_total,
            N_Model = n_model,
            N_Low = n_low,
            N_High = n_high,
            Median_Low = NA_real_,
            Median_High = NA_real_,
            Diff_HighMinusLow = NA_real_,
            WilcoxP = NA_real_,
            SpearmanRho = NA_real_,
            SpearmanP = NA_real_,
            stringsAsFactors = FALSE
          ))
        }
        
        cor_res <- safe_cor(df$CellLine_MO_DDRscore, df$LN_IC50)
        
        data.frame(
          Drug = dg,
          N = n_total,
          N_Model = n_model,
          N_Low = n_low,
          N_High = n_high,
          Median_Low = median(df$LN_IC50[df$CellLine_MO_Group == "Low"], na.rm = TRUE),
          Median_High = median(df$LN_IC50[df$CellLine_MO_Group == "High"], na.rm = TRUE),
          Diff_HighMinusLow =
            median(df$LN_IC50[df$CellLine_MO_Group == "High"], na.rm = TRUE) -
            median(df$LN_IC50[df$CellLine_MO_Group == "Low"], na.rm = TRUE),
          WilcoxP = tryCatch(wilcox.test(LN_IC50 ~ CellLine_MO_Group, data = df)$p.value, error = function(e) NA_real_),
          SpearmanRho = cor_res$Rho,
          SpearmanP = cor_res$CorP,
          stringsAsFactors = FALSE
        )
      })
    )
    
    drug_stats <- drug_stats %>%
      dplyr::mutate(
        WilcoxFDR = p.adjust(WilcoxP, method = "BH"),
        SpearmanFDR = p.adjust(SpearmanP, method = "BH"),
        SensitivityDirection = dplyr::case_when(
          is.na(Diff_HighMinusLow) ~ "NA",
          Diff_HighMinusLow < 0 ~ "High score more sensitive",
          Diff_HighMinusLow > 0 ~ "High score more resistant",
          TRUE ~ "No difference"
        ),
        WilcoxLabel = sig_label(WilcoxFDR)
      ) %>%
      dplyr::arrange(WilcoxFDR, SpearmanFDR, dplyr::desc(N))
    
    save_csv(
      drug_stats,
      file.path(GDSC_DIR, "GDSC_CellLine_MO_DDRscore_drug_sensitivity_stats.csv")
    )
    
    plot_drugs <- drug_stats %>%
      dplyr::filter(
        is.finite(WilcoxP),
        N >= 10,
        N_Low >= 3,
        N_High >= 3
      ) %>%
      dplyr::arrange(WilcoxFDR, dplyr::desc(abs(Diff_HighMinusLow))) %>%
      dplyr::slice_head(n = 12) %>%
      dplyr::pull(Drug)
    
    if (length(plot_drugs) == 0) {
      plot_drugs <- analysis_focus %>%
        dplyr::count(Drug, sort = TRUE) %>%
        dplyr::slice_head(n = 12) %>%
        dplyr::pull(Drug)
    }
    
    plot_df <- analysis_focus %>%
      dplyr::filter(Drug %in% plot_drugs) %>%
      dplyr::mutate(
        Drug = factor(Drug, levels = plot_drugs),
        CellLine_MO_Group = factor(CellLine_MO_Group, levels = c("Low", "High"))
      )
    
    label_df <- plot_df %>%
      dplyr::group_by(Drug) %>%
      dplyr::summarise(
        y_pos = max(LN_IC50, na.rm = TRUE) * 1.08,
        .groups = "drop"
      ) %>%
      dplyr::left_join(
        drug_stats %>%
          dplyr::select(Drug, WilcoxLabel, WilcoxFDR, SensitivityDirection),
        by = "Drug"
      ) %>%
      dplyr::mutate(
        x_pos = 1.5,
        Drug = factor(Drug, levels = plot_drugs)
      )
    
    p_gdsc_box <- ggplot(
      plot_df,
      aes(x = CellLine_MO_Group, y = LN_IC50, fill = CellLine_MO_Group)
    ) +
      geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.62, linewidth = 0.35) +
      geom_jitter(width = 0.14, size = 0.65, alpha = 0.35) +
      geom_text(
        data = label_df,
        aes(x = x_pos, y = y_pos, label = WilcoxLabel),
        inherit.aes = FALSE,
        size = 3.6,
        fontface = "bold"
      ) +
      facet_wrap(~Drug, scales = "free_y", ncol = 4) +
      scale_fill_manual(values = c(Low = LOW_COL, High = HIGH_COL)) +
      theme_pub(10.5) +
      theme(
        legend.position = "none",
        strip.text = element_text(face = "bold")
      ) +
      labs(
        x = "Cell-line MO-DDRscore group",
        y = "GDSC LN_IC50",
        title = "GDSC drug sensitivity associated with cell-line MO-DDRscore",
        subtitle = paste0(main_scope, "; lower LN_IC50 indicates higher sensitivity")
      )
    
    safe_ggsave(
      file.path(FIG3_DIR, "Fig3F1_GDSC_drug_sensitivity_boxplot.pdf"),
      p_gdsc_box,
      10.5,
      7.8
    )
    
    cor_plot_df <- drug_stats %>%
      dplyr::filter(
        Drug %in% plot_drugs,
        is.finite(SpearmanRho)
      ) %>%
      dplyr::mutate(
        Drug = factor(Drug, levels = rev(plot_drugs)),
        negLogFDR = -log10(SpearmanFDR + 1e-300),
        Association = dplyr::case_when(
          SpearmanRho < 0 ~ "Higher score -> lower LN_IC50",
          SpearmanRho > 0 ~ "Higher score -> higher LN_IC50",
          TRUE ~ "No association"
        )
      )
    
    if (nrow(cor_plot_df) > 0) {
      
      p_gdsc_cor <- ggplot(
        cor_plot_df,
        aes(x = SpearmanRho, y = Drug)
      ) +
        geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
        geom_point(
          aes(size = negLogFDR, color = Association),
          alpha = 0.9
        ) +
        scale_color_manual(
          values = c(
            "Higher score -> lower LN_IC50" = HIGH_COL,
            "Higher score -> higher LN_IC50" = LOW_COL,
            "No association" = "grey60"
          )
        ) +
        theme_pub(11) +
        theme(
          panel.grid.major.y = element_blank(),
          legend.title = element_blank()
        ) +
        labs(
          x = "Spearman correlation with LN_IC50",
          y = NULL,
          size = "-log10(FDR)",
          title = "Continuous association between cell-line MO-DDRscore and GDSC LN_IC50"
        )
      
      safe_ggsave(
        file.path(FIG3_DIR, "Fig3F2_GDSC_score_IC50_correlation_dotplot.pdf"),
        p_gdsc_cor,
        8.2,
        5.6
      )
    }
  }
  
} else {
  message("DepMap/GDSC files not complete. Fig3F drug sensitivity skipped.")
  message("Expected:")
  message("  ", DEPMAP_EXPR_FILE)
  message("  ", MODEL_MAP_FILE)
  message("  ", GDSC_FILE)
}

############################################################
## 13. Fig3G: KM survival by MO-DDRscore group
############################################################

cat("\nRunning Fig3G KM survival...\n")

surv_df <- score_clin %>%
  dplyr::filter(
    is.finite(time),
    time > 0,
    status %in% c(0, 1),
    !is.na(MO_DDRscore_group)
  ) %>%
  dplyr::mutate(
    MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High"))
  )

save_csv(surv_df, file.path(TAB3_DIR, "Fig3G_survival_input.csv"))

if (nrow(surv_df) >= 50 && sum(surv_df$status == 1) >= 10 && has_pkg("survminer")) {
  
  fit <- survival::survfit(
    survival::Surv(time, status) ~ MO_DDRscore_group,
    data = surv_df
  )
  
  p_km <- survminer::ggsurvplot(
    fit,
    data = surv_df,
    pval = TRUE,
    conf.int = FALSE,
    risk.table = TRUE,
    palette = c(LOW_COL, HIGH_COL),
    legend.title = "",
    legend.labs = c("Low", "High"),
    ggtheme = theme_bw(base_size = 12),
    title = paste0(FOCUS_CANCER, " MO-DDRscore survival")
  )$plot +
    labs(
      x = "Time (days)",
      y = "Overall survival probability"
    ) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0)
    )
  
  safe_ggsave(
    file.path(FIG3_DIR, "Fig3G_MO_DDRscore_KM_survival.pdf"),
    p_km,
    6.2,
    5.2
  )
  
  cox_score <- survival::coxph(
    survival::Surv(time, status) ~ MO_DDRscore_raw,
    data = surv_df
  )
  
  cox_group <- survival::coxph(
    survival::Surv(time, status) ~ MO_DDRscore_group,
    data = surv_df
  )
  
  save_csv(
    data.frame(
      Model = c("MO_DDRscore_raw", "MO_DDRscore_group"),
      HR = c(
        summary(cox_score)$coefficients[1, "exp(coef)"],
        summary(cox_group)$coefficients[1, "exp(coef)"]
      ),
      P = c(
        summary(cox_score)$coefficients[1, "Pr(>|z|)"],
        summary(cox_group)$coefficients[1, "Pr(>|z|)"]
      )
    ),
    file.path(TAB3_DIR, "Fig3G_survival_cox_summary.csv")
  )
}

############################################################
## 14. Supplementary: checkpoint expression
############################################################

cat("\nRunning supplementary checkpoint expression...\n")

checkpoint_genes <- intersect(
  c(
    "CD274", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2",
    "PDCD1LG2", "ICOS", "IDO1", "CD80", "CD86",
    "TNFRSF9", "TNFRSF4", "CD40", "CD40LG",
    "VSIR", "SIGLEC15", "BTLA", "ENTPD1", "CD70"
  ),
  rownames(expr_fig)
)

if (length(checkpoint_genes) >= 3) {
  
  chk_df <- as.data.frame(t(expr_fig[checkpoint_genes, , drop = FALSE]), check.names = FALSE)
  chk_df$Sample <- rownames(chk_df)
  
  chk_df <- chk_df %>%
    dplyr::left_join(
      sample_map %>% dplyr::select(Sample, MO_DDRscore_group),
      by = "Sample"
    )
  
  chk_long <- chk_df %>%
    tidyr::pivot_longer(
      cols = all_of(checkpoint_genes),
      names_to = "Gene",
      values_to = "Expression"
    ) %>%
    dplyr::filter(is.finite(Expression), !is.na(MO_DDRscore_group))
  
  chk_stat <- chk_long %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      Median_Low = median(Expression[MO_DDRscore_group == "Low"], na.rm = TRUE),
      Median_High = median(Expression[MO_DDRscore_group == "High"], na.rm = TRUE),
      Diff = Median_High - Median_Low,
      P = tryCatch(wilcox.test(Expression ~ MO_DDRscore_group)$p.value, error = function(e) NA_real_),
      y_pos = max(Expression, na.rm = TRUE) * 1.08,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      FDR = p.adjust(P, "BH"),
      Label = sig_label(FDR)
    ) %>%
    dplyr::arrange(FDR, dplyr::desc(abs(Diff)))
  
  save_csv(chk_stat, file.path(SUPP_TAB_DIR, "FigS1_checkpoint_expression_stats.csv"))
  
  showg <- chk_stat %>%
    dplyr::slice_head(n = min(TOP_CHECKPOINT_GENES, nrow(.))) %>%
    dplyr::pull(Gene)
  
  chk_label_df <- chk_stat %>%
    dplyr::filter(Gene %in% showg) %>%
    dplyr::mutate(
      Gene = factor(Gene, levels = showg),
      x_pos = 1.5
    )
  
  p_chk <- chk_long %>%
    dplyr::filter(Gene %in% showg) %>%
    dplyr::mutate(Gene = factor(Gene, levels = showg)) %>%
    ggplot(aes(MO_DDRscore_group, Expression, fill = MO_DDRscore_group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.62, linewidth = 0.32) +
    geom_jitter(width = 0.12, size = 0.45, alpha = 0.25) +
    geom_text(
      data = chk_label_df,
      aes(x = x_pos, y = y_pos, label = Label),
      inherit.aes = FALSE,
      size = 3.3,
      fontface = "bold"
    ) +
    facet_wrap(~Gene, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = c(Low = LOW_COL, High = HIGH_COL)) +
    theme_pub(10) +
    theme(legend.position = "none", strip.text = element_text(face = "bold")) +
    labs(
      x = NULL,
      y = "log2(TPM + 1)",
      title = "Immune checkpoint expression"
    )
  
  safe_ggsave(
    file.path(SUPP_DIR, "FigS1_checkpoint_expression_boxplot.pdf"),
    p_chk,
    8.2,
    6
  )
}

############################################################
## 15. Supplementary: clinical feature distribution
############################################################

cat("\nRunning supplementary clinical distribution...\n")

clinical_df <- score_clin %>%
  dplyr::filter(!is.na(MO_DDRscore_group)) %>%
  dplyr::mutate(
    MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High")),
    Age = suppressWarnings(as.numeric(age)),
    Gender = dplyr::case_when(
      toupper(as.character(gender)) %in% c("MALE", "M") ~ "Male",
      toupper(as.character(gender)) %in% c("FEMALE", "F") ~ "Female",
      TRUE ~ "Unknown"
    ),
    Stage = as.character(stage),
    Stage = ifelse(is.na(Stage) | Stage == "", "Unknown", Stage)
  )

save_csv(
  clinical_df,
  file.path(SUPP_TAB_DIR, "FigS2_clinical_feature_input.csv")
)

age_df <- clinical_df %>% dplyr::filter(is.finite(Age))

p_age_val <- tryCatch(
  wilcox.test(Age ~ MO_DDRscore_group, data = age_df)$p.value,
  error = function(e) NA_real_
)

gender_df <- clinical_df %>% dplyr::filter(Gender != "Unknown")
gender_tab <- table(gender_df$MO_DDRscore_group, gender_df$Gender)

p_gender_val <- tryCatch(
  fisher.test(gender_tab)$p.value,
  error = function(e) {
    tryCatch(chisq.test(gender_tab)$p.value, error = function(e2) NA_real_)
  }
)

stage_df <- clinical_df %>% dplyr::filter(Stage != "Unknown")
stage_tab <- table(stage_df$MO_DDRscore_group, stage_df$Stage)

p_stage_val <- tryCatch(
  fisher.test(stage_tab)$p.value,
  error = function(e) {
    tryCatch(chisq.test(stage_tab)$p.value, error = function(e2) NA_real_)
  }
)

p_age <- ggplot(age_df, aes(x = MO_DDRscore_group, y = Age, fill = MO_DDRscore_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.85, linewidth = 0.35) +
  geom_jitter(width = 0.14, size = 0.55, alpha = 0.35) +
  scale_fill_manual(values = c(Low = LOW_COL, High = HIGH_COL)) +
  theme_pub(11) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Age", title = paste0("Age\n", format_p(p_age_val)))

gender_plot_df <- gender_df %>%
  dplyr::count(MO_DDRscore_group, Gender) %>%
  dplyr::group_by(MO_DDRscore_group) %>%
  dplyr::mutate(Fraction = n / sum(n)) %>%
  dplyr::ungroup()

p_gender <- ggplot(gender_plot_df, aes(x = MO_DDRscore_group, y = Fraction, fill = Gender)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  theme_pub(11) +
  theme(legend.position = "right", legend.title = element_blank()) +
  labs(x = NULL, y = "Proportion", title = paste0("Gender\n", format_p(p_gender_val)))

stage_plot_df <- stage_df %>%
  dplyr::count(MO_DDRscore_group, Stage) %>%
  dplyr::group_by(MO_DDRscore_group) %>%
  dplyr::mutate(Fraction = n / sum(n)) %>%
  dplyr::ungroup()

p_stage <- ggplot(stage_plot_df, aes(x = MO_DDRscore_group, y = Fraction, fill = Stage)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  theme_pub(11) +
  theme(legend.position = "right", legend.title = element_blank()) +
  labs(x = NULL, y = "Proportion", title = paste0("Stage\n", format_p(p_stage_val)))

p_clinical <- p_age + p_gender + p_stage +
  patchwork::plot_layout(widths = c(1, 1.1, 1.5))

safe_ggsave(
  file.path(SUPP_DIR, "FigS2_clinical_feature_distribution.pdf"),
  p_clinical,
  11,
  3.8
)

############################################################
## 16. Supplementary: gene-level CNA altered fraction
############################################################

cat("\nRunning supplementary gene-level CNA altered fraction...\n")

is_gistic_cnv <- function(m) {
  vals <- unique(as.numeric(m[is.finite(m)]))
  vals <- vals[!is.na(vals)]
  length(vals) > 0 && all(vals %in% c(-2, -1, 0, 1, 2))
}

make_cnv_binary <- function(m, type = c("alter", "amp", "del")) {
  type <- match.arg(type)
  m <- as.matrix(m)
  storage.mode(m) <- "numeric"
  
  if (is_gistic_cnv(m)) {
    if (type == "alter") return(abs(m) >= 1)
    if (type == "amp") return(m >= 1)
    if (type == "del") return(m <= -1)
  }
  
  if (type == "alter") return(abs(m) > CNV_THRESHOLD)
  if (type == "amp") return(m > CNV_THRESHOLD)
  if (type == "del") return(m < -CNV_THRESHOLD)
}

collapse_cnv_to_patient <- function(cnv_mat, patients_keep) {
  if (is.null(cnv_mat)) return(NULL)
  
  cnv_pat <- patient_id(colnames(cnv_mat))
  keep_cols <- cnv_pat %in% patients_keep
  
  if (sum(keep_cols) == 0) return(NULL)
  
  m <- cnv_mat[, keep_cols, drop = FALSE]
  colnames(m) <- cnv_pat[keep_cols]
  
  if (any(duplicated(colnames(m)))) {
    sample_tab <- table(colnames(m))
    m_t_sum <- rowsum(t(m), group = colnames(m), reorder = FALSE)
    m_t_mean <- m_t_sum / as.numeric(sample_tab[rownames(m_t_sum)])
    m <- t(m_t_mean)
  }
  
  m
}

cnv_features <- data.frame(
  Patient = sample_map$Patient,
  AllGene_CNA_Fraction = NA_real_,
  AllGene_Amp_Fraction = NA_real_,
  AllGene_Del_Fraction = NA_real_,
  DDRGene_CNA_Fraction = NA_real_,
  DDRGene_Amp_Fraction = NA_real_,
  DDRGene_Del_Fraction = NA_real_
)

cnv_patient <- collapse_cnv_to_patient(cnv_mat, sample_map$Patient)

if (!is.null(cnv_patient)) {
  
  common_pat <- intersect(colnames(cnv_patient), sample_map$Patient)
  cnv_all <- cnv_patient[, common_pat, drop = FALSE]
  
  cnv_type <- ifelse(is_gistic_cnv(cnv_all), "GISTIC_discrete", "continuous_gene_level")
  message("CNV type detected: ", cnv_type)
  
  cnv_alter <- make_cnv_binary(cnv_all, "alter")
  cnv_amp <- make_cnv_binary(cnv_all, "amp")
  cnv_del <- make_cnv_binary(cnv_all, "del")
  
  ddr_cnv_genes <- intersect(DDR_GENES, rownames(cnv_all))
  
  tmp_cnv <- data.frame(
    Patient = common_pat,
    AllGene_CNA_Fraction = colMeans(cnv_alter, na.rm = TRUE),
    AllGene_Amp_Fraction = colMeans(cnv_amp, na.rm = TRUE),
    AllGene_Del_Fraction = colMeans(cnv_del, na.rm = TRUE),
    DDRGene_CNA_Fraction = if (length(ddr_cnv_genes) >= 5) {
      colMeans(cnv_alter[ddr_cnv_genes, , drop = FALSE], na.rm = TRUE)
    } else {
      NA_real_
    },
    DDRGene_Amp_Fraction = if (length(ddr_cnv_genes) >= 5) {
      colMeans(cnv_amp[ddr_cnv_genes, , drop = FALSE], na.rm = TRUE)
    } else {
      NA_real_
    },
    DDRGene_Del_Fraction = if (length(ddr_cnv_genes) >= 5) {
      colMeans(cnv_del[ddr_cnv_genes, , drop = FALSE], na.rm = TRUE)
    } else {
      NA_real_
    }
  )
  
  cnv_features <- sample_map %>%
    dplyr::select(Patient, Sample, MO_DDRscore_group) %>%
    dplyr::left_join(tmp_cnv, by = "Patient")
}

save_csv(
  cnv_features,
  file.path(SUPP_TAB_DIR, "FigS3_gene_level_CNA_fraction_table.csv")
)

cnv_long <- cnv_features %>%
  dplyr::filter(!is.na(MO_DDRscore_group)) %>%
  tidyr::pivot_longer(
    cols = c(
      AllGene_CNA_Fraction,
      AllGene_Amp_Fraction,
      AllGene_Del_Fraction,
      DDRGene_CNA_Fraction,
      DDRGene_Amp_Fraction,
      DDRGene_Del_Fraction
    ),
    names_to = "Feature",
    values_to = "Value"
  ) %>%
  dplyr::filter(is.finite(Value))

if (nrow(cnv_long) > 0) {
  
  cnv_long$Feature <- factor(
    cnv_long$Feature,
    levels = c(
      "AllGene_CNA_Fraction",
      "AllGene_Amp_Fraction",
      "AllGene_Del_Fraction",
      "DDRGene_CNA_Fraction",
      "DDRGene_Amp_Fraction",
      "DDRGene_Del_Fraction"
    ),
    labels = c(
      "All-gene CNA altered fraction",
      "All-gene amplification fraction",
      "All-gene deletion fraction",
      "DDR-gene CNA altered fraction",
      "DDR-gene amplification fraction",
      "DDR-gene deletion fraction"
    )
  )
  
  cnv_stat <- cnv_long %>%
    dplyr::group_by(Feature) %>%
    dplyr::summarise(
      Median_Low = median(Value[MO_DDRscore_group == "Low"], na.rm = TRUE),
      Median_High = median(Value[MO_DDRscore_group == "High"], na.rm = TRUE),
      Diff = Median_High - Median_Low,
      P = tryCatch(wilcox.test(Value ~ MO_DDRscore_group)$p.value, error = function(e) NA_real_),
      y_pos = max(Value, na.rm = TRUE) * 1.08,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      FDR = p.adjust(P, "BH"),
      Label = sig_label(FDR),
      x_pos = 1.5
    )
  
  save_csv(
    cnv_stat,
    file.path(SUPP_TAB_DIR, "FigS3_gene_level_CNA_fraction_stats.csv")
  )
  
  p_cnv <- ggplot(
    cnv_long,
    aes(x = MO_DDRscore_group, y = Value, fill = MO_DDRscore_group)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.62, linewidth = 0.35) +
    geom_jitter(width = 0.14, size = 0.5, alpha = 0.30) +
    geom_text(
      data = cnv_stat,
      aes(x = x_pos, y = y_pos, label = Label),
      inherit.aes = FALSE,
      size = 4,
      fontface = "bold"
    ) +
    facet_wrap(~Feature, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c(Low = LOW_COL, High = HIGH_COL)) +
    theme_pub(11.5) +
    theme(legend.position = "none", strip.text = element_text(face = "bold")) +
    labs(
      x = NULL,
      y = "Fraction of altered genes",
      title = "Exploratory gene-level CNA altered fraction"
    )
  
  safe_ggsave(
    file.path(SUPP_DIR, "FigS3_gene_level_CNA_fraction_boxplot.pdf"),
    p_cnv,
    9,
    5.5
  )
}

############################################################
## 17. Done
############################################################

cat("\nFig3 multi-omics + GDSC/DepMap pipeline finished.\n")
cat("Main outputs:\n")
cat(" - Fig3A1_DEG_volcano.pdf\n")
cat(" - Fig3A2_top_DEG_heatmap.pdf\n")
cat(" - Fig3B_focused_MSigDB_Reactome_GSEA.pdf\n")
cat(" - Fig3C_focused_pathway_ssGSEA_activity.pdf\n")
cat(" - Fig3D1_DDR_only_mutation_oncoplot.pdf\n")
cat(" - Fig3D2_DDR_mutation_burden.pdf\n")
cat(" - Fig3E1_ESTIMATE_tumor_microenvironment.pdf, if ESTIMATE available\n")
cat(" - Fig3E2_CIBERSORT_LM22_immune_cell_infiltration.pdf, if official CIBERSORT files exist\n")
cat(" - Fig3F1_GDSC_drug_sensitivity_boxplot.pdf, if DepMap/GDSC files exist\n")
cat(" - Fig3F2_GDSC_score_IC50_correlation_dotplot.pdf, if DepMap/GDSC files exist\n")
cat(" - Fig3G_MO_DDRscore_KM_survival.pdf\n")

cat("\nSupplementary outputs:\n")
cat(" - FigS1_checkpoint_expression_boxplot.pdf\n")
cat(" - FigS2_clinical_feature_distribution.pdf\n")
cat(" - FigS3_gene_level_CNA_fraction_boxplot.pdf\n")

cat("\nOutput directories:\n")
cat(" Main figures: ", FIG3_DIR, "\n")
cat(" Main tables:  ", TAB3_DIR, "\n")
cat(" Supp figures: ", SUPP_DIR, "\n")
cat(" Supp tables:  ", SUPP_TAB_DIR, "\n")