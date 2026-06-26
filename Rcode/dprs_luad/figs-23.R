############################################################
# MO-DDRscore clean pipeline
# Output only PDF figures to:
# D:/R_workspace/评分/AD_DDR_project/res
#
# Main outputs:
# Fig2: MO-DDRscore construction / ablation / SPIDR validation
# Fig3: DEG / Hallmark / mutation / CNV / immune / checkpoint / survival
#
# Fixed issues:
# 1) Use MO-DDRscore / MO-DDRweight naming
# 2) Raw MO-DDRscore used for statistics
# 3) 0-100 scaled MO-DDRscore used only for visualization
# 4) Fig2E uses tie-safe ranked dot plot, not histogram
# 5) High/Low rule: Low <= median cutoff; High > median cutoff
# 6) CNV burden calculated genome-wide from all genes in CNV matrix
# 7) PDF only, no PNG, no CSV outputs
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260510)

############################################################
# 0. Paths and parameters
############################################################

PROJECT_DIR <- "D:/R_workspace/评分/AD_DDR_project"
DATA_DIR <- file.path(PROJECT_DIR, "00_data")
RES_DIR <- file.path(PROJECT_DIR, "res")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

FOCUS_CANCER <- "LUAD"
TCGA_DIR <- "D:/R/R_workspace/梁老师文件/TCGA"

TCGA_EXPR_FILE <- file.path(
  TCGA_DIR,
  "mRNA_exp_TPM_only_TCGA/mRNA_exp_TPM_only_TCGA",
  paste0("TCGA-", FOCUS_CANCER, ".gene_expression_TPM.tsv")
)

TCGA_CLIN_FILE <- file.path(
  TCGA_DIR,
  "clinical",
  paste0("TCGA.", FOCUS_CANCER, ".sampleMap"),
  paste0(FOCUS_CANCER, "_clinicalMatrix")
)

TCGA_MUT_FILE <- file.path(TCGA_DIR, "mutation", paste0(FOCUS_CANCER, ".txt"))
TCGA_CNV_FILE <- file.path(TCGA_DIR, "cnv", paste0(FOCUS_CANCER, ".txt"))

DDR_GENE_FILE <- file.path(DATA_DIR, "DDR_236_genes.csv")
PPI_FILE <- file.path(DATA_DIR, "string_interactions_short_0.7.tsv")
SPIDR_FILE <- file.path(DATA_DIR, "spidr_supp_table3_clean.csv")

DEG_LOGFC <- 1.2
DEG_FDR <- 0.05
MAX_HEATMAP_GENES <- 50

LOW_COL <- "#4DBBD5"
HIGH_COL <- "#E64B35"
BLUE_COL <- "#2166AC"
RED_COL <- "#B2182B"

############################################################
# 1. Packages and helper functions
############################################################

cran_pkgs <- c(
  "data.table", "dplyr", "tidyr", "ggplot2",
  "survival", "pheatmap"
)

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
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

safe_pdf <- function(file, width = 7, height = 5) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
}

safe_ggsave <- function(file, p, width = 7, height = 5) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = file,
    plot = p,
    width = width,
    height = height,
    device = "pdf",
    useDingbats = FALSE
  )
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
  x[x %in% c(
    "", "NA", "NaN", "[Not Available]", "[Not Applicable]",
    "not reported", "Not Reported", "null", "NULL", "--"
  )] <- NA
  x
}

to_num <- function(x) suppressWarnings(as.numeric(clean_na(x)))

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
  m <- as.matrix(m)
  storage.mode(m) <- "numeric"
  mu <- rowMeans(m, na.rm = TRUE)
  sdv <- apply(m, 1, sd, na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- 1
  out <- sweep(sweep(m, 1, mu, "-"), 1, sdv, "/")
  out[!is.finite(out)] <- 0
  out
}

z_cols <- function(m) {
  m <- as.matrix(m)
  storage.mode(m) <- "numeric"
  mu <- colMeans(m, na.rm = TRUE)
  sdv <- apply(m, 2, sd, na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- 1
  out <- sweep(sweep(m, 2, mu, "-"), 2, sdv, "/")
  out[!is.finite(out)] <- 0
  out
}

safe_cor <- function(x, y, method = "spearman") {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10) return(c(cor = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method))
  c(cor = unname(ct$estimate), p = ct$p.value)
}

calc_cindex <- function(time, status, score) {
  df <- data.frame(
    time = as.numeric(time),
    status = as.numeric(status),
    score = as.numeric(score)
  )
  df <- df[is.finite(df$time) & is.finite(df$status) & is.finite(df$score), ]
  if (nrow(df) < 20 || sum(df$status == 1, na.rm = TRUE) < 5) return(NA_real_)
  suppressWarnings(
    tryCatch(
      survival::concordance(Surv(time, status) ~ score, data = df)$concordance,
      error = function(e) NA_real_
    )
  )
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
  
  list(
    weights = w,
    normalized = mat_norm,
    entropy = ent,
    divergence = div
  )
}

############################################################
# 2. DDR genes
############################################################

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
  tmp_gene <- safe_fread(DDR_GENE_FILE)
  gc <- find_col(tmp_gene, c("^Gene$", "gene", "symbol"))
  if (is.na(gc)) gc <- colnames(tmp_gene)[1]
  DDR_GENES <- unique(na.omit(clean_gene(tmp_gene[[gc]])))
} else {
  DDR_GENES <- unique(na.omit(clean_gene(fallback_genes)))
}

cat("DDR genes loaded:", length(DDR_GENES), "\n")

############################################################
# 3. Read TCGA data
############################################################

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
  mat
}

read_clin <- function(file) {
  x <- safe_fread(file)
  
  pc <- find_col(x, c("^_PATIENT$", "bcr_patient_barcode", "patient_id", "Patient", "sampleID"))
  if (is.na(pc)) pc <- colnames(x)[1]
  
  vc <- find_col(x, c("^vital_status$", "vital.status", "death", "dead"))
  dc <- find_col(x, c("^days_to_death$", "days.to.death"))
  fc <- find_col(x, c("^days_to_last_followup$", "days_to_last_follow_up", "days.to.last.followup", "days_to_followup"))
  
  if (is.na(vc) || (is.na(dc) && is.na(fc))) {
    stop("Cannot infer clinical survival columns.")
  }
  
  vital <- toupper(clean_na(x[[vc]]))
  death <- if (!is.na(dc)) to_num(x[[dc]]) else rep(NA_real_, nrow(x))
  follow <- if (!is.na(fc)) to_num(x[[fc]]) else rep(NA_real_, nrow(x))
  
  status <- ifelse(
    vital %in% c("DECEASED", "DEAD", "DIED"),
    1,
    ifelse(vital %in% c("LIVING", "ALIVE"), 0, NA_real_)
  )
  
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
  
  out <- out[
    !is.na(out$Patient) &
      out$Patient != "" &
      is.finite(out$time) &
      out$time > 0 &
      out$status %in% c(0, 1),
  ]
  
  out %>%
    dplyr::arrange(Patient, dplyr::desc(status), dplyr::desc(time)) %>%
    dplyr::distinct(Patient, .keep_all = TRUE)
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
  mat
}

cat("\nLoading TCGA data...\n")

tcga_expr <- read_expr(TCGA_EXPR_FILE)
tumor_samples <- colnames(tcga_expr)[sample_class(colnames(tcga_expr)) == "Tumor"]
normal_samples <- colnames(tcga_expr)[sample_class(colnames(tcga_expr)) == "Normal"]
tcga_clin <- read_clin(TCGA_CLIN_FILE)
mut <- read_mut(TCGA_MUT_FILE)
cnv_mat <- read_cnv(TCGA_CNV_FILE)

cat("Expression:", nrow(tcga_expr), "genes x", ncol(tcga_expr), "samples\n")
cat("Tumor:", length(tumor_samples), " Normal:", length(normal_samples), "\n")
cat("Clinical:", nrow(tcga_clin), " Events:", sum(tcga_clin$status == 1), "\n")

############################################################
# 4. MO-DDRweight construction
############################################################

cat("\nConstructing MO-DDRweight...\n")

ddr_genes <- intersect(DDR_GENES, rownames(tcga_expr))
if (length(ddr_genes) < 20) stop("Too few DDR genes in expression matrix.")

############################################################
# E1: Expression dysregulation
############################################################

expr_log <- log2(tcga_expr[ddr_genes, , drop = FALSE] + 1)

if (length(normal_samples) >= 5) {
  tumor_mean <- rowMeans(expr_log[, tumor_samples, drop = FALSE], na.rm = TRUE)
  normal_mean <- rowMeans(expr_log[, normal_samples, drop = FALSE], na.rm = TRUE)
  logFC <- tumor_mean - normal_mean
  
  pval <- sapply(ddr_genes, function(g) {
    tryCatch(
      wilcox.test(expr_log[g, tumor_samples], expr_log[g, normal_samples])$p.value,
      error = function(e) NA_real_
    )
  })
  
  fdr <- p.adjust(pval, "BH")
  ExprDysregulation <- abs(logFC) * (-log10(fdr + 1e-300))
} else {
  logFC <- rep(NA_real_, length(ddr_genes))
  pval <- fdr <- rep(NA_real_, length(ddr_genes))
  tumor_expr <- expr_log[, tumor_samples, drop = FALSE]
  ExprDysregulation <- minmax01(rowMeans(tumor_expr, na.rm = TRUE)) *
    minmax01(apply(tumor_expr, 1, sd, na.rm = TRUE))
}

expr_ev <- data.frame(
  Gene = ddr_genes,
  TumorNormal_logFC = logFC,
  TumorNormal_P = pval,
  TumorNormal_FDR = fdr,
  ExprDysregulation = ExprDysregulation
)

############################################################
# E2: Genomic alteration
############################################################

gen_ev <- data.frame(
  Gene = ddr_genes,
  MutationFreq = 0,
  AmpFreq = 0,
  DelFreq = 0,
  CNVAlterFreq = 0
)

if (!is.null(mut)) {
  pats <- unique(patient_id(tumor_samples))
  
  mf <- mut %>%
    dplyr::filter(Hugo_Symbol %in% ddr_genes, Patient %in% pats) %>%
    dplyr::distinct(Hugo_Symbol, Patient) %>%
    dplyr::group_by(Hugo_Symbol) %>%
    dplyr::summarise(
      MutationFreq = dplyr::n_distinct(Patient) / length(pats),
      .groups = "drop"
    )
  
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

############################################################
# E3: Network context
############################################################

read_ppi <- function(file) {
  if (!file.exists(file)) return(NULL)
  p <- tryCatch(safe_fread(file), error = function(e) NULL)
  if (is.null(p) || nrow(p) == 0) return(NULL)
  
  g1 <- find_col(p, c("protein1", "gene1", "Gene1", "preferredName_A", "node1", "from"))
  g2 <- find_col(p, c("protein2", "gene2", "Gene2", "preferredName_B", "node2", "to"))
  sc <- find_col(p, c("combined_score", "score", "confidence", "weight"))
  
  if (is.na(g1) || is.na(g2)) return(NULL)
  
  out <- data.frame(
    GeneA = clean_gene(p[[g1]]),
    GeneB = clean_gene(p[[g2]]),
    PPI_score = if (!is.na(sc)) as.numeric(p[[sc]]) else 1
  )
  
  out <- out[
    !is.na(out$GeneA) &
      !is.na(out$GeneB) &
      out$GeneA != out$GeneB,
  ]
  
  out$PPI_score[!is.finite(out$PPI_score)] <- 1
  out$PPI_score <- minmax01(out$PPI_score)
  out
}

ppi <- read_ppi(PPI_FILE)

net_ev <- data.frame(
  Gene = ddr_genes,
  PPI_Degree = 0,
  PPI_WeightedDegree = 0,
  NetworkContext = 0
)

if (!is.null(ppi)) {
  ppi_ddr <- ppi %>%
    dplyr::filter(GeneA %in% ddr_genes, GeneB %in% ddr_genes)
  
  if (nrow(ppi_ddr) > 0) {
    dl <- dplyr::bind_rows(
      ppi_ddr %>% dplyr::select(Gene = GeneA, PPI_score),
      ppi_ddr %>% dplyr::select(Gene = GeneB, PPI_score)
    )
    
    ds <- dl %>%
      dplyr::group_by(Gene) %>%
      dplyr::summarise(
        PPI_Degree = dplyr::n(),
        PPI_WeightedDegree = sum(PPI_score, na.rm = TRUE),
        .groups = "drop"
      )
    
    net_ev$PPI_Degree <- ds$PPI_Degree[match(net_ev$Gene, ds$Gene)]
    net_ev$PPI_WeightedDegree <- ds$PPI_WeightedDegree[match(net_ev$Gene, ds$Gene)]
    net_ev$PPI_Degree[is.na(net_ev$PPI_Degree)] <- 0
    net_ev$PPI_WeightedDegree[is.na(net_ev$PPI_WeightedDegree)] <- 0
    
    net_ev$NetworkContext <- minmax01(log1p(net_ev$PPI_Degree)) +
      minmax01(net_ev$PPI_WeightedDegree)
    net_ev$NetworkContext <- minmax01(net_ev$NetworkContext)
  }
}

############################################################
# E4: Co-expression coherence
############################################################

tumor_ddr <- log2(tcga_expr[ddr_genes, tumor_samples, drop = FALSE] + 1)
tumor_ddr_z <- t(z_rows(tumor_ddr))

cc <- suppressWarnings(
  cor(tumor_ddr_z, method = "spearman", use = "pairwise.complete.obs")
)

diag(cc) <- NA

mean_abs <- rowMeans(abs(cc), na.rm = TRUE)

top_abs <- apply(abs(cc), 1, function(x) {
  x <- sort(x[is.finite(x)], decreasing = TRUE)
  if (length(x) == 0) return(0)
  mean(head(x, max(1, ceiling(length(x) * 0.10))), na.rm = TRUE)
})

co_ev <- data.frame(
  Gene = ddr_genes,
  MeanAbsDDRCoexpression = mean_abs,
  Top10AbsDDRCoexpression = top_abs
)

co_ev$CoexpressionCoherence <- minmax01(co_ev$MeanAbsDDRCoexpression) +
  minmax01(co_ev$Top10AbsDDRCoexpression)
co_ev$CoexpressionCoherence <- minmax01(co_ev$CoexpressionCoherence)

############################################################
# Evidence fusion
############################################################

evidence <- expr_ev %>%
  dplyr::left_join(gen_ev, by = "Gene") %>%
  dplyr::left_join(net_ev, by = "Gene") %>%
  dplyr::left_join(co_ev, by = "Gene")

evidence_cols <- c(
  "ExprDysregulation",
  "GenomicAlteration",
  "NetworkContext",
  "CoexpressionCoherence"
)

for (v in evidence_cols) {
  evidence[[v]][!is.finite(evidence[[v]])] <- 0
}

ef <- entropy_fusion(evidence[, c("Gene", evidence_cols)], id_col = "Gene")
layer_weights <- ef$weights

norm_ev <- as.data.frame(ef$normalized, check.names = FALSE)
norm_ev$Gene <- rownames(ef$normalized)

MO_DDRweight <- evidence %>%
  dplyr::left_join(norm_ev, by = "Gene", suffix = c("", "_Norm"))

MO_DDRweight$Raw_MO_DDRweight <- 0

for (v in evidence_cols) {
  MO_DDRweight$Raw_MO_DDRweight <- MO_DDRweight$Raw_MO_DDRweight +
    MO_DDRweight[[paste0(v, "_Norm")]] * layer_weights[v]
}

MO_DDRweight$MO_DDRweight <- minmax01(MO_DDRweight$Raw_MO_DDRweight)

MO_DDRweight <- MO_DDRweight[
  order(MO_DDRweight$MO_DDRweight, decreasing = TRUE),
]

layer_df <- data.frame(
  EvidenceLayer = names(layer_weights),
  AdaptiveWeight = as.numeric(layer_weights),
  Entropy = as.numeric(ef$entropy[names(layer_weights)]),
  Divergence = as.numeric(ef$divergence[names(layer_weights)])
)

print(layer_df)

############################################################
# 5. Fig2A-D: MO-DDRweight construction figures
############################################################

cat("\nDrawing Fig2...\n")

############################################################
# Fig2A: evidence heatmap
############################################################

top_evi <- head(MO_DDRweight$Gene, min(40, nrow(MO_DDRweight)))

ev_mat <- as.matrix(MO_DDRweight[match(top_evi, MO_DDRweight$Gene), evidence_cols])
ev_mat <- apply(ev_mat, 2, minmax01)
rownames(ev_mat) <- top_evi

if (has_pkg("pheatmap")) {
  safe_pdf(file.path(RES_DIR, "Fig2A_gene_level_multiomics_evidence_heatmap.pdf"), 6.5, 8)
  pheatmap::pheatmap(
    ev_mat,
    cluster_cols = FALSE,
    color = colorRampPalette(c("#F7FBFF", "#6BAED6", "#08306B"))(100),
    fontsize_row = 7,
    main = "Gene-level multi-omics evidence"
  )
  dev.off()
}

############################################################
# Fig2B: entropy layer weights
############################################################

p_layer <- ggplot(
  layer_df,
  aes(x = reorder(EvidenceLayer, AdaptiveWeight), y = AdaptiveWeight)
) +
  geom_col(fill = "#4A5568", width = 0.72) +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    x = NULL,
    y = "Entropy-derived adaptive weight",
    title = "Adaptive evidence-layer contribution"
  )

safe_ggsave(
  file.path(RES_DIR, "Fig2B_entropy_adaptive_evidence_weights.pdf"),
  p_layer,
  5.5,
  4
)

############################################################
# Fig2C: MO-DDRweight distribution
############################################################

p_dist <- ggplot(MO_DDRweight, aes(x = MO_DDRweight)) +
  geom_histogram(bins = 30, fill = "#4575B4", color = "white") +
  theme_bw(base_size = 12) +
  labs(
    x = "MO-DDRweight",
    y = "Gene count",
    title = "MO-DDRweight distribution"
  )

safe_ggsave(
  file.path(RES_DIR, "Fig2C_MO_DDRweight_distribution.pdf"),
  p_dist,
  5,
  4
)

############################################################
# Fig2D: top weighted genes
############################################################

p_top <- ggplot(
  head(MO_DDRweight, 30),
  aes(x = reorder(Gene, MO_DDRweight), y = MO_DDRweight)
) +
  geom_segment(aes(xend = Gene, y = 0, yend = MO_DDRweight), color = "grey60") +
  geom_point(size = 2.5, color = HIGH_COL) +
  coord_flip() +
  theme_bw(base_size = 12) +
  labs(
    x = NULL,
    y = "MO-DDRweight",
    title = "Top MO-DDRweighted genes"
  )

safe_ggsave(
  file.path(RES_DIR, "Fig2D_top_MO_DDRweight_genes.pdf"),
  p_top,
  6,
  6
)

############################################################
# 6. Sample-level MO-DDRscore and baseline scores
############################################################

compute_score <- function(expr, samples, wt, score_name, ref_samples = tumor_samples) {
  genes <- intersect(wt$Gene, rownames(expr))
  
  expr_log <- log2(expr[genes, samples, drop = FALSE] + 1)
  ref_log <- log2(expr[genes, ref_samples, drop = FALSE] + 1)
  
  mu <- rowMeans(ref_log, na.rm = TRUE)
  sdv <- apply(ref_log, 1, sd, na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- 1
  
  z <- sweep(sweep(expr_log, 1, mu, "-"), 1, sdv, "/")
  z[!is.finite(z)] <- 0
  
  w <- wt$MO_DDRweight[match(genes, wt$Gene)]
  w[is.na(w)] <- 0
  
  sc <- colSums(z * w, na.rm = TRUE) / sum(abs(w), na.rm = TRUE)
  
  out <- data.frame(
    Sample = names(sc),
    Patient = patient_id(names(sc)),
    SampleType = sample_type(names(sc)),
    SampleClass = sample_class(names(sc))
  )
  
  out[[score_name]] <- as.numeric(sc)
  out
}

score_samples <- unique(c(tumor_samples, normal_samples))

mo_score <- compute_score(
  tcga_expr,
  score_samples,
  MO_DDRweight[, c("Gene", "MO_DDRweight")],
  "MO_DDRscore_raw"
)

eq <- compute_score(
  tcga_expr,
  score_samples,
  data.frame(Gene = ddr_genes, MO_DDRweight = 1),
  "Equal_DDRscore"
)

exprwt <- MO_DDRweight[, c("Gene", "ExprDysregulation_Norm")]
colnames(exprwt) <- c("Gene", "MO_DDRweight")

expronly <- compute_score(
  tcga_expr,
  score_samples,
  exprwt,
  "ExpressionOnly_DDRscore"
)

nonet <- data.frame(
  Gene = MO_DDRweight$Gene,
  MO_DDRweight = rowMeans(
    MO_DDRweight[, c("ExprDysregulation_Norm", "GenomicAlteration_Norm")],
    na.rm = TRUE
  )
)

nonets <- compute_score(
  tcga_expr,
  score_samples,
  nonet,
  "NoNetwork_DDRscore"
)

mo_score <- mo_score %>%
  dplyr::left_join(eq[, c("Sample", "Equal_DDRscore")], by = "Sample") %>%
  dplyr::left_join(expronly[, c("Sample", "ExpressionOnly_DDRscore")], by = "Sample") %>%
  dplyr::left_join(nonets[, c("Sample", "NoNetwork_DDRscore")], by = "Sample")

############################################################
# Tie-safe grouping and scaled score
############################################################

tumor_idx <- mo_score$SampleClass == "Tumor" & is.finite(mo_score$MO_DDRscore_raw)

score_min <- min(mo_score$MO_DDRscore_raw[tumor_idx], na.rm = TRUE)
score_max <- max(mo_score$MO_DDRscore_raw[tumor_idx], na.rm = TRUE)

mo_score$MO_DDRscore_0_100 <- NA_real_

if (is.finite(score_min) && is.finite(score_max) && score_max > score_min) {
  mo_score$MO_DDRscore_0_100 <- 100 *
    (mo_score$MO_DDRscore_raw - score_min) /
    (score_max - score_min)
}

mo_score$MO_DDRscore_0_100 <- pmax(pmin(mo_score$MO_DDRscore_0_100, 100), 0)

cut_raw <- median(mo_score$MO_DDRscore_raw[tumor_idx], na.rm = TRUE)
cut_scaled <- 100 * (cut_raw - score_min) / (score_max - score_min)

mo_score$MO_DDRscore_group <- NA_character_

mo_score$MO_DDRscore_group[tumor_idx] <- ifelse(
  mo_score$MO_DDRscore_raw[tumor_idx] > cut_raw,
  "High",
  "Low"
)

mo_score$MO_DDRscore_group <- factor(
  mo_score$MO_DDRscore_group,
  levels = c("Low", "High")
)

############################################################
# Fig2E: ranked sample-level MO-DDRscore
############################################################

plot_score <- mo_score %>%
  dplyr::filter(
    SampleClass == "Tumor",
    is.finite(MO_DDRscore_raw),
    !is.na(MO_DDRscore_group)
  ) %>%
  dplyr::arrange(MO_DDRscore_raw)

plot_score$Rank <- seq_len(nrow(plot_score))
gap_pos <- max(plot_score$Rank[plot_score$MO_DDRscore_group == "Low"], na.rm = TRUE)

p_rank <- ggplot(
  plot_score,
  aes(x = Rank, y = MO_DDRscore_0_100, color = MO_DDRscore_group)
) +
  geom_point(size = 1.2, alpha = 0.85) +
  geom_vline(
    xintercept = gap_pos + 0.5,
    linetype = 2,
    color = "grey40",
    linewidth = 0.5
  ) +
  scale_color_manual(values = c(Low = LOW_COL, High = HIGH_COL)) +
  theme_bw(base_size = 12) +
  labs(
    x = "LUAD tumor samples ranked by raw MO-DDRscore",
    y = "Scaled MO-DDRscore (0-100)",
    color = "Group",
    title = "Sample-level MO-DDRscore stratification"
  )

safe_ggsave(
  file.path(RES_DIR, "Fig2E_MO_DDRscore_ranked_stratification.pdf"),
  p_rank,
  6,
  4
)

############################################################
# 7. Fig2F: baseline / ablation comparison
############################################################

tumor_mo <- mo_score %>%
  dplyr::filter(SampleClass == "Tumor") %>%
  dplyr::arrange(Patient, SampleType) %>%
  dplyr::distinct(Patient, .keep_all = TRUE)

mo_score_clin <- tumor_mo %>%
  dplyr::inner_join(tcga_clin, by = "Patient")

tmb_df <- NULL

if (!is.null(mut)) {
  tmb_df <- mut %>%
    dplyr::filter(!is.na(Patient), !is.na(Hugo_Symbol)) %>%
    dplyr::group_by(Patient) %>%
    dplyr::summarise(
      MutationCount = dplyr::n(),
      MutatedGeneCount = dplyr::n_distinct(Hugo_Symbol),
      DDR_MutatedGeneCount = dplyr::n_distinct(Hugo_Symbol[Hugo_Symbol %in% DDR_GENES]),
      .groups = "drop"
    )
}

cnv_df <- NULL

if (!is.null(cnv_mat)) {
  cnv_df <- data.frame(
    Sample = colnames(cnv_mat),
    Patient = patient_id(colnames(cnv_mat)),
    CNVBurden = colMeans(abs(cnv_mat) > 0.2, na.rm = TRUE),
    AmpBurden = colMeans(cnv_mat > 0.2, na.rm = TRUE),
    DelBurden = colMeans(cnv_mat < -0.2, na.rm = TRUE)
  ) %>%
    dplyr::group_by(Patient) %>%
    dplyr::summarise(
      CNVBurden = mean(CNVBurden, na.rm = TRUE),
      AmpBurden = mean(AmpBurden, na.rm = TRUE),
      DelBurden = mean(DelBurden, na.rm = TRUE),
      .groups = "drop"
    )
}

multiomics_score <- mo_score_clin

if (!is.null(tmb_df)) {
  multiomics_score <- multiomics_score %>% dplyr::left_join(tmb_df, by = "Patient")
}

if (!is.null(cnv_df)) {
  multiomics_score <- multiomics_score %>% dplyr::left_join(cnv_df, by = "Patient")
}

for (v in c("MutationCount", "MutatedGeneCount", "DDR_MutatedGeneCount",
            "CNVBurden", "AmpBurden", "DelBurden")) {
  if (!v %in% colnames(multiomics_score)) multiomics_score[[v]] <- NA_real_
}

score_sets <- list(
  DNA_Repair = c(
    "BRCA1","BRCA2","RAD51","RAD50","MRE11","NBN","ATM","ATR",
    "CHEK1","CHEK2","ERCC1","ERCC2","MSH2","MSH6","MLH1","PMS2",
    "XRCC1","XRCC5","XRCC6","PRKDC","PARP1","LIG1","LIG3","LIG4",
    "FANCA","FANCD2","FANCI"
  ),
  Cell_Cycle_G2M = c(
    "AURKA","AURKB","BUB1","BUB1B","CCNB1","CCNB2","CDC20",
    "CDC25A","CDC25B","CDC25C","CDK1","MKI67","PLK1","TOP2A"
  ),
  DNA_Replication = c(
    "MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","PCNA","POLE",
    "POLD1","RFC1","RFC2","RPA1","RPA2","CDC6","CDC45"
  ),
  Checkpoint = c(
    "ATM","ATR","CHEK1","CHEK2","TP53","TP53BP1",
    "CLSPN","TOPBP1","HUS1","RAD9A","RAD17"
  )
)

expr_tumor_log <- log2(tcga_expr[, tumor_samples, drop = FALSE] + 1)
expr_tumor_z <- z_rows(expr_tumor_log)

ps <- list()

for (nm in names(score_sets)) {
  gs <- intersect(clean_gene(score_sets[[nm]]), rownames(expr_tumor_z))
  if (length(gs) >= 3) {
    ps[[nm]] <- colMeans(expr_tumor_z[gs, , drop = FALSE], na.rm = TRUE)
  }
}

path_score <- as.data.frame(ps, check.names = FALSE)
path_score$Sample <- rownames(path_score)

path_score <- path_score %>%
  dplyr::left_join(
    mo_score[, c(
      "Sample", "Patient", "MO_DDRscore_raw", "Equal_DDRscore",
      "ExpressionOnly_DDRscore", "NoNetwork_DDRscore", "MO_DDRscore_group"
    )],
    by = "Sample"
  ) %>%
  dplyr::left_join(
    multiomics_score[, c("Patient", "MutationCount", "CNVBurden", "AmpBurden", "DelBurden")],
    by = "Patient"
  )

bio_feats <- setdiff(
  colnames(path_score),
  c(
    "Sample", "Patient", "MO_DDRscore_raw", "Equal_DDRscore",
    "ExpressionOnly_DDRscore", "NoNetwork_DDRscore", "MO_DDRscore_group"
  )
)

bio_feats <- bio_feats[
  sapply(bio_feats, function(v) is.numeric(path_score[[v]]) || is.integer(path_score[[v]]))
]

score_cols <- c(
  "MO_DDRscore_raw",
  "Equal_DDRscore",
  "ExpressionOnly_DDRscore",
  "NoNetwork_DDRscore"
)

score_label <- c(
  MO_DDRscore_raw = "Full MO-DDRscore",
  Equal_DDRscore = "Equal-weight",
  ExpressionOnly_DDRscore = "Expression-only",
  NoNetwork_DDRscore = "No-network"
)

base_stats <- dplyr::bind_rows(lapply(score_cols, function(sc) {
  dplyr::bind_rows(lapply(bio_feats, function(ft) {
    cr <- safe_cor(path_score[[sc]], path_score[[ft]])
    data.frame(
      ScoreType = score_label[[sc]],
      Feature = ft,
      Correlation = cr["cor"],
      P = cr["p"],
      AbsCorrelation = abs(cr["cor"])
    )
  }))
}))

base_stats$FDR <- p.adjust(base_stats$P, "BH")

p_base <- ggplot(
  base_stats,
  aes(x = reorder(Feature, AbsCorrelation), y = Correlation, fill = ScoreType)
) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  geom_col(position = "dodge", width = 0.75) +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(
    x = NULL,
    y = "Spearman correlation",
    fill = "Score",
    title = "Ablation analysis of MO-DDRscore"
  )

safe_ggsave(
  file.path(RES_DIR, "Fig2F_MO_DDRscore_ablation_comparison.pdf"),
  p_base,
  8,
  5.8
)

############################################################
# 8. Fig2G: SPIDR external validation
############################################################

read_spidr <- function(file) {
  if (!file.exists(file)) return(NULL)
  
  sp <- tryCatch(safe_fread(file), error = function(e) NULL)
  if (is.null(sp)) return(NULL)
  
  g1 <- if ("gene_a_base" %in% colnames(sp)) {
    "gene_a_base"
  } else if ("gene_a_raw" %in% colnames(sp)) {
    "gene_a_raw"
  } else {
    find_col(sp, c("^GeneA$", "^gene_a$", "^Gene1$", "^gene1$"))
  }
  
  g2 <- if ("gene_b_base" %in% colnames(sp)) {
    "gene_b_base"
  } else if ("gene_b_raw" %in% colnames(sp)) {
    "gene_b_raw"
  } else {
    find_col(sp, c("^GeneB$", "^gene_b$", "^Gene2$", "^gene2$"))
  }
  
  if (is.na(g1) || is.na(g2) || !"gemini_sensitive" %in% colnames(sp)) {
    return(NULL)
  }
  
  lab <- if ("sl_label_le_1_0" %in% colnames(sp)) {
    "sl_label_le_1_0"
  } else if ("strong_sl_label_le_1_5" %in% colnames(sp)) {
    "strong_sl_label_le_1_5"
  } else {
    NA
  }
  
  out <- data.frame(
    GeneA = clean_gene(sp[[g1]]),
    GeneB = clean_gene(sp[[g2]]),
    GeminiScore = as.numeric(sp$gemini_sensitive)
  )
  
  out$SL_Label <- if (!is.na(lab)) {
    as.integer(as.character(sp[[lab]]))
  } else {
    ifelse(out$GeminiScore <= -1, 1, 0)
  }
  
  out <- out[
    !is.na(out$GeneA) &
      !is.na(out$GeneB) &
      out$GeneA != out$GeneB &
      is.finite(out$GeminiScore),
  ]
  
  out$PairKey <- ifelse(
    out$GeneA < out$GeneB,
    paste(out$GeneA, out$GeneB, sep = "__"),
    paste(out$GeneB, out$GeneA, sep = "__")
  )
  
  out <- out %>%
    dplyr::group_by(PairKey) %>%
    dplyr::summarise(
      GeneA = dplyr::first(GeneA),
      GeneB = dplyr::first(GeneB),
      GeminiScore = min(GeminiScore, na.rm = TRUE),
      SL_Label = max(SL_Label, na.rm = TRUE),
      .groups = "drop"
    )
  
  out$SLStrength <- minmax01(pmax(0, -out$GeminiScore))
  out
}

spidr <- read_spidr(SPIDR_FILE)

if (!is.null(spidr)) {
  spidr_ddr <- spidr %>%
    dplyr::filter(GeneA %in% ddr_genes, GeneB %in% ddr_genes)
  
  if (nrow(spidr_ddr) > 0) {
    sp_gene <- dplyr::bind_rows(
      spidr_ddr %>% dplyr::select(Gene = GeneA, SLStrength, SL_Label),
      spidr_ddr %>% dplyr::select(Gene = GeneB, SLStrength, SL_Label)
    ) %>%
      dplyr::group_by(Gene) %>%
      dplyr::summarise(
        SPIDR_PairCount = dplyr::n(),
        SPIDR_PositivePairCount = sum(SL_Label == 1, na.rm = TRUE),
        SPIDR_PositivePairFraction = mean(SL_Label == 1, na.rm = TRUE),
        SPIDR_MeanSLStrength = mean(SLStrength, na.rm = TRUE),
        SPIDR_Top10SLStrength = mean(
          head(
            sort(SLStrength, decreasing = TRUE),
            max(1, ceiling(length(SLStrength) * 0.10))
          ),
          na.rm = TRUE
        ),
        .groups = "drop"
      )
    
    spv <- MO_DDRweight %>%
      dplyr::select(Gene, MO_DDRweight) %>%
      dplyr::left_join(sp_gene, by = "Gene")
    
    for (v in setdiff(colnames(spv), c("Gene", "MO_DDRweight"))) {
      spv[[v]][is.na(spv[[v]])] <- 0
    }
    
    sp_cor <- dplyr::bind_rows(lapply(c(
      "SPIDR_PositivePairCount",
      "SPIDR_PositivePairFraction",
      "SPIDR_MeanSLStrength",
      "SPIDR_Top10SLStrength"
    ), function(v) {
      cr <- safe_cor(spv$MO_DDRweight, spv[[v]])
      data.frame(
        Feature = v,
        Cor = cr["cor"],
        P = cr["p"]
      )
    }))
    
    sp_cor$FDR <- p.adjust(sp_cor$P, "BH")
    
    p_sp <- ggplot(
      sp_cor,
      aes(x = reorder(Feature, Cor), y = Cor)
    ) +
      geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
      geom_col(fill = "#805AD5", width = 0.72) +
      coord_flip() +
      theme_bw(base_size = 12) +
      labs(
        x = NULL,
        y = "Spearman correlation with MO-DDRweight",
        title = "SPIDR external functional validation"
      )
    
    safe_ggsave(
      file.path(RES_DIR, "Fig2G_SPIDR_external_functional_validation.pdf"),
      p_sp,
      6,
      4.5
    )
  }
} else {
  message("SPIDR file not found or unreadable; Fig2G skipped.")
}

############################################################
# 9. Fig3A: DEG analysis
############################################################

cat("\nRunning Fig3 multi-omics analyses...\n")

de_samples <- intersect(
  mo_score$Sample[mo_score$SampleClass == "Tumor" & !is.na(mo_score$MO_DDRscore_group)],
  colnames(tcga_expr)
)

group <- factor(
  mo_score$MO_DDRscore_group[match(de_samples, mo_score$Sample)],
  levels = c("Low", "High")
)

expr_de <- log2(tcga_expr[, de_samples, drop = FALSE] + 1)

if (has_pkg("limma")) {
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  
  fit <- limma::lmFit(expr_de, design)
  cont <- limma::makeContrasts(High - Low, levels = design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont))
  deg <- limma::topTable(fit2, number = Inf, adjust.method = "BH")
  deg$Gene <- rownames(deg)
  deg <- deg[, c("Gene", setdiff(colnames(deg), "Gene"))]
  deg$FDR <- deg$adj.P.Val
} else {
  deg <- dplyr::bind_rows(lapply(rownames(expr_de), function(g) {
    x <- expr_de[g, group == "High"]
    y <- expr_de[g, group == "Low"]
    
    data.frame(
      Gene = g,
      logFC = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
      P.Value = tryCatch(t.test(x, y)$p.value, error = function(e) NA_real_)
    )
  }))
  
  deg$FDR <- p.adjust(deg$P.Value, "BH")
  deg$adj.P.Val <- deg$FDR
}

deg$Significance <- "NS"
deg$Significance[deg$FDR < DEG_FDR & deg$logFC > DEG_LOGFC] <- "Up"
deg$Significance[deg$FDR < DEG_FDR & deg$logFC < -DEG_LOGFC] <- "Down"
deg <- deg[order(deg$FDR, -abs(deg$logFC)), ]

p_vol <- ggplot(
  deg,
  aes(logFC, -log10(FDR + 1e-300), color = Significance)
) +
  geom_point(size = 0.8, alpha = 0.8) +
  scale_color_manual(
    values = c(Down = LOW_COL, NS = "grey70", Up = HIGH_COL),
    breaks = c("Down", "NS", "Up")
  ) +
  geom_vline(xintercept = c(-DEG_LOGFC, DEG_LOGFC), linetype = 2, color = "grey50") +
  geom_hline(yintercept = -log10(DEG_FDR), linetype = 2, color = "grey50") +
  theme_bw(base_size = 12) +
  labs(
    x = "log2 fold change (High vs Low)",
    y = "-log10(FDR)",
    color = NULL,
    title = "DEGs by MO-DDRscore group"
  )

safe_ggsave(
  file.path(RES_DIR, "Fig3A_DEG_volcano.pdf"),
  p_vol,
  6,
  5
)

topg <- unique(c(
  head(deg$Gene[deg$Significance == "Up"], MAX_HEATMAP_GENES / 2),
  head(deg$Gene[deg$Significance == "Down"], MAX_HEATMAP_GENES / 2)
))

topg <- intersect(topg, rownames(expr_de))

if (length(topg) >= 5 && has_pkg("pheatmap")) {
  sample_order_df <- data.frame(
    Sample = de_samples,
    Group = group,
    MO_DDRscore_raw = mo_score$MO_DDRscore_raw[match(de_samples, mo_score$Sample)]
  ) %>%
    dplyr::arrange(Group, MO_DDRscore_raw)
  
  hm <- z_rows(expr_de[topg, sample_order_df$Sample, drop = FALSE])
  hm[hm > 2.5] <- 2.5
  hm[hm < -2.5] <- -2.5
  
  anno <- data.frame(
    MO_DDRscore_group = sample_order_df$Group,
    MO_DDRscore = sample_order_df$MO_DDRscore_raw
  )
  rownames(anno) <- sample_order_df$Sample
  
  ann_colors <- list(
    MO_DDRscore_group = c(Low = LOW_COL, High = HIGH_COL)
  )
  
  gap_pos <- sum(sample_order_df$Group == "Low", na.rm = TRUE)
  if (gap_pos <= 0 || gap_pos >= nrow(sample_order_df)) gap_pos <- NULL
  
  safe_pdf(file.path(RES_DIR, "Fig3A_top_DEG_heatmap_ordered_by_group.pdf"), 8, 7)
  pheatmap::pheatmap(
    hm,
    annotation_col = anno,
    annotation_colors = ann_colors,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    gaps_col = gap_pos,
    show_colnames = FALSE,
    fontsize_row = 7,
    color = colorRampPalette(c(BLUE_COL, "white", RED_COL))(100),
    main = "Top DEGs ordered by MO-DDRscore group"
  )
  dev.off()
}

############################################################
# 10. Fig3B: Hallmark pathway analysis
############################################################

get_sets <- function() {
  if (has_pkg("msigdbr")) {
    ms <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", category = "H"), error = function(e) NULL)
    if (!is.null(ms) && nrow(ms) > 0) {
      return(lapply(split(ms$gene_symbol, ms$gs_name), clean_gene))
    }
  }
  
  list(
    HALLMARK_DNA_REPAIR = score_sets$DNA_Repair,
    HALLMARK_G2M_CHECKPOINT = score_sets$Cell_Cycle_G2M,
    HALLMARK_E2F_TARGETS = c(
      "E2F1","E2F2","E2F3","MCM2","MCM3","MCM4","MCM5",
      "MCM6","MCM7","PCNA","TYMS","TK1","CDC6","CDK2","CCNE1"
    ),
    HALLMARK_MYC_TARGETS_V1 = c(
      "MYC","NPM1","NCL","RPLP0","RPS3","LDHA","ODC1",
      "CAD","HSPD1","HSPE1","MCM4","MCM5","MCM6"
    ),
    HALLMARK_P53_PATHWAY = c(
      "TP53","CDKN1A","MDM2","GADD45A","BAX","BBC3",
      "PMAIP1","DDB2","RRM2B","SESN1","SESN2"
    ),
    HALLMARK_HYPOXIA = c(
      "HIF1A","VEGFA","CA9","LDHA","SLC2A1","ENO1",
      "PGK1","BNIP3","NDRG1","EGLN3"
    ),
    HALLMARK_INTERFERON_GAMMA_RESPONSE = c(
      "STAT1","IRF1","CXCL9","CXCL10","CXCL11","GBP1",
      "GBP5","IDO1","HLA-DRA","HLA-DRB1"
    ),
    HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = c(
      "VIM","CDH2","SNAI1","SNAI2","TWIST1","ZEB1","ZEB2",
      "COL1A1","COL1A2","FN1","MMP2","MMP9"
    )
  )
}

sets <- get_sets()
expr_z <- z_rows(expr_tumor_log)

hs <- list()

for (nm in names(sets)) {
  gs <- intersect(clean_gene(sets[[nm]]), rownames(expr_z))
  if (length(gs) >= 3) {
    hs[[nm]] <- colMeans(expr_z[gs, , drop = FALSE], na.rm = TRUE)
  }
}

hall_scores <- as.data.frame(hs, check.names = FALSE)

if (ncol(hall_scores) > 0) {
  hall_scores$Sample <- rownames(hall_scores)
  
  hall_scores <- hall_scores %>%
    dplyr::left_join(
      mo_score[, c("Sample", "MO_DDRscore_raw", "MO_DDRscore_group")],
      by = "Sample"
    )
  
  hcols <- setdiff(colnames(hall_scores), c("Sample", "MO_DDRscore_raw", "MO_DDRscore_group"))
  
  hstat <- dplyr::bind_rows(lapply(hcols, function(pw) {
    dd <- hall_scores[
      is.finite(hall_scores[[pw]]) &
        !is.na(hall_scores$MO_DDRscore_group),
    ]
    
    cr <- safe_cor(dd$MO_DDRscore_raw, dd[[pw]])
    
    data.frame(
      Pathway = pw,
      Cor = cr["cor"],
      CorP = cr["p"],
      Diff = median(dd[[pw]][dd$MO_DDRscore_group == "High"], na.rm = TRUE) -
        median(dd[[pw]][dd$MO_DDRscore_group == "Low"], na.rm = TRUE),
      WilcoxP = tryCatch(
        wilcox.test(dd[[pw]] ~ dd$MO_DDRscore_group)$p.value,
        error = function(e) NA_real_
      )
    )
  }))
  
  hstat$FDR <- p.adjust(hstat$WilcoxP, "BH")
  hstat <- hstat[order(hstat$FDR, -abs(hstat$Diff)), ]
  
  hstat$Pathway_label <- gsub("^HALLMARK_", "", hstat$Pathway)
  hstat$Pathway_label <- gsub("_", " ", hstat$Pathway_label)
  
  p_h <- ggplot(
    head(hstat, 25),
    aes(x = reorder(Pathway_label, Diff), y = Diff)
  ) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
    geom_point(
      aes(size = -log10(WilcoxP + 1e-300), color = Cor),
      alpha = 0.9
    ) +
    scale_color_gradient2(low = BLUE_COL, mid = "white", high = RED_COL, midpoint = 0) +
    coord_flip() +
    theme_bw(base_size = 11) +
    labs(
      x = NULL,
      y = "Median difference (High - Low)",
      size = "-log10(P)",
      color = "Spearman r",
      title = "Hallmark pathway alteration"
    )
  
  safe_ggsave(
    file.path(RES_DIR, "Fig3B_Hallmark_pathway_bubble.pdf"),
    p_h,
    7.2,
    6
  )
  
  if (has_pkg("pheatmap")) {
    top_pw <- head(hstat$Pathway, min(20, nrow(hstat)))
    
    hall_order_df <- hall_scores %>%
      dplyr::filter(!is.na(MO_DDRscore_group), is.finite(MO_DDRscore_raw)) %>%
      dplyr::arrange(MO_DDRscore_group, MO_DDRscore_raw)
    
    hm <- as.matrix(t(hall_order_df[, top_pw, drop = FALSE]))
    colnames(hm) <- hall_order_df$Sample
    hm <- z_rows(hm)
    hm[hm > 2.5] <- 2.5
    hm[hm < -2.5] <- -2.5
    
    rownames(hm) <- gsub("^HALLMARK_", "", rownames(hm))
    rownames(hm) <- gsub("_", " ", rownames(hm))
    
    anno <- data.frame(
      MO_DDRscore_group = hall_order_df$MO_DDRscore_group,
      MO_DDRscore = hall_order_df$MO_DDRscore_raw
    )
    rownames(anno) <- hall_order_df$Sample
    
    ann_colors <- list(
      MO_DDRscore_group = c(Low = LOW_COL, High = HIGH_COL)
    )
    
    gap_pos <- sum(hall_order_df$MO_DDRscore_group == "Low", na.rm = TRUE)
    if (gap_pos <= 0 || gap_pos >= nrow(hall_order_df)) gap_pos <- NULL
    
    safe_pdf(file.path(RES_DIR, "Fig3B_Hallmark_pathway_heatmap.pdf"), 9, 6.5)
    pheatmap::pheatmap(
      hm,
      annotation_col = anno,
      annotation_colors = ann_colors,
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      gaps_col = gap_pos,
      show_colnames = FALSE,
      color = colorRampPalette(c(BLUE_COL, "white", RED_COL))(100),
      fontsize_row = 7,
      main = "Hallmark pathway scores ordered by MO-DDRscore group"
    )
    dev.off()
  }
}

############################################################
# 11. Fig3C: Driver mutation frequency
############################################################

driver_genes <- c(
  "TP53","KRAS","EGFR","STK11","KEAP1","BRAF","MET","ERBB2",
  "ALK","ROS1","RET","NF1","PIK3CA","RBM10","SMARCA4"
)

if (!is.null(mut)) {
  pg <- tumor_mo %>%
    dplyr::select(Patient, MO_DDRscore_group) %>%
    dplyr::distinct()
  
  mut2 <- mut %>%
    dplyr::left_join(pg, by = "Patient") %>%
    dplyr::filter(!is.na(MO_DDRscore_group))
  
  ngrp <- table(pg$MO_DDRscore_group)
  
  dstat <- dplyr::bind_rows(lapply(driver_genes, function(g) {
    pats <- unique(mut2$Patient[mut2$Hugo_Symbol == g])
    
    ht <- as.numeric(ngrp["High"])
    lt <- as.numeric(ngrp["Low"])
    
    hm <- sum(pg$Patient[pg$MO_DDRscore_group == "High"] %in% pats)
    lm <- sum(pg$Patient[pg$MO_DDRscore_group == "Low"] %in% pats)
    
    p <- tryCatch(
      fisher.test(matrix(c(hm, ht - hm, lm, lt - lm), nrow = 2))$p.value,
      error = function(e) NA_real_
    )
    
    data.frame(
      Gene = g,
      High_Freq = hm / ht,
      Low_Freq = lm / lt,
      Diff = hm / ht - lm / lt,
      P = p
    )
  }))
  
  dstat$FDR <- p.adjust(dstat$P, "BH")
  
  dlong <- dstat %>%
    dplyr::select(Gene, High_Freq, Low_Freq) %>%
    tidyr::pivot_longer(
      cols = c(High_Freq, Low_Freq),
      names_to = "Group",
      values_to = "Frequency"
    )
  
  dlong$Group <- ifelse(dlong$Group == "High_Freq", "High", "Low")
  dlong$Group <- factor(dlong$Group, levels = c("Low", "High"))
  
  gene_order <- dstat %>%
    dplyr::mutate(TotalFreq = High_Freq + Low_Freq) %>%
    dplyr::arrange(TotalFreq) %>%
    dplyr::pull(Gene)
  
  dlong$Gene <- factor(dlong$Gene, levels = gene_order)
  
  p_drv <- ggplot(
    dlong,
    aes(x = Gene, y = Frequency, fill = Group)
  ) +
    geom_col(position = position_dodge(width = 0.75), width = 0.68) +
    coord_flip() +
    scale_fill_manual(values = c(Low = LOW_COL, High = HIGH_COL)) +
    theme_bw(base_size = 12) +
    labs(
      x = NULL,
      y = "Mutation frequency",
      fill = "MO-DDRscore group",
      title = "LUAD driver mutation frequency"
    )
  
  safe_ggsave(
    file.path(RES_DIR, "Fig3C_driver_mutation_frequency_barplot.pdf"),
    p_drv,
    6,
    5
  )
}

############################################################
# 12. Fig3D: Genome-wide CNV burden
############################################################

if (!is.null(cnv_df)) {
  cnv_plot <- cnv_df %>%
    dplyr::left_join(
      tumor_mo[, c("Patient", "MO_DDRscore_group")],
      by = "Patient"
    ) %>%
    dplyr::filter(!is.na(MO_DDRscore_group))
  
  cnv_long <- cnv_plot %>%
    tidyr::pivot_longer(
      cols = dplyr::any_of(c("CNVBurden", "AmpBurden", "DelBurden")),
      names_to = "Feature",
      values_to = "Value"
    )
  
  cnv_long$Feature <- factor(
    cnv_long$Feature,
    levels = c("CNVBurden", "AmpBurden", "DelBurden"),
    labels = c("CNV burden", "Amplification burden", "Deletion burden")
  )
  
  cnv_long$MO_DDRscore_group <- factor(
    cnv_long$MO_DDRscore_group,
    levels = c("Low", "High")
  )
  
  p_cnv <- ggplot(
    cnv_long,
    aes(MO_DDRscore_group, Value, fill = MO_DDRscore_group)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 0.55, alpha = 0.35) +
    facet_wrap(~ Feature, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c(Low = LOW_COL, High = HIGH_COL)) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    labs(
      x = NULL,
      y = "Burden",
      title = "Genome-wide CNV burden by MO-DDRscore group"
    )
  
  safe_ggsave(
    file.path(RES_DIR, "Fig3D_genome_wide_CNV_burden_boxplot.pdf"),
    p_cnv,
    7,
    4
  )
}

############################################################
# 13. Fig3E: Immune / HLA / CYT feature scores
############################################################

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
  if (length(gs) >= 2) {
    im[[nm]] <- colMeans(expr_imm_z[gs, , drop = FALSE], na.rm = TRUE)
  }
}

immune_score_df <- as.data.frame(im, check.names = FALSE)

if (ncol(immune_score_df) > 0) {
  immune_score_df$Sample <- rownames(immune_score_df)
  
  immune_score_df <- immune_score_df %>%
    dplyr::left_join(
      mo_score[, c("Sample", "MO_DDRscore_raw", "MO_DDRscore_group")],
      by = "Sample"
    )
  
  imcols <- setdiff(
    colnames(immune_score_df),
    c("Sample", "MO_DDRscore_raw", "MO_DDRscore_group")
  )
  
  istat <- dplyr::bind_rows(lapply(imcols, function(v) {
    dd <- immune_score_df[
      is.finite(immune_score_df[[v]]) &
        !is.na(immune_score_df$MO_DDRscore_group),
    ]
    
    cr <- safe_cor(dd$MO_DDRscore_raw, dd[[v]])
    
    data.frame(
      Feature = v,
      Cor = cr["cor"],
      CorP = cr["p"],
      Diff = median(dd[[v]][dd$MO_DDRscore_group == "High"], na.rm = TRUE) -
        median(dd[[v]][dd$MO_DDRscore_group == "Low"], na.rm = TRUE),
      WilcoxP = tryCatch(
        wilcox.test(dd[[v]] ~ dd$MO_DDRscore_group)$p.value,
        error = function(e) NA_real_
      )
    )
  }))
  
  istat$FDR <- p.adjust(istat$WilcoxP, "BH")
  
  istat$Feature_label <- gsub("_", " ", istat$Feature)
  istat$Feature_label <- gsub("T cell CD8", "CD8 T cell", istat$Feature_label)
  
  p_im <- ggplot(
    istat,
    aes(x = reorder(Feature_label, Diff), y = Diff)
  ) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
    geom_point(
      aes(size = -log10(WilcoxP + 1e-300), color = Cor),
      alpha = 0.9
    ) +
    scale_color_gradient2(low = BLUE_COL, mid = "white", high = RED_COL, midpoint = 0) +
    coord_flip() +
    theme_bw(base_size = 12) +
    labs(
      x = NULL,
      y = "Median difference (High - Low)",
      size = "-log10(P)",
      color = "Spearman r",
      title = "Immune, HLA and cytolytic feature differences"
    )
  
  safe_ggsave(
    file.path(RES_DIR, "Fig3E_immune_HLA_CYT_feature_differences.pdf"),
    p_im,
    6.5,
    5
  )
}

############################################################
# 14. Fig3F: Immune checkpoint expression
############################################################

checkpoint_genes <- intersect(
  c(
    "CD274","PDCD1","CTLA4","LAG3","TIGIT","HAVCR2",
    "PDCD1LG2","ICOS","IDO1","CD80","CD86","TNFRSF9",
    "TNFRSF4","CD40","CD40LG","VSIR","SIGLEC15"
  ),
  rownames(tcga_expr)
)

if (length(checkpoint_genes) >= 3) {
  chk_df <- as.data.frame(
    t(log2(tcga_expr[checkpoint_genes, tumor_samples, drop = FALSE] + 1)),
    check.names = FALSE
  )
  
  chk_df$Sample <- rownames(chk_df)
  
  chk_df <- chk_df %>%
    dplyr::left_join(
      mo_score[, c("Sample", "MO_DDRscore_group")],
      by = "Sample"
    ) %>%
    dplyr::filter(!is.na(MO_DDRscore_group))
  
  chk_long <- chk_df %>%
    tidyr::pivot_longer(
      cols = all_of(checkpoint_genes),
      names_to = "Gene",
      values_to = "Expression"
    )
  
  chk_long$MO_DDRscore_group <- factor(
    chk_long$MO_DDRscore_group,
    levels = c("Low", "High")
  )
  
  chk_stat <- chk_long %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      Diff = median(Expression[MO_DDRscore_group == "High"], na.rm = TRUE) -
        median(Expression[MO_DDRscore_group == "Low"], na.rm = TRUE),
      P = tryCatch(
        wilcox.test(Expression ~ MO_DDRscore_group)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    )
  
  chk_stat$FDR <- p.adjust(chk_stat$P, "BH")
  
  showg <- head(chk_stat$Gene[order(chk_stat$FDR)], min(12, nrow(chk_stat)))
  
  p_chk <- ggplot(
    chk_long[chk_long$Gene %in% showg, ],
    aes(MO_DDRscore_group, Expression, fill = MO_DDRscore_group)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.12, size = 0.45, alpha = 0.25) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = c(Low = LOW_COL, High = HIGH_COL)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none") +
    labs(
      x = NULL,
      y = "log2(TPM+1)",
      title = "Immune checkpoint expression"
    )
  
  safe_ggsave(
    file.path(RES_DIR, "Fig3F_checkpoint_expression_boxplot.pdf"),
    p_chk,
    8,
    6
  )
}

############################################################
# 15. Fig3G: Survival analysis
############################################################

surv_df <- tumor_mo %>%
  dplyr::inner_join(tcga_clin, by = "Patient") %>%
  dplyr::filter(
    is.finite(time),
    time > 0,
    status %in% c(0, 1),
    !is.na(MO_DDRscore_group)
  )

surv_df$MO_DDRscore_group <- factor(
  surv_df$MO_DDRscore_group,
  levels = c("Low", "High")
)

if (nrow(surv_df) >= 30 && sum(surv_df$status == 1, na.rm = TRUE) >= 5 && has_pkg("survminer")) {
  fit <- survival::survfit(
    Surv(time, status) ~ MO_DDRscore_group,
    data = surv_df
  )
  
  p_km <- survminer::ggsurvplot(
    fit,
    data = surv_df,
    pval = TRUE,
    conf.int = FALSE,
    risk.table = FALSE,
    palette = c(LOW_COL, HIGH_COL),
    legend.title = "",
    legend.labs = c("Low", "High"),
    ggtheme = theme_bw(base_size = 12),
    title = paste0(FOCUS_CANCER, " MO-DDRscore survival")
  )$plot +
    labs(
      x = "Time (days)",
      y = "Overall survival probability"
    )
  
  safe_ggsave(
    file.path(RES_DIR, "Fig3G_MO_DDRscore_KM_survival.pdf"),
    p_km,
    6,
    5
  )
} else {
  message("survminer not installed or survival samples/events insufficient; KM skipped.")
}

if (has_pkg("timeROC")) {
  roc_obj <- tryCatch(
    timeROC::timeROC(
      T = surv_df$time,
      delta = surv_df$status,
      marker = surv_df$MO_DDRscore_raw,
      cause = 1,
      times = c(365, 1095, 1825),
      iid = TRUE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(roc_obj)) {
    safe_pdf(file.path(RES_DIR, "Fig3G_MO_DDRscore_timeROC.pdf"), 5.5, 5)
    plot(roc_obj, time = 365, col = "#E64B35", title = FALSE)
    plot(roc_obj, time = 1095, add = TRUE, col = "#4DBBD5")
    plot(roc_obj, time = 1825, add = TRUE, col = "#2F855A")
    legend(
      "bottomright",
      legend = paste0(
        c("1-year", "3-year", "5-year"),
        " AUC=",
        sprintf("%.3f", as.numeric(roc_obj$AUC))
      ),
      col = c("#E64B35", "#4DBBD5", "#2F855A"),
      lwd = 2,
      bty = "n"
    )
    title(paste0(FOCUS_CANCER, " MO-DDRscore time-dependent ROC"))
    dev.off()
  }
} else {
  message("timeROC not installed; timeROC skipped.")
}

############################################################
# 16. Final report
############################################################

cat("\nMO-DDRscore pipeline finished.\n")
cat("All PDF figures saved to:\n")
cat(RES_DIR, "\n\n")

cat("Key figure files:\n")
cat("Fig2A_gene_level_multiomics_evidence_heatmap.pdf\n")
cat("Fig2B_entropy_adaptive_evidence_weights.pdf\n")
cat("Fig2C_MO_DDRweight_distribution.pdf\n")
cat("Fig2D_top_MO_DDRweight_genes.pdf\n")
cat("Fig2E_MO_DDRscore_ranked_stratification.pdf\n")
cat("Fig2F_MO_DDRscore_ablation_comparison.pdf\n")
cat("Fig2G_SPIDR_external_functional_validation.pdf\n")
cat("Fig3A_DEG_volcano.pdf\n")
cat("Fig3A_top_DEG_heatmap_ordered_by_group.pdf\n")
cat("Fig3B_Hallmark_pathway_bubble.pdf\n")
cat("Fig3B_Hallmark_pathway_heatmap.pdf\n")
cat("Fig3C_driver_mutation_frequency_barplot.pdf\n")
cat("Fig3D_genome_wide_CNV_burden_boxplot.pdf\n")
cat("Fig3E_immune_HLA_CYT_feature_differences.pdf\n")
cat("Fig3F_checkpoint_expression_boxplot.pdf\n")
cat("Fig3G_MO_DDRscore_KM_survival.pdf\n")
cat("Fig3G_MO_DDRscore_timeROC.pdf\n")