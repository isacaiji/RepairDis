############################################################
# Figure 2: official-style multi-omics analysis for MO-DDRscore
#
# Reference-style modules:
#   A DDR pathway activity: GSVA/ssGSEA with curated DDR pathway file
#   B DEG heatmap: limma + ComplexHeatmap with clinical annotations
#   C Hallmark pathways: MSigDB Hallmark GSVA and GSEA
#   D Somatic mutation: maftools oncoplot for DDR genes
#   E TMB: maftools::tmb
#   F CNV: GISTIC2 / FGA-FGG-FGL official outputs if available
#   G Immune microenvironment: ESTIMATE official package if installed;
#      optional ssGSEA using external standard immune GMT
#   H Immune checkpoints: expression comparison
#
# No ranking landscape, no custom CNV burden as main result.
############################################################

options(stringsAsFactors = FALSE)
set.seed(20260513)

############################
# 0. Parameters
############################

PROJECT_DIR <- "D:/R_workspace/评分/AD_DDR_project"
PROC_DIR    <- file.path(PROJECT_DIR, "01_processed")
DATA_DIR    <- file.path(PROJECT_DIR, "00_data")
OUT_DIR     <- file.path(PROJECT_DIR, "02_Figure2_official_multiomics")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FOCUS_CANCER <- "LUAD"
TCGA_DIR <- "D:/R/R_workspace/梁老师文件/TCGA"

TCGA_EXPR_FILE <- file.path(TCGA_DIR, "mRNA_exp_TPM_only_TCGA/mRNA_exp_TPM_only_TCGA",
                            paste0("TCGA-", FOCUS_CANCER, ".gene_expression_TPM.tsv"))
TCGA_CLIN_PROCESSED_FILE <- file.path(PROC_DIR, paste0(FOCUS_CANCER, "_clinical_processed.csv"))
MO_SCORE_FILE <- file.path(PROC_DIR, paste0(FOCUS_CANCER, "_MO_DDRscore.csv"))

# Required for DDR pathway analysis
DDR_GENE_FILE <- file.path(DATA_DIR, "DDR_236_genes.csv")   # Gene, Primary_Pathway

# Required for maftools mutation/TMB
MAF_FILE <- file.path(TCGA_DIR, "mutation", paste0(FOCUS_CANCER, ".txt"))

# Optional official CNV outputs
GISTIC_DIR <- file.path(DATA_DIR, "GISTIC2_LUAD")
GISTIC_ALL_LESIONS <- file.path(GISTIC_DIR, "all_lesions.conf_99.txt")
GISTIC_AMP_GENES   <- file.path(GISTIC_DIR, "amp_genes.conf_99.txt")
GISTIC_DEL_GENES   <- file.path(GISTIC_DIR, "del_genes.conf_99.txt")
GISTIC_SCORES      <- file.path(GISTIC_DIR, "scores.gistic")
FGA_FILE <- file.path(DATA_DIR, "LUAD_FGA_FGG_FGL.csv")  # Patient/Sample + FGA/FGG/FGL

# Optional standard immune signature GMT, e.g. Charoentong 28 immune signatures
IMMUNE_GMT_FILE <- file.path(DATA_DIR, "Charoentong_28immune.gmt")

# Optional checkpoint gene file, one gene per row
ICI_GENE_FILE <- file.path(DATA_DIR, "ICI_genes_35.txt")

DEG_ADJ_P <- 0.05
DEG_LOGFC <- 1.00
TOP_HEAT_GENES <- 60
TOP_SHOW_N <- 20
TOP_MUT_GENES <- 25
TOP_CHECKPOINTS <- 20

############################
# 1. Packages
############################

cran_pkgs <- c("data.table","dplyr","tidyr","ggplot2","ggrepel","ggpubr",
               "pheatmap","ComplexHeatmap","circlize","grid","patchwork","survival")
bioc_pkgs <- c("limma","GSVA","GSEABase","clusterProfiler","msigdbr","maftools")

install_if_missing <- function(pkgs, bioc = FALSE) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
        BiocManager::install(p, ask = FALSE, update = FALSE)
      } else {
        install.packages(p, repos = "https://cloud.r-project.org")
      }
    }
  }
}

install_if_missing(cran_pkgs, FALSE)
install_if_missing(bioc_pkgs, TRUE)

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr); library(ggplot2)
  library(ggrepel); library(ggpubr); library(pheatmap); library(ComplexHeatmap)
  library(circlize); library(grid); library(patchwork); library(survival)
  library(limma); library(GSVA); library(GSEABase); library(clusterProfiler)
  library(msigdbr); library(maftools)
})

############################
# 2. Helper functions
############################

save_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file)
}

safe_ggsave <- function(file, p, w = 7, h = 5) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file, p, width = w, height = h, device = "pdf", useDingbats = FALSE)
}

clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub(" ///.*$| //.*$|;.*$|,.*$|\\s*\\([^\\)]*\\)$", "", x)
  x[x %in% c("", "NA", "---", "NULL", "N/A", "NAN")] <- NA
  x
}

patient_id <- function(x) substr(gsub("\\.", "-", as.character(x)), 1, 12)
sample_type <- function(x) substr(gsub("\\.", "-", as.character(x)), 14, 15)

sample_class <- function(x) {
  code <- sample_type(x)
  ifelse(code %in% c("01","02","03","05","06","07"), "Tumor",
         ifelse(code %in% c("10","11","12","13","14"), "Normal", "Other"))
}

pick_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) NA_character_ else hit[1]
}

simple_gender <- function(x) {
  x <- toupper(as.character(x))
  dplyr::case_when(
    x %in% c("MALE","M","1") ~ "Male",
    x %in% c("FEMALE","F","0") ~ "Female",
    TRUE ~ NA_character_
  )
}

simple_stage <- function(x) {
  x <- toupper(as.character(x))
  dplyr::case_when(
    grepl("STAGE I[^I]|^I$|^IA|^IB", x) ~ "I",
    grepl("STAGE II[^I]|^II$|^IIA|^IIB", x) ~ "II",
    grepl("STAGE III|^III|^IIIA|^IIIB|^IIIC", x) ~ "III",
    grepl("STAGE IV|^IV|^IVA|^IVB", x) ~ "IV",
    TRUE ~ NA_character_
  )
}

simple_status <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  dplyr::case_when(x == 1 ~ "Dead", x == 0 ~ "Alive", TRUE ~ NA_character_)
}

safe_ssgsea <- function(expr_mat, gene_sets) {
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "numeric"
  gene_sets <- lapply(gene_sets, function(x) intersect(unique(clean_gene(x)), rownames(expr_mat)))
  gene_sets <- gene_sets[sapply(gene_sets, length) >= 2]
  if (length(gene_sets) == 0) stop("No valid gene sets after matching expression matrix.")
  
  res <- tryCatch({
    GSVA::gsva(expr_mat, gene_sets, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE, verbose = FALSE)
  }, error = function(e) {
    if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
      param <- GSVA::ssgseaParam(expr_mat, gene_sets, normalize = TRUE)
      GSVA::gsva(param, verbose = FALSE)
    } else stop(e)
  })
  as.matrix(res)
}

box_compare <- function(df, group_col, value_col, title = NULL, ylab = NULL) {
  df <- df %>% filter(!is.na(.data[[group_col]]), is.finite(.data[[value_col]]))
  ggplot(df, aes(x = .data[[group_col]], y = .data[[value_col]], fill = .data[[group_col]])) +
    geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.16, size = 0.65, alpha = 0.45) +
    ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", size = 3.2) +
    scale_fill_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(), axis.text = element_text(color = "black"),
          legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(x = NULL, y = ylab, title = title)
}

############################
# 3. Load data
############################

TCGA_EXPR_RDS <- file.path(PROC_DIR, paste0(FOCUS_CANCER, "_tcga_expr_tpm_matrix.rds"))
if (exists("tcga_expr")) {
  cat("tcga_expr exists in environment.\n")
} else if (file.exists(TCGA_EXPR_RDS)) {
  tcga_expr <- readRDS(TCGA_EXPR_RDS)
} else {
  expr_raw <- fread(TCGA_EXPR_FILE, data.table = FALSE, check.names = FALSE)
  gene_col <- colnames(expr_raw)[1]
  expr_raw[[gene_col]] <- clean_gene(expr_raw[[gene_col]])
  expr_raw <- expr_raw %>%
    filter(!is.na(.data[[gene_col]]), .data[[gene_col]] != "") %>%
    group_by(.data[[gene_col]]) %>%
    summarise(across(everything(), ~ suppressWarnings(mean(as.numeric(.x), na.rm = TRUE))), .groups = "drop")
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

tcga_clin <- fread(TCGA_CLIN_PROCESSED_FILE, data.table = FALSE, check.names = FALSE) %>%
  mutate(Patient = as.character(Patient),
         time = suppressWarnings(as.numeric(time)),
         status = suppressWarnings(as.numeric(status)))

mo_score <- fread(MO_SCORE_FILE, data.table = FALSE, check.names = FALSE) %>%
  mutate(Sample = gsub("\\.", "-", as.character(Sample)),
         Patient = as.character(Patient),
         SampleType = sample_type(Sample),
         MO_DDRscore_raw = suppressWarnings(as.numeric(MO_DDRscore_raw)),
         MO_DDRscore_group = as.character(MO_DDRscore_group))

tumor_samples <- colnames(tcga_expr)[sample_class(colnames(tcga_expr)) == "Tumor"]

score_df <- mo_score %>%
  filter(Sample %in% tumor_samples, Sample %in% colnames(tcga_expr),
         MO_DDRscore_group %in% c("Low","High"), is.finite(MO_DDRscore_raw)) %>%
  arrange(Patient, SampleType) %>%
  distinct(Patient, .keep_all = TRUE)

score_df$MO_DDRscore_group <- factor(score_df$MO_DDRscore_group, levels = c("Low","High"))
expr_log <- log2(tcga_expr[, score_df$Sample, drop = FALSE] + 1)

cat("Figure 2 samples:", nrow(score_df), "\n")
print(table(score_df$MO_DDRscore_group))
save_csv(score_df, file.path(OUT_DIR, "Fig2_sample_group_table.csv"))

############################
# 4. Fig2A DDR pathway activity: GSVA/ssGSEA
############################

if (!file.exists(DDR_GENE_FILE)) stop("DDR_GENE_FILE not found: ", DDR_GENE_FILE)

ddr_anno <- data.table::fread(
  DDR_GENE_FILE,
  data.table = FALSE,
  check.names = FALSE
)

# 兼容 Gene / gene / Symbol 等列名
if (!"Gene" %in% colnames(ddr_anno)) {
  gene_col <- intersect(
    c("gene", "GeneSymbol", "Gene_Symbol", "Symbol", "symbol", "Hugo_Symbol"),
    colnames(ddr_anno)
  )[1]
  
  if (is.na(gene_col)) {
    colnames(ddr_anno)[1] <- "Gene"
  } else {
    colnames(ddr_anno)[colnames(ddr_anno) == gene_col] <- "Gene"
  }
}

# 兼容 Pathway / pathway / Primary_Pathway
if (!"Primary_Pathway" %in% colnames(ddr_anno)) {
  pathway_col <- intersect(
    c("Pathway", "pathway", "Primary_pathway", "primary_pathway",
      "PrimaryPathway", "Category", "category", "Type", "type"),
    colnames(ddr_anno)
  )[1]
  
  if (is.na(pathway_col)) {
    stop("DDR_GENE_FILE must contain a pathway column, such as Pathway or Primary_Pathway.")
  } else {
    colnames(ddr_anno)[colnames(ddr_anno) == pathway_col] <- "Primary_Pathway"
  }
}

ddr_anno <- ddr_anno %>%
  dplyr::mutate(
    Gene = clean_gene(Gene),
    Primary_Pathway = trimws(as.character(Primary_Pathway)),
    Primary_Pathway = dplyr::case_when(
      Primary_Pathway %in% c("DNA damage checkpoint", "DNA damage checkpoint ") ~ "Checkpoint",
      Primary_Pathway %in% c("Replication stress", "Replication Stress") ~ "Replication_stress",
      TRUE ~ Primary_Pathway
    )
  ) %>%
  dplyr::filter(
    !is.na(Gene),
    Gene != "",
    !is.na(Primary_Pathway),
    Primary_Pathway != ""
  ) %>%
  dplyr::distinct(Gene, Primary_Pathway)

ddr_sets <- split(ddr_anno$Gene, ddr_anno$Primary_Pathway)

ddr_sets <- lapply(ddr_sets, function(x) {
  intersect(unique(clean_gene(x)), rownames(expr_log))
})

ddr_sets <- ddr_sets[sapply(ddr_sets, length) >= 2]

if (length(ddr_sets) < 2) {
  stop("Too few DDR pathway gene sets matched to expression matrix.")
}

cat("\nDDR pathway sets used:\n")
print(sapply(ddr_sets, length))
ddr_sets <- lapply(ddr_sets, function(x) intersect(unique(x), rownames(expr_log)))
ddr_sets <- ddr_sets[sapply(ddr_sets, length) >= 2]
if (length(ddr_sets) < 2) stop("Too few DDR pathway gene sets matched.")

ddr_ssgsea <- safe_ssgsea(expr_log, ddr_sets)

ddr_long <- as.data.frame(t(ddr_ssgsea), check.names = FALSE) %>%
  mutate(Sample = rownames(.)) %>%
  pivot_longer(cols = -Sample, names_to = "Pathway", values_to = "Activity") %>%
  left_join(score_df %>% select(Sample, Patient, MO_DDRscore_group, MO_DDRscore_raw), by = "Sample")

ddr_stat <- ddr_long %>%
  group_by(Pathway) %>%
  summarise(P = wilcox.test(Activity ~ MO_DDRscore_group)$p.value,
            Median_High = median(Activity[MO_DDRscore_group == "High"], na.rm = TRUE),
            Median_Low = median(Activity[MO_DDRscore_group == "Low"], na.rm = TRUE),
            Diff = Median_High - Median_Low, .groups = "drop") %>%
  mutate(FDR = p.adjust(P, method = "BH")) %>% arrange(P)

save_csv(ddr_long, file.path(OUT_DIR, "Fig2A_DDR_pathway_ssGSEA_scores.csv"))
save_csv(ddr_stat, file.path(OUT_DIR, "Fig2A_DDR_pathway_statistics.csv"))

p_ddr <- ggplot(ddr_long, aes(MO_DDRscore_group, Activity, fill = MO_DDRscore_group)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.16, size = 0.5, alpha = 0.45) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 3) +
  facet_wrap(~Pathway, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black"),
        legend.position = "none", strip.text = element_text(face = "bold")) +
  labs(x = NULL, y = "ssGSEA score", title = "DDR pathway activity")
safe_ggsave(file.path(OUT_DIR, "Fig2A_DDR_pathway_activity_boxplot.pdf"), p_ddr, 9, 7)

ddr_mat <- ddr_ssgsea[, score_df$Sample, drop = FALSE]
ddr_mat <- t(scale(t(ddr_mat))); ddr_mat[ddr_mat > 2] <- 2; ddr_mat[ddr_mat < -2] <- -2
anno_col <- data.frame(Group = score_df$MO_DDRscore_group)
rownames(anno_col) <- score_df$Sample
pdf(file.path(OUT_DIR, "Fig2A_DDR_pathway_activity_heatmap.pdf"), width = 8, height = 5)
pheatmap(ddr_mat, annotation_col = anno_col, show_colnames = FALSE,
         color = colorRampPalette(c("#2B6CB0","white","#C53030"))(100),
         main = "DDR pathway activity")
dev.off()

############################
# 5. Fig2B DEG heatmap: limma + ComplexHeatmap
############################

group_factor <- factor(score_df$MO_DDRscore_group, levels = c("Low","High"))
design <- model.matrix(~0 + group_factor)
colnames(design) <- c("Low","High")
fit <- lmFit(expr_log, design)
contrast <- makeContrasts(High - Low, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contrast))
deg <- topTable(fit2, number = Inf, adjust.method = "BH")
deg$Gene <- rownames(deg)
deg <- deg %>%
  mutate(Regulation = case_when(
    adj.P.Val < DEG_ADJ_P & logFC >= DEG_LOGFC ~ "Up in High",
    adj.P.Val < DEG_ADJ_P & logFC <= -DEG_LOGFC ~ "Down in High",
    TRUE ~ "NS")) %>%
  arrange(adj.P.Val)
save_csv(deg, file.path(OUT_DIR, "Fig2B_DEG_High_vs_Low_limma.csv"))

top_up <- deg %>% filter(adj.P.Val < DEG_ADJ_P, logFC >= DEG_LOGFC) %>%
  arrange(adj.P.Val, desc(abs(logFC))) %>% slice_head(n = ceiling(TOP_HEAT_GENES/2)) %>% pull(Gene)
top_down <- deg %>% filter(adj.P.Val < DEG_ADJ_P, logFC <= -DEG_LOGFC) %>%
  arrange(adj.P.Val, desc(abs(logFC))) %>% slice_head(n = floor(TOP_HEAT_GENES/2)) %>% pull(Gene)
heat_genes <- unique(c(top_up, top_down))
if (length(heat_genes) < 10) {
  heat_genes <- deg %>% arrange(adj.P.Val) %>% slice_head(n = TOP_HEAT_GENES) %>% pull(Gene)
}
heat_genes <- intersect(heat_genes, rownames(expr_log))
save_csv(deg %>% filter(Gene %in% heat_genes), file.path(OUT_DIR, "Fig2B_DEG_heatmap_genes.csv"))

age_col <- pick_col(tcga_clin, c("age","age_at_initial_pathologic_diagnosis","age_at_diagnosis"))
gender_col <- pick_col(tcga_clin, c("gender","sex"))
stage_col <- pick_col(tcga_clin, c("stage","pathologic_stage","ajcc_pathologic_stage","clinical_stage"))

plot_meta <- score_df %>%
  left_join(tcga_clin, by = "Patient") %>%
  mutate(Age = if (!is.na(age_col)) suppressWarnings(as.numeric(.data[[age_col]])) else NA_real_,
         Gender = if (!is.na(gender_col)) simple_gender(.data[[gender_col]]) else NA_character_,
         Stage = if (!is.na(stage_col)) simple_stage(.data[[stage_col]]) else NA_character_,
         Status = simple_status(status)) %>%
  filter(is.finite(Age), Gender %in% c("Female","Male"),
         Stage %in% c("I","II","III","IV"), Status %in% c("Alive","Dead")) %>%
  arrange(MO_DDRscore_group, MO_DDRscore_raw)

hm_mat <- expr_log[heat_genes, plot_meta$Sample, drop = FALSE]
hm_mat <- t(scale(t(hm_mat))); hm_mat[hm_mat > 2] <- 2; hm_mat[hm_mat < -2] <- -2; hm_mat[is.na(hm_mat)] <- 0

anno_df <- data.frame(
  Group = factor(as.character(plot_meta$MO_DDRscore_group), levels = c("Low","High")),
  Score = plot_meta$MO_DDRscore_raw,
  Age = plot_meta$Age,
  Gender = factor(plot_meta$Gender, levels = c("Female","Male")),
  Stage = factor(plot_meta$Stage, levels = intersect(c("I","II","III","IV"), unique(plot_meta$Stage))),
  Status = factor(plot_meta$Status, levels = c("Alive","Dead")),
  row.names = plot_meta$Sample
)

top_anno <- HeatmapAnnotation(
  Group = anno_df$Group, Score = anno_df$Score, Age = anno_df$Age,
  Gender = anno_df$Gender, Stage = anno_df$Stage, Status = anno_df$Status,
  col = list(
    Group = c(Low = "#3B75AF", High = "#C84630"),
    Score = circlize::colorRamp2(quantile(anno_df$Score, c(0.05,0.5,0.95), na.rm = TRUE), c("#3B75AF","white","#C84630")),
    Age = circlize::colorRamp2(quantile(anno_df$Age, c(0.05,0.5,0.95), na.rm = TRUE), c("#E8F1FA","white","#D6604D")),
    Gender = c(Female = "#DDA0DD", Male = "#7FB3D5"),
    Stage = c(I = "#D9F0D3", II = "#ADDD8E", III = "#41AB5D", IV = "#005A32")[levels(anno_df$Stage)],
    Status = c(Alive = "#4DBBD5", Dead = "#E64B35")
  ),
  annotation_name_side = "left", show_annotation_name = TRUE
)

ht <- Heatmap(hm_mat, name = "Expression\nZ-score",
              col = circlize::colorRamp2(c(-2,0,2), c("#2B6CB0","white","#C53030")),
              top_annotation = top_anno,
              cluster_rows = TRUE, cluster_columns = FALSE,
              show_column_names = FALSE, show_row_names = TRUE,
              row_names_gp = grid::gpar(fontsize = 6),
              column_split = anno_df$Group,
              column_title = "Differentially expressed genes between High and Low MO-DDRscore groups",
              column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"))
pdf(file.path(OUT_DIR, "Fig2B_DEG_heatmap_with_clinical_annotation.pdf"), width = 10, height = 9, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE)
dev.off()

############################
# 6. Fig2C Hallmark GSVA/GSEA
############################

hallmark <- tryCatch(msigdbr(species = "Homo sapiens", category = "H"),
                     error = function(e) msigdbr(species = "Homo sapiens", collection = "H"))
term2gene <- hallmark %>% select(gs_name, gene_symbol) %>%
  mutate(gene_symbol = clean_gene(gene_symbol)) %>%
  filter(!is.na(gene_symbol), gene_symbol != "")

hallmark_sets <- split(term2gene$gene_symbol, term2gene$gs_name)
hallmark_sets <- lapply(hallmark_sets, function(x) intersect(unique(x), rownames(expr_log)))
hallmark_sets <- hallmark_sets[sapply(hallmark_sets, length) >= 5]

hallmark_gsva <- safe_ssgsea(expr_log, hallmark_sets)
hallmark_long <- as.data.frame(t(hallmark_gsva), check.names = FALSE) %>%
  mutate(Sample = rownames(.)) %>%
  pivot_longer(cols = -Sample, names_to = "Hallmark", values_to = "Score") %>%
  left_join(score_df %>% select(Sample, MO_DDRscore_group), by = "Sample")

hallmark_stat <- hallmark_long %>%
  group_by(Hallmark) %>%
  summarise(P = wilcox.test(Score ~ MO_DDRscore_group)$p.value,
            Median_High = median(Score[MO_DDRscore_group=="High"], na.rm = TRUE),
            Median_Low = median(Score[MO_DDRscore_group=="Low"], na.rm = TRUE),
            Diff = Median_High - Median_Low, .groups = "drop") %>%
  mutate(FDR = p.adjust(P, "BH")) %>% arrange(P)
save_csv(hallmark_long, file.path(OUT_DIR, "Fig2C_Hallmark_GSVA_scores.csv"))
save_csv(hallmark_stat, file.path(OUT_DIR, "Fig2C_Hallmark_GSVA_statistics.csv"))

gsva_show <- hallmark_stat %>% arrange(P) %>% slice_head(n = TOP_SHOW_N) %>%
  mutate(Pathway = gsub("^HALLMARK_", "", Hallmark),
         Pathway = gsub("_", " ", Pathway),
         Pathway = factor(Pathway, levels = rev(Pathway)),
         Direction = ifelse(Diff > 0, "High", "Low"))

p_gsva <- ggplot(gsva_show, aes(Diff, Pathway, color = Direction, size = -log10(FDR + 1e-300))) +
  geom_point(alpha = 0.9) +
  scale_color_manual(values = c(High = "#C84630", Low = "#3B75AF")) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text = element_text(color = "black")) +
  labs(x = "Median difference (High - Low)", y = NULL, size = "-log10(FDR)", title = "Hallmark GSVA")
safe_ggsave(file.path(OUT_DIR, "Fig2C_Hallmark_GSVA_dotplot.pdf"), p_gsva, 7.5, 6)

gene_rank <- deg$t; names(gene_rank) <- deg$Gene
gene_rank <- sort(gene_rank[is.finite(gene_rank)], decreasing = TRUE)
gsea_res <- clusterProfiler::GSEA(geneList = gene_rank, TERM2GENE = term2gene, pvalueCutoff = 1, verbose = FALSE)
gsea_df <- as.data.frame(gsea_res)
save_csv(gsea_df, file.path(OUT_DIR, "Fig2C_Hallmark_GSEA_results.csv"))
if (nrow(gsea_df) > 0) {
  gsea_show <- gsea_df %>% arrange(p.adjust) %>% slice_head(n = TOP_SHOW_N) %>%
    mutate(Pathway = gsub("^HALLMARK_", "", ID),
           Pathway = gsub("_", " ", Pathway),
           Pathway = factor(Pathway, levels = rev(Pathway)),
           Direction = ifelse(NES > 0, "High", "Low"))
  p_gsea <- ggplot(gsea_show, aes(NES, Pathway, size = -log10(p.adjust + 1e-300), color = Direction)) +
    geom_point(alpha = 0.9) +
    scale_color_manual(values = c(High = "#C84630", Low = "#3B75AF")) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(), axis.text = element_text(color = "black")) +
    labs(x = "NES", y = NULL, size = "-log10(FDR)", title = "Hallmark GSEA")
  safe_ggsave(file.path(OUT_DIR, "Fig2C_Hallmark_GSEA_dotplot.pdf"), p_gsea, 7.5, 6)
}

############################
# 7. Fig2D/E mutation and TMB: maftools
############################

make_maf_like <- function(file) {
  x <- fread(file, data.table = FALSE, check.names = FALSE)
  if (all(c("Hugo_Symbol","Tumor_Sample_Barcode") %in% colnames(x))) {
    x$Hugo_Symbol <- clean_gene(x$Hugo_Symbol)
    x$Tumor_Sample_Barcode <- patient_id(x$Tumor_Sample_Barcode)
    return(x)
  }
  stop("MAF_FILE must contain Hugo_Symbol and Tumor_Sample_Barcode for maftools.")
}

if (file.exists(MAF_FILE)) {
  maf_df <- make_maf_like(MAF_FILE) %>% filter(!is.na(Hugo_Symbol), Hugo_Symbol != "")
  maf_group <- score_df %>% transmute(Tumor_Sample_Barcode = Patient, MO_DDRscore_group = as.character(MO_DDRscore_group))
  maf_obj <- read.maf(maf = maf_df, clinicalData = maf_group, verbose = FALSE)
  saveRDS(maf_obj, file.path(OUT_DIR, "Fig2D_maftools_maf_object.rds"))
  
  all_ddr_genes <- unique(ddr_anno$Gene)
  ddr_maf_df <- maf_df %>% filter(Hugo_Symbol %in% all_ddr_genes)
  if (nrow(ddr_maf_df) > 0) {
    ddr_maf_obj <- read.maf(maf = ddr_maf_df, clinicalData = maf_group, verbose = FALSE)
    gene_summary <- getGeneSummary(ddr_maf_obj)
    save_csv(gene_summary, file.path(OUT_DIR, "Fig2D_DDR_mutation_gene_summary.csv"))
    top_mut_genes <- gene_summary %>% arrange(desc(MutatedSamples)) %>% slice_head(n = TOP_MUT_GENES) %>% pull(Hugo_Symbol)
    pdf(file.path(OUT_DIR, "Fig2D_DDR_mutation_oncoplot_maftools.pdf"), width = 10, height = 7, useDingbats = FALSE)
    oncoplot(maf = ddr_maf_obj, genes = top_mut_genes, clinicalFeatures = "MO_DDRscore_group",
             sortByAnnotation = TRUE, removeNonMutated = FALSE)
    dev.off()
  }
  
  tmb_df <- tmb(maf = maf_obj)
  tmb_df$Tumor_Sample_Barcode <- patient_id(tmb_df$Tumor_Sample_Barcode)
  tmb_df <- tmb_df %>%
    left_join(score_df %>% transmute(Tumor_Sample_Barcode = Patient, MO_DDRscore_group), by = "Tumor_Sample_Barcode") %>%
    filter(!is.na(MO_DDRscore_group))
  save_csv(tmb_df, file.path(OUT_DIR, "Fig2E_TMB_maftools.csv"))
  tmb_col <- intersect(c("total_perMB","total"), colnames(tmb_df))[1]
  if (!is.na(tmb_col)) {
    p_tmb <- box_compare(tmb_df, "MO_DDRscore_group", tmb_col, "Tumor mutational burden", tmb_col)
    safe_ggsave(file.path(OUT_DIR, "Fig2E_TMB_boxplot_maftools.pdf"), p_tmb, 4.5, 4.8)
  }
} else {
  message("MAF_FILE not found. Skip mutation/TMB.")
}

############################
# 8. Fig2F CNV: GISTIC2 and FGA/FGG/FGL official outputs
############################

if (all(file.exists(c(GISTIC_ALL_LESIONS, GISTIC_AMP_GENES, GISTIC_DEL_GENES)))) {
  gistic_obj <- tryCatch({
    readGistic(gisticAllLesionsFile = GISTIC_ALL_LESIONS,
               gisticAmpGenesFile = GISTIC_AMP_GENES,
               gisticDelGenesFile = GISTIC_DEL_GENES,
               gisticScoresFile = if (file.exists(GISTIC_SCORES)) GISTIC_SCORES else NULL,
               isTCGA = TRUE)
  }, error = function(e) {
    message("readGistic failed: ", e$message); NULL
  })
  if (!is.null(gistic_obj)) {
    saveRDS(gistic_obj, file.path(OUT_DIR, "Fig2F_GISTIC_object.rds"))
    pdf(file.path(OUT_DIR, "Fig2F_GISTIC_chromosome_plot.pdf"), width = 10, height = 5, useDingbats = FALSE)
    gisticChromPlot(gistic = gistic_obj, markBands = "all")
    dev.off()
    pdf(file.path(OUT_DIR, "Fig2F_GISTIC_bubble_plot.pdf"), width = 8, height = 6, useDingbats = FALSE)
    gisticBubblePlot(gistic = gistic_obj)
    dev.off()
  }
} else {
  message("GISTIC2 files not found. Skip official GISTIC CNV plots.")
}

if (file.exists(FGA_FILE)) {
  fga_df <- fread(FGA_FILE, data.table = FALSE, check.names = FALSE)
  if ("Sample" %in% colnames(fga_df)) fga_df$Patient <- patient_id(fga_df$Sample)
  if (!"Patient" %in% colnames(fga_df)) stop("FGA_FILE must contain Patient or Sample column.")
  fga_cols <- intersect(c("FGA","FGG","FGL"), colnames(fga_df))
  if (length(fga_cols) > 0) {
    fga_long <- fga_df %>%
      mutate(Patient = as.character(Patient)) %>%
      left_join(score_df %>% select(Patient, MO_DDRscore_group), by = "Patient") %>%
      filter(!is.na(MO_DDRscore_group)) %>%
      pivot_longer(cols = all_of(fga_cols), names_to = "Metric", values_to = "Value")
    save_csv(fga_long, file.path(OUT_DIR, "Fig2F_FGA_FGG_FGL_long.csv"))
    p_fga <- ggplot(fga_long, aes(MO_DDRscore_group, Value, fill = MO_DDRscore_group)) +
      geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
      geom_jitter(width = 0.16, size = 0.6, alpha = 0.45) +
      stat_compare_means(method = "wilcox.test", label = "p.signif", size = 3.2) +
      facet_wrap(~Metric, scales = "free_y", nrow = 1) +
      scale_fill_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
      theme_bw(base_size = 11) +
      theme(panel.grid = element_blank(), axis.text = element_text(color = "black"), legend.position = "none") +
      labs(x = NULL, y = "Fraction", title = "Genome alteration fractions")
    safe_ggsave(file.path(OUT_DIR, "Fig2F_FGA_FGG_FGL_boxplot.pdf"), p_fga, 8.5, 4.2)
  }
} else {
  message("FGA/FGG/FGL file not found. Skip FGA comparison.")
}

############################################################
# 9. ESTIMATE official workflow
############################################################

estimate_ok <- requireNamespace("estimate", quietly = TRUE)

write_estimate_gct <- function(expr_mat, file) {
  
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "numeric"
  
  expr_mat <- expr_mat[!is.na(rownames(expr_mat)) & rownames(expr_mat) != "", , drop = FALSE]
  expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), , drop = FALSE]
  
  expr_mat[!is.finite(expr_mat)] <- 0
  
  gct_df <- data.frame(
    Name = rownames(expr_mat),
    Description = rownames(expr_mat),
    expr_mat,
    check.names = FALSE
  )
  
  con <- file(file, open = "wt")
  writeLines("#1.2", con)
  writeLines(paste(nrow(gct_df), ncol(gct_df) - 2, sep = "\t"), con)
  close(con)
  
  write.table(
    gct_df,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    append = TRUE
  )
}

if (!estimate_ok) {
  
  message("Package 'estimate' is not installed. Skip official ESTIMATE workflow.")
  
} else {
  
  estimate_res <- tryCatch({
    
    estimate_in <- file.path(OUT_DIR, "ESTIMATE_input_expression.gct")
    estimate_common <- file.path(OUT_DIR, "ESTIMATE_common_genes.gct")
    estimate_out <- file.path(OUT_DIR, "ESTIMATE_scores.gct")
    
    # ESTIMATE 用表达矩阵。这里用 log2(TPM+1)，保持和前面分析一致。
    est_mat <- expr_log[, score_df$Sample, drop = FALSE]
    
    write_estimate_gct(est_mat, estimate_in)
    
    estimate::filterCommonGenes(
      input.f = estimate_in,
      output.f = estimate_common,
      id = "GeneSymbol"
    )
    
    estimate::estimateScore(
      input.ds = estimate_common,
      output.ds = estimate_out,
      platform = "illumina"
    )
    
    est_score <- data.table::fread(
      estimate_out,
      skip = 2,
      data.table = FALSE,
      check.names = FALSE
    )
    
    colnames(est_score)[1] <- "Metric"
    
    if ("Description" %in% colnames(est_score)) {
      est_score$Description <- NULL
    }
    
    rownames(est_score) <- est_score$Metric
    est_score$Metric <- NULL
    
    est_score_t <- as.data.frame(t(est_score), check.names = FALSE)
    est_score_t$Sample <- rownames(est_score_t)
    
    est_long <- est_score_t %>%
      dplyr::left_join(
        score_df %>% dplyr::select(Sample, Patient, MO_DDRscore_group),
        by = "Sample"
      ) %>%
      dplyr::filter(!is.na(MO_DDRscore_group))
    
    estimate_cols <- intersect(
      c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity"),
      colnames(est_long)
    )
    
    if (length(estimate_cols) == 0) {
      stop("No ESTIMATE score columns found in output.")
    }
    
    est_long <- est_long %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(estimate_cols),
        names_to = "Metric",
        values_to = "Value"
      )
    
    save_csv(
      est_long,
      file.path(OUT_DIR, "Fig2G_ESTIMATE_scores_long.csv")
    )
    
    p_est <- ggplot(est_long, aes(MO_DDRscore_group, Value, fill = MO_DDRscore_group)) +
      geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
      geom_jitter(width = 0.16, size = 0.6, alpha = 0.45) +
      ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", size = 3.2) +
      facet_wrap(~Metric, scales = "free_y", nrow = 1) +
      scale_fill_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none"
      ) +
      labs(
        x = NULL,
        y = "Score",
        title = "ESTIMATE scores"
      )
    
    safe_ggsave(
      file.path(OUT_DIR, "Fig2G_ESTIMATE_boxplot.pdf"),
      p_est,
      9,
      4.2
    )
    
    TRUE
    
  }, error = function(e) {
    
    message("Skip ESTIMATE workflow because of error: ", e$message)
    FALSE
  })
}
############################
# 10. Fig2H immune checkpoint expression
############################

if (file.exists(ICI_GENE_FILE)) {
  ici_raw <- fread(ICI_GENE_FILE, data.table = FALSE, check.names = FALSE)
  checkpoint_genes <- unique(na.omit(clean_gene(ici_raw[[1]])))
} else {
  checkpoint_genes <- c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","TIGIT","HAVCR2",
                        "IDO1","IDO2","CD276","CD80","CD86","CD70","CD27","CD28",
                        "CD40","CD40LG","TNFRSF4","TNFSF4","TNFRSF9","TNFSF9",
                        "BTLA","VSIR","ADORA2A","CD160","ICOS","ICOSLG")
  message("ICI_GENE_FILE not found. Use common checkpoint gene list.")
}
checkpoint_genes <- intersect(unique(clean_gene(checkpoint_genes)), rownames(expr_log))

if (length(checkpoint_genes) >= 2) {
  ck_long <- as.data.frame(t(expr_log[checkpoint_genes, score_df$Sample, drop = FALSE]), check.names = FALSE) %>%
    mutate(Sample = rownames(.)) %>%
    pivot_longer(cols = -Sample, names_to = "Checkpoint", values_to = "Expression") %>%
    left_join(score_df %>% select(Sample, Patient, MO_DDRscore_group), by = "Sample")
  ck_stat <- ck_long %>%
    group_by(Checkpoint) %>%
    summarise(P = wilcox.test(Expression ~ MO_DDRscore_group)$p.value,
              Median_High = median(Expression[MO_DDRscore_group=="High"], na.rm = TRUE),
              Median_Low = median(Expression[MO_DDRscore_group=="Low"], na.rm = TRUE),
              Diff = Median_High - Median_Low, .groups = "drop") %>%
    mutate(FDR = p.adjust(P, "BH")) %>% arrange(P)
  save_csv(ck_long, file.path(OUT_DIR, "Fig2H_Checkpoint_expression_long.csv"))
  save_csv(ck_stat, file.path(OUT_DIR, "Fig2H_Checkpoint_expression_statistics.csv"))
  ck_show <- ck_stat %>% arrange(P) %>% slice_head(n = min(TOP_CHECKPOINTS, n())) %>% pull(Checkpoint)
  p_ck <- ck_long %>%
    filter(Checkpoint %in% ck_show) %>%
    mutate(Checkpoint = factor(Checkpoint, levels = rev(ck_show))) %>%
    ggplot(aes(MO_DDRscore_group, Expression, fill = MO_DDRscore_group)) +
    geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.16, size = 0.45, alpha = 0.45) +
    stat_compare_means(method = "wilcox.test", label = "p.signif", size = 3) +
    facet_wrap(~Checkpoint, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), axis.text = element_text(color = "black"), legend.position = "none") +
    labs(x = NULL, y = "log2(TPM+1)", title = "Immune checkpoint expression")
  safe_ggsave(file.path(OUT_DIR, "Fig2H_Checkpoint_expression_boxplot.pdf"), p_ck, 11, 7)
}

############################
# 11. Summary
############################

summary_df <- data.frame(
  Item = c("Samples","Low_group","High_group","DDR_pathway_count","DDR_gene_count",
           "DEG_up_in_high","DEG_down_in_high","MAF_file_exists",
           "GISTIC_files_exist","FGA_file_exists","Immune_GMT_exists","Checkpoint_genes_matched"),
  Value = c(nrow(score_df), sum(score_df$MO_DDRscore_group == "Low"), sum(score_df$MO_DDRscore_group == "High"),
            length(ddr_sets), length(unique(ddr_anno$Gene)),
            sum(deg$Regulation == "Up in High"), sum(deg$Regulation == "Down in High"),
            file.exists(MAF_FILE), all(file.exists(c(GISTIC_ALL_LESIONS, GISTIC_AMP_GENES, GISTIC_DEL_GENES))),
            file.exists(FGA_FILE), file.exists(IMMUNE_GMT_FILE), length(checkpoint_genes))
)

save_csv(summary_df, file.path(OUT_DIR, "Fig2_analysis_summary.csv"))

cat("\nDone. Outputs saved to:\n", OUT_DIR, "\n")
cat("Recommended main panels:\n")
cat("A: Fig2A_DDR_pathway_activity_boxplot.pdf or heatmap\n")
cat("B: Fig2B_DEG_heatmap_with_clinical_annotation.pdf\n")
cat("C: Fig2C_Hallmark_GSEA_dotplot.pdf / Fig2C_Hallmark_GSVA_dotplot.pdf\n")
cat("D: Fig2D_DDR_mutation_oncoplot_maftools.pdf\n")
cat("E: Fig2E_TMB_boxplot_maftools.pdf\n")
cat("F: Fig2F_GISTIC_* or Fig2F_FGA_FGG_FGL_boxplot.pdf\n")
cat("G: Fig2G_ESTIMATE_boxplot.pdf / Fig2G_Immune_ssGSEA_boxplot.pdf\n")
cat("H: Fig2H_Checkpoint_expression_boxplot.pdf\n")
