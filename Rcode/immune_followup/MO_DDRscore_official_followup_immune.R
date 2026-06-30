############################################################
# Official follow-up immune analysis after basic IOBR results
# MO-DDRscore high/low groups
#
# This script avoids homemade gene sets:
# - Immune deconvolution results are read from IOBR outputs.
# - Pathway/module scores are computed only from MSigDB gene sets via msigdbr.
# - ssGSEA is run with Bioconductor GSVA.
# - TIDE is handled as an official web-server workflow: export input and merge
#   downloaded TIDE output if available.
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260513)

############################
# 0. Packages
############################

cran_pkgs <- c(
  "data.table", "dplyr", "tidyr", "ggplot2",
  "tibble", "stringr", "msigdbr", "pheatmap"
)

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_pkgs <- c("GSVA")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(stringr)
  library(msigdbr)
  library(GSVA)
  library(pheatmap)
})

############################
# 1. Paths
############################

PROJECT_DIR <- "D:/R_workspace/评分/AD_DDR_project"
PROC_DIR <- file.path(PROJECT_DIR, "01_processed")
FIG2_DIR <- file.path(PROJECT_DIR, "02_Figure2_MultiOmics_Immune")
BASIC_DIR <- file.path(FIG2_DIR, "Immune_basic_MO_DDRscore")
OUT_DIR <- file.path(FIG2_DIR, "Immune_official_followup_MO_DDRscore")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

EXPR_FILE <- file.path(PROC_DIR, "LUAD_tcga_expr_tpm_matrix.rds")
ANNOT_FILE <- file.path(BASIC_DIR, "Immune_sample_annotation.csv")
IMMUNE_LONG_FILE <- file.path(BASIC_DIR, "Immune_deconvolution_long.csv")
IMMUNE_STAT_FILE <- file.path(BASIC_DIR, "Immune_deconvolution_group_comparison.csv")
MARKER_LONG_FILE <- file.path(BASIC_DIR, "Immune_marker_expression_long.csv")
MARKER_STAT_FILE <- file.path(BASIC_DIR, "Immune_marker_group_comparison.csv")

TMB_CANDIDATES <- c(
  file.path(FIG2_DIR, "Fig2E_TMB_maftools_table.csv"),
  file.path(PROJECT_DIR, "02_Figure2_official_multiomics", "Fig2E_TMB_maftools.csv"),
  file.path(PROJECT_DIR, "03_Figure3_multiomics", "TMB_mutation_count_table.csv"),
  file.path(PROJECT_DIR, "table", "Fig3_multiomics_GDSC", "Fig3D2_DDR_mutation_burden_table.csv")
)

MAF_RDS_CANDIDATES <- c(
  file.path(PROJECT_DIR, "02_Figure2_official_multiomics", "Fig2D_maftools_maf_object.rds")
)

STEMNESS_CANDIDATES <- c(
  file.path(PROC_DIR, "LUAD_stemness_scores.csv"),
  file.path(FIG2_DIR, "LUAD_stemness_scores.csv"),
  file.path(OUT_DIR, "LUAD_stemness_scores.csv")
)

TIDE_RESULT_CANDIDATES <- c(
  file.path(OUT_DIR, "TIDE_result.csv"),
  file.path(OUT_DIR, "TIDE_result.txt"),
  file.path(FIG2_DIR, "TIDE_result.csv"),
  file.path(FIG2_DIR, "TIDE_result.txt")
)

############################
# 2. Helper functions
############################

save_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file)
}

clean_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub(" ///.*$| //.*$|;.*$|,.*$|\\s*\\([^\\)]*\\)$", "", x)
  x[x %in% c("", "NA", "---", "NULL", "N/A")] <- NA
  x
}

patient_id <- function(x) substr(gsub("\\.", "-", x), 1, 12)
safe_num <- function(x) suppressWarnings(as.numeric(x))

read_csv_if_exists <- function(files) {
  hit <- files[file.exists(files)][1]
  if (length(hit) == 0 || is.na(hit)) return(NULL)
  message("Reading: ", hit)
  data.table::fread(hit, data.table = FALSE, check.names = FALSE)
}

find_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

normalize_sample_col <- function(df) {
  id_col <- find_col(df, c(
    "Sample", "sample", "ID", "id", "Tumor_Sample_Barcode",
    "Tumor.Sample.Barcode", "Patient", "patient", "submitter_id",
    "V1", "X", "...1"
  ))
  if (!is.na(id_col)) {
    colnames(df)[colnames(df) == id_col] <- "Sample"
    df$Sample <- gsub("\\.", "-", df$Sample)
    df$Patient <- patient_id(df$Sample)
  }
  df
}

safe_wilcox <- function(x, g) {
  keep <- is.finite(x) & !is.na(g)
  x <- x[keep]
  g <- droplevels(factor(g[keep], levels = c("Low", "High")))
  if (length(unique(g)) < 2 || length(unique(x)) < 2) return(NA_real_)
  tryCatch(wilcox.test(x ~ g)$p.value, error = function(e) NA_real_)
}

group_compare <- function(df, value_col, group_col = "MO_DDRscore_group") {
  x <- safe_num(df[[value_col]])
  g <- factor(df[[group_col]], levels = c("Low", "High"))
  low <- x[g == "Low"]
  high <- x[g == "High"]
  data.frame(
    N_Low = sum(is.finite(low)),
    N_High = sum(is.finite(high)),
    Median_Low = median(low, na.rm = TRUE),
    Median_High = median(high, na.rm = TRUE),
    Delta_High_minus_Low = median(high, na.rm = TRUE) - median(low, na.rm = TRUE),
    P = safe_wilcox(x, g)
  )
}

safe_spearman <- function(x, y) {
  x <- safe_num(x)
  y <- safe_num(y)
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 10 || length(unique(x)) < 3 || length(unique(y)) < 3) {
    return(c(Rho = NA_real_, P = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
  c(Rho = unname(ct$estimate), P = ct$p.value)
}

zscore <- function(x) {
  x <- safe_num(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

plot_box_facet <- function(df, feature_col, value_col, title, file, ncol = 4,
                           width = 10, height = 7) {
  if (nrow(df) == 0) return(NULL)
  p <- ggplot(df, aes(x = MO_DDRscore_group, y = .data[[value_col]],
                      fill = MO_DDRscore_group)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 0.45, alpha = 0.35) +
    scale_fill_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
    facet_wrap(stats::as.formula(paste("~", feature_col)), scales = "free_y", ncol = ncol) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text = element_text(color = "black"),
      strip.text = element_text(face = "bold", size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(x = NULL, y = value_col, title = title)
  ggsave(file, p, width = width, height = height, useDingbats = FALSE)
}

run_ssgsea <- function(expr_log, gene_sets) {
  # Compatible with both newer and older GSVA APIs.
  if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    param <- GSVA::ssgseaParam(expr_log, gene_sets, minSize = 5, maxSize = 500,
                               normalize = TRUE)
    out <- GSVA::gsva(param, verbose = FALSE)
  } else {
    out <- GSVA::gsva(
      expr = expr_log,
      gset.idx.list = gene_sets,
      method = "ssgsea",
      kcdf = "Gaussian",
      min.sz = 5,
      max.sz = 500,
      ssgsea.norm = TRUE,
      verbose = FALSE
    )
  }
  as.matrix(out)
}

############################
# 3. Load expression and annotation
############################

anno <- data.table::fread(ANNOT_FILE, data.table = FALSE, check.names = FALSE)
anno$Sample <- gsub("\\.", "-", anno$Sample)
if (!"Patient" %in% colnames(anno)) {
  anno$Patient <- patient_id(anno$Sample)
}
anno <- anno %>%
  mutate(
    Patient = as.character(Patient),
    MO_DDRscore_raw = safe_num(MO_DDRscore_raw),
    MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High"))
  )

expr <- readRDS(EXPR_FILE)
expr <- as.matrix(expr)
storage.mode(expr) <- "numeric"
rownames(expr) <- clean_gene(rownames(expr))
colnames(expr) <- gsub("\\.", "-", colnames(expr))
expr <- expr[!is.na(rownames(expr)) & rownames(expr) != "", , drop = FALSE]

# Average duplicated gene symbols.
gene_count <- table(rownames(expr))
expr <- rowsum(expr, group = rownames(expr), reorder = FALSE)
expr <- sweep(expr, 1, as.numeric(gene_count[rownames(expr)]), "/")
expr[!is.finite(expr)] <- 0
expr[expr < 0] <- 0

samples <- intersect(anno$Sample, colnames(expr))
anno <- anno %>% filter(Sample %in% samples) %>% arrange(match(Sample, samples))
expr_tpm <- expr[, anno$Sample, drop = FALSE]
expr_log <- log2(expr_tpm + 1)

cat("Samples:", ncol(expr_tpm), "\n")
print(table(anno$MO_DDRscore_group))

############################
# 4. Read basic IOBR outputs
############################

immune_long <- data.table::fread(IMMUNE_LONG_FILE, data.table = FALSE, check.names = FALSE) %>%
  mutate(Sample = gsub("\\.", "-", Sample), Score = safe_num(Score)) %>%
  left_join(anno %>% select(Sample, Patient, MO_DDRscore_raw, MO_DDRscore_group),
            by = "Sample", suffix = c("", ".anno")) %>%
  mutate(
    Patient = dplyr::coalesce(Patient, Patient.anno),
    MO_DDRscore_raw = dplyr::coalesce(MO_DDRscore_raw, MO_DDRscore_raw.anno),
    MO_DDRscore_group = dplyr::coalesce(MO_DDRscore_group, MO_DDRscore_group.anno)
  ) %>%
  select(-matches("\\.anno$"))

immune_stat <- data.table::fread(IMMUNE_STAT_FILE, data.table = FALSE, check.names = FALSE) %>%
  mutate(
    P = safe_num(P),
    FDR = safe_num(FDR),
    Delta_High_minus_Low = safe_num(Delta_High_minus_Low),
    Direction = ifelse(Delta_High_minus_Low > 0, "High_up", "High_down")
  )

marker_long <- data.table::fread(MARKER_LONG_FILE, data.table = FALSE, check.names = FALSE) %>%
  mutate(Sample = gsub("\\.", "-", Sample), Expression = safe_num(Expression)) %>%
  left_join(anno %>% select(Sample, Patient, MO_DDRscore_raw, MO_DDRscore_group),
            by = "Sample", suffix = c("", ".anno")) %>%
  mutate(
    Patient = dplyr::coalesce(Patient, Patient.anno),
    MO_DDRscore_raw = dplyr::coalesce(MO_DDRscore_raw, MO_DDRscore_raw.anno),
    MO_DDRscore_group = dplyr::coalesce(MO_DDRscore_group, MO_DDRscore_group.anno)
  ) %>%
  select(-matches("\\.anno$"))

marker_stat <- data.table::fread(MARKER_STAT_FILE, data.table = FALSE, check.names = FALSE) %>%
  mutate(
    P = safe_num(P),
    FDR = safe_num(FDR),
    Delta_High_minus_Low = safe_num(Delta_High_minus_Low),
    Direction = ifelse(Delta_High_minus_Low > 0, "High_up", "High_down")
  )

############################
# 5. Official MSigDB gene sets
############################

cat("Loading MSigDB gene sets from msigdbr...\n")

msig_all <- tryCatch(
  msigdbr::msigdbr(species = "Homo sapiens"),
  error = function(e) {
    msigdbr::msigdbr(db_species = "HS", species = "Homo sapiens")
  }
)

gene_col <- if ("gene_symbol" %in% colnames(msig_all)) {
  "gene_symbol"
} else if ("human_gene_symbol" %in% colnames(msig_all)) {
  "human_gene_symbol"
} else {
  stop("Cannot identify gene symbol column in msigdbr output.")
}

msig_all <- msig_all %>%
  mutate(
    gs_name = as.character(gs_name),
    Gene = clean_gene(.data[[gene_col]])
  ) %>%
  filter(!is.na(Gene), Gene %in% rownames(expr_log))

# These are official MSigDB sets selected because your basic analysis showed:
# B-cell/plasma-cell depletion, myeloid/neutrophil remodeling,
# antigen-presentation/checkpoint/CXCL9-11 activation, and lower IPS.
module_patterns <- c(
  HALLMARK_IFNG = "^HALLMARK_INTERFERON_GAMMA_RESPONSE$",
  HALLMARK_IFNA = "^HALLMARK_INTERFERON_ALPHA_RESPONSE$",
  HALLMARK_INFLAMMATION = "^HALLMARK_INFLAMMATORY_RESPONSE$",
  HALLMARK_TGF_BETA = "^HALLMARK_TGF_BETA_SIGNALING$",
  HALLMARK_EMT = "^HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION$",
  HALLMARK_IL6_JAK_STAT3 = "^HALLMARK_IL6_JAK_STAT3_SIGNALING$",
  HALLMARK_COMPLEMENT = "^HALLMARK_COMPLEMENT$",
  GO_B_CELL_ACTIVATION = "^(GO|GOBP)_B_CELL_ACTIVATION$",
  GO_BCR_SIGNALING = "^(GO|GOBP)_B_CELL_RECEPTOR_SIGNALING_PATHWAY$",
  GO_PLASMA_CELL_DIFFERENTIATION = "^(GO|GOBP)_PLASMA_CELL_DIFFERENTIATION$",
  GO_ANTIGEN_PRESENTATION = "^(GO|GOBP)_ANTIGEN_PROCESSING_AND_PRESENTATION$",
  GO_MHC_I_PRESENTATION = "MHC_CLASS_I.*ANTIGEN_PROCESSING|ANTIGEN_PROCESSING.*MHC_CLASS_I",
  GO_RESPONSE_TO_IFNG = "^(GO|GOBP)_RESPONSE_TO_INTERFERON_GAMMA$",
  GO_NEUTROPHIL_ACTIVATION = "^(GO|GOBP)_NEUTROPHIL_ACTIVATION$",
  GO_MYELOID_LEUKOCYTE_ACTIVATION = "^(GO|GOBP)_MYELOID_LEUKOCYTE_ACTIVATION$",
  GO_MONONUCLEAR_CELL_MIGRATION = "^(GO|GOBP)_MONONUCLEAR_CELL_MIGRATION$",
  GO_CYTOKINE_MEDIATED_SIGNALING = "^(GO|GOBP)_CYTOKINE_MEDIATED_SIGNALING_PATHWAY$",
  REACTOME_INTERFERON_GAMMA_SIGNALING = "^REACTOME_INTERFERON_GAMMA_SIGNALING$",
  REACTOME_ANTIGEN_PRESENTATION = "REACTOME_.*ANTIGEN.*PRESENTATION",
  REACTOME_BCR_SIGNALING = "REACTOME_.*B_CELL_RECEPTOR|REACTOME_.*BCR",
  REACTOME_PD1_SIGNALING = "REACTOME_.*PD_1|REACTOME_.*PD1"
)

selected_sets <- bind_rows(lapply(names(module_patterns), function(label) {
  pat <- module_patterns[[label]]
  msig_all %>%
    filter(grepl(pat, gs_name, ignore.case = FALSE)) %>%
    mutate(Module = label)
}))

if (nrow(selected_sets) == 0) {
  stop("No selected MSigDB gene sets were matched. Check msigdbr version / set names.")
}

matched_set_summary <- selected_sets %>%
  distinct(Module, gs_name, Gene) %>%
  count(Module, gs_name, name = "N_genes_matched") %>%
  arrange(Module, gs_name)

save_csv(matched_set_summary, file.path(OUT_DIR, "Official_MSigDB_selected_sets.csv"))

missing_modules <- setdiff(names(module_patterns), unique(selected_sets$Module))
if (length(missing_modules) > 0) {
  writeLines(missing_modules, file.path(OUT_DIR, "Official_MSigDB_missing_modules.txt"))
}

msig_list <- split(selected_sets$Gene, paste(selected_sets$Module, selected_sets$gs_name, sep = "__"))
msig_list <- lapply(msig_list, unique)
msig_list <- msig_list[lengths(msig_list) >= 5]

cat("MSigDB sets used:", length(msig_list), "\n")

############################
# 6. ssGSEA for official MSigDB sets
############################

ssgsea_mat <- run_ssgsea(expr_log, msig_list)

ssgsea_long <- as.data.frame(ssgsea_mat, check.names = FALSE) %>%
  rownames_to_column("SetID") %>%
  tidyr::separate(SetID, into = c("Module", "MSigDB_Set"), sep = "__", remove = FALSE) %>%
  pivot_longer(
    cols = -c(SetID, Module, MSigDB_Set),
    names_to = "Sample",
    values_to = "Score"
  ) %>%
  left_join(anno %>% select(Sample, Patient, MO_DDRscore_raw, MO_DDRscore_group), by = "Sample") %>%
  mutate(Score = safe_num(Score))

save_csv(ssgsea_long, file.path(OUT_DIR, "Official_MSigDB_ssGSEA_long.csv"))

ssgsea_stat <- ssgsea_long %>%
  group_by(Module, MSigDB_Set) %>%
  group_modify(~ group_compare(.x, "Score")) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(P, method = "BH"),
         Direction = ifelse(Delta_High_minus_Low > 0, "High_up", "High_down")) %>%
  arrange(FDR, P)

save_csv(ssgsea_stat, file.path(OUT_DIR, "Official_MSigDB_ssGSEA_group_comparison.csv"))

ssgsea_cor <- ssgsea_long %>%
  group_by(Module, MSigDB_Set) %>%
  summarise(
    N = n(),
    Rho = safe_spearman(Score, MO_DDRscore_raw)["Rho"],
    P = safe_spearman(Score, MO_DDRscore_raw)["P"],
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(P, method = "BH"),
         Direction = ifelse(Rho > 0, "Positive", "Negative")) %>%
  arrange(FDR, P)

save_csv(ssgsea_cor, file.path(OUT_DIR, "Official_MSigDB_ssGSEA_vs_MO_DDRscore_spearman.csv"))

plot_box_facet(
  ssgsea_long %>% mutate(Label = Module),
  "Label", "Score",
  "Official MSigDB ssGSEA modules by MO-DDRscore group",
  file.path(OUT_DIR, "Fig_official_MSigDB_ssGSEA_boxplots.pdf"),
  ncol = 4, width = 12, height = 8
)

dot_df <- ssgsea_stat %>%
  filter(is.finite(FDR)) %>%
  mutate(
    Label = paste(Module, MSigDB_Set, sep = "\n"),
    Label = factor(Label, levels = rev(Label)),
    LogFDR = -log10(pmax(FDR, 1e-300))
  )

if (nrow(dot_df) > 0) {
  p <- ggplot(dot_df, aes(x = Delta_High_minus_Low, y = Label)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
    geom_point(aes(size = LogFDR, color = Direction), alpha = 0.9) +
    scale_color_manual(values = c(High_down = "#3B75AF", High_up = "#C84630")) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text = element_text(color = "black")) +
    labs(
      x = "Median difference: High - Low",
      y = NULL,
      color = "Direction",
      size = "-log10(FDR)",
      title = "Official MSigDB immune pathway differences"
    )
  ggsave(file.path(OUT_DIR, "Fig_official_MSigDB_ssGSEA_delta_dotplot.pdf"),
         p, width = 9, height = max(5, 0.24 * nrow(dot_df)), useDingbats = FALSE)
}

############################
# 7. Official deconvolution-derived axes from your basic results
############################

# These are not new gene sets. They summarize official deconvolution outputs
# already generated by IOBR.
axis_features <- list(
  Global_TME = c("ESTIMATEScore_estimate", "ImmuneScore_estimate",
                 "StromalScore_estimate", "TumorPurity_estimate"),
  B_Plasma_axis = c("Bcells_EPIC", "B_lineage_MCPcounter",
                    "B_cells_quantiseq", "B_cell_TIMER",
                    "Plasma_cells_xCell", "Class-switched_memory_B-cells_xCell"),
  Myeloid_Neutrophil_axis = c("Monocytic_lineage_MCPcounter",
                              "Neutrophils_MCPcounter", "Neutrophils_quantiseq",
                              "Neutrophil_TIMER", "Myeloid_dendritic_cells_MCPcounter"),
  Macrophage_axis = c("Macrophages_M1_quantiseq", "Macrophages_M2_quantiseq",
                      "Macrophages_EPIC", "Macrophage_TIMER"),
  IPS_axis = c("IPS_IPS", "AZ_IPS", "MHC_IPS", "CP_IPS", "EC_IPS", "SC_IPS")
)

immune_wide <- immune_long %>%
  mutate(FeatureID = Feature) %>%
  select(Sample, Patient, MO_DDRscore_raw, MO_DDRscore_group, FeatureID, Score) %>%
  pivot_wider(names_from = FeatureID, values_from = Score)

axis_df <- anno %>% select(Sample, Patient, MO_DDRscore_raw, MO_DDRscore_group)

for (ax in names(axis_features)) {
  feats <- intersect(axis_features[[ax]], colnames(immune_wide))
  if (length(feats) == 0) next
  tmp <- immune_wide %>%
    select(Sample, all_of(feats))
  tmp[, feats] <- lapply(tmp[, feats, drop = FALSE], zscore)
  tmp[[ax]] <- rowMeans(tmp[, feats, drop = FALSE], na.rm = TRUE)
  axis_df <- axis_df %>% left_join(tmp %>% select(Sample, all_of(ax)), by = "Sample")
}

if (all(c("Macrophages_M1_quantiseq", "Macrophages_M2_quantiseq") %in% colnames(immune_wide))) {
  axis_df <- axis_df %>%
    left_join(
      immune_wide %>%
        transmute(
          Sample,
          quanTIseq_M1_M2_ratio =
            (safe_num(Macrophages_M1_quantiseq) + 1e-6) /
            (safe_num(Macrophages_M2_quantiseq) + 1e-6)
        ),
      by = "Sample"
    )
}

axis_long <- axis_df %>%
  pivot_longer(
    cols = -c(Sample, Patient, MO_DDRscore_raw, MO_DDRscore_group),
    names_to = "Axis",
    values_to = "Score"
  ) %>%
  filter(is.finite(Score))

axis_stat <- axis_long %>%
  group_by(Axis) %>%
  group_modify(~ group_compare(.x, "Score")) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(P, method = "BH"),
         Direction = ifelse(Delta_High_minus_Low > 0, "High_up", "High_down")) %>%
  arrange(FDR, P)

axis_cor <- axis_long %>%
  group_by(Axis) %>%
  summarise(
    N = n(),
    Rho = safe_spearman(Score, MO_DDRscore_raw)["Rho"],
    P = safe_spearman(Score, MO_DDRscore_raw)["P"],
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(P, method = "BH"),
         Direction = ifelse(Rho > 0, "Positive", "Negative")) %>%
  arrange(FDR, P)

save_csv(axis_long, file.path(OUT_DIR, "Official_deconvolution_result_driven_axes_long.csv"))
save_csv(axis_stat, file.path(OUT_DIR, "Official_deconvolution_result_driven_axes_group_comparison.csv"))
save_csv(axis_cor, file.path(OUT_DIR, "Official_deconvolution_result_driven_axes_vs_MO_DDRscore_spearman.csv"))

plot_box_facet(
  axis_long, "Axis", "Score",
  "Result-driven axes from official deconvolution outputs",
  file.path(OUT_DIR, "Fig_official_deconvolution_axes_boxplots.pdf"),
  ncol = 3, width = 9, height = 5.5
)

############################
# 8. TMB integration
############################

tmb_df <- read_csv_if_exists(TMB_CANDIDATES)
if (!is.null(tmb_df)) {
  tmb_df <- normalize_sample_col(tmb_df)
  tmb_col <- find_col(tmb_df, c("total_perMB", "TMB", "tmb", "nonsynonymous_perMB",
                                "Mutation_Burden", "total"))
  if (!is.na(tmb_col)) {
    tmb_df <- tmb_df %>%
      mutate(
        TMB_value = safe_num(.data[[tmb_col]]),
        TMB_log1p = log10(TMB_value + 1)
      ) %>%
      left_join(
        anno %>% select(Patient, MO_DDRscore_raw, MO_DDRscore_group),
        by = "Patient",
        suffix = c("", ".anno")
      ) %>%
      mutate(
        MO_DDRscore_raw = dplyr::coalesce(
          if ("MO_DDRscore_raw" %in% colnames(.)) safe_num(MO_DDRscore_raw) else NA_real_,
          if ("MO_DDRscore_raw.anno" %in% colnames(.)) safe_num(MO_DDRscore_raw.anno) else NA_real_
        ),
        MO_DDRscore_group = dplyr::coalesce(
          if ("MO_DDRscore_group" %in% colnames(.)) as.character(MO_DDRscore_group) else NA_character_,
          if ("MO_DDRscore_group.anno" %in% colnames(.)) as.character(MO_DDRscore_group.anno) else NA_character_
        ),
        MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High"))
      ) %>%
      select(-matches("\\.anno$"))

    tmb_stat <- group_compare(tmb_df, "TMB_value")
    tmb_cor <- safe_spearman(tmb_df$TMB_value, tmb_df$MO_DDRscore_raw)

    save_csv(tmb_df, file.path(OUT_DIR, "Official_TMB_merged.csv"))
    save_csv(tmb_stat, file.path(OUT_DIR, "Official_TMB_group_comparison.csv"))
    save_csv(data.frame(Rho = tmb_cor["Rho"], P = tmb_cor["P"]),
             file.path(OUT_DIR, "Official_TMB_vs_MO_DDRscore_spearman.csv"))

    p <- ggplot(tmb_df, aes(x = MO_DDRscore_group, y = TMB_log1p,
                            fill = MO_DDRscore_group)) +
      geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.85) +
      geom_jitter(width = 0.15, size = 0.55, alpha = 0.35) +
      scale_fill_manual(values = c(Low = "#3B75AF", High = "#C84630")) +
      theme_bw(base_size = 11) +
      theme(panel.grid = element_blank(), legend.position = "none",
            axis.text = element_text(color = "black"),
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      labs(x = NULL, y = "log10(TMB + 1)",
           title = paste0("TMB, Wilcoxon P = ", signif(tmb_stat$P, 3)))
    ggsave(file.path(OUT_DIR, "Fig_official_TMB_by_group.pdf"),
           p, width = 4.6, height = 4.1, useDingbats = FALSE)
  }
} else {
  writeLines("No TMB file found. TMB integration skipped.",
             file.path(OUT_DIR, "NOTE_TMB_missing.txt"))
}

############################
# 9. MATH from MAF object if VAF columns are available
############################

maf_rds <- MAF_RDS_CANDIDATES[file.exists(MAF_RDS_CANDIDATES)][1]
if (!is.na(maf_rds) && length(maf_rds) > 0) {
  maf_obj <- readRDS(maf_rds)
  maf_data <- NULL
  if (isS4(maf_obj) && "data" %in% slotNames(maf_obj)) {
    maf_data <- as.data.frame(maf_obj@data, check.names = FALSE)
  } else if (is.data.frame(maf_obj)) {
    maf_data <- maf_obj
  }

  if (!is.null(maf_data) && nrow(maf_data) > 0) {
    sample_col <- find_col(maf_data, c("Tumor_Sample_Barcode", "Sample", "sample"))
    if (!is.na(sample_col)) {
      colnames(maf_data)[colnames(maf_data) == sample_col] <- "Sample"
      maf_data$Patient <- patient_id(maf_data$Sample)

      if (all(c("t_alt_count", "t_ref_count") %in% colnames(maf_data))) {
        maf_data$VAF_for_MATH <- safe_num(maf_data$t_alt_count) /
          (safe_num(maf_data$t_alt_count) + safe_num(maf_data$t_ref_count))
      } else {
        vaf_col <- find_col(maf_data, c("VAF", "vaf", "TumorVAF", "i_TumorVAF_WU", "tumor_f"))
        if (!is.na(vaf_col)) maf_data$VAF_for_MATH <- safe_num(maf_data[[vaf_col]])
      }

      if ("VAF_for_MATH" %in% colnames(maf_data)) {
        maf_data$VAF_for_MATH[maf_data$VAF_for_MATH > 1] <-
          maf_data$VAF_for_MATH[maf_data$VAF_for_MATH > 1] / 100

        math_df <- maf_data %>%
          filter(is.finite(VAF_for_MATH), VAF_for_MATH > 0) %>%
          group_by(Patient) %>%
          summarise(
            N_mut_for_MATH = n(),
            Median_VAF = median(VAF_for_MATH, na.rm = TRUE),
            MAD_VAF = mad(VAF_for_MATH, na.rm = TRUE),
            MATH = 100 * MAD_VAF / Median_VAF,
            .groups = "drop"
          ) %>%
          filter(N_mut_for_MATH >= 5, is.finite(MATH)) %>%
          left_join(anno %>% select(Patient, MO_DDRscore_raw, MO_DDRscore_group), by = "Patient")

        if (nrow(math_df) > 0) {
          math_stat <- group_compare(math_df, "MATH")
          math_cor <- safe_spearman(math_df$MATH, math_df$MO_DDRscore_raw)

          save_csv(math_df, file.path(OUT_DIR, "Official_MATH_score.csv"))
          save_csv(math_stat, file.path(OUT_DIR, "Official_MATH_group_comparison.csv"))
          save_csv(data.frame(Rho = math_cor["Rho"], P = math_cor["P"]),
                   file.path(OUT_DIR, "Official_MATH_vs_MO_DDRscore_spearman.csv"))
        }
      }
    }
  }
} else {
  writeLines("No MAF RDS found. MATH calculation skipped.",
             file.path(OUT_DIR, "NOTE_MATH_missing.txt"))
}

############################
# 10. Stemness integration if public stemness table is available
############################

stem_df <- read_csv_if_exists(STEMNESS_CANDIDATES)
if (!is.null(stem_df)) {
  stem_df <- normalize_sample_col(stem_df)
  stem_col <- find_col(stem_df, c("mRNAsi", "RNAss", "DNAss", "Stemness",
                                  "stemness", "StemnessScore"))
  if (!is.na(stem_col)) {
    stem_df <- stem_df %>%
      mutate(Stemness_value = safe_num(.data[[stem_col]])) %>%
      left_join(
        anno %>% select(Patient, MO_DDRscore_raw, MO_DDRscore_group),
        by = "Patient",
        suffix = c("", ".anno")
      ) %>%
      mutate(
        MO_DDRscore_raw = dplyr::coalesce(
          if ("MO_DDRscore_raw" %in% colnames(.)) safe_num(MO_DDRscore_raw) else NA_real_,
          if ("MO_DDRscore_raw.anno" %in% colnames(.)) safe_num(MO_DDRscore_raw.anno) else NA_real_
        ),
        MO_DDRscore_group = dplyr::coalesce(
          if ("MO_DDRscore_group" %in% colnames(.)) as.character(MO_DDRscore_group) else NA_character_,
          if ("MO_DDRscore_group.anno" %in% colnames(.)) as.character(MO_DDRscore_group.anno) else NA_character_
        ),
        MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High"))
      ) %>%
      select(-matches("\\.anno$"))

    stem_stat <- group_compare(stem_df, "Stemness_value")
    stem_cor <- safe_spearman(stem_df$Stemness_value, stem_df$MO_DDRscore_raw)

    save_csv(stem_df, file.path(OUT_DIR, "Official_stemness_merged.csv"))
    save_csv(stem_stat, file.path(OUT_DIR, "Official_stemness_group_comparison.csv"))
    save_csv(data.frame(Rho = stem_cor["Rho"], P = stem_cor["P"]),
             file.path(OUT_DIR, "Official_stemness_vs_MO_DDRscore_spearman.csv"))
  }
} else {
  writeLines(
    c(
      "No public stemness table found.",
      "If you want this module, download mRNAsi/RNAss/DNAss from a public source,",
      "save Sample/Patient + mRNAsi/RNAss/DNAss columns as LUAD_stemness_scores.csv, and rerun."
    ),
    file.path(OUT_DIR, "NOTE_stemness_missing.txt")
  )
}

############################
# 11. TIDE official workflow support
############################

tide_tpm <- data.frame(Gene = rownames(expr_tpm), expr_tpm, check.names = FALSE)
save_csv(tide_tpm, file.path(OUT_DIR, "TIDE_upload_gene_by_sample_TPM.csv"))

tide_log <- data.frame(Gene = rownames(expr_log), expr_log, check.names = FALSE)
save_csv(tide_log, file.path(OUT_DIR, "TIDE_upload_gene_by_sample_log2TPMplus1.csv"))

writeLines(
  c(
    "Use the official TIDE web workflow. Do not approximate TIDE with a homemade score.",
    "Upload one of:",
    "1. TIDE_upload_gene_by_sample_TPM.csv",
    "2. TIDE_upload_gene_by_sample_log2TPMplus1.csv",
    "After downloading TIDE output, save it as TIDE_result.csv or TIDE_result.txt in this folder and rerun this script.",
    "The script will merge all numeric TIDE columns and compare High vs Low."
  ),
  file.path(OUT_DIR, "README_TIDE_official_upload.txt")
)

tide_res <- read_csv_if_exists(TIDE_RESULT_CANDIDATES)
if (!is.null(tide_res)) {
  tide_res <- normalize_sample_col(tide_res) %>%
    left_join(
      anno %>% select(Patient, MO_DDRscore_raw, MO_DDRscore_group),
      by = "Patient",
      suffix = c("", ".anno")
    ) %>%
    mutate(
      MO_DDRscore_raw = dplyr::coalesce(
        if ("MO_DDRscore_raw" %in% colnames(.)) safe_num(MO_DDRscore_raw) else NA_real_,
        if ("MO_DDRscore_raw.anno" %in% colnames(.)) safe_num(MO_DDRscore_raw.anno) else NA_real_
      ),
      MO_DDRscore_group = dplyr::coalesce(
        if ("MO_DDRscore_group" %in% colnames(.)) as.character(MO_DDRscore_group) else NA_character_,
        if ("MO_DDRscore_group.anno" %in% colnames(.)) as.character(MO_DDRscore_group.anno) else NA_character_
      ),
      MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High"))
    ) %>%
    select(-matches("\\.anno$"))

  numeric_cols <- setdiff(
    names(tide_res)[sapply(tide_res, function(x) any(is.finite(safe_num(x))))],
    c(
      "Sample", "Patient", "MO_DDRscore_raw", "MO_DDRscore_group",
      "No benefits", "Responder", "CTL.flag"
    )
  )

  tide_stat <- bind_rows(lapply(numeric_cols, function(cc) {
    tmp <- tide_res %>% mutate(Value = safe_num(.data[[cc]]))
    out <- group_compare(tmp, "Value")
    out$Feature <- cc
    out
  })) %>%
    mutate(FDR = p.adjust(P, method = "BH")) %>%
    arrange(FDR, P)

  tide_cor <- bind_rows(lapply(numeric_cols, function(cc) {
    cr <- safe_spearman(safe_num(tide_res[[cc]]), tide_res$MO_DDRscore_raw)
    data.frame(Feature = cc, Rho = cr["Rho"], P = cr["P"])
  })) %>%
    mutate(FDR = p.adjust(P, method = "BH")) %>%
    arrange(FDR, P)

  save_csv(tide_res, file.path(OUT_DIR, "Official_TIDE_merged.csv"))
  save_csv(tide_stat, file.path(OUT_DIR, "Official_TIDE_group_comparison.csv"))
  save_csv(tide_cor, file.path(OUT_DIR, "Official_TIDE_vs_MO_DDRscore_spearman.csv"))

  tide_cat_cols <- intersect(c("Responder", "No benefits", "CTL.flag"), colnames(tide_res))
  if (length(tide_cat_cols) > 0) {
    tide_cat_summary <- bind_rows(lapply(tide_cat_cols, function(cc) {
      tab <- table(tide_res$MO_DDRscore_group, tide_res[[cc]], useNA = "ifany")
      out <- as.data.frame(tab, stringsAsFactors = FALSE)
      colnames(out) <- c("MO_DDRscore_group", "Category", "N")
      out$Feature <- cc
      out
    })) %>%
      select(Feature, MO_DDRscore_group, Category, N)

    tide_cat_test <- bind_rows(lapply(tide_cat_cols, function(cc) {
      tab <- table(tide_res$MO_DDRscore_group, tide_res[[cc]], useNA = "no")
      p <- if (nrow(tab) >= 2 && ncol(tab) >= 2) {
        fisher.test(tab)$p.value
      } else {
        NA_real_
      }
      data.frame(Feature = cc, Fisher_P = p)
    })) %>%
      mutate(FDR = p.adjust(Fisher_P, method = "BH"))

    save_csv(tide_cat_summary, file.path(OUT_DIR, "Official_TIDE_categorical_summary.csv"))
    save_csv(tide_cat_test, file.path(OUT_DIR, "Official_TIDE_categorical_fisher.csv"))
  }
}

############################
# 12. Integrated heatmap
############################

heat_ssgsea <- ssgsea_long %>%
  mutate(ID = paste0("MSigDB:", Module)) %>%
  group_by(ID, Sample) %>%
  summarise(Value = mean(Score, na.rm = TRUE), .groups = "drop")

heat_axis <- axis_long %>%
  transmute(ID = paste0("IOBR_axis:", Axis), Sample, Value = Score)

heat_df <- bind_rows(heat_ssgsea, heat_axis) %>%
  group_by(ID) %>%
  mutate(Value_z = zscore(Value)) %>%
  ungroup()

if (nrow(heat_df) > 0) {
  heat_mat <- heat_df %>%
    select(ID, Sample, Value_z) %>%
    pivot_wider(names_from = Sample, values_from = Value_z) %>%
    column_to_rownames("ID") %>%
    as.matrix()

  sample_order <- anno %>% arrange(MO_DDRscore_group, MO_DDRscore_raw) %>% pull(Sample)
  sample_order <- intersect(sample_order, colnames(heat_mat))
  heat_mat <- heat_mat[, sample_order, drop = FALSE]

  ann_col <- anno %>%
    select(Sample, MO_DDRscore_group, MO_DDRscore_raw) %>%
    distinct(Sample, .keep_all = TRUE) %>%
    filter(Sample %in% sample_order) %>%
    arrange(match(Sample, sample_order)) %>%
    column_to_rownames("Sample")

  pdf(file.path(OUT_DIR, "Fig_official_integrated_followup_heatmap.pdf"),
      width = 10, height = 8)
  pheatmap::pheatmap(
    heat_mat,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_colnames = FALSE,
    annotation_col = ann_col,
    color = colorRampPalette(c("#3B75AF", "white", "#C84630"))(100),
    breaks = seq(-2.5, 2.5, length.out = 101),
    main = "Official immune follow-up modules"
  )
  dev.off()
}

############################
# 13. Auto summary
############################

summary_df <- data.frame(
  Item = c(
    "N_samples",
    "N_MSigDB_sets_used",
    "N_MSigDB_group_FDR_lt_0.05",
    "N_MSigDB_correlation_FDR_lt_0.05",
    "N_official_deconvolution_axis",
    "TMB_integrated",
    "TIDE_result_integrated",
    "Stemness_integrated"
  ),
  Value = c(
    nrow(anno),
    length(msig_list),
    sum(ssgsea_stat$FDR < 0.05, na.rm = TRUE),
    sum(ssgsea_cor$FDR < 0.05, na.rm = TRUE),
    length(unique(axis_long$Axis)),
    ifelse(!is.null(tmb_df), "YES", "NO"),
    ifelse(exists("tide_res") && !is.null(tide_res), "YES", "NO"),
    ifelse(!is.null(stem_df), "YES", "NO")
  )
)

save_csv(summary_df, file.path(OUT_DIR, "Official_followup_summary.csv"))

cat("\nDone.\n")
cat("Output directory:\n", OUT_DIR, "\n")
print(summary_df)
