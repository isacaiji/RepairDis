############################################################
# Standard publication-style immune figure suite for MO-DDRscore
# Uses standard packages:
# ComplexHeatmap for landscape/marker heatmaps
# ggplot2/ggpubr-style panels for dotplots and boxplots
# linkET is handled by the separate Figure 6I-like script
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260513)

############################
# 0. Packages
############################

pkgs <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", "scales", "stringr",
  "ComplexHeatmap", "circlize", "grid", "patchwork"
)

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p %in% c("ComplexHeatmap")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(p, update = FALSE, ask = FALSE)
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
  library(scales)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(patchwork)
})

############################
# 1. Paths
############################

PROJECT_DIR <- "D:/R_workspace/\u8bc4\u5206/AD_DDR_project"
FIG2_DIR <- file.path(PROJECT_DIR, "02_Figure2_MultiOmics_Immune")
BASIC_DIR <- file.path(FIG2_DIR, "Immune_basic_MO_DDRscore")
FOLLOW_DIR <- file.path(FIG2_DIR, "Immune_official_followup_MO_DDRscore")
IOBR_DIR <- file.path(FIG2_DIR, "Immune_IOBR_signature_MO_DDRscore")
NETWORK_DIR <- file.path(FIG2_DIR, "Gene_Immune_Correlation_MO_DDRweight")
OUT_DIR <- file.path(FIG2_DIR, "Immune_publication_figures_MO_DDRscore")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

############################
# 2. Helpers
############################

read_csv <- function(file) {
  data.table::fread(file, data.table = FALSE, check.names = FALSE)
}

save_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file)
}

num <- function(x) suppressWarnings(as.numeric(x))

safe_z <- function(x, cap = 2.5) {
  x <- num(x)
  if (sum(is.finite(x)) < 3 || sd(x, na.rm = TRUE) == 0) {
    z <- rep(0, length(x))
  } else {
    z <- as.numeric(scale(x))
  }
  pmax(pmin(z, cap), -cap)
}

sig_star <- function(p) {
  p <- num(p)
  ifelse(
    is.na(p), "",
    ifelse(p < 0.001, "***",
           ifelse(p < 0.01, "**",
                  ifelse(p < 0.05, "*", "ns")))
  )
}

neglog <- function(p, cap = 35) {
  p <- pmax(num(p), .Machine$double.xmin)
  pmin(-log10(p), cap)
}

clean_feature <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("Homologous recombination", "HR", x)
  x <- gsub("DNA replication", "DNA replication", x)
  x <- gsub("Mismatch Repair", "MMR", x)
  x <- gsub("Base excision repair", "BER", x)
  x <- gsub("Nucleotide excision repair", "NER", x)
  x <- gsub("B Plasma axis", "B/plasma axis", x)
  x <- gsub("Myeloid Neutrophil axis", "Myeloid/neutrophil axis", x)
  x <- gsub("quanTIseq M1 M2 ratio", "M1/M2 ratio", x)
  x <- gsub("TIP Release of cancer cell antigens", "Cancer antigen release", x)
  x <- gsub("TIP Infiltration of immune cells into tumors.*", "Immune-cell infiltration TIP", x)
  x <- gsub("Antigen Processing and Presentation.*", "Antigen processing", x)
  x <- gsub("TCR signaling Pathway.*", "TCR signaling", x)
  x <- gsub("Natural Killer Cell Cytotoxicity.*", "NK cytotoxicity", x)
  x
}

theme_pub <- function(base_size = 9) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      panel.grid.major = element_line(color = "#E7ECF3", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "#CBD4E1", linewidth = 0.35),
      strip.background = element_rect(fill = "#EAF0F7", color = "#CBD4E1", linewidth = 0.35),
      strip.text = element_text(face = "bold", color = "#17213A"),
      axis.text = element_text(color = "#263242"),
      axis.title = element_text(face = "bold", color = "#17213A"),
      plot.title = element_text(face = "bold", color = "#17213A", size = base_size + 2),
      plot.subtitle = element_text(color = "#6A7588"),
      legend.title = element_text(face = "bold", color = "#17213A"),
      legend.text = element_text(color = "#263242"),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

save_plot <- function(plot, name, width, height) {
  pdf_file <- file.path(OUT_DIR, paste0(name, ".pdf"))
  png_file <- file.path(OUT_DIR, paste0(name, ".png"))
  tryCatch(
    ggsave(pdf_file, plot, width = width, height = height, device = cairo_pdf, bg = "white"),
    error = function(e) ggsave(pdf_file, plot, width = width, height = height, device = "pdf", bg = "white")
  )
  ggsave(png_file, plot, width = width, height = height, dpi = 500, bg = "white")
}

save_heatmap <- function(ht, name, width, height) {
  pdf_file <- file.path(OUT_DIR, paste0(name, ".pdf"))
  png_file <- file.path(OUT_DIR, paste0(name, ".png"))
  pdf(pdf_file, width = width, height = height, useDingbats = FALSE)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  png(png_file, width = width, height = height, units = "in", res = 500)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}

low_col <- "#2E6FAD"
high_col <- "#C94B3A"
ink <- "#17213A"
muted <- "#6A7588"
group_cols <- c(Low = low_col, High = high_col)

############################
# 3. Load shared data
############################

axis_long <- read_csv(file.path(FOLLOW_DIR, "Official_deconvolution_result_driven_axes_long.csv")) %>%
  mutate(
    Sample = gsub("\\.", "-", Sample),
    Patient = as.character(Patient),
    MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High")),
    MO_DDRscore_raw = num(MO_DDRscore_raw),
    Axis = as.character(Axis),
    Score = num(Score)
  )

sample_anno <- axis_long %>%
  distinct(Sample, Patient, MO_DDRscore_raw, MO_DDRscore_group) %>%
  filter(!is.na(MO_DDRscore_group), is.finite(MO_DDRscore_raw)) %>%
  arrange(MO_DDRscore_group, MO_DDRscore_raw)

sample_order <- sample_anno$Sample

iobr_long <- read_csv(file.path(IOBR_DIR, "IOBR_official_signature_scores_long.csv")) %>%
  mutate(
    Sample = gsub("\\.", "-", Sample),
    Patient = as.character(Patient),
    MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High")),
    MO_DDRscore_raw = num(MO_DDRscore_raw),
    Signature = as.character(Signature),
    Score = num(Score)
  )

iobr_gc <- read_csv(file.path(IOBR_DIR, "IOBR_official_signature_group_comparison.csv")) %>%
  mutate(
    Signature = as.character(Signature),
    Delta_High_minus_Low = num(Delta_High_minus_Low),
    FDR = num(FDR)
  )

msig_gc <- read_csv(file.path(FOLLOW_DIR, "Official_MSigDB_ssGSEA_group_comparison.csv")) %>%
  mutate(
    Module = as.character(Module),
    MSigDB_Set = as.character(MSigDB_Set),
    Delta_High_minus_Low = num(Delta_High_minus_Low),
    FDR = num(FDR)
  )

deconv_gc <- read_csv(file.path(BASIC_DIR, "Immune_deconvolution_group_comparison.csv")) %>%
  mutate(
    Method = as.character(Method),
    Feature = as.character(Feature),
    Delta_High_minus_Low = num(Delta_High_minus_Low),
    FDR = num(FDR)
  )

tide <- read_csv(file.path(FOLLOW_DIR, "Official_TIDE_merged.csv")) %>%
  mutate(
    Sample = gsub("\\.", "-", Sample),
    Patient = as.character(Patient),
    MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High")),
    MO_DDRscore_raw = num(MO_DDRscore_raw)
  )

tide_gc <- read_csv(file.path(FOLLOW_DIR, "Official_TIDE_group_comparison.csv")) %>%
  mutate(
    Feature = as.character(Feature),
    Delta_High_minus_Low = num(Delta_High_minus_Low),
    FDR = num(FDR)
  )

tide_cat <- read_csv(file.path(FOLLOW_DIR, "Official_TIDE_categorical_summary.csv")) %>%
  mutate(
    MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High")),
    N = num(N)
  )

tide_fisher <- read_csv(file.path(FOLLOW_DIR, "Official_TIDE_categorical_fisher.csv")) %>%
  mutate(FDR = num(FDR), Fisher_P = num(Fisher_P))

tmb <- read_csv(file.path(FOLLOW_DIR, "Official_TMB_merged.csv")) %>%
  mutate(
    Sample = gsub("\\.", "-", Sample),
    Patient = as.character(Patient),
    MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High")),
    MO_DDRscore_raw = num(MO_DDRscore_raw),
    TMB_value = num(TMB_value),
    TMB_log1p = log1p(TMB_value)
  )

tmb_gc <- read_csv(file.path(FOLLOW_DIR, "Official_TMB_group_comparison.csv")) %>%
  mutate(Feature = "TMB", Delta_High_minus_Low = num(Delta_High_minus_Low), P = num(P), FDR = num(P))

math <- read_csv(file.path(FOLLOW_DIR, "Official_MATH_score.csv")) %>%
  mutate(
    Patient = as.character(Patient),
    MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High")),
    MO_DDRscore_raw = num(MO_DDRscore_raw),
    MATH = num(MATH)
  )

math_gc <- read_csv(file.path(FOLLOW_DIR, "Official_MATH_group_comparison.csv")) %>%
  mutate(Feature = "MATH", Delta_High_minus_Low = num(Delta_High_minus_Low), P = num(P), FDR = num(P))

marker_long <- read_csv(file.path(BASIC_DIR, "Immune_marker_expression_long.csv")) %>%
  mutate(
    Sample = gsub("\\.", "-", Sample),
    Patient = as.character(Patient),
    MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High")),
    MO_DDRscore_raw = num(MO_DDRscore_raw),
    Panel = as.character(Panel),
    Gene = as.character(Gene),
    Expression = num(Expression)
  )

marker_gc <- read_csv(file.path(BASIC_DIR, "Immune_marker_group_comparison.csv")) %>%
  mutate(
    Panel = as.character(Panel),
    Gene = as.character(Gene),
    Delta_High_minus_Low = num(Delta_High_minus_Low),
    FDR = num(FDR)
  )

############################
# 4. Figure A: ComplexHeatmap immune landscape
############################

axis_keep <- c("Global_TME", "B_Plasma_axis", "Myeloid_Neutrophil_axis", "quanTIseq_M1_M2_ratio", "IPS_axis")
iobr_keep <- c(
  "DDR", "Homologous_recombination", "DNA_replication", "Mismatch_Repair",
  "Base_excision_repair", "Nucleotide_excision_repair",
  "MDSC_Peng_et_al", "TMEscoreA_plus",
  "TIP_Release_of_cancer_cell_antigens",
  "Antigen_Processing_and_Presentation_Li_et_al"
)
tide_keep <- c("TIDE", "Exclusion", "MDSC", "CAF", "CD274")

feature_landscape <- bind_rows(
  iobr_long %>%
    filter(Signature %in% iobr_keep) %>%
    transmute(Sample, Patient, MO_DDRscore_group, MO_DDRscore_raw,
              Category = "DDR / immune programs",
              FeatureID = Signature,
              Feature = clean_feature(Signature),
              Value = Score),
  axis_long %>%
    filter(Axis %in% axis_keep) %>%
    transmute(Sample, Patient, MO_DDRscore_group, MO_DDRscore_raw,
              Category = "Immune axes",
              FeatureID = Axis,
              Feature = clean_feature(Axis),
              Value = ifelse(Axis == "quanTIseq_M1_M2_ratio", log10(pmax(Score, 0) + 1), Score)),
  tide %>%
    select(Sample, Patient, MO_DDRscore_group, MO_DDRscore_raw, all_of(tide_keep)) %>%
    pivot_longer(cols = all_of(tide_keep), names_to = "FeatureID", values_to = "Value") %>%
    transmute(Sample, Patient, MO_DDRscore_group, MO_DDRscore_raw,
              Category = "TIDE / ICB",
              FeatureID,
              Feature = FeatureID,
              Value = num(Value)),
  tmb %>%
    transmute(Sample = sample_anno$Sample[match(Patient, sample_anno$Patient)],
              Patient, MO_DDRscore_group, MO_DDRscore_raw,
              Category = "Genomic instability",
              FeatureID = "TMB",
              Feature = "TMB",
              Value = TMB_log1p),
  math %>%
    transmute(Sample = sample_anno$Sample[match(Patient, sample_anno$Patient)],
              Patient, MO_DDRscore_group, MO_DDRscore_raw,
              Category = "Genomic instability",
              FeatureID = "MATH",
              Feature = "MATH",
              Value = MATH)
) %>%
  filter(Sample %in% sample_order, is.finite(Value)) %>%
  mutate(
    FeatureKey = paste(Category, FeatureID, sep = "__")
  ) %>%
  distinct(FeatureKey, Sample, .keep_all = TRUE)

land_info <- feature_landscape %>%
  distinct(FeatureKey, Category, Feature) %>%
  mutate(
    Category = factor(Category, levels = c("DDR / immune programs", "Genomic instability", "Immune axes", "TIDE / ICB")),
    Feature = factor(Feature, levels = Feature)
  ) %>%
  arrange(Category)

land_wide <- feature_landscape %>%
  select(FeatureKey, Sample, Value) %>%
  pivot_wider(names_from = Sample, values_from = Value)

land_mat <- as.data.frame(land_wide)
rownames(land_mat) <- land_mat$FeatureKey
land_mat$FeatureKey <- NULL
land_mat <- as.matrix(land_mat[land_info$FeatureKey, sample_order, drop = FALSE])
land_z <- t(apply(land_mat, 1, safe_z))
rownames(land_z) <- as.character(land_info$Feature)

score_col_fun <- colorRamp2(
  quantile(sample_anno$MO_DDRscore_raw, probs = c(0.02, 0.5, 0.98), na.rm = TRUE),
  c("#2E6FAD", "white", "#C94B3A")
)

ha_land <- HeatmapAnnotation(
  `MO-DDRscore group` = sample_anno$MO_DDRscore_group,
  `MO-DDRscore` = sample_anno$MO_DDRscore_raw,
  col = list(
    `MO-DDRscore group` = group_cols,
    `MO-DDRscore` = score_col_fun
  ),
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
  simple_anno_size = unit(3.5, "mm")
)

ht_land <- Heatmap(
  land_z,
  name = "Z-score",
  col = colorRamp2(c(-2.2, 0, 2.2), c("#2E6FAD", "white", "#C94B3A")),
  top_annotation = ha_land,
  row_split = land_info$Category,
  column_split = sample_anno$MO_DDRscore_group,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8),
  row_title_gp = gpar(fontsize = 7, fontface = "bold"),
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  heatmap_legend_param = list(title_gp = gpar(fontface = "bold"))
)

save_heatmap(ht_land, "Std_A_ComplexHeatmap_immune_landscape", 8.6, 6.2)
save_csv(feature_landscape, file.path(OUT_DIR, "Std_A_immune_landscape_long.csv"))

############################
# 5. Figure B: standard deconvolution dotplot
############################

deconv_map <- tibble::tribble(
  ~Method, ~Feature, ~Display, ~Class,
  "estimate", "ImmuneScore_estimate", "ImmuneScore", "Global TME",
  "estimate", "StromalScore_estimate", "StromalScore", "Global TME",
  "estimate", "TumorPurity_estimate", "Tumor purity", "Global TME",
  "epic", "Bcells_EPIC", "B cells", "B lineage",
  "timer", "B_cell_TIMER", "B cells", "B lineage",
  "mcpcounter", "B_lineage_MCPcounter", "B lineage", "B lineage",
  "quantiseq", "B_cells_quantiseq", "B cells", "B lineage",
  "xcell", "Plasma_cells_xCell", "Plasma cells", "B lineage",
  "epic", "CD8_Tcells_EPIC", "CD8 T cells", "T / NK lineage",
  "quantiseq", "NK_cells_quantiseq", "NK cells", "T / NK lineage",
  "mcpcounter", "NK_cells_MCPcounter", "NK cells", "T / NK lineage",
  "quantiseq", "Macrophages_M1_quantiseq", "M1 macrophages", "Myeloid",
  "quantiseq", "Macrophages_M2_quantiseq", "M2 macrophages", "Myeloid",
  "timer", "Neutrophil_TIMER", "Neutrophils", "Myeloid",
  "quantiseq", "Neutrophils_quantiseq", "Neutrophils", "Myeloid",
  "mcpcounter", "Monocytic_lineage_MCPcounter", "Monocytic lineage", "Myeloid",
  "mcpcounter", "Fibroblasts_MCPcounter", "Fibroblasts", "Stroma",
  "epic", "CAFs_EPIC", "CAFs", "Stroma",
  "ips", "IPS_IPS", "IPS", "Immunophenotype",
  "ips", "AZ_IPS", "AZ", "Immunophenotype",
  "ips", "MHC_IPS", "MHC", "Immunophenotype"
)

deconv_dot <- deconv_gc %>%
  inner_join(deconv_map, by = c("Method", "Feature")) %>%
  mutate(
    Class = factor(Class, levels = c("Global TME", "B lineage", "T / NK lineage", "Myeloid", "Stroma", "Immunophenotype")),
    Method = factor(Method, levels = c("estimate", "epic", "timer", "mcpcounter", "quantiseq", "xcell", "ips")),
    Direction = ifelse(Delta_High_minus_Low >= 0, "High-score enriched", "Low-score enriched"),
    Stars = sig_star(FDR),
    NegLogFDR = neglog(FDR, cap = 25),
    SignedScore = sign(Delta_High_minus_Low) * NegLogFDR,
    DisplayMethod = paste0(Display, "  [", Method, "]")
  ) %>%
  arrange(Class, SignedScore) %>%
  mutate(
    DisplayMethod = factor(DisplayMethod, levels = rev(unique(DisplayMethod))),
    StarX = SignedScore + ifelse(SignedScore >= 0, 0.75, -0.75)
  )

pB <- ggplot(deconv_dot, aes(SignedScore, DisplayMethod)) +
  geom_vline(xintercept = 0, color = "#9BA7B7", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = SignedScore, yend = DisplayMethod, color = Direction),
               linewidth = 0.55, alpha = 0.70) +
  geom_point(aes(size = NegLogFDR, fill = Direction),
             shape = 21, color = "white", stroke = 0.35, alpha = 0.98) +
  geom_text(aes(x = StarX, label = ifelse(Stars == "ns", "", Stars)),
            size = 2.7, color = ink, fontface = "bold") +
  facet_grid(Class ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("High-score enriched" = high_col, "Low-score enriched" = low_col), name = NULL) +
  scale_color_manual(values = c("High-score enriched" = high_col, "Low-score enriched" = low_col), name = NULL) +
  scale_size_continuous(range = c(2.0, 7.2), name = "-log10(FDR)") +
  labs(
    title = "Differential immune infiltration across deconvolution algorithms",
    subtitle = "Direction is High-minus-Low MO-DDRscore; distance from zero represents -log10(FDR)",
    x = "Signed significance, -log10(FDR)",
    y = NULL
  ) +
  theme_pub(8) +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0),
    panel.grid.major = element_line(color = "#EDF1F6", linewidth = 0.25),
    panel.grid.major.y = element_blank(),
    legend.position = "right"
  )

save_plot(pB, "Std_B_Deconvolution_delta_dotplot", 7.8, 8.4)
save_csv(deconv_dot, file.path(OUT_DIR, "Std_B_deconvolution_delta_dotplot_data.csv"))

############################
# 6. Figure C: immune function and pathway dotplot
############################

iobr_select_patterns <- paste(
  c(
    "DDR$", "Homologous_recombination", "DNA_replication", "Mismatch_Repair",
    "Base_excision_repair", "Nucleotide_excision_repair",
    "MDSC_Peng", "TMEscoreA", "TMEscore_CIR", "TMEscore_plus",
    "TIP_Release", "TIP_Infiltration", "Antigen_Processing",
    "TCR_signaling", "Natural_Killer_Cell_Cytotoxicity",
    "CD8_T_cells_Bindea", "B_cells_Danaher", "B_cells_Bindea",
    "TLS_Nature", "MHC_Class_II", "WNT_target"
  ),
  collapse = "|"
)

iobr_func <- iobr_gc %>%
  filter(grepl(iobr_select_patterns, Signature)) %>%
  mutate(
    Source = "IOBR signature",
    Category = case_when(
      grepl("DDR|Homologous|DNA_replication|Mismatch|excision", Signature) ~ "DDR / repair",
      grepl("B_cells|TLS|MHC_Class_II", Signature) ~ "B/TLS/MHC-II",
      grepl("MDSC|TMEscore|WNT", Signature) ~ "Suppression / exclusion",
      TRUE ~ "Antigen / effector programs"
    ),
    Label = clean_feature(Signature)
  ) %>%
  select(Source, Category, Label, Delta_High_minus_Low, FDR)

msig_func <- msig_gc %>%
  group_by(Module) %>%
  summarise(
    Delta_High_minus_Low = median(Delta_High_minus_Low, na.rm = TRUE),
    FDR = min(FDR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(Module %in% c(
    "REACTOME_ANTIGEN_PRESENTATION", "REACTOME_BCR_SIGNALING",
    "GO_BCR_SIGNALING", "GO_B_CELL_ACTIVATION", "GO_PLASMA_CELL_DIFFERENTIATION",
    "REACTOME_PD1_SIGNALING", "HALLMARK_TGF_BETA", "HALLMARK_INFLAMMATION",
    "GO_MYELOID_LEUKOCYTE_ACTIVATION", "GO_MONONUCLEAR_CELL_MIGRATION",
    "HALLMARK_COMPLEMENT", "HALLMARK_IL6_JAK_STAT3"
  )) %>%
  mutate(
    Source = "MSigDB ssGSEA",
    Category = case_when(
      grepl("BCR|B_CELL|PLASMA", Module) ~ "B/TLS/MHC-II",
      grepl("PD1|ANTIGEN", Module) ~ "Antigen / effector programs",
      grepl("MYELOID|MIGRATION|INFLAMMATION|COMPLEMENT|IL6|TGF", Module) ~ "Suppression / exclusion",
      TRUE ~ "Immune programs"
    ),
    Label = clean_feature(Module)
  ) %>%
  select(Source, Category, Label, Delta_High_minus_Low, FDR)

func_dot <- bind_rows(iobr_func, msig_func) %>%
  filter(is.finite(FDR), is.finite(Delta_High_minus_Low)) %>%
  group_by(Category) %>%
  arrange(FDR, .by_group = TRUE) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  mutate(
    Category = factor(Category, levels = c("DDR / repair", "Antigen / effector programs", "B/TLS/MHC-II", "Suppression / exclusion")),
    Label = factor(Label, levels = rev(unique(Label))),
    Stars = sig_star(FDR),
    NegLogFDR = neglog(FDR, cap = 45),
    LabelX = Delta_High_minus_Low + ifelse(Delta_High_minus_Low >= 0, 0.04, -0.04)
  )

pC <- ggplot(func_dot, aes(Delta_High_minus_Low, Label)) +
  geom_vline(xintercept = 0, color = "#9BA7B7", linewidth = 0.35) +
  geom_point(aes(size = NegLogFDR, fill = Delta_High_minus_Low),
             shape = 21, color = "white", stroke = 0.45, alpha = 0.98) +
  geom_text(aes(x = LabelX, label = ifelse(Stars == "ns", "", Stars)),
            size = 2.5, color = ink, fontface = "bold") +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = low_col, mid = "white", high = high_col, midpoint = 0,
    name = "Median delta\nHigh - Low"
  ) +
  scale_size_continuous(range = c(2.2, 7.8), name = "-log10(FDR)") +
  labs(
    title = "Official immune and DDR program differences",
    subtitle = "IOBR published signatures and MSigDB ssGSEA modules",
    x = "Median difference, High - Low MO-DDRscore",
    y = NULL
  ) +
  theme_pub(8) +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0),
    panel.grid.major.y = element_blank(),
    legend.position = "right"
  )

save_plot(pC, "Std_C_Official_functional_program_dotplot", 8.5, 8.2)
save_csv(func_dot, file.path(OUT_DIR, "Std_C_functional_program_dotplot_data.csv"))

############################
# 7. Figure D: TMB/MATH/TIDE violin-box panels
############################

tide_box_features <- c("TIDE", "Exclusion", "MDSC", "CAF", "CD274")

box_long <- bind_rows(
  tmb %>%
    filter(!is.na(MO_DDRscore_group)) %>%
    transmute(MO_DDRscore_group, Feature = "TMB\n(log1p)", Value = TMB_log1p),
  math %>%
    filter(!is.na(MO_DDRscore_group)) %>%
    transmute(MO_DDRscore_group, Feature = "MATH", Value = MATH),
  tide %>%
    filter(!is.na(MO_DDRscore_group)) %>%
    select(MO_DDRscore_group, all_of(tide_box_features)) %>%
    pivot_longer(cols = all_of(tide_box_features), names_to = "Feature", values_to = "Value") %>%
    mutate(Value = num(Value))
) %>%
  filter(is.finite(Value)) %>%
  mutate(
    Feature = factor(Feature, levels = c("TMB\n(log1p)", "MATH", "TIDE", "Exclusion", "MDSC", "CAF", "CD274"))
  )

box_sig <- bind_rows(
  tmb_gc %>% transmute(Feature = "TMB\n(log1p)", FDR = FDR),
  math_gc %>% transmute(Feature = "MATH", FDR = FDR),
  tide_gc %>% filter(Feature %in% tide_box_features) %>% transmute(Feature, FDR = FDR)
) %>%
  mutate(
    Feature = factor(Feature, levels = levels(box_long$Feature)),
    Stars = sig_star(FDR)
  )

pD1 <- ggplot(box_long, aes(MO_DDRscore_group, Value, fill = MO_DDRscore_group)) +
  geom_violin(width = 0.88, alpha = 0.18, color = NA, trim = TRUE) +
  geom_boxplot(width = 0.23, outlier.shape = NA, linewidth = 0.35, color = "#273246", alpha = 0.82) +
  geom_jitter(aes(color = MO_DDRscore_group), width = 0.12, size = 0.45, alpha = 0.25, show.legend = FALSE) +
  geom_text(data = box_sig, aes(x = 1.5, y = Inf, label = Stars),
            inherit.aes = FALSE, vjust = 1.2, size = 4, fontface = "bold", color = ink) +
  facet_wrap(~Feature, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = group_cols, drop = FALSE) +
  scale_color_manual(values = group_cols, drop = FALSE) +
  labs(
    title = "Genomic instability and TIDE-derived immunotherapy-related features",
    subtitle = "Violin/box plots compare Low and High MO-DDRscore groups",
    x = NULL,
    y = "Value"
  ) +
  theme_pub(8) +
  theme(legend.position = "none", panel.grid.major.x = element_blank())

responder_plot <- tide_cat %>%
  filter(Feature == "Responder") %>%
  group_by(MO_DDRscore_group) %>%
  mutate(
    Proportion = N / sum(N),
    Category = ifelse(as.character(Category) == "TRUE", "Predicted responder", "Non-responder")
  ) %>%
  ungroup()

responder_p <- tide_fisher %>%
  filter(Feature == "Responder") %>%
  pull(FDR)

pD2 <- ggplot(responder_plot, aes(MO_DDRscore_group, Proportion, fill = Category)) +
  geom_col(width = 0.58, color = "white", linewidth = 0.35) +
  geom_text(aes(label = percent(Proportion, accuracy = 0.1)),
            position = position_stack(vjust = 0.5), color = "white",
            fontface = "bold", size = 3.1) +
  annotate("text", x = 1.5, y = 1.05, label = sig_star(responder_p),
           size = 5, fontface = "bold", color = ink) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1.1), expand = c(0, 0)) +
  scale_fill_manual(values = c("Predicted responder" = high_col, "Non-responder" = "#AAB4C3")) +
  labs(
    title = "TIDE-predicted ICB response",
    x = NULL,
    y = "Proportion",
    fill = NULL
  ) +
  theme_pub(8) +
  theme(panel.grid.major.x = element_blank(), legend.position = "bottom")

pD <- pD1 / pD2 + plot_layout(heights = c(1, 0.88))

save_plot(pD, "Std_D_TMB_MATH_TIDE_violin_and_response", 10.2, 6.2)
save_csv(box_long, file.path(OUT_DIR, "Std_D_TMB_MATH_TIDE_boxplot_long.csv"))
save_csv(responder_plot, file.path(OUT_DIR, "Std_D_TIDE_responder_fraction.csv"))

############################
# 8. Figure E: ComplexHeatmap marker expression
############################

marker_keep <- c(
  "CD274", "PDCD1LG2", "LAG3", "TNFRSF9",
  "TAP1", "TAP2", "NLRC5", "PSMB8", "PSMB9",
  "CXCL9", "CXCL10", "CXCL11"
)

marker_info <- marker_long %>%
  filter(Gene %in% marker_keep) %>%
  distinct(Gene, Panel) %>%
  mutate(
    Panel = factor(Panel, levels = c("Checkpoint", "HLA_APM", "Chemokine_Receptor")),
    Gene = factor(Gene, levels = marker_keep)
  ) %>%
  arrange(Panel, Gene)

marker_wide <- marker_long %>%
  filter(Gene %in% marker_keep, Sample %in% sample_order) %>%
  select(Gene, Sample, Expression) %>%
  distinct(Gene, Sample, .keep_all = TRUE) %>%
  pivot_wider(names_from = Sample, values_from = Expression)

marker_mat <- as.data.frame(marker_wide)
rownames(marker_mat) <- marker_mat$Gene
marker_mat$Gene <- NULL
marker_mat <- as.matrix(marker_mat[as.character(marker_info$Gene), sample_order, drop = FALSE])
marker_z <- t(apply(marker_mat, 1, safe_z))
rownames(marker_z) <- as.character(marker_info$Gene)

ha_marker <- HeatmapAnnotation(
  `MO-DDRscore group` = sample_anno$MO_DDRscore_group,
  `MO-DDRscore` = sample_anno$MO_DDRscore_raw,
  col = list(
    `MO-DDRscore group` = group_cols,
    `MO-DDRscore` = score_col_fun
  ),
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
  simple_anno_size = unit(3.5, "mm")
)

ht_marker <- Heatmap(
  marker_z,
  name = "Z-score",
  col = colorRamp2(c(-2.2, 0, 2.2), c("#2E6FAD", "white", "#C94B3A")),
  top_annotation = ha_marker,
  row_split = marker_info$Panel,
  column_split = sample_anno$MO_DDRscore_group,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 8, fontface = "bold"),
  row_title_gp = gpar(fontsize = 7, fontface = "bold"),
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  heatmap_legend_param = list(title_gp = gpar(fontface = "bold"))
)

save_heatmap(ht_marker, "Std_E_ComplexHeatmap_checkpoint_APM_CXCL_markers", 7.5, 4.8)
save_csv(marker_long %>% filter(Gene %in% marker_keep), file.path(OUT_DIR, "Std_E_marker_expression_heatmap_long.csv"))

############################
# 9. Figure F: linkET regulatory map
############################

linket_script <- "C:/Users/nova/Documents/New project/MO_DDRweight_gene_immune_cell_linkET_plot.R"
if (file.exists(linket_script)) {
  source(linket_script, local = TRUE)
}

############################
# 10. Summary
############################

summary_df <- data.frame(
  Panel = c(
    "A", "B", "C", "D", "E", "F"
  ),
  Figure = c(
    "Std_A_ComplexHeatmap_immune_landscape.pdf",
    "Std_B_Deconvolution_delta_dotplot.pdf",
    "Std_C_Official_functional_program_dotplot.pdf",
    "Std_D_TMB_MATH_TIDE_violin_and_response.pdf",
    "Std_E_ComplexHeatmap_checkpoint_APM_CXCL_markers.pdf",
    "Figure6I_linkET_MO_DDRweight_gene_immune_cell_network.pdf"
  ),
  Package = c(
    "ComplexHeatmap",
    "ggplot2",
    "ggplot2",
    "ggplot2 / patchwork",
    "ComplexHeatmap",
    "linkET"
  ),
  Meaning = c(
    "Integrated MO-DDRscore immune landscape",
    "Differential immune infiltration across deconvolution methods",
    "Official IOBR/MSigDB immune and DDR program differences",
    "TMB, MATH, TIDE and predicted response distributions",
    "Checkpoint, antigen-presentation, and chemokine marker expression",
    "Gene-immune cell regulatory association map"
  )
)

save_csv(summary_df, file.path(OUT_DIR, "Std_immune_publication_figure_suite_summary.csv"))

cat("\nDone.\n")
cat("Output directory:\n", OUT_DIR, "\n")
print(summary_df)
