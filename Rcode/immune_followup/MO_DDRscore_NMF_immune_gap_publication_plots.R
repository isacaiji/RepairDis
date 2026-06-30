############################################################
# Publication plots for MO-DDRscore / NMF immune gap analyses
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(20260513)

pkgs <- c(
  "data.table", "dplyr", "tidyr", "ggplot2", "scales",
  "stringr", "patchwork", "survival"
)

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
  library(scales)
  library(stringr)
  library(patchwork)
  library(survival)
})

############################
# 1. Paths
############################

PROJECT_DIR <- "D:/R_workspace/评分/AD_DDR_project"
FIG2_DIR <- file.path(PROJECT_DIR, "02_Figure2_MultiOmics_Immune")
GAP_DIR <- Sys.getenv(
  "IMMUNE_GAP_DIR",
  unset = file.path(FIG2_DIR, "Immune_reference_gap_fill_allinone")
)
OUT_DIR <- Sys.getenv(
  "IMMUNE_GAP_PLOT_OUT_DIR",
  unset = file.path(GAP_DIR, "publication_plots")
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

############################
# 2. Helpers
############################

read_csv <- function(file) data.table::fread(file, data.table = FALSE, check.names = FALSE)

save_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file)
}

num <- function(x) suppressWarnings(as.numeric(x))

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

clean_label <- function(x) {
  x <- gsub("_", " ", x)
  x <- gsub("Homologous recombination", "HR", x)
  x <- gsub("Mismatch Repair", "MMR", x)
  x <- gsub("Base excision repair", "BER", x)
  x <- gsub("Nucleotide excision repair", "NER", x)
  x <- gsub("DNA replication", "DNA replication", x)
  x <- gsub("MDSC Peng et al", "MDSC", x)
  x <- gsub("TIP Release of cancer cell antigens", "Cancer antigen release", x)
  x <- gsub("TIP Infiltration of immune cells into tumors.*", "Immune-cell infiltration TIP", x)
  x <- gsub("Antigen Processing and Presentation.*", "Antigen processing", x)
  x <- gsub("Natural Killer Cell Cytotoxicity.*", "NK cytotoxicity", x)
  x <- gsub("TCR signaling Pathway.*", "TCR signaling", x)
  x <- gsub("quantiseq", "quanTIseq", x, ignore.case = TRUE)
  x <- gsub("mcpcounter", "MCPcounter", x, ignore.case = TRUE)
  x <- gsub("estimate", "ESTIMATE", x, ignore.case = TRUE)
  x <- gsub("timer", "TIMER", x, ignore.case = TRUE)
  x <- gsub("epic", "EPIC", x, ignore.case = TRUE)
  x <- gsub("xcell", "xCell", x, ignore.case = TRUE)
  x <- gsub("ips", "IPS", x, ignore.case = TRUE)
  x
}

theme_pub <- function(base_size = 8.5) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      panel.grid.major = element_line(color = "#E8EDF4", linewidth = 0.25),
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

format_p <- function(p) {
  p <- num(p)
  ifelse(is.na(p), "P = NA", ifelse(p < 0.001, paste0("P = ", format(p, scientific = TRUE, digits = 2)),
                                    paste0("P = ", signif(p, 3))))
}

low_col <- "#2E6FAD"
high_col <- "#C94B3A"
c1_col <- "#2E6FAD"
c2_col <- "#C94B3A"
ink <- "#17213A"
group_cols <- c(Low = low_col, High = high_col, C1 = c1_col, C2 = c2_col)

############################
# 3. Load result tables
############################

anno <- read_csv(file.path(GAP_DIR, "Master_MO_DDRscore_NMF_annotation.csv")) %>%
  mutate(
    MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High")),
    NMF_subtype = factor(NMF_subtype, levels = c("C1", "C2")),
    MO_DDRscore_raw = num(MO_DDRscore_raw)
  )

nmf_test <- read_csv(file.path(GAP_DIR, "NMF_subtype_vs_MO_DDRscore_group_test.csv"))
deconv_cmp <- read_csv(file.path(GAP_DIR, "NMF_immune_deconvolution_group_comparison.csv"))
iobr_cmp <- read_csv(file.path(GAP_DIR, "NMF_IOBR_signature_group_comparison.csv"))
tide_cmp <- read_csv(file.path(GAP_DIR, "NMF_TIDE_features_group_comparison.csv"))
tide_long <- read_csv(file.path(GAP_DIR, "NMF_TIDE_features_long.csv"))
tide_resp <- read_csv(file.path(GAP_DIR, "NMF_TIDE_Responder_counts.csv"))
tide_cat_test <- read_csv(file.path(GAP_DIR, "NMF_TIDE_categorical_tests.csv"))
tmb <- read_csv(file.path(GAP_DIR, "NMF_TMB_merged.csv"))
math <- read_csv(file.path(GAP_DIR, "NMF_MATH_merged.csv"))
tmb_cmp <- read_csv(file.path(GAP_DIR, "NMF_TMB_group_comparison.csv"))
math_cmp <- read_csv(file.path(GAP_DIR, "NMF_MATH_group_comparison.csv"))
cp_mo <- read_csv(file.path(GAP_DIR, "Official_checkpoint_MO_DDRscore_group_comparison.csv"))
cp_nmf <- read_csv(file.path(GAP_DIR, "Official_checkpoint_NMF_subtype_group_comparison.csv"))
surv_df <- read_csv(file.path(GAP_DIR, "Survival_MO_NMF_TMB_MATH_merged.csv"))
surv_tests <- read_csv(file.path(GAP_DIR, "Survival_MO_NMF_TMB_MATH_tests.csv"))

############################
# 4. Figure A: NMF concordance with MO-DDRscore
############################

prop_df <- anno %>%
  filter(!is.na(NMF_subtype), !is.na(MO_DDRscore_group)) %>%
  count(NMF_subtype, MO_DDRscore_group) %>%
  group_by(NMF_subtype) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup()

fisher_p <- nmf_test$Value[nmf_test$Item == "Fisher_P"][1]
cramer_v <- nmf_test$Value[nmf_test$Item == "Cramers_V"][1]

pA1 <- ggplot(prop_df, aes(NMF_subtype, Proportion, fill = MO_DDRscore_group)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.35) +
  geom_text(aes(label = paste0(n, "\n", percent(Proportion, accuracy = 0.1))),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 3.1) +
  scale_fill_manual(values = c(Low = low_col, High = high_col), drop = FALSE) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "NMF subtypes recapitulate MO-DDRscore states",
    subtitle = paste0(format_p(fisher_p), "; Cramer's V = ", signif(num(cramer_v), 3)),
    x = "NMF-derived DDR subtype",
    y = "Proportion",
    fill = "MO-DDRscore"
  ) +
  theme_pub(9) +
  theme(panel.grid.major.x = element_blank(), legend.position = "right")

pA2 <- ggplot(anno, aes(NMF_subtype, MO_DDRscore_raw, fill = NMF_subtype)) +
  geom_violin(width = 0.78, color = NA, alpha = 0.18, trim = TRUE) +
  geom_boxplot(width = 0.22, outlier.shape = NA, linewidth = 0.35, alpha = 0.9) +
  geom_jitter(aes(color = NMF_subtype), width = 0.12, size = 0.55, alpha = 0.28, show.legend = FALSE) +
  scale_fill_manual(values = c(C1 = c1_col, C2 = c2_col), drop = FALSE) +
  scale_color_manual(values = c(C1 = c1_col, C2 = c2_col), drop = FALSE) +
  labs(
    title = "MO-DDRscore gradient across NMF subtypes",
    x = "NMF-derived DDR subtype",
    y = "MO-DDRscore"
  ) +
  theme_pub(9) +
  theme(legend.position = "none", panel.grid.major.x = element_blank())

pA <- pA1 + pA2 + plot_layout(widths = c(1.05, 1))
save_plot(pA, "Gap_A_NMF_MO_DDRscore_concordance", 9.4, 4.2)

############################
# 5. Figure B: NMF immune deconvolution lollipop
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
  "quantiseq", "T_cells_CD8_quantiseq", "CD8 T cells", "T / NK lineage",
  "mcpcounter", "CD8_T_cells_MCPcounter", "CD8 T cells", "T / NK lineage",
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

deconv_dot <- deconv_cmp %>%
  inner_join(deconv_map, by = c("Method", "Feature")) %>%
  mutate(
    Class = factor(Class, levels = c("Global TME", "B lineage", "T / NK lineage", "Myeloid", "Stroma", "Immunophenotype")),
    Method = clean_label(Method),
    Direction = ifelse(num(Delta_Group2_minus_Group1) >= 0, "C2 enriched", "C1 enriched"),
    NegLogFDR = neglog(FDR, cap = 25),
    SignedScore = sign(num(Delta_Group2_minus_Group1)) * NegLogFDR,
    Label = paste0(Display, "  [", Method, "]"),
    Stars = ifelse(Significance == "ns", "", Significance)
  ) %>%
  arrange(Class, SignedScore) %>%
  mutate(
    Label = factor(Label, levels = rev(unique(Label))),
    StarX = SignedScore + ifelse(SignedScore >= 0, 0.75, -0.75)
  )

pB <- ggplot(deconv_dot, aes(SignedScore, Label)) +
  geom_vline(xintercept = 0, color = "#9BA7B7", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = SignedScore, yend = Label, color = Direction),
               linewidth = 0.55, alpha = 0.70) +
  geom_point(aes(size = NegLogFDR, fill = Direction),
             shape = 21, color = "white", stroke = 0.35, alpha = 0.98) +
  geom_text(aes(x = StarX, label = Stars), size = 2.6, color = ink, fontface = "bold") +
  facet_grid(Class ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("C2 enriched" = c2_col, "C1 enriched" = c1_col), name = NULL) +
  scale_color_manual(values = c("C2 enriched" = c2_col, "C1 enriched" = c1_col), name = NULL) +
  scale_size_continuous(range = c(2.0, 7.2), name = "-log10(FDR)") +
  labs(
    title = "Immune infiltration differences between NMF-derived DDR subtypes",
    subtitle = "Direction is C2-minus-C1; distance from zero represents -log10(FDR)",
    x = "Signed significance, -log10(FDR)",
    y = NULL
  ) +
  theme_pub(8) +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0),
    panel.grid.major.y = element_blank(),
    legend.position = "right"
  )

save_plot(pB, "Gap_B_NMF_immune_deconvolution_lollipop", 7.8, 8.3)
save_csv(deconv_dot, file.path(OUT_DIR, "Gap_B_plot_data.csv"))

############################
# 6. Figure C: NMF IOBR immune programs
############################

iobr_keep_pattern <- paste(
  c(
    "DDR$", "Homologous_recombination", "DNA_replication", "Mismatch_Repair",
    "Base_excision_repair", "Nucleotide_excision_repair",
    "MDSC_Peng", "TMEscoreA", "TMEscore_CIR", "TMEscore_plus",
    "TIP_Release", "TIP_Infiltration", "Antigen_Processing",
    "TCR_signaling", "Natural_Killer_Cell_Cytotoxicity",
    "CD8_T_cells_Bindea", "B_cells_Danaher", "B_cells_Bindea",
    "TLS_Nature", "MHC_Class_II", "Type_II_IFN", "Mast_cells"
  ),
  collapse = "|"
)

iobr_dot <- iobr_cmp %>%
  filter(grepl(iobr_keep_pattern, Signature)) %>%
  mutate(
    Category = case_when(
      grepl("DDR|Homologous|DNA_replication|Mismatch|excision", Signature) ~ "DDR / repair",
      grepl("B_cells|TLS|MHC_Class_II|Type_II_IFN", Signature) ~ "B/TLS/MHC-II",
      grepl("MDSC|TMEscore|WNT|Mast_cells", Signature) ~ "Suppression / stroma",
      TRUE ~ "Antigen / effector"
    ),
    Category = factor(Category, levels = c("DDR / repair", "Antigen / effector", "B/TLS/MHC-II", "Suppression / stroma")),
    Label = clean_label(Signature),
    Delta = num(Delta_Group2_minus_Group1),
    NegLogFDR = neglog(FDR, cap = 45),
    Stars = ifelse(Significance == "ns", "", Significance)
  ) %>%
  group_by(Category) %>%
  arrange(FDR, .by_group = TRUE) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  mutate(
    Label = factor(Label, levels = rev(unique(Label))),
    StarX = Delta + ifelse(Delta >= 0, 0.045, -0.045)
  )

pC <- ggplot(iobr_dot, aes(Delta, Label)) +
  geom_vline(xintercept = 0, color = "#9BA7B7", linewidth = 0.35) +
  geom_point(aes(size = NegLogFDR, fill = Delta),
             shape = 21, color = "white", stroke = 0.45, alpha = 0.98) +
  geom_text(aes(x = StarX, label = Stars), size = 2.5, color = ink, fontface = "bold") +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = c1_col, mid = "white", high = c2_col, midpoint = 0,
                       name = "Median delta\nC2 - C1") +
  scale_size_continuous(range = c(2.2, 7.8), name = "-log10(FDR)") +
  labs(
    title = "Official IOBR immune and DDR programs across NMF subtypes",
    subtitle = "C2 recapitulates high-score DDR repair and immunosuppressive program activation",
    x = "Median difference, C2 - C1",
    y = NULL
  ) +
  theme_pub(8) +
  theme(
    strip.text.y = element_text(angle = 0, hjust = 0),
    panel.grid.major.y = element_blank(),
    legend.position = "right"
  )

save_plot(pC, "Gap_C_NMF_IOBR_program_dotplot", 8.5, 8.2)
save_csv(iobr_dot, file.path(OUT_DIR, "Gap_C_plot_data.csv"))

############################
# 7. Figure D: NMF TIDE / TMB / MATH and responder
############################

tide_features <- c("TIDE", "Dysfunction", "Exclusion", "MDSC", "CAF", "CD274", "CD8", "CTL")
tide_box <- tide_long %>%
  mutate(
    NMF_subtype = factor(NMF_subtype, levels = c("C1", "C2")),
    Score = num(Score)
  ) %>%
  filter(Feature %in% tide_features, is.finite(Score)) %>%
  mutate(Feature = factor(Feature, levels = tide_features))

box_long <- bind_rows(
  tmb %>% transmute(NMF_subtype = factor(NMF_subtype, levels = c("C1", "C2")), Feature = "TMB\n(log1p)", Value = num(TMB_log1p)),
  math %>% transmute(NMF_subtype = factor(NMF_subtype, levels = c("C1", "C2")), Feature = "MATH", Value = num(MATH)),
  tide_box %>% transmute(NMF_subtype, Feature = as.character(Feature), Value = Score)
) %>%
  filter(!is.na(NMF_subtype), is.finite(Value)) %>%
  mutate(Feature = factor(Feature, levels = c("TMB\n(log1p)", "MATH", tide_features)))

sig_df <- bind_rows(
  tmb_cmp %>% transmute(Feature = "TMB\n(log1p)", FDR = num(FDR), Stars = sig_star(FDR)),
  math_cmp %>% transmute(Feature = "MATH", FDR = num(FDR), Stars = sig_star(FDR)),
  tide_cmp %>% filter(Feature %in% tide_features) %>% transmute(Feature, FDR = num(FDR), Stars = sig_star(FDR))
) %>%
  mutate(Feature = factor(Feature, levels = levels(box_long$Feature)))

pD1 <- ggplot(box_long, aes(NMF_subtype, Value, fill = NMF_subtype)) +
  geom_violin(width = 0.86, alpha = 0.18, color = NA, trim = TRUE) +
  geom_boxplot(width = 0.22, outlier.shape = NA, linewidth = 0.35, color = "#273246", alpha = 0.85) +
  geom_jitter(aes(color = NMF_subtype), width = 0.11, size = 0.45, alpha = 0.24, show.legend = FALSE) +
  geom_text(data = sig_df, aes(x = 1.5, y = Inf, label = Stars),
            inherit.aes = FALSE, vjust = 1.18, size = 3.8, fontface = "bold", color = ink) +
  facet_wrap(~Feature, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c(C1 = c1_col, C2 = c2_col), drop = FALSE) +
  scale_color_manual(values = c(C1 = c1_col, C2 = c2_col), drop = FALSE) +
  labs(
    title = "Genomic instability and TIDE-derived features across NMF subtypes",
    x = NULL,
    y = "Value"
  ) +
  theme_pub(8) +
  theme(legend.position = "none", panel.grid.major.x = element_blank())

resp_p <- tide_cat_test %>% filter(CategoryVariable == "Responder") %>% pull(Fisher_P)
resp_plot <- tide_resp %>%
  mutate(
    NMF_subtype = factor(Group, levels = c("C1", "C2")),
    Category = ifelse(as.character(Category) == "TRUE", "Predicted responder", "Non-responder")
  )

pD2 <- ggplot(resp_plot, aes(NMF_subtype, Proportion, fill = Category)) +
  geom_col(width = 0.58, color = "white", linewidth = 0.35) +
  geom_text(aes(label = percent(Proportion, accuracy = 0.1)),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 3.1) +
  annotate("text", x = 1.5, y = 1.05, label = sig_star(resp_p),
           size = 5, fontface = "bold", color = ink) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1.1), expand = c(0, 0)) +
  scale_fill_manual(values = c("Predicted responder" = c2_col, "Non-responder" = "#AAB4C3")) +
  labs(
    title = "TIDE-predicted ICB response by NMF subtype",
    x = NULL,
    y = "Proportion",
    fill = NULL
  ) +
  theme_pub(8) +
  theme(panel.grid.major.x = element_blank(), legend.position = "bottom")

pD <- pD1 / pD2 + plot_layout(heights = c(1.4, 0.85))
save_plot(pD, "Gap_D_NMF_TIDE_TMB_MATH_response", 10.2, 7.4)
save_csv(box_long, file.path(OUT_DIR, "Gap_D_boxplot_data.csv"))

############################
# 8. Figure E: checkpoint panel effect dotplot
############################

cp_dot <- bind_rows(
  cp_mo %>%
    transmute(
      Comparison = "MO-DDRscore\nHigh vs Low",
      Gene,
      Delta = num(Delta_Group2_minus_Group1),
      FDR = num(FDR),
      Stars = sig_star(FDR)
    ),
  cp_nmf %>%
    transmute(
      Comparison = "NMF\nC2 vs C1",
      Gene,
      Delta = num(Delta_Group2_minus_Group1),
      FDR = num(FDR),
      Stars = sig_star(FDR)
    )
) %>%
  mutate(
    Comparison = factor(Comparison, levels = c("MO-DDRscore\nHigh vs Low", "NMF\nC2 vs C1")),
    Gene = factor(Gene, levels = rev(unique(cp_nmf$Gene))),
    NegLogFDR = neglog(FDR, cap = 18),
    Stars = ifelse(Stars == "ns", "", Stars)
  )

pE <- ggplot(cp_dot, aes(Comparison, Gene)) +
  geom_point(aes(size = NegLogFDR, fill = Delta),
             shape = 21, color = "white", stroke = 0.45, alpha = 0.98) +
  geom_text(aes(label = Stars), nudge_x = 0.24, size = 2.8, fontface = "bold", color = ink) +
  scale_fill_gradient2(low = c1_col, mid = "white", high = c2_col, midpoint = 0,
                       name = "Median delta") +
  scale_size_continuous(range = c(2.2, 8.2), name = "-log10(FDR)") +
  labs(
    title = "Immune checkpoint expression supports the high-score/C2 phenotype",
    subtitle = "Official IOBR immune checkpoint panel",
    x = NULL,
    y = NULL
  ) +
  theme_pub(9) +
  theme(panel.grid.major = element_line(color = "#EDF1F6", linewidth = 0.25))

save_plot(pE, "Gap_E_checkpoint_panel_effect_dotplot", 5.8, 4.8)
save_csv(cp_dot, file.path(OUT_DIR, "Gap_E_plot_data.csv"))

############################
# 9. Figure F: survival linkage
############################

make_km_df <- function(df, group_col) {
  d <- df %>%
    mutate(
      time = num(time),
      status = num(status),
      .group_chr = as.character(.data[[group_col]])
    ) %>%
    filter(is.finite(time), time > 0, status %in% c(0, 1), !is.na(.group_chr), .group_chr != "") %>%
    mutate(.group = as.factor(.group_chr))
  if (nrow(d) < 20 || nlevels(d$.group) < 2) return(NULL)
  fit <- survfit(Surv(time, status) ~ .group, data = d)
  s <- summary(fit)
  out <- data.frame(
    time = s$time,
    surv = s$surv,
    lower = s$lower,
    upper = s$upper,
    strata = s$strata,
    stringsAsFactors = FALSE
  )
  out$Group <- sub("^\\.group=", "", out$strata)
  out$GroupVariable <- group_col
  out
}

make_km_plot <- function(df, group_col, title, palette = NULL) {
  km <- make_km_df(df, group_col)
  if (is.null(km)) return(ggplot() + theme_void() + labs(title = title))
  pval <- surv_tests %>% filter(GroupVariable == group_col) %>% pull(Logrank_P)
  ptxt <- if (length(pval) == 0) "" else format_p(pval[1])
  if (is.null(palette)) {
    lev <- unique(km$Group)
    palette <- setNames(hue_pal()(length(lev)), lev)
  }
  ggplot(km, aes(time / 365, surv, color = Group)) +
    geom_step(linewidth = 0.75) +
    scale_color_manual(values = palette, drop = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = c(0.01, 0.01)) +
    coord_cartesian(xlim = c(0, max(km$time, na.rm = TRUE) / 365)) +
    annotate("text", x = Inf, y = 0.08, label = ptxt, hjust = 1.05, color = ink, fontface = "bold", size = 3.0) +
    labs(title = title, x = "Time, years", y = "Overall survival", color = NULL) +
    theme_pub(8) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

pal_mo <- c(Low = low_col, High = high_col)
pal_nmf <- c(C1 = c1_col, C2 = c2_col)
pal4_mo_tmb <- c(
  "Low_TMB-low" = "#8DB9DD",
  "Low_TMB-high" = "#2E6FAD",
  "High_TMB-low" = "#E69A8C",
  "High_TMB-high" = "#C94B3A"
)
pal4_mo_math <- c(
  "Low_MATH-low" = "#8DB9DD",
  "Low_MATH-high" = "#2E6FAD",
  "High_MATH-low" = "#E69A8C",
  "High_MATH-high" = "#C94B3A"
)

pF1 <- make_km_plot(surv_df, "MO_DDRscore_group", "MO-DDRscore group", pal_mo)
pF2 <- make_km_plot(surv_df, "NMF_subtype", "NMF subtype", pal_nmf)
pF3 <- make_km_plot(surv_df, "MO_TMB_group", "MO-DDRscore + TMB", pal4_mo_tmb)
pF4 <- make_km_plot(surv_df, "MO_MATH_group", "MO-DDRscore + MATH", pal4_mo_math)

pF <- (pF1 + pF2) / (pF3 + pF4) +
  plot_annotation(title = "Survival linkage of score-defined DDR state and genomic instability")

save_plot(pF, "Gap_F_survival_MO_NMF_TMB_MATH", 10.2, 7.2)

############################
# 10. Summary
############################

summary_df <- data.frame(
  Panel = c("A", "B", "C", "D", "E", "F"),
  Figure = c(
    "Gap_A_NMF_MO_DDRscore_concordance.pdf",
    "Gap_B_NMF_immune_deconvolution_lollipop.pdf",
    "Gap_C_NMF_IOBR_program_dotplot.pdf",
    "Gap_D_NMF_TIDE_TMB_MATH_response.pdf",
    "Gap_E_checkpoint_panel_effect_dotplot.pdf",
    "Gap_F_survival_MO_NMF_TMB_MATH.pdf"
  ),
  Meaning = c(
    "NMF subtype and MO-DDRscore group concordance",
    "NMF C2-minus-C1 immune infiltration differences",
    "NMF C2-minus-C1 IOBR immune and DDR program differences",
    "NMF subtype TIDE/TMB/MATH and predicted responder comparison",
    "Checkpoint expression effects in MO-DDRscore and NMF comparisons",
    "Survival linkage for MO-DDRscore, NMF, and TMB/MATH combinations"
  ),
  stringsAsFactors = FALSE
)

save_csv(summary_df, file.path(OUT_DIR, "Gap_publication_plot_summary.csv"))

cat("\nDone.\n")
cat("Output directory:\n", OUT_DIR, "\n")
print(summary_df)
