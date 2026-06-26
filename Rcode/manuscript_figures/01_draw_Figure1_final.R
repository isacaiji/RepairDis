############################################################
# RepairDis Figure 1 final panels
# - all asterisks are generated from FDR, not raw P
# - Figure 1C data are exported into per-cancer folders
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(grid)
})

BASE_DIR <- "D:/R_workspace/summary/Repairium_pancancer_analysis_plan"
PLOT_DIR <- file.path(BASE_DIR, "03-res/plots/张岩-图")
FIG_DIR <- file.path(PLOT_DIR, "Figure1")
DATA_DIR <- file.path(BASE_DIR, "03-res/plots_data")
SEL_DIR <- file.path(DATA_DIR, "main_figure_selected_data")
RAW_A_DIR <- file.path(BASE_DIR, "03-res/Figure2_pan_cancer/A_DDR_pathway_activity")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

pal_low <- "#20AEB3"
pal_high <- "#D45F5F"
pal_dark <- "#10243C"
pal_grid <- "#E7EDF3"

cancer_order <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
  "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
  "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
  "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
)

sig_from_fdr <- function(fdr, ns_blank = FALSE) {
  out <- dplyr::case_when(
    is.na(fdr) ~ "",
    fdr < 1e-4 ~ "****",
    fdr < 1e-3 ~ "***",
    fdr < 1e-2 ~ "**",
    fdr < 5e-2 ~ "*",
    TRUE ~ "ns"
  )
  if (ns_blank) out[out == "ns"] <- ""
  out
}

save_plot <- function(p, stem, width, height) {
  pdf_file <- file.path(FIG_DIR, paste0(stem, ".pdf"))
  png_file <- file.path(FIG_DIR, paste0(stem, ".png"))
  tryCatch(
    ggsave(pdf_file, p, width = width, height = height, useDingbats = FALSE),
    error = function(e) {
      fallback <- file.path(FIG_DIR, paste0(stem, "_FDR_FIXED.pdf"))
      message("PDF locked, writing fallback: ", fallback)
      ggsave(fallback, p, width = width, height = height, useDingbats = FALSE)
    }
  )
  tryCatch(
    ggsave(png_file, p, width = width, height = height, dpi = 450, bg = "white"),
    error = function(e) {
      fallback <- file.path(FIG_DIR, paste0(stem, "_FDR_FIXED.png"))
      message("PNG locked, writing fallback: ", fallback)
      ggsave(fallback, p, width = width, height = height, dpi = 450, bg = "white")
    }
  )
}

theme_repair <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(color = pal_dark),
      plot.title = element_text(face = "bold", size = base_size + 5, hjust = 0),
      plot.subtitle = element_text(size = base_size + 0.5, color = "#617083", hjust = 0),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = pal_grid, linewidth = 0.28),
      panel.border = element_rect(color = "#CAD5E0", fill = NA, linewidth = 0.45),
      axis.text = element_text(color = "#263849"),
      axis.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "#EEF4F8", color = "#D7E1EA", linewidth = 0.45),
      strip.text = element_text(face = "bold", color = pal_dark),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(8, 12, 8, 8)
    )
}

############################################################
# Figure 1C. Overall DDR pathway activity
############################################################

ddr_long <- fread(file.path(DATA_DIR, "Figure1D_overall_DDR_activity_long.csv"), data.table = FALSE)
ddr_stat <- fread(file.path(DATA_DIR, "Figure1D_overall_DDR_activity_group_comparison.csv"), data.table = FALSE)

ddr_long <- ddr_long %>%
  mutate(
    Cancer = factor(Cancer, levels = cancer_order),
    MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High")),
    Overall_DDR_activity = as.numeric(Overall_DDR_activity)
  ) %>%
  filter(!is.na(Cancer), MO_DDRscore_group %in% c("Low", "High"), is.finite(Overall_DDR_activity))

ddr_stat <- ddr_stat %>%
  mutate(
    Cancer = factor(Cancer, levels = cancer_order),
    FDR = as.numeric(FDR),
    Significance = sig_from_fdr(FDR, ns_blank = FALSE)
  ) %>%
  filter(!is.na(Cancer))

for (ca in levels(droplevels(ddr_long$Cancer))) {
  ca_dir <- file.path(RAW_A_DIR, ca)
  dir.create(ca_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(
    ddr_long %>% filter(Cancer == ca),
    file.path(ca_dir, paste0(ca, "_overall_DDR_activity_for_Figure1C.csv"))
  )
  fwrite(
    ddr_stat %>% filter(Cancer == ca),
    file.path(ca_dir, paste0(ca, "_overall_DDR_activity_for_Figure1C_stats.csv"))
  )
}

y_rng <- range(ddr_long$Overall_DDR_activity, na.rm = TRUE)
y_span <- diff(y_rng)
label_df <- ddr_long %>%
  group_by(Cancer) %>%
  summarise(y = max(Overall_DDR_activity, na.rm = TRUE), .groups = "drop") %>%
  left_join(ddr_stat %>% select(Cancer, Significance), by = "Cancer") %>%
  mutate(y = y + 0.055 * y_span)

p1c <- ggplot(ddr_long, aes(Cancer, Overall_DDR_activity, fill = MO_DDRscore_group)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.34, color = "#8EA1AF") +
  geom_violin(aes(color = MO_DDRscore_group),
              position = position_dodge(width = 0.78), width = 0.72,
              trim = TRUE, scale = "width", alpha = 0.78, linewidth = 0.32) +
  geom_boxplot(aes(color = MO_DDRscore_group),
               position = position_dodge(width = 0.78), width = 0.13,
               outlier.shape = NA, alpha = 0.86, linewidth = 0.28) +
  geom_point(aes(group = MO_DDRscore_group),
             position = position_jitterdodge(jitter.width = 0.11, dodge.width = 0.78),
             size = 0.45, alpha = 0.35, color = "#202020") +
  geom_text(data = label_df, aes(x = Cancer, y = y, label = Significance),
            inherit.aes = FALSE, angle = 90, size = 2.7, fontface = "bold", color = pal_dark) +
  scale_fill_manual(values = c(Low = pal_low, High = pal_high), name = "MO-DDRscore") +
  scale_color_manual(values = c(Low = "#08777A", High = "#983F42"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.16))) +
  labs(
    title = "C  Pan-cancer overall DDR pathway activity",
    subtitle = "Mean standardized activity across DDR pathway modules in MO-DDRscore-low and -high tumors",
    x = NULL,
    y = "Overall DDR pathway activity"
  ) +
  theme_repair(10) +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 8.4, face = "bold"),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    legend.justification = "left"
  )

save_plot(p1c, "Figure1C_overall_DDR_activity", 14.5, 5.8)

############################################################
# Figure 1D. Cancer hallmark landscape
############################################################

hall <- fread(file.path(SEL_DIR, "Figure1D_cancer_hallmark_selected.csv"), data.table = FALSE) %>%
  mutate(
    Cancer = factor(Cancer, levels = cancer_order),
    NES = as.numeric(NES),
    FDR = as.numeric(FDR),
    Significance = sig_from_fdr(FDR, ns_blank = TRUE),
    Feature_Display = gsub("Dna", "DNA", Feature_Display),
    Feature_Display = gsub("E2f", "E2F", Feature_Display),
    Feature_Display = gsub("G2m", "G2M", Feature_Display),
    Feature_Display = gsub("Tgf", "TGF", Feature_Display),
    Feature_Display = gsub("Mtorc1", "MTORC1", Feature_Display),
    Feature_Display = gsub("Myc", "MYC", Feature_Display),
    Feature_Display = gsub("Pi3k Akt Mtor", "PI3K-AKT-MTOR", Feature_Display),
    Feature_Display = factor(Feature_Display, levels = rev(unique(Feature_Display)))
  ) %>%
  filter(!is.na(Cancer), is.finite(NES), is.finite(FDR))

p1d <- ggplot(hall, aes(Cancer, Feature_Display, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.52, width = 0.94, height = 0.90) +
  geom_text(aes(label = Significance), size = 2.35, fontface = "bold", color = pal_dark) +
  scale_fill_gradient(low = "#FFF2EC", high = pal_high, name = "NES") +
  labs(
    title = "D  Cancer hallmark programs",
    subtitle = "MSigDB Hallmark GSEA in MO-DDRscore-high versus -low tumors",
    x = "Cancer type",
    y = NULL
  ) +
  theme_repair(9.5) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 8.2, face = "bold"),
    axis.text.y = element_text(size = 8.2, face = "bold"),
    legend.key.height = unit(18, "pt")
  )

save_plot(p1d, "Figure1D_cancer_hallmark_landscape", 9.6, 5.2)

############################################################
# Figure 1E. Aging-related hallmark landscape
############################################################

aging_order <- c(
  "Genomic instability", "Telomere attrition", "Epigenetic alterations",
  "Loss of proteostasis", "Disabled macroautophagy",
  "Deregulated nutrient sensing", "Mitochondrial dysfunction",
  "Cellular senescence", "Stem cell exhaustion",
  "Altered intercellular communication", "Chronic inflammation",
  "Extracellular matrix changes"
)

aging <- fread(file.path(SEL_DIR, "Figure1E_aging_hallmark_selected.csv"), data.table = FALSE) %>%
  mutate(
    Cancer = factor(Cancer, levels = cancer_order),
    Feature = factor(Feature, levels = rev(aging_order)),
    FDR = as.numeric(FDR),
    EffectScaledWithinFeature = as.numeric(EffectScaledWithinFeature),
    Significance = sig_from_fdr(FDR, ns_blank = TRUE),
    PlotValue = pmax(pmin(EffectScaledWithinFeature, 2.8), -2.8)
  ) %>%
  filter(!is.na(Cancer), !is.na(Feature), is.finite(PlotValue))

p1e <- ggplot(aging, aes(Cancer, Feature, fill = PlotValue)) +
  geom_tile(color = "white", linewidth = 0.50, width = 0.94, height = 0.90) +
  geom_text(aes(label = Significance), size = 2.2, fontface = "bold", color = pal_dark) +
  scale_fill_gradient2(low = pal_low, mid = "white", high = pal_high,
                       midpoint = 0, limits = c(-2.8, 2.8), oob = scales::squish,
                       name = "High-Low\nscaled effect") +
  labs(
    title = "E  Aging-related hallmark landscape",
    subtitle = "Cancer-specific high-low differences in aging-related pathway scores",
    x = "Cancer type",
    y = NULL
  ) +
  theme_repair(9.3) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 8.0, face = "bold"),
    axis.text.y = element_text(size = 8.1, face = "bold"),
    legend.key.height = unit(18, "pt")
  )

save_plot(p1e, "Figure1E_aging_hallmark_landscape", 10.2, 5.4)

############################################################
# Figure 1F. TMB/MATH, significant results only
############################################################

tmb <- fread(file.path(SEL_DIR, "Figure1F_TMB_MATH_selected_stats.csv"), data.table = FALSE) %>%
  mutate(
    Cancer = factor(Cancer, levels = cancer_order),
    Feature_Display = factor(Feature_Display, levels = c("TMB", "MATH")),
    Effect = as.numeric(Effect),
    FDR = as.numeric(FDR),
    Significance = sig_from_fdr(FDR, ns_blank = TRUE),
    DirectionCol = ifelse(Effect >= 0, "Higher in High", "Higher in Low")
  ) %>%
  filter(!is.na(Cancer), !is.na(Feature_Display), is.finite(Effect), is.finite(FDR), FDR < 0.05)

keep_cancers <- tmb %>%
  count(Cancer, name = "n_sig") %>%
  arrange(Cancer) %>%
  pull(Cancer) %>%
  as.character()

tmb <- tmb %>%
  mutate(Cancer = factor(as.character(Cancer), levels = keep_cancers))

p1f <- ggplot(tmb, aes(y = Cancer, x = Effect, color = DirectionCol)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, color = "#AAB7C2") +
  geom_segment(aes(x = 0, xend = Effect, yend = Cancer), linewidth = 0.85, alpha = 0.78) +
  geom_point(size = 2.5, alpha = 0.95) +
  geom_text(aes(label = Significance), nudge_x = 0.08, color = pal_dark,
            fontface = "bold", size = 2.8) +
  facet_grid(Feature_Display ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = c("Higher in High" = pal_high, "Higher in Low" = pal_low), name = NULL) +
  labs(
    title = "F  Tumor mutation burden and MATH",
    subtitle = "Only FDR-significant TMB/MATH differences are shown",
    x = "Median difference (High - Low)",
    y = NULL
  ) +
  theme_repair(9.5) +
  theme(
    axis.text.y = element_text(size = 8.5, face = "bold"),
    legend.position = "top",
    strip.text.y = element_text(angle = 0)
  )

save_plot(p1f, "Figure1F_TMB_MATH", 7.8, 4.8)

fwrite(
  data.frame(
    Panel = c("Figure1C", "Figure1D", "Figure1E", "Figure1F"),
    Significance_rule = "Asterisks were recalculated from FDR: * <0.05, ** <0.01, *** <0.001, **** <0.0001.",
    stringsAsFactors = FALSE
  ),
  file.path(FIG_DIR, "Figure1_significance_rule.csv")
)

cat("Figure 1 final panels regenerated in:\n", FIG_DIR, "\n")
