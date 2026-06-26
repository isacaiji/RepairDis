############################################################
# Figure 1C. Pan-cancer overall DDR pathway activity
# - Asterisks are recalculated from FDR, not raw P.
# - Two versions are exported:
#   1) complete data
#   2) data with Overall_DDR_activity < -5 removed
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(grid)
})

BASE_DIR <- "D:/R_workspace/summary/Repairium_pancancer_analysis_plan"
PLOT_ROOT <- file.path(BASE_DIR, "03-res", "plots", "张岩-图")
FIG_DIR <- file.path(PLOT_ROOT, "Figure1")
DATA_IN_DIR <- file.path(BASE_DIR, "03-res", "plots_data")
SELECTED_DIR <- file.path(DATA_IN_DIR, "main_figure_selected_data")
RAW_A_DIR <- file.path(BASE_DIR, "03-res", "Figure2_pan_cancer", "A_DDR_pathway_activity")
FIG_DATA_DIR <- file.path(FIG_DIR, "plot_data")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DATA_DIR, recursive = TRUE, showWarnings = FALSE)

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

safe_wilcox <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]
  g <- droplevels(g[ok])
  if (length(unique(g)) < 2) return(NA_real_)
  if (any(table(g) < 2)) return(NA_real_)
  suppressWarnings(stats::wilcox.test(x ~ g)$p.value)
}

read_activity <- function() {
  candidates <- c(
    file.path(DATA_IN_DIR, "Figure1D_overall_DDR_activity_long.csv"),
    file.path(SELECTED_DIR, "Figure1C_overall_DDR_activity_selected_long.csv")
  )
  in_file <- candidates[file.exists(candidates)][1]
  if (is.na(in_file)) {
    stop("Cannot find Figure1C overall DDR activity input file.")
  }

  x <- fread(in_file, data.table = FALSE, check.names = FALSE)
  required <- c("Cancer", "Sample", "MO_DDRscore_group", "Overall_DDR_activity")
  miss <- setdiff(required, colnames(x))
  if (length(miss) > 0) {
    stop("Input file lacks required columns: ", paste(miss, collapse = ", "))
  }

  x %>%
    mutate(
      Cancer = factor(Cancer, levels = cancer_order),
      MO_DDRscore_group = factor(MO_DDRscore_group, levels = c("Low", "High")),
      Overall_DDR_activity = as.numeric(Overall_DDR_activity)
    ) %>%
    filter(
      !is.na(Cancer),
      MO_DDRscore_group %in% c("Low", "High"),
      is.finite(Overall_DDR_activity)
    )
}

calc_stats <- function(x) {
  res <- lapply(levels(droplevels(x$Cancer)), function(ca) {
    d <- x %>% filter(Cancer == ca)
    g <- droplevels(d$MO_DDRscore_group)
    low <- d$Overall_DDR_activity[g == "Low"]
    high <- d$Overall_DDR_activity[g == "High"]
    data.frame(
      Cancer = ca,
      N_Low = sum(is.finite(low)),
      N_High = sum(is.finite(high)),
      Median_Low = median(low, na.rm = TRUE),
      Median_High = median(high, na.rm = TRUE),
      Effect = median(high, na.rm = TRUE) - median(low, na.rm = TRUE),
      P = safe_wilcox(d$Overall_DDR_activity, g),
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows() %>%
    mutate(
      FDR = p.adjust(P, method = "BH"),
      Significance = sig_from_fdr(FDR, ns_blank = FALSE),
      Direction = case_when(
        is.na(Effect) ~ "NA",
        Effect > 0 ~ "Higher in High",
        Effect < 0 ~ "Lower in High",
        TRUE ~ "No median difference"
      )
    )
  res
}

export_per_cancer <- function(x, stat, suffix) {
  for (ca in levels(droplevels(x$Cancer))) {
    ca_dir <- file.path(RAW_A_DIR, ca)
    dir.create(ca_dir, recursive = TRUE, showWarnings = FALSE)

    d_ca <- x %>% filter(Cancer == ca)
    s_ca <- stat %>% filter(Cancer == ca)

    fwrite(
      d_ca,
      file.path(ca_dir, paste0(ca, "_overall_DDR_activity_for_Figure1C_", suffix, ".csv"))
    )
    fwrite(
      s_ca,
      file.path(ca_dir, paste0(ca, "_overall_DDR_activity_for_Figure1C_", suffix, "_stats.csv"))
    )

    if (suffix == "complete") {
      fwrite(
        d_ca,
        file.path(ca_dir, paste0(ca, "_overall_DDR_activity_for_Figure1C.csv"))
      )
      fwrite(
        s_ca,
        file.path(ca_dir, paste0(ca, "_overall_DDR_activity_for_Figure1C_stats.csv"))
      )
    }
  }
}

save_plot <- function(p, stem, width = 14.5, height = 5.9) {
  pdf_file <- file.path(FIG_DIR, paste0(stem, ".pdf"))
  png_file <- file.path(FIG_DIR, paste0(stem, ".png"))

  tryCatch(
    ggsave(pdf_file, p, width = width, height = height, useDingbats = FALSE),
    error = function(e) {
      fallback <- file.path(FIG_DIR, paste0(stem, "_new.pdf"))
      message("PDF locked, writing fallback: ", fallback)
      ggsave(fallback, p, width = width, height = height, useDingbats = FALSE)
    }
  )

  tryCatch(
    ggsave(png_file, p, width = width, height = height, dpi = 600),
    error = function(e) {
      fallback <- file.path(FIG_DIR, paste0(stem, "_new.png"))
      message("PNG locked, writing fallback: ", fallback)
      ggsave(fallback, p, width = width, height = height, dpi = 600)
    }
  )
}

theme_figure1c <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", color = pal_dark, size = base_size + 5),
      plot.subtitle = element_text(color = "#6C7A89", size = base_size + 1),
      axis.title.x = element_text(face = "bold", color = pal_dark, margin = margin(t = 8)),
      axis.title.y = element_text(face = "bold", color = pal_dark, margin = margin(r = 8)),
      axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 8.4, face = "bold", color = pal_dark),
      axis.text.y = element_text(color = pal_dark),
      axis.line = element_line(color = pal_dark, linewidth = 0.45),
      axis.ticks = element_line(color = pal_dark, linewidth = 0.35),
      panel.grid.major.y = element_line(color = pal_grid, linewidth = 0.30),
      panel.grid.major.x = element_blank(),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_text(face = "bold", color = pal_dark),
      legend.text = element_text(color = pal_dark),
      plot.margin = margin(8, 10, 8, 12)
    )
}

draw_activity_plot <- function(x, stat, stem, subtitle_extra = NULL) {
  cancer_order_panel <- stat %>%
    mutate(Effect = as.numeric(Effect)) %>%
    filter(is.finite(Effect)) %>%
    arrange(Effect) %>%
    pull(Cancer) %>%
    as.character()

  x <- x %>%
    mutate(
      Cancer = factor(as.character(Cancer), levels = cancer_order_panel),
      MO_DDRscore_group = factor(as.character(MO_DDRscore_group), levels = c("Low", "High"))
    ) %>%
    filter(!is.na(Cancer))

  stat <- stat %>%
    mutate(Cancer = factor(Cancer, levels = cancer_order_panel)) %>%
    filter(!is.na(Cancer))

  y_rng <- range(x$Overall_DDR_activity, na.rm = TRUE)
  y_span <- diff(y_rng)
  if (!is.finite(y_span) || y_span == 0) y_span <- 1
  label_y <- y_rng[2] + 0.08 * y_span

  label_df <- stat %>%
    mutate(
      y = label_y,
      Significance = sig_from_fdr(FDR, ns_blank = FALSE)
    )

  subtitle <- "Mean standardized activity across 9 DDR pathway modules in MO-DDRscore-low and -high tumors"
  if (!is.null(subtitle_extra)) subtitle <- paste0(subtitle, "; ", subtitle_extra)

  p <- ggplot(x, aes(Cancer, Overall_DDR_activity, fill = MO_DDRscore_group)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.34, color = "#8EA1AF") +
    geom_violin(
      aes(color = MO_DDRscore_group),
      position = position_dodge(width = 0.78),
      width = 0.72,
      trim = TRUE,
      scale = "width",
      alpha = 0.78,
      linewidth = 0.32
    ) +
    geom_boxplot(
      aes(color = MO_DDRscore_group),
      position = position_dodge(width = 0.78),
      width = 0.13,
      outlier.shape = NA,
      alpha = 0.86,
      linewidth = 0.28
    ) +
    geom_text(
      data = label_df,
      aes(x = Cancer, y = y, label = Significance),
      inherit.aes = FALSE,
      angle = 0,
      size = 2.55,
      fontface = "bold",
      color = pal_dark
    ) +
    scale_fill_manual(values = c(Low = pal_low, High = pal_high), name = "MO-DDRscore") +
    scale_color_manual(values = c(Low = "#08777A", High = "#983F42"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.13))) +
    labs(
      title = "C  Pan-cancer overall DDR pathway activity",
      subtitle = subtitle,
      x = "Cancer type",
      y = "Overall DDR pathway activity score"
    ) +
    coord_cartesian(clip = "off") +
    theme_figure1c(10)

  save_plot(p, stem)
  invisible(p)
}

ddr_complete <- read_activity()
stat_complete <- calc_stats(ddr_complete)

ddr_no_below_minus5 <- ddr_complete %>%
  filter(Overall_DDR_activity >= -5)
stat_no_below_minus5 <- calc_stats(ddr_no_below_minus5)

fwrite(ddr_complete, file.path(FIG_DATA_DIR, "Figure1C_overall_DDR_activity_complete.csv"))
fwrite(stat_complete, file.path(FIG_DATA_DIR, "Figure1C_overall_DDR_activity_complete_stats.csv"))
fwrite(ddr_no_below_minus5, file.path(FIG_DATA_DIR, "Figure1C_overall_DDR_activity_noBelowMinus5.csv"))
fwrite(stat_no_below_minus5, file.path(FIG_DATA_DIR, "Figure1C_overall_DDR_activity_noBelowMinus5_stats.csv"))

export_per_cancer(ddr_complete, stat_complete, "complete")
export_per_cancer(ddr_no_below_minus5, stat_no_below_minus5, "noBelowMinus5")

writeLines(
  c(
    "Figure1C overall DDR pathway activity calculation:",
    "1. For each tumor sample, DDR pathway module activities were calculated upstream and summarized as Overall_DDR_activity.",
    "2. Overall_DDR_activity represents the mean standardized activity across 9 DDR pathway modules.",
    "3. Within each cancer type, tumors were grouped by MO-DDRscore median into Low and High groups.",
    "4. Group differences were tested using Wilcoxon rank-sum tests.",
    "5. Asterisks are based on BH-adjusted FDR: * <0.05, ** <0.01, *** <0.001, **** <0.0001.",
    "6. The noBelowMinus5 version removes samples with Overall_DDR_activity < -5 and recalculates Wilcoxon P values and FDR."
  ),
  con = file.path(FIG_DATA_DIR, "Figure1C_overall_DDR_activity_README.txt")
)

draw_activity_plot(
  ddr_complete,
  stat_complete,
  "Figure1C_overall_DDR_activity",
  subtitle_extra = "complete data"
)

draw_activity_plot(
  ddr_no_below_minus5,
  stat_no_below_minus5,
  "Figure1C_overall_DDR_activity_noBelowMinus5",
  subtitle_extra = "samples with activity < -5 removed"
)

fwrite(
  data.frame(
    Version = c("complete", "noBelowMinus5"),
    N_samples = c(nrow(ddr_complete), nrow(ddr_no_below_minus5)),
    N_removed_below_minus5 = c(0, nrow(ddr_complete) - nrow(ddr_no_below_minus5)),
    N_cancers = c(dplyr::n_distinct(ddr_complete$Cancer), dplyr::n_distinct(ddr_no_below_minus5$Cancer)),
    Significance_rule = "Asterisks are recalculated from BH-adjusted FDR: * <0.05, ** <0.01, *** <0.001, **** <0.0001.",
    stringsAsFactors = FALSE
  ),
  file.path(FIG_DATA_DIR, "Figure1C_overall_DDR_activity_run_summary.csv")
)

cat("Figure1C regenerated.\n")
cat("Script:\n", file.path(PLOT_ROOT, "01_C.R"), "\n")
cat("Figure output:\n", FIG_DIR, "\n")
cat("Data output:\n", FIG_DATA_DIR, "\n")
