############################################################
# Figure 3F: TIDE-predicted immunotherapy response
# Purpose:
#   Replot TIDE predicted responder proportions with cancers
#   ordered by High - Low responder proportion.
#   Significance stars are based on FDR, not raw P.
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

BASE_DIR <- "D:/R_workspace/summary/Repairium_pancancer_analysis_plan"
SRC_DIR <- file.path(BASE_DIR, "03-res", "Figure2_pan_cancer", "K_TIDE_response")
OUT_ROOT <- file.path(BASE_DIR, "03-res", "plots", "张岩-图")
FIG_DIR <- file.path(OUT_ROOT, "Figure3")
DATA_DIR <- file.path(OUT_ROOT, "plot_data")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

SUMMARY_FILE <- file.path(SRC_DIR, "K_TIDE_response_summary_all_cancers.csv")
STAT_FILE <- file.path(SRC_DIR, "K_TIDE_response_group_comparison_all_cancers.csv")

pal_low <- "#2CB7B8"
pal_high <- "#D75B68"
pal_dark <- "#102A43"
pal_grid <- "#E5ECF2"

sig_label_fdr <- function(fdr) {
  dplyr::case_when(
    is.na(fdr) ~ "",
    fdr < 1e-4 ~ "****",
    fdr < 1e-3 ~ "***",
    fdr < 1e-2 ~ "**",
    fdr < 5e-2 ~ "*",
    TRUE ~ ""
  )
}

resp_df <- data.table::fread(SUMMARY_FILE, data.table = FALSE) %>%
  dplyr::filter(StatusLabel == "Predicted responder") %>%
  dplyr::mutate(
    Group = factor(Group, levels = c("Low", "High")),
    Percent = as.numeric(Percent)
  ) %>%
  dplyr::select(Cancer, Group, N, Total, Percent)

stat_df <- data.table::fread(STAT_FILE, data.table = FALSE) %>%
  dplyr::mutate(
    P_value = as.numeric(P_value),
    FDR = p.adjust(P_value, method = "BH"),
    FDR_label = sig_label_fdr(FDR),
    Low_Responder_Percent = as.numeric(Low_Responder_Percent),
    High_Responder_Percent = as.numeric(High_Responder_Percent),
    Delta_High_minus_Low = High_Responder_Percent - Low_Responder_Percent
  ) %>%
  dplyr::select(
    Cancer, P_value, FDR, FDR_label,
    Low_Responder, Low_Total, Low_Responder_Percent,
    High_Responder, High_Total, High_Responder_Percent,
    Delta_High_minus_Low
  )

order_df <- stat_df %>%
  dplyr::arrange(dplyr::desc(High_Responder_Percent), dplyr::desc(Delta_High_minus_Low)) %>%
  dplyr::mutate(Cancer = as.character(Cancer))

cancer_order <- order_df$Cancer

plot_df <- resp_df %>%
  dplyr::left_join(stat_df, by = "Cancer") %>%
  dplyr::mutate(Cancer = factor(Cancer, levels = rev(cancer_order)))

seg_df <- stat_df %>%
  dplyr::mutate(Cancer = factor(Cancer, levels = rev(cancer_order)))

star_df <- stat_df %>%
  dplyr::filter(FDR_label != "") %>%
  dplyr::mutate(
    Cancer = factor(Cancer, levels = rev(cancer_order)),
    x_star = 103.5
  )

data.table::fwrite(
  plot_df,
  file.path(DATA_DIR, "Figure3F_TIDE_predicted_response_plot_data.csv")
)
data.table::fwrite(
  stat_df %>% dplyr::arrange(dplyr::desc(High_Responder_Percent), dplyr::desc(Delta_High_minus_Low)),
  file.path(DATA_DIR, "Figure3F_TIDE_predicted_response_statistics_FDR.csv")
)

p <- ggplot() +
  geom_segment(
    data = seg_df,
    aes(
      x = Low_Responder_Percent,
      xend = High_Responder_Percent,
      y = Cancer,
      yend = Cancer
    ),
    color = "#8D99A6",
    linewidth = 0.55,
    alpha = 0.90
  ) +
  geom_point(
    data = plot_df,
    aes(x = Percent, y = Cancer, color = Group),
    size = 2.35,
    alpha = 0.98
  ) +
  geom_text(
    data = star_df,
    aes(x = x_star, y = Cancer, label = FDR_label),
    color = pal_dark,
    fontface = "bold",
    size = 3.4,
    hjust = 0
  ) +
  scale_color_manual(
    values = c(Low = pal_low, High = pal_high),
    name = "MO-DDRscore"
  ) +
  scale_x_continuous(
    limits = c(0, 108),
    breaks = c(0, 25, 50, 75, 100),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "TIDE-predicted response",
    x = "Predicted responder proportion (%)",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = pal_grid, linewidth = 0.45),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "#68798A", linewidth = 0.55),
    axis.text = element_text(color = pal_dark, face = "bold"),
    axis.text.y = element_text(size = 8.6),
    axis.title.x = element_text(color = pal_dark, face = "bold", size = 12),
    legend.position = "top",
    legend.title = element_text(color = pal_dark, face = "bold"),
    legend.text = element_text(color = pal_dark),
    plot.title = element_text(color = pal_dark, face = "bold", size = 18, hjust = 0),
    plot.margin = margin(8, 48, 8, 8)
  )

ggsave(
  file.path(FIG_DIR, "Figure3F_TIDE_predicted_response.png"),
  p, width = 5.2, height = 6.2, dpi = 600
)
ggsave(
  file.path(FIG_DIR, "Figure3F_TIDE_predicted_response.pdf"),
  p, width = 5.2, height = 6.2, useDingbats = FALSE
)
ggsave(
  file.path(FIG_DIR, "Figure3F_TIDE_predicted_response.tiff"),
  p, width = 5.2, height = 6.2, dpi = 600, compression = "lzw"
)

legend_txt <- c(
  "Figure 3F. TIDE-predicted immunotherapy response.",
  "Predicted responder proportions were compared between MO-DDRscore-low and -high tumors within each cancer type.",
  "Dots indicate the predicted responder proportion in each group, and grey lines connect the paired low and high groups from the same cancer type.",
  "Cancer types are ordered by the predicted responder proportion in the MO-DDRscore-high group.",
  "Asterisks indicate FDR-adjusted significance from Fisher's exact test."
)
writeLines(
  legend_txt,
  con = file.path(DATA_DIR, "Figure3F_TIDE_predicted_response_legend.txt"),
  useBytes = TRUE
)

message("Done.")
message("Figure output: ", FIG_DIR)
message("Plot data output: ", DATA_DIR)
