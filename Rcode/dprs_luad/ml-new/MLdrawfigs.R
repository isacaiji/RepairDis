############################################################
# 02_ML_DPRS_Figure4_plotting.R
#
# Purpose:
#   Replot Figure 4 from ML-DPRS modeling outputs.
#
# Input directory:
#   D:/R_workspace/评分/AD_DDR_project/04-ML
#
# Required input files:
#   Fig4B_all_model_Cindex_summary.csv
#   Fig4B_final_selected_model_info.csv
#   Fig4C_ML_DDRscore_all_sets.csv
#   Fig4C_final_signature_genes.csv
#
# Output:
#   04-ML/02_Figure4_plots/figures
#   04-ML/02_Figure4_plots/tables
############################################################

options(stringsAsFactors = FALSE)

############################
# 0. Paths
############################

ML_DIR <- "D:/R_workspace/评分/AD_DDR_project/04-ML"

OUT_DIR <- file.path(ML_DIR, "02_Figure4_plots")
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

############################
# 1. Packages
############################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(survival)
  library(grid)
  library(ComplexHeatmap)
  library(circlize)
})

############################
# 2. Helper functions
############################

save_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file)
}

save_pdf <- function(file, p, w = 7, h = 5) {
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

calc_cindex <- function(time, status, score) {
  df <- data.frame(
    time = as.numeric(time),
    status = as.numeric(status),
    score = as.numeric(score)
  )
  
  df <- df[
    is.finite(df$time) &
      is.finite(df$status) &
      is.finite(df$score),
  ]
  
  if (nrow(df) < 30 || sum(df$status == 1) < 5) {
    return(NA_real_)
  }
  
  as.numeric(
    survival::concordance(
      survival::Surv(time, status) ~ score,
      data = df,
      reverse = TRUE
    )$concordance
  )
}

get_time_points <- function(time) {
  if (median(time, na.rm = TRUE) > 100) {
    c(365, 1095, 1825)
  } else {
    c(1, 3, 5)
  }
}

############################
# 3. Read modeling outputs
############################

cindex_file <- file.path(ML_DIR, "Fig4B_all_model_Cindex_summary.csv")
selected_file <- file.path(ML_DIR, "Fig4B_final_selected_model_info.csv")
score_file <- file.path(ML_DIR, "Fig4C_ML_DDRscore_all_sets.csv")
gene_file <- file.path(ML_DIR, "Fig4C_final_signature_genes.csv")

if (!file.exists(cindex_file)) stop("Cannot find: ", cindex_file)
if (!file.exists(selected_file)) stop("Cannot find: ", selected_file)
if (!file.exists(score_file)) stop("Cannot find: ", score_file)
if (!file.exists(gene_file)) stop("Cannot find: ", gene_file)

cindex_summary <- data.table::fread(cindex_file, data.table = FALSE, check.names = FALSE)
selected_info <- data.table::fread(selected_file, data.table = FALSE, check.names = FALSE)
final_score <- data.table::fread(score_file, data.table = FALSE, check.names = FALSE)
signature_genes <- data.table::fread(gene_file, data.table = FALSE, check.names = FALSE)

best_model <- as.character(selected_info$Model[1])

cat("\nSelected model:\n")
print(best_model)

dataset_cols <- intersect(
  c("Training", "Testing", "GSE72094", "GSE68465"),
  colnames(cindex_summary)
)

if (length(dataset_cols) < 2) {
  stop("Cannot find enough dataset columns in C-index summary.")
}

############################
# 4. Prepare C-index summary
############################

clean_cindex_summary <- function(cindex_summary, best_model, dataset_cols) {
  
  cindex_summary %>%
    dplyr::mutate(
      dplyr::across(dplyr::all_of(dataset_cols), as.numeric),
      Average = rowMeans(
        dplyr::across(dplyr::all_of(dataset_cols)),
        na.rm = TRUE
      ),
      Train_Test_Gap = if (all(c("Training", "Testing") %in% dataset_cols)) {
        abs(Training - Testing)
      } else {
        NA_real_
      },
      IsSelected = Model == best_model
    ) %>%
    dplyr::group_by(Model) %>%
    dplyr::arrange(
      dplyr::desc(Average),
      dplyr::desc(IsSelected),
      Train_Test_Gap,
      .by_group = TRUE
    ) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(
      dplyr::desc(Average),
      dplyr::desc(IsSelected),
      Train_Test_Gap
    )
}

cindex_clean <- clean_cindex_summary(
  cindex_summary = cindex_summary,
  best_model = best_model,
  dataset_cols = dataset_cols
)

save_csv(
  cindex_clean,
  file.path(TAB_DIR, "Fig4B_clean_Cindex_summary.csv")
)

############################
# 5. C-index heatmap
############################

plot_cindex_heatmap_final <- function(cindex_summary,
                                      best_model,
                                      dataset_cols,
                                      top_n = 30,
                                      out_pdf,
                                      out_csv) {
  
  dat <- clean_cindex_summary(
    cindex_summary = cindex_summary,
    best_model = best_model,
    dataset_cols = dataset_cols
  )
  
  if (!is.null(top_n)) {
    dat <- dat %>% dplyr::slice_head(n = top_n)
  }
  
  save_csv(dat, out_csv)
  
  mat <- dat %>%
    dplyr::select(dplyr::all_of(dataset_cols)) %>%
    as.matrix()
  
  rownames(mat) <- dat$Model
  colnames(mat) <- dataset_cols
  
  n_model <- nrow(mat)
  
  if (n_model <= 15) {
    row_fs <- 8.6
    cell_fs <- 8.2
    col_fs <- 8.5
    row_h_cm <- 0.46
    pdf_h <- 7.2
    pdf_w <- 10.6
    heat_w_cm <- 6.2
    bar_w_cm <- 2.4
    bar_num_fs <- 7.2
  } else if (n_model <= 30) {
    row_fs <- 7.3
    cell_fs <- 7.0
    col_fs <- 8.2
    row_h_cm <- 0.36
    pdf_h <- 9.4
    pdf_w <- 11.0
    heat_w_cm <- 6.3
    bar_w_cm <- 2.4
    bar_num_fs <- 6.6
  } else if (n_model <= 60) {
    row_fs <- 6.0
    cell_fs <- 5.6
    col_fs <- 8.0
    row_h_cm <- 0.28
    pdf_h <- max(10.5, 0.12 * n_model + 4.2)
    pdf_w <- 11.4
    heat_w_cm <- 6.4
    bar_w_cm <- 2.5
    bar_num_fs <- 5.4
  } else {
    row_fs <- 5.1
    cell_fs <- 4.8
    col_fs <- 7.6
    row_h_cm <- 0.24
    pdf_h <- max(16, 0.115 * n_model + 4.5)
    pdf_w <- 11.8
    heat_w_cm <- 6.4
    bar_w_cm <- 2.5
    bar_num_fs <- 4.8
  }
  
  col_fun <- circlize::colorRamp2(
    c(0.50, 0.65, 0.85),
    c("#4195C1", "#FFFFFF", "#CB5746")
  )
  
  row_label_col <- ifelse(dat$IsSelected, "#C53030", "black")
  row_label_face <- ifelse(dat$IsSelected, "bold", "plain")
  
  avg_for_bar <- as.numeric(dat$Average)
  avg_label <- sprintf("%.3f", avg_for_bar)
  bar_fill <- ifelse(dat$IsSelected, "#C53030", "steelblue")
  
  avg_bar_anno <- ComplexHeatmap::AnnotationFunction(
    which = "row",
    n = nrow(dat),
    width = grid::unit(bar_w_cm, "cm"),
    fun = function(index, k, n) {
      
      vals <- avg_for_bar[index]
      labs <- avg_label[index]
      fills <- bar_fill[index]
      nr <- length(index)
      
      grid::pushViewport(
        grid::viewport(
          xscale = c(0, 0.85),
          yscale = c(0.5, nr + 0.5),
          clip = "on"
        )
      )
      
      y_pos <- rev(seq_len(nr))
      
      grid::grid.rect(
        x = grid::unit(0, "native"),
        y = grid::unit(y_pos, "native"),
        width = grid::unit(vals, "native"),
        height = grid::unit(0.68, "native"),
        just = c("left", "center"),
        gp = grid::gpar(fill = fills, col = NA)
      )
      
      text_x <- pmax(vals - 0.012, 0.05)
      
      grid::grid.text(
        label = labs,
        x = grid::unit(text_x, "native"),
        y = grid::unit(y_pos, "native"),
        just = "right",
        gp = grid::gpar(
          fontsize = bar_num_fs,
          col = "white",
          fontface = "bold"
        )
      )
      
      grid::grid.xaxis(
        at = c(0, 0.5, 0.7),
        label = c("0", "0.5", "0.7"),
        gp = grid::gpar(fontsize = 6)
      )
      
      grid::popViewport()
    }
  )
  
  right_ha <- ComplexHeatmap::rowAnnotation(
    `Average C-index` = avg_bar_anno,
    annotation_name_side = "top",
    annotation_name_rot = 0,
    annotation_name_gp = grid::gpar(
      fontsize = 6.5,
      fontface = "bold"
    )
  )
  
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = "C-index",
    col = col_fun,
    right_annotation = right_ha,
    
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    rect_gp = grid::gpar(col = "black", lwd = 0.65),
    
    show_row_names = TRUE,
    show_column_names = TRUE,
    
    row_names_side = "left",
    row_names_max_width = grid::unit(6.8, "cm"),
    row_names_gp = grid::gpar(
      fontsize = row_fs,
      col = row_label_col,
      fontface = row_label_face
    ),
    
    column_names_side = "top",
    column_names_rot = 0,
    column_names_centered = TRUE,
    column_names_gp = grid::gpar(
      fontsize = col_fs,
      fontface = "bold",
      col = "black"
    ),
    
    heatmap_legend_param = list(
      title = "C-index",
      title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 9),
      at = c(0.50, 0.60, 0.70, 0.80),
      labels = c("0.50", "0.60", "0.70", "0.80"),
      legend_height = grid::unit(3.6, "cm")
    ),
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      v <- mat[i, j]
      
      if (dat$IsSelected[i]) {
        grid::grid.rect(
          x = x,
          y = y,
          width = width,
          height = height,
          gp = grid::gpar(fill = NA, col = "#C53030", lwd = 1.35)
        )
      }
      
      grid::grid.text(
        sprintf("%.3f", v),
        x,
        y,
        gp = grid::gpar(
          fontsize = cell_fs,
          col = "black",
          fontface = ifelse(dat$IsSelected[i], "bold", "plain")
        )
      )
    },
    
    column_title = "Machine-learning model performance",
    column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
    
    width = grid::unit(heat_w_cm, "cm"),
    height = grid::unit(row_h_cm * n_model, "cm")
  )
  
  pdf(
    out_pdf,
    width = pdf_w,
    height = pdf_h,
    useDingbats = FALSE
  )
  
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = grid::unit(c(5, 5, 5, 5), "mm")
  )
  
  dev.off()
  
  invisible(dat)
}

cindex_top15 <- plot_cindex_heatmap_final(
  cindex_summary = cindex_summary,
  best_model = best_model,
  dataset_cols = dataset_cols,
  top_n = 15,
  out_pdf = file.path(FIG_DIR, "Fig4B_Cindex_heatmap_top15.pdf"),
  out_csv = file.path(TAB_DIR, "Fig4B_Cindex_heatmap_top15.csv")
)

cindex_top30 <- plot_cindex_heatmap_final(
  cindex_summary = cindex_summary,
  best_model = best_model,
  dataset_cols = dataset_cols,
  top_n = 30,
  out_pdf = file.path(FIG_DIR, "Fig4B_Cindex_heatmap_top30.pdf"),
  out_csv = file.path(TAB_DIR, "Fig4B_Cindex_heatmap_top30.csv")
)

cindex_all <- plot_cindex_heatmap_final(
  cindex_summary = cindex_summary,
  best_model = best_model,
  dataset_cols = dataset_cols,
  top_n = NULL,
  out_pdf = file.path(FIG_DIR, "Fig4B_Cindex_heatmap_all_models.pdf"),
  out_csv = file.path(TAB_DIR, "Fig4B_Cindex_heatmap_all_models.csv")
)

############################
# 6. Average C-index barplot
############################

top_bar <- cindex_clean %>%
  dplyr::slice_head(n = 30) %>%
  dplyr::mutate(
    DisplayModel = Model,
    DisplayModel = factor(DisplayModel, levels = rev(DisplayModel)),
    Group = ifelse(IsSelected, "Selected model", "Other models")
  )

p_avg <- ggplot(top_bar, aes(x = Average, y = DisplayModel, fill = Group)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = sprintf("%.3f", Average)),
    hjust = -0.08,
    size = 3.0,
    fontface = "bold"
  ) +
  scale_fill_manual(values = c("Selected model" = "#C53030", "Other models" = "grey75")) +
  geom_vline(xintercept = 0.65, linetype = 2, color = "grey45") +
  coord_cartesian(xlim = c(0.50, max(top_bar$Average, na.rm = TRUE) + 0.07)) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(color = "black", size = 7.5),
    axis.text.x = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "Average C-index",
    y = NULL,
    title = "Top machine-learning models ranked by average C-index"
  )

save_pdf(
  file.path(FIG_DIR, "Fig4B_average_Cindex_barplot_top30.pdf"),
  p_avg,
  8.5,
  7
)

############################
# 7. Prepare DPRS score table
############################

final_score <- final_score %>%
  dplyr::mutate(
    Dataset = as.character(Dataset),
    time = as.numeric(time),
    status = as.numeric(status),
    ML_DDRscore = as.numeric(ML_DDRscore)
  ) %>%
  dplyr::filter(
    is.finite(time),
    time > 0,
    status %in% c(0, 1),
    is.finite(ML_DDRscore)
  )

if (!"RiskGroup" %in% colnames(final_score)) {
  final_score <- final_score %>%
    dplyr::group_by(Dataset) %>%
    dplyr::mutate(
      Cutoff = median(ML_DDRscore, na.rm = TRUE),
      RiskGroup = ifelse(ML_DDRscore >= Cutoff, "High", "Low")
    ) %>%
    dplyr::ungroup()
}

final_score$RiskGroup <- factor(final_score$RiskGroup, levels = c("Low", "High"))
final_score$Dataset <- factor(
  final_score$Dataset,
  levels = c("Training", "Testing", "GSE72094", "GSE68465")
)

save_csv(
  final_score,
  file.path(TAB_DIR, "Fig4_replot_DPRS_scores.csv")
)

score_cindex <- final_score %>%
  dplyr::group_by(Dataset) %>%
  dplyr::summarise(
    N = dplyr::n(),
    Events = sum(status == 1),
    Cindex = calc_cindex(time, status, ML_DDRscore),
    .groups = "drop"
  )

save_csv(
  score_cindex,
  file.path(TAB_DIR, "Fig4D_selected_model_Cindex_by_dataset.csv")
)

############################
# 8. Risk score and survival status distribution
############################

risk_rank <- final_score %>%
  dplyr::filter(is.finite(ML_DDRscore)) %>%
  dplyr::arrange(Dataset, ML_DDRscore) %>%
  dplyr::group_by(Dataset) %>%
  dplyr::mutate(Rank = dplyr::row_number()) %>%
  dplyr::ungroup()

p_risk <- ggplot(risk_rank, aes(x = Rank, y = ML_DDRscore, color = RiskGroup)) +
  geom_point(size = 0.75, alpha = 0.9) +
  facet_wrap(~Dataset, scales = "free_x", nrow = 1) +
  scale_color_manual(values = c("Low" = "#2B6CB0", "High" = "#C53030")) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    legend.title = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = "Patient rank",
    y = "DPRS",
    title = paste0("Risk score distribution: ", best_model)
  )

save_pdf(
  file.path(FIG_DIR, "Fig4C_DPRS_risk_score_distribution.pdf"),
  p_risk,
  12,
  4
)

p_status <- ggplot(risk_rank, aes(x = Rank, y = time, color = factor(status))) +
  geom_point(size = 0.75, alpha = 0.9) +
  facet_wrap(~Dataset, scales = "free_x", nrow = 1) +
  scale_color_manual(
    values = c("0" = "#2B6CB0", "1" = "#C53030"),
    labels = c("Alive/Censored", "Dead")
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    legend.title = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = "Patient rank",
    y = "Survival time",
    color = "Status",
    title = "Survival status distribution"
  )

save_pdf(
  file.path(FIG_DIR, "Fig4C_DPRS_survival_status_distribution.pdf"),
  p_status,
  12,
  4
)

############################
# 9. KM plots
############################

if (requireNamespace("survminer", quietly = TRUE)) {
  
  for (ds in levels(final_score$Dataset)) {
    
    df <- final_score %>% dplyr::filter(Dataset == ds)
    if (nrow(df) < 40 || sum(df$status == 1) < 8) next
    
    fit <- survival::survfit(
      survival::Surv(time, status) ~ RiskGroup,
      data = df
    )
    
    p_km <- survminer::ggsurvplot(
      fit,
      data = df,
      pval = TRUE,
      risk.table = TRUE,
      risk.table.height = 0.23,
      palette = c("#2B6CB0", "#C53030"),
      legend.title = "",
      legend.labs = c("Low DPRS", "High DPRS"),
      title = paste0(ds, " cohort"),
      xlab = "Time",
      ylab = "Overall survival probability",
      ggtheme = theme_bw(base_size = 12)
    )
    
    pdf(
      file.path(FIG_DIR, paste0("Fig4D_DPRS_KM_", ds, ".pdf")),
      width = 5.8,
      height = 6.3,
      useDingbats = FALSE
    )
    
    print(p_km)
    dev.off()
  }
} else {
  message("Package survminer is not installed. Skip KM plots.")
}

############################
# 10. timeROC plots
############################

if (requireNamespace("timeROC", quietly = TRUE)) {
  
  auc_all <- list()
  
  for (ds in levels(final_score$Dataset)) {
    
    df <- final_score %>% dplyr::filter(Dataset == ds)
    if (nrow(df) < 40 || sum(df$status == 1) < 8) next
    
    times <- get_time_points(df$time)
    time_labels <- c("1-year", "3-year", "5-year")
    
    roc <- timeROC::timeROC(
      T = df$time,
      delta = df$status,
      marker = df$ML_DDRscore,
      cause = 1,
      weighting = "marginal",
      times = times,
      ROC = TRUE
    )
    
    auc_df <- data.frame(
      Dataset = ds,
      Model = best_model,
      Time = time_labels,
      TimeValue = times,
      AUC = as.numeric(roc$AUC),
      stringsAsFactors = FALSE
    )
    
    auc_all[[ds]] <- auc_df
    
    save_csv(
      auc_df,
      file.path(TAB_DIR, paste0("Fig4E_DPRS_timeROC_AUC_", ds, ".csv"))
    )
    
    roc_df <- data.frame(
      FP_1y = roc$FP[, 1],
      TP_1y = roc$TP[, 1],
      FP_3y = roc$FP[, 2],
      TP_3y = roc$TP[, 2],
      FP_5y = roc$FP[, 3],
      TP_5y = roc$TP[, 3]
    )
    
    p_roc <- ggplot() +
      geom_line(data = roc_df, aes(x = FP_1y, y = TP_1y), linewidth = 1.1, color = "#EA6433") +
      geom_line(data = roc_df, aes(x = FP_3y, y = TP_3y), linewidth = 1.1, color = "#16A5D9") +
      geom_line(data = roc_df, aes(x = FP_5y, y = TP_5y), linewidth = 1.1, color = "#FFC21A") +
      geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey65", linewidth = 0.8) +
      annotate(
        "text",
        x = 0.60,
        y = 0.24,
        size = 4.1,
        hjust = 0,
        label = paste0("1-year AUC = ", sprintf("%.3f", auc_df$AUC[1])),
        color = "#EA6433",
        fontface = "bold"
      ) +
      annotate(
        "text",
        x = 0.60,
        y = 0.15,
        size = 4.1,
        hjust = 0,
        label = paste0("3-year AUC = ", sprintf("%.3f", auc_df$AUC[2])),
        color = "#16A5D9",
        fontface = "bold"
      ) +
      annotate(
        "text",
        x = 0.60,
        y = 0.06,
        size = 4.1,
        hjust = 0,
        label = paste0("5-year AUC = ", sprintf("%.3f", auc_df$AUC[3])),
        color = "#D99A00",
        fontface = "bold"
      ) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black", face = "bold"),
        panel.border = element_rect(color = "black", linewidth = 0.8),
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      labs(
        x = "False positive rate",
        y = "True positive rate",
        title = paste0(ds, " DPRS time-dependent ROC")
      )
    
    save_pdf(
      file.path(FIG_DIR, paste0("Fig4E_DPRS_timeROC_", ds, ".pdf")),
      p_roc,
      5.5,
      5.2
    )
  }
  
  auc_all <- dplyr::bind_rows(auc_all)
  
  save_csv(
    auc_all,
    file.path(TAB_DIR, "Fig4E_DPRS_timeROC_AUC_all.csv")
  )
  
  if (nrow(auc_all) > 0) {
    p_auc_sum <- ggplot(auc_all, aes(x = Time, y = AUC, fill = Dataset)) +
      geom_col(
        position = position_dodge(width = 0.75),
        width = 0.65,
        color = "black",
        linewidth = 0.2
      ) +
      geom_text(
        aes(label = sprintf("%.3f", AUC)),
        position = position_dodge(width = 0.75),
        vjust = -0.25,
        size = 3
      ) +
      coord_cartesian(ylim = c(0.45, max(auc_all$AUC, na.rm = TRUE) + 0.08)) +
      theme_classic(base_size = 12) +
      theme(
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black", face = "bold"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      labs(
        x = NULL,
        y = "AUC",
        title = "DPRS time-dependent AUC across cohorts"
      )
    
    save_pdf(
      file.path(FIG_DIR, "Fig4E_DPRS_timeROC_AUC_summary.pdf"),
      p_auc_sum,
      7,
      5
    )
  }
  
} else {
  message("Package timeROC is not installed. Skip time-dependent ROC.")
}

############################
# 11. Final signature gene coefficient plot
############################

coef_col <- intersect(
  c("StepCoxCoef", "CoxCoef", "coef", "Coefficient", "coefficients"),
  colnames(signature_genes)
)[1]

if (!is.na(coef_col) && "Gene" %in% colnames(signature_genes)) {
  
  gene_plot <- signature_genes %>%
    dplyr::mutate(
      Gene = as.character(Gene),
      Coef = as.numeric(.data[[coef_col]]),
      Direction = ifelse(Coef >= 0, "Risk", "Protective")
    ) %>%
    dplyr::filter(is.finite(Coef)) %>%
    dplyr::arrange(Coef)
  
  p_gene <- ggplot(
    gene_plot,
    aes(x = reorder(Gene, Coef), y = Coef, color = Direction)
  ) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_segment(
      aes(xend = Gene, y = 0, yend = Coef),
      linewidth = 0.8,
      color = "grey65"
    ) +
    geom_point(size = 2.8) +
    scale_color_manual(values = c("Risk" = "#C53030", "Protective" = "#2B6CB0")) +
    coord_flip() +
    theme_classic(base_size = 12) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold"),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      x = NULL,
      y = ifelse(coef_col == "StepCoxCoef", "StepCox coefficient", "Training Cox coefficient"),
      title = paste0("Final signature genes of ", best_model)
    )
  
  save_pdf(
    file.path(FIG_DIR, "Fig4F_DPRS_signature_gene_lollipop.pdf"),
    p_gene,
    6.5,
    max(4.8, 0.28 * nrow(gene_plot) + 2)
  )
  
  save_csv(
    gene_plot,
    file.path(TAB_DIR, "Fig4F_DPRS_signature_gene_lollipop_data.csv")
  )
}

############################
# 12. Output manifest
############################

plot_manifest <- data.frame(
  File = c(
    "Fig4B_Cindex_heatmap_top15.pdf",
    "Fig4B_Cindex_heatmap_top30.pdf",
    "Fig4B_Cindex_heatmap_all_models.pdf",
    "Fig4B_average_Cindex_barplot_top30.pdf",
    "Fig4C_DPRS_risk_score_distribution.pdf",
    "Fig4C_DPRS_survival_status_distribution.pdf",
    "Fig4D_DPRS_KM_*.pdf",
    "Fig4E_DPRS_timeROC_*.pdf",
    "Fig4E_DPRS_timeROC_AUC_summary.pdf",
    "Fig4F_DPRS_signature_gene_lollipop.pdf"
  ),
  Description = c(
    "C-index heatmap for top 15 machine-learning models",
    "C-index heatmap for top 30 machine-learning models",
    "C-index heatmap for all machine-learning models",
    "Average C-index barplot for top 30 machine-learning models",
    "Risk score distribution across cohorts",
    "Survival status distribution across cohorts",
    "Kaplan-Meier curves by DPRS risk group",
    "Time-dependent ROC curves for DPRS",
    "Summary barplot of 1-, 3-, and 5-year AUC values",
    "Final DPRS signature gene coefficient plot"
  ),
  Path = c(
    file.path(FIG_DIR, "Fig4B_Cindex_heatmap_top15.pdf"),
    file.path(FIG_DIR, "Fig4B_Cindex_heatmap_top30.pdf"),
    file.path(FIG_DIR, "Fig4B_Cindex_heatmap_all_models.pdf"),
    file.path(FIG_DIR, "Fig4B_average_Cindex_barplot_top30.pdf"),
    file.path(FIG_DIR, "Fig4C_DPRS_risk_score_distribution.pdf"),
    file.path(FIG_DIR, "Fig4C_DPRS_survival_status_distribution.pdf"),
    file.path(FIG_DIR, "Fig4D_DPRS_KM_*.pdf"),
    file.path(FIG_DIR, "Fig4E_DPRS_timeROC_*.pdf"),
    file.path(FIG_DIR, "Fig4E_DPRS_timeROC_AUC_summary.pdf"),
    file.path(FIG_DIR, "Fig4F_DPRS_signature_gene_lollipop.pdf")
  ),
  stringsAsFactors = FALSE
)

save_csv(
  plot_manifest,
  file.path(OUT_DIR, "ML_DPRS_Figure4_plot_manifest.csv")
)

sink(file.path(OUT_DIR, "ML_DPRS_Figure4_plotting_session_info.txt"))
cat("ML-DPRS Figure 4 plotting finished.\n")
cat("Input directory:", ML_DIR, "\n")
cat("Output directory:", OUT_DIR, "\n")
cat("Selected model:", best_model, "\n\n")
cat("Generated figures:\n")
print(plot_manifest)
cat("\nSession info:\n")
print(sessionInfo())
sink()

cat("\nML-DPRS Figure 4 plotting finished.\n")
cat("Input directory:", ML_DIR, "\n")
cat("Output directory:", OUT_DIR, "\n")
cat("\nMain output figures:\n")
cat(" - figures/Fig4B_Cindex_heatmap_top15.pdf\n")
cat(" - figures/Fig4B_Cindex_heatmap_top30.pdf\n")
cat(" - figures/Fig4B_Cindex_heatmap_all_models.pdf\n")
cat(" - figures/Fig4B_average_Cindex_barplot_top30.pdf\n")
cat(" - figures/Fig4C_DPRS_risk_score_distribution.pdf\n")
cat(" - figures/Fig4C_DPRS_survival_status_distribution.pdf\n")
cat(" - figures/Fig4D_DPRS_KM_*.pdf\n")
cat(" - figures/Fig4E_DPRS_timeROC_*.pdf\n")
cat(" - figures/Fig4E_DPRS_timeROC_AUC_summary.pdf\n")
cat(" - figures/Fig4F_DPRS_signature_gene_lollipop.pdf\n")