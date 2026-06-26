############################################################
# RepairDis Figure 3 final panel collector
# Figure 3 panels are finalized DPRS outputs; this script
# records and refreshes the panel files without redrawing models.
############################################################

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

BASE_DIR <- "D:/R_workspace/summary/Repairium_pancancer_analysis_plan"
PLOT_DIR <- file.path(BASE_DIR, "03-res/plots/张岩-图")
FIG_DIR <- file.path(PLOT_DIR, "Figure3")
SRC_DIRS <- c(
  "D:/研究生/wen/图/fig3",
  file.path(PLOT_DIR, "Figure3")
)

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

expected <- c(
  "Figure3A_Cindex_heatmap_top30.pdf",
  "Figure3B_DPRS_signature_gene_lollipop.pdf",
  "Figure3C_DPRS_risk_score_distribution.pdf",
  "Figure3D_DPRS_survival_status_distribution.pdf",
  "Figure3E_KM_Training.pdf",
  "Figure3F_KM_Testing.pdf",
  "Figure3G_KM_GSE72094.pdf",
  "Figure3H_KM_GSE68465.pdf",
  "Figure3I_timeROC_Training.pdf",
  "Figure3J_timeROC_Testing.pdf",
  "Figure3K_timeROC_GSE72094.pdf",
  "Figure3L_timeROC_GSE68465.pdf",
  "Figure3M_published_signature_benchmark.pdf",
  "Figure3N_published_signature_timeROC_GSE68465.pdf",
  "Figure3O_DPRS_stratified_predicted_drug_sensitivity.pdf"
)

find_source <- function(fname) {
  for (d in SRC_DIRS) {
    p <- file.path(d, fname)
    if (file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = TRUE))
  }
  NA_character_
}

manifest <- data.frame(
  File = expected,
  Source = vapply(expected, find_source, character(1)),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(manifest))) {
  if (!is.na(manifest$Source[i]) && file.exists(manifest$Source[i])) {
    target_i <- file.path(FIG_DIR, manifest$File[i])
    if (normalizePath(manifest$Source[i], winslash = "/", mustWork = TRUE) !=
        normalizePath(target_i, winslash = "/", mustWork = FALSE)) {
      file.copy(manifest$Source[i], target_i, overwrite = TRUE)
    }
  }
}

manifest$Target <- file.path(FIG_DIR, manifest$File)
manifest$Exists <- file.exists(manifest$Target)
manifest$Note <- ifelse(
  manifest$Exists,
  "Final DPRS panel copied/confirmed. Regeneration should use the original DPRS modeling scripts.",
  "Missing; check D:/研究生/wen/图/fig3 or original DPRS output directories."
)

data.table::fwrite(manifest, file.path(FIG_DIR, "Figure3_panel_inventory.csv"))

cat("Figure 3 panel inventory written to:\n")
cat(file.path(FIG_DIR, "Figure3_panel_inventory.csv"), "\n")
print(manifest)
