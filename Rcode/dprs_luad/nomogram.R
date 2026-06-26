############################################################
# Nomogram / Calibration / TimeAUC / DCA / Multi-Cox for DPRS
# Fixed version
############################################################

options(stringsAsFactors = FALSE)

############################
# 0. Parameters
############################

# Keep the script ASCII-safe: "\u8bc4\u5206" is "评分".
# You can override these paths before running:
#   Sys.setenv(DPRS_PROJECT_DIR = "D:/R_workspace/评分/AD_DDR_project")
#   Sys.setenv(DPRS_OUT_DIR = "...")
resolve_project_dir <- function() {
  env_dir <- Sys.getenv("DPRS_PROJECT_DIR", unset = NA_character_)
  if (!is.na(env_dir) && dir.exists(env_dir)) {
    return(normalizePath(env_dir, winslash = "/", mustWork = TRUE))
  }
  
  direct_dir <- file.path("D:/R_workspace", "\u8bc4\u5206", "AD_DDR_project")
  direct_dir <- enc2native(direct_dir)
  if (dir.exists(direct_dir)) {
    return(normalizePath(direct_dir, winslash = "/", mustWork = TRUE))
  }
  
  candidates <- list.dirs("D:/R_workspace", recursive = TRUE, full.names = TRUE)
  candidates <- candidates[basename(candidates) == "AD_DDR_project"]
  if (length(candidates) > 0) {
    return(normalizePath(candidates[1], winslash = "/", mustWork = TRUE))
  }
  
  stop("Cannot find AD_DDR_project. Set DPRS_PROJECT_DIR to the project path.")
}

PROJECT_DIR <- resolve_project_dir()
FIG4_DIR <- file.path(PROJECT_DIR, "04-5_Mime1Matched_StepCoxFixed")

OUT_DIR <- Sys.getenv("DPRS_OUT_DIR", unset = NA_character_)
if (is.na(OUT_DIR) || !nzchar(OUT_DIR)) {
  OUT_DIR <- file.path(FIG4_DIR, "Clinical_Nomogram_DCA")
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DPRS_FILE <- file.path(FIG4_DIR, "Fig4C_DPRS_all_sets.csv")
TCGA_CLIN_PROCESSED_FILE <- file.path(
  PROJECT_DIR,
  "01_processed",
  "LUAD_clinical_processed.csv"
)

CAL_BOOT_B <- as.integer(Sys.getenv("CAL_BOOT_B", unset = "200"))
if (!is.finite(CAL_BOOT_B) || CAL_BOOT_B < 1) CAL_BOOT_B <- 200L

############################
# 1. Packages
############################

pkgs <- c(
  "data.table", "dplyr", "ggplot2", "survival",
  "rms", "timeROC", "dcurves", "survminer"
)

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(survival)
  library(rms)
  library(timeROC)
  library(dcurves)
  library(survminer)
})

############################
# 2. Helper functions
############################

save_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file)
}

first_existing_col <- function(df, candidates, required = TRUE) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    if (required) {
      stop("Missing required column. Tried: ", paste(candidates, collapse = ", "))
    }
    return(NA_character_)
  }
  hit[1]
}

stage_to_num <- function(x) {
  raw <- toupper(trimws(as.character(x)))
  raw[raw %in% c("", "NA", "N/A", "NONE", "NULL")] <- NA_character_
  
  x <- gsub("^PATHOLOGIC", "", raw)
  x <- gsub("^CLINICAL", "", x)
  x <- gsub("^STAGE", "", x)
  x <- gsub("[^A-Z0-9]", "", x)
  x[grepl("DISCREP|UNKNOWN|NOTREPORTED|NOTAVAILABLE", x)] <- NA_character_
  
  dplyr::case_when(
    grepl("^4", x) | grepl("^IV", x) ~ 4,
    grepl("^3", x) | grepl("^III", x) ~ 3,
    grepl("^2", x) | grepl("^II", x) ~ 2,
    grepl("^1", x) | grepl("^I", x) ~ 1,
    TRUE ~ NA_real_
  )
}

tn_to_num <- function(x, prefix = "T") {
  x <- toupper(trimws(as.character(x)))
  x <- gsub(paste0("^", prefix), "", x)
  x <- gsub("[^0-9]", "", x)
  
  dplyr::case_when(
    grepl("^0", x) ~ 0,
    grepl("^1", x) ~ 1,
    grepl("^2", x) ~ 2,
    grepl("^3", x) ~ 3,
    grepl("^4", x) ~ 4,
    TRUE ~ NA_real_
  )
}

get_score_col <- function(df) {
  cand <- c("DPRS", "RiskScore", "riskScore", "RS", "ML_DDRscore")
  first_existing_col(df, cand, required = TRUE)
}

cal_to_df <- function(cal, time_label) {
  # calibrate.cph objects keep useful columns in the underlying matrix.
  # as.data.frame(cal) dispatches to a method that returns only one column.
  x <- as.data.frame(unclass(cal))
  
  pred_col <- intersect(c("mean.predicted", "predy"), colnames(x))[1]
  obs_col <- intersect(
    c("KM.corrected", "calibrated.corrected", "KM", "calibrated.orig", "actual"),
    colnames(x)
  )[1]
  
  if (is.na(pred_col) || is.na(obs_col)) {
    stop(
      "Cannot identify calibration columns. Available columns: ",
      paste(colnames(x), collapse = ", ")
    )
  }
  
  out <- data.frame(
    Predicted = as.numeric(x[[pred_col]]),
    Observed = as.numeric(x[[obs_col]]),
    Time = time_label
  )
  
  out[is.finite(out$Predicted) & is.finite(out$Observed), , drop = FALSE]
}

cox_table <- function(fit, model_name) {
  s <- summary(fit)
  conf <- as.data.frame(s$conf.int)
  coef <- as.data.frame(s$coefficients)
  
  data.frame(
    Model = model_name,
    Variable = rownames(coef),
    HR = conf[["exp(coef)"]],
    Lower95CI = conf[["lower .95"]],
    Upper95CI = conf[["upper .95"]],
    Pvalue = coef[["Pr(>|z|)"]],
    row.names = NULL,
    check.names = FALSE
  )
}

surv_risk <- function(fit, newdata, time) {
  pmin(pmax(as.numeric(1 - rms::survest(fit, newdata = newdata, times = time)$surv), 0), 1)
}

run_calibrate <- function(...) {
  warning_messages <- character()
  fit <- withCallingHandlers(
    rms::calibrate(...),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  attr(fit, "captured_warnings") <- unique(warning_messages)
  fit
}

############################
# 3. Load and merge data
############################

if (!file.exists(DPRS_FILE)) stop("DPRS file not found: ", DPRS_FILE)
if (!file.exists(TCGA_CLIN_PROCESSED_FILE)) {
  stop("Clinical file not found: ", TCGA_CLIN_PROCESSED_FILE)
}

dprs <- data.table::fread(DPRS_FILE, data.table = FALSE, check.names = FALSE)
clin <- data.table::fread(TCGA_CLIN_PROCESSED_FILE, data.table = FALSE, check.names = FALSE)

score_col <- get_score_col(dprs)
id_col <- first_existing_col(dprs, c("ID", "Patient", "Sample"))
dataset_col <- first_existing_col(dprs, c("Dataset", "dataset", "Set"))
risk_group_col <- first_existing_col(dprs, c("RiskGroup", "risk_group", "Group"), required = FALSE)
clin_patient_col <- first_existing_col(clin, c("Patient", "ID", "sample", "Sample"))
stage_col <- first_existing_col(clin, c("stage", "Stage", "pathologic_stage"))
age_col <- first_existing_col(clin, c("age", "Age", "age_at_diagnosis"), required = FALSE)
gender_col <- first_existing_col(clin, c("gender", "Gender", "sex", "Sex"), required = FALSE)
t_col <- first_existing_col(clin, c("T_stage", "T", "t_stage"), required = FALSE)
n_col <- first_existing_col(clin, c("N_stage", "N", "n_stage"), required = FALSE)

dprs_tcga <- dprs %>%
  dplyr::filter(.data[[dataset_col]] %in% c("Training", "Testing")) %>%
  dplyr::mutate(
    Patient = as.character(.data[[id_col]]),
    Dataset = as.character(.data[[dataset_col]]),
    DPRS = as.numeric(.data[[score_col]]),
    RiskGroup = if (!is.na(risk_group_col)) as.character(.data[[risk_group_col]]) else NA_character_,
    time = as.numeric(.data[["time"]]),
    status = as.numeric(.data[["status"]])
  ) %>%
  dplyr::select(Patient, Dataset, time, status, DPRS, RiskGroup)

clin_use <- clin %>%
  dplyr::mutate(
    Patient = as.character(.data[[clin_patient_col]]),
    stage_num = stage_to_num(.data[[stage_col]])
  )

clin_use$age <- if (!is.na(age_col)) as.numeric(clin_use[[age_col]]) else NA_real_
clin_use$gender <- if (!is.na(gender_col)) as.character(clin_use[[gender_col]]) else NA_character_
clin_use$T_num <- if (!is.na(t_col)) tn_to_num(clin_use[[t_col]], "T") else NA_real_
clin_use$N_num <- if (!is.na(n_col)) tn_to_num(clin_use[[n_col]], "N") else NA_real_

dat <- dprs_tcga %>%
  dplyr::left_join(
    clin_use %>% dplyr::select(Patient, age, gender, stage_num, T_num, N_num),
    by = "Patient"
  ) %>%
  dplyr::filter(
    is.finite(time),
    time > 0,
    status %in% c(0, 1),
    is.finite(DPRS),
    is.finite(stage_num)
  )

dat$Stage <- factor(
  dat$stage_num,
  levels = c(1, 2, 3, 4),
  labels = c("I", "II", "III", "IV")
)
dat$RiskGroup <- factor(dat$RiskGroup, levels = c("Low", "High"))

cat("Project directory:", PROJECT_DIR, "\n")
cat("Output directory:", OUT_DIR, "\n")
cat("Total input samples:", nrow(dat), "\n")
cat("Total events:", sum(dat$status == 1), "\n")
print(table(dat$Stage, useNA = "ifany"))

save_csv(dat, file.path(OUT_DIR, "Master_Analysis_Data.csv"))

############################
# 4. Fit Cox models
############################

fit_combined <- survival::coxph(
  survival::Surv(time, status) ~ DPRS + Stage,
  data = dat, x = TRUE, y = TRUE
)
fit_dprs <- survival::coxph(
  survival::Surv(time, status) ~ DPRS,
  data = dat, x = TRUE, y = TRUE
)
fit_stage <- survival::coxph(
  survival::Surv(time, status) ~ Stage,
  data = dat, x = TRUE, y = TRUE
)

dat$Combined_LP <- as.numeric(predict(fit_combined, type = "lp"))
dat$DPRS_LP <- as.numeric(predict(fit_dprs, type = "lp"))
dat$Stage_LP <- as.numeric(predict(fit_stage, type = "lp"))

############################
# 5. Figure C: Nomogram
############################

dat_rms <- dat %>% dplyr::select(time, status, DPRS, Stage)

dd <- rms::datadist(dat_rms)
options(datadist = "dd")

fit_nom <- rms::cph(
  survival::Surv(time, status) ~ DPRS + Stage,
  data = dat_rms,
  x = TRUE,
  y = TRUE,
  surv = TRUE,
  time.inc = 1095
)

surv_fun <- rms::Survival(fit_nom)

nom <- rms::nomogram(
  fit_nom,
  fun = list(
    function(x) surv_fun(365, x),
    function(x) surv_fun(1095, x),
    function(x) surv_fun(1825, x)
  ),
  funlabel = c("Pr(OS, 1-year)", "Pr(OS, 3-year)", "Pr(OS, 5-year)"),
  lp = FALSE
)

pdf(file.path(OUT_DIR, "Figure3C_nomogram.pdf"), width = 9, height = 6)
plot(nom, xfrac = 0.35)
dev.off()

############################
# 6. Figure D: Calibration curves
############################

set.seed(20260513)
m_use <- max(30, floor(nrow(dat_rms) / 5))

cal1 <- run_calibrate(
  fit_nom, cmethod = "KM", method = "boot", u = 365, m = m_use, B = CAL_BOOT_B
)
cal3 <- run_calibrate(
  fit_nom, cmethod = "KM", method = "boot", u = 1095, m = m_use, B = CAL_BOOT_B
)
cal5 <- run_calibrate(
  fit_nom, cmethod = "KM", method = "boot", u = 1825, m = m_use, B = CAL_BOOT_B
)

cal_warnings <- unique(c(
  attr(cal1, "captured_warnings"),
  attr(cal3, "captured_warnings"),
  attr(cal5, "captured_warnings")
))
if (length(cal_warnings) > 0) {
  save_csv(
    data.frame(Warning = cal_warnings),
    file.path(OUT_DIR, "Figure3D_calibration_warnings.csv")
  )
  cat("Calibration warnings captured:", length(cal_warnings), "\n")
}

cal_df <- dplyr::bind_rows(
  cal_to_df(cal1, "1-year"),
  cal_to_df(cal3, "3-year"),
  cal_to_df(cal5, "5-year")
)

save_csv(cal_df, file.path(OUT_DIR, "Figure3D_calibration_data.csv"))

p_cal <- ggplot(cal_df, aes(Predicted, Observed, color = Time)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_point(size = 1.8, alpha = 0.8) +
  geom_line(linewidth = 0.9) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "Predicted Survival Probability",
    y = "Observed Survival Probability",
    title = "Calibration Curves for Nomogram Predicted Survival"
  )

ggsave(file.path(OUT_DIR, "Figure3D_calibration_curves.pdf"), p_cal, width = 6, height = 5)

############################
# 7. Figure E: Time-dependent AUC from 1 to 5 years
############################

times_auc <- seq(365, 1825, by = 365)

marker_list <- list(
  Combined_model = dat$Combined_LP,
  DPRS = dat$DPRS_LP,
  Stage = dat$Stage_LP
)
if (sum(is.finite(dat$T_num)) >= 30) marker_list$T_stage <- dat$T_num
if (sum(is.finite(dat$N_num)) >= 30) marker_list$N_stage <- dat$N_num

auc_df <- dplyr::bind_rows(lapply(names(marker_list), function(nm) {
  marker <- marker_list[[nm]]
  keep <- is.finite(dat$time) &
    dat$time > 0 &
    dat$status %in% c(0, 1) &
    is.finite(marker)
  
  roc <- timeROC::timeROC(
    T = dat$time[keep],
    delta = dat$status[keep],
    marker = marker[keep],
    cause = 1,
    weighting = "marginal",
    times = times_auc,
    ROC = TRUE
  )
  
  data.frame(
    Variable = nm,
    Year = seq_along(times_auc),
    Time = times_auc,
    AUC = as.numeric(roc$AUC)
  )
}))

save_csv(auc_df, file.path(OUT_DIR, "Figure3E_time_dependent_AUC_1to5.csv"))

p_auc <- ggplot(auc_df, aes(Year, AUC, color = Variable)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.8) +
  scale_x_continuous(breaks = seq_along(times_auc)) +
  coord_cartesian(ylim = c(0.45, 0.95)) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(x = "Time, years", y = "AUC", title = "Time-dependent AUC curves")

ggsave(file.path(OUT_DIR, "Figure3E_time_dependent_AUC_curves.pdf"), p_auc, width = 6, height = 5)

############################
# 8. Figure F: DCA at 3 years
############################

fit_dprs_cph <- rms::cph(
  survival::Surv(time, status) ~ DPRS,
  data = dat_rms, x = TRUE, y = TRUE, surv = TRUE
)
fit_stage_cph <- rms::cph(
  survival::Surv(time, status) ~ Stage,
  data = dat_rms, x = TRUE, y = TRUE, surv = TRUE
)

dat_dca <- dat_rms
dat_dca$Risk_Combined <- surv_risk(fit_nom, dat_dca, 1095)
dat_dca$Risk_DPRS <- surv_risk(fit_dprs_cph, dat_dca, 1095)
dat_dca$Risk_Stage <- surv_risk(fit_stage_cph, dat_dca, 1095)

dca_res <- dcurves::dca(
  survival::Surv(time, status) ~ Risk_Combined + Risk_DPRS + Risk_Stage,
  data = dat_dca,
  time = 1095,
  thresholds = seq(0.01, 0.80, by = 0.01)
)

p_dca <- plot(dca_res, smooth = TRUE) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "Decision curve analysis at 3 years",
    x = "High Risk Threshold",
    y = "Net Benefit"
  )

ggsave(file.path(OUT_DIR, "Figure3F_DCA_3year.pdf"), p_dca, width = 6, height = 5)

############################
# 9. Univariable and multivariable Cox regression
############################

cox_dat <- dat %>%
  dplyr::filter(!is.na(age), !is.na(gender), !is.na(Stage), !is.na(RiskGroup)) %>%
  dplyr::mutate(
    gender = as.factor(tools::toTitleCase(tolower(gender)))
  )

fit_uni_risk <- survival::coxph(survival::Surv(time, status) ~ RiskGroup, data = cox_dat)
fit_uni_age <- survival::coxph(survival::Surv(time, status) ~ age, data = cox_dat)
fit_uni_gender <- survival::coxph(survival::Surv(time, status) ~ gender, data = cox_dat)
fit_uni_stage <- survival::coxph(survival::Surv(time, status) ~ Stage, data = cox_dat)

fit_multi <- survival::coxph(
  survival::Surv(time, status) ~ RiskGroup + age + gender + Stage,
  data = cox_dat
)

fit_multi_dprs_cont <- survival::coxph(
  survival::Surv(time, status) ~ DPRS + age + gender + Stage,
  data = cox_dat
)

cox_uni_table <- dplyr::bind_rows(
  cox_table(fit_uni_risk, "Univariable: RiskGroup"),
  cox_table(fit_uni_age, "Univariable: age"),
  cox_table(fit_uni_gender, "Univariable: gender"),
  cox_table(fit_uni_stage, "Univariable: Stage")
)
cox_multi_table <- cox_table(fit_multi, "Multivariable: RiskGroup + clinical")
cox_multi_dprs_cont_table <- cox_table(
  fit_multi_dprs_cont,
  "Multivariable: continuous DPRS + clinical"
)

save_csv(cox_uni_table, file.path(OUT_DIR, "Figure3X_Univariate_Cox.csv"))
save_csv(cox_multi_table, file.path(OUT_DIR, "Figure3X_Multivariate_Cox_RiskGroup.csv"))
save_csv(
  cox_multi_dprs_cont_table,
  file.path(OUT_DIR, "Figure3X_Multivariate_Cox_DPRScontinuous.csv")
)

p_forest <- survminer::ggforest(
  fit_multi,
  data = cox_dat,
  main = "Multivariate Cox Regression Analysis",
  fontsize = 1.0,
  noDigits = 2
)

ggsave(file.path(OUT_DIR, "Figure3X_Multivariate_Forest.pdf"), plot = p_forest, width = 8, height = 6)

p_forest_dprs_cont <- survminer::ggforest(
  fit_multi_dprs_cont,
  data = cox_dat,
  main = "Multivariate Cox Regression Analysis (Continuous DPRS)",
  fontsize = 1.0,
  noDigits = 2
)

ggsave(
  file.path(OUT_DIR, "Figure3X_Multivariate_Forest_DPRScontinuous.pdf"),
  plot = p_forest_dprs_cont,
  width = 8,
  height = 6
)

############################
# 10. Save model summary
############################

summary_file <- file.path(OUT_DIR, "Cox_model_summary.txt")
sink(summary_file)
on.exit({
  while (sink.number() > 0) sink()
}, add = TRUE)

cat("Combined model:\n")
print(summary(fit_combined))
cat("\nDPRS model:\n")
print(summary(fit_dprs))
cat("\nStage model:\n")
print(summary(fit_stage))
cat("\nUnivariable Cox table:\n")
print(cox_uni_table)
cat("\nMultivariate model (RiskGroup + clinical):\n")
print(summary(fit_multi))
cat("\nMultivariate model (continuous DPRS + clinical):\n")
print(summary(fit_multi_dprs_cont))
sink()

cat("\nDone.\n")
cat("Outputs saved to:\n", OUT_DIR, "\n")
