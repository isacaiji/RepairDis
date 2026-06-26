############################################################
# GEO preprocessing for GSE72094 and GSE68465
############################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("Biobase", quietly = TRUE)) BiocManager::install("Biobase")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(GEOquery)
library(Biobase)
library(data.table)
library(dplyr)

PROJECT_DIR <- "D:/R_workspace/评分/AD_DDR_project"
DATA_DIR <- file.path(PROJECT_DIR, "00_data")

clean_gene_symbol_geo <- function(x) {
  x <- as.character(x)
  x <- gsub(" ///.*$", "", x)
  x <- gsub(" //.*$", "", x)
  x <- gsub(";.*$", "", x)
  x <- gsub(",.*$", "", x)
  x <- gsub("\\s*\\([^\\)]*\\)$", "", x)
  x <- toupper(trimws(x))
  x[x %in% c("", "NA", "---", "NULL", "N/A", "NA_NA")] <- NA
  x
}

safe_num <- function(x) {
  x <- as.character(x)
  x <- gsub(",", "", x)
  x <- gsub("[^0-9.\\-]", "", x)
  suppressWarnings(as.numeric(x))
}

flatten_pheno_characteristics <- function(pheno) {
  pheno <- as.data.frame(pheno, check.names = FALSE)
  out <- data.frame(Sample = rownames(pheno), stringsAsFactors = FALSE)
  
  if ("geo_accession" %in% colnames(pheno)) {
    out$Sample <- as.character(pheno$geo_accession)
  }
  
  char_cols <- grep("^characteristics_ch1", colnames(pheno), value = TRUE)
  
  if (length(char_cols) > 0) {
    for (cc in char_cols) {
      vals <- as.character(pheno[[cc]])
      has_key <- grepl(":", vals)
      
      keys <- ifelse(has_key, sub(":.*$", "", vals), NA)
      values <- ifelse(has_key, sub("^[^:]+:\\s*", "", vals), vals)
      
      keys <- tolower(keys)
      keys <- gsub("[^A-Za-z0-9]+", "_", keys)
      keys <- gsub("_+$", "", keys)
      keys <- gsub("^_+", "", keys)
      
      for (k in unique(keys[!is.na(keys) & keys != ""])) {
        new_col <- k
        if (new_col %in% colnames(out)) {
          new_col <- paste0(k, "_", cc)
        }
        out[[new_col]] <- ifelse(keys == k, values, NA)
      }
    }
  }
  
  cbind(out, pheno)
}

guess_symbol_col <- function(feature) {
  cn <- colnames(feature)
  pats <- c(
    "^gene_symbol$",
    "gene.?symbol",
    "^symbol$",
    "symbol",
    "gene_assignment",
    "gene.?name",
    "Gene Symbol",
    "GENE_SYMBOL"
  )
  
  for (p in pats) {
    hit <- grep(p, cn, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  
  NA_character_
}

get_platform_annotation <- function(eset) {
  fd <- Biobase::fData(eset)
  if (!is.null(fd) && nrow(fd) > 0 && ncol(fd) > 0) {
    symbol_col <- guess_symbol_col(fd)
    if (!is.na(symbol_col)) return(fd)
  }
  
  platform_id <- annotation(eset)
  message("Feature annotation not found in series matrix. Trying to download GPL: ", platform_id)
  
  gpl <- tryCatch({
    GEOquery::getGEO(platform_id)
  }, error = function(e) NULL)
  
  if (is.null(gpl)) {
    warning("Cannot download GPL annotation: ", platform_id)
    return(fd)
  }
  
  gpl_tab <- GEOquery::Table(gpl)
  rownames(gpl_tab) <- gpl_tab$ID
  gpl_tab
}

aggregate_probe_to_gene_max_iqr <- function(expr, feature) {
  symbol_col <- guess_symbol_col(feature)
  
  if (is.na(symbol_col)) {
    stop("Cannot identify gene symbol column. Feature columns are: ",
         paste(colnames(feature), collapse = ", "))
  }
  
  symbols <- clean_gene_symbol_geo(feature[rownames(expr), symbol_col])
  
  # 如果 rownames 不能直接匹配，尝试用 Probe/ID 列匹配
  if (all(is.na(symbols))) {
    if ("ID" %in% colnames(feature)) {
      idx <- match(rownames(expr), feature$ID)
      symbols <- clean_gene_symbol_geo(feature[idx, symbol_col])
    }
  }
  
  valid <- !is.na(symbols) & symbols != ""
  expr2 <- expr[valid, , drop = FALSE]
  symbols2 <- symbols[valid]
  
  storage.mode(expr2) <- "numeric"
  
  q99 <- suppressWarnings(quantile(expr2, 0.99, na.rm = TRUE))
  if (is.finite(q99) && q99 > 100) {
    expr2 <- log2(expr2 + 1)
  }
  
  iqr_val <- apply(expr2, 1, IQR, na.rm = TRUE)
  
  keep_idx <- tapply(seq_along(symbols2), symbols2, function(ii) {
    ii[which.max(iqr_val[ii])]
  })
  
  keep_idx <- as.integer(keep_idx)
  gene_expr <- expr2[keep_idx, , drop = FALSE]
  rownames(gene_expr) <- names(keep_idx)
  gene_expr <- gene_expr[!duplicated(rownames(gene_expr)), , drop = FALSE]
  
  gene_expr
}

find_col <- function(df, patterns) {
  cn <- colnames(df)
  for (p in patterns) {
    hit <- grep(p, cn, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

parse_status <- function(x) {
  y <- tolower(as.character(x))
  y <- trimws(y)
  
  out <- rep(NA_real_, length(y))
  
  out[y %in% c("1", "dead", "deceased", "death", "died", "yes", "true", "event", "relapse", "progression")] <- 1
  out[y %in% c("0", "alive", "living", "censored", "no", "false", "none", "non-event", "non_event")] <- 0
  
  # 如果是 numeric 0/1
  yn <- suppressWarnings(as.numeric(y))
  idx <- is.na(out) & yn %in% c(0, 1)
  out[idx] <- yn[idx]
  
  out
}

make_geo_clinical_auto <- function(pheno_full, gse_id) {
  sample_col <- if ("geo_accession" %in% colnames(pheno_full)) {
    "geo_accession"
  } else if ("Sample" %in% colnames(pheno_full)) {
    "Sample"
  } else {
    colnames(pheno_full)[1]
  }
  
  time_col <- find_col(
    pheno_full,
    c(
      "os.*time",
      "overall.*survival.*time",
      "survival.*time",
      "follow.*up.*time",
      "followup.*time",
      "follow.*up.*day",
      "days.*follow",
      "days.*survival",
      "survival.*days",
      "time.*death",
      "months.*survival",
      "survival.*month",
      "os.*month",
      "os.*day"
    )
  )
  
  status_col <- find_col(
    pheno_full,
    c(
      "os.*event",
      "overall.*survival.*event",
      "survival.*event",
      "vital.*status",
      "death.*status",
      "dead",
      "deceased",
      "status"
    )
  )
  
  candidate_cols <- grep(
    "os|surv|death|dead|vital|status|follow|time|month|day|event",
    colnames(pheno_full),
    ignore.case = TRUE,
    value = TRUE
  )
  
  fwrite(
    data.frame(CandidateClinicalColumns = candidate_cols),
    file.path(DATA_DIR, paste0(gse_id, "_candidate_clinical_columns.csv"))
  )
  
  if (is.na(time_col) || is.na(status_col)) {
    warning("Cannot automatically infer clinical time/status for ", gse_id,
            ". Please inspect candidate clinical columns.")
    return(NULL)
  }
  
  time <- safe_num(pheno_full[[time_col]])
  
  # 如果列名里出现 month，就转成 days
  if (grepl("month", time_col, ignore.case = TRUE)) {
    time <- time * 30.44
  }
  
  status <- parse_status(pheno_full[[status_col]])
  
  clin <- data.frame(
    Sample = as.character(pheno_full[[sample_col]]),
    time = as.numeric(time),
    status = as.numeric(status),
    RawTimeCol = time_col,
    RawStatusCol = status_col,
    stringsAsFactors = FALSE
  )
  
  clin <- clin[!is.na(clin$Sample) & clin$Sample != "", ]
  clin <- clin[is.finite(clin$time) & clin$time > 0 & clin$status %in% c(0, 1), ]
  clin <- clin[!duplicated(clin$Sample), ]
  
  clin
}

process_gse_series <- function(gse_id) {
  message("\n==============================")
  message("Processing ", gse_id)
  message("==============================")
  
  series_file <- file.path(DATA_DIR, paste0(gse_id, "_series_matrix.txt.gz"))
  if (!file.exists(series_file)) {
    stop("Series matrix not found: ", series_file)
  }
  
  eset <- GEOquery::getGEO(filename = series_file)
  
  if (is.list(eset)) {
    eset <- eset[[1]]
  }
  
  expr <- Biobase::exprs(eset)
  pheno <- Biobase::pData(eset)
  feature <- get_platform_annotation(eset)
  
  message("Probe expression matrix: ", paste(dim(expr), collapse = " x "))
  message("Pheno table: ", paste(dim(pheno), collapse = " x "))
  message("Feature table: ", paste(dim(feature), collapse = " x "))
  
  # 保存原始 pheno 和 feature
  pheno_full <- flatten_pheno_characteristics(pheno)
  
  fwrite(
    pheno_full,
    file.path(DATA_DIR, paste0(gse_id, "_pheno_full.csv"))
  )
  
  feature_out <- data.frame(
    Probe = rownames(feature),
    feature,
    check.names = FALSE
  )
  
  fwrite(
    feature_out,
    file.path(DATA_DIR, paste0(gse_id, "_feature_annotation.csv"))
  )
  
  # probe -> gene
  gene_expr <- aggregate_probe_to_gene_max_iqr(expr, feature)
  
  gene_expr_df <- data.frame(
    Gene = rownames(gene_expr),
    gene_expr,
    check.names = FALSE
  )
  
  # 最终表达文件
  fwrite(
    gene_expr_df,
    file.path(DATA_DIR, paste0(gse_id, "_expression.csv"))
  )
  
  # 自动整理 clinical
  clin <- make_geo_clinical_auto(pheno_full, gse_id)
  
  if (!is.null(clin)) {
    fwrite(
      clin,
      file.path(DATA_DIR, paste0(gse_id, "_clinical.csv"))
    )
    
    message("Clinical auto parsed:")
    message("n = ", nrow(clin), "; events = ", sum(clin$status == 1, na.rm = TRUE))
    print(table(clin$status, useNA = "ifany"))
    message("time summary:")
    print(summary(clin$time))
    message("time col = ", unique(clin$RawTimeCol))
    message("status col = ", unique(clin$RawStatusCol))
  } else {
    message("Clinical auto parsing failed. Check candidate clinical columns file.")
  }
  
  message("Gene expression matrix: ", paste(dim(gene_expr), collapse = " x "))
  
  invisible(list(
    expr_gene = gene_expr,
    pheno_full = pheno_full,
    clinical = clin,
    feature = feature
  ))
}

res72094 <- process_gse_series("GSE72094")
res68465 <- process_gse_series("GSE68465")

############################################################
# Quick check
############################################################

cat("\nFinal generated files:\n")
print(file.exists(file.path(DATA_DIR, "GSE72094_expression.csv")))
print(file.exists(file.path(DATA_DIR, "GSE72094_clinical.csv")))
print(file.exists(file.path(DATA_DIR, "GSE68465_expression.csv")))
print(file.exists(file.path(DATA_DIR, "GSE68465_clinical.csv")))

if (file.exists(file.path(DATA_DIR, "GSE72094_clinical.csv"))) {
  cat("\nGSE72094 clinical:\n")
  x <- fread(file.path(DATA_DIR, "GSE72094_clinical.csv"), data.table = FALSE)
  print(dim(x))
  print(table(x$status))
  print(summary(x$time))
}

if (file.exists(file.path(DATA_DIR, "GSE68465_clinical.csv"))) {
  cat("\nGSE68465 clinical:\n")
  x <- fread(file.path(DATA_DIR, "GSE68465_clinical.csv"), data.table = FALSE)
  print(dim(x))
  print(table(x$status))
  print(summary(x$time))
}



























############################################################
# Fix GSE72094 / GSE68465 expression and clinical
############################################################

library(GEOquery)
library(Biobase)
library(data.table)
library(dplyr)

DATA_DIR <- "D:/R_workspace/评分/AD_DDR_project/00_data"

clean_gene_symbol_geo2 <- function(x) {
  x <- as.character(x)
  
  out <- sapply(x, function(xx) {
    if (is.na(xx) || xx == "") return(NA_character_)
    
    xx <- trimws(xx)
    
    # Affymetrix GPL96 sometimes: "GENE1 /// GENE2"
    if (grepl("///", xx, fixed = TRUE)) {
      xx <- strsplit(xx, "///", fixed = TRUE)[[1]][1]
    }
    
    # Sometimes: "NM_XXX // SYMBOL // description"
    if (grepl(" // ", xx, fixed = TRUE)) {
      parts <- strsplit(xx, " // ", fixed = TRUE)[[1]]
      parts <- trimws(parts)
      parts <- parts[parts != ""]
      if (length(parts) >= 2) {
        xx <- parts[2]
      } else if (length(parts) == 1) {
        xx <- parts[1]
      }
    }
    
    xx <- gsub(";.*$", "", xx)
    xx <- gsub(",.*$", "", xx)
    xx <- gsub("\\s*\\([^\\)]*\\)$", "", xx)
    xx <- toupper(trimws(xx))
    
    if (xx %in% c("", "NA", "---", "NULL", "N/A", "NAN")) return(NA_character_)
    
    xx
  })
  
  as.character(out)
}

guess_symbol_col <- function(feature) {
  cn <- colnames(feature)
  
  pats <- c(
    "^GeneSymbol$",
    "^Gene Symbol$",
    "^gene_symbol$",
    "gene.?symbol",
    "^symbol$",
    "symbol",
    "gene_assignment"
  )
  
  for (p in pats) {
    hit <- grep(p, cn, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  
  NA_character_
}

guess_probe_col <- function(feature) {
  cn <- colnames(feature)
  
  pats <- c("^ID$", "^Probe$", "^ID_REF$", "probe")
  
  for (p in pats) {
    hit <- grep(p, cn, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  
  NA_character_
}

make_gene_expression_fixed <- function(gse_id) {
  message("\n==============================")
  message("Fixing expression: ", gse_id)
  message("==============================")
  
  series_file <- file.path(DATA_DIR, paste0(gse_id, "_series_matrix.txt.gz"))
  if (!file.exists(series_file)) stop("Missing file: ", series_file)
  
  eset <- GEOquery::getGEO(filename = series_file)
  if (is.list(eset)) eset <- eset[[1]]
  
  expr <- Biobase::exprs(eset)
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  
  feature_file1 <- file.path(DATA_DIR, paste0(gse_id, "_feature_annotation_fixed.csv"))
  feature_file2 <- file.path(DATA_DIR, paste0(gse_id, "_feature_annotation.csv"))
  
  if (file.exists(feature_file1)) {
    feature <- fread(feature_file1, data.table = FALSE, check.names = FALSE)
  } else {
    feature <- fread(feature_file2, data.table = FALSE, check.names = FALSE)
  }
  
  sym_col <- guess_symbol_col(feature)
  probe_col <- guess_probe_col(feature)
  
  message("symbol col = ", sym_col)
  message("probe col = ", probe_col)
  
  if (is.na(sym_col) || is.na(probe_col)) {
    stop("Cannot identify probe or symbol column for ", gse_id)
  }
  
  map <- data.frame(
    Probe = as.character(feature[[probe_col]]),
    Gene = clean_gene_symbol_geo2(feature[[sym_col]]),
    stringsAsFactors = FALSE
  )
  
  map <- map[!is.na(map$Gene) & map$Gene != "" & !is.na(map$Probe) & map$Probe != "", ]
  map <- map[!duplicated(map$Probe), ]
  
  common_probe <- intersect(rownames(expr), map$Probe)
  
  message("expr probes = ", nrow(expr))
  message("mapped probes = ", nrow(map))
  message("common probes = ", length(common_probe))
  
  if (length(common_probe) < 1000) {
    debug <- data.frame(
      expr_probe_head = head(rownames(expr), 20),
      map_probe_head = head(map$Probe, 20),
      map_gene_head = head(map$Gene, 20)
    )
    fwrite(debug, file.path(DATA_DIR, paste0(gse_id, "_probe_mapping_debug.csv")))
    stop("Too few common probes for ", gse_id)
  }
  
  expr2 <- expr[common_probe, , drop = FALSE]
  map2 <- map[match(common_probe, map$Probe), , drop = FALSE]
  
  q99 <- suppressWarnings(quantile(expr2, 0.99, na.rm = TRUE))
  if (is.finite(q99) && q99 > 100) {
    expr2 <- log2(expr2 + 1)
  }
  
  iqr_val <- apply(expr2, 1, IQR, na.rm = TRUE)
  
  probe_rank <- data.frame(
    Probe = common_probe,
    Gene = map2$Gene,
    IQR = iqr_val,
    stringsAsFactors = FALSE
  )
  
  probe_rank <- probe_rank %>%
    dplyr::filter(!is.na(Gene), Gene != "") %>%
    dplyr::arrange(Gene, dplyr::desc(IQR)) %>%
    dplyr::distinct(Gene, .keep_all = TRUE)
  
  selected_probes <- probe_rank$Probe
  selected_genes <- probe_rank$Gene
  
  gene_expr <- expr2[selected_probes, , drop = FALSE]
  rownames(gene_expr) <- selected_genes
  
  gene_expr <- gene_expr[!duplicated(rownames(gene_expr)), , drop = FALSE]
  
  out <- data.frame(
    Gene = rownames(gene_expr),
    gene_expr,
    check.names = FALSE
  )
  
  fwrite(out, file.path(DATA_DIR, paste0(gse_id, "_expression.csv")))
  
  message(gse_id, " final gene expression dim: ", paste(dim(gene_expr), collapse = " x "))
  
  invisible(gene_expr)
}

expr72094 <- make_gene_expression_fixed("GSE72094")
expr68465 <- make_gene_expression_fixed("GSE68465")

############################################################
# Make GSE68465 clinical manually
############################################################

safe_num <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "--", "NA", "N/A", "NULL")] <- NA
  suppressWarnings(as.numeric(x))
}

make_gse68465_clinical <- function() {
  ph <- fread(
    file.path(DATA_DIR, "GSE68465_pheno_full.csv"),
    data.table = FALSE,
    check.names = FALSE
  )
  
  sample_col <- if ("geo_accession" %in% colnames(ph)) "geo_accession" else "Sample"
  
  time_month <- safe_num(ph$months_to_last_contact_or_death)
  status_raw <- tolower(trimws(as.character(ph$vital_status)))
  
  status <- ifelse(
    status_raw == "dead", 1,
    ifelse(status_raw == "alive", 0, NA_real_)
  )
  
  clin <- data.frame(
    Sample = as.character(ph[[sample_col]]),
    time = time_month * 30.44,
    status = status,
    RawTimeCol = "months_to_last_contact_or_death",
    RawStatusCol = "vital_status",
    stringsAsFactors = FALSE
  )
  
  clin <- clin[!is.na(clin$Sample) & clin$Sample != "", ]
  clin <- clin[is.finite(clin$time) & clin$time > 0 & clin$status %in% c(0, 1), ]
  clin <- clin[!duplicated(clin$Sample), ]
  
  fwrite(clin, file.path(DATA_DIR, "GSE68465_clinical.csv"))
  
  cat("\nGSE68465 clinical generated:\n")
  print(dim(clin))
  print(table(clin$status))
  print(summary(clin$time))
  
  invisible(clin)
}

clin68465 <- make_gse68465_clinical()

############################################################
# Final check
############################################################

cat("\nFinal check:\n")

x1 <- fread(file.path(DATA_DIR, "GSE72094_expression.csv"), data.table = FALSE)
x2 <- fread(file.path(DATA_DIR, "GSE68465_expression.csv"), data.table = FALSE)
c1 <- fread(file.path(DATA_DIR, "GSE72094_clinical.csv"), data.table = FALSE)
c2 <- fread(file.path(DATA_DIR, "GSE68465_clinical.csv"), data.table = FALSE)

cat("GSE72094 expression:", paste(dim(x1), collapse = " x "), "\n")
cat("GSE72094 clinical:", paste(dim(c1), collapse = " x "), "; events=", sum(c1$status == 1), "\n")
cat("Matched GSE72094 samples:", length(intersect(colnames(x1)[-1], c1$Sample)), "\n\n")

cat("GSE68465 expression:", paste(dim(x2), collapse = " x "), "\n")
cat("GSE68465 clinical:", paste(dim(c2), collapse = " x "), "; events=", sum(c2$status == 1), "\n")
cat("Matched GSE68465 samples:", length(intersect(colnames(x2)[-1], c2$Sample)), "\n")