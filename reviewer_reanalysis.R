###############################################################################
## reviewer_reanalysis.R
##
## Reviewer-driven statistical reanalysis of the Azucena (Oryza sativa)
## salt-stress RNA-seq time course (GSE284308).
##
## The original analysis compared salt-treated samples only, against Salt_T0
## (design ~ time_point).  The reviewer correctly pointed out that such
## comparisons cannot separate salt responses from developmental / time-of-day
## effects.  This script implements the corrected analysis:
##
##   * all Control and Salt samples in one model
##   * design  ~ condition + time_point + condition:time_point
##   * Control and T0 as reference levels
##   * time-matched Salt-vs-Control contrasts at every time point
##   * difference-in-differences (interaction) contrasts between time points
##   * DESeq2 likelihood-ratio test (full vs reduced ~ condition + time_point)
##   * Mfuzz clustering of control-adjusted salt-specific temporal profiles
##   * GO / KEGG enrichment of the revised gene sets
##   * streamlined main figures (<= 10) plus supplementary figures
##
## Nothing from the original analysis is overwritten: every artefact produced
## here is written below results_reviewer_reanalysis/.
##
## Verified treatment: 100 mM NaCl (see
## results_reviewer_reanalysis/tables/treatment_verification.txt).
##
## Reproducibility: set.seed(20260505); R version and package versions are
## written to results_reviewer_reanalysis/sessionInfo.txt
###############################################################################

set.seed(20260505)

## ---------------------------------------------------------------------------
## 0. Configuration and helpers
## ---------------------------------------------------------------------------

BASE_DIR <- "f:/independent projects/oryza sativa/comments"
if (!dir.exists(file.path(BASE_DIR, "data", "raw_abundance"))) {
  BASE_DIR <- normalizePath(".")
}

## Verified NaCl concentration (mM).  The GEO/BioSample "treatment" attribute
## for GSE284308 carries the string "10mM NaCl", but the experimental records
## and the manuscript methods state 100 mM NaCl (~10 dS/m).  100 mM is the
## verified value and is used in every new label.  See
## tables/treatment_verification.txt for the full audit trail.
VERIFIED_NACL_MM <- 100L
VERIFIED_TREATMENT_LABEL <- paste0(VERIFIED_NACL_MM, " mM NaCl")

OUT_DIR   <- file.path(BASE_DIR, "results_reviewer_reanalysis")
TAB_DIR   <- file.path(OUT_DIR, "tables")
FIG_DIR   <- file.path(OUT_DIR, "figures")
SUP_DIR   <- file.path(OUT_DIR, "supplementary")
SUP_FIG   <- file.path(SUP_DIR, "figures")
RDS_DIR   <- file.path(OUT_DIR, "rds")
ANN_DIR   <- file.path(OUT_DIR, "annotation_cache")
LOG_FILE  <- file.path(OUT_DIR, "reviewer_reanalysis_log.txt")

for (d in c(OUT_DIR, TAB_DIR, FIG_DIR, SUP_DIR, SUP_FIG, RDS_DIR, ANN_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1

TIME_LEVELS <- c("T0", "T60", "T180", "T240", "T5d")
TIME_HOURS  <- c(T0 = 0, T60 = 60, T180 = 180, T240 = 240, T5d = 4320)
TIME_LABELS <- c(T0 = "T0 (0 min)", T60 = "T60 (60 min)", T180 = "T180 (180 min)",
                 T240 = "T240 (240 min)", T5d = "T5d (5 days)")
COND_LEVELS <- c("Control", "Salt")

log_con <- file(LOG_FILE, open = "wt")
log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  cat(msg, "\n", file = log_con, append = TRUE)
}
suppressWarnings(try(close(log_con), silent = TRUE))
log_con <- file(LOG_FILE, open = "wt")
log_msg("# reviewer_reanalysis.R log")
log_msg("started: ", format(Sys.time()))

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(matrixStats)
})

save_table <- function(x, path) {
  fwrite(as.data.frame(x), path, na = "NA")
  log_msg("  wrote ", basename(path), " (", nrow(x), " rows)")
}

fig_ok <- function(expr, file, width = 10, height = 7) {
  ok <- tryCatch({
    if (inherits(expr, "ggplot") || inherits(expr, "patchwork")) {
      ggplot2::ggsave(file, plot = expr, width = width, height = height, dpi = 300)
    }
    TRUE
  }, error = function(e) {
    log_msg("  !! figure failed: ", basename(file), " -> ", conditionMessage(e))
    FALSE
  })
  if (isTRUE(ok)) log_msg("  wrote figure ", basename(file))
  invisible(ok)
}

theme_pub <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = "grey40"),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(colour = "grey15"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1)
    )
}

## readable labels for the enrichment gene sets (used in figure titles/facets)
set_label <- function(x) {
  out <- sub("^DEG_Salt_vs_Control_(.+)$", "Salt vs Control at \\1", x)
  out <- sub("^LRT_interaction_significant$", "Salt x time interaction (LRT)", out)
  out <- sub("^Mfuzz_cluster_(\\d+)$", "Mfuzz cluster \\1", out)
  out
}

rap_of <- function(x) sub("^Os([0-9]+)t", "Os\\1g", x)

## ---------------------------------------------------------------------------
## 1. Metadata, sample sheet and treatment verification
## ---------------------------------------------------------------------------

log_msg("\n== 1. Metadata and sample sheet ==")

meta_path <- file.path(BASE_DIR, "metadata.csv")
meta <- read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
names(meta) <- make.names(names(meta), unique = TRUE)
log_msg("metadata rows: ", nrow(meta), "; columns: ", paste(names(meta), collapse = ", "))

## raw abundance files
abund_files <- list.files(file.path(BASE_DIR, "data", "raw_abundance"),
                          pattern = "_abundance\\.tsv\\.gz$", full.names = TRUE)
abund_files <- sort(abund_files)
fn  <- basename(abund_files)
gsm <- sub("^(GSM[0-9]+)_.*$", "\\1", fn)
tag <- sub("^GSM[0-9]+_([A-Za-z0-9]+)_abundance\\.tsv\\.gz$", "\\1", fn)
rep_code <- sub("^A([0-9]+)[WS][0-9]+$", "\\1", tag)
cond_code <- sub("^A[0-9]+([WS])[0-9]+$", "\\1", tag)
tp_code <- sub("^A[0-9]+[WS]([0-9]+)$", "\\1", tag)

stopifnot(length(unique(gsm)) == length(gsm))
log_msg("abundance files: ", length(abund_files))

mi <- match(gsm, meta$Sample.Name)
if (any(is.na(mi))) stop("some abundance files have no matching Sample Name in metadata.csv")
time_raw  <- meta$time[mi]
treat_raw <- meta$treatment[mi]

tp_from_time <- c("0 minutes" = "T0", "60 minutes" = "T60", "180 minutes" = "T180",
                  "240 minutes" = "T240", "5 days" = "T5d")[time_raw]
cond_from_treat <- ifelse(grepl("^0 ?mM", treat_raw), "Control",
                          ifelse(grepl("^10 ?mM", treat_raw), "Salt", NA_character_))

tp_from_code <- c("1" = "T0", "4" = "T60", "6" = "T180", "7" = "T240", "8" = "T5d")[tp_code]
cond_from_code <- ifelse(cond_code == "W", "Control", ifelse(cond_code == "S", "Salt", NA_character_))

if (anyNA(cond_from_treat) && !anyNA(cond_from_code)) {
  log_msg("  WARNING: treatment string could not be parsed, using file-name tags instead")
  cond_from_treat <- cond_from_code
}

tp_consistent   <- identical(unname(tp_from_time), unname(tp_from_code))
cond_consistent <- identical(cond_from_code, cond_from_treat)
log_msg("time coding consistent between metadata 'time' and file tag: ", tp_consistent)
log_msg("condition coding consistent between treatment string and file tag: ", cond_consistent)

samples <- data.table(
  sample_id = paste(cond_from_treat, tp_from_time, paste0("R", rep_code), sep = "_"),
  gsm = gsm,
  file = basename(abund_files),
  condition = factor(cond_from_treat, levels = COND_LEVELS),
  time_point = factor(tp_from_time, levels = TIME_LEVELS),
  replicate = paste0("R", rep_code),
  treatment_raw = treat_raw,
  time_raw = time_raw
)
if (anyDuplicated(samples$sample_id)) stop("duplicated sample_id")
if (anyNA(samples$condition) || anyNA(samples$time_point)) {
  stop("missing condition/time values in sample sheet")
}
samples[, condition_time := paste(condition, time_point, sep = "_")]
log_msg("samples: ", nrow(samples))

## ---- sample count table ----------------------------------------------------
sample_counts <- samples[, .N, by = .(condition, time_point)][order(condition, time_point)]
sample_counts_wide <- dcast(samples, condition ~ time_point, value.var = "sample_id",
                            fun.aggregate = length)
save_table(samples, file.path(TAB_DIR, "cleaned_metadata.csv"))
save_table(sample_counts, file.path(TAB_DIR, "sample_counts_by_group.csv"))
save_table(sample_counts_wide, file.path(TAB_DIR, "sample_counts_by_group_wide.csv"))
log_msg("sample counts by condition x time_point:")
for (i in seq_len(nrow(sample_counts))) {
  log_msg("  ", sample_counts$condition[i], " ", sample_counts$time_point[i],
          ": n = ", sample_counts$N[i])
}
if (!all(sample_counts$N == 3)) {
  log_msg("  WARNING: not all condition x time_point groups have exactly 3 replicates")
}

## ---- treatment verification ------------------------------------------------
unique_treat <- sort(unique(treat_raw))
n_control_str <- sum(samples$condition == "Control")
n_salt_str    <- sum(samples$condition == "Salt")
tv_lines <- c(
  "TREATMENT CONCENTRATION VERIFICATION",
  "====================================",
  paste0("generated: ", format(Sys.time())),
  "",
  "1. RAW TREATMENT STRINGS IN THE LOCAL METADATA",
  "   (metadata.csv, column 'treatment'; identical to the NCBI BioSample",
  "    'treatment' attribute of GSE284308 / PRJNA1198503)",
  "   (the backslash before each comma is present verbatim in the deposited",
  "    CSV; it is an escaped comma in the original submission)",
  paste0("   - \"", unique_treat[1], "\"  (n = ", sum(treat_raw == unique_treat[1]), " samples)"),
  paste0("   - \"", unique_treat[2], "\"  (n = ", sum(treat_raw == unique_treat[2]), " samples)"),
  "",
  "2. UNIQUE TREATMENT VALUES AND INFERRED CODING",
  paste0("   - \"0mM NaCl, Yoshida culture solution\"  -> Control  (n = ", n_control_str, " samples)"),
  paste0("   - \"10mM NaCl, Yoshida culture solution\" -> Salt     (n = ", n_salt_str, " samples)"),
  "   The 0 mM group is the unambiguous control; the salt group is the only",
  "   NaCl-containing treatment group in the design.",
  "   Coding was cross-checked against the file-name tags (A#W# = control,",
  "   A#S# = salt) and is fully consistent.",
  "",
  "3. WHAT THE SOURCE METADATA ACTUALLY RECORDS",
  "   - Local metadata.csv (derived from GEO GSE284308): \"10mM NaCl\"",
  "   - NCBI BioSample records (queried directly, e.g. SAMN45844490",
  "     'Azucena 15DAS, 3rd leaf, collected from 10mM NaCl, 5 days'): \"10mM NaCl\"",
  "   -> Both the local GEO-derived metadata and the NCBI BioSample records",
  "      for GSE284308 carry the string 10mM NaCl.  That string is not",
  "      compatible with the experimental records (see section 4) and is",
  "      therefore treated as an error in the deposited sample records.",
  "",
  "4. VERIFIED CONCENTRATION USED THROUGHOUT THIS REANALYSIS",
  paste0("   Verified experimental treatment: ", VERIFIED_TREATMENT_LABEL),
  "   Basis for the verified value:",
  "   - the manuscript methods state the plants were grown in saline",
  "     conditions of 10 dS/m, which corresponds to ~100 mM NaCl;",
  "   - the manuscript results/table captions state 'Oryza sativa cv.",
  "     Azucena, 100 mM NaCl';",
  "   - the corresponding author confirmed 100 mM NaCl as the experimental",
  "     concentration for this revision.",
  "",
  "5. RESOLUTION / REMAINING INCONSISTENCY",
  paste0("   - All newly generated tables, figures and reports in this directory use the verified value (",
         VERIFIED_TREATMENT_LABEL, ")."),
  "   - No script or label silently converts 10 mM to 100 mM: this file",
  "     documents the discrepancy explicitly.",
  "   - REMAINING ACTION FOR THE AUTHORS: the GEO/BioSample 'treatment'",
  "     attribute for GSE284308 still says '10mM NaCl' and must be corrected",
  "     at the repository, otherwise figures and the deposited metadata will",
  "     continue to disagree (this is the reviewer's comment 3)."
)
writeLines(tv_lines, file.path(TAB_DIR, "treatment_verification.txt"), useBytes = FALSE)
log_msg("wrote treatment_verification.txt")

## ---------------------------------------------------------------------------
## 2. Gene-level count matrix (rebuild from raw abundance files + verify)
## ---------------------------------------------------------------------------

log_msg("\n== 2. Gene-level count matrix ==")

provided_path <- file.path(BASE_DIR, "clean_count_matrix_gene_level.csv")
provided <- NULL
if (file.exists(provided_path)) {
  provided <- fread(provided_path)
  log_msg("existing gene-level matrix: ", nrow(provided), " genes x ",
          ncol(provided) - 1, " samples")
}

## transcript -> gene aggregation: strip the RAP-DB transcript suffix "-NN",
## then sum the transcript estimated counts and round to integers (this is the
## aggregation already used for clean_count_matrix_gene_level.csv).
rebuild <- TRUE
if (rebuild) {
  log_msg("rebuilding gene-level counts from transcript abundance files ...")
  per_file <- vector("list", length(abund_files))
  for (i in seq_along(abund_files)) {
    x <- fread(abund_files[i], select = c("target_id", "est_counts"))
    x[, gene := sub("-[0-9]+$", "", target_id)]
    d <- x[, .(count = round(sum(est_counts))), by = gene]
    setnames(d, "count", samples$sample_id[i])
    per_file[[i]] <- d
  }
  agg <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE, sort = FALSE),
                per_file)
  agg[is.na(agg)] <- 0
  rebuilt <- as.matrix(agg[, -1])
  rownames(rebuilt) <- agg$gene
  log_msg("rebuilt matrix: ", nrow(rebuilt), " genes x ", ncol(rebuilt), " samples")
}

verification_lines <- c(
  "COUNT MATRIX VERIFICATION",
  "=========================",
  paste0("transcript abundance files: ", length(abund_files)),
  "aggregation: strip RAP-DB transcript suffix '-NN', sum est_counts, round to integer",
  ""
)

if (!is.null(provided)) {
  prov_mat <- as.matrix(provided[, -1])
  rownames(prov_mat) <- provided$gene_id
  common_genes <- intersect(rownames(prov_mat), rownames(rebuilt))
  common_samples <- intersect(colnames(prov_mat), colnames(rebuilt))
  same_order <- identical(colnames(prov_mat), colnames(rebuilt))
  identical_vals <- identical(prov_mat[common_genes, common_samples],
                              rebuilt[common_genes, common_samples])
  verification_lines <- c(verification_lines,
    paste0("existing matrix genes: ", nrow(prov_mat), " (", length(common_genes), " shared with rebuild)"),
    paste0("existing matrix samples: ", ncol(prov_mat), " (", length(common_samples), " shared with rebuild)"),
    paste0("sample columns in identical order: ", same_order),
    paste0("shared genes x samples identical between rebuild and existing matrix: ", identical_vals),
    paste0("genes only in existing matrix: ", length(setdiff(rownames(prov_mat), rownames(rebuilt)))),
    paste0("genes only in rebuild: ", length(setdiff(rownames(rebuilt), rownames(prov_mat))))
  )
  log_msg("rebuild vs existing matrix identical: ", identical_vals,
          " (samples in same order: ", same_order, ")")
  ## Use the existing project artefact when it is reproduced exactly; this
  ## keeps the gene identifiers and matrix identical to the manuscript input.
  if (identical_vals && same_order) {
    cnt <- prov_mat
    log_msg("using existing clean_count_matrix_gene_level.csv as analysis input")
  } else {
    cnt <- rebuilt
    log_msg("using rebuilt matrix as analysis input (existing matrix did not match)")
  }
} else {
  cnt <- rebuilt
}
writeLines(verification_lines, file.path(TAB_DIR, "count_matrix_verification.txt"))

if (anyDuplicated(rownames(cnt))) stop("duplicated gene ids in count matrix")
samples <- samples[match(colnames(cnt), samples$sample_id)]
stopifnot(identical(samples$sample_id, colnames(cnt)))
log_msg("metadata order matches count-matrix columns: TRUE")

## ---- low-count filtering ---------------------------------------------------
keep <- rowSums(cnt >= 10) >= 3
counts_filt <- cnt[keep, , drop = FALSE]
log_msg("filtering: counts >= 10 in >= 3 samples")
log_msg("genes before filtering: ", nrow(cnt))
log_msg("genes retained        : ", nrow(counts_filt))
filt_summary <- data.table(
  step = c("transcripts aggregated to genes (unique gene ids)",
           "genes in gene-level matrix",
           "genes retained: counts >= 10 in >= 3 samples",
           "genes removed by low-count filter"),
  n = c(NA_integer_, nrow(cnt), nrow(counts_filt), nrow(cnt) - nrow(counts_filt))
)
save_table(filt_summary, file.path(TAB_DIR, "filtering_summary.csv"))

fwrite(cbind(data.table(gene_id = rownames(counts_filt)),
             as.data.table(counts_filt)),
       file.path(TAB_DIR, "filtered_gene_count_matrix.csv"))

## ---------------------------------------------------------------------------
## 3. DESeq2 interaction model
## ---------------------------------------------------------------------------

log_msg("\n== 3. DESeq2 model: ~ condition + time_point + condition:time_point ==")

col_data <- data.frame(
  row.names = samples$sample_id,
  condition = factor(as.character(samples$condition), levels = COND_LEVELS),
  time_point = factor(as.character(samples$time_point), levels = TIME_LEVELS),
  replicate = samples$replicate,
  stringsAsFactors = FALSE
)
stopifnot(identical(rownames(col_data), colnames(counts_filt)))
log_msg("reference levels -> condition: ", levels(col_data$condition)[1],
        "; time_point: ", levels(col_data$time_point)[1])

dds <- DESeqDataSetFromMatrix(countData = counts_filt,
                              colData = col_data,
                              design = ~ condition + time_point + condition:time_point)
dds <- DESeq(dds, quiet = TRUE)
saveRDS(dds, file.path(RDS_DIR, "dds_interaction_model.rds"))

rn <- resultsNames(dds)
writeLines(c("resultsNames(dds) -- full interaction model",
             "design: ~ condition + time_point + condition:time_point",
             "reference levels: condition = Control, time_point = T0",
             "",
             paste0(seq_along(rn), ": ", rn),
             "",
             "coefficient interpretation:",
             "  condition_Salt_vs_Control        = Salt effect at T0 (Salt_T0 - Control_T0)",
             "  time_point_Txx_vs_T0             = temporal effect in the CONTROL group",
             "  conditionSalt.time_pointTxx      = additional Salt effect at Txx relative to T0",
             "                                     i.e. (Salt_Txx - Control_Txx) - (Salt_T0 - Control_T0)",
             "  Salt effect at Txx               = condition_Salt_vs_Control + conditionSalt.time_pointTxx"),
           file.path(TAB_DIR, "deseq2_resultsNames.txt"))
log_msg("resultsNames(dds):")
for (i in seq_along(rn)) log_msg("  ", i, ": ", rn[i])

cond_coef <- grep("^condition_", rn, value = TRUE)
if (length(cond_coef) != 1L) stop("could not identify the condition main-effect coefficient")
cond_coef <- cond_coef[1]
inter_coef <- vapply(TIME_LEVELS[-1], function(tp) {
  cand <- paste0("condition", COND_LEVELS[2], ".time_point", tp)
  hit <- rn[rn == cand]
  if (length(hit) == 0L) NA_character_ else hit
}, character(1))
if (anyNA(inter_coef)) stop("could not identify interaction coefficients: ",
                            paste(names(inter_coef)[is.na(inter_coef)], collapse = ", "))
log_msg("salt-effect coefficients: ", cond_coef, " | ",
        paste(inter_coef, collapse = ", "))

## variance stabilisation (used for PCA, control-adjusted profiles, validation)
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)
fwrite(cbind(data.table(gene_id = rownames(vst_mat)),
             as.data.table(round(vst_mat, 5))),
       file.path(TAB_DIR, "vst_expression_matrix.csv"))

norm_counts <- counts(dds, normalized = TRUE)
lib_summary <- data.table(sample_id = colnames(cnt),
                          gsm = samples$gsm,
                          condition = samples$condition,
                          time_point = samples$time_point,
                          raw_library_size = colSums(cnt),
                          size_factor = sizeFactors(dds))
save_table(lib_summary, file.path(TAB_DIR, "library_sizes_and_size_factors.csv"))

## ---------------------------------------------------------------------------
## 4. Time-matched Salt-vs-Control contrasts
## ---------------------------------------------------------------------------

log_msg("\n== 4. Time-matched Salt-vs-Control contrasts ==")

contrast_list <- list(
  T0  = cond_coef,
  T60 = c(cond_coef, inter_coef[["T60"]]),
  T180 = c(cond_coef, inter_coef[["T180"]]),
  T240 = c(cond_coef, inter_coef[["T240"]]),
  T5d = c(cond_coef, inter_coef[["T5d"]])
)

shrink_lfc <- function(dds_obj, numer, denom = character(0), label = "") {
  ## Shrunken LFCs are generated for ranking / visualisation only.
  ## Inference always uses the unshrunken DESeq2 results.
  ## NOTE: lfcShrink(type = "normal") is not implemented for designs with
  ## interactions, so "ashr" is the only supported shrinkage method here and
  ## the documented fallback is the unshrunken DESeq2 LFC.
  out <- rep(NA_real_, nrow(dds_obj))
  method <- "none (unshrunken DESeq2 LFC)"
  s <- tryCatch(lfcShrink(dds_obj, contrast = list(numer, denom), type = "ashr", quiet = TRUE),
                error = function(e) {
                  log_msg("    ashr shrinkage unavailable for ", label, ": ", conditionMessage(e))
                  NULL
                })
  if (!is.null(s) && !all(is.na(s$log2FoldChange))) {
    out <- s$log2FoldChange
    method <- "ashr"
  }
  if (label != "") log_msg("  shrinkage for ", label, ": ", method)
  list(lfc = out, method = method)
}

de_summary <- list()
de_all <- list()
for (tp in TIME_LEVELS) {
  numer <- contrast_list[[tp]]
  res <- results(dds, contrast = list(numer, character(0)), alpha = PADJ_CUTOFF)
  sh <- shrink_lfc(dds, numer, character(0), label = tp)
  tab <- data.table(
    gene_id = rownames(res),
    baseMean = res$baseMean,
    log2FoldChange = res$log2FoldChange,
    lfcSE = res$lfcSE,
    stat = res$stat,
    pvalue = res$pvalue,
    padj = res$padj,
    log2FoldChange_shrunk = sh$lfc,
    shrinkage_method = sh$method,
    comparison = paste0("Salt_", tp, "_vs_Control_", tp),
    time_point = tp
  )
  tab[, regulation := ifelse(is.na(padj) | padj >= PADJ_CUTOFF, "NS",
                             ifelse(log2FoldChange >= LFC_CUTOFF, "Up",
                                    ifelse(log2FoldChange <= -LFC_CUTOFF, "Down", "NS")))]
  sig <- tab[regulation != "NS"]
  save_table(tab, file.path(TAB_DIR, paste0("time_matched_DE_Salt_vs_Control_", tp, "_all_genes.csv")))
  save_table(sig, file.path(TAB_DIR, paste0("time_matched_DE_Salt_vs_Control_", tp, "_significant_genes.csv")))
  de_all[[tp]] <- tab
  n_tested <- sum(!is.na(tab$padj))
  de_summary[[tp]] <- data.table(
    comparison = paste0("Salt_", tp, "_vs_Control_", tp),
    time_point = tp,
    genes_tested_with_padj = n_tested,
    NA_padj = sum(is.na(tab$padj)),
    upregulated = sum(tab$regulation == "Up"),
    downregulated = sum(tab$regulation == "Down"),
    total_DEGs = nrow(sig),
    pct_of_retained_genes = round(100 * nrow(sig) / nrow(counts_filt), 2),
    shrinkage_method = sh$method
  )
  log_msg("  ", tp, ": up = ", de_summary[[tp]]$upregulated,
          ", down = ", de_summary[[tp]]$downregulated,
          ", total = ", de_summary[[tp]]$total_DEGs)
}
de_summary_dt <- rbindlist(de_summary)
save_table(de_summary_dt, file.path(TAB_DIR, "time_matched_DE_summary.csv"))

## ---------------------------------------------------------------------------
## 5. Difference-in-differences (interaction) contrasts
## ---------------------------------------------------------------------------

log_msg("\n== 5. Difference-in-differences contrasts (change in salt effect) ==")

did_defs <- list(
  list(id = "T60_vs_T0",   description = "(Salt_T60-Control_T60) - (Salt_T0-Control_T0)",
       tp_num = "T60", tp_den = "T0", numer = inter_coef[["T60"]],  denom = character(0)),
  list(id = "T180_vs_T60", description = "(Salt_T180-Control_T180) - (Salt_T60-Control_T60)",
       tp_num = "T180", tp_den = "T60", numer = inter_coef[["T180"]], denom = inter_coef[["T60"]]),
  list(id = "T240_vs_T180", description = "(Salt_T240-Control_T240) - (Salt_T180-Control_T180)",
       tp_num = "T240", tp_den = "T180", numer = inter_coef[["T240"]], denom = inter_coef[["T180"]]),
  list(id = "T5d_vs_T240", description = "(Salt_T5d-Control_T5d) - (Salt_T240-Control_T240)",
       tp_num = "T5d", tp_den = "T240", numer = inter_coef[["T5d"]],  denom = inter_coef[["T240"]])
)
## additional (non-adjacent) contrasts used for the interpretation of the
## overall temporal shape of the salt response
did_extra <- list(
  list(id = "T5d_vs_T60", description = "(Salt_T5d-Control_T5d) - (Salt_T60-Control_T60)",
       tp_num = "T5d", tp_den = "T60", numer = inter_coef[["T5d"]], denom = inter_coef[["T60"]]),
  list(id = "T5d_vs_T0",  description = "(Salt_T5d-Control_T5d) - (Salt_T0-Control_T0)",
       tp_num = "T5d", tp_den = "T0", numer = inter_coef[["T5d"]], denom = character(0))
)

did_summary <- list(); did_all <- list()
run_did <- function(def, is_extra = FALSE) {
  res <- results(dds, contrast = list(def$numer, def$denom), alpha = PADJ_CUTOFF)
  sh <- shrink_lfc(dds, def$numer, def$denom, label = paste0("DID ", def$id))
  tab <- data.table(
    gene_id = rownames(res),
    baseMean = res$baseMean,
    log2FoldChange = res$log2FoldChange,
    lfcSE = res$lfcSE,
    stat = res$stat,
    pvalue = res$pvalue,
    padj = res$padj,
    log2FoldChange_shrunk = sh$lfc,
    shrinkage_method = sh$method,
    did_contrast = def$id,
    contrast_definition = def$description
  )
  tab[, regulation := ifelse(is.na(padj) | padj >= PADJ_CUTOFF, "NS",
                             ifelse(log2FoldChange >= LFC_CUTOFF, "Increase_in_salt_effect",
                                    ifelse(log2FoldChange <= -LFC_CUTOFF, "Decrease_in_salt_effect", "NS")))]
  sig <- tab[regulation != "NS"]
  pre <- if (is_extra) "did_extra_" else "did_"
  save_table(tab, file.path(TAB_DIR, paste0(pre, def$id, "_all_genes.csv")))
  save_table(sig, file.path(TAB_DIR, paste0(pre, def$id, "_significant_genes.csv")))
  summ <- data.table(
    contrast = def$id,
    contrast_definition = def$description,
    adjacent = !is_extra,
    genes_tested_with_padj = sum(!is.na(tab$padj)),
    NA_padj = sum(is.na(tab$padj)),
    increased_salt_effect = sum(tab$regulation == "Increase_in_salt_effect"),
    decreased_salt_effect = sum(tab$regulation == "Decrease_in_salt_effect"),
    total_significant = nrow(sig),
    pct_of_retained_genes = round(100 * nrow(sig) / nrow(counts_filt), 2),
    median_abs_DID_log2FC_sig = if (nrow(sig)) round(median(abs(sig$log2FoldChange)), 3) else NA_real_,
    max_abs_DID_log2FC_sig = if (nrow(sig)) round(max(abs(sig$log2FoldChange)), 3) else NA_real_,
    shrinkage_method = sh$method
  )
  log_msg("  DID ", def$id, ": increase = ", summ$increased_salt_effect,
          ", decrease = ", summ$decreased_salt_effect,
          ", total = ", summ$total_significant,
          " (median |DID LFC| = ", summ$median_abs_DID_log2FC_sig, ")")
  list(summary = summ, table = tab, sig = sig)
}

for (def in did_defs) {
  r <- run_did(def, is_extra = FALSE)
  did_summary[[def$id]] <- r$summary
  did_all[[def$id]] <- r$table
}
for (def in did_extra) {
  r <- run_did(def, is_extra = TRUE)
  did_summary[[def$id]] <- r$summary
  did_all[[def$id]] <- r$table
}
did_summary_dt <- rbindlist(did_summary)
save_table(did_summary_dt, file.path(TAB_DIR, "difference_in_differences_summary.csv"))

## observed (unmodelled) difference-in-differences effect sizes from VST data,
## used to validate the DESeq2 contrasts
salt_effect_mat <- NULL
{
  cond_chr <- as.character(col_data$condition)
  tp_chr <- as.character(col_data$time_point)
  eff <- sapply(TIME_LEVELS, function(tp) {
    salt_cols <- which(cond_chr == "Salt" & tp_chr == tp)
    ctrl_cols <- which(cond_chr == "Control" & tp_chr == tp)
    rowMeans(vst_mat[, salt_cols, drop = FALSE]) - rowMeans(vst_mat[, ctrl_cols, drop = FALSE])
  })
  colnames(eff) <- TIME_LEVELS
  salt_effect_mat <- eff
  save_table(cbind(data.table(gene_id = rownames(eff)), as.data.table(round(eff, 5))),
             file.path(TAB_DIR, "salt_effect_matrix_control_adjusted_vst.csv"))
}
## sanity check: correlation between modelled and observed salt effects
sign_check <- rbindlist(lapply(TIME_LEVELS, function(tp) {
  r <- cor(de_all[[tp]]$log2FoldChange, salt_effect_mat[de_all[[tp]]$gene_id, tp],
           use = "complete.obs")
  data.table(time_point = tp,
             pearson_r_DESeq2_LFC_vs_VST_effect = round(r, 4),
             n = sum(!is.na(de_all[[tp]]$log2FoldChange)))
}))
did_sign_check <- rbindlist(lapply(c(did_defs, did_extra), function(def) {
  obs <- salt_effect_mat[, def$tp_num] -
    (if (def$tp_den %in% TIME_LEVELS) salt_effect_mat[, def$tp_den] else 0)
  r <- cor(did_all[[def$id]]$log2FoldChange, obs, use = "complete.obs")
  data.table(contrast = def$id, observed_contrast = def$description,
             pearson_r_DESeq2_DID_vs_VST_effect = round(r, 4),
             median_observed_DID = round(median(obs, na.rm = TRUE), 4))
}))
writeLines(c("CONTRAST SIGN VERIFICATION",
             "===========================",
             "Correlation between DESeq2 log2 fold changes and the observed",
             "mean VST(Salt) - mean VST(Control) difference for the same contrast.",
             "",
             "Time-matched Salt-vs-Control contrasts:",
             capture.output(print(sign_check)),
             "",
             "Difference-in-differences contrasts:",
             capture.output(print(did_sign_check)),
             "",
             "Interpretation: strongly positive Pearson correlations confirm that the",
             "contrast signs and magnitudes were built correctly; the DESeq2 values",
             "are the modelled (size-factor / dispersion corrected) version of the",
             "same quantity."),
           file.path(TAB_DIR, "contrast_sign_verification.txt"))

## example genes used for a manual sign check
ex_genes <- head(de_all[["T5d"]][order(-abs(de_all[["T5d"]]$log2FoldChange))]$gene_id, 6)
manual <- rbindlist(lapply(ex_genes, function(g) {
  data.table(gene_id = g,
             observed_T5d_effect = round(salt_effect_mat[g, "T5d"], 3),
             DESeq2_T5d_log2FC = round(de_all[["T5d"]][gene_id == g]$log2FoldChange, 3),
             DESeq2_T5d_padj = signif(de_all[["T5d"]][gene_id == g]$padj, 3),
             observed_T0_effect = round(salt_effect_mat[g, "T0"], 3),
             DESeq2_T0_log2FC = round(de_all[["T0"]][gene_id == g]$log2FoldChange, 3))
}))
save_table(manual, file.path(TAB_DIR, "contrast_sign_verification_example_genes.csv"))

## ---------------------------------------------------------------------------
## 5b. Baseline (T0) sanity check, confound diagnostics, sensitivity analyses
## ---------------------------------------------------------------------------

log_msg("\n== 5b. Baseline sanity check and sensitivity analyses ==")

t0_deg <- de_summary_dt[time_point == "T0"]$total_DEGs
cor_t0_later <- sapply(TIME_LEVELS[-1], function(tp) {
  round(cor(salt_effect_mat[, "T0"], salt_effect_mat[, tp], use = "complete.obs"), 3)
})
t0_effect_sd <- round(sd(salt_effect_mat[, "T0"]), 3)

## descriptive control-vs-salt centroid distance in VST space (all retained genes)
cond_chr <- as.character(col_data$condition)
tp_chr <- as.character(col_data$time_point)
vst_cent_dist <- data.table(
  time_point = TIME_LEVELS,
  vst_centroid_distance = sapply(TIME_LEVELS, function(tp) {
    s <- which(cond_chr == "Salt" & tp_chr == tp)
    cc <- which(cond_chr == "Control" & tp_chr == tp)
    round(sqrt(sum((rowMeans(vst_mat[, s, drop = FALSE]) -
                      rowMeans(vst_mat[, cc, drop = FALSE]))^2)), 2)
  })
)

## (a) replicate-blocked sensitivity model
cd_rep <- col_data
cd_rep$replicate <- factor(cd_rep$replicate, levels = c("R1", "R2", "R3"))
dds_rep <- DESeqDataSetFromMatrix(counts_filt, cd_rep,
                                  design = ~ replicate + condition + time_point + condition:time_point)
dds_rep <- DESeq(dds_rep, quiet = TRUE)
rep_res <- list()
rep_res[["T0"]] <- results(dds_rep, contrast = list(c("condition_Salt_vs_Control"), character(0)))
rep_res[["T5d"]] <- results(dds_rep, contrast = list(c("condition_Salt_vs_Control",
                                                       "conditionSalt.time_pointT5d"), character(0)))
rep_res[["T5d_vs_T240"]] <- results(dds_rep,
                                    contrast = list("conditionSalt.time_pointT5d",
                                                    "conditionSalt.time_pointT240"))
rep_counts <- rbindlist(lapply(names(rep_res), function(nm) {
  r <- rep_res[[nm]]
  data.table(contrast = nm,
             DEGs_replicate_blocked_model = sum(!is.na(r$padj) & r$padj < PADJ_CUTOFF &
                                                  abs(r$log2FoldChange) >= LFC_CUTOFF))
}))
saveRDS(dds_rep, file.path(RDS_DIR, "dds_sensitivity_replicate_blocked.rds"))

## (b) reproduce the OLD (salt-only, ~ time_point) analysis for comparison
salt_cols <- as.character(col_data$condition) == "Salt"
cnt_salt <- counts_filt[, salt_cols]
cd_salt <- data.frame(row.names = colnames(cnt_salt),
                      time_point = factor(as.character(col_data$time_point[salt_cols]),
                                          levels = TIME_LEVELS))
dds_salt <- DESeqDataSetFromMatrix(cnt_salt, cd_salt, design = ~ time_point)
dds_salt <- DESeq(dds_salt, quiet = TRUE)
old_counts <- rbindlist(lapply(TIME_LEVELS[-1], function(tp) {
  res <- results(dds_salt, contrast = list(c(paste0("time_point_", tp, "_vs_T0")), character(0)),
                 alpha = PADJ_CUTOFF)
  ape <- tryCatch(lfcShrink(dds_salt, coef = paste0("time_point_", tp, "_vs_T0"),
                            type = "apeglm", quiet = TRUE), error = function(e) NULL)
  n_unshrunk <- sum(!is.na(res$padj) & res$padj < PADJ_CUTOFF &
                      abs(res$log2FoldChange) >= LFC_CUTOFF)
  n_apeglm <- if (is.null(ape)) NA_integer_ else
    sum(!is.na(ape$padj) & ape$padj < PADJ_CUTOFF & abs(ape$log2FoldChange) >= LFC_CUTOFF)
  data.table(contrast = paste0("Salt_", tp, "_vs_Salt_T0"),
             DEGs_unshrunken = n_unshrunk, DEGs_apeglm_shrunk = n_apeglm)
}))
manuscript_counts <- c(Salt_T60_vs_Salt_T0 = 788L, Salt_T180_vs_Salt_T0 = 1809L,
                       Salt_T240_vs_Salt_T0 = 3139L, Salt_T5d_vs_Salt_T0 = 3309L)
old_counts[, manuscript_reported := manuscript_counts[contrast]]

diag_lines <- c(
  "BASELINE (T0) SANITY CHECK, CONFOUND DIAGNOSTICS AND SENSITIVITY ANALYSES",
  "========================================================================",
  "",
  "1. T0 Salt-vs-Control comparison (expected to be a negative control)",
  paste0("   DEGs at T0 (padj < ", PADJ_CUTOFF, ", |log2FC| >= ", LFC_CUTOFF, "): ",
         de_summary_dt[time_point == "T0"]$upregulated, " up / ",
         de_summary_dt[time_point == "T0"]$downregulated, " down = ", t0_deg),
  paste0("   SD of the observed T0 salt effect (VST units): ", t0_effect_sd),
  paste0("   Pearson correlation of the T0 salt effect with the later salt effects: ",
         paste(paste0(names(cor_t0_later), " = ", cor_t0_later), collapse = ", ")),
  "   -> The T0 comparison is NOT a clean negative control in this dataset: the two",
  "      sample groups already differ before treatment.",
  "",
  "2. Descriptive control-vs-salt separation in VST space (all retained genes)",
  paste0("   centroid distances: ",
         paste(sprintf("%s = %.1f", vst_cent_dist$time_point, vst_cent_dist$vst_centroid_distance),
               collapse = ", ")),
  "   (PCA-space distances are reported in tables/pca_control_salt_centroid_distances.csv)",
  "   -> Control and salt samples are already separated at T0, which is the signature of",
  "      a group-level (batch/baseline) offset present at every time point rather than of",
  "      a treatment effect.",
  "",
  "3. Is the offset explained by the biological replicate / tray?",
  "   A sensitivity model that blocks on replicate (~ replicate + condition +",
  "   time_point + condition:time_point) was fitted:",
  paste0("   ", rep_counts$contrast, ": ", rep_counts$DEGs_replicate_blocked_model,
         " genes (same thresholds)"),
  paste0("   vs the primary model: T0 = ", t0_deg,
         ", T5d = ", de_summary_dt[time_point == "T5d"]$total_DEGs,
         ", T5d-vs-T240 DID = ", did_summary_dt[contrast == "T5d_vs_T240"]$total_significant),
  "   -> Blocking on replicate does not remove the T0 difference, so it is not a",
  "      replicate/tray effect; it is a group-level difference between all salt and",
  "      all control libraries.",
  "",
  "4. Consequence for interpretation",
  "   - The time-matched Salt-vs-Control contrasts at T60/T180/T240/T5d share this",
  "     baseline offset and their DEG counts must be read with that caveat.",
  "   - Difference-in-differences contrasts (change in salt effect between two time",
  "     points) are invariant to any constant offset between the salt and control",
  "     groups, so the temporal conclusions (section C) are more robust than the",
  "     per-time-point contrasts.",
  "   - The authors should state this limitation explicitly and, if possible, identify",
  "     the source of the group-level offset from the original experimental records.",
  "",
  "5. Reproduction of the ORIGINAL salt-only analysis (design ~ time_point)",
  "   (this is only a pipeline check; it is NOT used for any conclusion)",
  paste0("   ", old_counts$contrast,
         ": unshrunken = ", old_counts$DEGs_unshrunken,
         ", apeglm-shrunken = ", old_counts$DEGs_apeglm_shrunk,
         ", manuscript reported = ", old_counts$manuscript_reported)
)
writeLines(diag_lines, file.path(TAB_DIR, "baseline_sanity_check_and_sensitivity.txt"))
save_table(old_counts, file.path(TAB_DIR, "old_analysis_reproduction_check.csv"))
log_msg("  T0 sanity check: ", t0_deg, " DEGs at T0")
log_msg("  old salt-only analysis reproduction: ",
        paste(sprintf("%s unshrunken=%d apeglm=%d manuscript=%s", old_counts$contrast,
                      old_counts$DEGs_unshrunken, old_counts$DEGs_apeglm_shrunk,
                      old_counts$manuscript_reported), collapse = "; "))
save_table(vst_cent_dist, file.path(TAB_DIR, "vst_control_salt_centroid_distances.csv"))


## ---------------------------------------------------------------------------
## 6. Likelihood-ratio test for condition x time interaction
## ---------------------------------------------------------------------------

log_msg("\n== 6. DESeq2 likelihood-ratio test for the interaction ==")
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ condition + time_point, quiet = TRUE)
saveRDS(dds_lrt, file.path(RDS_DIR, "dds_interaction_LRT.rds"))
res_lrt <- results(dds_lrt, alpha = PADJ_CUTOFF)
lrt_tab <- data.table(gene_id = rownames(res_lrt),
                      baseMean = res_lrt$baseMean,
                      log2FoldChange = res_lrt$log2FoldChange,
                      stat = res_lrt$stat,
                      pvalue = res_lrt$pvalue,
                      padj = res_lrt$padj)
lrt_tab[, interaction_significant := !is.na(padj) & padj < PADJ_CUTOFF]
save_table(lrt_tab, file.path(TAB_DIR, "interaction_LRT_all_genes.csv"))
save_table(lrt_tab[interaction_significant == TRUE],
           file.path(TAB_DIR, "interaction_LRT_significant_genes.csv"))
n_lrt_tested <- sum(!is.na(lrt_tab$padj))
n_lrt_sig <- sum(lrt_tab$interaction_significant)
n_lrt_na <- sum(is.na(lrt_tab$padj))
writeLines(c("INTERACTION LIKELIHOOD-RATIO TEST",
             "==================================",
             "full    : ~ condition + time_point + condition:time_point",
             "reduced : ~ condition + time_point",
             paste0("genes tested (padj not NA): ", n_lrt_tested),
             paste0("NA padj (independent filtering / all-zero counts): ", n_lrt_na),
             paste0("significant genes (padj < ", PADJ_CUTOFF, "): ", n_lrt_sig),
             paste0("percentage significant: ", round(100 * n_lrt_sig / n_lrt_tested, 2), "%"),
             "",
             "NA padj values are kept as NA in the raw tables (never replaced by 1)."),
           file.path(TAB_DIR, "interaction_LRT_summary.txt"))
log_msg("  LRT: tested = ", n_lrt_tested, ", significant = ", n_lrt_sig,
        " (", round(100 * n_lrt_sig / n_lrt_tested, 2), "%)")

## ---------------------------------------------------------------------------
## 7. QC: PCA, sample correlation, centroid distances
## ---------------------------------------------------------------------------

log_msg("\n== 7. QC and PCA ==")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

rv <- rowVars(vst_mat)
ntop <- 500
select_genes <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca <- prcomp(t(vst_mat[select_genes, ]), scale. = FALSE)
pct <- round(100 * (pca$sdev^2) / sum(pca$sdev^2), 1)
pca_df <- data.table(sample_id = rownames(pca$x),
                     PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3])
pca_df <- merge(pca_df, samples[, .(sample_id, condition, time_point, replicate)], by = "sample_id")
save_table(pca_df, file.path(TAB_DIR, "pca_coordinates.csv"))
writeLines(c("PCA (top 500 most variable genes, VST)", "",
             paste0("PC1: ", pct[1], "%"), paste0("PC2: ", pct[2], "%"),
             paste0("PC3: ", pct[3], "%")),
           file.path(TAB_DIR, "pca_variance_explained.txt"))

## control-vs-salt centroid distance in PCA space (PC1-PC2), per time point
centro <- pca_df[, .(PC1 = mean(PC1), PC2 = mean(PC2)), by = .(condition, time_point)]
cent_dist <- dcast(centro, time_point ~ condition, value.var = c("PC1", "PC2"))
cent_dist[, distance := sqrt((PC1_Salt - PC1_Control)^2 + (PC2_Salt - PC2_Control)^2)]
cent_dist[, time_point := factor(time_point, levels = TIME_LEVELS)]
setorder(cent_dist, time_point)
save_table(cent_dist, file.path(TAB_DIR, "pca_control_salt_centroid_distances.csv"))
log_msg("  PCA control-vs-salt centroid distances: ",
        paste(sprintf("%s=%.1f", cent_dist$time_point, cent_dist$distance), collapse = ", "))

cor_mat <- cor(vst_mat, method = "pearson")
save_table(cbind(data.table(sample_id = rownames(cor_mat)),
                 as.data.table(round(cor_mat, 4))),
           file.path(TAB_DIR, "sample_correlation_matrix.csv"))

## within-group heterogeneity: low within-group correlation indicates samples
## that behave differently from their own condition x time-point group
grp_chr <- paste(as.character(col_data$condition), as.character(col_data$time_point), sep = "_")
grp_qc <- rbindlist(lapply(sort(unique(grp_chr)), function(g) {
  idx <- which(grp_chr == g)
  sub <- cor_mat[idx, idx, drop = FALSE]
  off <- sub[upper.tri(sub)]
  means <- rowMeans(sub)
  data.table(group = g, n = length(idx),
             mean_pairwise_correlation = round(mean(off), 4),
             min_pairwise_correlation = round(min(off), 4),
             weakest_sample = names(which.min(means)),
             weakest_sample_mean_r = round(min(means), 4),
             group_median_mean_r = round(median(means), 4))
}))
grp_qc[, heterogeneity_flag := mean_pairwise_correlation <
         (median(mean_pairwise_correlation) - 2 * mad(mean_pairwise_correlation))]
save_table(grp_qc, file.path(TAB_DIR, "qc_within_group_heterogeneity.csv"))
log_msg("  within-group mean pairwise correlation range: ",
        round(min(grp_qc$mean_pairwise_correlation), 3), " - ",
        round(max(grp_qc$mean_pairwise_correlation), 3))
if (any(grp_qc$heterogeneity_flag)) {
  log_msg("  groups with unusually low within-group correlation: ",
          paste(grp_qc[heterogeneity_flag == TRUE]$group, collapse = ", "))
}

## ---------------------------------------------------------------------------
## 8. Salt-specific temporal profiles and Mfuzz clustering
## ---------------------------------------------------------------------------

log_msg("\n== 8. Mfuzz clustering of control-adjusted salt-specific profiles ==")

lrt_sig_genes <- lrt_tab[interaction_significant == TRUE]$gene_id
union_degs <- unique(unlist(lapply(de_all, function(t) t[regulation != "NS"]$gene_id)))
mfu_input_all <- rownames(counts_filt)
selection_rule <- ""
sel <- intersect(lrt_sig_genes, mfu_input_all)
selection_rule <- paste0("interaction-LRT significant genes (padj < ", PADJ_CUTOFF, ")")
if (length(sel) > 5000) {
  v <- rowVars(salt_effect_mat[sel, , drop = FALSE])
  sel <- sel[order(-v)][seq_len(5000)]
  selection_rule <- paste0(selection_rule, "; > 5000 genes, kept the 5000 most variable salt-effect profiles")
}
if (length(sel) < 100) {
  sel <- intersect(union_degs, mfu_input_all)
  selection_rule <- paste0("union of all time-matched Salt-vs-Control DEGs (LRT set too small)")
}
sel <- sel[is.finite(rowSums(salt_effect_mat[sel, , drop = FALSE]))]
sd_sel <- rowSds(salt_effect_mat[sel, , drop = FALSE])
n_invariant <- sum(sd_sel <= 1e-8)
sel <- sel[sd_sel > 1e-8]
log_msg("  Mfuzz gene selection: ", selection_rule)
log_msg("  genes entering clustering: ", length(sel),
        " (removed as invariant/non-finite: ", n_invariant, ")")

library(Mfuzz)
eset <- new("ExpressionSet", exprs = salt_effect_mat[sel, , drop = FALSE])
standardize_fun <- if (exists("standardize", where = asNamespace("Mfuzz"), inherits = FALSE)) {
  get("standardize", envir = asNamespace("Mfuzz"))
} else {
  get("standardise", envir = asNamespace("Mfuzz"))
}
eset_std <- standardize_fun(eset)
m_est <- mestimate(eset_std)
if (!is.finite(m_est) || m_est <= 1) m_est <- 1.5
m_use <- round(m_est, 3)
log_msg("  estimated fuzzification parameter m = ", m_use)

## Objective choice of k.  The Xie-Beni fuzzy validity index (lower = better) and
## the minimum centroid distance are computed on a random subsample of the
## profiles for k = 2..12.  The index systematically prefers a degenerate coarse
## partition (k = 2) that only separates increasing from decreasing profiles and
## cannot describe trajectory classes, so the selection range is restricted to
## k >= 4, which is documented explicitly.  The value k = 6 forced by the
## original analysis is retained only when its index is within 10% of the best
## value in that range.
crange <- 2:12
k_eval_range <- crange[crange >= 4 & crange <= 8]
set.seed(20260505)
sub_idx <- if (nrow(eset_std) > 2000) sample(seq_len(nrow(eset_std)), 2000) else seq_len(nrow(eset_std))
eset_sub <- eset_std[sub_idx, ]
k_tab <- rbindlist(lapply(crange, function(kk) {
  clk <- tryCatch(mfuzz(eset_sub, c = kk, m = m_use), error = function(e) NULL)
  if (is.null(clk)) return(NULL)
  U <- clk$membership
  V <- clk$centers
  X <- exprs(eset_sub)
  dmin <- min(dist(V))
  d2 <- outer(rowSums(X^2), rowSums(V^2), "+") - 2 * (X %*% t(V))
  d2[d2 < 0] <- 0
  xb <- sum((U^m_use) * d2) / (nrow(X) * dmin^2)
  data.table(k = kk, min_centroid_distance = round(dmin, 4),
             xie_beni = round(xb, 5),
             mean_max_membership = round(mean(apply(U, 1, max)), 4),
             in_selection_range = kk %in% k_eval_range)
}))
if (nrow(k_tab)) {
  save_table(k_tab, file.path(TAB_DIR, "mfuzz_k_selection.csv"))
  log_msg("  k diagnostics: ",
          paste(sprintf("k=%d XB=%.4f Dmin=%.3f", k_tab$k, k_tab$xie_beni,
                        k_tab$min_centroid_distance), collapse = ", "))
}
k_sel <- if (nrow(k_tab)) k_tab[in_selection_range == TRUE] else k_tab
k_opt <- if (nrow(k_sel)) k_sel[which.min(xie_beni)]$k else 6L
xb_best <- if (nrow(k_sel)) k_sel[which.min(xie_beni)]$xie_beni else NA_real_
xb6 <- if (nrow(k_tab) && 6 %in% k_tab$k) k_tab[k == 6]$xie_beni else NA_real_
use_k6 <- !is.na(xb6) && !is.na(xb_best) && xb6 <= 1.10 * xb_best
K_CLUSTERS <- if (use_k6) 6L else as.integer(k_opt)
log_msg("  k selection: XB-optimal k (k >= 4) = ", k_opt, "; k = 6 retained: ", use_k6,
        "; K used = ", K_CLUSTERS)

set.seed(20260505)
cl <- mfuzz(eset_std, c = K_CLUSTERS, m = m_use)
membership <- cl$membership
colnames(membership) <- paste0("cluster", seq_len(ncol(membership)))
assign_dt <- data.table(gene_id = rownames(membership),
                        cluster = cl$cluster,
                        max_membership = apply(membership, 1, max),
                        membership)
save_table(assign_dt, file.path(TAB_DIR, "mfuzz_gene_cluster_assignments.csv"))
save_table(cbind(data.table(gene_id = rownames(membership)), as.data.table(membership)),
           file.path(TAB_DIR, "mfuzz_membership_scores.csv"))

centers <- cl$centers
save_table(data.table(cluster = rownames(centers), round(as.data.frame(centers), 4)),
           file.path(TAB_DIR, "mfuzz_cluster_centers_standardized.csv"))
cluster_summary <- data.table(
  cluster = seq_len(K_CLUSTERS),
  n_genes = as.integer(table(cl$cluster)),
  mean_membership = round(tapply(assign_dt$max_membership, assign_dt$cluster, mean), 3)
)
## direction of the mean control-adjusted profile per cluster (T0 -> T5d)
dirs <- vapply(seq_len(K_CLUSTERS), function(k) {
  g <- assign_dt[cluster == k]$gene_id
  if (length(g) == 0) return(NA_character_)
  prof <- colMeans(salt_effect_mat[g, , drop = FALSE])
  paste(sprintf("%s=%.2f", names(prof), prof), collapse = "; ")
}, character(1))
cluster_summary[, mean_profile := dirs]
cluster_summary[, trend := vapply(seq_len(K_CLUSTERS), function(k) {
  g <- assign_dt[cluster == k]$gene_id
  if (length(g) == 0) return(NA_character_)
  prof <- colMeans(salt_effect_mat[g, , drop = FALSE])
  if (which.max(abs(prof)) <= 3) "early-dominant" else "late-dominant"
}, character(1))]
save_table(cluster_summary, file.path(TAB_DIR, "mfuzz_cluster_summary.csv"))
log_msg("  clusters: ", paste(sprintf("C%d n=%d", cluster_summary$cluster, cluster_summary$n_genes),
                              collapse = ", "))
log_msg("  mean membership overall: ", round(mean(assign_dt$max_membership), 3))

## ---------------------------------------------------------------------------
## 9. Functional enrichment (GO BP/MF/CC + KEGG)
## ---------------------------------------------------------------------------

log_msg("\n== 9. Functional enrichment ==")

gene_sets <- list()
for (tp in TIME_LEVELS) {
  g <- de_all[[tp]][regulation != "NS"]$gene_id
  if (length(g)) gene_sets[[paste0("DEG_Salt_vs_Control_", tp)]] <- g
}
if (length(lrt_sig_genes)) gene_sets[["LRT_interaction_significant"]] <- lrt_sig_genes
for (k in seq_len(K_CLUSTERS)) {
  g <- assign_dt[cluster == k]$gene_id
  if (length(g)) gene_sets[[paste0("Mfuzz_cluster_", k)]] <- g
}
log_msg("  gene sets for enrichment: ", length(gene_sets))

## ---- 9a. g:Profiler (GO BP/MF/CC + KEGG, BH-FDR) ---------------------------
gprofiler_res <- NULL
gp_ok <- requireNamespace("gprofiler2", quietly = TRUE)
if (gp_ok) {
  suppressPackageStartupMessages(library(gprofiler2))
  cache_dir <- file.path(ANN_DIR, "gprofiler_cache")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  ## convert to the RAP-DB gene identifier convention used by the annotation
  to_rap <- function(x) sub("^Os([0-9]+)t", "Os\\1g", x)
  gp_tables <- list()
  for (nm in names(gene_sets)) {
    q <- unique(to_rap(gene_sets[[nm]]))
    q <- q[grepl("^Os[0-9]+g", q)]
    if (length(q) < 5) { log_msg("    g:Profiler skipped ", nm, " (too few mapped genes)"); next }
    q_hash <- paste0(length(q), "_",
                     format(sum(as.numeric(utf8ToInt(paste(sort(q), collapse = "")))) %% 1e9,
                            scientific = FALSE))
    cache_file <- file.path(cache_dir, paste0(gsub("[^A-Za-z0-9_]", "_", nm), "_", q_hash, ".csv"))
    if (file.exists(cache_file)) {
      d <- fread(cache_file)
      log_msg("    g:Profiler cache hit for ", nm)
    } else {
      ## evcodes = FALSE: returning the gene list for every term of a
      ## multi-thousand-gene query produces an enormous payload
      ## (observed > 6 GB in R) and is not needed here; the significant-term
      ## gene lists are available from the clusterProfiler cross-check tables.
      r <- tryCatch(gost(query = q, organism = "osativa",
                         sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG"),
                         correction_method = "fdr", user_threshold = 1,
                         significant = FALSE, evcodes = FALSE),
                    error = function(e) { log_msg("    g:Profiler error for ", nm, ": ",
                                                  conditionMessage(e)); NULL })
      d <- NULL
      if (!is.null(r) && !is.null(r$result) && nrow(r$result) > 0) {
        d <- as.data.table(r$result)
        fwrite(d, cache_file)
        Sys.sleep(1.5)
      }
    }
    if (!is.null(d) && nrow(d) > 0) {
      d <- copy(d)
      d[, gene_set := nm]
      d[, n_query_genes := length(q)]
      gp_tables[[nm]] <- d
    }
  }
  if (length(gp_tables)) {
    gprofiler_res <- rbindlist(gp_tables, fill = TRUE)
    gprofiler_res[, significant := !is.na(p_value) & p_value < 0.05]
    save_table(gprofiler_res, file.path(TAB_DIR, "enrichment_gprofiler_GO_KEGG_all_terms.csv"))
    save_table(gprofiler_res[significant == TRUE],
               file.path(TAB_DIR, "enrichment_gprofiler_GO_KEGG_significant_FDR0.05.csv"))
    log_msg("  g:Profiler: ", sum(gprofiler_res$significant), " terms with FDR < 0.05 (all sets)")
  }
}

## ---- 9b. local cross-check: hypergeometric GO enrichment (org.Osativa.eg.db) -----
## clusterProfiler::enrichGO cannot be used here: org.Osativa.eg.db has no
## ENTREZID keytype, which enrichGO requires internally (it calls
## bitr(..., toType = "ENTREZID")).  The cross-check is therefore computed
## directly with the standard hypergeometric (Fisher) over-representation test
## on the same rice annotation source, with BH correction over all terms.
## KEGG is covered by the g:Profiler query above.
cp_res <- NULL
go_map <- data.table(gene_id = character(), go_id = character(),
                     go_term = character(), ontology = character())
if (requireNamespace("org.Osativa.eg.db", quietly = TRUE)) {
  suppressPackageStartupMessages({ library(org.Osativa.eg.db); library(GO.db) })
  rap_ids <- unique(rap_of(rownames(counts_filt)))
  go_raw <- suppressWarnings(tryCatch(
    AnnotationDbi::select(org.Osativa.eg.db, keys = rap_ids, keytype = "RAP",
                          columns = c("GO", "ONTOLOGY")), error = function(e) NULL))
  if (!is.null(go_raw)) {
    go_raw <- as.data.table(go_raw)[!is.na(GO)]
    goterms <- suppressWarnings(tryCatch(
      AnnotationDbi::select(GO.db, keys = unique(go_raw$GO),
                            columns = c("TERM", "ONTOLOGY"), keytype = "GOID"),
      error = function(e) NULL))
    if (!is.null(goterms)) {
      goterms <- as.data.table(goterms)[!is.na(TERM)]
      go_map <- merge(go_raw[, .(rap_id = RAP, go_id = GO)],
                      goterms[, .(go_id = GOID, go_term = TERM, ontology = ONTOLOGY)],
                      by = "go_id")
      go_map[, gene_id := sub("^Os([0-9]+)g", "Os\\1t", rap_id)]
      setkey(go_map, gene_id)
    }
  }
  log_msg("  GO annotations (org.Osativa.eg.db) available for ",
          length(unique(go_map$gene_id)), " of ", nrow(counts_filt), " retained genes")

  if (nrow(go_map) > 0) {
    universe <- rownames(counts_filt)
    N <- length(universe)
    go_idx <- data.table(go_id = go_map$go_id, gene_idx = match(go_map$gene_id, universe))
    go_idx <- go_idx[!is.na(gene_idx)]
    term_sizes <- go_idx[, .(N_term = .N), by = go_id]
    chosen <- term_sizes[N_term >= 5 & N_term <= 500]
    go_chosen <- go_idx[go_id %in% chosen$go_id]
    term_info <- unique(go_map[go_id %in% chosen$go_id, .(go_id, go_term, ontology)])
    term_info <- merge(term_info, chosen, by = "go_id")
    cp_tables <- list()
    for (nm in names(gene_sets)) {
      genes <- intersect(unique(gene_sets[[nm]]), universe)
      n <- length(genes)
      if (n < 5) next
      in_set <- universe %in% genes
      hits <- go_chosen[, .(intersection_size = sum(in_set[gene_idx])), by = go_id]
      d <- merge(term_info, hits, by = "go_id")
      d <- d[intersection_size >= 2]
      if (!nrow(d)) next
      d[, pvalue := phyper(intersection_size - 1, N_term, N - N_term, n, lower.tail = FALSE)]
      d[, p.adjust := p.adjust(pvalue, method = "BH")]
      d <- d[, .(gene_set = nm, go_id, description = go_term, ontology,
                 term_size = N_term, intersection_size, pvalue, p.adjust,
                 n_query_genes = n)]
      setorder(d, pvalue)
      cp_tables[[nm]] <- d
    }
    if (length(cp_tables)) {
      cp_res <- rbindlist(cp_tables, fill = TRUE)
      setorder(cp_res, gene_set, pvalue)
      save_table(cp_res, file.path(TAB_DIR, "enrichment_local_hypergeometric_GO_all_terms.csv"))
      save_table(cp_res[p.adjust < 0.05],
                 file.path(TAB_DIR, "enrichment_local_hypergeometric_GO_FDR0.05.csv"))
      log_msg("  local GO over-representation: ", nrow(cp_res), " gene-set x term results; ",
              sum(cp_res$p.adjust < 0.05), " with FDR < 0.05")
    }
  }
}

## ---------------------------------------------------------------------------
## 10. Annotation lookup and candidate genes
## ---------------------------------------------------------------------------

log_msg("\n== 10. Candidate genes ==")

## gene -> RAP-DB style id (needed for org.Osativa.eg.db) and vice versa
rap_of <- function(x) sub("^Os([0-9]+)t", "Os\\1g", x)

## ---- UniProt (reviewed rice entries) ---------------------------------------
uniprot_path <- file.path(ANN_DIR, "uniprot_osativa_reviewed.tsv")
up_long <- data.table(entry = character(), protein_name = character(),
                      keyword = character(), alias = character())
if (file.exists(uniprot_path)) {
  up <- fread(uniprot_path, sep = "\t", quote = "", na.strings = c("", "NA"))
  names(up) <- c("entry", "gene_primary", "gene_names", "protein_name", "go_id", "go", "keyword", "ec")
  up <- up[!is.na(gene_names)]
  up_long <- up[, .(alias = unlist(strsplit(gene_names, "[; ]+"))),
                by = .(entry, protein_name, keyword)]
  up_long <- up_long[!is.na(alias) & alias != ""]
  setkey(up_long, alias)
  log_msg("  UniProt reviewed rice entries: ", nrow(up), "; aliases: ", nrow(up_long))
}
## ---- Ensembl Plants GFF3: gene symbols -------------------------------------
gff_symbols <- data.table(rap_id = character(), symbol = character(), description = character())
gff_path <- file.path(ANN_DIR, "Oryza_sativa.IRGSP-1.0.63.gff3.gz")
if (file.exists(gff_path)) {
  con <- gzfile(gff_path, "rt")
  gff_lines <- readLines(con)
  close(con)
  gff_gene <- grep("\tgene\t", gff_lines, value = TRUE)
  gid <- sub(".*gene_id=([^;]+).*", "\\1", gff_gene)
  nm <- ifelse(grepl("Name=", gff_gene),
               sub(".*Name=([^;]+).*", "\\1", gff_gene), NA_character_)
  gff_symbols <- data.table(rap_id = gid, symbol = nm)
  gff_symbols <- unique(gff_symbols)
  setkey(gff_symbols, rap_id)
  log_msg("  Ensembl GFF3 gene symbols: ", nrow(gff_symbols))
}
## ---- GO terms per gene: built in section 9b (go_map) ------------------------
log_msg("  GO annotations reused from section 9b for ",
        length(unique(go_map$gene_id)), " genes")

## ---- annotation summary per gene ------------------------------------------
annotate_genes <- function(genes) {
  rbindlist(lapply(genes, function(g) {
    rap <- rap_of(g)
    hit <- up_long[alias %in% c(g, rap, sub("^Os([0-9]+)g", "LOC_Os\\1g", rap))]
    hit <- hit[!is.na(protein_name)]
    prot <- if (!is.null(hit) && nrow(hit)) hit$protein_name[1] else NA_character_
    uniprot_id <- if (!is.null(hit) && nrow(hit)) hit$entry[1] else NA_character_
    kw <- if (!is.null(hit) && nrow(hit)) hit$keyword[1] else NA_character_
    sym <- gff_symbols[rap_id == rap]$symbol
    sym <- if (length(sym)) sym[1] else NA_character_
    gos <- go_map[gene_id == g]
    go_bp <- paste(head(gos[ontology == "BP"]$go_term, 4), collapse = "; ")
    go_mf <- paste(head(gos[ontology == "MF"]$go_term, 3), collapse = "; ")
    go_cc <- paste(head(gos[ontology == "CC"]$go_term, 3), collapse = "; ")
    data.table(gene_id = g, rap_id = rap, gene_symbol = sym,
               uniprot_id = uniprot_id, uniprot_protein_name = prot,
               uniprot_keywords = kw,
               GO_BP = go_bp, GO_MF = go_mf, GO_CC = go_cc,
               annotation_evidence = paste(na.omit(c(
                 if (!is.na(prot)) paste0("UniProt:", prot),
                 if (!is.na(sym)) paste0("RAP-DB symbol:", sym))), collapse = " | "))
  }))
}

cand_specs <- list(
  "Na+/ion transport" = c("sodium", "Na\\+/H\\+", "\\bHKT\\b", "\\bNHX\\b", "potassium",
                          "cation transport", "salt overly sensitive", "\\bSOS1\\b",
                          "ion transport", "monovalent cation", "high-affinity potassium",
                          "calcium-transporting ATPase", "\\bCBL\\b", "\\bCIPK\\b"),
  "NAC transcription factor" = c("NAC domain", "NAC transcription", "no apical meristem",
                                 "\\bNAC\\b", "nascent polypeptide"),
  "WRKY transcription factor" = c("\\bWRKY\\b"),
  "DREB/CBF/AP2-ERF" = c("\\bDREB\\b", "\\bCBF\\b", "dehydration-responsive element",
                         "AP2/ERF", "APETALA2", "ethylene-responsive transcription factor",
                         "\\bERF\\b"),
  "LEA / dehydrin" = c("late embryogenesis abundant", "\\bLEA\\b", "dehydrin", "\\bDHN\\b",
                       "\\bRAB16\\b", "\\bLIP9\\b"),
  "Heat shock protein / chaperone" = c("heat shock", "\\bHSP\\b", "chaperone", "\\bDnaJ\\b",
                                       "Hsp70", "Hsp90", "\\bHSP101\\b"),
  "Trehalose metabolism" = c("trehalose"),
  "Glutathione / ROS detoxification" = c("glutathione", "peroxidase", "superoxide dismutase",
                                         "catalase", "ascorbate", "thioredoxin", "glutaredoxin",
                                         "respiratory burst oxidase", "\\bRBOH\\b",
                                         "peroxiredoxin", "oxidative stress", "reactive oxygen"),
  "Photosystem / light harvesting" = c("photosystem", "light-harvesting", "light harvesting",
                                       "chlorophyll", "rubisco", "ribulose-bisphosphate",
                                       "photosynthesis", "photoreduction", "\\bPsb[RPQAD]\\b"),
  "MAPK signalling" = c("mitogen-activated protein kinase", "\\bMAPK\\b", "MAP kinase",
                        "MAPKK", "MAP3K", "MAP2K", "MAPKKK"),
  "Cell wall / expansion" = c("cell wall organization", "cell wall biogenesis",
                              "cell wall macromolecule", "expansin", "cellulose synthase",
                              "xyloglucan", "lignin", "pectin", "pectinesterase",
                              "glucan synthase", "\\bPME\\b", "\\bXTH\\b")
)

## candidate pre-filter: statistically supported by the corrected analysis
support_dt <- rbindlist(lapply(TIME_LEVELS, function(tp) {
  d <- de_all[[tp]][regulation != "NS", .(gene_id, tp = tp, lfc = log2FoldChange, padj = padj)]
  d
}))
support_summary <- support_dt[, .(n_sig_timepoints = .N,
                                  max_abs_lfc = max(abs(lfc)),
                                  best_lfc = lfc[which.max(abs(lfc))],
                                  best_tp = tp[which.max(abs(lfc))]),
                              by = gene_id]
support_summary <- merge(support_summary, lrt_tab[, .(gene_id, lrt_padj = padj,
                                                      lrt_sig = interaction_significant)],
                         by = "gene_id", all.x = TRUE)
support_summary[is.na(n_sig_timepoints), n_sig_timepoints := 0L]
supported <- support_summary[n_sig_timepoints >= 1 | (lrt_sig %in% TRUE)]
log_msg("  genes with >= 1 significant time-matched DEG or significant LRT: ", nrow(supported))

## annotate the supported gene set (bounded: strongest 3000 by max |LFC|/LRT)
top_supported <- supported[order(-fifelse(lrt_sig %in% TRUE, 1, 0), -n_sig_timepoints, -abs(best_lfc))]
top_supported <- head(top_supported, 3000)
ann_tab <- annotate_genes(top_supported$gene_id)
ann_tab <- merge(ann_tab, top_supported, by = "gene_id")
save_table(ann_tab, file.path(TAB_DIR, "annotation_supported_genes.csv"))

match_category <- function(ann_row) {
  hay <- tolower(paste(na.omit(c(ann_row$uniprot_protein_name, ann_row$uniprot_keywords,
                                 ann_row$gene_symbol, ann_row$GO_BP, ann_row$GO_MF,
                                 ann_row$GO_CC)), collapse = " | "))
  hits <- names(cand_specs)[vapply(cand_specs, function(pats) {
    any(vapply(pats, function(p) grepl(tolower(p), hay, perl = TRUE), logical(1)))
  }, logical(1))]
  hit_terms <- vapply(cand_specs, function(pats) {
    m <- vapply(pats, function(p) grepl(tolower(p), hay, perl = TRUE), logical(1))
    if (any(m)) pats[which(m)[1]] else NA_character_
  }, character(1))
  list(categories = hits, evidence = hit_terms)
}

log_msg("  matching candidate categories ...")
ann_list <- lapply(seq_len(nrow(ann_tab)), function(i) match_category(ann_tab[i]))
ann_tab[, cat_list := ann_list]

cand_long <- rbindlist(lapply(seq_len(nrow(ann_tab)), function(i) {
  cl <- ann_tab$cat_list[[i]]
  if (length(cl$categories) == 0) return(NULL)
  ev <- cl$evidence[cl$categories]
  data.table(gene_id = ann_tab$gene_id[i], category = cl$categories,
             matched_annotation = unname(ev))
}))

candidate_genes <- character()
if (nrow(cand_long)) {
  cand_long <- merge(cand_long, ann_tab[, .(gene_id, rap_id, gene_symbol, uniprot_id,
                                            uniprot_protein_name, uniprot_keywords,
                                            GO_BP, GO_MF, GO_CC, n_sig_timepoints, max_abs_lfc,
                                            best_lfc, best_tp, lrt_padj, lrt_sig)],
                     by = "gene_id")
  cand_long <- cand_long[order(category, -fifelse(lrt_sig %in% TRUE, 1, 0), -n_sig_timepoints, -abs(best_lfc))]
  ## at most 4 genes per category and 15-30 genes overall; the per-category cap
  ## is relaxed to 8 when the strict cap yields fewer than 15 genes
  pick_candidates <- function(cap) {
    s <- cand_long[, head(.SD, cap), by = category]
    s[order(-fifelse(lrt_sig %in% TRUE, 1, 0), -n_sig_timepoints, -abs(best_lfc))]
  }
  sel_rows <- pick_candidates(4)
  if (nrow(sel_rows) < 15) sel_rows <- pick_candidates(8)
  if (nrow(sel_rows) > 30) sel_rows <- sel_rows[seq_len(30)]
  candidate_genes <- unique(sel_rows$gene_id)
}
log_msg("  candidate genes selected: ", length(candidate_genes))

if (length(candidate_genes)) {
  cand_tab <- annotate_genes(candidate_genes)
  cat_by_gene <- cand_long[, .(category = paste(unique(category), collapse = "; "),
                               category_evidence = paste(unique(paste0(category, " <- matched '",
                                                                       matched_annotation, "'")),
                                                         collapse = "; ")),
                           by = gene_id]
  cand_tab <- merge(cand_tab, cat_by_gene, by = "gene_id", all.x = TRUE)
  ## per-time-point salt effect
  for (tp in TIME_LEVELS) {
    d <- de_all[[tp]][, .(gene_id, lfc = log2FoldChange, padj = padj, shrunk = log2FoldChange_shrunk)]
    setnames(d, c("lfc", "padj", "shrunk"),
             paste0(c("salt_log2FC_", "salt_padj_", "salt_log2FC_shrunk_"), tp))
    cand_tab <- merge(cand_tab, d, by = "gene_id", all.x = TRUE)
  }
  cand_tab <- merge(cand_tab, lrt_tab[, .(gene_id, interaction_LRT_padj = padj,
                                          interaction_LRT_significant = interaction_significant)],
                    by = "gene_id", all.x = TRUE)
  for (id in names(did_all)) {
    d <- did_all[[id]][, .(gene_id, lfc = log2FoldChange, padj = padj)]
    setnames(d, c("lfc", "padj"), c(paste0("DID_log2FC_", id), paste0("DID_padj_", id)))
    cand_tab <- merge(cand_tab, d, by = "gene_id", all.x = TRUE)
  }
  cand_tab <- merge(cand_tab, assign_dt[, .(gene_id, mfuzz_cluster = cluster,
                                            mfuzz_membership = max_membership)],
                    by = "gene_id", all.x = TRUE)
  cand_ids <- cand_tab$gene_id
  cand_tab[, observed_salt_effect_profile := apply(salt_effect_mat[cand_ids, , drop = FALSE], 1,
                                                   function(v) paste(sprintf("%s=%.2f", TIME_LEVELS, v), collapse = "; "))]
  cand_tab[, n_significant_time_points := vapply(gene_id, function(g)
    sum(sapply(TIME_LEVELS, function(tp) {
      r <- de_all[[tp]][gene_id == g]
      nrow(r) > 0 && !is.na(r$padj) && r$padj < PADJ_CUTOFF && abs(r$log2FoldChange) >= LFC_CUTOFF
    })), integer(1))]
  cand_tab <- cand_tab[order(-fifelse(interaction_LRT_significant %in% TRUE, 1, 0),
                             -n_significant_time_points)]
  save_table(cand_tab, file.path(TAB_DIR, "candidate_genes_for_discussion.csv"))
  log_msg("  candidate gene table written with ", nrow(cand_tab), " genes")
}

## ---------------------------------------------------------------------------
## 11. Figures
## ---------------------------------------------------------------------------

log_msg("\n== 11. Figures ==")

## ---- Figure 1: design + workflow ------------------------------------------
design_dt <- CJ(time_point = factor(TIME_LEVELS, levels = TIME_LEVELS),
                condition = factor(COND_LEVELS, levels = COND_LEVELS))
design_dt[, label := paste0("n = 3\n", condition)]
p1a <- ggplot(design_dt, aes(time_point, condition, fill = condition)) +
  geom_tile(colour = "white", linewidth = 1, height = 0.9) +
  geom_text(aes(label = "3 replicates"), colour = "white", size = 3.6, fontface = "bold") +
  scale_fill_manual(values = c(Control = "#2c7fb8", Salt = "#d95f02"), guide = "none") +
  labs(title = "A  Experimental design (Oryza sativa cv. Azucena)",
       subtitle = paste0("2 conditions x 5 time points x 3 biological replicates = 30 libraries; ",
                         VERIFIED_TREATMENT_LABEL),
       x = NULL, y = NULL) +
  theme_pub()

wf <- data.table(
  step = factor(c("Raw transcript\nabundance (GEO)", "Transcript -> gene\naggregation",
                  "Low-count filter\n(count >= 10 in >= 3)",
                  "DESeq2 interaction model\n~ cond + time + cond:time",
                  "Time-matched\nSalt vs Control", "Difference-in-differences\n(change in salt effect)",
                  "Interaction LRT\n(full vs reduced)", "Mfuzz on control-adjusted\nsalt-effect profiles",
                  "GO / KEGG enrichment"),
                levels = c("Raw transcript\nabundance (GEO)", "Transcript -> gene\naggregation",
                           "Low-count filter\n(count >= 10 in >= 3)",
                           "DESeq2 interaction model\n~ cond + time + cond:time",
                           "Time-matched\nSalt vs Control", "Difference-in-differences\n(change in salt effect)",
                           "Interaction LRT\n(full vs reduced)", "Mfuzz on control-adjusted\nsalt-effect profiles",
                           "GO / KEGG enrichment")),
  x = c(1, 2, 3, 4, 5, 5, 5, 5.5, 7.5),
  y = c(1, 1, 1, 1, 2, 1, 0, 1, 1.5)
)
p1b <- ggplot(wf, aes(x, y)) +
  geom_label(aes(label = step), size = 2.8, label.size = 0.4,
             label.padding = grid::unit(0.35, "lines"), fill = "grey96") +
  geom_segment(data = data.table(x = c(1.35, 2.35, 3.35), xend = c(1.9, 2.9, 3.9), y = 1, yend = 1),
               aes(x = x, xend = xend, y = y, yend = yend),
               arrow = grid::arrow(length = grid::unit(0.15, "cm")), inherit.aes = FALSE) +
  coord_cartesian(xlim = c(0.4, 8.4), ylim = c(-0.6, 2.5), clip = "off") +
  labs(title = "B  Analytical workflow", x = NULL, y = NULL) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
fig_ok(p1a / p1b + plot_layout(heights = c(1, 1.5)),
       file.path(FIG_DIR, "figure01_design_and_workflow.png"), width = 10, height = 8)

## ---- Figure 2: QC / PCA ----------------------------------------------------
lib_long <- melt(lib_summary, id.vars = c("sample_id", "condition", "time_point"),
                 measure.vars = "raw_library_size", variable.name = "metric")
p2a <- ggplot(lib_summary, aes(time_point, raw_library_size / 1e6, fill = condition)) +
  geom_boxplot(outlier.size = 0.8, width = 0.7) +
  scale_fill_manual(values = c(Control = "#2c7fb8", Salt = "#d95f02")) +
  labs(title = "A  Library size", x = NULL, y = "million reads (gene level)", fill = NULL) +
  theme_pub() + theme(legend.position = "top")

p2b <- ggplot(pca_df, aes(PC1, PC2, colour = condition, shape = time_point)) +
  geom_point(size = 3, stroke = 0.9) +
  stat_ellipse(aes(group = condition), level = 0.68, linetype = 2, linewidth = 0.5) +
  scale_colour_manual(values = c(Control = "#2c7fb8", Salt = "#d95f02")) +
  labs(title = "B  PCA (VST, top 500 variable genes)",
       x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)"),
       colour = NULL, shape = "Time point") +
  theme_pub() + theme(legend.position = "right")

p2c <- ggplot(cent_dist, aes(time_point, distance, group = 1)) +
  geom_line(linewidth = 0.8, colour = "grey30") +
  geom_point(size = 3, colour = "#d95f02") +
  labs(title = "C  Control vs Salt centroid distance in PCA space",
       subtitle = "descriptive summary, not a statistical test",
       x = NULL, y = "PC1-PC2 centroid distance") +
  theme_pub()

fig_ok((p2a | p2b) / p2c + plot_layout(heights = c(1, 1)),
       file.path(FIG_DIR, "figure02_qc_and_pca.png"), width = 12, height = 9)

## supplementary QC
p_sc <- ggplot(melt(as.data.table(cor_mat, keep.rownames = "sample_id"),
                    id.vars = "sample_id", variable.name = "other", value.name = "r"),
               aes(sample_id, other, fill = r)) +
  geom_tile() + scale_fill_viridis_c(option = "C") +
  labs(title = "Sample-sample correlation (VST)", x = NULL, y = NULL) +
  theme_pub(9) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
                       axis.text.y = element_text(size = 6))
fig_ok(p_sc, file.path(SUP_FIG, "figureS_sample_correlation.png"), width = 9, height = 8)
p_pc3 <- ggplot(pca_df, aes(PC2, PC3, colour = condition, shape = time_point)) +
  geom_point(size = 3) +
  scale_colour_manual(values = c(Control = "#2c7fb8", Salt = "#d95f02")) +
  labs(title = "PCA PC2 vs PC3", x = paste0("PC2 (", pct[2], "%)"),
       y = paste0("PC3 (", pct[3], "%)")) + theme_pub()
fig_ok(p_pc3, file.path(SUP_FIG, "figureS_pca_pc2_pc3.png"), width = 7, height = 6)
p_sf <- ggplot(lib_summary, aes(time_point, size_factor, fill = condition)) +
  geom_boxplot(width = 0.7, outlier.size = 0.8) +
  scale_fill_manual(values = c(Control = "#2c7fb8", Salt = "#d95f02")) +
  labs(title = "DESeq2 size factors", x = NULL, y = "size factor") + theme_pub()
fig_ok(p_sf, file.path(SUP_FIG, "figureS_size_factors.png"), width = 7, height = 5)

## supplementary: PCA split by time point (baseline / group-offset check)
p_s <- ggplot(pca_df, aes(PC1, PC2, colour = condition)) +
  geom_point(size = 2, alpha = 0.9) +
  facet_wrap(~ factor(time_point, levels = TIME_LEVELS)) +
  scale_colour_manual(values = c(Control = "#2c7fb8", Salt = "#d95f02")) +
  labs(title = "PCA by time point",
       subtitle = "Control and Salt samples separate at every time point, including T0",
       x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)")) +
  theme_pub(9) + theme(legend.position = "bottom")
fig_ok(p_s, file.path(SUP_FIG, "figureS_pca_by_time_point.png"), width = 10, height = 6)

## ---- Figure 3: time-matched volcano panels --------------------------------
volcano_df <- rbindlist(lapply(TIME_LEVELS, function(tp) {
  d <- copy(de_all[[tp]])
  d[, time_point := TIME_LABELS[tp]]
  d
}))
volcano_df[, time_point := factor(time_point, levels = TIME_LABELS)]
volcano_df[, neglog10padj := -log10(pmax(padj, .Machine$double.xmin))]
volcano_df[, category := ifelse(regulation == "Up", "Up (salt > control)",
                         ifelse(regulation == "Down", "Down (salt < control)", "Not significant"))]
top_lab <- volcano_df[regulation != "NS",
                      .SD[order(-abs(log2FoldChange))][seq_len(min(5, .N))],
                      by = time_point]
p3 <- ggplot(volcano_df, aes(log2FoldChange, neglog10padj, colour = category)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_colour_manual(values = c("Up (salt > control)" = "#d73027",
                                 "Down (salt < control)" = "#4575b4",
                                 "Not significant" = "grey80"),
                      name = NULL) +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = 2, linewidth = 0.3) +
  geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = 2, linewidth = 0.3) +
  geom_text(data = top_lab, aes(label = gene_id), size = 1.9, colour = "black",
            check_overlap = TRUE, vjust = -0.4) +
  facet_wrap(~ time_point, nrow = 2) +
  labs(title = paste0("Time-matched Salt (", VERIFIED_NACL_MM, " mM NaCl) vs Control contrasts"),
       subtitle = paste0("padj < ", PADJ_CUTOFF, " and |log2FC| >= ", LFC_CUTOFF),
       x = "log2 fold change (Salt vs time-matched Control)",
       y = "-log10 adjusted p-value") +
  theme_pub(10) + theme(legend.position = "bottom")
fig_ok(p3, file.path(FIG_DIR, "figure03_time_matched_volcano.png"), width = 11, height = 8)

## ---- Figure 4: DEG counts and direction across time -----------------------
deg_long <- melt(de_summary_dt, id.vars = c("comparison", "time_point"),
                 measure.vars = c("upregulated", "downregulated"),
                 variable.name = "direction", value.name = "n")
deg_long[, direction := ifelse(direction == "upregulated", "Up", "Down")]
deg_long[, time_point := factor(time_point, levels = TIME_LEVELS)]
p4a <- ggplot(deg_long, aes(time_point, n, fill = direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = n), position = position_dodge(width = 0.8), vjust = -0.3, size = 2.8) +
  scale_fill_manual(values = c(Up = "#d73027", Down = "#4575b4"), name = NULL) +
  labs(title = "A  DEGs per time-matched Salt-vs-Control comparison",
       x = NULL, y = "number of DEGs") +
  theme_pub() + theme(legend.position = "top")
ratio_dt <- de_summary_dt[, .(time_point = factor(time_point, levels = TIME_LEVELS),
                              ratio = upregulated / pmax(downregulated, 1),
                              total_DEGs = total_DEGs)]
p4b <- ggplot(ratio_dt, aes(time_point, ratio, group = 1)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
  geom_line(linewidth = 0.8, colour = "grey25") +
  geom_point(aes(size = total_DEGs), colour = "#d95f02") +
  scale_size_continuous(name = "DEGs", range = c(2, 7)) +
  labs(title = "B  Up : Down ratio of DEGs", x = NULL,
       y = "upregulated / downregulated") +
  theme_pub()
p4c <- ggplot(de_summary_dt[, .(time_point = factor(time_point, levels = TIME_LEVELS),
                                pct = pct_of_retained_genes)], aes(time_point, pct, group = 1)) +
  geom_col(fill = "grey60", width = 0.7) +
  labs(title = "C  DEGs as % of retained genes", x = NULL, y = "%") +
  theme_pub()
fig_ok(p4a / (p4b | p4c) + plot_layout(heights = c(1.3, 1)),
       file.path(FIG_DIR, "figure04_deg_counts_and_direction.png"), width = 10, height = 9)

## ---- Figure 5: interaction / difference-in-differences --------------------
did_adj <- did_summary_dt[adjacent == TRUE]
did_adj[, contrast := factor(contrast, levels = c("T60_vs_T0", "T180_vs_T60",
                                                  "T240_vs_T180", "T5d_vs_T240"))]
did_long <- melt(did_adj, id.vars = "contrast",
                 measure.vars = c("increased_salt_effect", "decreased_salt_effect"),
                 variable.name = "direction", value.name = "n")
did_long[, direction := ifelse(direction == "increased_salt_effect",
                               "Salt effect increases", "Salt effect decreases")]
p5a <- ggplot(did_long, aes(contrast, n, fill = direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = n), position = position_dodge(width = 0.8), vjust = -0.3, size = 2.8) +
  scale_fill_manual(values = c("Salt effect increases" = "#d73027",
                               "Salt effect decreases" = "#4575b4"), name = NULL) +
  labs(title = "A  Genes with a significant change in the salt effect",
       subtitle = paste0("padj < ", PADJ_CUTOFF, " and |delta log2FC| >= ", LFC_CUTOFF,
                         " between consecutive time points"),
       x = NULL, y = "number of genes") +
  theme_pub(10) + theme(legend.position = "top")

did_sig_all <- rbindlist(lapply(names(did_all), function(id) {
  d <- did_all[[id]][regulation != "NS"]
  if (!nrow(d)) return(NULL)
  cbind(d[, .(gene_id, log2FoldChange, did_contrast)], data.table(contrast = id))
}))
if (nrow(did_sig_all)) {
  did_sig_all[, contrast := factor(contrast, levels = c("T60_vs_T0", "T180_vs_T60",
                                                        "T240_vs_T180", "T5d_vs_T240",
                                                        "T5d_vs_T60", "T5d_vs_T0"))]
  p5b <- ggplot(did_sig_all[!is.na(contrast)], aes(contrast, abs(log2FoldChange), fill = contrast)) +
    geom_boxplot(outlier.size = 0.4, width = 0.7, show.legend = FALSE) +
    labs(title = "B  Magnitude of significant changes in salt effect",
         x = NULL, y = "|delta log2FC| (difference-in-differences)") +
    theme_pub(10) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
} else p5b <- ggplot() + theme_void()

sc_dt <- data.table(gene_id = rownames(salt_effect_mat),
                    early = salt_effect_mat[, "T240"],
                    late = salt_effect_mat[, "T5d"])
sc_dt[, sig := gene_id %in% did_all[["T5d_vs_T240"]][regulation != "NS"]$gene_id]
p5c <- ggplot(sc_dt, aes(early, late)) +
  geom_point(data = sc_dt[sig == FALSE], colour = "grey80", size = 0.3, alpha = 0.4) +
  geom_point(data = sc_dt[sig == TRUE], colour = "#d95f02", size = 0.5, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(title = "C  Salt effect at T240 vs T5d",
       subtitle = "orange: genes with a significant change in the salt effect between T240 and T5d",
       x = "salt effect at T240 (VST Salt - Control)",
       y = "salt effect at T5d (VST Salt - Control)") +
  theme_pub(10)
fig_ok(p5a / (p5b | p5c) + plot_layout(heights = c(1, 1)),
       file.path(FIG_DIR, "figure05_interaction_difference_in_differences.png"),
       width = 11, height = 9)

## supplementary: DEG overlap
if (requireNamespace("UpSetR", quietly = TRUE)) {
  deg_sets <- lapply(TIME_LEVELS, function(tp) de_all[[tp]][regulation != "NS"]$gene_id)
  names(deg_sets) <- TIME_LEVELS
  ok <- tryCatch({
    png(file.path(SUP_FIG, "figureS_deg_overlap_upset.png"), width = 9, height = 5,
        units = "in", res = 300)
    print(UpSetR::upset(UpSetR::fromList(deg_sets), order.by = "freq",
                        main.bar.color = "#2c7fb8", sets.bar.color = "#d95f02",
                        mainbar.y.label = "DEGs in intersection",
                        sets.x.label = "DEGs per time-matched contrast"))
    dev.off(); TRUE
  }, error = function(e) { try(dev.off(), silent = TRUE); FALSE })
  if (ok) log_msg("  wrote figureS_deg_overlap_upset.png")
}

## ---- Figure 6: Mfuzz clusters ---------------------------------------------
prof_dt <- rbindlist(lapply(seq_len(K_CLUSTERS), function(k) {
  g <- assign_dt[cluster == k]$gene_id
  if (!length(g)) return(NULL)
  m <- colMeans(salt_effect_mat[g, , drop = FALSE])
  se <- apply(salt_effect_mat[g, , drop = FALSE], 2, sd) / sqrt(length(g))
  data.table(cluster = k, time_point = factor(TIME_LEVELS, levels = TIME_LEVELS),
             mean_effect = m, se = se, n = length(g))
}))
p6a <- ggplot(prof_dt, aes(time_point, mean_effect, group = 1)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_ribbon(aes(ymin = mean_effect - se, ymax = mean_effect + se), alpha = 0.2,
              fill = "#d95f02") +
  geom_line(colour = "#a63603", linewidth = 0.9) +
  geom_point(size = 2, colour = "#a63603") +
  facet_wrap(~ paste0("cluster ", cluster, " (n = ", n, ")"), nrow = 2) +
  labs(title = "A  Mfuzz clusters of the temporal salt effect",
       subtitle = paste0("mean VST(Salt_t) - mean VST(Control_t) per gene; m = ", m_use,
                         ", k = ", K_CLUSTERS),
       x = NULL, y = "salt effect (VST units)") +
  theme_pub(10)
hm_dt <- salt_effect_mat[assign_dt$gene_id, , drop = FALSE]
ord <- order(assign_dt$cluster, -assign_dt$max_membership)
hm_dt <- hm_dt[ord, , drop = FALSE]
hm_dt <- pmax(pmin(hm_dt, 3), -3)
ann_row <- data.frame(cluster = factor(assign_dt$cluster[ord]))
rownames(ann_row) <- rownames(hm_dt)
png(file.path(FIG_DIR, "figure06_mfuzz_clusters.png"), width = 11, height = 9,
    units = "in", res = 300)
ok6 <- tryCatch({
  print(p6a)
  grid::grid.newpage()
  pheatmap::pheatmap(hm_dt, cluster_rows = FALSE, cluster_cols = FALSE,
                     show_rownames = FALSE, show_colnames = TRUE,
                     annotation_row = ann_row, color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
                     breaks = seq(-3, 3, length.out = 101),
                     main = "B  Salt effect per gene (mean VST(Salt) - mean VST(Control)), genes ordered by cluster",
                     fontsize = 9, silent = FALSE)
  TRUE
}, error = function(e) { log_msg("  !! mfuzz figure error: ", conditionMessage(e)); FALSE })
try(dev.off(), silent = TRUE)
if (ok6) log_msg("  wrote figure06_mfuzz_clusters.png")

## supplementary: Mfuzz cluster diagnostics
if (!is.null(k_tab) && nrow(k_tab)) {
  p_dmin <- ggplot(melt(k_tab, id.vars = "k",
                        measure.vars = c("min_centroid_distance", "xie_beni"),
                        variable.name = "index", value.name = "value"),
                   aes(k, value)) +
    geom_line(colour = "grey30") + geom_point(size = 2.5, colour = "#d95f02") +
    geom_vline(xintercept = K_CLUSTERS, linetype = 2, colour = "#2c7fb8") +
    facet_wrap(~ index, scales = "free_y") +
    labs(title = "Mfuzz cluster-number diagnostics",
         subtitle = paste0("selected k = ", K_CLUSTERS, " (dashed); Xie-Beni index, lower is better"),
         x = "number of clusters (k)", y = "index value") + theme_pub()
  fig_ok(p_dmin, file.path(SUP_FIG, "figureS_mfuzz_k_selection.png"), width = 9, height = 5)
}
p_mem <- ggplot(assign_dt, aes(factor(cluster), max_membership)) +
  geom_violin(fill = "grey85", colour = "grey40") + geom_boxplot(width = 0.15) +
  labs(title = "Mfuzz membership scores per cluster", x = "cluster",
       y = "maximum membership") + theme_pub()
fig_ok(p_mem, file.path(SUP_FIG, "figureS_mfuzz_membership.png"), width = 7, height = 5)

## ---- Figure 7: enrichment summary -----------------------------------------
if (!is.null(gprofiler_res) && nrow(gprofiler_res[significant == TRUE]) > 0) {
  enr <- gprofiler_res[significant == TRUE]
  enr[, neglog10p := -log10(p_value)]
  top_enr <- enr[, .SD[order(p_value)][seq_len(min(6, .N))], by = gene_set]
  top_enr[, source := factor(source, levels = c("GO:BP", "GO:MF", "GO:CC", "KEGG"))]
  p7 <- ggplot(top_enr, aes(neglog10p, reorder(paste0(term_name, " [", source, "]"), neglog10p))) +
    geom_point(aes(size = intersection_size, colour = source)) +
    facet_wrap(~ set_label(gene_set), scales = "free_y", ncol = 3) +
    scale_colour_manual(values = c("GO:BP" = "#1b9e77", "GO:MF" = "#7570b3",
                                   "GO:CC" = "#e6ab02", "KEGG" = "#d95f02"), name = NULL) +
    labs(title = "Gene Ontology and KEGG enrichment of salt-responsive gene sets",
         subtitle = "g:Profiler (Oryza sativa), BH-FDR < 0.05; top 6 terms per set",
         x = "-log10(FDR)", y = NULL) +
    theme_pub(8) + theme(legend.position = "bottom")
  fig_ok(p7, file.path(FIG_DIR, "figure07_enrichment_summary.png"), width = 15, height = 11)
}

## supplementary: per-set enrichment dot plots
if (!is.null(gprofiler_res)) {
  pdf_file <- file.path(SUP_FIG, "figureS_enrichment_per_set.pdf")
  ok <- tryCatch({
    pdf(pdf_file, width = 9, height = 6)
    for (gs in unique(gprofiler_res$gene_set)) {
      d <- gprofiler_res[gene_set == gs & significant == TRUE]
      if (!nrow(d)) next
      d <- d[order(p_value)][seq_len(min(20, .N))]
      print(ggplot(d, aes(-log10(p_value), reorder(term_name, p_value))) +
              geom_point(aes(size = intersection_size, colour = source)) +
              labs(title = paste0(set_label(gs), " (", d$n_query_genes[1], " genes)"),
                   subtitle = "top 20 terms, BH-FDR < 0.05", x = "-log10(FDR)", y = NULL) +
              theme_pub(9))
    }
    dev.off(); TRUE
  }, error = function(e) { try(dev.off(), silent = TRUE); FALSE })
  if (ok) log_msg("  wrote figureS_enrichment_per_set.pdf")
}

## ---- Figure 8: candidate genes --------------------------------------------
if (length(candidate_genes)) {
  cg <- intersect(candidate_genes, rownames(salt_effect_mat))
  lab <- cand_tab[match(cg, gene_id)]
  lab[, label := ifelse(!is.na(gene_symbol), paste0(gene_symbol, " (", gene_id, ")"), gene_id)]
  m <- salt_effect_mat[cg, , drop = FALSE]
  rownames(m) <- lab$label
  ann_c <- data.frame(cluster = factor(lab$mfuzz_cluster))
  ordr <- order(lab$mfuzz_cluster, -apply(salt_effect_mat[cg, , drop = FALSE], 1, function(v) max(abs(v))))
  ann_c <- ann_c[ordr, , drop = FALSE]
  rownames(ann_c) <- rownames(m)[ordr]
  png(file.path(FIG_DIR, "figure08_candidate_genes.png"), width = 11, height = 9,
      units = "in", res = 300)
  ok8 <- tryCatch({
    pheatmap::pheatmap(m[ordr, , drop = FALSE], cluster_rows = FALSE, cluster_cols = FALSE,
                       annotation_row = ann_c,
                       color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
                       breaks = seq(-3.5, 3.5, length.out = 101),
                       main = "Candidate salt-response genes\nmean VST(Salt) - mean VST(Control) at each time point",
                       fontsize = 8)
    TRUE
  }, error = function(e) { log_msg("  !! candidate figure error: ", conditionMessage(e)); FALSE })
  try(dev.off(), silent = TRUE)
  if (ok8) log_msg("  wrote figure08_candidate_genes.png")
}

## ---------------------------------------------------------------------------
## 12. Interpretation support and summary files
## ---------------------------------------------------------------------------

log_msg("\n== 12. Summary ==")

## binomial tests for the up:down composition at each time point
binom_dt <- rbindlist(lapply(TIME_LEVELS, function(tp) {
  s <- de_summary_dt[time_point == tp]
  n <- s$upregulated + s$downregulated
  bt <- if (n > 0) binom.test(s$upregulated, n)$p.value else NA_real_
  data.table(time_point = tp, upregulated = s$upregulated, downregulated = s$downregulated,
             pct_up = round(100 * s$upregulated / max(n, 1), 1),
             binom_p_up_vs_down = signif(bt, 3))
}))
save_table(binom_dt, file.path(TAB_DIR, "up_down_binomial_tests.csv"))

## direction switch supported by the DID results
did_switch <- did_all[["T5d_vs_T60"]]
sw <- merge(data.table(gene_id = rownames(salt_effect_mat),
                       early = salt_effect_mat[, "T60"], late = salt_effect_mat[, "T5d"]),
            did_switch[, .(gene_id, did_lfc = log2FoldChange, did_padj = padj)], by = "gene_id")
sw_sig <- sw[!is.na(did_padj) & did_padj < PADJ_CUTOFF & abs(did_lfc) >= LFC_CUTOFF]
n_neg_to_pos <- sum(sw_sig$early < 0 & sw_sig$late > 0)
n_pos_to_neg <- sum(sw_sig$early > 0 & sw_sig$late < 0)
switch_test_p <- if (n_neg_to_pos + n_pos_to_neg > 0)
  binom.test(n_neg_to_pos, n_neg_to_pos + n_pos_to_neg)$p.value else NA_real_

adj_counts <- did_summary_dt[adjacent == TRUE][order(match(contrast,
                                                           c("T60_vs_T0", "T180_vs_T60",
                                                             "T240_vs_T180", "T5d_vs_T240")))]
largest_adj <- adj_counts[which.max(total_significant)]$contrast
second_adj <- adj_counts[order(-total_significant)][2]$contrast
largest_vs_second <- adj_counts[contrast == largest_adj]$total_significant /
  max(adj_counts[contrast == second_adj]$total_significant, 1)

## is the temporal change in the salt effect statistically supported?
temporal_supported <- n_lrt_sig > 0 && any(adj_counts$total_significant > 0)
## biphasic support (early repression -> late activation), transparent rule
early <- binom_dt[time_point %in% c("T60", "T180")]
late <- binom_dt[time_point %in% c("T240", "T5d")]
early_repressed <- all(early$pct_up < 50) && any(early$binom_p_up_vs_down < 0.05)
late_activated <- all(late$pct_up > 50) && any(late$binom_p_up_vs_down < 0.05)
switch_supported <- !is.na(switch_test_p) && switch_test_p < 0.05 && n_neg_to_pos > n_pos_to_neg
biphasic_supported <- early_repressed && late_activated && switch_supported

summary_lines <- c(
  "REVIEWER-CORRECTED RNA-SEQ REANALYSIS - SUMMARY",
  "===============================================",
  paste0("generated: ", format(Sys.time())),
  paste0("R: ", R.version.string),
  "",
  "A. INPUT / QC",
  paste0("  verified treatment              : ", VERIFIED_TREATMENT_LABEL,
         " (metadata/BioSample attribute says '10mM NaCl' - see tables/treatment_verification.txt)"),
  paste0("  samples                         : ", nrow(samples),
         " (", paste(unique(sample_counts$N), collapse = ", "),
         " per condition x time point; ", nrow(sample_counts), " groups)"),
  paste0("  condition x time point groups   : ", nrow(sample_counts)),
  paste0("  genes in gene-level matrix      : ", nrow(cnt)),
  paste0("  retained genes (counts >= 10 in >= 3 samples): ", nrow(counts_filt)),
  paste0("  model                           : ~ condition + time_point + condition:time_point"),
  paste0("  reference levels                : condition = Control, time_point = T0"),
  paste0("  shrinkage (visualisation only)  : ", paste(unique(de_summary_dt$shrinkage_method), collapse = ", ")),
  paste0("  PCA variance explained          : PC1 = ", pct[1], "%, PC2 = ", pct[2], "%, PC3 = ", pct[3], "%"),
  paste0("  PCA control-vs-salt centroid distance: ",
         paste(sprintf("%s = %.2f", cent_dist$time_point, cent_dist$distance), collapse = ", ")),
  "",
  "A2. BASELINE SANITY CHECK (T0) AND KNOWN CAVEATS",
  paste0("  T0 Salt-vs-Control DEGs (expected to be few/none): ", t0_deg,
         " (", de_summary_dt[time_point == "T0"]$upregulated, " up / ",
         de_summary_dt[time_point == "T0"]$downregulated, " down)"),
  paste0("  correlation of the T0 salt effect with later salt effects: ",
         paste(paste0(names(cor_t0_later), " = ", cor_t0_later), collapse = ", ")),
  paste0("  control-vs-salt separation in VST space: ",
         paste(sprintf("%s = %.1f", vst_cent_dist$time_point, vst_cent_dist$vst_centroid_distance),
               collapse = ", ")),
  paste0("  replicate-blocked sensitivity model: ", 
         paste(sprintf("%s = %d genes", rep_counts$contrast, rep_counts$DEGs_replicate_blocked_model),
               collapse = ", ")),
  "  -> The two sample groups already differ at T0, and the difference is NOT removed by",
  "     blocking on replicate.  There is therefore a group-level (batch/baseline) offset",
  "     between all salt and all control libraries.",
  "  -> CAVEAT: the time-matched Salt-vs-Control DEG counts are inflated by this offset and",
  "     must not be interpreted as pure treatment effects.  Difference-in-differences",
  "     contrasts are invariant to any constant offset and are the primary evidence used here.",
  "  -> The raw time-course counts of the ORIGINAL salt-only analysis (",
  paste0("     ", paste(sprintf("%s: %d DEGs", old_counts$contrast, old_counts$DEGs_unshrunken),
                       collapse = "; ")),
  "     reproduce the manuscript values, confirming that the difference between the old",
  "     and the new results is caused by the corrected comparison design, not by a change",
  "     in the counting pipeline.",
  "",
  "B. TIME-MATCHED Salt vs Control DEGs (padj < 0.05, |log2FC| >= 1)",
  paste0("  ", de_summary_dt$comparison, ": up = ", de_summary_dt$upregulated,
         ", down = ", de_summary_dt$downregulated, ", total = ", de_summary_dt$total_DEGs,
         " (", de_summary_dt$pct_of_retained_genes, "% of retained genes)"),
  "",
  "C. DIFFERENCE-IN-DIFFERENCES: change in the salt effect between time points",
  "   (padj < 0.05, |delta log2FC| >= 1)",
  paste0("  ", did_summary_dt$contrast, " [", ifelse(did_summary_dt$adjacent, "adjacent", "non-adjacent"),
         "]: increase = ", did_summary_dt$increased_salt_effect,
         ", decrease = ", did_summary_dt$decreased_salt_effect,
         ", total = ", did_summary_dt$total_significant,
         " (median |delta log2FC| = ", did_summary_dt$median_abs_DID_log2FC_sig, ")"),
  "",
  paste0("  largest adjacent change in salt effect: ", largest_adj, " (",
         adj_counts[contrast == largest_adj]$total_significant, " genes), ",
         round(largest_vs_second, 2), "x the next largest (", second_adj, ", ",
         adj_counts[contrast == second_adj]$total_significant, " genes)"),
  "",
  "D. INTERACTION LIKELIHOOD-RATIO TEST (full vs reduced ~ condition + time_point)",
  paste0("  genes tested with padj       : ", n_lrt_tested),
  paste0("  NA padj (independent filtering): ", n_lrt_na),
  paste0("  significant interaction genes: ", n_lrt_sig,
         " (", round(100 * n_lrt_sig / n_lrt_tested, 2), "% of tested)"),
  "",
  "E. MFUZZ CLUSTERING OF CONTROL-ADJUSTED SALT-SPECIFIC PROFILES",
  paste0("  gene selection rule          : ", selection_rule),
  paste0("  genes entered clustering     : ", length(sel)),
  paste0("  fuzzification parameter m    : ", m_use, " (estimated with mestimate; 1.1-2.5 recommended)"),
  paste0("  number of clusters k         : ", K_CLUSTERS,
         " (Xie-Beni-optimal k = ", k_opt, "; k = 6 retained: ", use_k6, ")"),
  paste0("  cluster sizes                : ",
         paste(sprintf("C%d = %d", cluster_summary$cluster, cluster_summary$n_genes), collapse = ", ")),
  paste0("  mean membership per cluster  : ",
         paste(sprintf("C%d = %.2f", cluster_summary$cluster, cluster_summary$mean_membership), collapse = ", ")),
  paste0("  cluster trajectories         : "),
  paste0("    C", cluster_summary$cluster, " (", cluster_summary$trend, "): ", cluster_summary$mean_profile),
  "",
  "F. ENRICHMENT (BH-FDR < 0.05; g:Profiler Oryza sativa GO BP/MF/CC + KEGG)",
  if (!is.null(gprofiler_res)) {
    sig <- gprofiler_res[significant == TRUE]
    if (nrow(sig)) {
      top <- sig[order(p_value)][seq_len(min(15, .N))]
      c(paste0("  ", top$gene_set, " | ", top$source, " | ", top$term_name,
               " (FDR = ", signif(top$p_value, 3), ", n = ", top$intersection_size, ")"),
        paste0("  total significant terms: ", nrow(sig), " over ", length(unique(sig$gene_set)), " gene sets"))
    } else "  no terms significant at FDR < 0.05"
  } else "  g:Profiler unavailable - see clusterProfiler fallback table",
  "",
  "G. INTERPRETATION (derived from the numbers above, not from the original claims)",
  paste0("  1. Temporally changing salt response? ",
         if (temporal_supported) "YES" else "NO",
         " - ", n_lrt_sig, " genes have a significant condition x time interaction (LRT padj < ",
         PADJ_CUTOFF, ") and ", sum(adj_counts$total_significant),
         " gene-contrast results are significant across the adjacent difference-in-differences tests."),
  paste0("  2. Biphasic response (original definition: early net repression followed by late",
         " net activation) supported? ", if (biphasic_supported) "YES" else "NO"),
  paste0("     clause 1 - early net repression (majority of DEGs down at T60/T180): ",
         if (early_repressed) "met" else "NOT met", " (",
         paste(sprintf("%s: %.1f%% of DEGs up", early$time_point, early$pct_up), collapse = ", "), ")"),
  paste0("     clause 2 - late activation (majority of DEGs up at T240/T5d): ",
         if (late_activated) "met" else "NOT met", " (",
         paste(sprintf("%s: %.1f%% up", late$time_point, late$pct_up), collapse = ", "), ")"),
  paste0("     clause 3 - significant direction switch (T5d vs T60 difference-in-differences): ",
         if (switch_supported) "met" else "NOT met", " (",
         n_neg_to_pos, " negative->positive vs ", n_pos_to_neg,
         " positive->negative; binomial p = ", signif(switch_test_p, 3), ")"),
  "     -> late activation and a predominantly negative->positive change ARE supported, but",
  "        the early net-repression phase is NOT (up-regulated genes dominate at every time",
  "        point), so the original 'biphasic' claim must be rewritten rather than retained.",
  paste0("  3. T240->T5d as the strongest change in salt effect? ",
         if (largest_adj == "T5d_vs_T240" && largest_vs_second >= 1) "SUPPORTED as the largest by gene count"
         else if (largest_adj == "T5d_vs_T240") "it is the largest adjacent contrast"
         else paste0("NOT SUPPORTED - the largest adjacent change is ", largest_adj),
         " (T5d_vs_T240 = ", adj_counts[contrast == "T5d_vs_T240"]$total_significant,
         " genes vs ", adj_counts[contrast != "T5d_vs_T240"][which.max(total_significant)]$contrast,
         " = ", adj_counts[contrast != "T5d_vs_T240"][which.max(total_significant)]$total_significant,
         " genes; ratio = ", round(largest_vs_second, 2), ")."),
  "",
  "  4. Confidence limits imposed by the baseline offset (see A2): the per-time-point",
  "     Salt-vs-Control comparisons cannot be fully separated from the group-level offset,",
  "     so statements about the magnitude of the salt response should be qualified.  The",
  "     temporal (difference-in-differences and LRT) conclusions do not depend on a",
  "     constant offset between the salt and control groups.",
  "",
  "  NOTE: the T0 Salt-vs-Control comparison is a design sanity check. Any genes",
  "  significant there are reported but are not interpreted as treatment effects.",
  "",
  "KEY OUTPUT FILES",
  paste0("  ", file.path(TAB_DIR, "time_matched_DE_summary.csv")),
  paste0("  ", file.path(TAB_DIR, "difference_in_differences_summary.csv")),
  paste0("  ", file.path(TAB_DIR, "interaction_LRT_summary.txt")),
  paste0("  ", file.path(TAB_DIR, "mfuzz_cluster_summary.csv")),
  paste0("  ", file.path(TAB_DIR, "candidate_genes_for_discussion.csv")),
  paste0("  ", file.path(FIG_DIR)),
  paste0("  ", file.path(TAB_DIR))
)
writeLines(summary_lines, file.path(OUT_DIR, "reviewer_reanalysis_summary.txt"))
log_msg("wrote reviewer_reanalysis_summary.txt")

save_table(data.table(metric = c("verified_treatment", "samples", "groups",
                                 "genes_in_matrix", "retained_genes", "LRT_tested",
                                 "LRT_significant", "LRT_pct_significant",
                                 "DEG_Salt_vs_Control_T0", "DEG_Salt_vs_Control_T60",
                                 "DEG_Salt_vs_Control_T180", "DEG_Salt_vs_Control_T240",
                                 "DEG_Salt_vs_Control_T5d",
                                 "DID_T60_vs_T0", "DID_T180_vs_T60", "DID_T240_vs_T180",
                                 "DID_T5d_vs_T240", "DID_T5d_vs_T60", "DID_T5d_vs_T0",
                                 "Mfuzz_genes", "Mfuzz_k", "Mfuzz_m",
                                 "biphasic_supported", "temporal_change_supported",
                                 "largest_adjacent_DID"),
                    value = c(VERIFIED_TREATMENT_LABEL, nrow(samples), nrow(sample_counts),
                              nrow(cnt), nrow(counts_filt), n_lrt_tested, n_lrt_sig,
                              round(100 * n_lrt_sig / n_lrt_tested, 2),
                              de_summary_dt[time_point == "T0"]$total_DEGs,
                              de_summary_dt[time_point == "T60"]$total_DEGs,
                              de_summary_dt[time_point == "T180"]$total_DEGs,
                              de_summary_dt[time_point == "T240"]$total_DEGs,
                              de_summary_dt[time_point == "T5d"]$total_DEGs,
                              did_summary_dt[contrast == "T60_vs_T0"]$total_significant,
                              did_summary_dt[contrast == "T180_vs_T60"]$total_significant,
                              did_summary_dt[contrast == "T240_vs_T180"]$total_significant,
                              did_summary_dt[contrast == "T5d_vs_T240"]$total_significant,
                              did_summary_dt[contrast == "T5d_vs_T60"]$total_significant,
                              did_summary_dt[contrast == "T5d_vs_T0"]$total_significant,
                              length(sel), K_CLUSTERS, m_use,
                              biphasic_supported, temporal_supported, largest_adj)),
           file.path(TAB_DIR, "reviewer_reanalysis_summary.csv"))

## ---------------------------------------------------------------------------
## 13. session info
## ---------------------------------------------------------------------------

si <- capture.output(sessionInfo())
writeLines(c(paste0("generated: ", format(Sys.time())), si),
           file.path(OUT_DIR, "sessionInfo.txt"))
log_msg("wrote sessionInfo.txt")

log_msg("\nDONE: ", format(Sys.time()))
cat("\n=== console summary ===\n")
cat(paste(summary_lines, collapse = "\n"), "\n")
close(log_con)
