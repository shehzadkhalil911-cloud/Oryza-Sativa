#!/usr/bin/env Rscript

# Complete RNA-seq time-series DEG pipeline for Oryza sativa salt stress.
# Inputs expected in the project root:
#   - metadata.csv
#   - GSE284308_RAW.tar containing *_abundance.tsv.gz files
# Outputs are written under results/.

options(stringsAsFactors = FALSE)
set.seed(20260505)

# ------------------------------------------------------------------------------
# STEP 0: Configuration, packages, and output directories
# ------------------------------------------------------------------------------

cfg <- list(
  project_name = "rice_salt_timeseries_10mM_NaCl",
  metadata_file = "metadata.csv",
  raw_tar = "GSE284308_RAW.tar",
  extract_dir = file.path("data", "raw_abundance"),
  results_dir = "results",
  alpha = 0.05,
  lfc_cutoff = 1,
  min_count = 10,
  min_samples = 3,
  heatmap_top_n = 100,
  top_labels = 10,
  baseline_time = "T0",
  time_levels = c("T0", "T60", "T180", "T240", "T5d"),
  time_labels = c(T0 = "0 min", T60 = "60 min", T180 = "180 min",
                  T240 = "240 min", T5d = "5 days"),
  auto_install = TRUE
)

dirs <- file.path(
  cfg$results_dir,
  c("tables", "plots", "plots/pca", "plots/volcano", "plots/heatmaps",
    "plots/enrichment", "plots/overlap", "plots/qc", "reports")
)
invisible(lapply(c(cfg$extract_dir, dirs), dir.create, recursive = TRUE, showWarnings = FALSE))

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = "")))
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

cran_pkgs <- c(
  "tidyverse", "ggplot2", "pheatmap", "VennDiagram", "RColorBrewer",
  "ggpubr", "dplyr", "data.table", "ggrepel", "BiocManager"
)
bioc_pkgs_required <- c("DESeq2", "EnhancedVolcano")
bioc_pkgs_optional <- c("apeglm", "ComplexHeatmap", "circlize", "clusterProfiler", "enrichplot")
bioc_pkgs_core <- c(
  "S4Vectors", "IRanges", "BiocGenerics", "GenomeInfoDb", "XVector",
  "DelayedArray", "SummarizedExperiment", "Biobase", "GenomicRanges",
  "AnnotationDbi", "BiocParallel"
)

install_missing_packages <- function() {
  missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_cran) > 0 && isTRUE(cfg$auto_install)) {
    install.packages(missing_cran, repos = "https://cloud.r-project.org")
  }
  
  if (!requireNamespace("BiocManager", quietly = TRUE) && isTRUE(cfg$auto_install)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  
  bioc_pkgs <- unique(c(bioc_pkgs_core, bioc_pkgs_required, bioc_pkgs_optional))
  installed_names <- rownames(installed.packages())
  missing_bioc <- setdiff(bioc_pkgs, installed_names)
  s4vectors_old <- "S4Vectors" %in% installed_names &&
    package_version(as.character(utils::packageVersion("S4Vectors"))) < package_version("0.47.6")
  
  if ((length(missing_bioc) > 0 || s4vectors_old) &&
      isTRUE(cfg$auto_install) && requireNamespace("BiocManager", quietly = TRUE)) {
    BiocManager::install(unique(c(bioc_pkgs_core, missing_bioc)), ask = FALSE, update = TRUE)
  }
  
  still_missing <- c(
    cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)],
    bioc_pkgs_required[!vapply(bioc_pkgs_required, requireNamespace, logical(1), quietly = TRUE)]
  )
  if (length(still_missing) > 0) {
    stop("Required packages could not be loaded or installed: ", paste(still_missing, collapse = ", "))
  }
}

install_missing_packages()

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(ggplot2)
  library(pheatmap)
  library(EnhancedVolcano)
  library(VennDiagram)
  library(RColorBrewer)
  library(ggpubr)
  library(dplyr)
  library(data.table)
  library(ggrepel)
})

has_apeglm <- requireNamespace("apeglm", quietly = TRUE)
has_complex_heatmap <- requireNamespace("ComplexHeatmap", quietly = TRUE) &&
  requireNamespace("circlize", quietly = TRUE)
has_clusterprofiler <- requireNamespace("clusterProfiler", quietly = TRUE)
has_enrichplot <- requireNamespace("enrichplot", quietly = TRUE)

pub_theme <- function(base_size = 12) {
  theme_pubr(base_size = base_size, border = FALSE) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2),
      plot.subtitle = element_text(color = "grey25"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold")
    )
}

time_cols <- c(
  T0 = "#4B5563",
  T60 = "#0072B2",
  T180 = "#009E73",
  T240 = "#E69F00",
  T5d = "#D55E00"
)
condition_cols <- c(Control = "#6B7280", Salt = "#DC2626")
functional_cols <- c(
  Signaling = "#0072B2",
  Hormone = "#CC79A7",
  Transcription = "#6A3D9A",
  `Osmotic stress` = "#56B4E9",
  `Ion transport` = "#009E73",
  Redox = "#E69F00",
  Metabolism = "#D55E00",
  Photosynthesis = "#66A61E",
  `Cell structure` = "#8C564B",
  `Protein homeostasis` = "#1F78B4",
  Development = "#7570B3",
  KEGG = "#4B5563"
)

save_plot <- function(plot, basename, width = 8, height = 6, dpi = 1200) {
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
  ggsave(paste0(basename, ".pdf"), plot, width = width, height = height, device = pdf_device)
  ggsave(paste0(basename, ".tiff"), plot, width = width, height = height, dpi = dpi,
         compression = "lzw", device = "tiff")
}

write_csv_safe <- function(x, path) {
  data.table::fwrite(as.data.frame(x), path)
}

# ------------------------------------------------------------------------------
# STEP 1: Load data
# ------------------------------------------------------------------------------

log_msg("STEP 1: Loading metadata and abundance files")

if (!file.exists(cfg$metadata_file)) {
  stop("Missing metadata file: ", cfg$metadata_file)
}
if (!file.exists(cfg$raw_tar) && length(list.files(cfg$extract_dir, pattern = "_abundance\\.tsv\\.gz$", full.names = TRUE)) == 0) {
  stop("Missing raw abundance archive or extracted abundance files.")
}

metadata_raw <- data.table::fread(cfg$metadata_file)

if (file.exists(cfg$raw_tar)) {
  tar_files <- utils::untar(cfg$raw_tar, list = TRUE)
  abundance_members <- tar_files[grepl("_abundance\\.tsv\\.gz$", tar_files)]
  if (length(abundance_members) == 0) {
    stop("No *_abundance.tsv.gz files found inside ", cfg$raw_tar)
  }
  existing <- list.files(cfg$extract_dir, pattern = "_abundance\\.tsv\\.gz$", full.names = FALSE)
  missing_members <- abundance_members[!basename(abundance_members) %in% existing]
  if (length(missing_members) > 0) {
    log_msg("Extracting ", length(missing_members), " abundance files to ", cfg$extract_dir)
    utils::untar(cfg$raw_tar, files = missing_members, exdir = cfg$extract_dir)
  }
}

abundance_files <- list.files(cfg$extract_dir, pattern = "_abundance\\.tsv\\.gz$", full.names = TRUE)
if (length(abundance_files) == 0) {
  stop("No abundance files available after extraction.")
}

file_map <- tibble(
  abundance_file = abundance_files,
  file_basename = basename(abundance_files),
  sample_id = stringr::str_extract(file_basename, "GSM\\d+"),
  library_code = stringr::str_match(file_basename, "GSM\\d+_([^_]+)_abundance\\.tsv\\.gz$")[, 2],
  replicate = stringr::str_match(library_code, "^A(\\d+)")[, 2]
)

# ------------------------------------------------------------------------------
# STEP 2: Metadata cleaning
# ------------------------------------------------------------------------------

log_msg("STEP 2: Cleaning metadata")

metadata_clean <- metadata_raw %>%
  as_tibble() %>%
  transmute(
    run = .data[["Run"]],
    sample_id = .data[["Sample Name"]],
    genotype = .data[["genotype"]],
    organism = .data[["Organism"]],
    tissue = .data[["tissue"]],
    treatment_raw = .data[["treatment"]],
    time_raw = .data[["time"]],
    condition = case_when(
      str_detect(treatment_raw, regex("^\\s*10\\s*mM|^\\s*10mM\\s+NaCl", ignore_case = TRUE)) ~ "Salt",
      str_detect(treatment_raw, regex("^\\s*0\\s*mM|^\\s*0mM\\s+NaCl", ignore_case = TRUE)) ~ "Control",
      TRUE ~ NA_character_
    ),
    time_point = case_when(
      str_detect(time_raw, regex("^\\s*0\\s+minutes", ignore_case = TRUE)) ~ "T0",
      str_detect(time_raw, regex("^\\s*60\\s+minutes", ignore_case = TRUE)) ~ "T60",
      str_detect(time_raw, regex("^\\s*180\\s+minutes", ignore_case = TRUE)) ~ "T180",
      str_detect(time_raw, regex("^\\s*240\\s+minutes", ignore_case = TRUE)) ~ "T240",
      str_detect(time_raw, regex("5\\s+days", ignore_case = TRUE)) ~ "T5d",
      TRUE ~ NA_character_
    ),
    time_minutes = case_when(
      time_point == "T0" ~ 0,
      time_point == "T60" ~ 60,
      time_point == "T180" ~ 180,
      time_point == "T240" ~ 240,
      time_point == "T5d" ~ 5 * 24 * 60,
      TRUE ~ NA_real_
    )
  ) %>%
  left_join(file_map, by = "sample_id") %>%
  mutate(
    replicate = if_else(is.na(replicate), as.character(row_number()), replicate),
    condition = factor(condition, levels = c("Control", "Salt")),
    time_point = factor(time_point, levels = cfg$time_levels),
    time_label = factor(unname(cfg$time_labels[as.character(time_point)]), levels = unname(cfg$time_labels)),
    sample_name = paste(condition, time_point, paste0("R", replicate), sep = "_")
  ) %>%
  arrange(condition, time_point, as.integer(replicate))

if (anyNA(metadata_clean$condition) || anyNA(metadata_clean$time_point) || anyNA(metadata_clean$abundance_file)) {
  problem_rows <- metadata_clean %>%
    filter(is.na(condition) | is.na(time_point) | is.na(abundance_file))
  write_csv_safe(problem_rows, file.path(cfg$results_dir, "tables", "metadata_rows_needing_attention.csv"))
  stop("Metadata cleaning produced missing condition/time/file values. See results/tables/metadata_rows_needing_attention.csv")
}

if (anyDuplicated(metadata_clean$sample_name)) {
  metadata_clean <- metadata_clean %>%
    group_by(sample_name) %>%
    mutate(sample_name = if (n() > 1) paste0(sample_name, "_", row_number()) else sample_name) %>%
    ungroup()
}

write_csv_safe(metadata_clean, file.path(cfg$results_dir, "tables", "metadata_clean.csv"))

read_abundance_counts <- function(path, sample_name) {
  dt <- data.table::fread(path, select = c("target_id", "est_counts", "tpm"))
  dt[, gene_id := sub("-[0-9]+$", "", target_id)]
  counts_dt <- dt[, .(
    est_counts = sum(est_counts, na.rm = TRUE),
    tpm = sum(tpm, na.rm = TRUE)
  ), by = gene_id]
  setnames(counts_dt, c("est_counts", "tpm"), c(sample_name, paste0(sample_name, "_TPM")))
  counts_dt
}

log_msg("Importing and aggregating transcript-level estimated counts to gene IDs")
count_tables <- purrr::map2(metadata_clean$abundance_file, metadata_clean$sample_name, read_abundance_counts)
count_dt <- purrr::reduce(count_tables, full_join, by = "gene_id")
count_dt[is.na(count_dt)] <- 0

tpm_cols <- grep("_TPM$", names(count_dt), value = TRUE)
tpm_matrix <- count_dt %>% select(gene_id, all_of(tpm_cols))
names(tpm_matrix) <- sub("_TPM$", "", names(tpm_matrix))

count_matrix <- count_dt %>%
  select(gene_id, all_of(metadata_clean$sample_name)) %>%
  as.data.frame()
rownames(count_matrix) <- count_matrix$gene_id
count_matrix$gene_id <- NULL
count_matrix <- round(as.matrix(count_matrix))
storage.mode(count_matrix) <- "integer"

write_csv_safe(
  tibble(gene_id = rownames(count_matrix)) %>% bind_cols(as.data.frame(count_matrix)),
  file.path(cfg$results_dir, "tables", "clean_count_matrix_gene_level.csv")
)
write_csv_safe(tpm_matrix, file.path(cfg$results_dir, "tables", "clean_tpm_matrix_gene_level.csv"))

# ------------------------------------------------------------------------------
# STEP 3: Filter low-expression genes
# ------------------------------------------------------------------------------

log_msg("STEP 3: Filtering low-expression genes")

keep <- rowSums(count_matrix >= cfg$min_count) >= cfg$min_samples
count_matrix_filtered <- count_matrix[keep, , drop = FALSE]

filter_summary <- tibble(
  total_genes = nrow(count_matrix),
  retained_genes = nrow(count_matrix_filtered),
  removed_low_expression = nrow(count_matrix) - nrow(count_matrix_filtered),
  min_count = cfg$min_count,
  min_samples = cfg$min_samples
)
write_csv_safe(filter_summary, file.path(cfg$results_dir, "tables", "filtering_summary.csv"))

metadata_deseq <- metadata_clean %>%
  select(sample_name, run, sample_id, condition, time_point, time_label, time_minutes,
         replicate, genotype, organism, tissue, treatment_raw, time_raw) %>%
  as.data.frame()
rownames(metadata_deseq) <- metadata_deseq$sample_name
metadata_deseq <- metadata_deseq[colnames(count_matrix_filtered), , drop = FALSE]

# ------------------------------------------------------------------------------
# STEP 4: DESeq2 object creation
# ------------------------------------------------------------------------------

log_msg("STEP 4: Creating DESeq2 objects")

dds_all <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = metadata_deseq,
  design = ~ condition + time_point
)
dds_all$condition <- relevel(dds_all$condition, ref = "Control")
dds_all$time_point <- relevel(dds_all$time_point, ref = cfg$baseline_time)

salt_samples <- rownames(metadata_deseq)[metadata_deseq$condition == "Salt"]
dds_salt <- dds_all[, salt_samples]
design(dds_salt) <- ~ time_point
dds_salt$time_point <- relevel(droplevels(dds_salt$time_point), ref = cfg$baseline_time)

# ------------------------------------------------------------------------------
# STEP 5: Run DESeq2
# ------------------------------------------------------------------------------

log_msg("STEP 5: Running DESeq2")

dds_all <- DESeq(dds_all)
dds_salt <- DESeq(dds_salt)

norm_counts <- counts(dds_all, normalized = TRUE)
write_csv_safe(
  tibble(gene_id = rownames(norm_counts)) %>% bind_cols(as.data.frame(norm_counts)),
  file.path(cfg$results_dir, "tables", "deseq2_normalized_counts.csv")
)

vsd_all <- vst(dds_all, blind = FALSE)
vst_matrix <- assay(vsd_all)
salt_vst <- vst_matrix[, salt_samples]
write_csv_safe(
  tibble(gene_id = rownames(vst_matrix)) %>% bind_cols(as.data.frame(vst_matrix)),
  file.path(cfg$results_dir, "tables", "deseq2_vst_expression.csv")
)

rld_matrix <- NULL
if (ncol(dds_all) <= 50) {
  rld_matrix <- tryCatch({
    log_msg("Computing rlog transformation")
    assay(rlog(dds_all, blind = FALSE))
  }, error = function(e) {
    log_msg("rlog skipped: ", conditionMessage(e))
    NULL
  })
  if (!is.null(rld_matrix)) {
    write_csv_safe(
      tibble(gene_id = rownames(rld_matrix)) %>% bind_cols(as.data.frame(rld_matrix)),
      file.path(cfg$results_dir, "tables", "deseq2_rlog_expression.csv")
    )
  }
}

saveRDS(dds_all, file.path(cfg$results_dir, "tables", "dds_all_samples.rds"))
saveRDS(dds_salt, file.path(cfg$results_dir, "tables", "dds_salt_timecourse.rds"))
saveRDS(vsd_all, file.path(cfg$results_dir, "tables", "vst_all_samples.rds"))

# ------------------------------------------------------------------------------
# Publication QC plots
# ------------------------------------------------------------------------------

log_msg("Generating QC plots")

lib_sizes <- metadata_clean %>%
  mutate(total_counts = colSums(count_matrix)[sample_name]) %>%
  arrange(condition, time_point, replicate)

p_lib <- ggplot(lib_sizes, aes(x = time_label, y = total_counts / 1e6, fill = condition)) +
  geom_col(position = position_dodge(width = 0.75), color = "grey25", linewidth = 0.15, width = 0.65) +
  scale_fill_manual(values = condition_cols) +
  labs(
    title = "RNA-seq library sizes",
    x = NULL,
    y = "Total estimated counts (millions)",
    fill = "Condition"
  ) +
  pub_theme() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_plot(p_lib, file.path(cfg$results_dir, "plots", "qc", "library_size_qc"), width = 7, height = 4.5)

sample_cor <- cor(vst_matrix, method = "pearson")
ann_col <- metadata_deseq[, c("condition", "time_point", "replicate"), drop = FALSE]
ann_col$time_point <- factor(ann_col$time_point, levels = cfg$time_levels)
annotation_colors <- list(
  condition = condition_cols,
  time_point = time_cols
)
pheatmap(
  sample_cor,
  annotation_col = ann_col,
  annotation_row = ann_col,
  annotation_colors = annotation_colors,
  color = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),
  border_color = NA,
  fontsize = 8,
  main = "Sample correlation using VST expression",
  filename = file.path(cfg$results_dir, "plots", "qc", "sample_correlation_heatmap.pdf"),
  width = 8,
  height = 7.5
)
tiff(file.path(cfg$results_dir, "plots", "qc", "sample_correlation_heatmap.tiff"),
    width = 8, height = 7.5, units = "in", res = 1200, compression = "lzw")
pheatmap(
  sample_cor,
  annotation_col = ann_col,
  annotation_row = ann_col,
  annotation_colors = annotation_colors,
  color = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),
  border_color = NA,
  fontsize = 8,
  main = "Sample correlation using VST expression"
)
dev.off()

# ------------------------------------------------------------------------------
# STEP 6: Apply shrinkage and STEP 7: Extract all contrasts
# ------------------------------------------------------------------------------

log_msg("STEP 6-7: Extracting shrunken contrasts for salt-treated time course")

contrast_table <- tribble(
  ~contrast_id, ~time_point, ~coef_name, ~label,
  "T60_vs_T0", "T60", "time_point_T60_vs_T0", "60 min vs 0 min",
  "T180_vs_T0", "T180", "time_point_T180_vs_T0", "180 min vs 0 min",
  "T240_vs_T0", "T240", "time_point_T240_vs_T0", "240 min vs 0 min",
  "T5d_vs_T0", "T5d", "time_point_T5d_vs_T0", "5 days vs 0 min"
)

available_results_names <- resultsNames(dds_salt)
missing_coef <- setdiff(contrast_table$coef_name, available_results_names)
if (length(missing_coef) > 0) {
  stop("Expected DESeq2 coefficient(s) not found: ", paste(missing_coef, collapse = ", "),
       "\nAvailable: ", paste(available_results_names, collapse = ", "))
}

extract_shrunken_result <- function(dds, contrast_id, time_point, coef_name, label) {
  shrink_type <- if (has_apeglm) "apeglm" else "normal"
  res <- lfcShrink(dds, coef = coef_name, type = shrink_type)
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    as_tibble() %>%
    mutate(
      contrast_id = contrast_id,
      comparison = label,
      time_point = time_point,
      shrinkage = shrink_type,
      baseMean = if_else(is.na(baseMean), 0, baseMean),
      log2FoldChange = if_else(is.na(log2FoldChange), 0, log2FoldChange),
      padj = if_else(is.na(padj), 1, padj),
      pvalue = if_else(is.na(pvalue), 1, pvalue),
      regulation = case_when(
        padj < cfg$alpha & log2FoldChange >= cfg$lfc_cutoff ~ "Up",
        padj < cfg$alpha & log2FoldChange <= -cfg$lfc_cutoff ~ "Down",
        TRUE ~ "Not significant"
      ),
      significant = regulation != "Not significant"
    ) %>%
    arrange(padj, desc(abs(log2FoldChange)))
  write_csv_safe(res_df, file.path(cfg$results_dir, "tables", paste0("DESeq2_", contrast_id, "_all_genes.csv")))
  write_csv_safe(
    filter(res_df, significant),
    file.path(cfg$results_dir, "tables", paste0("DESeq2_", contrast_id, "_significant_DEGs.csv"))
  )
  res_df
}

deg_results <- pmap_dfr(contrast_table, extract_shrunken_result, dds = dds_salt)
write_csv_safe(deg_results, file.path(cfg$results_dir, "tables", "DESeq2_all_timepoint_contrasts.csv"))

deg_summary <- deg_results %>%
  count(comparison, time_point, regulation) %>%
  group_by(comparison, time_point) %>%
  complete(regulation = c("Up", "Down", "Not significant"), fill = list(n = 0)) %>%
  ungroup() %>%
  mutate(time_point = factor(time_point, levels = cfg$time_levels))
write_csv_safe(deg_summary, file.path(cfg$results_dir, "tables", "DEG_summary_by_timepoint.csv"))

p_deg_summary <- deg_summary %>%
  filter(regulation != "Not significant") %>%
  mutate(direction_n = if_else(regulation == "Down", -n, n)) %>%
  ggplot(aes(x = factor(time_point, levels = cfg$time_levels), y = direction_n, fill = regulation)) +
  geom_col(width = 0.65, color = "grey20", linewidth = 0.15) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  scale_fill_manual(values = c(Up = "#DC2626", Down = "#2563EB")) +
  scale_x_discrete(labels = cfg$time_labels) +
  labs(
    title = "Differential expression across the salt-stress time course",
    x = NULL,
    y = "Number of significant DEGs",
    fill = "Regulation"
  ) +
  pub_theme() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_plot(p_deg_summary, file.path(cfg$results_dir, "plots", "qc", "deg_summary_barplot"), width = 7, height = 4.8)

# ------------------------------------------------------------------------------
# PCA plot using VST expression
# ------------------------------------------------------------------------------

log_msg("Generating publication-quality PCA plot")

pca <- prcomp(t(vst_matrix), center = TRUE, scale. = FALSE)
pca_percent <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
pca_df <- as_tibble(pca$x[, 1:4], rownames = "sample_name") %>%
  left_join(metadata_clean, by = "sample_name")

p_pca <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = time_point, shape = condition)
) +
  stat_ellipse(
    aes(group = time_point, color = time_point),
    type = "norm",
    level = 0.68,
    linewidth = 0.35,
    alpha = 0.45,
    show.legend = FALSE
  ) +
  geom_point(size = 3.4, alpha = 0.95) +
  geom_text_repel(
    aes(label = paste0(condition, "_", cfg$time_labels[as.character(time_point)], "_R", replicate)),
    size = 2.5,
    max.overlaps = 25,
    show.legend = FALSE,
    box.padding = 0.3,
    min.segment.length = 0
  ) +
  scale_color_manual(values = time_cols, labels = cfg$time_labels, name = "Time point") +
  scale_shape_manual(values = c(Control = 16, Salt = 17), name = "Condition") +
  labs(
    title = "PCA of VST-normalized rice leaf transcriptomes",
    subtitle = "Oryza sativa under control and 10 mM NaCl time-series treatment",
    x = paste0("PC1 (", pca_percent[1], "%)"),
    y = paste0("PC2 (", pca_percent[2], "%)")
  ) +
  pub_theme(base_size = 12) +
  guides(color = guide_legend(override.aes = list(size = 4)))

save_plot(p_pca, file.path(cfg$results_dir, "plots", "pca", "PCA_VST_time_condition"), width = 8.2, height = 6.2)
write_csv_safe(pca_df, file.path(cfg$results_dir, "tables", "PCA_coordinates_VST.csv"))

# ------------------------------------------------------------------------------
# Enhanced volcano plots for every time point
# ------------------------------------------------------------------------------

log_msg("Generating EnhancedVolcano plots")

make_volcano <- function(res_df, contrast_id, label) {
  top_genes <- res_df %>%
    filter(padj < cfg$alpha, abs(log2FoldChange) > cfg$lfc_cutoff) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = cfg$top_labels) %>%
    pull(gene_id)
  
  plot_obj <- EnhancedVolcano::EnhancedVolcano(
    as.data.frame(res_df),
    lab = res_df$gene_id,
    x = "log2FoldChange",
    y = "padj",
    selectLab = top_genes,
    xlab = expression(Log[2] ~ "fold change"),
    ylab = expression("-" ~ Log[10] ~ "adjusted P"),
    title = paste0(label, " in salt-treated rice"),
    subtitle = paste0("DESeq2 with ", unique(res_df$shrinkage), " shrinkage; padj < ",
                      cfg$alpha, ", |log2FC| > ", cfg$lfc_cutoff),
    pCutoff = cfg$alpha,
    FCcutoff = cfg$lfc_cutoff,
    pointSize = 1.7,
    labSize = 3.2,
    boxedLabels = FALSE,
    drawConnectors = TRUE,
    widthConnectors = 0.45,
    colConnectors = "grey35",
    max.overlaps = Inf,
    col = c("grey72", "#6B7280", "#2563EB", "#DC2626"),
    legendLabels = c("NS", "Log2FC", "padj", "padj and Log2FC"),
    legendPosition = "right",
    legendLabSize = 10,
    legendIconSize = 3,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    border = "partial",
    borderWidth = 0.4,
    borderColour = "grey35"
  ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey25")
    )
  
  save_plot(plot_obj, file.path(cfg$results_dir, "plots", "volcano", paste0("EnhancedVolcano_", contrast_id)), width = 8, height = 6.5)
  invisible(plot_obj)
}

pwalk(contrast_table, function(contrast_id, time_point, coef_name, label) {
  make_volcano(filter(deg_results, .data$contrast_id == contrast_id), contrast_id, label)
})




#!/usr/bin/env Rscript

# ==============================================================================
#  SECTION 1 : RNA-Seq Preprocessing Quality Control & PCA
#  Project    : rice_salt_timeseries_10mM_NaCl  (GSE284308)
#  Organism   : Oryza sativa
#  Design     : Control vs 10 mM NaCl, five time points (T0–T5d)
#
#  Requires   : outputs from my_new_analysis.R (Steps 0-7)
#                 results/tables/dds_all_samples.rds
#                 results/tables/vst_all_samples.rds
#                 results/tables/metadata_clean.csv
#                 results/tables/clean_count_matrix_gene_level.csv
#                 results/tables/deseq2_normalized_counts.csv
#                 results/tables/DESeq2_all_timepoint_contrasts.csv
#
#  Outputs    : results/plots/qc/          – all QC figures (PDF + PNG)
#               results/plots/pca/         – PCA figures   (PDF + PNG)
#               results/tables/            – QC metric CSVs
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(20260505)

# ── Re-use the same cfg block so paths always agree ──────────────────────────
cfg <- list(
  project_name  = "rice_salt_timeseries_10mM_NaCl",
  metadata_file = "metadata.csv",
  results_dir   = "results",
  alpha         = 0.05,
  lfc_cutoff    = 1,
  min_count     = 10,
  min_samples   = 3,
  baseline_time = "T0",
  time_levels   = c("T0", "T60", "T180", "T240", "T5d"),
  time_labels   = c(T0 = "0 min", T60 = "60 min", T180 = "180 min",
                    T240 = "240 min", T5d = "5 days"),
  auto_install  = TRUE
)

# ── Output subdirectories (same as parent script) ────────────────────────────
dirs <- file.path(cfg$results_dir,
                  c("tables", "plots/qc", "plots/pca"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# ── Logging helper (matches parent script) ───────────────────────────────────
log_msg <- function(...) {
  cat(sprintf("[%s] %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              paste0(..., collapse = "")))
}

# ── Package installation (same logic as parent script) ───────────────────────
cran_pkgs <- c("tidyverse", "ggplot2", "pheatmap", "RColorBrewer",
               "ggpubr", "dplyr", "data.table", "ggrepel",
               "patchwork", "scales", "FactoMineR", "factoextra",
               "BiocManager")
bioc_pkgs <- c("DESeq2", "ComplexHeatmap", "circlize", "BiocParallel")

if (isTRUE(cfg$auto_install)) {
  missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace,
                                    logical(1), quietly = TRUE)]
  if (length(missing_cran) > 0)
    install.packages(missing_cran, repos = "https://cloud.r-project.org")
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  
  missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace,
                                    logical(1), quietly = TRUE)]
  if (length(missing_bioc) > 0)
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
  library(ggpubr)
  library(data.table)
  library(scales)
  library(FactoMineR)
  library(factoextra)
})
has_complex_heatmap <- requireNamespace("ComplexHeatmap", quietly = TRUE) &&
  requireNamespace("circlize", quietly = TRUE)
if (has_complex_heatmap) {
  suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize) })
}

# ── Colour palettes (identical to parent script) ─────────────────────────────
time_cols <- c(
  T0   = "#4B5563",
  T60  = "#0072B2",
  T180 = "#009E73",
  T240 = "#E69F00",
  T5d  = "#D55E00"
)
condition_cols <- c(Control = "#6B7280", Salt = "#DC2626")

# ── Publication theme (identical to parent script) ───────────────────────────
pub_theme <- function(base_size = 12) {
  theme_pubr(base_size = base_size, border = FALSE) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0, size = base_size + 2),
      plot.subtitle = element_text(color = "grey25"),
      axis.title    = element_text(face = "bold"),
      legend.title  = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text    = element_text(face = "bold")
    )
}

# ── save_plot / write_csv_safe helpers (identical to parent script) ───────────
save_plot <- function(plot, basename, width = 8, height = 6, dpi = 1200) {
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
  ggsave(paste0(basename, ".pdf"), plot, width = width, height = height, device = pdf_device)
  ggsave(paste0(basename, ".tiff"), plot, width = width, height = height, dpi = dpi,
         compression = "lzw", device = "tiff")
}

write_csv_safe <- function(x, path) {
  data.table::fwrite(as.data.frame(x), path)
}


# ==============================================================================
#  LOAD OBJECTS produced by my_new_analysis.R
# ==============================================================================

log_msg("Loading RDS and CSV objects from my_new_analysis.R")

rds_dds <- file.path(cfg$results_dir, "tables", "dds_all_samples.rds")
rds_vsd <- file.path(cfg$results_dir, "tables", "vst_all_samples.rds")

if (!file.exists(rds_dds) || !file.exists(rds_vsd)) {
  stop("Run my_new_analysis.R first to generate RDS objects in results/tables/")
}

dds_all  <- readRDS(rds_dds)
vsd_all  <- readRDS(rds_vsd)
vst_matrix <- assay(vsd_all)

metadata_clean <- data.table::fread(
  file.path(cfg$results_dir, "tables", "metadata_clean.csv")
) %>%
  as_tibble() %>%
  mutate(
    condition  = factor(condition,  levels = c("Control", "Salt")),
    time_point = factor(time_point, levels = cfg$time_levels),
    time_label = factor(unname(cfg$time_labels[as.character(time_point)]),
                         levels = unname(cfg$time_labels))
  )

salt_samples <- metadata_clean %>%
  filter(condition == "Salt") %>%
  pull(sample_name)
salt_vst <- vst_matrix[, salt_samples]

count_matrix_raw <- data.table::fread(
  file.path(cfg$results_dir, "tables", "clean_count_matrix_gene_level.csv")
) %>% as.data.frame() %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()

norm_counts <- data.table::fread(
  file.path(cfg$results_dir, "tables", "deseq2_normalized_counts.csv")
) %>% as.data.frame() %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()

deg_results <- data.table::fread(
  file.path(cfg$results_dir, "tables", "DESeq2_all_timepoint_contrasts.csv")
) %>% as_tibble() %>%
  mutate(time_point = factor(time_point, levels = cfg$time_levels))

log_msg(sprintf("  Loaded: %d genes × %d samples", nrow(vst_matrix), ncol(vst_matrix)))


# ==============================================================================
#  QC-1 : Per-Sample Library Sizes  (raw counts)
# ==============================================================================

log_msg("QC-1: Library size distribution")

lib_size_df <- metadata_clean %>%
  mutate(
    raw_counts       = colSums(count_matrix_raw)[sample_name],
    norm_counts_sum  = colSums(norm_counts)[sample_name]
  ) %>%
  arrange(condition, time_point, replicate)

write_csv_safe(lib_size_df,
               file.path(cfg$results_dir, "tables", "qc_library_sizes.csv"))

# ── Grouped bar: raw counts ───────────────────────────────────────────────────
p_lib_raw <- ggplot(lib_size_df,
                    aes(x = time_label, y = raw_counts / 1e6,
                        fill = condition,
                        group = interaction(time_label, condition))) +
  geom_col(position = position_dodge(width = 0.75),
           colour = "grey25", linewidth = 0.15, width = 0.65) +
  geom_hline(yintercept = 10, linetype = "dashed",
             colour = "firebrick", linewidth = 0.4) +
  annotate("text", x = 0.65, y = 10.4,
           label = "10 M threshold", size = 2.8, colour = "firebrick", hjust = 0) +
  scale_fill_manual(values = condition_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title    = "RNA-seq library sizes – raw estimated counts",
    subtitle = "Dashed line: 10 million reads (quality threshold)",
    x        = NULL,
    y        = "Total estimated counts (millions)",
    fill     = "Condition"
  ) +
  pub_theme() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_plot(p_lib_raw,
          file.path(cfg$results_dir, "plots", "qc", "qc1_library_sizes_raw"),
          width = 7.5, height = 4.8)

# ── Box plot: per-gene count distribution (log10) ────────────────────────────
count_long <- as.data.frame(log10(norm_counts + 1)) %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = "sample_name", values_to = "log10_norm") %>%
  left_join(metadata_clean %>% select(sample_name, condition, time_point,
                                      time_label, replicate),
            by = "sample_name")

p_count_box <- ggplot(count_long,
                      aes(x = sample_name, y = log10_norm,
                          fill = condition, colour = time_point)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.35, alpha = 0.7) +
  scale_fill_manual(values = condition_cols) +
  scale_colour_manual(values = time_cols) +
  scale_x_discrete(labels = function(x) {
    meta <- metadata_clean[match(x, metadata_clean$sample_name), ]
    paste0(meta$condition, "\n", meta$time_label, "\nR", meta$replicate)
  }) +
  labs(
    title    = "Normalised count distribution across samples",
    subtitle = "DESeq2 size-factor normalised counts, log₁₀(count + 1)",
    x        = NULL,
    y        = expression(log[10]~(normalised~count + 1)),
    fill     = "Condition",
    colour   = "Time point"
  ) +
  pub_theme(base_size = 10) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "right")

save_plot(p_count_box,
          file.path(cfg$results_dir, "plots", "qc", "qc1_count_distribution_boxplot"),
          width = 12, height = 5.5)

log_msg("  ✔ Library size plots saved")


# ==============================================================================
#  QC-2 : Dispersion & Size Factor diagnostics
# ==============================================================================

log_msg("QC-2: Dispersion & size factor plots")

# ── Size factors ─────────────────────────────────────────────────────────────
sf_df <- tibble(
  sample_name  = colnames(dds_all),
  size_factor  = sizeFactors(dds_all)
) %>%
  left_join(metadata_clean %>% select(sample_name, condition,
                                      time_point, time_label, replicate),
            by = "sample_name")

write_csv_safe(sf_df,
               file.path(cfg$results_dir, "tables", "qc_size_factors.csv"))

p_sf <- ggplot(sf_df, aes(x = sample_name, y = size_factor,
                          fill = condition, colour = time_point)) +
  geom_col(linewidth = 0.3, alpha = 0.85, width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  scale_fill_manual(values = condition_cols) +
  scale_colour_manual(values = time_cols) +
  scale_x_discrete(labels = function(x) {
    meta <- metadata_clean[match(x, metadata_clean$sample_name), ]
    paste0(meta$condition, "\n", meta$time_label, "\nR", meta$replicate)
  }) +
  labs(
    title    = "DESeq2 size factors per sample",
    subtitle = "Values near 1 indicate uniform sequencing depth",
    x        = NULL,
    y        = "Size factor",
    fill     = "Condition",
    colour   = "Time point"
  ) +
  pub_theme(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

save_plot(p_sf,
          file.path(cfg$results_dir, "plots", "qc", "qc2_size_factors"),
          width = 11, height = 5)

# ── Dispersion plot ───────────────────────────────────────────────────────────
pdf(file.path(cfg$results_dir, "plots", "qc", "qc2_dispersion_estimates.pdf"),
    width = 6.5, height = 5.5)
plotDispEsts(dds_all,
             main = "DESeq2 dispersion estimates\nOryza sativa – salt stress time course",
             cex = 0.4, cex.lab = 1.1)
dev.off()
tiff(file.path(cfg$results_dir, "plots", "qc", "qc2_dispersion_estimates.tiff"),
    width = 6.5, height = 5.5, units = "in", res = 1200, compression = "lzw")
plotDispEsts(dds_all,
             main = "DESeq2 dispersion estimates\nOryza sativa – salt stress time course",
             cex = 0.4, cex.lab = 1.1)
dev.off()

log_msg("  ✔ Dispersion and size factor plots saved")


# ==============================================================================
#  QC-3 : Sample-to-Sample Correlation Heatmap
# ==============================================================================

log_msg("QC-3: Sample correlation heatmap")

sample_cor <- cor(vst_matrix, method = "pearson")
write_csv_safe(
  as.data.frame(sample_cor) %>% rownames_to_column("sample"),
  file.path(cfg$results_dir, "tables", "qc_sample_pearson_correlation.csv")
)

ann_df <- metadata_clean %>%
  select(sample_name, condition, time_point, replicate) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("sample_name")
ann_df$replicate <- factor(ann_df$replicate)

annotation_colors <- list(
  condition  = condition_cols,
  time_point = time_cols,
  replicate  = setNames(
    RColorBrewer::brewer.pal(max(3, length(unique(ann_df$replicate))), "Set2"),
    levels(ann_df$replicate)
  )
)

# ── pheatmap version (always available) ──────────────────────────────────────
pheatmap(
  sample_cor,
  annotation_col    = ann_df,
  annotation_row    = ann_df,
  annotation_colors = annotation_colors,
  color             = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),
  border_color      = NA,
  fontsize          = 8,
  display_numbers   = TRUE,
  number_format     = "%.2f",
  number_color      = "grey20",
  fontsize_number   = 5.5,
  clustering_method = "ward.D2",
  main              = "Sample Pearson correlation — VST expression\nOryza sativa 10 mM NaCl time course",
  filename          = file.path(cfg$results_dir, "plots", "qc",
                                "qc3_sample_correlation_heatmap.pdf"),
  width = 9, height = 8
)

# ── TIFF version via pheatmap ──────────────────────────────────────────────────
tiff(file.path(cfg$results_dir, "plots", "qc",
               "qc3_sample_correlation_heatmap.tiff"),
    width = 9, height = 8, units = "in", res = 1200, compression = "lzw")
pheatmap(
  sample_cor,
  annotation_col    = ann_df,
  annotation_row    = ann_df,
  annotation_colors = annotation_colors,
  color             = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),
  border_color      = NA,
  fontsize          = 8,
  display_numbers   = TRUE,
  number_format     = "%.2f",
  number_color      = "grey20",
  fontsize_number   = 5.5,
  clustering_method = "ward.D2",
  main              = "Sample Pearson correlation — VST expression"
)
dev.off()

# ── Spearman distance heatmap as complement ───────────────────────────────────
sample_spearman <- cor(vst_matrix, method = "spearman")
pheatmap(
  sample_spearman,
  annotation_col    = ann_df,
  annotation_row    = ann_df,
  annotation_colors = annotation_colors,
  color             = colorRampPalette(rev(brewer.pal(9, "PuOr")))(100),
  border_color      = NA,
  fontsize          = 8,
  clustering_method = "ward.D2",
  main              = "Sample Spearman correlation — VST expression",
  filename          = file.path(cfg$results_dir, "plots", "qc",
                                "qc3_sample_spearman_heatmap.pdf"),
  width = 9, height = 8
)
tiff(file.path(cfg$results_dir, "plots", "qc",
              "qc3_sample_spearman_heatmap.tiff"),
    width = 9, height = 8, units = "in", res = 1200, compression = "lzw")
pheatmap(
  sample_spearman,
  annotation_col    = ann_df,
  annotation_row    = ann_df,
  annotation_colors = annotation_colors,
  color             = colorRampPalette(rev(brewer.pal(9, "PuOr")))(100),
  border_color      = NA,
  fontsize          = 8,
  clustering_method = "ward.D2",
  main              = "Sample Spearman correlation — VST expression"
)
dev.off()

log_msg("  ✔ Sample correlation heatmaps saved")


# ==============================================================================
#  QC-4 : Mean–Variance & Cook's Distance
# ==============================================================================

log_msg("QC-4: Mean–variance trend and Cook's distance")

# ── Mean–variance (log-normalised) ───────────────────────────────────────────
mv_df <- tibble(
  mean_expr = rowMeans(log2(norm_counts + 1)),
  var_expr  = apply(log2(norm_counts + 1), 1, var)
)

p_mv <- ggplot(mv_df, aes(x = mean_expr, y = var_expr)) +
  geom_hex(bins = 80) +
  geom_smooth(method = "loess", se = FALSE,
              colour = "firebrick", linewidth = 0.7) +
  scale_fill_viridis_c(option = "plasma", trans = "log10",
                       name = "Gene count") +
  labs(
    title    = "Mean–variance relationship (normalised counts)",
    subtitle = "Hex bin density; loess trend in red",
    x        = expression(log[2]~(mean~normalised~count + 1)),
    y        = expression(Variance~(log[2]~counts))
  ) +
  pub_theme()

save_plot(p_mv,
          file.path(cfg$results_dir, "plots", "qc", "qc4_mean_variance"),
          width = 6.5, height = 5)

# ── Cook's distance outlier check ────────────────────────────────────────────
cooks_max <- apply(assays(dds_all)[["cooks"]], 2, max, na.rm = TRUE)
cooks_df  <- tibble(
  sample_name = names(cooks_max),
  max_cooks   = cooks_max
) %>%
  left_join(metadata_clean %>% select(sample_name, condition,
                                      time_point, time_label, replicate),
            by = "sample_name")

write_csv_safe(cooks_df,
               file.path(cfg$results_dir, "tables", "qc_cooks_distance.csv"))

p_cooks <- ggplot(cooks_df,
                  aes(x = sample_name, y = max_cooks,
                      fill = condition, colour = time_point)) +
  geom_col(linewidth = 0.3, alpha = 0.85, width = 0.7) +
  geom_hline(yintercept = qf(0.99, p = ncol(dds_all), df2 = ncol(dds_all) - ncol(dds_all)),
             linetype = "dashed", colour = "firebrick", linewidth = 0.4) +
  scale_fill_manual(values = condition_cols) +
  scale_colour_manual(values = time_cols) +
  scale_x_discrete(labels = function(x) {
    meta <- metadata_clean[match(x, metadata_clean$sample_name), ]
    paste0(meta$condition, "\n", meta$time_label, "\nR", meta$replicate)
  }) +
  labs(
    title    = "Maximum Cook's distance per sample",
    subtitle = "Potential outlier samples have unusually high Cook's D",
    x        = NULL,
    y        = "Maximum Cook's distance",
    fill     = "Condition", colour = "Time point"
  ) +
  pub_theme(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

save_plot(p_cooks,
          file.path(cfg$results_dir, "plots", "qc", "qc4_cooks_distance"),
          width = 11, height = 5)

log_msg("  ✔ Mean–variance and Cook's distance plots saved")


# ==============================================================================
#  PCA-1 : Full PCA — all samples, all timepoints
# ==============================================================================

log_msg("PCA-1: Full PCA with all samples")

pca_obj     <- prcomp(t(vst_matrix), center = TRUE, scale. = FALSE)
pca_percent <- round(100 * (pca_obj$sdev^2 / sum(pca_obj$sdev^2)), 2)

pca_df <- as_tibble(pca_obj$x[, 1:6], rownames = "sample_name") %>%
  left_join(metadata_clean, by = "sample_name")

write_csv_safe(pca_df,
               file.path(cfg$results_dir, "tables", "PCA_coordinates_VST.csv"))

# ── PC1 vs PC2 ───────────────────────────────────────────────────────────────
p_pca12 <- ggplot(pca_df,
                  aes(x = PC1, y = PC2,
                      colour = time_point, shape = condition)) +
  stat_ellipse(aes(group = time_point, colour = time_point),
               type = "norm", level = 0.68,
               linewidth = 0.4, alpha = 0.5, show.legend = FALSE) +
  geom_point(size = 3.8, alpha = 0.95) +
  geom_text_repel(
    aes(label = paste0(condition, "_",
                       cfg$time_labels[as.character(time_point)],
                       "_R", replicate)),
    size = 2.5, max.overlaps = 30, show.legend = FALSE,
    box.padding = 0.35, min.segment.length = 0
  ) +
  scale_colour_manual(values = time_cols,
                      labels = cfg$time_labels,
                      name   = "Time point") +
  scale_shape_manual(values = c(Control = 16, Salt = 17),
                     name   = "Condition") +
  labs(
    title    = "PCA of VST-normalised rice leaf transcriptomes",
    subtitle = "Oryza sativa — control and 10 mM NaCl salt stress time course",
    x        = sprintf("PC1 (%.1f%%)", pca_percent[1]),
    y        = sprintf("PC2 (%.1f%%)", pca_percent[2])
  ) +
  pub_theme(base_size = 12) +
  guides(colour = guide_legend(override.aes = list(size = 4)))

save_plot(p_pca12,
          file.path(cfg$results_dir, "plots", "pca", "pca1_PC1vsPC2"),
          width = 8.5, height = 6.5)

# ── PC2 vs PC3 (stress progression axis) ────────────────────────────────────
p_pca23 <- ggplot(pca_df,
                  aes(x = PC2, y = PC3,
                      colour = time_point, shape = condition)) +
  geom_point(size = 3.8, alpha = 0.95) +
  geom_text_repel(
    aes(label = paste0(condition, "_R", replicate)),
    size = 2.5, max.overlaps = 25, show.legend = FALSE,
    box.padding = 0.3
  ) +
  scale_colour_manual(values = time_cols,
                      labels = cfg$time_labels,
                      name   = "Time point") +
  scale_shape_manual(values = c(Control = 16, Salt = 17),
                     name   = "Condition") +
  labs(
    title = "PCA — PC2 vs PC3 (stress progression axis)",
    x     = sprintf("PC2 (%.1f%%)", pca_percent[2]),
    y     = sprintf("PC3 (%.1f%%)", pca_percent[3])
  ) +
  pub_theme(base_size = 12)

save_plot(p_pca23,
          file.path(cfg$results_dir, "plots", "pca", "pca1_PC2vsPC3"),
          width = 7.5, height = 6)

log_msg("  ✔ Full PCA plots saved")


# ==============================================================================
#  PCA-2 : Scree Plot — variance explained & cumulative
# ==============================================================================

log_msg("PCA-2: Scree plot")

n_pcs    <- min(20, length(pca_percent))
scree_df <- tibble(
  PC         = 1:n_pcs,
  Variance   = pca_percent[1:n_pcs],
  Cumulative = cumsum(pca_percent[1:n_pcs])
)

write_csv_safe(scree_df,
               file.path(cfg$results_dir, "tables", "PCA_variance_explained.csv"))

p_scree <- ggplot(scree_df, aes(x = PC)) +
  geom_bar(aes(y = Variance), stat = "identity",
           fill = "#0072B2", alpha = 0.8, width = 0.7) +
  geom_line(aes(y = Cumulative), colour = "#D55E00",
            linewidth = 0.9, linetype = "solid") +
  geom_point(aes(y = Cumulative), colour = "#D55E00", size = 2.8) +
  geom_hline(yintercept = 80, linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  annotate("text", x = n_pcs * 0.75, y = 81.5,
           label = "80% cumulative variance", size = 3, colour = "grey30") +
  scale_x_continuous(breaks = 1:n_pcs) +
  scale_y_continuous(
    name     = "Variance explained (%)",
    sec.axis = sec_axis(~ ., name = "Cumulative variance (%)")
  ) +
  labs(
    title    = "PCA scree plot — variance explained per PC",
    subtitle = "Bars: individual variance; orange line: cumulative variance",
    x        = "Principal component"
  ) +
  pub_theme(base_size = 12)

save_plot(p_scree,
          file.path(cfg$results_dir, "plots", "pca", "pca2_scree_variance"),
          width = 7.5, height = 5)

log_msg("  ✔ Scree plot saved")


# ==============================================================================
#  PCA-3 : Replicate Consistency — inter-replicate distance per group
# ==============================================================================

log_msg("PCA-3: Replicate consistency assessment")

rep_dist_df <- pca_df %>%
  select(sample_name, condition, time_point, time_label,
         replicate, PC1, PC2, PC3) %>%
  group_by(condition, time_point) %>%
  summarise(
    n_rep          = n(),
    sd_PC1         = sd(PC1),
    sd_PC2         = sd(PC2),
    mean_euclidean = {
      coords <- cbind(PC1, PC2, PC3)
      dists  <- as.vector(dist(coords, method = "euclidean"))
      mean(dists)
    },
    .groups = "drop"
  ) %>%
  mutate(time_label = factor(unname(cfg$time_labels[as.character(time_point)]),
                             levels = unname(cfg$time_labels)))

write_csv_safe(rep_dist_df,
               file.path(cfg$results_dir, "tables", "qc_replicate_pc_distances.csv"))

p_rep_dist <- ggplot(rep_dist_df,
                     aes(x = time_label, y = mean_euclidean,
                         fill = condition, group = condition)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6, colour = "grey25", linewidth = 0.15) +
  geom_point(position = position_dodge(width = 0.7), size = 2) +
  scale_fill_manual(values = condition_cols) +
  labs(
    title    = "Replicate consistency — mean pairwise PC distance",
    subtitle = "Euclidean distance in PC1–3 space; lower = more consistent replicates",
    x        = "Time point",
    y        = "Mean pairwise Euclidean distance (PC1–3)",
    fill     = "Condition"
  ) +
  pub_theme() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_plot(p_rep_dist,
          file.path(cfg$results_dir, "plots", "pca", "pca3_replicate_consistency"),
          width = 7.5, height = 5)

log_msg("  ✔ Replicate consistency plot saved")


# ==============================================================================
#  PCA-4 : Temporal Trajectory — mean PC centroid per phase
# ==============================================================================

log_msg("PCA-4: Temporal trajectory (centroid path)")

centroid_df <- pca_df %>%
  group_by(condition, time_point) %>%
  summarise(
    PC1 = mean(PC1), PC2 = mean(PC2), PC3 = mean(PC3),
    .groups = "drop"
  ) %>%
  mutate(
    time_order = match(as.character(time_point), cfg$time_levels),
    time_label = factor(unname(cfg$time_labels[as.character(time_point)]),
                        levels = unname(cfg$time_labels))
  ) %>%
  arrange(condition, time_order)

write_csv_safe(centroid_df,
               file.path(cfg$results_dir, "tables", "PCA_temporal_centroids.csv"))

p_traj <- ggplot(centroid_df,
                 aes(x = PC1, y = PC2, colour = condition)) +
  geom_path(aes(group = condition),
            arrow = arrow(length = unit(0.25, "cm"),
                          type = "closed", ends = "last"),
            linewidth = 0.9, alpha = 0.7) +
  geom_point(aes(shape = time_point), size = 5, alpha = 0.95) +
  geom_label_repel(
    aes(label = time_label),
    size = 3, box.padding = 0.45, max.overlaps = 20,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = condition_cols, name = "Condition") +
  scale_shape_manual(values = c(T0 = 15, T60 = 16, T180 = 17,
                                T240 = 18, T5d = 8),
                     labels = cfg$time_labels,
                     name   = "Time point") +
  labs(
    title    = "Temporal transcriptomic trajectory",
    subtitle = "Mean PC1/PC2 centroid per condition × time point",
    x        = sprintf("PC1 (%.1f%%)", pca_percent[1]),
    y        = sprintf("PC2 (%.1f%%)", pca_percent[2])
  ) +
  pub_theme(base_size = 12)

save_plot(p_traj,
          file.path(cfg$results_dir, "plots", "pca", "pca4_temporal_trajectory"),
          width = 8, height = 6.5)

log_msg("  ✔ Temporal trajectory plot saved")


# ==============================================================================
#  PCA-5 : PC Loading contributions (top drivers of PC1 & PC2)
# ==============================================================================

log_msg("PCA-5: PC loading contributions")

loadings_df <- as.data.frame(pca_obj$rotation[, 1:4]) %>%
  rownames_to_column("gene_id")

write_csv_safe(loadings_df,
               file.path(cfg$results_dir, "tables", "PCA_gene_loadings_PC1_to_PC4.csv"))

# Top 30 contributors for PC1
top_pc1 <- loadings_df %>%
  slice_max(abs(PC1), n = 30) %>%
  arrange(PC1)

p_load1 <- ggplot(top_pc1,
                  aes(x = reorder(gene_id, PC1), y = PC1,
                      fill = PC1 > 0)) +
  geom_col(width = 0.75, colour = NA) +
  scale_fill_manual(values = c("TRUE" = "#DC2626", "FALSE" = "#2563EB"),
                    labels = c("TRUE" = "Positive (salt-induced)",
                               "FALSE" = "Negative (salt-repressed)"),
                    name = "Direction") +
  coord_flip() +
  labs(
    title    = sprintf("Top 30 gene loadings on PC1 (%.1f%% variance)",
                       pca_percent[1]),
    subtitle = "Largest contributors to sample separation along PC1",
    x        = NULL,
    y        = "PC1 loading"
  ) +
  pub_theme(base_size = 11)

save_plot(p_load1,
          file.path(cfg$results_dir, "plots", "pca", "pca5_PC1_gene_loadings"),
          width = 7.5, height = 7)

# Top 30 contributors for PC2
top_pc2 <- loadings_df %>%
  slice_max(abs(PC2), n = 30) %>%
  arrange(PC2)

p_load2 <- ggplot(top_pc2,
                  aes(x = reorder(gene_id, PC2), y = PC2,
                      fill = PC2 > 0)) +
  geom_col(width = 0.75, colour = NA) +
  scale_fill_manual(values = c("TRUE" = "#009E73", "FALSE" = "#E69F00"),
                    name = "Direction") +
  coord_flip() +
  labs(
    title    = sprintf("Top 30 gene loadings on PC2 (%.1f%% variance)",
                       pca_percent[2]),
    x        = NULL,
    y        = "PC2 loading"
  ) +
  pub_theme(base_size = 11)

save_plot(p_load2,
          file.path(cfg$results_dir, "plots", "pca", "pca5_PC2_gene_loadings"),
          width = 7.5, height = 7)

log_msg("  ✔ PC loading plots saved")


# ==============================================================================
#  PCA-6 : Condition separation per time point (PC1 strip plot)
# ==============================================================================

log_msg("PCA-6: Per-time-point PC1 separation between conditions")

p_pc1_strip <- ggplot(pca_df,
                      aes(x = time_label, y = PC1,
                          colour = condition, shape = condition)) +
  geom_jitter(width = 0.07, height = 0, size = 3.5, alpha = 0.9) +
  stat_summary(aes(group = condition),
               fun = mean, geom = "crossbar",
               fatten = 2, width = 0.35,
               position = position_dodge(width = 0)) +
  scale_colour_manual(values = condition_cols) +
  scale_shape_manual(values = c(Control = 16, Salt = 17)) +
  labs(
    title    = "PC1 score per sample — condition separation over time",
    subtitle = "Crossbar marks group mean; circles = Control, triangles = Salt",
    x        = "Time point",
    y        = sprintf("PC1 (%.1f%%)", pca_percent[1]),
    colour   = "Condition", shape = "Condition"
  ) +
  pub_theme()

save_plot(p_pc1_strip,
          file.path(cfg$results_dir, "plots", "pca", "pca6_PC1_by_timepoint"),
          width = 7.5, height = 5)

log_msg("  ✔ PC1 strip plot saved")


# ==============================================================================
#  SUMMARY REPORT — print to console
# ==============================================================================

log_msg("Generating summary report to console")

cat("\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat("  QC & PCA SUMMARY — Oryza sativa salt stress (GSE284308)\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

cat("  DATASET\n")
cat(sprintf("  ├─ Genes in VST matrix       : %d\n", nrow(vst_matrix)))
cat(sprintf("  ├─ Samples                   : %d\n", ncol(vst_matrix)))
cat(sprintf("  └─ Time points               : %s\n\n",
            paste(cfg$time_levels, collapse = " → ")))

cat("  LIBRARY SIZES (raw estimated counts)\n")
lib_tbl <- lib_size_df %>%
  select(sample_name, condition, time_label, raw_counts) %>%
  mutate(raw_millions = round(raw_counts / 1e6, 2))
print(as.data.frame(lib_tbl[, c("sample_name","condition","time_label","raw_millions")]),
      row.names = FALSE)

cat("\n  SIZE FACTORS\n")
print(as.data.frame(sf_df[, c("sample_name","condition","time_label","size_factor")]),
      row.names = FALSE)

cat("\n  SAMPLE CORRELATIONS (Pearson, VST) — min / mean / max\n")
cat(sprintf("  ├─ Min  : %.4f\n", min(sample_cor[lower.tri(sample_cor)])))
cat(sprintf("  ├─ Mean : %.4f\n", mean(sample_cor[lower.tri(sample_cor)])))
cat(sprintf("  └─ Max  : %.4f\n\n", max(sample_cor[lower.tri(sample_cor)])))

cat("  PCA — VARIANCE EXPLAINED\n")
print(scree_df, n = n_pcs)

cat("\n  REPLICATE CONSISTENCY (mean Euclidean distance in PC1–3)\n")
print(as.data.frame(
  rep_dist_df[, c("condition","time_label","mean_euclidean")]
), row.names = FALSE)

cat("\n  OUTPUT FILES\n")
all_outputs <- list.files(
  path       = file.path(cfg$results_dir, c("plots/qc","plots/pca","tables")),
  pattern    = "qc|PCA|pca",
  full.names = TRUE,
  recursive  = FALSE
)
invisible(lapply(sort(all_outputs), function(f) cat("  ", f, "\n")))

cat("\n══════════════════════════════════════════════════════════════════════\n")
log_msg("Section 1 (QC & PCA) complete.")


#!/usr/bin/env Rscript

# ==============================================================================
#  SECTION 2 : Differential Expression Analysis
#  Project    : rice_salt_timeseries_10mM_NaCl  (GSE284308)
#  Organism   : Oryza sativa
#  Design     : Salt time course — T60 / T180 / T240 / T5d vs T0
#               + sequential contrasts  T180 vs T60  |  T5d vs T240
#
#  Requires   : my_new_analysis.R outputs
#                 results/tables/dds_salt_timecourse.rds
#                 results/tables/vst_all_samples.rds
#                 results/tables/metadata_clean.csv
#                 results/tables/DESeq2_all_timepoint_contrasts.csv
#
#  Outputs
#    results/tables/
#      DEG_summary_by_timepoint.csv
#      DEG_updown_proportions.csv
#      DEG_overlap_matrix.csv
#      DEG_sequential_*.csv
#      DEG_master_annotated.csv
#    results/plots/volcano/
#      EnhancedVolcano_*  (PDF + PNG)
#      volcano_panel_all_timepoints.*
#    results/plots/heatmaps/
#      heatmap_top100_DEGs.*
#      heatmap_zscored_clusters.*
#    results/plots/overlap/
#      venn_up_DEGs.*  |  venn_down_DEGs.*
#      upset_DEG_overlap.*
#      barplot_updown_proportions.*
#      barplot_deg_counts_timeline.*
# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(20260505)

# ── cfg (identical to parent) ─────────────────────────────────────────────────
cfg <- list(
  project_name  = "rice_salt_timeseries_10mM_NaCl",
  results_dir   = "results",
  alpha         = 0.05,
  lfc_cutoff    = 1,
  min_count     = 10,
  min_samples   = 3,
  heatmap_top_n = 100,
  top_labels    = 15,
  baseline_time = "T0",
  time_levels   = c("T0", "T60", "T180", "T240", "T5d"),
  time_labels   = c(T0 = "0 min", T60 = "60 min", T180 = "180 min",
                    T240 = "240 min", T5d = "5 days"),
  auto_install  = TRUE
)

dirs <- file.path(cfg$results_dir,
                  c("tables", "plots/volcano", "plots/heatmaps",
                    "plots/overlap", "plots/qc"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              paste0(..., collapse = "")))
}

# ── packages ──────────────────────────────────────────────────────────────────
cran_pkgs <- c("tidyverse", "ggplot2", "ggrepel", "patchwork", "pheatmap",
               "RColorBrewer", "ggpubr", "data.table", "scales",
               "ggvenn", "UpSetR", "BiocManager")
bioc_pkgs <- c("DESeq2", "EnhancedVolcano", "ComplexHeatmap", "circlize")

if (isTRUE(cfg$auto_install)) {
  miss_c <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss_c)) install.packages(miss_c, repos = "https://cloud.r-project.org")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  miss_b <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss_b)) BiocManager::install(miss_b, ask = FALSE, update = FALSE)
  # ggvenn from GitHub if needed
  if (!requireNamespace("ggvenn", quietly = TRUE))
    tryCatch(remotes::install_github("yanlinlin82/ggvenn"), error = function(e) NULL)
}

suppressPackageStartupMessages({
  library(DESeq2); library(tidyverse); library(ggplot2); library(ggrepel)
  library(patchwork); library(pheatmap); library(RColorBrewer)
  library(ggpubr); library(data.table); library(scales); library(EnhancedVolcano)
})
has_ch  <- requireNamespace("ComplexHeatmap", quietly = TRUE) &&
  requireNamespace("circlize",       quietly = TRUE)
has_ggvenn <- requireNamespace("ggvenn",  quietly = TRUE)
has_upset  <- requireNamespace("UpSetR",  quietly = TRUE)
if (has_ch) suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize) })
if (has_ggvenn) library(ggvenn)
if (has_upset)  library(UpSetR)

# ── shared aesthetics (matches parent script) ─────────────────────────────────
time_cols <- c(T0="#4B5563", T60="#0072B2", T180="#009E73",
               T240="#E69F00", T5d="#D55E00")
condition_cols <- c(Control="#6B7280", Salt="#DC2626")

pub_theme <- function(base_size = 12) {
  theme_pubr(base_size = base_size, border = FALSE) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0, size = base_size + 2),
      plot.subtitle = element_text(color = "grey25"),
      axis.title    = element_text(face = "bold"),
      legend.title  = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text    = element_text(face = "bold")
    )
}

save_plot <- function(plot, basename, width = 8, height = 6, dpi = 1200) {
  pdf_dev <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
  ggsave(paste0(basename, ".pdf"), plot, width = width, height = height, device = pdf_dev)
  ggsave(paste0(basename, ".tiff"), plot, width = width, height = height, dpi = dpi,
         compression = "lzw", device = "tiff")
}
write_csv_safe <- function(x, path) data.table::fwrite(as.data.frame(x), path)


# ==============================================================================
#  LOAD OBJECTS
# ==============================================================================
log_msg("Loading objects from my_new_analysis.R")

dds_salt   <- readRDS(file.path(cfg$results_dir, "tables", "dds_salt_timecourse.rds"))
vsd_all    <- readRDS(file.path(cfg$results_dir, "tables", "vst_all_samples.rds"))
vst_matrix <- assay(vsd_all)

metadata_clean <- data.table::fread(
  file.path(cfg$results_dir, "tables", "metadata_clean.csv")) %>%
  as_tibble() %>%
  mutate(condition  = factor(condition,  levels = c("Control","Salt")),
          time_point = factor(time_point, levels = cfg$time_levels))

salt_samples <- metadata_clean %>%
  filter(condition == "Salt") %>%
  pull(sample_name)
salt_vst <- vst_matrix[, salt_samples]

deg_all <- data.table::fread(
  file.path(cfg$results_dir, "tables", "DESeq2_all_timepoint_contrasts.csv")) %>%
  as_tibble() %>%
  mutate(time_point = factor(time_point, levels = cfg$time_levels[-1]))  # T60–T5d

log_msg(sprintf("  Loaded: %d gene×contrast rows | %d samples",
                nrow(deg_all), ncol(vst_matrix)))


# ==============================================================================
#  STEP 2-A : Sequential contrasts  (Early vs Mid, Mid vs Late proxy)
#             T180_vs_T60  |  T5d_vs_T240
# ==============================================================================
log_msg("Step 2-A: Sequential contrasts")

has_apeglm <- requireNamespace("apeglm", quietly = TRUE)
shrink_type <- if (has_apeglm) "apeglm" else "normal"

# Re-level to T60 → extract T180 vs T60
dds_seq60 <- dds_salt
dds_seq60$time_point <- relevel(dds_seq60$time_point, ref = "T60")
dds_seq60 <- suppressMessages(DESeq(dds_seq60))

res_180v60 <- lfcShrink(dds_seq60,
                        coef = "time_point_T180_vs_T60",
                        type = shrink_type) %>%
  as.data.frame() %>% rownames_to_column("gene_id") %>% as_tibble() %>%
  mutate(contrast_id = "T180_vs_T60",
         comparison  = "180 min vs 60 min",
         time_point  = "T180",
         shrinkage   = shrink_type,
         padj        = replace_na(padj, 1),
         log2FoldChange = replace_na(log2FoldChange, 0),
         regulation  = case_when(
           padj < cfg$alpha & log2FoldChange >=  cfg$lfc_cutoff ~ "Up",
           padj < cfg$alpha & log2FoldChange <= -cfg$lfc_cutoff ~ "Down",
           TRUE ~ "Not significant"))

# Re-level to T240 → extract T5d vs T240
dds_seq240 <- dds_salt
dds_seq240$time_point <- relevel(dds_seq240$time_point, ref = "T240")
dds_seq240 <- suppressMessages(DESeq(dds_seq240))

res_5dv240 <- lfcShrink(dds_seq240,
                        coef = "time_point_T5d_vs_T240",
                        type = shrink_type) %>%
  as.data.frame() %>% rownames_to_column("gene_id") %>% as_tibble() %>%
  mutate(contrast_id = "T5d_vs_T240",
         comparison  = "5 days vs 240 min",
         time_point  = "T5d",
         shrinkage   = shrink_type,
         padj        = replace_na(padj, 1),
         log2FoldChange = replace_na(log2FoldChange, 0),
         regulation  = case_when(
           padj < cfg$alpha & log2FoldChange >=  cfg$lfc_cutoff ~ "Up",
           padj < cfg$alpha & log2FoldChange <= -cfg$lfc_cutoff ~ "Down",
           TRUE ~ "Not significant"))

write_csv_safe(res_180v60,
               file.path(cfg$results_dir, "tables", "DEG_sequential_T180_vs_T60.csv"))
write_csv_safe(res_5dv240,
               file.path(cfg$results_dir, "tables", "DEG_sequential_T5d_vs_T240.csv"))

seq_summary <- bind_rows(res_180v60, res_5dv240) %>%
  count(comparison, regulation) %>%
  pivot_wider(names_from = regulation, values_from = n, values_fill = 0)
write_csv_safe(seq_summary,
               file.path(cfg$results_dir, "tables", "DEG_sequential_summary.csv"))

log_msg("  ✔ Sequential contrasts done")
cat("\n  Sequential contrast summary:\n"); print(seq_summary)


# ==============================================================================
#  STEP 2-B : DEG count summary & timeline bar chart
# ==============================================================================
log_msg("Step 2-B: DEG count summary")

deg_summary <- deg_all %>%
  count(comparison, time_point, regulation) %>%
  complete(nesting(comparison, time_point),
           regulation = c("Up","Down","Not significant"),
           fill = list(n = 0)) %>%
  mutate(time_point = factor(time_point, levels = cfg$time_levels))

write_csv_safe(deg_summary,
               file.path(cfg$results_dir, "tables", "DEG_summary_by_timepoint.csv"))

# ── Bidirectional bar chart ───────────────────────────────────────────────────
deg_bar_df <- deg_summary %>%
  filter(regulation != "Not significant") %>%
  mutate(
    n_signed   = if_else(regulation == "Down", -n, n),
    time_label = factor(unname(cfg$time_labels[as.character(time_point)]),
                        levels = unname(cfg$time_labels))
  )

p_deg_bar <- ggplot(deg_bar_df,
                    aes(x = time_label, y = n_signed, fill = regulation)) +
  geom_col(width = 0.65, colour = "grey20", linewidth = 0.15) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text(aes(label = abs(n_signed),
                vjust  = if_else(n_signed > 0, -0.3, 1.1)),
            size = 3.2, colour = "grey20") +
  scale_fill_manual(values = c(Up = "#DC2626", Down = "#2563EB")) +
  scale_y_continuous(labels = abs,
                     expand  = expansion(mult = c(0.12, 0.12))) +
  labs(
    title    = "DEG counts across the salt-stress time course",
    subtitle = sprintf("DESeq2 | padj < %.2f | |log₂FC| > %d | %s shrinkage",
                       cfg$alpha, cfg$lfc_cutoff, shrink_type),
    x        = "Time point (salt vs T0)",
    y        = "Number of significant DEGs",
    fill     = "Regulation"
  ) +
  pub_theme()

save_plot(p_deg_bar,
          file.path(cfg$results_dir, "plots", "overlap", "barplot_deg_counts_timeline"),
          width = 7.5, height = 5)

# ── Up/Down proportion chart ─────────────────────────────────────────────────
deg_prop <- deg_summary %>%
  filter(regulation != "Not significant") %>%
  group_by(time_point) %>%
  mutate(total = sum(n),
         prop  = n / total,
         time_label = factor(unname(cfg$time_labels[as.character(time_point)]),
                             levels = unname(cfg$time_labels)))

write_csv_safe(deg_prop,
               file.path(cfg$results_dir, "tables", "DEG_updown_proportions.csv"))

p_prop <- ggplot(deg_prop,
                 aes(x = time_label, y = prop, fill = regulation)) +
  geom_col(width = 0.65, colour = "grey20", linewidth = 0.15) +
  geom_text(aes(label = paste0(round(prop * 100, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 3.2, colour = "white", fontface = "bold") +
  scale_fill_manual(values = c(Up = "#DC2626", Down = "#2563EB")) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title    = "Proportion of up- vs down-regulated DEGs per time point",
    x        = NULL, y = "Proportion", fill = "Regulation"
  ) +
  pub_theme()

save_plot(p_prop,
          file.path(cfg$results_dir, "plots", "overlap", "barplot_updown_proportions"),
          width = 7.5, height = 5)

log_msg("  ✔ DEG count and proportion charts saved")


# ==============================================================================
#  STEP 2-C : Volcano plots  — per time point + combined panel
# ==============================================================================
log_msg("Step 2-C: Volcano plots")

contrast_ids <- unique(deg_all$contrast_id)

make_volcano <- function(res_df, cid, label) {
  top_genes <- res_df %>%
    filter(padj < cfg$alpha, abs(log2FoldChange) > cfg$lfc_cutoff) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = cfg$top_labels) %>%
    pull(gene_id)
  
  p <- EnhancedVolcano(
    as.data.frame(res_df),
    lab        = res_df$gene_id,
    x          = "log2FoldChange",
    y          = "padj",
    selectLab  = top_genes,
    xlab       = expression(Log[2]~"fold change"),
    ylab       = expression("-"~Log[10]~"adjusted P"),
    title      = paste0(label, " — Oryza sativa salt stress"),
    subtitle   = sprintf("DESeq2 | %s shrinkage | padj < %.2f | |log₂FC| > %d",
                         unique(res_df$shrinkage), cfg$alpha, cfg$lfc_cutoff),
    pCutoff    = cfg$alpha,
    FCcutoff   = cfg$lfc_cutoff,
    pointSize  = 1.7,
    labSize    = 3.2,
    boxedLabels   = FALSE,
    drawConnectors = TRUE,
    widthConnectors = 0.45,
    colConnectors   = "grey35",
    max.overlaps    = Inf,
    col = c("grey72","#6B7280","#2563EB","#DC2626"),
    legendLabels    = c("NS","Log2FC","padj","padj & Log2FC"),
    legendPosition  = "right",
    legendLabSize   = 10,
    legendIconSize  = 3,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    border          = "partial",
    borderWidth     = 0.4,
    borderColour    = "grey35"
  ) + theme(plot.title    = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size  = 9,     colour = "grey25"))
  
  save_plot(p,
            file.path(cfg$results_dir, "plots", "volcano",
                      paste0("EnhancedVolcano_", cid)),
            width = 8.5, height = 7)
  invisible(p)
}

volcano_list <- lapply(contrast_ids, function(cid) {
  sub <- filter(deg_all, contrast_id == cid)
  lbl <- unique(sub$comparison)
  make_volcano(sub, cid, lbl)
})
names(volcano_list) <- contrast_ids

# ── Combined 2×2 panel volcano ────────────────────────────────────────────────
panel_list <- lapply(contrast_ids, function(cid) {
  df  <- filter(deg_all, contrast_id == cid)
  lbl <- unique(df$comparison)
  up  <- sum(df$regulation == "Up")
  dn  <- sum(df$regulation == "Down")
  
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj),
                 colour = regulation)) +
    geom_point(alpha = 0.45, size = 0.9) +
    geom_point(data = filter(df, regulation != "Not significant"),
               alpha = 0.85, size = 1.4) +
    geom_vline(xintercept = c(-cfg$lfc_cutoff, cfg$lfc_cutoff),
               linetype = "dashed", colour = "grey50", linewidth = 0.45) +
    geom_hline(yintercept = -log10(cfg$alpha),
               linetype = "dashed", colour = "grey50", linewidth = 0.45) +
    geom_text_repel(
      data = filter(df, regulation != "Not significant") %>%
        arrange(padj) %>% slice_head(n = 10),
      aes(label = gene_id), size = 2.2, max.overlaps = 15,
      segment.alpha = 0.5, show.legend = FALSE
    ) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("↑%d  ↓%d", up, dn),
             hjust = 1.1, vjust = 1.4, size = 3.5) +
    scale_colour_manual(
      values = c(Up = "#DC2626", Down = "#2563EB", `Not significant` = "grey72"),
      guide  = "none") +
    labs(title = lbl,
         x     = expression(log[2]~FC),
         y     = expression(-log[10]~padj)) +
    pub_theme(base_size = 10)
})

p_panel <- wrap_plots(panel_list, ncol = 2) +
  plot_annotation(
    title    = "Differential expression — all time points vs T0",
    subtitle = "Oryza sativa 10 mM NaCl salt stress (GSE284308)",
    theme    = theme(plot.title    = element_text(face = "bold", size = 14),
                     plot.subtitle = element_text(colour = "grey25"))
  )

save_plot(p_panel,
          file.path(cfg$results_dir, "plots", "volcano",
                    "volcano_panel_all_timepoints"),
          width = 13, height = 11)

log_msg("  ✔ Volcano plots saved")


# ==============================================================================
#  STEP 2-D : DEG overlap — Venn & UpSet
# ==============================================================================
log_msg("Step 2-D: DEG overlap (Venn + UpSet)")

# Gene sets per time point
up_sets <- lapply(setNames(contrast_ids, contrast_ids), function(cid)
  filter(deg_all, contrast_id == cid, regulation == "Up")$gene_id)
dn_sets <- lapply(setNames(contrast_ids, contrast_ids), function(cid)
  filter(deg_all, contrast_id == cid, regulation == "Down")$gene_id)

names(up_sets) <- unname(cfg$time_labels[
  sub("_vs_T0","", contrast_ids)])
names(dn_sets) <- names(up_sets)

# ── Overlap matrix ────────────────────────────────────────────────────────────
make_overlap_matrix <- function(sets) {
  nms <- names(sets)
  mat <- outer(nms, nms, Vectorize(function(a, b)
    length(intersect(sets[[a]], sets[[b]]))))
  rownames(mat) <- colnames(mat) <- nms
  mat
}
ov_up <- make_overlap_matrix(up_sets)
ov_dn <- make_overlap_matrix(dn_sets)
write_csv_safe(as.data.frame(ov_up) %>% rownames_to_column("time"),
               file.path(cfg$results_dir, "tables", "DEG_overlap_matrix_up.csv"))
write_csv_safe(as.data.frame(ov_dn) %>% rownames_to_column("time"),
               file.path(cfg$results_dir, "tables", "DEG_overlap_matrix_down.csv"))

# ── ggvenn (up-regulated) ─────────────────────────────────────────────────────
if (has_ggvenn && length(up_sets) <= 4) {
  p_venn_up <- ggvenn(up_sets,
                      fill_color = unname(time_cols[-1]),
                      stroke_size = 0.5, set_name_size = 3.5,
                      text_size   = 3) +
    labs(title = "Up-regulated DEG overlap across time points") +
    pub_theme()
  save_plot(p_venn_up,
            file.path(cfg$results_dir, "plots", "overlap", "venn_up_DEGs"),
            width = 8, height = 7)
  
  p_venn_dn <- ggvenn(dn_sets,
                      fill_color = c("#93C5FD","#60A5FA","#3B82F6","#1D4ED8"),
                      stroke_size = 0.5, set_name_size = 3.5,
                      text_size   = 3) +
    labs(title = "Down-regulated DEG overlap across time points") +
    pub_theme()
  save_plot(p_venn_dn,
            file.path(cfg$results_dir, "plots", "overlap", "venn_down_DEGs"),
            width = 8, height = 7)
  log_msg("  ✔ Venn diagrams saved")
}

# ── UpSet plot (all time points) ──────────────────────────────────────────────
if (has_upset) {
  make_upset_df <- function(sets) {
    all_genes <- unique(unlist(sets))
    as.data.frame(setNames(
      lapply(sets, function(s) as.integer(all_genes %in% s)),
      names(sets)
    ), row.names = all_genes)
  }
  
  upset_up <- make_upset_df(up_sets)
  upset_dn <- make_upset_df(dn_sets)
  
  pdf(file.path(cfg$results_dir, "plots", "overlap", "upset_up_DEGs.pdf"),
      width = 10, height = 6)
  print(upset(upset_up, order.by = "freq", nsets = length(up_sets),
              main.bar.color = "#DC2626", sets.bar.color = time_cols[-1],
              text.scale = 1.3,
              mainbar.y.label = "Intersection size",
              sets.x.label   = "Up-DEGs per time point",
              mb.ratio = c(0.65, 0.35)))
  dev.off()
  tiff(file.path(cfg$results_dir, "plots", "overlap", "upset_up_DEGs.tiff"),
      width = 10, height = 6, units = "in", res = 1200, compression = "lzw")
  print(upset(upset_up, order.by = "freq", nsets = length(up_sets),
              main.bar.color = "#DC2626", sets.bar.color = time_cols[-1],
              text.scale = 1.3,
              mainbar.y.label = "Intersection size",
              sets.x.label   = "Up-DEGs per time point",
              mb.ratio = c(0.65, 0.35)))
  dev.off()
  
  pdf(file.path(cfg$results_dir, "plots", "overlap", "upset_down_DEGs.pdf"),
      width = 10, height = 6)
  print(upset(upset_dn, order.by = "freq", nsets = length(dn_sets),
              main.bar.color = "#2563EB", sets.bar.color = time_cols[-1],
              text.scale = 1.3,
              mainbar.y.label = "Intersection size",
              sets.x.label   = "Down-DEGs per time point"))
  dev.off()
  tiff(file.path(cfg$results_dir, "plots", "overlap", "upset_down_DEGs.tiff"),
      width = 10, height = 6, units = "in", res = 1200, compression = "lzw")
  print(upset(upset_dn, order.by = "freq", nsets = length(dn_sets),
              main.bar.color = "#2563EB", sets.bar.color = time_cols[-1],
              text.scale = 1.3,
              mainbar.y.label = "Intersection size",
              sets.x.label   = "Down-DEGs per time point"))
  dev.off()
  log_msg("  ✔ UpSet plots saved")
}


# ==============================================================================
#  STEP 2-E : DEG heatmaps
# ==============================================================================
log_msg("Step 2-E: DEG expression heatmaps")

# Top N DEGs by minimum adjusted p-value across all contrasts
top_degs <- deg_all %>%
  filter(regulation != "Not significant") %>%
  group_by(gene_id) %>%
  summarise(min_padj = min(padj), max_lfc = max(abs(log2FoldChange)),
            .groups = "drop") %>%
  arrange(min_padj, desc(max_lfc)) %>%
  slice_head(n = cfg$heatmap_top_n) %>%
  pull(gene_id)

# Subset and z-score VST matrix
salt_samples <- metadata_clean %>%
  filter(condition == "Salt") %>%
  pull(sample_name)
vst_sub  <- vst_matrix[top_degs, salt_samples, drop = FALSE]
vst_z    <- t(scale(t(vst_sub)))
vst_z[is.nan(vst_z)] <- 0

# Annotation bars
ann_col <- metadata_clean %>%
  filter(condition == "Salt") %>%
  arrange(time_point) %>%
  select(sample_name, time_point) %>%
  as.data.frame() %>%
  column_to_rownames("sample_name")
ann_colors <- list(time_point = time_cols[cfg$time_levels[-1]])

col_order <- metadata_clean %>%
  filter(condition == "Salt") %>%
  arrange(time_point) %>%
  pull(sample_name)
vst_z <- vst_z[, col_order, drop = FALSE]

# ── pheatmap: top 100 DEGs ────────────────────────────────────────────────────
ann_colors <- list(time_point = time_cols[cfg$time_levels])  # Include T0

# Or if time_cols is not defined, define it:
time_cols <- c(T0 = "#4B5563", T60 = "#0072B2", T180 = "#009E73", 
               T240 = "#E69F00", T5d = "#D55E00")

ann_colors <- list(time_point = time_cols)

pheatmap(
  vst_z,
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  clustering_method = "ward.D2",
  cluster_cols      = FALSE,
  show_rownames     = nrow(vst_z) <= 60,
  border_color      = NA,
  fontsize          = 8,
  fontsize_row      = 6,
  main              = sprintf("Top %d DEGs — z-scored VST expression (salt samples)",
                              cfg$heatmap_top_n),
  filename = file.path(cfg$results_dir, "plots", "heatmaps",
                       "heatmap_top100_DEGs.pdf"),
  width = 9, height = 11
)

tiff(file.path(cfg$results_dir, "plots", "heatmaps", "heatmap_top100_DEGs.tiff"),
    width = 9, height = 11, units = "in", res = 1200, compression = "lzw")
pheatmap(vst_z,
         annotation_col = ann_col, annotation_colors = ann_colors,
         color          = colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
         clustering_method = "ward.D2", cluster_cols = FALSE,
         show_rownames  = nrow(vst_z) <= 60,
         border_color   = NA, fontsize = 8, fontsize_row = 6,
         main = sprintf("Top %d DEGs — z-scored VST expression", cfg$heatmap_top_n))
dev.off()

# ── k-means clustered heatmap (6 clusters) ───────────────────────────────────
k_clust <- 6
set.seed(20260505)
km     <- kmeans(vst_z, centers = k_clust, nstart = 25, iter.max = 100)
k_labs <- paste0("Cluster ", km$cluster)

if (has_ch) {
  col_fun <- colorRamp2(c(-2, 0, 2),
                        c("#2563EB", "white", "#DC2626"))
  ht <- Heatmap(
    vst_z,
    name               = "z-score",
    col                = col_fun,
    cluster_rows       = TRUE,
    cluster_columns    = FALSE,
    split              = km$cluster,
    column_title       = "Salt-stressed time course →",
    row_title_rot      = 0,
    show_row_names     = FALSE,
    show_column_names  = TRUE,
    column_names_gp    = gpar(fontsize = 7),
    top_annotation     = HeatmapAnnotation(
      Time = ann_col$time_point,
      col  = ann_colors,
      annotation_name_side = "left"
    ),
    heatmap_legend_param = list(title = "z-score", at = c(-2,0,2))
  )
  pdf(file.path(cfg$results_dir, "plots", "heatmaps",
                "heatmap_zscored_clusters_ComplexHeatmap.pdf"),
      width = 9, height = 12)
  draw(ht, heatmap_legend_side = "right",
       column_title = "Top DEGs — temporal expression clusters")
  dev.off()
  tiff(file.path(cfg$results_dir, "plots", "heatmaps",
                "heatmap_zscored_clusters_ComplexHeatmap.tiff"),
      width = 9, height = 12, units = "in", res = 1200, compression = "lzw")
  draw(ht, heatmap_legend_side = "right",
       column_title = "Top DEGs — temporal expression clusters")
  dev.off()
}

# Save cluster assignments
clust_df <- tibble(gene_id = names(km$cluster), cluster = km$cluster) %>%
  arrange(cluster)
write_csv_safe(clust_df,
               file.path(cfg$results_dir, "tables", "DEG_kmeans_cluster_assignments.csv"))

# ── Cluster mean trajectory plot ──────────────────────────────────────────────
traj_df <- as.data.frame(vst_z) %>%
  rownames_to_column("gene_id") %>%
  left_join(clust_df, by = "gene_id") %>%
  pivot_longer(-c(gene_id, cluster), names_to = "sample_name", values_to = "zscore") %>%
  left_join(metadata_clean %>% select(sample_name, time_point), by = "sample_name") %>%
  group_by(cluster, time_point) %>%
  summarise(mean_z = mean(zscore), se_z = sd(zscore)/sqrt(n()), .groups = "drop") %>%
  mutate(time_label = factor(unname(cfg$time_labels[as.character(time_point)]),
                             levels = unname(cfg$time_labels[-1])))

p_traj <- ggplot(traj_df,
                 aes(x = time_label, y = mean_z,
                     colour = factor(cluster), group = cluster)) +
  geom_ribbon(aes(ymin = mean_z - se_z, ymax = mean_z + se_z,
                  fill = factor(cluster)), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  facet_wrap(~paste0("Cluster ", cluster), ncol = 3) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  scale_colour_brewer(palette = "Set1", guide = "none") +
  scale_fill_brewer(palette   = "Set1", guide = "none") +
  labs(
    title    = "DEG cluster mean expression trajectories (salt samples)",
    subtitle = sprintf("k-means k = %d; ribbon = ±1 SE", k_clust),
    x        = "Time point",
    y        = "Mean z-score (VST)"
  ) +
  pub_theme(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_plot(p_traj,
          file.path(cfg$results_dir, "plots", "heatmaps",
                    "heatmap_cluster_trajectories"),
          width = 11, height = 8)

log_msg("  ✔ DEG heatmaps and cluster trajectories saved")


# ==============================================================================
#  STEP 2-F : Master annotated DEG table
# ==============================================================================
log_msg("Step 2-F: Master annotated DEG table")

# Pivot wide: one row per gene, columns for each contrast
deg_wide <- deg_all %>%
  select(gene_id, contrast_id, log2FoldChange, padj, regulation) %>%
  pivot_wider(names_from  = contrast_id,
              values_from = c(log2FoldChange, padj, regulation),
              names_glue  = "{contrast_id}_{.value}")

master_df <- deg_wide %>%
  mutate(
    n_contrasts_sig = rowSums(across(ends_with("_regulation"),
                                     ~ . != "Not significant"), na.rm = TRUE),
    ever_up    = if_any(ends_with("_regulation"), ~ . == "Up"),
    ever_down  = if_any(ends_with("_regulation"), ~ . == "Down"),
    phase      = case_when(
      !ever_up & !ever_down ~ "Not significant",
      ever_up  & !ever_down ~ "Up only",
      !ever_up &  ever_down ~ "Down only",
      TRUE                  ~ "Mixed"
    )
  ) %>%
  arrange(desc(n_contrasts_sig))

write_csv_safe(master_df,
               file.path(cfg$results_dir, "tables", "DEG_master_annotated.csv"))

log_msg("  ✔ Master DEG table saved")


# ==============================================================================
#  CONSOLE SUMMARY
# ==============================================================================
cat("\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat("  SECTION 2 SUMMARY — DEG Analysis\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

cat("  DEG COUNTS (salt vs T0)\n")
print(as.data.frame(
  deg_summary %>%
    filter(regulation != "Not significant") %>%
    select(time_point, regulation, n) %>%
    pivot_wider(names_from = regulation, values_from = n, values_fill = 0)
), row.names = FALSE)

cat("\n  SEQUENTIAL CONTRAST SUMMARY\n")
print(seq_summary, row.names = FALSE)

cat("\n  TOP 20 MOST CONSISTENTLY DE GENES\n")
print(master_df %>%
        filter(n_contrasts_sig > 0) %>%
        select(gene_id, n_contrasts_sig, phase) %>%
        slice_head(n = 20), n = 20)

cat("\n  k-MEANS CLUSTER SIZES\n")
print(table(clust_df$cluster))

cat("\n  OUTPUT FILES\n")
out_files <- list.files(
  path = file.path(cfg$results_dir, c("plots/volcano","plots/heatmaps",
                                      "plots/overlap","tables")),
  pattern = "DEG|volcano|heatmap|venn|upset|barplot",
  full.names = TRUE, recursive = FALSE)
invisible(lapply(sort(out_files), function(f) cat("  ", f, "\n")))
cat("\n══════════════════════════════════════════════════════════════════════\n")
log_msg("Section 2 (DEG Analysis) complete.")


# ==============================================================================

options(stringsAsFactors = FALSE)
set.seed(20260505)

# ── Configuration (same as parent scripts) ──────────────────────────────────
cfg <- list(
  project_name  = "rice_salt_timeseries_10mM_NaCl",
  results_dir   = "results",
  alpha         = 0.05,
  lfc_cutoff    = 1,
  baseline_time = "T0",
  time_levels   = c("T0", "T60", "T180", "T240", "T5d"),
  time_labels   = c(T0 = "0 min", T60 = "60 min", T180 = "180 min",
                    T240 = "240 min", T5d = "5 days"),
  time_minutes  = c(T0 = 0, T60 = 60, T180 = 180, T240 = 240, T5d = 7200),
  auto_install  = TRUE
)

# ── Output directories ──────────────────────────────────────────────────────
dirs <- file.path(cfg$results_dir,
                  c("tables", "plots/temporal", "plots/wgcna", "plots/tf"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# ── Logging & helpers ───────────────────────────────────────────────────────
log_msg <- function(...) {
  cat(sprintf("[%s] %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              paste0(..., collapse = "")))
}

write_csv_safe <- function(x, path) {
  data.table::fwrite(as.data.frame(x), path)
}

save_plot <- function(plot, basename, width = 8, height = 6, dpi = 1200) {
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
  ggsave(paste0(basename, ".pdf"), plot, width = width, height = height, device = pdf_device)
  ggsave(paste0(basename, ".tiff"), plot, width = width, height = height, dpi = dpi,
         compression = "lzw", device = "tiff")
}

# ── Publication theme ───────────────────────────────────────────────────────
pub_theme <- function(base_size = 12) {
  ggpubr::theme_pubr(base_size = base_size, border = FALSE) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0, size = base_size + 2),
      plot.subtitle = ggplot2::element_text(color = "grey25"),
      axis.title    = ggplot2::element_text(face = "bold"),
      legend.title  = ggplot2::element_text(face = "bold"),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.25),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      strip.text    = ggplot2::element_text(face = "bold")
    )
}

# ── Colour palettes ─────────────────────────────────────────────────────────
time_cols <- c(T0 = "#4B5563", T60 = "#0072B2", T180 = "#009E73",
               T240 = "#E69F00", T5d = "#D55E00")
condition_cols <- c(Control = "#6B7280", Salt = "#DC2626")

# ── Package installation & loading ──────────────────────────────────────────
cran_pkgs <- c("tidyverse", "ggplot2", "pheatmap", "RColorBrewer",
               "ggpubr", "dplyr", "data.table", "ggrepel", "igraph",
               "ggnetwork", "reshape2", "scales", "patchwork",
               "BiocManager")
bioc_pkgs <- c("DESeq2", "Mfuzz", "maSigPro", "WGCNA", "ComplexHeatmap", 
               "circlize", "BiocParallel")

if (isTRUE(cfg$auto_install)) {
  missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace,
                                    logical(1), quietly = TRUE)]
  if (length(missing_cran) > 0) {
    install.packages(missing_cran, repos = "https://cloud.r-project.org")
  }
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  
  missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace,
                                    logical(1), quietly = TRUE)]
  if (length(missing_bioc) > 0) {
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  }
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(Mfuzz)
  library(maSigPro)
  library(WGCNA)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(ggpubr)
  library(data.table)
  library(scales)
  library(igraph)
  library(ggnetwork)
  library(patchwork)
})

has_complex_heatmap <- requireNamespace("ComplexHeatmap", quietly = TRUE) &&
  requireNamespace("circlize", quietly = TRUE)
if (has_complex_heatmap) {
  suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize) })
}

# ── WGCNA settings ──────────────────────────────────────────────────────────
allowWGCNAThreads(nThreads = parallel::detectCores() - 1)

# ==============================================================================
#  LOAD REQUIRED OBJECTS FROM PREVIOUS SECTIONS
# ==============================================================================

metadata_salt <- as.data.frame(colData(dds_salt))

# ── Robust Mfuzz clustering ────────────────────────────────────────────────

# Ensure clean data
log_msg("Cleaning data for Mfuzz...")

# Remove rows with any NA/Inf
keep_rows <- apply(salt_vst, 1, function(x) {
  all(!is.na(x) & !is.infinite(x) & is.finite(x))
})
salt_vst_clean <- salt_vst[keep_rows, ]
log_msg("  Kept ", sum(keep_rows), "/", nrow(salt_vst), " genes after NA/Inf removal")

# Remove constant rows
row_sd <- apply(salt_vst_clean, 1, sd, na.rm = TRUE)
keep_rows2 <- row_sd > 0.1  # Keep genes with some variation
salt_vst_clean <- salt_vst_clean[keep_rows2, ]
log_msg("  Kept ", sum(keep_rows2), "/", length(row_sd), " genes after low-variance filtering")

# Take top 5000 most variable genes if still too many
if(nrow(salt_vst_clean) > 5000) {
  row_var <- apply(salt_vst_clean, 1, var)
  top_genes <- order(row_var, decreasing = TRUE)[1:5000]
  salt_vst_clean <- salt_vst_clean[top_genes, ]
  log_msg("  Reduced to top 5000 most variable genes")
}

# Create ExpressionSet
expr_set <- Biobase::ExpressionSet(
  assayData = as.matrix(salt_vst_clean),
  phenoData = Biobase::AnnotatedDataFrame(metadata_salt)
)

# Standardize
expr_set_std <- standardise(expr_set)

# Use fixed m value (2.0 is standard)
m_value <- 2.0
log_msg("  Using m = ", m_value)

# Direct clustering with c=6 (skip optimization)
optimal_c <- 6
log_msg("Running Mfuzz with c = ", optimal_c, "...")

set.seed(20260505)
mfuzz_fit <- mfuzz(expr_set_std, centers = optimal_c, m = m_value, iter.max = 200)

log_msg("  Clustering complete!")

# Extract results
mfuzz_clusters <- mfuzz_fit$cluster
mfuzz_membership <- mfuzz_fit$membership

log_msg("  Cluster sizes:")
print(table(mfuzz_clusters))

# Create output dataframe
mfuzz_df <- tibble(
  gene_id = rownames(expr_set_std),
  mfuzz_cluster = as.numeric(mfuzz_clusters),
  membership_score = apply(mfuzz_membership, 1, max)
)

# Save
write_csv_safe(mfuzz_df, file.path(cfg$results_dir, "tables", "temporal_clusters_mfuzz.csv"))

log_msg("  ✔ Mfuzz clustering complete: ", optimal_c, " clusters identified")
# Extract cluster trajectories for visualization
cluster_trajectories <- as.data.frame(exprs(expr_set_std)) %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "expression") %>%
  left_join(mfuzz_df, by = "gene_id") %>%
  left_join(metadata_salt %>% rownames_to_column("sample"), by = "sample") %>%
  group_by(mfuzz_cluster, time_point) %>%
  summarise(
    mean_expr = mean(expression),
    sd_expr = sd(expression),
    se_expr = sd_expr / sqrt(n()),
    n_genes = n(),
    .groups = "drop"
  )

# Plot trajectories
library(ggplot2)

cluster_plot <- ggplot(cluster_trajectories, 
                       aes(x = factor(time_point, levels = cfg$time_levels), 
                           y = mean_expr, 
                           group = mfuzz_cluster,
                           color = factor(mfuzz_cluster))) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  geom_ribbon(aes(ymin = mean_expr - se_expr, 
                  ymax = mean_expr + se_expr, 
                  fill = factor(mfuzz_cluster)),
              alpha = 0.2, color = NA) +
  facet_wrap(~mfuzz_cluster, scales = "free_y") +
  scale_color_brewer(palette = "Set1", name = "Cluster") +
  scale_fill_brewer(palette = "Set1", name = "Cluster") +
  labs(title = "Mfuzz Temporal Clustering Patterns",
       x = "Time Point",
       y = "Mean Z-score Expression") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(cluster_plot, 
          file.path(cfg$results_dir, "plots/temporal_clustering", "cluster_trajectories"))


# Summary statistics for each cluster
cluster_summary <- mfuzz_df %>%
  group_by(mfuzz_cluster) %>%
  summarise(
    n_genes = n(),
    mean_membership = mean(membership_score),
    high_confidence_genes = sum(membership_score > cfg$membership_threshold),
    pct_high_confidence = 100 * high_confidence_genes / n(),
    .groups = "drop"
  )

print(cluster_summary)
write_csv_safe(cluster_summary, 
               file.path(cfg$results_dir, "tables", "cluster_summary_stats.csv"))
# Order genes by cluster and membership score
gene_order <- mfuzz_df %>%
  arrange(mfuzz_cluster, desc(membership_score)) %>%
  pull(gene_id)

# Order samples by time
sample_order <- metadata_salt %>%
  arrange(time_point) %>%
  rownames()

# Create heatmap matrix
heatmap_mat <- exprs(expr_set_std)[gene_order, sample_order]

# Annotations
col_anno <- data.frame(
  Time = metadata_salt[sample_order, "time_point"],
  row.names = sample_order
)

row_anno <- data.frame(
  Cluster = factor(mfuzz_df$mfuzz_cluster[match(gene_order, mfuzz_df$gene_id)]),
  Membership = mfuzz_df$membership_score[match(gene_order, mfuzz_df$gene_id)],
  row.names = gene_order
)

# Plot heatmap
# Create the heatmap and save directly
pdf(file.path(cfg$results_dir, "plots/temporal_clustering", "cluster_heatmap.pdf"),
    width = 8, height = 10)

pheatmap(heatmap_mat,
         annotation_col = col_anno,
         annotation_row = row_anno,
         show_rownames = FALSE,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         main = "Mfuzz Temporal Clusters (ordered by membership score)",
         fontsize = 8,
         fontsize_row = 4,
         annotation_legend = TRUE)

dev.off()

# Also save as TIFF
tiff(file.path(cfg$results_dir, "plots/temporal_clustering", "cluster_heatmap.tiff"),
    width = 8, height = 10, units = "in", res = 1200, compression = "lzw")

pheatmap(heatmap_mat,
         annotation_col = col_anno,
         annotation_row = row_anno,
         show_rownames = FALSE,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         main = "Mfuzz Temporal Clusters (ordered by membership score)",
         fontsize = 8,
         fontsize_row = 4,
         annotation_legend = TRUE)

dev.off()

log_msg("Heatmap saved to: ", file.path(cfg$results_dir, "plots/temporal_clustering", "cluster_heatmap.pdf"))

log_msg("=== MFUZZ CLUSTERING SUMMARY ===")
log_msg("Total clusters: ", length(unique(mfuzz_df$mfuzz_cluster)))
log_msg("Total genes clustered: ", nrow(mfuzz_df))
log_msg("Average membership score: ", round(mean(mfuzz_df$membership_score), 3))
log_msg("Genes with high confidence (membership > ", cfg$membership_threshold, "): ", 
        sum(mfuzz_df$membership_score > cfg$membership_threshold))
log_msg("")
log_msg("Cluster sizes:")
mfuzz_df %>% count(mfuzz_cluster) %>% arrange(mfuzz_cluster) %>% 
  mutate(percentage = round(100 * n / sum(n), 1)) %>% 
  print()

log_msg("  ✔ Mfuzz clustering complete: ", optimal_c, " clusters identified")



# ── Limma temporal trend analysis ──────────────────────────────────────────
log_msg("3.1.5: Performing temporal trend analysis with limma...")

library(limma)

# Prepare data
salt_vst_matrix <- salt_vst
time_points <- metadata_salt[colnames(salt_vst_matrix), "time_minutes"]

# Create design matrix with natural splines for flexible temporal patterns
library(splines)
design_limma <- model.matrix(~ ns(time_points, df = 3))  # 3 degrees of freedom for flexibility

# Fit linear model
fit_limma <- lmFit(salt_vst_matrix, design_limma)
fit_limma <- eBayes(fit_limma)

# Extract significant temporal genes
temporal_results <- data.frame(
  gene_id = rownames(fit_limma),
  F_stat = fit_limma$F,
  F_pvalue = fit_limma$F.p.value,
  adj_pvalue = p.adjust(fit_limma$F.p.value, method = "BH")
) %>%
  mutate(
    significant = adj_pvalue < 0.05,
    temporal_pattern = case_when(
      significant ~ "Significant temporal change",
      TRUE ~ "No significant change"
    )
  )

log_msg("  Found ", sum(temporal_results$significant), " genes with significant temporal changes")

# Get coefficients for each gene
coefficients <- fit_limma$coefficients
colnames(coefficients) <- c("Intercept", "ns1", "ns2", "ns3")

# Add to results
temporal_results <- cbind(temporal_results, coefficients)

# Save results
write_csv_safe(temporal_results,
               file.path(cfg$results_dir, "tables", "limma_temporal_results.csv"))

# Identify genes with specific patterns
temporal_results <- temporal_results %>%
  mutate(
    pattern_direction = case_when(
      ns1 > 0 & ns2 > 0 & ns3 > 0 ~ "Increasing",
      ns1 < 0 & ns2 < 0 & ns3 < 0 ~ "Decreasing",
      (ns1 > 0 & ns2 < 0) | (ns1 < 0 & ns2 > 0) ~ "Peak/Nadir at middle",
      TRUE ~ "Complex pattern"
    )
  )

log_msg("  ✔ Limma temporal analysis complete")

# ── Integrate Limma temporal results with Mfuzz clusters ────────────────────
log_msg("3.1.6: Integrating Limma results with Mfuzz clusters...")

# Ensure limma_results is available
if(exists("temporal_results")) {
  
  # Combine Mfuzz clusters with limma results
  integrated_results <- mfuzz_df %>%
    left_join(temporal_results %>% 
                select(gene_id, F_stat, F_pvalue, adj_pvalue, significant, pattern_direction, 
                       Intercept, ns1, ns2, ns3),
              by = "gene_id")
  
  # For each cluster, summarize temporal patterns
  cluster_pattern_summary <- integrated_results %>%
    group_by(mfuzz_cluster) %>%
    summarise(
      n_genes = n(),
      n_significant = sum(significant, na.rm = TRUE),
      pct_significant = 100 * n_significant / n(),
      mean_F_stat = mean(F_stat, na.rm = TRUE),
      mean_membership = mean(membership_score),
      .groups = "drop"
    )
  
  log_msg("\n=== Cluster Summary with Temporal Significance ===")
  print(cluster_pattern_summary)
  
  # Pattern distribution per cluster
  pattern_distribution <- integrated_results %>%
    filter(significant) %>%
    group_by(mfuzz_cluster, pattern_direction) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(mfuzz_cluster) %>%
    mutate(percentage = 100 * count / sum(count))
  
  log_msg("\n=== Pattern Direction per Cluster ===")
  print(pattern_distribution)
  
  # Save integrated results
  write_csv_safe(integrated_results,
                 file.path(cfg$results_dir, "tables", "mfuzz_limma_integrated.csv"))
  
  write_csv_safe(cluster_pattern_summary,
                 file.path(cfg$results_dir, "tables", "cluster_temporal_summary.csv"))
  
  write_csv_safe(pattern_distribution,
                 file.path(cfg$results_dir, "tables", "cluster_pattern_directions.csv"))
  
} else {
  log_msg("  Warning: temporal_results not found. Using ANOVA results instead.")
}

# Plot: Pattern direction distribution across clusters
pattern_plot <- pattern_distribution %>%
  ggplot(aes(x = factor(mfuzz_cluster), y = percentage, fill = pattern_direction)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "RdYlBu", name = "Temporal Pattern") +
  labs(title = "Distribution of Temporal Patterns by Cluster",
       subtitle = "For significantly changing genes only",
       x = "Mfuzz Cluster",
       y = "Percentage of significant genes") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

save_plot(pattern_plot,
          file.path(cfg$results_dir, "plots/temporal_clustering", 
                    "cluster_pattern_distribution"))
# Create volcano plots for each cluster
for(clust in unique(integrated_results$mfuzz_cluster)) {
  cluster_data <- integrated_results %>%
    filter(mfuzz_cluster == clust)
  
  # Define significance and logFC-like (using F-stat)
  cluster_data <- cluster_data %>%
    mutate(
      significance = ifelse(adj_pvalue < 0.05, "Significant", "Not Significant"),
      logF = log10(F_stat + 1)  # Similar to logFC but for F-statistic
    )
  
  volcano <- ggplot(cluster_data, aes(x = logF, y = -log10(adj_pvalue + 1e-10),
                                      color = significance)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    labs(title = paste("Cluster", clust, "- Temporal Significance"),
         subtitle = paste("n =", nrow(cluster_data), "genes"),
         x = "log10(F-statistic + 1)",
         y = "-log10(adjusted p-value)") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  save_plot(volcano,
            file.path(cfg$results_dir, "plots/temporal_clustering", 
                      paste0("cluster_", clust, "_volcano")))
}

log_msg("  Volcano plots saved for each cluster")
# Create heatmap of cluster expression patterns
library(pheatmap)

# Get top 50 most representative genes per cluster (highest membership score)
top_genes <- integrated_results %>%
  group_by(mfuzz_cluster) %>%
  arrange(desc(membership_score)) %>%
  slice_head(n = 50) %>%
  pull(gene_id)

# Prepare heatmap matrix
heatmap_mat <- salt_vst[top_genes, ]
rownames(heatmap_mat) <- top_genes

# Order samples by time
sample_order <- metadata_salt %>%
  arrange(time_point) %>%
  rownames()
heatmap_mat <- heatmap_mat[, sample_order]

# Create annotations
col_annotation <- data.frame(
  Time = metadata_salt[sample_order, "time_point"],
  row.names = sample_order
)

row_annotation <- integrated_results %>%
  filter(gene_id %in% top_genes) %>%
  arrange(mfuzz_cluster, desc(membership_score)) %>%
  select(gene_id, mfuzz_cluster, membership_score, significant) %>%
  column_to_rownames("gene_id")

# Annotation colors
ann_colors <- list(
  Time = c(T0 = "#4B5563", T60 = "#0072B2", T180 = "#009E73", 
           T240 = "#E69F00", T5d = "#D55E00"),
  mfuzz_cluster = setNames(brewer.pal(length(unique(row_annotation$mfuzz_cluster)), "Set1"),
                           sort(unique(row_annotation$mfuzz_cluster))),
  significant = c("TRUE" = "red", "FALSE" = "gray")
)

# Save heatmap
pdf(file.path(cfg$results_dir, "plots/temporal_clustering", 
              "representative_genes_heatmap.pdf"),
    width = 10, height = 12)

pheatmap(heatmap_mat,
         annotation_col = col_annotation,
         annotation_row = row_annotation[, c("mfuzz_cluster", "significant")],
         annotation_colors = ann_colors,
         show_rownames = FALSE,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale = "row",
         main = "Top Representative Genes per Mfuzz Cluster",
         fontsize = 8,
         border_color = NA)

dev.off()

tiff(file.path(cfg$results_dir, "plots/temporal_clustering", 
              "representative_genes_heatmap.tiff"),
    width = 10, height = 12, units = "in", res = 1200, compression = "lzw")
pheatmap(heatmap_mat,
         annotation_col = col_annotation,
         annotation_row = row_annotation[, c("mfuzz_cluster", "significant")],
         annotation_colors = ann_colors,
         show_rownames = FALSE,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale = "row",
         main = "Top Representative Genes per Mfuzz Cluster",
         fontsize = 8,
         border_color = NA)

dev.off()

log_msg("  Representative genes heatmap saved")
# Plot mean expression trajectories per cluster with confidence bands
cluster_trajectories <- integrated_results %>%
  # Get expression data for each gene
  left_join(salt_vst %>% as.data.frame() %>% rownames_to_column("gene_id"),
            by = "gene_id") %>%
  pivot_longer(cols = starts_with("SALT_"),  # Adjust column pattern as needed
               names_to = "sample", values_to = "expression") %>%
  left_join(metadata_salt %>% rownames_to_column("sample"), by = "sample") %>%
  group_by(mfuzz_cluster, time_point) %>%
  summarise(
    mean_expr = mean(expression, na.rm = TRUE),
    se_expr = sd(expression, na.rm = TRUE) / sqrt(n()),
    n_genes = n_distinct(gene_id),
    n_significant = sum(significant, na.rm = TRUE),
    .groups = "drop"
  )

# Create trajectory plot
trajectory_plot <- cluster_trajectories %>%
  ggplot(aes(x = factor(time_point, levels = cfg$time_levels), 
             y = mean_expr, 
             color = factor(mfuzz_cluster),
             group = mfuzz_cluster)) +
  geom_line(size = 1.2) +
  geom_point(size = 2.5) +
  geom_ribbon(aes(ymin = mean_expr - se_expr, 
                  ymax = mean_expr + se_expr,
                  fill = factor(mfuzz_cluster)),
              alpha = 0.2, color = NA) +
  facet_wrap(~mfuzz_cluster, scales = "free_y", ncol = 3) +
  scale_color_brewer(palette = "Set1", name = "Cluster") +
  scale_fill_brewer(palette = "Set1", name = "Cluster") +
  labs(title = "Mfuzz Cluster Expression Trajectories",
       subtitle = paste("Mean ± SE (", nrow(integrated_results), " genes)"),
       x = "Time Point",
       y = "Mean VST Expression") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(trajectory_plot,
          file.path(cfg$results_dir, "plots/temporal_clustering", 
                    "cluster_trajectories_with_CI"),
          width = 12, height = 8)

log_msg("  Trajectory plot saved")
# Generate comprehensive summary
sink(file.path(cfg$results_dir, "tables", "temporal_analysis_summary.txt"))

cat("=== TEMPORAL CLUSTERING ANALYSIS SUMMARY ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("MFUZZ CLUSTERING:\n")
cat("  Total genes clustered:", nrow(mfuzz_df), "\n")
cat("  Number of clusters:", length(unique(mfuzz_df$mfuzz_cluster)), "\n")
cat("  Mean membership score:", round(mean(mfuzz_df$membership_score), 3), "\n")
cat("  High confidence genes (membership >", cfg$membership_threshold, "):", 
    sum(mfuzz_df$membership_score > cfg$membership_threshold), "\n\n")

cat("CLUSTER SIZES:\n")
print(table(mfuzz_df$mfuzz_cluster))

cat("\nTEMPORAL SIGNIFICANCE (Limma):\n")
cat("  Total significant genes:", sum(temporal_results$significant, na.rm = TRUE), "\n")
cat("  Percentage of total:", 
    round(100 * sum(temporal_results$significant, na.rm = TRUE) / nrow(temporal_results), 1), "%\n\n")

cat("SIGNIFICANCE PER CLUSTER:\n")
print(cluster_pattern_summary)

cat("\nPATTERN DISTRIBUTION:\n")
print(pattern_distribution)

sink()

log_msg("  Summary report saved to: temporal_analysis_summary.txt")
log_msg("  ✔ Integration complete!")

# ── Save all genes in each temporal cluster ──────────────────────────────────
log_msg("Saving gene lists for each Mfuzz cluster...")

# Ensure mfuzz_df exists
if(!exists("mfuzz_df")) {
  log_msg("  Loading Mfuzz results...")
  mfuzz_df <- read_csv(file.path(cfg$results_dir, "tables", "temporal_clusters_mfuzz.csv"))
}

# Create cluster directory
cluster_dir <- file.path(cfg$results_dir, "tables", "cluster_genes")
dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)

# Get unique clusters
clusters <- sort(unique(mfuzz_df$mfuzz_cluster))

# Save genes for each cluster
cluster_summary <- data.frame()

for(clust in clusters) {
  # Get genes in this cluster
  cluster_genes <- mfuzz_df %>%
    filter(mfuzz_cluster == clust) %>%
    arrange(desc(membership_score))
  
  # Extract just gene IDs
  gene_ids <- cluster_genes$gene_id
  
  # Save full data frame with membership scores
  write_csv_safe(cluster_genes,
                 file.path(cluster_dir, paste0("cluster_", clust, "_all_genes.csv")))
  
  # Save simple gene list (one per line)
  writeLines(gene_ids,
             file.path(cluster_dir, paste0("cluster_", clust, "_genes.txt")))
  
  # Save high-confidence genes (membership > 0.7)
  high_conf <- cluster_genes %>%
    filter(membership_score > 0.7)
  
  if(nrow(high_conf) > 0) {
    write_csv_safe(high_conf,
                   file.path(cluster_dir, paste0("cluster_", clust, "_high_confidence_genes.csv")))
    
    writeLines(high_conf$gene_id,
               file.path(cluster_dir, paste0("cluster_", clust, "_high_confidence_genes.txt")))
  }
  
  # Summary for this cluster
  cluster_summary <- rbind(cluster_summary, data.frame(
    cluster = clust,
    total_genes = nrow(cluster_genes),
    high_confidence_genes = nrow(high_conf),
    min_membership = round(min(cluster_genes$membership_score), 3),
    max_membership = round(max(cluster_genes$membership_score), 3),
    mean_membership = round(mean(cluster_genes$membership_score), 3)
  ))
  
  log_msg("  Cluster ", clust, ": ", nrow(cluster_genes), 
          " genes (", nrow(high_conf), " high confidence)")
}

# Save cluster summary
write_csv_safe(cluster_summary,
               file.path(cluster_dir, "cluster_summary.csv"))

log_msg("\nCluster Summary:")
print(cluster_summary)

# Save all genes in a single file with cluster annotation
all_genes_with_clusters <- mfuzz_df %>%
  arrange(mfuzz_cluster, desc(membership_score))

write_csv_safe(all_genes_with_clusters,
               file.path(cfg$results_dir, "tables", "all_genes_with_cluster_assignments.csv"))

log_msg("\n  All gene lists saved to: ", cluster_dir)


# ── Save cluster statistics and metadata ────────────────────────────────────
log_msg("Saving cluster statistics...")

# Calculate cluster-specific statistics
cluster_stats <- mfuzz_df %>%
  group_by(mfuzz_cluster) %>%
  summarise(
    n_genes = n(),
    mean_membership = mean(membership_score),
    sd_membership = sd(membership_score),
    median_membership = median(membership_score),
    min_membership = min(membership_score),
    max_membership = max(membership_score),
    q25_membership = quantile(membership_score, 0.25),
    q75_membership = quantile(membership_score, 0.75),
    .groups = "drop"
  )

write_csv_safe(cluster_stats,
               file.path(cfg$results_dir, "tables", "cluster_statistics.csv"))

# For each cluster, calculate expression summary
if(exists("salt_vst")) {
  cluster_expr_summary <- data.frame()
  
  for(clust in clusters) {
    cluster_genes <- mfuzz_df %>% filter(mfuzz_cluster == clust) %>% pull(gene_id)
    cluster_expr <- salt_vst[cluster_genes, ]
    
    # Calculate mean expression per time point
    time_points <- colnames(cluster_expr)
    mean_expr <- rowMeans(cluster_expr)
    mean_profile <- colMeans(cluster_expr)
    
    # Create summary
    cluster_expr_summary <- rbind(cluster_expr_summary, data.frame(
      cluster = clust,
      n_genes = length(cluster_genes),
      mean_expression_all = mean(mean_expr),
      expression_range_min = min(cluster_expr),
      expression_range_max = max(cluster_expr),
      mean_expression_T0 = mean(cluster_expr[, grep("T0", colnames(cluster_expr))]),
      mean_expression_T60 = mean(cluster_expr[, grep("T60", colnames(cluster_expr))]),
      mean_expression_T180 = mean(cluster_expr[, grep("T180", colnames(cluster_expr))]),
      mean_expression_T240 = mean(cluster_expr[, grep("T240", colnames(cluster_expr))]),
      mean_expression_T5d = mean(cluster_expr[, grep("T5d", colnames(cluster_expr))]),
      fold_change_T5d_vs_T0 = mean(cluster_expr[, grep("T5d", colnames(cluster_expr))]) / 
        mean(cluster_expr[, grep("T0", colnames(cluster_expr))])
    ))
  }
  
  write_csv_safe(cluster_expr_summary,
                 file.path(cfg$results_dir, "tables", "cluster_expression_summary.csv"))
}

log_msg("  Cluster statistics saved")



