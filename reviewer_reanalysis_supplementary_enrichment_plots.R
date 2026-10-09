###############################################################################
## reviewer_reanalysis_supplementary_enrichment_plots.R
##
## Supplementary enrichment figures (PNG) for the reviewer-corrected reanalysis.
##
## This script is purely presentational: it reads the enrichment tables that
## were produced by reviewer_reanalysis.R and renders publication-quality PNG
## figures into results_reviewer_reanalysis/supplementary/figures/.
##
## No statistical result is recalculated or modified here.  Significance is
## always the Benjamini-Hochberg FDR reported in the input tables (FDR < 0.05).
##
## Inputs
##   tables/enrichment_gprofiler_GO_KEGG_significant_FDR0.05.csv
##   tables/enrichment_gprofiler_GO_KEGG_all_terms.csv
##   tables/enrichment_local_hypergeometric_GO_all_terms.csv
##   tables/enrichment_local_hypergeometric_GO_FDR0.05.csv
##
## Outputs (PNG, 300 dpi)
##   supplementary/figures/figureS_enrichment_overview_per_set.png
##   supplementary/figures/figureS_enrichment_<gene set>.png          (11 files)
##   supplementary/figures/figureS_enrichment_term_heatmap.png
##   supplementary/figures/figureS_enrichment_significant_terms_by_source.png
##   supplementary/figures/figureS_enrichment_GO_by_ontology.png
##   supplementary/figures/figureS_enrichment_KEGG_pathways.png
##   supplementary/figures/figureS_enrichment_local_GO_crosscheck.png
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

## ---------------------------------------------------------------------------
## 0. Paths, palette and helpers
## ---------------------------------------------------------------------------

BASE_DIR <- "f:/independent projects/oryza sativa/comments"
if (!dir.exists(file.path(BASE_DIR, "results_reviewer_reanalysis"))) {
  BASE_DIR <- normalizePath(".")
}
OUT_DIR <- file.path(BASE_DIR, "results_reviewer_reanalysis")
TAB_DIR <- file.path(OUT_DIR, "tables")
SUP_FIG <- file.path(OUT_DIR, "supplementary", "figures")
dir.create(SUP_FIG, showWarnings = FALSE, recursive = TRUE)

LOG <- file.path(SUP_FIG, "supplementary_enrichment_plots_log.txt")
log_con <- file(LOG, open = "wt")
log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  cat(msg, "\n", file = log_con, append = TRUE)
}
log_msg("# supplementary enrichment PNG generation")
log_msg("started: ", format(Sys.time()))

SRC_COLOURS <- c("GO:BP" = "#1b9e77", "GO:MF" = "#7570b3",
                 "GO:CC" = "#e6ab02", "KEGG" = "#d95f02")
SRC_LABELS <- c("GO:BP" = "GO biological process", "GO:MF" = "GO molecular function",
                "GO:CC" = "GO cellular component", "KEGG" = "KEGG pathway")
FDR_CUT <- 0.05

## gene-set ordering and human-readable titles
SET_ORDER <- c("DEG_Salt_vs_Control_T0", "DEG_Salt_vs_Control_T60",
               "DEG_Salt_vs_Control_T180", "DEG_Salt_vs_Control_T240",
               "DEG_Salt_vs_Control_T5d", "LRT_interaction_significant",
               "Mfuzz_cluster_1", "Mfuzz_cluster_2", "Mfuzz_cluster_3",
               "Mfuzz_cluster_4", "Mfuzz_cluster_5")

set_title <- function(s) {
  out <- sub("^DEG_Salt_vs_Control_(.+)$", "Salt vs Control at \\1 (time-matched DEGs)", s)
  out <- sub("^LRT_interaction_significant$",
             "Genes with a significant salt x time interaction (LRT)", out)
  out <- sub("^Mfuzz_cluster_(\\d+)$", "Mfuzz cluster \\1 (salt-effect profiles)", out)
  out
}

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
          strip.background = element_rect(fill = "grey92", colour = "grey40"),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", size = base_size + 1),
          plot.subtitle = element_text(size = base_size - 1, colour = "grey25"),
          plot.caption = element_text(size = base_size - 2, colour = "grey35", hjust = 0),
          legend.position = "bottom")
}

wrap_term <- function(x, width = 52) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

save_png <- function(plot, file, width = 9, height = 7, dpi = 300) {
  ok <- tryCatch({
    ggsave(file, plot = plot, width = width, height = height, dpi = dpi,
           bg = "white", limitsize = FALSE)
    TRUE
  }, error = function(e) {
    log_msg("  !! failed: ", basename(file), " -> ", conditionMessage(e))
    FALSE
  })
  if (ok) log_msg("  wrote ", basename(file),
                  " (", round(file.info(file)$size / 1024), " kB)")
  invisible(ok)
}

## ---------------------------------------------------------------------------
## 1. Load enrichment results
## ---------------------------------------------------------------------------

sig_file <- file.path(TAB_DIR, "enrichment_gprofiler_GO_KEGG_significant_FDR0.05.csv")
all_file <- file.path(TAB_DIR, "enrichment_gprofiler_GO_KEGG_all_terms.csv")
if (!file.exists(sig_file) || !file.exists(all_file)) {
  stop("g:Profiler enrichment tables not found; run reviewer_reanalysis.R first")
}

sig <- fread(sig_file)
allt <- fread(all_file)
sig[, gene_set := factor(gene_set, levels = SET_ORDER)]
allt[, gene_set := factor(gene_set, levels = SET_ORDER)]
sig <- sig[!is.na(p_value)]
setorder(sig, p_value)
log_msg("significant terms (FDR < ", FDR_CUT, "): ", nrow(sig),
        " over ", length(unique(sig$gene_set)), " gene sets")

## number of query genes per set (from the tables themselves)
nq <- unique(allt[, .(gene_set, n_query_genes)])[order(gene_set)]
log_msg("query sizes: ", paste(sprintf("%s=%d", nq$gene_set, nq$n_query_genes), collapse = ", "))

## ---------------------------------------------------------------------------
## 2. Individual per-gene-set dot plots (top 15 terms by FDR, significance shown)
## ---------------------------------------------------------------------------

make_set_plot <- function(gs, n_show = 15) {
  d_all <- allt[gene_set == gs][order(p_value)][seq_len(min(n_show, .N))]
  if (!nrow(d_all)) return(NULL)
  d_all <- copy(d_all)
  d_all[, sig_flag := ifelse(!is.na(p_value) & p_value < FDR_CUT, "FDR < 0.05", "not significant")]
  d_all[, term_lab := wrap_term(term_name)]
  d_all[, term_lab := factor(term_lab, levels = rev(unique(term_lab)))]
  d_all[, neglog10 := -log10(pmax(p_value, .Machine$double.xmin))]
  n_sig_set <- nrow(sig[gene_set == gs])
  ggplot(d_all, aes(neglog10, term_lab)) +
    geom_segment(aes(x = 0, xend = neglog10, y = term_lab, yend = term_lab),
                 colour = "grey85", linewidth = 0.4) +
    geom_point(aes(colour = source, size = intersection_size, shape = sig_flag)) +
    geom_text(aes(label = paste0(intersection_size, "/", term_size)),
              hjust = -0.25, size = 2.5, colour = "grey25") +
    geom_vline(xintercept = -log10(FDR_CUT), linetype = 2, colour = "grey40", linewidth = 0.4) +
    scale_colour_manual(values = SRC_COLOURS, labels = SRC_LABELS, name = NULL) +
    scale_shape_manual(values = c("FDR < 0.05" = 16, "not significant" = 1), name = NULL) +
    scale_size_continuous(range = c(1.5, 6), name = "genes in term") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.18))) +
    labs(title = set_title(as.character(gs)),
         subtitle = paste0(nq[gene_set == gs]$n_query_genes, " query genes; ",
                           n_sig_set, " terms with FDR < 0.05 (of ",
                           nrow(allt[gene_set == gs]), " tested); ",
                           "top ", nrow(d_all), " terms by FDR shown"),
         x = expression(-log[10]~"(BH-adjusted P value)"), y = NULL,
         caption = paste0("labels: genes in term / genes annotated to the term; ",
                          "dashed line = FDR 0.05; g:Profiler (Oryza sativa)")) +
    theme_pub(10) +
    theme(legend.position = "right",
          axis.text.y = element_text(size = 7.5, lineheight = 0.95))
}

for (gs in SET_ORDER) {
  p <- make_set_plot(gs)
  if (is.null(p)) { log_msg("  no enrichment rows for ", gs); next }
  h <- max(5, 0.36 * min(15, nrow(allt[gene_set == gs])) + 2.2)
  save_png(p, file.path(SUP_FIG, paste0("figureS_enrichment_", gs, ".png")),
           width = 9.5, height = h)
}

## ---------------------------------------------------------------------------
## 3. Combined overview: top 6 terms per gene set in one multi-panel figure
## ---------------------------------------------------------------------------

top_ov <- sig[, .SD[order(p_value)][seq_len(min(6, .N))], by = gene_set]
top_ov[, term_lab := wrap_term(term_name, 44)]
top_ov[, panel := factor(as.character(gene_set), levels = SET_ORDER)]
top_ov[, term_lab := factor(paste(as.character(panel), term_lab, sep = "___"),
                            levels = rev(unique(paste(as.character(panel), term_lab, sep = "___"))))]
ov <- ggplot(top_ov, aes(-log10(p_value), term_lab)) +
  geom_col(aes(fill = source), width = 0.72) +
  geom_text(aes(label = paste0(intersection_size, " genes")), hjust = -0.15, size = 2.3,
            colour = "grey20") +
  geom_vline(xintercept = -log10(FDR_CUT), linetype = 2, colour = "grey30", linewidth = 0.4) +
  facet_wrap(~ panel, scales = "free_y", ncol = 2,
             labeller = labeller(panel = function(x) set_title(x))) +
  scale_fill_manual(values = SRC_COLOURS, labels = SRC_LABELS, name = NULL) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.22))) +
  labs(title = "Gene Ontology and KEGG enrichment of salt-responsive gene sets (top 6 terms per set)",
       subtitle = paste0("g:Profiler, Oryza sativa, GO BP/MF/CC + KEGG, BH-FDR; ",
                         "dashed line = FDR 0.05; number of genes in the intersection shown next to each bar"),
       x = expression(-log[10]~"(BH-adjusted P value)"), y = NULL) +
  theme_pub(9) +
  theme(axis.text.y = element_text(size = 7, lineheight = 0.95),
        legend.position = "bottom")
save_png(ov, file.path(SUP_FIG, "figureS_enrichment_overview_per_set.png"),
         width = 15, height = 16)

## ---------------------------------------------------------------------------
## 4. Heatmap: terms significant in several gene sets
## ---------------------------------------------------------------------------

## terms that are significant in >= 2 gene sets (union of the top 5 terms per set)
key_terms <- unique(sig[, .SD[order(p_value)][seq_len(min(5, .N))], by = gene_set]$term_id)
hm <- allt[term_id %in% key_terms,
           .(gene_set, term_id, term_name, source, p_value, intersection_size, term_size)]
n_sets_sig <- hm[, .(n_sig = sum(!is.na(p_value) & p_value < FDR_CUT)), by = term_id]
keep <- n_sets_sig[n_sig >= 2]$term_id
hm <- hm[term_id %in% keep]
hm[, neglog10 := ifelse(is.na(p_value), NA_real_,
                        -log10(pmax(pmin(p_value, 1), .Machine$double.xmin)))]
hm[, neglog10 := pmin(neglog10, 15)]
term_order <- hm[, .(best = min(p_value, na.rm = TRUE)), by = .(term_id, term_name, source)][order(best)]
hm[, term_lab := factor(wrap_term(term_name, 48),
                        levels = rev(wrap_term(term_order$term_name, 48)))]
hm[, lab := ifelse(!is.na(neglog10) & !is.na(p_value) & p_value < FDR_CUT,
                   as.character(intersection_size), "")]
log_msg("heatmap terms (significant in >= 2 gene sets): ", length(unique(hm$term_id)))

p_hm <- ggplot(hm, aes(gene_set, term_lab, fill = neglog10)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = lab), size = 2.2, colour = "grey10") +
  scale_fill_viridis_c(option = "D", na.value = "grey90", name = expression(-log[10]~"FDR"),
                       limits = c(0, 15)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", sub("^DEG Salt vs Control ", "DEG ", x))) +
  labs(title = "Enrichment across gene sets: terms significant in at least two gene sets",
       subtitle = "cell values = number of genes in the intersection; grey = not significant (FDR >= 0.05) or not tested; colour capped at 15",
       x = NULL, y = NULL) +
  theme_pub(9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(size = 7, lineheight = 0.95),
        legend.position = "right",
        panel.grid.major.y = element_blank())
save_png(p_hm, file.path(SUP_FIG, "figureS_enrichment_term_heatmap.png"),
         width = 12, height = max(6, 0.32 * length(unique(hm$term_id)) + 3.5))

## ---------------------------------------------------------------------------
## 5. Significant terms per gene set and source
## ---------------------------------------------------------------------------

ct <- sig[, .N, by = .(gene_set, source)]
ct[, gene_set := factor(as.character(gene_set), levels = rev(SET_ORDER))]
p_ct <- ggplot(ct, aes(N, gene_set, fill = source)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = N), hjust = -0.2, size = 2.7, colour = "grey20") +
  scale_fill_manual(values = SRC_COLOURS, labels = SRC_LABELS, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  scale_y_discrete(labels = function(x) set_title(x)) +
  labs(title = "Number of significantly enriched terms per gene set",
       subtitle = paste0("g:Profiler (Oryza sativa), BH-FDR < ", FDR_CUT,
                         "; total = ", nrow(sig), " significant gene-set x term results"),
       x = "significant terms", y = NULL) +
  theme_pub(10) + theme(legend.position = "bottom")
save_png(p_ct, file.path(SUP_FIG, "figureS_enrichment_significant_terms_by_source.png"),
         width = 11, height = 6.5)

## ---------------------------------------------------------------------------
## 6. GO terms by ontology (top terms per ontology across all gene sets)
## ---------------------------------------------------------------------------

go <- copy(sig[source %in% c("GO:BP", "GO:MF", "GO:CC")])
go[, ontology := factor(source, levels = c("GO:BP", "GO:MF", "GO:CC"))]
top_go <- go[, .SD[order(p_value)][seq_len(min(10, .N))], by = ontology]
top_go[, term_lab := wrap_term(term_name, 46)]
top_go[, key := paste(ontology, term_lab, sep = "___")]
top_go[, key := factor(key, levels = rev(unique(key)))]
p_go <- ggplot(top_go, aes(-log10(p_value), key, colour = gene_set)) +
  geom_segment(aes(x = 0, xend = -log10(p_value), y = key, yend = key),
               colour = "grey88", linewidth = 0.4) +
  geom_point(aes(size = intersection_size)) +
  geom_vline(xintercept = -log10(FDR_CUT), linetype = 2, colour = "grey35", linewidth = 0.4) +
  facet_wrap(~ ontology, scales = "free_y", ncol = 3,
             labeller = labeller(ontology = c("GO:BP" = "GO biological process",
                                              "GO:MF" = "GO molecular function",
                                              "GO:CC" = "GO cellular component"))) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_size_continuous(range = c(1.5, 6), name = "genes in term") +
  labs(title = "Top Gene Ontology terms across salt-responsive gene sets",
       subtitle = paste0("g:Profiler, Oryza sativa, BH-FDR; the 10 most significant terms per ",
                         "ontology are shown, coloured by the gene set in which they were found"),
       x = expression(-log[10]~"(BH-adjusted P value)"), y = NULL, colour = "gene set") +
  theme_pub(9) +
  theme(axis.text.y = element_text(size = 7, lineheight = 0.95),
        legend.position = "bottom", legend.box = "vertical")
save_png(p_go, file.path(SUP_FIG, "figureS_enrichment_GO_by_ontology.png"),
         width = 16, height = 9)

## ---------------------------------------------------------------------------
## 7. All significant KEGG pathways across gene sets
## ---------------------------------------------------------------------------

kegg <- copy(sig[source == "KEGG"])
if (nrow(kegg)) {
  kegg[, term_lab := wrap_term(term_name, 48)]
  agg <- kegg[, .(neglog10 = max(-log10(pmax(p_value, .Machine$double.xmin))),
                  n_genes = max(intersection_size),
                  n_sets = uniqueN(gene_set),
                  sets = paste(sort(unique(as.character(gene_set))), collapse = ", ")),
              by = term_lab]
  agg[, term_lab := factor(term_lab, levels = agg[order(neglog10)]$term_lab)]
  p_kegg <- ggplot(agg, aes(neglog10, term_lab)) +
    geom_col(aes(fill = neglog10), width = 0.72) +
    geom_text(aes(label = paste0(n_genes, " genes | ", n_sets, " set",
                                 ifelse(n_sets > 1, "s", ""))),
              hjust = -0.12, size = 2.5, colour = "grey20") +
    geom_vline(xintercept = -log10(FDR_CUT), linetype = 2, colour = "grey30", linewidth = 0.4) +
    scale_fill_viridis_c(option = "C", guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.30))) +
    labs(title = "Significantly enriched KEGG pathways in salt-responsive gene sets",
         subtitle = paste0(nrow(kegg), " significant gene-set x pathway results over ",
                           uniqueN(kegg$term_id), " distinct pathways; g:Profiler, BH-FDR < ",
                           FDR_CUT),
         x = expression(-log[10]~"(BH-adjusted P value)"), y = NULL,
         caption = "labels: largest number of query genes in the pathway across gene sets | number of gene sets") +
    theme_pub(10) +
    theme(axis.text.y = element_text(size = 8, lineheight = 0.95))
  save_png(p_kegg, file.path(SUP_FIG, "figureS_enrichment_KEGG_pathways.png"),
           width = 11, height = max(5, 0.34 * nrow(agg) + 3))
} else {
  log_msg("  no significant KEGG pathways")
}

## ---------------------------------------------------------------------------
## 8. Offline GO cross-check (org.Osativa.eg.db) with coverage caveat
## ---------------------------------------------------------------------------

loc_file <- file.path(TAB_DIR, "enrichment_local_hypergeometric_GO_all_terms.csv")
if (file.exists(loc_file)) {
  loc <- fread(loc_file)
  loc[, gene_set := factor(gene_set, levels = SET_ORDER)]
  loc[, neglog10 := -log10(pmax(p.adjust, .Machine$double.xmin))]
  loc[, sig_flag := ifelse(p.adjust < FDR_CUT, "FDR < 0.05", "not significant")]
  loc_top <- loc[order(pvalue)][seq_len(min(20, .N))]
  loc_top[, lab := paste0(wrap_term(description, 40), " [", ontology, "]")]
  loc_top[, lab := factor(lab, levels = rev(unique(lab)))]
  p_loc <- ggplot(loc_top, aes(neglog10, lab)) +
    geom_segment(aes(x = 0, xend = neglog10, y = lab, yend = lab),
                 colour = "grey86", linewidth = 0.4) +
    geom_point(aes(size = intersection_size, shape = sig_flag, colour = gene_set)) +
    geom_vline(xintercept = -log10(FDR_CUT), linetype = 2, colour = "grey40", linewidth = 0.4) +
    scale_shape_manual(values = c("FDR < 0.05" = 16, "not significant" = 1), name = NULL) +
    scale_size_continuous(range = c(1.5, 6), name = "genes in term") +
    labs(title = "GO over-representation test with the local org.Osativa.eg.db annotation",
         subtitle = paste0("hypergeometric test, BH-FDR, universe = all 25,078 ",
                           "retained genes; top 20 gene-set x term results shown"),
         x = expression(-log[10]~"(BH-adjusted P value)"), y = NULL, colour = "gene set",
         caption = paste0("Note: the locally installed org.Osativa.eg.db (v0.1) contains only 95 GO terms for ",
                          "18,098 of 25,078 genes, so this test has lower resolution than the g:Profiler ",
                          "analysis in the other enrichment figures. It is shown as a control using the ",
                          "annotation that is available locally.")) +
    theme_pub(9) +
    theme(axis.text.y = element_text(size = 7, lineheight = 0.95),
          legend.position = "bottom", legend.box = "vertical",
          plot.caption = element_text(size = 7, colour = "grey35", hjust = 0))
  save_png(p_loc, file.path(SUP_FIG, "figureS_enrichment_local_GO_crosscheck.png"),
           width = 11, height = 8)
} else {
  log_msg("  local hypergeometric table not found - skipped")
}

## ---------------------------------------------------------------------------
## 9. Session info of this plotting run
## ---------------------------------------------------------------------------

writeLines(c(paste0("generated: ", format(Sys.time())),
             "presentational script: reviewer_reanalysis_supplementary_enrichment_plots.R",
             "input tables: enrichment_gprofiler_GO_KEGG_{all_terms,significant_FDR0.05}.csv,",
             "              enrichment_local_hypergeometric_GO_*.csv",
             capture.output(sessionInfo())),
           file.path(SUP_FIG, "supplementary_enrichment_plots_sessionInfo.txt"))

## manifest of the PNG files produced
pngs <- list.files(SUP_FIG, pattern = "^figureS_enrichment.*\\.png$", full.names = TRUE)
manifest <- data.table(file = basename(pngs),
                       size_kB = round(file.info(pngs)$size / 1024, 1),
                       generated = format(file.info(pngs)$mtime, "%Y-%m-%d %H:%M:%S"))
fwrite(manifest, file.path(SUP_FIG, "supplementary_enrichment_figures_manifest.csv"))
log_msg("\n", nrow(manifest), " PNG enrichment figures written:")
for (i in seq_len(nrow(manifest))) {
  log_msg("  ", manifest$file[i], " (", manifest$size_kB[i], " kB)")
}
log_msg("DONE: ", format(Sys.time()))
close(log_con)
