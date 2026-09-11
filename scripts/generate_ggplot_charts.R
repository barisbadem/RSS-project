# =============================================================================
# D-REGION allele frequency by geography - publication-style bar charts
#
# DATA SOURCE: raw observed haplotype counts from our own analysis of
# KIARVA's real, public 1000 Genomes IGHD genotype data (CC BY-NC 4.0). NOT
# copied from any published article. Same counts as the accompanying Excel
# workbook and GraphPad Prism export (single source, cross-checked).
#
# Everything below this one raw table is COMPUTED, not restated: R itself
# builds the per-gene contingency tables, runs the chi-square test of
# independence and the Benjamini-Hochberg correction across all 16 tested
# genes, derives percentages and 95% Wilson confidence intervals, and only
# then plots - nothing is pre-summarized or hardcoded per gene/allele.
# =============================================================================

library(ggplot2)

raw_csv <- "gene,allele,region,count
IGHD2-2,IGHD2-2*01,Africa,454
IGHD2-2,IGHD2-2*01,Europe,327
IGHD2-2,IGHD2-2*01,Asia,603
IGHD2-2,IGHD2-2*01_S7032,Africa,4
IGHD2-2,IGHD2-2*01_S7032,Europe,0
IGHD2-2,IGHD2-2*01_S7032,Asia,0
IGHD2-2,IGHD2-2*02,Africa,309
IGHD2-2,IGHD2-2*02,Europe,409
IGHD2-2,IGHD2-2*02,Asia,200
IGHD2-2,IGHD2-2*04_S0329,Africa,1
IGHD2-2,IGHD2-2*04_S0329,Europe,11
IGHD2-2,IGHD2-2*04_S0329,Asia,5
IGHD2-2,IGHD2-2*05_S3769,Africa,28
IGHD2-2,IGHD2-2*05_S3769,Europe,1
IGHD2-2,IGHD2-2*05_S3769,Asia,0
IGHD2-21,IGHD2-21*01,Africa,279
IGHD2-21,IGHD2-21*01,Europe,92
IGHD2-21,IGHD2-21*01,Asia,140
IGHD2-21,IGHD2-21*02,Africa,299
IGHD2-21,IGHD2-21*02,Europe,318
IGHD2-21,IGHD2-21*02,Asia,444
IGHD2-8,IGHD2-8*01,Africa,531
IGHD2-8,IGHD2-8*01,Europe,482
IGHD2-8,IGHD2-8*01,Asia,448
IGHD2-8,IGHD2-8*02,Africa,79
IGHD2-8,IGHD2-8*02,Europe,80
IGHD2-8,IGHD2-8*02,Asia,82
IGHD3-10,IGHD3-10*01,Africa,653
IGHD3-10,IGHD3-10*01,Europe,599
IGHD3-10,IGHD3-10*01,Asia,682
IGHD3-10,IGHD3-10*01_S3902,Africa,0
IGHD3-10,IGHD3-10*01_S3902,Europe,4
IGHD3-10,IGHD3-10*01_S3902,Asia,0
IGHD3-10,IGHD3-10*03,Africa,51
IGHD3-10,IGHD3-10*03,Europe,81
IGHD3-10,IGHD3-10*03,Asia,54
IGHD3-16,IGHD3-16*02,Africa,688
IGHD3-16,IGHD3-16*02,Europe,671
IGHD3-16,IGHD3-16*02,Asia,713
IGHD3-16,IGHD3-16*02_S3738,Africa,2
IGHD3-16,IGHD3-16*02_S3738,Europe,2
IGHD3-16,IGHD3-16*02_S3738,Asia,0
IGHD3-16,IGHD3-16*03,Africa,58
IGHD3-16,IGHD3-16*03,Europe,93
IGHD3-16,IGHD3-16*03,Asia,63
IGHD3-3,IGHD3-3*01,Africa,617
IGHD3-3,IGHD3-3*01,Europe,762
IGHD3-3,IGHD3-3*01,Asia,624
IGHD3-3,IGHD3-3*03_S4150,Africa,81
IGHD3-3,IGHD3-3*03_S4150,Europe,2
IGHD3-3,IGHD3-3*03_S4150,Asia,0
IGHD4-17,IGHD4-17*01,Africa,339
IGHD4-17,IGHD4-17*01,Europe,336
IGHD4-17,IGHD4-17*01,Asia,352
IGHD4-17,IGHD4-17*01/IGHD4-4*01_S0251,Africa,351
IGHD4-17,IGHD4-17*01/IGHD4-4*01_S0251,Europe,346
IGHD4-17,IGHD4-17*01/IGHD4-4*01_S0251,Asia,362
IGHD4-23,IGHD4-23*01,Africa,438
IGHD4-23,IGHD4-23*01,Europe,434
IGHD4-23,IGHD4-23*01,Asia,492
IGHD4-23,IGHD4-23*01_S7487,Africa,54
IGHD4-23,IGHD4-23*01_S7487,Europe,2
IGHD4-23,IGHD4-23*01_S7487,Asia,0
IGHD4-4,IGHD4-11*01/IGHD4-4*01,Africa,467
IGHD4-4,IGHD4-11*01/IGHD4-4*01,Europe,469
IGHD4-4,IGHD4-11*01/IGHD4-4*01,Asia,543
IGHD4-4,IGHD4-4*01,Africa,301
IGHD4-4,IGHD4-4*01,Europe,295
IGHD4-4,IGHD4-4*01,Asia,247
IGHD4-4,IGHD4-4*01_S0251,Africa,5
IGHD4-4,IGHD4-4*01_S0251,Europe,0
IGHD4-4,IGHD4-4*01_S0251,Asia,0
IGHD4-4,IGHD4-4*01_S6581,Africa,1
IGHD4-4,IGHD4-4*01_S6581,Europe,0
IGHD4-4,IGHD4-4*01_S6581,Asia,0
IGHD5-18,IGHD5-18*01,Africa,533
IGHD5-18,IGHD5-18*01,Europe,547
IGHD5-18,IGHD5-18*01,Asia,509
IGHD5-18,IGHD5-18*02,Africa,43
IGHD5-18,IGHD5-18*02,Europe,57
IGHD5-18,IGHD5-18*02,Asia,47
IGHD5-18/5-5,IGHD5-18*01/IGHD5-5*01,Africa,704
IGHD5-18/5-5,IGHD5-18*01/IGHD5-5*01,Europe,720
IGHD5-18/5-5,IGHD5-18*01/IGHD5-5*01,Asia,677
IGHD5-18/5-5,IGHD5-18*02,Africa,28
IGHD5-18/5-5,IGHD5-18*02,Europe,34
IGHD5-18/5-5,IGHD5-18*02,Asia,27
IGHD5-18/5-5,IGHD5-5*01_S6943,Africa,2
IGHD5-18/5-5,IGHD5-5*01_S6943,Europe,0
IGHD5-18/5-5,IGHD5-5*01_S6943,Asia,0
IGHD5-5,IGHD5-5*01,Africa,584
IGHD5-5,IGHD5-5*01,Europe,662
IGHD5-5,IGHD5-5*01,Asia,504
IGHD5-5,IGHD5-5*01_S6943,Africa,2
IGHD5-5,IGHD5-5*01_S6943,Europe,0
IGHD5-5,IGHD5-5*01_S6943,Asia,0
IGHD6-13,IGHD6-13*01,Africa,769
IGHD6-13,IGHD6-13*01,Europe,765
IGHD6-13,IGHD6-13*01,Asia,789
IGHD6-13,IGHD6-13*01_S0744,Africa,2
IGHD6-13,IGHD6-13*01_S0744,Europe,1
IGHD6-13,IGHD6-13*01_S0744,Asia,1
IGHD6-13,IGHD6-13*01_S5237,Africa,1
IGHD6-13,IGHD6-13*01_S5237,Europe,0
IGHD6-13,IGHD6-13*01_S5237,Asia,2
IGHD6-19,IGHD6-19*01,Africa,757
IGHD6-19,IGHD6-19*01,Europe,752
IGHD6-19,IGHD6-19*01,Asia,771
IGHD6-19,IGHD6-19*01_S2183,Africa,1
IGHD6-19,IGHD6-19*01_S2183,Europe,0
IGHD6-19,IGHD6-19*01_S2183,Asia,1
IGHD6-25,IGHD6-25*01,Africa,546
IGHD6-25,IGHD6-25*01,Europe,529
IGHD6-25,IGHD6-25*01,Asia,586
IGHD6-25,IGHD6-25*01_S2758,Africa,0
IGHD6-25,IGHD6-25*01_S2758,Europe,1
IGHD6-25,IGHD6-25*01_S2758,Asia,0
IGHD6-6,IGHD6-6*01,Africa,753
IGHD6-6,IGHD6-6*01,Europe,786
IGHD6-6,IGHD6-6*01,Asia,640
IGHD6-6,IGHD6-6*01_S1911,Africa,1
IGHD6-6,IGHD6-6*01_S1911,Europe,0
IGHD6-6,IGHD6-6*01_S1911,Asia,0
"

# ---- 1) parse the raw table - this is the ONLY hardcoded input; every
#         number below is computed from it, nothing is pre-summarized ------
raw_long <- read.csv(text = raw_csv, stringsAsFactors = FALSE)
gene_order <- unique(raw_long$gene)
raw_long$gene   <- factor(raw_long$gene, levels = gene_order)
raw_long$region <- factor(raw_long$region, levels = c("Africa", "Europe", "Asia"))

# ---- 2) percentages + 95% Wilson confidence intervals (base R, no extra pkg) ----
totals <- aggregate(count ~ gene + region, data = raw_long, sum)
names(totals)[3] <- "total"
raw_long <- merge(raw_long, totals, by = c("gene", "region"))
raw_long$pct <- 100 * raw_long$count / raw_long$total

ci <- t(mapply(function(x, n) {
  if (n == 0) return(c(NA_real_, NA_real_))
  100 * prop.test(x, n, correct = TRUE)$conf.int
}, raw_long$count, raw_long$total))
raw_long$ci_low  <- ci[, 1]
raw_long$ci_high <- ci[, 2]

# ---- 3) chi-square test of independence, per gene, built straight from
#         raw_long via xtabs (allele x region contingency table) - then ONE
#         Benjamini-Hochberg correction across all 16 tested genes (the
#         correct family of tests, not just the ones that end up significant)
chi_results <- do.call(rbind, lapply(split(raw_long, raw_long$gene), function(sub) {
  tab <- xtabs(count ~ allele + region, data = sub)
  test <- suppressWarnings(chisq.test(tab))
  data.frame(gene = as.character(unique(sub$gene)), n_allele = nrow(tab),
             chi2 = unname(test$statistic), df = unname(test$parameter),
             p_raw = test$p.value, stringsAsFactors = FALSE)
}))
chi_results$p_BH <- p.adjust(chi_results$p_raw, method = "BH")
chi_results <- chi_results[order(chi_results$p_BH), ]
cat("\n==== Chi-square test, ALL 16 tested genes (BH correction over this full set) ====\n")
print(chi_results, row.names = FALSE)

sig_genes <- as.character(chi_results$gene[chi_results$p_BH < 0.05])
cat(sprintf("\n%d / %d genes statistically significant (BH p<0.05): %s\n",
            length(sig_genes), nrow(chi_results), paste(sig_genes, collapse = ", ")))

# for PLOTTING ORDER only (the stats table above stays sorted by
# significance): group genes by IGHD family number (the digit right after
# "IGHD", e.g. IGHD4-4 and IGHD4-23 are both family 4) so family-mates are
# plotted next to each other, then by position number within a family
family_num   <- as.integer(sub("^IGHD(\\d+)-.*$", "\\1", sig_genes))
position_num <- as.integer(sub("^IGHD\\d+-(\\d+).*$", "\\1", sig_genes))
sig_genes_plot_order <- sig_genes[order(family_num, position_num)]

# ---- 4) GraphPad Prism-style theme ------------------------------------------
theme_prism <- theme_bw(base_size = 12) +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black", size = 0.6),  # `size` (not linewidth) for compat with older ggplot2
    axis.ticks       = element_line(color = "black"),
    axis.text        = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    legend.position  = "top",
    legend.title     = element_blank(),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 9.5, color = "grey30"),
    text             = element_text(family = "sans")
  )

region_colors <- c(Africa = "#2A78D6", Europe = "#EB6834", Asia = "#1BAF7A")

# ---- 4b) where to save the HD figures. Saves into your current R working
#          directory (see getwd() in the console, or Session > Set Working
#          Directory > Choose Directory in RStudio to pick a folder first).
out_dir <- "figures"
dir.create(out_dir, showWarnings = FALSE)
cat(sprintf("\nSaving figures into: %s\n", normalizePath(out_dir)))

# ---- 5) one figure per SIGNIFICANT gene, one facet panel per allele (never
#         combined into a single set of bars - each allele keeps its own
#         separate panel). Printed in family order (see sig_genes_plot_order
#         above) so family-mates (e.g. IGHD4-4, IGHD4-23) appear back to back.
for (g in sig_genes_plot_order) {
  sub  <- subset(raw_long, gene == g)
  stat <- subset(chi_results, gene == g)
  subtitle <- sprintf("Chi-square test (recomputed from this data in R): chi2=%.2f, df=%d\nraw p=%.2e, BH-adjusted p=%.2e (SIGNIFICANT *)",
                       stat$chi2, stat$df, stat$p_raw, stat$p_BH)

  p <- ggplot(sub, aes(x = region, y = pct, fill = region)) +
    geom_col(width = 0.65, color = "black", size = 0.3) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, size = 0.4, color = "black") +
    facet_wrap(~ allele) +
    scale_fill_manual(values = region_colors) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    labs(title = as.character(g), subtitle = subtitle, x = NULL, y = "Allele frequency (%)") +
    theme_prism

  print(p)

  n_allele_panels <- nlevels(factor(sub$allele))
  fig_width <- max(6, 3 * n_allele_panels)  # widen for genes with more allele panels

  # vector PDF - what most journals actually want, scales to any size with no quality loss
  ggsave(file.path(out_dir, paste0(gsub("[/*]", "-", g), ".pdf")), plot = p,
         width = fig_width, height = 5, units = "in")
  # high-resolution TIFF as well, in case a raster format is required instead
  ggsave(file.path(out_dir, paste0(gsub("[/*]", "-", g), ".tiff")), plot = p,
         width = fig_width, height = 5, units = "in", dpi = 600, compression = "lzw")
}

cat(sprintf("\nDone - %d PDF + %d TIFF file(s) written to: %s\n",
            length(sig_genes_plot_order), length(sig_genes_plot_order), normalizePath(out_dir)))
