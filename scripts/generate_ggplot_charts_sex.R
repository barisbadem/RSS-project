# =============================================================================
# D-REGION allele frequency by BIOLOGICAL SEX - publication-style bar charts
#
# Companion to the geography figures, run with a deliberately identical
# method so the two are directly comparable: per gene, an allele x group
# contingency table, a chi-square test of independence, and ONE
# Benjamini-Hochberg correction across the full family of tested genes.
#
# DATA SOURCE: raw observed haplotype counts from our own analysis of
# KIARVA's real, public 1000 Genomes IGHD genotype data (CC BY-NC 4.0),
# joined to recorded sex from the 1000 Genomes pedigree file. NOT copied
# from any published article. Same counts as the accompanying Excel workbook.
#
# Unlike the geography analysis, nothing is subsampled: the sexes are already
# near-balanced (1203 male / 1227 female), so every individual with a
# recorded sex is used.
#
# Sex is also independent of superpopulation in this cohort
# (chi2 = 4.76, df = 4, p = 0.313), so real geographic allele differences
# cannot leak in disguised as sex differences.
#
# Everything below the one raw table is COMPUTED, not restated: R itself
# builds the contingency tables, runs the tests, applies the correction,
# derives percentages and 95% Wilson confidence intervals, and only then
# plots. Nothing is pre-summarized or hardcoded per gene/allele.
#
# IGH sits on chromosome 14 - an autosome - so the expectation going in is
# NO difference. Every tested gene is plotted, because here the informative
# result is the absence of an effect, and that has to be shown, not asserted.
# =============================================================================

library(ggplot2)

raw_csv <- "gene,allele,sex,count
IGHD2-2,IGHD2-2*01,Male,1444
IGHD2-2,IGHD2-2*01,Female,1473
IGHD2-2,IGHD2-2*01_S7032,Male,1
IGHD2-2,IGHD2-2*01_S7032,Female,4
IGHD2-2,IGHD2-2*02,Male,835
IGHD2-2,IGHD2-2*02,Female,842
IGHD2-2,IGHD2-2*04_S0329,Male,6
IGHD2-2,IGHD2-2*04_S0329,Female,16
IGHD2-2,IGHD2-2*05_S3769,Male,26
IGHD2-2,IGHD2-2*05_S3769,Female,33
IGHD2-8,IGHD2-8*01,Male,1460
IGHD2-8,IGHD2-8*01,Female,1450
IGHD2-8,IGHD2-8*01_S7032,Male,3
IGHD2-8,IGHD2-8*01_S7032,Female,0
IGHD2-8,IGHD2-8*02,Male,225
IGHD2-8,IGHD2-8*02,Female,224
IGHD2-21,IGHD2-21*01,Male,535
IGHD2-21,IGHD2-21*01,Female,512
IGHD2-21,IGHD2-21*02,Male,1115
IGHD2-21,IGHD2-21*02,Female,1188
IGHD3-3,IGHD3-3*01,Male,1957
IGHD3-3,IGHD3-3*01,Female,1982
IGHD3-3,IGHD3-3*03_S4150,Male,61
IGHD3-3,IGHD3-3*03_S4150,Female,60
IGHD3-10,IGHD3-10*01,Male,1960
IGHD3-10,IGHD3-10*01,Female,1983
IGHD3-10,IGHD3-10*01_S3902,Male,1
IGHD3-10,IGHD3-10*01_S3902,Female,3
IGHD3-10,IGHD3-10*03,Male,173
IGHD3-10,IGHD3-10*03,Female,174
IGHD3-16,IGHD3-16*02,Male,2072
IGHD3-16,IGHD3-16*02,Female,2106
IGHD3-16,IGHD3-16*02_S0652,Male,1
IGHD3-16,IGHD3-16*02_S0652,Female,1
IGHD3-16,IGHD3-16*02_S3738,Male,2
IGHD3-16,IGHD3-16*02_S3738,Female,4
IGHD3-16,IGHD3-16*03,Male,189
IGHD3-16,IGHD3-16*03,Female,193
IGHD3-22,IGHD3-22*01,Male,2056
IGHD3-22,IGHD3-22*01,Female,2085
IGHD3-22,IGHD3-22*01_S7002,Male,0
IGHD3-22,IGHD3-22*01_S7002,Female,3
IGHD4-4,IGHD4-11*01/IGHD4-4*01,Male,1450
IGHD4-4,IGHD4-11*01/IGHD4-4*01,Female,1481
IGHD4-4,IGHD4-4*01,Male,844
IGHD4-4,IGHD4-4*01,Female,831
IGHD4-4,IGHD4-4*01_S0251,Male,3
IGHD4-4,IGHD4-4*01_S0251,Female,3
IGHD4-4,IGHD4-4*01_S6581,Male,1
IGHD4-4,IGHD4-4*01_S6581,Female,1
IGHD4-17,IGHD4-17*01,Male,1011
IGHD4-17,IGHD4-17*01,Female,1034
IGHD4-17,IGHD4-17*01/IGHD4-4*01_S0251,Male,1043
IGHD4-17,IGHD4-17*01/IGHD4-4*01_S0251,Female,1062
IGHD4-23,IGHD4-23*01,Male,1363
IGHD4-23,IGHD4-23*01,Female,1431
IGHD4-23,IGHD4-23*01_S7487,Male,51
IGHD4-23,IGHD4-23*01_S7487,Female,57
IGHD5-5,IGHD5-5*01,Male,1724
IGHD5-5,IGHD5-5*01,Female,1700
IGHD5-5,IGHD5-5*01_S6943,Male,0
IGHD5-5,IGHD5-5*01_S6943,Female,2
IGHD5-18,IGHD5-18*01,Male,1601
IGHD5-18,IGHD5-18*01,Female,1621
IGHD5-18,IGHD5-18*02,Male,119
IGHD5-18,IGHD5-18*02,Female,153
IGHD5-18/5-5,IGHD5-18*01/IGHD5-5*01,Male,2085
IGHD5-18/5-5,IGHD5-18*01/IGHD5-5*01,Female,2079
IGHD5-18/5-5,IGHD5-18*02,Male,75
IGHD5-18/5-5,IGHD5-18*02,Female,91
IGHD5-18/5-5,IGHD5-5*01_S6943,Male,0
IGHD5-18/5-5,IGHD5-5*01_S6943,Female,2
IGHD6-6,IGHD6-6*01,Male,2101
IGHD6-6,IGHD6-6*01,Female,2139
IGHD6-6,IGHD6-6*01_S1911,Male,1
IGHD6-6,IGHD6-6*01_S1911,Female,1
IGHD6-6,IGHD6-6*01_S4962,Male,2
IGHD6-6,IGHD6-6*01_S4962,Female,2
IGHD6-13,IGHD6-13*01,Male,2283
IGHD6-13,IGHD6-13*01,Female,2326
IGHD6-13,IGHD6-13*01_S0744,Male,1
IGHD6-13,IGHD6-13*01_S0744,Female,4
IGHD6-13,IGHD6-13*01_S2484,Male,2
IGHD6-13,IGHD6-13*01_S2484,Female,2
IGHD6-13,IGHD6-13*01_S5237,Male,7
IGHD6-13,IGHD6-13*01_S5237,Female,2
IGHD6-13,IGHD6-13*01_S5932,Male,3
IGHD6-13,IGHD6-13*01_S5932,Female,0
IGHD6-19,IGHD6-19*01,Male,2245
IGHD6-19,IGHD6-19*01,Female,2289
IGHD6-19,IGHD6-19*01_S2183,Male,3
IGHD6-19,IGHD6-19*01_S2183,Female,1
IGHD6-25,IGHD6-25*01,Male,1612
IGHD6-25,IGHD6-25*01,Female,1732
IGHD6-25,IGHD6-25*01_S2758,Male,0
IGHD6-25,IGHD6-25*01_S2758,Female,4
"

# ---- 1) parse the raw table - the ONLY hardcoded input ---------------------
raw_long <- read.csv(text = raw_csv, stringsAsFactors = FALSE)
gene_order <- unique(raw_long$gene)
raw_long$gene <- factor(raw_long$gene, levels = gene_order)
raw_long$sex  <- factor(raw_long$sex,  levels = c("Male", "Female"))

# ---- 2) percentages + 95% Wilson confidence intervals (base R) -------------
totals <- aggregate(count ~ gene + sex, data = raw_long, sum)
names(totals)[3] <- "total"
raw_long <- merge(raw_long, totals, by = c("gene", "sex"))
raw_long$pct <- 100 * raw_long$count / raw_long$total

ci <- t(mapply(function(x, n) {
  if (n == 0) return(c(NA_real_, NA_real_))
  100 * prop.test(x, n, correct = TRUE)$conf.int
}, raw_long$count, raw_long$total))
raw_long$ci_low  <- ci[, 1]
raw_long$ci_high <- ci[, 2]

# ---- 3) chi-square per gene from xtabs, then ONE BH correction across the
#         full set of tested genes (the correct family, not a selected subset)
chi_results <- do.call(rbind, lapply(split(raw_long, raw_long$gene), function(sub) {
  tab  <- xtabs(count ~ allele + sex, data = sub)
  test <- suppressWarnings(chisq.test(tab))
  data.frame(gene = as.character(unique(sub$gene)), n_allele = nrow(tab),
             chi2 = unname(test$statistic), df = unname(test$parameter),
             p_raw = test$p.value, stringsAsFactors = FALSE)
}))
chi_results$p_BH <- p.adjust(chi_results$p_raw, method = "BH")
chi_results <- chi_results[order(chi_results$p_BH), ]
cat("\n==== Chi-square test by sex, ALL tested genes (BH over this full set) ====\n")
print(chi_results, row.names = FALSE)

sig_genes <- as.character(chi_results$gene[chi_results$p_BH < 0.05])
cat(sprintf("\n%d / %d genes statistically significant (BH p<0.05)%s\n",
            length(sig_genes), nrow(chi_results),
            if (length(sig_genes)) paste0(": ", paste(sig_genes, collapse = ", ")) else " - NONE"))
cat(sprintf("Smallest uncorrected p across all genes: %.4f\n", min(chi_results$p_raw)))

# plot every tested gene, in IGHD family order (family number, then position)
all_genes    <- as.character(chi_results$gene)
family_num   <- as.integer(sub("^IGHD(\\d+)-.*$", "\\1", all_genes))
position_num <- as.integer(sub("^IGHD\\d+-(\\d+).*$", "\\1", all_genes))
genes_plot_order <- all_genes[order(family_num, position_num)]

# ---- 4) GraphPad Prism-style theme -----------------------------------------
theme_prism <- theme_bw(base_size = 12) +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black", size = 0.6),  # `size` not `linewidth`, for older ggplot2
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

sex_colors <- c(Male = "#3A7BD5", Female = "#E8698C")   # blue / pink

# ---- 4b) where the HD figures are saved: a "figures_sex" folder inside your
#          current R working directory (getwd(), or Session > Set Working
#          Directory in RStudio to choose one first)
out_dir <- "figures_sex"
dir.create(out_dir, showWarnings = FALSE)
# (no normalizePath() - on Windows it can choke on non-ASCII characters in a
# path, e.g. a OneDrive folder with Turkish letters)
cat(sprintf("\nSaving figures into: %s/%s\n", getwd(), out_dir))

# ---- 5) one figure per gene, one facet panel per allele (never all alleles
#         collapsed into a single set of bars), in family order -------------
for (g in genes_plot_order) {
  sub  <- subset(raw_long, gene == g)
  stat <- subset(chi_results, gene == g)
  verdict <- if (stat$p_BH < 0.05) "SIGNIFICANT *" else "not significant"
  subtitle <- sprintf("Chi-square test (recomputed from this data in R): chi2=%.2f, df=%d\nraw p=%.3f, BH-adjusted p=%.3f (%s)",
                      stat$chi2, stat$df, stat$p_raw, stat$p_BH, verdict)

  p <- ggplot(sub, aes(x = sex, y = pct, fill = sex)) +
    geom_col(width = 0.6, color = "black", size = 0.3) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.18, size = 0.4, color = "black") +
    facet_wrap(~ allele, nrow = 1) +   # all alleles of a gene side by side, one row
    scale_fill_manual(values = sex_colors) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    labs(title = as.character(g), subtitle = subtitle, x = NULL, y = "Allele frequency (%)") +
    theme_prism

  print(p)

  n_allele_panels <- nlevels(factor(sub$allele))
  fig_width <- max(5.5, 2.6 * n_allele_panels)   # one row, so width grows with panel count

  # vector PDF - what most journals want, scales with no quality loss
  ggsave(file.path(out_dir, paste0(gsub("[/*]", "-", g), "_sex.pdf")), plot = p,
         width = fig_width, height = 5, units = "in")
  # 600 dpi TIFF as well, in case a raster format is required
  ggsave(file.path(out_dir, paste0(gsub("[/*]", "-", g), "_sex.tiff")), plot = p,
         width = fig_width, height = 5, units = "in", dpi = 600, compression = "lzw")
}

cat(sprintf("\nDone - %d PDF + %d TIFF file(s) written to: %s/%s\n",
            length(genes_plot_order), length(genes_plot_order), getwd(), out_dir))
