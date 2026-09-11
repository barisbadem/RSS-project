# =============================================================================
# D-REGION allele frequency by geography - publication-style bar charts
#
# DATA SOURCE (read this before using these numbers anywhere): these raw
# haplotype counts are NOT copied from a published journal article. They are
# our own analysis output, computed directly from KIARVA's real, public
# 1000 Genomes IGHD genotype file (ScilifelabDataCentre/kiarva-backend,
# CC BY-NC 4.0), using real ~410-person-per-superpopulation cohorts
# (Africa=AFR, Europe=EUR, Asia=EAS+SAS), counted per real observed
# haplotype. The same underlying counts were already cross-checked in two
# other formats in this project: an Excel workbook with the full
# observed/expected/chi-square-contribution tables, and a GraphPad Prism
# (.pzfx) contingency table per gene. All three must show identical numbers.
#
# ALL 16 real D genes that have 2+ observed alleles are included below (not
# just the ones that turn out significant) because the Benjamini-Hochberg
# correction must be computed across the SAME family of tests as the
# original analysis - correcting across only a pre-selected subset of
# "already significant" genes would silently understate the adjusted
# p-values. This script derives chi-square + BH from these raw counts in R
# itself, then plots only the genes that come out significant - nothing
# about "which genes are significant" is hardcoded.
# =============================================================================

library(ggplot2)

# ---- 1) raw observed haplotype counts, ALL 16 multi-allele D genes ---------
GENE_VEC <- c("IGHD2-2", "IGHD2-2", "IGHD2-2", "IGHD2-2", "IGHD2-2", "IGHD2-21", "IGHD2-21",
  "IGHD2-8", "IGHD2-8", "IGHD3-10", "IGHD3-10", "IGHD3-10", "IGHD3-16", "IGHD3-16", "IGHD3-16",
  "IGHD3-3", "IGHD3-3", "IGHD4-17", "IGHD4-17", "IGHD4-23", "IGHD4-23", "IGHD4-4", "IGHD4-4",
  "IGHD4-4", "IGHD4-4", "IGHD5-18", "IGHD5-18", "IGHD5-18/5-5", "IGHD5-18/5-5", "IGHD5-18/5-5",
  "IGHD5-5", "IGHD5-5", "IGHD6-13", "IGHD6-13", "IGHD6-13", "IGHD6-19", "IGHD6-19", "IGHD6-25",
  "IGHD6-25", "IGHD6-6", "IGHD6-6")

ALLELE_VEC <- c("IGHD2-2*01", "IGHD2-2*01_S7032", "IGHD2-2*02", "IGHD2-2*04_S0329",
  "IGHD2-2*05_S3769", "IGHD2-21*01", "IGHD2-21*02", "IGHD2-8*01", "IGHD2-8*02", "IGHD3-10*01",
  "IGHD3-10*01_S3902", "IGHD3-10*03", "IGHD3-16*02", "IGHD3-16*02_S3738", "IGHD3-16*03",
  "IGHD3-3*01", "IGHD3-3*03_S4150", "IGHD4-17*01", "IGHD4-17*01/IGHD4-4*01_S0251", "IGHD4-23*01",
  "IGHD4-23*01_S7487", "IGHD4-11*01/IGHD4-4*01", "IGHD4-4*01", "IGHD4-4*01_S0251",
  "IGHD4-4*01_S6581", "IGHD5-18*01", "IGHD5-18*02", "IGHD5-18*01/IGHD5-5*01", "IGHD5-18*02",
  "IGHD5-5*01_S6943", "IGHD5-5*01", "IGHD5-5*01_S6943", "IGHD6-13*01", "IGHD6-13*01_S0744",
  "IGHD6-13*01_S5237", "IGHD6-19*01", "IGHD6-19*01_S2183", "IGHD6-25*01", "IGHD6-25*01_S2758",
  "IGHD6-6*01", "IGHD6-6*01_S1911")

AFRICA_VEC <- c(454, 4, 309, 1, 28, 279, 299, 531, 79, 653, 0, 51, 688, 2, 58, 617, 81, 339, 351,
  438, 54, 467, 301, 5, 1, 533, 43, 704, 28, 2, 584, 2, 769, 2, 1, 757, 1, 546, 0, 753, 1)

EUROPE_VEC <- c(327, 0, 409, 11, 1, 92, 318, 482, 80, 599, 4, 81, 671, 2, 93, 762, 2, 336, 346,
  434, 2, 469, 295, 0, 0, 547, 57, 720, 34, 0, 662, 0, 765, 1, 0, 752, 0, 529, 1, 786, 0)

ASIA_VEC <- c(603, 0, 200, 5, 0, 140, 444, 448, 82, 682, 0, 54, 713, 0, 63, 624, 0, 352, 362,
  492, 0, 543, 247, 0, 0, 509, 47, 677, 27, 0, 504, 0, 789, 1, 2, 771, 1, 586, 0, 640, 0)

wide <- data.frame(gene = GENE_VEC, allele = ALLELE_VEC,
                    Africa = AFRICA_VEC, Europe = EUROPE_VEC, Asia = ASIA_VEC,
                    stringsAsFactors = FALSE)
gene_order <- unique(wide$gene)

# long format: one row per gene x allele x region
raw <- do.call(rbind, lapply(c("Africa", "Europe", "Asia"), function(reg) {
  data.frame(gene = wide$gene, allele = wide$allele, region = reg, count = wide[[reg]],
             stringsAsFactors = FALSE)
}))
raw$gene   <- factor(raw$gene, levels = gene_order)
raw$region <- factor(raw$region, levels = c("Africa", "Europe", "Asia"))
region_tr  <- c(Africa = "Afrika", Europe = "Avrupa", Asia = "Asya")
raw$region_label <- factor(region_tr[as.character(raw$region)], levels = region_tr)

# ---- 2) percentages + 95% Wilson confidence intervals (base R, no extra pkg) ----
totals <- aggregate(count ~ gene + region, data = raw, sum)
names(totals)[3] <- "total"
raw <- merge(raw, totals, by = c("gene", "region"))
raw$pct <- 100 * raw$count / raw$total

ci <- t(mapply(function(x, n) {
  if (n == 0) return(c(NA_real_, NA_real_))
  100 * prop.test(x, n, correct = TRUE)$conf.int
}, raw$count, raw$total))
raw$ci_low  <- ci[, 1]
raw$ci_high <- ci[, 2]

# ---- 3) chi-square test of independence, per gene, across ALL 16 tested
#         genes, then ONE Benjamini-Hochberg correction across all of them
#         (the correct family of tests - not just the ones that end up
#         significant) ---------------------------------------------------
chi_results <- do.call(rbind, lapply(split(wide, wide$gene), function(sub) {
  tab <- as.matrix(sub[, c("Africa", "Europe", "Asia")])
  rownames(tab) <- sub$allele
  test <- suppressWarnings(chisq.test(tab))
  data.frame(gene = unique(sub$gene), n_allele = nrow(sub),
             chi2 = unname(test$statistic), df = unname(test$parameter),
             p_raw = test$p.value, stringsAsFactors = FALSE)
}))
chi_results$p_BH <- p.adjust(chi_results$p_raw, method = "BH")
chi_results <- chi_results[order(chi_results$p_BH), ]
cat("\n==== Ki-kare testi, TUM 16 test edilen gen (BH duzeltmesi bu tam kumeye gore) ====\n")
print(chi_results, row.names = FALSE)

sig_genes <- chi_results$gene[chi_results$p_BH < 0.05]
cat(sprintf("\n%d / %d gen istatistiksel olarak anlamli (BH p<0.05): %s\n",
            length(sig_genes), nrow(chi_results), paste(sig_genes, collapse = ", ")))

# ---- 4) GraphPad Prism-style theme ------------------------------------------
theme_prism <- theme_bw(base_size = 12) +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black", linewidth = 0.6),
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

region_colors <- c(Afrika = "#2A78D6", Avrupa = "#EB6834", Asya = "#1BAF7A")

# ---- 5) one figure per SIGNIFICANT gene, one facet panel per allele (never
#         combined into a single set of bars - each allele keeps its own
#         separate panel) ----------------------------------------------------
for (g in sig_genes) {
  sub  <- subset(raw, gene == g)
  stat <- subset(chi_results, gene == g)
  subtitle <- sprintf("Ki-kare testi (bu veriden R'de yeniden hesaplandi): chi2=%.2f, df=%d\nham p=%.2e, BH-duzeltmeli p=%.2e (ANLAMLI *)",
                       stat$chi2, stat$df, stat$p_raw, stat$p_BH)

  p <- ggplot(sub, aes(x = region_label, y = pct, fill = region_label)) +
    geom_col(width = 0.65, color = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, linewidth = 0.4, color = "black") +
    facet_wrap(~ allele) +
    scale_fill_manual(values = region_colors) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    labs(title = as.character(g), subtitle = subtitle, x = NULL, y = "Alel sikligi (%)") +
    theme_prism

  print(p)
}
