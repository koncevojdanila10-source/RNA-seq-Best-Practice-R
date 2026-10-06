#!/usr/bin/env Rscript
#
# x4_ase_crosscheck.R — is MBASED's allelic imbalance real, or is it the
# genotype calls?
#
# Step 16 called 54% of gene-sample tests imbalanced. Published estimates for
# genuine ASE sit in the low single digits, so something is wrong, and the
# question is what. Two hypotheses, and they make opposite predictions:
#
#   biology  — the imbalance is real. Then it should be a property of genes,
#              reproducible across samples, and unrelated to how each library
#              was sequenced.
#
#   genotypes — the heterozygous calls are wrong. The upstream pipeline found
#              that only 29% of its 30 059 "heterozygous" sites carried an
#              alternate fraction above 0.4. A homozygous site miscalled as
#              heterozygous looks perfectly monoallelic to MBASED, which does
#              not question the call. Then the per-sample rate of significant
#              ASE should track the per-sample miscall rate, and the significant
#              genes should be the ones sitting at extreme allele fractions.
#
# Three tests are run, each able to separate the two:
#
#   1. per-sample ASE rate against the miscall rate and against library
#      composition;
#   2. allele fractions of the genes MBASED calls significant, against the rest;
#   3. agreement with the upstream per-site binomial test, which used the same
#      reads but a different model.
#
suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(GenomicRanges)
})
source(file.path(Sys.getenv("R_DIR", unset = "R"), "00_common.R"))

banner("x4 - is the ASE signal biology or genotyping?")

load(file.path(RDATA_DIR, "16_ase_mbased.RData"))

st <- load_sample_table()
rd <- read_tsv(file.path(PROJECT_DIR, "inputs", "read_distribution.tsv"),
               show_col_types = FALSE) |> as.data.frame()
rownames(rd) <- rd$sample

genes_gr <- readRDS(file.path(RESULTS_DIR, "chr19_txdb_genes.rds"))

MIN_DEPTH <- 10
# The upstream threshold: a real heterozygote should carry well over this share
# of alternate reads. Below it, the call itself is in doubt.
HET_MIN_ALT <- 0.4

# ---------------------------------------------------------------------------
# 1. Per-sample miscall rate, straight from the ASEReadCounter tables
# ---------------------------------------------------------------------------
banner("how many heterozygous calls look like miscalls?")

read_sites <- function(sid) {
    f <- file.path(ASE_DIR, paste0(sid, ".ase.tsv"))
    if (!file.exists(f)) return(NULL)
    x <- read_tsv(f, show_col_types = FALSE) |>
        filter(totalCount >= MIN_DEPTH) |>
        mutate(sample = sid, alt_fraction = altCount / totalCount)
    x
}

sites <- bind_rows(lapply(st$sample, read_sites))
message("sites at depth >= ", MIN_DEPTH, ": ", nrow(sites))

miscall <- sites |>
    group_by(sample) |>
    summarise(
        sites            = n(),
        median_alt       = round(median(alt_fraction), 4),
        pct_below_0.4    = round(100 * mean(alt_fraction < HET_MIN_ALT), 1),
        pct_genuine_het  = round(100 * mean(alt_fraction >= HET_MIN_ALT), 1),
        .groups = "drop"
    )

per_sample <- summary_tab |>
    select(sample, group, genes, significant, pct_significant, median_MAF) |>
    left_join(miscall, by = "sample") |>
    mutate(intronic_pct = rd[sample, "intronic_pct"])
print(as.data.frame(per_sample), row.names = FALSE)
save_tab(per_sample, "x4_per_sample")

# ---------------------------------------------------------------------------
# 2. What predicts the per-sample ASE rate?
# ---------------------------------------------------------------------------
banner("what the per-sample ASE rate tracks")

for (v in c("pct_below_0.4", "intronic_pct", "sites")) {
    s <- spearman_test(per_sample$pct_significant, per_sample[[v]])
    message(sprintf("  %%significant vs %-14s %s", v, s$text))
}
s_gene <- spearman_test(per_sample$median_MAF, per_sample$pct_below_0.4)
message("  median MAF     vs miscall rate:  ", s_gene$text)

save_fig(
    ggplot(per_sample, aes(pct_below_0.4, pct_significant,
                           colour = group, label = sample)) +
        geom_smooth(aes(group = 1), method = "lm", se = FALSE,
                    colour = "grey40", linewidth = 0.5, formula = y ~ x) +
        geom_point(size = 3) +
        ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
        scale_colour_manual(values = GROUP_COLOURS) +
        labs(x = paste0("sites with alternate fraction < ", HET_MIN_ALT, " (%)"),
             y = "genes called imbalanced by MBASED (%)",
             title = "Apparent ASE tracks the rate of doubtful heterozygous calls") +
        theme_bw(),
    "x4_ase_vs_miscall", width = 8, height = 6
)

# ---------------------------------------------------------------------------
# 3. Are the significant genes the ones at extreme allele fractions?
# ---------------------------------------------------------------------------
# If MBASED were finding regulation, its significant genes would sit at
# moderate imbalance. If it is finding miscalled homozygotes, they sit at the
# extremes, where one allele is barely observed at all.
banner("allele fractions of significant against non-significant genes")

snv <- GRanges(seqnames = sites$contig, ranges = IRanges(sites$position, width = 1))
hits <- findOverlaps(snv, genes_gr, ignore.strand = TRUE)
sites$geneID <- NA_character_
sites$geneID[queryHits(hits)] <- mcols(genes_gr)$gene_id[subjectHits(hits)]

gene_sites <- sites |>
    filter(!is.na(geneID)) |>
    group_by(sample, geneID) |>
    summarise(n_sites = n(),
              mean_alt = mean(alt_fraction),
              min_alt  = min(alt_fraction),
              .groups = "drop")

joined <- ase_all |>
    select(sample, geneID, MAF, imbalance, pASE_BH, sig) |>
    inner_join(gene_sites, by = c("sample", "geneID"))
message("genes matched between MBASED and the site tables: ", nrow(joined))

cmp <- joined |>
    group_by(sig) |>
    summarise(genes = n(),
              median_mean_alt = round(median(mean_alt), 4),
              median_min_alt  = round(median(min_alt), 4),
              pct_with_a_site_below_0.1 = round(100 * mean(min_alt < 0.1), 1),
              .groups = "drop")
print(as.data.frame(cmp), row.names = FALSE)
save_tab(cmp, "x4_significant_vs_allele_fraction")

w <- wilcox.test(joined$min_alt[joined$sig], joined$min_alt[!joined$sig])
message("minimum alternate fraction, significant vs not: Wilcoxon p = ",
        format.pval(w$p.value, digits = 3, eps = 1e-300))

save_fig(
    ggplot(joined, aes(sig, min_alt, fill = sig)) +
        geom_violin(alpha = 0.7) +
        geom_boxplot(width = 0.15, outlier.size = 0.4) +
        scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick")) +
        scale_x_discrete(labels = c(`FALSE` = "not called", `TRUE` = "called imbalanced")) +
        labs(x = NULL, y = "lowest alternate fraction in the gene",
             title = "Genes called imbalanced are the ones with a near-absent allele") +
        theme_bw() + theme(legend.position = "none"),
    "x4_min_alt_by_significance", width = 7, height = 6
)

# ---------------------------------------------------------------------------
# 4. Against the upstream binomial test
# ---------------------------------------------------------------------------
# A second model on the same reads. It shares MBASED's blind spot — it also
# trusts the genotype calls — so agreement is not validation. Disagreement
# would be informative, and the size of each list is worth recording.
banner("MBASED against the upstream per-site binomial test")

bin_path <- file.path(ASE_DIR, "ase_tested.tsv")
if (!file.exists(bin_path)) {
    message("upstream ase_tested.tsv not found - skipping this comparison")
} else {
    bin <- read_tsv(bin_path, show_col_types = FALSE)
    bsnv <- GRanges(seqnames = bin$contig, ranges = IRanges(bin$position, width = 1))
    bh <- findOverlaps(bsnv, genes_gr, ignore.strand = TRUE)
    bin$geneID <- NA_character_
    bin$geneID[queryHits(bh)] <- mcols(genes_gr)$gene_id[subjectHits(bh)]

    # A gene counts as imbalanced by the binomial test if any of its sites is.
    bin_gene <- bin |>
        filter(!is.na(geneID)) |>
        group_by(sample, geneID) |>
        summarise(bin_sig = any(fdr < FDR_CUT, na.rm = TRUE), .groups = "drop")

    both <- ase_all |>
        select(sample, geneID, mbased_sig = sig) |>
        inner_join(bin_gene, by = c("sample", "geneID"))

    tb <- table(MBASED = both$mbased_sig, binomial = both$bin_sig)
    print(tb)
    jac <- sum(both$mbased_sig & both$bin_sig) /
           sum(both$mbased_sig | both$bin_sig)
    message("gene-sample pairs compared: ", nrow(both))
    message("MBASED significant:   ", sum(both$mbased_sig),
            "  binomial significant: ", sum(both$bin_sig))
    message("Jaccard index: ", round(jac, 3))
    message("Both methods take the heterozygous calls on trust, so this ",
            "agreement measures how similarly they respond to the same bad ",
            "input - it is not independent confirmation.")
    save_tab(as.data.frame(tb), "x4_mbased_vs_binomial")
}

# ---------------------------------------------------------------------------
# 5. Verdict
# ---------------------------------------------------------------------------
banner("verdict")

s_mis <- spearman_test(per_sample$pct_significant, per_sample$pct_below_0.4)
message("Per-sample ASE rate against miscall rate: ", s_mis$text)
message("Significant genes carry a lowest alternate fraction of ",
        cmp$median_min_alt[cmp$sig], " against ",
        cmp$median_min_alt[!cmp$sig], " for the rest.")
if (s_mis$rho > 0.6 && cmp$median_min_alt[cmp$sig] < cmp$median_min_alt[!cmp$sig]) {
    message("")
    message("Both predictions of the genotyping hypothesis hold. The imbalance ",
            "MBASED reports is a property of the variant calls, not of allelic")
    message("regulation, and it cannot be filtered away: the genotypes were ",
            "derived from the same reads the expression is measured in, so any")
    message("filter that removes the doubtful calls also removes the evidence. ",
            "ASE is NOT assessable on this dataset - the same verdict the")
    message("upstream pipeline reached, now reproduced by a second method with ",
            "a different statistical model.")
    message("")
    message("Assessing it properly needs genotypes from an orthogonal source: ",
            "DNA from the same individuals, or at minimum an external panel.")
} else {
    message("The genotyping hypothesis is not supported in this form - ",
            "the results need a closer look before any verdict.")
}

save(per_sample, joined, cmp, file = file.path(RDATA_DIR, "x4_ase_crosscheck.RData"))
banner("x4 done")
