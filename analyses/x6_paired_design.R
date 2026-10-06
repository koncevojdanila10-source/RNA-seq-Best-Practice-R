#!/usr/bin/env Rscript
#
# x6_paired_design.R — the within-patient contrast, treated as a paired design.
#
# Step 11 fitted ~ 0 + CONDITION on all twelve samples. That model does not
# know that acute (D2_x) and subacute (D1_x) are the SAME four patients sampled
# twice. A model that does know can subtract each person's own baseline, so
# between-person variability stops counting against the time effect:
#
#     ~ PATIENT + CONDITION
#
# This is the one place in the project where the design itself can buy back
# power, so it is worth measuring rather than assuming. Three fits, on the
# eight patient samples only (controls excluded: they differ technically from
# every patient library, see analyses/x2_confounder.R):
#
#   unpaired_8   ~ CONDITION               what step 11 effectively did
#   paired_4     ~ PATIENT + CONDITION     all four pairs
#   paired_3     ~ PATIENT + CONDITION     pair P3 dropped
#
# Why paired_3 exists: the upstream pipeline compared genotypes called from the
# RNA reads and found that D1_3 and D2_3 agree at 0.7345, below every pair the
# study says is one person and in the range of unrelated individuals (0.7300 to
# 0.8150). If they are different people, calling them a pair makes the model
# attribute a between-person difference to the time effect. A sample that
# appears in no complete pair carries no information about the contrast, so
# dropping the pair is equivalent to giving each of its samples its own patient
# term; this fit shows what that costs and what it changes.
#
# The same script also checks the claim made earlier that this contrast is
# "free of the composition confounder". It is freer than the contrasts against
# controls, not free: library composition still differs inside the pairs.
#
suppressPackageStartupMessages({
    library(DESeq2)
    library(RNASeqPower)
})
source(file.path(Sys.getenv("R_DIR", unset = "R"), "00_common.R"))

banner("x6 - paired design, within-patient contrast")

load(file.path(RDATA_DIR, "11_de.RData"))
load(file.path(RDATA_DIR, "12_annot.RData"))

rd <- read_tsv(file.path(PROJECT_DIR, "inputs", "read_distribution.tsv"),
               show_col_types = FALSE) |> as.data.frame()
rownames(rd) <- rd$sample
rel <- read_tsv(file.path(PROJECT_DIR, "inputs", "relatedness.tsv"),
                show_col_types = FALSE) |> as.data.frame()

raw <- counts(ddsHTSeq, normalized = FALSE)

pat <- sampleTable[sampleTable$CONDITION != "control", ]
pat$CONDITION <- factor(as.character(pat$CONDITION), levels = c("subacute", "acute"))
pat$PATIENT   <- droplevels(factor(pat$PATIENT))
message("patient samples: ", nrow(pat), " from ", nlevels(pat$PATIENT), " patients")

# ---------------------------------------------------------------------------
# 1. Is the pairing real?
# ---------------------------------------------------------------------------
banner("genotype concordance of the claimed pairs")
claimed <- rel[rel$claimed_same == "YES", c("pair", "concordance_nonref", "called", "agrees")]
print(claimed, row.names = FALSE)
best_other <- max(rel$concordance_nonref[rel$claimed_same == "no"])
message("highest concordance among pairs NOT claimed to be one person: ",
        best_other)
conflicts <- claimed$pair[claimed$agrees != "OK"]
message("pairs in conflict with the study's claim: ",
        if (length(conflicts)) paste(conflicts, collapse = ", ") else "none")
save_tab(claimed, "x6_pair_identity")

# ---------------------------------------------------------------------------
# 2. Does composition still differ inside the pairs?
# ---------------------------------------------------------------------------
banner("library composition inside the pairs")

rp_ids <- BIO$ensembl_gene_id[
    grepl("^RP[LS][0-9]", BIO$external_gene_name) & BIO$gene_biotype == "protein_coding"
]
rp_rows <- strip_version(rownames(raw)) %in% rp_ids
rp_fraction <- colSums(raw[rp_rows, , drop = FALSE]) / colSums(raw) * 100

pairs <- do.call(rbind, lapply(levels(pat$PATIENT), function(p) {
    a <- rownames(pat)[pat$PATIENT == p & pat$CONDITION == "acute"]
    s <- rownames(pat)[pat$PATIENT == p & pat$CONDITION == "subacute"]
    data.frame(
        patient    = p, acute = a, subacute = s,
        d_intronic = round(rd[a, "intronic_pct"] - rd[s, "intronic_pct"], 3),
        d_rp       = round(unname(rp_fraction[a] - rp_fraction[s]), 3)
    )
}))
print(pairs, row.names = FALSE)
save_tab(pairs, "x6_pair_composition")

message("pairs where the acute sample has MORE intronic reads: ",
        sum(pairs$d_intronic > 0), " of ", nrow(pairs))
message("pairs where the acute sample has FEWER ribosomal reads: ",
        sum(pairs$d_rp < 0), " of ", nrow(pairs))
message("With four pairs the smallest two-sided sign-test p-value is ",
        signif(2 * 0.5^4, 3), ", so this cannot reach significance however ",
        "consistent it is. The direction is reported, not tested.")

# ---------------------------------------------------------------------------
# 3. Three fits
# ---------------------------------------------------------------------------
banner("fitting")

fit_one <- function(samples, design, label) {
    cd <- pat[samples, , drop = FALSE]
    cd$PATIENT <- droplevels(cd$PATIENT)
    m <- raw[, samples, drop = FALSE]
    keep <- rowSums(m >= MIN_COUNT) > ncol(m) * MIN_FRACTION
    dds <- DESeqDataSetFromMatrix(m[keep, ], colData = cd, design = design)
    dds <- DESeq(dds, quiet = TRUE)
    res <- as.data.frame(results(dds, name = "CONDITION_acute_vs_subacute"))
    res <- res[!is.na(res$padj), , drop = FALSE]
    message("  ", label, ": ", ncol(m), " samples, ", nrow(dds), " genes")
    list(label = label, dds = dds, res = res, n_samples = ncol(m))
}

all8 <- rownames(pat)
no_p3 <- rownames(pat)[pat$PATIENT != "P3"]

fits <- list(
    unpaired_8 = fit_one(all8,  ~ CONDITION,            "unpaired_8"),
    paired_4   = fit_one(all8,  ~ PATIENT + CONDITION,  "paired_4"),
    paired_3   = fit_one(no_p3, ~ PATIENT + CONDITION,  "paired_3")
)

# ---------------------------------------------------------------------------
# 4. What the pairing changes
# ---------------------------------------------------------------------------
banner("comparison")

summarise_fit <- function(f) {
    r <- f$res
    # DESeq2's dispersion uses the same variance model as edgeR's, so its
    # square root is the biological coefficient of variation.
    bcv <- sqrt(median(dispersions(f$dds), na.rm = TRUE))
    data.frame(
        fit            = f$label,
        samples        = f$n_samples,
        genes          = nrow(r),
        median_BCV     = round(bcv, 3),
        padj_lt_0.05   = sum(r$padj < FDR_CUT),
        and_twofold    = sum(r$padj < FDR_CUT & abs(r$log2FoldChange) > LFC_CUT),
        min_padj       = signif(min(r$padj), 3),
        raw_p_lt_0.05_pct = round(100 * mean(r$pvalue < 0.05, na.rm = TRUE), 1)
    )
}
cmp <- do.call(rbind, lapply(fits, summarise_fit))
rownames(cmp) <- NULL
print(cmp, row.names = FALSE)
save_tab(cmp, "x6_fit_comparison")

message("\nunder the null, about 5% of raw p-values fall below 0.05; ",
        "a much larger share means real signal or a miscalibrated model")

ratio <- cmp$median_BCV[cmp$fit == "paired_4"] / cmp$median_BCV[cmp$fit == "unpaired_8"]
message("paired / unpaired BCV ratio (all four pairs): ", round(ratio, 3))
message(if (ratio < 0.8) {
    "  -> pairing removes a substantial share of the variability"
} else if (ratio < 1) {
    "  -> pairing helps a little"
} else {
    paste0("  -> pairing does not reduce variability here: the patients are not ",
           "more alike across time than different from each other")
})

# ---------------------------------------------------------------------------
# 5. Power again, at the variability the paired model actually leaves
# ---------------------------------------------------------------------------
# In a balanced paired design the variance of the effect estimate is 2*s^2/n,
# with s the residual variability after removing each patient's baseline. That
# is the same form rnapower assumes for two independent groups of n, so the
# paired model's residual BCV can be dropped into it. Residual degrees of
# freedom here are very small (3 for paired_4, 2 for paired_3), so the
# dispersion itself is imprecise and these figures carry that uncertainty.
banner("power at the paired model's residual variability")

depth_used <- round(median(rowMeans(counts(ddsHTSeq, normalized = TRUE))[
    rowMeans(counts(ddsHTSeq, normalized = TRUE)) >= MIN_COUNT]))

pw <- do.call(rbind, lapply(fits, function(f) {
    bcv <- sqrt(median(dispersions(f$dds), na.rm = TRUE))
    n   <- if (f$label == "paired_3") 3L else 4L
    n_tests <- nrow(f$res)
    data.frame(
        fit = f$label, pairs_or_per_group = n, BCV = round(bcv, 3),
        power_no_correction = round(rnapower(depth = depth_used, n = n, cv = bcv,
                                             effect = 2.0, alpha = FDR_CUT), 3),
        power_bonferroni    = round(rnapower(depth = depth_used, n = n, cv = bcv,
                                             effect = 2.0, alpha = FDR_CUT / n_tests), 3)
    )
}))
rownames(pw) <- NULL
print(pw, row.names = FALSE)
save_tab(pw, "x6_power")

# ---------------------------------------------------------------------------
# 6. Which genes, if any
# ---------------------------------------------------------------------------
banner("top genes in each fit")

sym <- setNames(BIO$external_gene_name, BIO$ensembl_gene_id)
tops <- list()
for (f in fits) {
    r <- f$res[order(f$res$pvalue), ][seq_len(min(8, nrow(f$res))), ]
    r$symbol <- unname(sym[strip_version(rownames(r))])
    r$ribosomal <- strip_version(rownames(r)) %in% rp_ids
    r$fit <- f$label
    tops[[f$label]] <- r
    message("\n", f$label, ":")
    print(data.frame(symbol = r$symbol, log2FC = round(r$log2FoldChange, 2),
                     p = signif(r$pvalue, 3), padj = signif(r$padj, 3),
                     ribosomal = r$ribosomal), row.names = FALSE)
}
save_tab(do.call(rbind, lapply(tops, function(t) t[, c("fit", "symbol", "log2FoldChange",
                                                       "pvalue", "padj", "ribosomal")])),
         "x6_top_genes")

# The composition check: if library composition leaks into this contrast, the
# ribosomal genes should lean the same way as in the contrasts against
# controls — lower in the sample with more intronic reads, i.e. in acute.
rp_lfc <- sapply(fits, function(f) {
    ids <- rownames(f$res)[strip_version(rownames(f$res)) %in% rp_ids]
    median(f$res[ids, "log2FoldChange"], na.rm = TRUE)
})
message("\nmedian log2FC of the ", sum(rp_rows), " ribosomal protein genes, acute vs subacute:")
print(round(rp_lfc, 3))
message("negative = lower in acute, the direction composition predicts")

# ---------------------------------------------------------------------------
# 7. p-value histograms
# ---------------------------------------------------------------------------
ph <- do.call(rbind, lapply(fits, function(f)
    data.frame(fit = f$label, p = f$res$pvalue[!is.na(f$res$pvalue)])))
save_fig(
    ggplot(ph, aes(p)) +
        geom_histogram(bins = 20, boundary = 0, fill = "grey70", colour = "white") +
        facet_wrap(~ fit, nrow = 1) +
        labs(x = "raw p-value, acute vs subacute", y = "genes",
             title = "Calibration of the within-patient contrast",
             subtitle = "flat = no signal; a spike near 0 = signal; a U shape = miscalibrated model") +
        theme_bw(),
    "x6_pvalue_histograms", width = 12, height = 4
)

fit_results <- lapply(fits, function(f) f$res)
save(fit_results, cmp, pw, pairs, rp_lfc,
     file = file.path(RDATA_DIR, "x6_paired.RData"))
banner("x6 done")
