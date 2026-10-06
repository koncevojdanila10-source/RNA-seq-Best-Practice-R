#!/usr/bin/env Rscript
#
# x2_confounder.R — is the patient-vs-control signal biology or library
# composition?
#
# Step 12 produced an unambiguous enrichment result: in both patient-vs-control
# contrasts the top terms are Translation, Ribosome, rRNA processing and
# SRP-dependent targeting to the ER. Those are ribosomal protein genes. The
# within-patient contrast, on the same twelve libraries, produced nothing.
#
# There is a mundane explanation. The upstream QC found that every control
# library carries a much higher intronic fraction (21-27%) than every patient
# library (7-14%). A library with more pre-mRNA and degraded material spends a
# smaller share of its reads on the mature transcripts of short, very highly
# expressed genes — and ribosomal protein genes are exactly that. So RP genes
# would look systematically lower in controls for reasons that have nothing to
# do with stroke.
#
# This script tests that explanation instead of asserting it:
#
#   1. how completely intronic fraction separates the groups;
#   2. whether the per-sample ribosomal-protein share tracks intronic fraction,
#      including *within* the patients alone, where group membership cannot
#      explain it;
#   3. whether RP genes are formally over-represented among the DEGs;
#   4. what happens to the DEG count when intronic fraction is put into the
#      model as a covariate.
#
# Point 4 comes with a warning attached, made explicit in the output: when a
# covariate is nearly collinear with the group, "correcting" for it removes the
# group effect by construction. The honest reading of a collapse there is that
# the two cannot be separated in this design, not that the biology was refuted.
#
suppressPackageStartupMessages({
    library(DESeq2)
    library(limma)
})
source(file.path(Sys.getenv("R_DIR", unset = "R"), "00_common.R"))

banner("x2 - library composition as a confounder")

load(file.path(RDATA_DIR, "11_de.RData"))
load(file.path(RDATA_DIR, "12_annot.RData"))

rd_path <- file.path(PROJECT_DIR, "inputs", "read_distribution.tsv")
rd <- read_tsv(rd_path, show_col_types = FALSE) |> as.data.frame()
rownames(rd) <- rd$sample
rd <- rd[colnames(ddsHTSeq), ]
stopifnot(!any(is.na(rd$intronic_pct)))

# ---------------------------------------------------------------------------
# 1. How completely does intronic fraction separate the groups?
# ---------------------------------------------------------------------------
banner("group separation by intronic fraction")

is_control  <- rd$group == "control"
max_patient <- max(rd$intronic_pct[!is_control])
min_control <- min(rd$intronic_pct[is_control])

message("patients: ", sprintf("%.1f-%.1f%%", min(rd$intronic_pct[!is_control]), max_patient))
message("controls: ", sprintf("%.1f-%.1f%%", min_control, max(rd$intronic_pct[is_control])))
message("separation is ",
        if (max_patient < min_control) "COMPLETE - no overlap" else "incomplete",
        " (gap ", sprintf("%.1f", min_control - max_patient), " percentage points)")

w <- wilcox.test(rd$intronic_pct[is_control], rd$intronic_pct[!is_control])
message("Wilcoxon control vs patient: p = ", signif(w$p.value, 3),
        "  (n = 4 vs 8, so ", signif(2 / choose(12, 4), 3),
        " is the smallest p this test can return)")

fit_cov <- lm(intronic_pct ~ CONDITION, data = cbind(rd, CONDITION = sampleTable[rownames(rd), "CONDITION"]))
r2 <- summary(fit_cov)$r.squared
message("variance in intronic fraction explained by group: R2 = ", round(r2, 3))

# ---------------------------------------------------------------------------
# 2. Does the ribosomal-protein share track intronic fraction?
# ---------------------------------------------------------------------------
banner("ribosomal protein share against intronic fraction")

# Cytoplasmic ribosomal protein genes: symbols RPL*/RPS* that are protein
# coding. The biotype filter matters — chr19 carries a lot of RP pseudogenes
# with the same symbol stem, and they are not what the enrichment picked up.
rp_symbols <- BIO$external_gene_name[
    grepl("^RP[LS][0-9]", BIO$external_gene_name) &
    BIO$gene_biotype == "protein_coding"
]
rp_ids_clean <- BIO$ensembl_gene_id[BIO$external_gene_name %in% rp_symbols]
raw <- counts(ddsHTSeq, normalized = FALSE)
rp_rows <- strip_version(rownames(raw)) %in% rp_ids_clean

message("ribosomal protein genes among the ", nrow(raw),
        " tested: ", sum(rp_rows),
        " (", paste(head(sort(rp_symbols), 8), collapse = ", "), " ...)")

rd$rp_fraction <- colSums(raw[rp_rows, , drop = FALSE]) / colSums(raw) * 100
rd$CONDITION   <- sampleTable[rownames(rd), "CONDITION"]

print(rd[, c("sample", "group", "intronic_pct", "rp_fraction")], row.names = FALSE)

sp_all <- spearman_test(rd$rp_fraction, rd$intronic_pct)
message("\nall 12 samples:      ", sp_all$text)

# The decisive test. Across all twelve, group membership could produce this
# correlation on its own. Within the patients alone there is no control group
# to drive it, so a correlation surviving here is a continuous technical
# gradient rather than a group difference.
pat <- rd[!is_control, ]
sp_pat <- spearman_test(pat$rp_fraction, pat$intronic_pct)
message("patients only (n = ", nrow(pat), "): ", sp_pat$text)
if (sp_pat$p >= 0.05)
    message("  -> not significant at n = ", nrow(pat),
            "; with this many samples that is weak evidence either way, ",
            "not evidence of absence")

save_tab(rd, "x2_composition")
save_fig(
    # `label` stays inside geom_text_repel rather than the top-level aes:
    # inherited, geom_smooth also receives it and drops it with a warning.
    ggplot(rd, aes(intronic_pct, rp_fraction, colour = CONDITION)) +
        # group = 1 forces a single fit: the colour aesthetic would otherwise
        # split the data by group and draw three separate lines.
        geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "grey40",
                    linewidth = 0.5, formula = y ~ x) +
        geom_point(size = 3) +
        ggrepel::geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE) +
        scale_colour_brewer(palette = "Set1") +
        labs(x = "intronic reads (%)", y = "counts on ribosomal protein genes (%)",
             title = "Ribosomal protein share tracks library composition",
             subtitle = paste0("Spearman rho = ", round(sp_all$rho, 3), ", p = ",
                               format.pval(sp_all$p, digits = 3, eps = 1e-300),
                               "; all 12 samples")) +
        theme_bw(),
    "x2_rp_vs_intronic", width = 8, height = 6
)

# ---------------------------------------------------------------------------
# 3. Are RP genes over-represented among the DEGs?
# ---------------------------------------------------------------------------
banner("over-representation of RP genes among DEGs")

universe <- rownames(deseq_res[[1]])
rp_universe <- strip_version(universe) %in% rp_ids_clean

fisher_tab <- list()
for (nm in names(CONTRASTS)) {
    degs <- deg_lists[[nm]]$DESeq2
    if (!length(degs)) { message("  ", nm, ": no DEGs, skipped"); next }
    is_deg <- universe %in% degs
    tb <- table(RP = rp_universe, DEG = is_deg)
    ft <- fisher.test(tb)
    lfc_rp <- deseq_res[[nm]][universe[rp_universe], "log2FoldChange"]
    fisher_tab[[nm]] <- data.frame(
        contrast     = nm,
        rp_genes     = sum(rp_universe),
        rp_among_deg = sum(rp_universe & is_deg),
        deg_total    = sum(is_deg),
        odds_ratio   = unname(ft$estimate),
        p_value      = ft$p.value,
        median_lfc_rp = median(lfc_rp, na.rm = TRUE)
    )
    message("  ", nm, ": ", sum(rp_universe & is_deg), " of ", sum(rp_universe),
            " RP genes are DEGs, OR = ", round(ft$estimate, 2),
            ", p = ", signif(ft$p.value, 3),
            ", median log2FC = ", round(median(lfc_rp, na.rm = TRUE), 2))
}
if (length(fisher_tab)) {
    rp_enrich <- do.call(rbind, fisher_tab)
    rownames(rp_enrich) <- NULL
    save_tab(rp_enrich, "x2_rp_enrichment")
}

# ---------------------------------------------------------------------------
# 4. Refit with composition in the model
# ---------------------------------------------------------------------------
banner("refitting with intronic fraction as a covariate")

message("NOTE: intronic fraction is ", round(100 * r2),
        "% explained by group. A covariate that close to the grouping removes")
message("      the group effect largely by construction, so a collapse below is")
message("      evidence that the two are inseparable here - not that the")
message("      biological difference has been refuted.")

dds_adj <- ddsHTSeq
colData(dds_adj)$intronic <- scale(rd[colnames(dds_adj), "intronic_pct"])[, 1]
design(dds_adj) <- ~ 0 + CONDITION + intronic
dds_adj <- DESeq(dds_adj, quiet = TRUE)

adj_tab <- list()
for (nm in names(CONTRASTS)) {
    cv <- makeContrasts(contrasts = paste(CONTRASTS[[nm]], collapse = " - "),
                        levels = resultsNames(dds_adj))
    r  <- as.data.frame(results(dds_adj, contrast = as.numeric(cv)))
    r  <- r[!is.na(r$padj), ]
    n_adj <- sum(abs(r$log2FoldChange) > LFC_CUT & r$padj < FDR_CUT)
    n_raw <- length(deg_lists[[nm]]$DESeq2)
    adj_tab[[nm]] <- data.frame(
        contrast = nm, deg_uncorrected = n_raw, deg_with_covariate = n_adj,
        retained_pct = if (n_raw) round(100 * n_adj / n_raw, 1) else NA_real_
    )
    message("  ", nm, ": ", n_raw, " -> ", n_adj, " DEGs")
}
deg_adjusted <- do.call(rbind, adj_tab)
rownames(deg_adjusted) <- NULL
print(deg_adjusted)
save_tab(deg_adjusted, "x2_deg_with_covariate")

save(rd, rp_symbols, deg_adjusted, file = file.path(RDATA_DIR, "x2_confounder.RData"))
banner("x2 done")
