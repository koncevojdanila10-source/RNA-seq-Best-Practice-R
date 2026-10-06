#!/usr/bin/env Rscript
#
# 11_de_three_tools.R — differential expression with DESeq2, edgeR and
# limma-voom, and the comparison between them.
#
# The course ran one contrast on four samples. This dataset has twelve samples
# in three groups, so every pairwise contrast is run, and each tool is fitted
# once and queried three times — which is how these packages are meant to be
# used, since the dispersion estimate is shared across contrasts.
#
#   acute    vs control     stroke at 24 h against healthy donors
#   subacute vs control     the same patients at day 7 against healthy donors
#   acute    vs subacute    within-patient, the time course
#
# A caveat that belongs beside the results rather than after them: the control
# libraries differ technically from the patient libraries. Every control
# carries a higher intronic fraction (23.6% against 8-10%) and a much lower
# duplication rate, and that separation is complete. Both control contrasts are
# therefore confounded by library composition; analyses/x2_confounder.R
# quantifies it. The within-patient contrast is the one to trust.
#
suppressPackageStartupMessages({
    library(DESeq2)
    library(edgeR)
    library(limma)
    library(vsn)
    library(ggpubr)
    library(pheatmap)
    library(EnhancedVolcano)
    library(ggVennDiagram)
})
source(file.path(Sys.getenv("R_DIR", unset = "R"), "00_common.R"))

banner("11 - differential expression, three tools")

sampleTable <- load_sample_table()
message("samples: ", nrow(sampleTable),
        " in groups: ", paste(levels(sampleTable$CONDITION), collapse = ", "))

# Every pairwise comparison between the three groups, named so the name can be
# used directly as a file-name stem.
CONTRASTS <- list(
    acute_vs_control    = c("CONDITIONacute",    "CONDITIONcontrol"),
    subacute_vs_control = c("CONDITIONsubacute", "CONDITIONcontrol"),
    acute_vs_subacute   = c("CONDITIONacute",    "CONDITIONsubacute")
)

contrast_string <- function(nm) paste(CONTRASTS[[nm]], collapse = " - ")

# ---------------------------------------------------------------------------
# 1. DESeq2
# ---------------------------------------------------------------------------
banner("DESeq2")

# ~ 0 + CONDITION drops the intercept. Keeping it would silently make whichever
# level sorts first alphabetically the reference for everything, and with three
# groups every contrast would then be read against that one level.
ddsHTSeq <- DESeqDataSetFromHTSeqCount(
    sampleTable = sampleTable,
    directory   = COUNTS_DIR,
    design      = ~ 0 + CONDITION
)
message("genes before prefiltering: ", nrow(ddsHTSeq))

# Prefilter: at least MIN_COUNT counts in more than MIN_FRACTION of samples.
# The threshold is empirical and the trade-off runs both ways. Set the fraction
# too high and a gene expressed in one group only is discarded before any test
# sees it; set it too low and almost nothing is gained in statistical power.
keep <- rowSums(counts(ddsHTSeq) >= MIN_COUNT) > ncol(ddsHTSeq) * MIN_FRACTION
ddsHTSeq <- ddsHTSeq[keep, ]
message("genes after prefiltering:  ", nrow(ddsHTSeq))

ddsHTSeq <- DESeq(ddsHTSeq)
COU <- counts(ddsHTSeq, normalized = TRUE)
message("resultsNames: ", paste(resultsNames(ddsHTSeq), collapse = ", "))

deseq_res <- list()
for (nm in names(CONTRASTS)) {
    cv <- makeContrasts(contrasts = contrast_string(nm),
                        levels = resultsNames(ddsHTSeq))
    r  <- as.data.frame(results(ddsHTSeq, contrast = as.numeric(cv)))
    # Genes DESeq2 filtered internally carry NA in padj and never entered the
    # multiple-testing correction, so they are dropped rather than counted.
    deseq_res[[nm]] <- r[!is.na(r$padj), ]
    message("  ", nm, ": ", nrow(deseq_res[[nm]]), " genes with a padj")
}

# ---------------------------------------------------------------------------
# 2. Transformations, and whether rlog earns its cost
# ---------------------------------------------------------------------------
banner("variance-stabilising transformations")

common     <- rownames(deseq_res[[1]])
rld        <- assay(rlog(ddsHTSeq))[common, ]
ntd        <- assay(normTransform(ddsHTSeq))[common, ]
log.counts <- log2(counts(ddsHTSeq, normalized = FALSE) + 1)[common, ]

p1 <- meanSdPlot(log.counts, ranks = FALSE, plot = FALSE)
p2 <- meanSdPlot(ntd,        ranks = FALSE, plot = FALSE)
p3 <- meanSdPlot(rld,        ranks = FALSE, plot = FALSE)
save_fig(
    ggarrange(p1$gg, p2$gg, p3$gg, nrow = 1,
              labels = c("log2", "log2 + size factors", "rlog")),
    "11_normalisation_comparison", width = 14, height = 5
)

# Library sizes by group: the plot that shows whether the groups were sequenced
# comparably, before any model is fitted.
libsize <- data.frame(
    sample = colnames(ddsHTSeq),
    reads  = colSums(counts(ddsHTSeq)) / 1e6,
    group  = sampleTable[colnames(ddsHTSeq), "CONDITION"]
)
save_fig(
    ggplot(libsize, aes(sample, reads, fill = group)) +
        geom_col() +
        scale_fill_brewer(palette = "Set1") +
        labs(y = "assigned reads (millions)", x = NULL,
             title = "Library size by sample") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)),
    "11_library_sizes"
)

# ---------------------------------------------------------------------------
# 3. edgeR — quasi-likelihood F test
# ---------------------------------------------------------------------------
banner("edgeR")

design_edger <- model.matrix(~ 0 + CONDITION, data = sampleTable)
colnames(design_edger) <- gsub(" ", "_", colnames(design_edger))

y <- DGEList(counts(ddsHTSeq))
y <- y[filterByExpr(y, design = design_edger), , keep.lib.sizes = FALSE]
message("genes after filterByExpr: ", nrow(y))

# edgeR does not carry library sizes along the way DESeq2 does: they have to be
# recomputed after filtering, or TMM normalises against the pre-filter totals.
y <- calcNormFactors(y)
y$samples$lib.size <- colSums(y$counts)
y <- estimateDisp(y, design_edger, tagwise = TRUE)

fit_edger <- glmQLFit(y, design_edger)
edger_res <- list()
for (nm in names(CONTRASTS)) {
    cv  <- makeContrasts(contrasts = contrast_string(nm),
                         levels = colnames(design_edger))
    lrt <- glmQLFTest(fit_edger, contrast = cv)
    edger_res[[nm]] <- as.data.frame(topTags(lrt, sort.by = "none", n = nrow(lrt)))
}

# ---------------------------------------------------------------------------
# 4. limma-voom
# ---------------------------------------------------------------------------
banner("limma-voom")

# voom models the mean-variance relationship explicitly. That is what lets a
# linear model built for microarray intensities work on integer counts.
png(file.path(FIG_DIR, "11_voom_mean_variance.png"),
    width = 1600, height = 1200, res = 200)
y_voom <- voom(y, design_edger, plot = TRUE)
invisible(dev.off())

fit_limma <- lmFit(y_voom, design_edger)
limma_res <- list()
for (nm in names(CONTRASTS)) {
    cv  <- makeContrasts(contrasts = contrast_string(nm),
                         levels = colnames(design_edger))
    tmp <- eBayes(contrasts.fit(fit_limma, cv))
    limma_res[[nm]] <- topTable(tmp, sort.by = "none", n = Inf)
}

# ---------------------------------------------------------------------------
# 5. Significant genes, and the Venn diagrams
# ---------------------------------------------------------------------------
banner("agreement between the three tools")

sig_genes <- function(nm, adjusted = TRUE) {
    d <- deseq_res[[nm]]; e <- edger_res[[nm]]; l <- limma_res[[nm]]
    if (adjusted) {
        list(DESeq2 = rownames(d)[abs(d$log2FoldChange) > LFC_CUT & d$padj      < FDR_CUT],
             edgeR  = rownames(e)[abs(e$logFC)          > LFC_CUT & e$FDR       < FDR_CUT],
             limma  = rownames(l)[abs(l$logFC)          > LFC_CUT & l$adj.P.Val < FDR_CUT])
    } else {
        list(DESeq2 = rownames(d)[abs(d$log2FoldChange) > LFC_CUT & d$pvalue  < FDR_CUT],
             edgeR  = rownames(e)[abs(e$logFC)          > LFC_CUT & e$PValue  < FDR_CUT],
             limma  = rownames(l)[abs(l$logFC)          > LFC_CUT & l$P.Value < FDR_CUT])
    }
}

deg_lists <- list(); deg_lists_raw <- list(); counts_tbl <- list()
for (nm in names(CONTRASTS)) {
    deg_lists[[nm]]     <- sig_genes(nm, adjusted = TRUE)
    deg_lists_raw[[nm]] <- sig_genes(nm, adjusted = FALSE)

    v_fdr <- ggVennDiagram(deg_lists[[nm]], label_alpha = 0) +
        scale_fill_gradient(low = "white", high = "red") +
        ggtitle(paste0(nm, " - FDR < ", FDR_CUT)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    v_raw <- ggVennDiagram(deg_lists_raw[[nm]], label_alpha = 0) +
        scale_fill_gradient(low = "white", high = "orange") +
        ggtitle(paste0(nm, " - raw p < ", FDR_CUT, ", no correction")) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    save_fig(ggarrange(v_fdr, v_raw, ncol = 2), paste0("11_venn_", nm),
             width = 12, height = 6)

    counts_tbl[[nm]] <- data.frame(
        contrast = nm,
        tool     = names(deg_lists[[nm]]),
        fdr      = lengths(deg_lists[[nm]]),
        raw_p    = lengths(deg_lists_raw[[nm]])
    )
}
deg_counts <- do.call(rbind, counts_tbl)
rownames(deg_counts) <- NULL
print(deg_counts)
save_tab(deg_counts, "11_deg_counts")

# ---------------------------------------------------------------------------
# 6. Sample-level views
# ---------------------------------------------------------------------------
banner("MDS, PCA, heatmap, volcano")

png(file.path(FIG_DIR, "11_mds.png"), width = 1600, height = 1400, res = 200)
plotMDS(rld, col = mypal[sampleTable[colnames(rld), "CONDITION"]],
        main = "MDS (base R)")
legend("topleft", fill = mypal[seq_along(levels(sampleTable$CONDITION))],
       legend = levels(sampleTable$CONDITION))
invisible(dev.off())

pcaData <- plotPCA(rlog(ddsHTSeq), intgroup = "CONDITION", returnData = TRUE)
pv <- round(100 * attr(pcaData, "percentVar"))
save_fig(
    ggplot(pcaData, aes(PC1, PC2, color = CONDITION)) +
        geom_point(size = 4, alpha = 0.8) +
        xlab(paste0("PC1: ", pv[1], "% variance")) +
        ylab(paste0("PC2: ", pv[2], "% variance")) +
        scale_color_brewer(palette = "Set1") +
        theme_bw() + ggtitle("PCA"),
    "11_pca"
)

# The heatmap is drawn on genes already selected for differing between the
# groups, so the clustering it shows follows from that selection and is not
# independent evidence for it.
for (nm in names(CONTRASTS)) {
    g <- deg_lists[[nm]]$DESeq2
    if (length(g) < 2) { message("  ", nm, ": too few DEGs for a heatmap"); next }
    anno <- as.data.frame(colData(ddsHTSeq)[, "CONDITION", drop = FALSE])
    png(file.path(FIG_DIR, paste0("11_heatmap_", nm, ".png")),
        width = 1600, height = 1800, res = 200)
    pheatmap(rld[g, ], cluster_rows = TRUE, cluster_cols = TRUE,
             show_rownames = FALSE, annotation_col = anno,
             main = paste0("DEGs (DESeq2) - ", nm))
    invisible(dev.off())
}

for (nm in names(CONTRASTS)) {
    save_fig(
        EnhancedVolcano(deseq_res[[nm]],
                        lab = rownames(deseq_res[[nm]]),
                        x = "log2FoldChange", y = "padj",
                        pCutoff = FDR_CUT, FCcutoff = LFC_CUT,
                        title = paste0("DESeq2 - ", nm), subtitle = NULL),
        paste0("11_volcano_", nm), width = 9, height = 8
    )
}

# ---------------------------------------------------------------------------
# 7. Persist
# ---------------------------------------------------------------------------
for (nm in names(CONTRASTS)) {
    save_tab(deseq_res[[nm]], paste0("11_deseq2_", nm), rowname_col = "gene_id")
    save_tab(edger_res[[nm]], paste0("11_edger_",  nm), rowname_col = "gene_id")
    save_tab(limma_res[[nm]], paste0("11_limma_",  nm), rowname_col = "gene_id")
}

save(sampleTable, ddsHTSeq, COU, rld, ntd, log.counts,
     deseq_res, edger_res, limma_res, deg_lists, deg_lists_raw, deg_counts,
     CONTRASTS, design_edger, y, y_voom,
     file = file.path(RDATA_DIR, "11_de.RData"))
message("\nenvironment -> ", file.path(RDATA_DIR, "11_de.RData"))
banner("11 done")
