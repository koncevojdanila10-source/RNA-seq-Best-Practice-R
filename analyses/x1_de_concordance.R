#!/usr/bin/env Rscript
#
# x1_de_concordance.R — how well do DESeq2, edgeR and limma-voom actually agree?
#
# The Venn diagrams in step 11 answer this by eye. A Venn shows the size of an
# overlap but not whether the genes outside it were near misses or genuine
# disagreements, and it says nothing at all about the genes neither tool called.
# This script puts numbers on both:
#
#   1. Spearman correlation of log2 fold change over every gene the tools share,
#      not just the significant ones. This is the honest measure of agreement,
#      because it does not depend on where the threshold happens to fall.
#   2. Jaccard index of the significant sets, per contrast.
#   3. Whether the disagreements concentrate at low expression, which is the
#      usual explanation and one worth confirming rather than assuming.
#
suppressPackageStartupMessages({
    library(ggpubr)
})
source(file.path(Sys.getenv("R_DIR", unset = "R"), "00_common.R"))

banner("x1 - formal agreement between the three tools")

load(file.path(RDATA_DIR, "11_de.RData"))

# The three tools filter slightly differently (1292 against 1284 genes here),
# so every comparison runs on the intersection. Comparing a tool's result to
# another tool's absence would inflate disagreement for free.
lfc_of <- function(nm) {
    d <- deseq_res[[nm]]; e <- edger_res[[nm]]; l <- limma_res[[nm]]
    g <- Reduce(intersect, list(rownames(d), rownames(e), rownames(l)))
    data.frame(
        gene     = g,
        DESeq2   = d[g, "log2FoldChange"],
        edgeR    = e[g, "logFC"],
        limma    = l[g, "logFC"],
        baseMean = d[g, "baseMean"]
    )
}

jaccard <- function(a, b) {
    if (!length(a) && !length(b)) return(NA_real_)
    length(intersect(a, b)) / length(union(a, b))
}

PAIRS <- list(c("DESeq2", "edgeR"), c("DESeq2", "limma"), c("edgeR", "limma"))

conc <- list()
for (nm in names(CONTRASTS)) {
    tab <- lfc_of(nm)
    for (p in PAIRS) {
        conc[[length(conc) + 1]] <- data.frame(
            contrast     = nm,
            pair         = paste(p, collapse = " vs "),
            genes        = nrow(tab),
            spearman_lfc = cor(tab[[p[1]]], tab[[p[2]]],
                               method = "spearman", use = "complete.obs"),
            jaccard_deg  = jaccard(deg_lists[[nm]][[p[1]]], deg_lists[[nm]][[p[2]]]),
            n_deg_1      = length(deg_lists[[nm]][[p[1]]]),
            n_deg_2      = length(deg_lists[[nm]][[p[2]]])
        )
    }
}
concordance <- do.call(rbind, conc)
rownames(concordance) <- NULL
print(concordance, digits = 4)
save_tab(concordance, "x1_concordance")

# --- where the disagreements sit -------------------------------------------
# A gene called by one tool and not by the others is either a borderline case
# or a real difference in how the tools model the data. If the borderline
# explanation is right, those genes should sit at lower expression than the
# genes all three agree on.
banner("expression of agreed vs disputed genes")

agree_tab <- list()
for (nm in names(CONTRASTS)) {
    tab  <- lfc_of(nm)
    sets <- deg_lists[[nm]]
    n_called <- rowSums(vapply(sets, function(s) tab$gene %in% s, logical(nrow(tab))))

    if (max(n_called) == 0) {
        message("  ", nm, ": no significant gene in any tool - nothing to compare")
        next
    }

    d <- data.frame(
        contrast = nm,
        called_by = factor(n_called, levels = 0:3,
                           labels = c("none", "one tool", "two tools", "all three")),
        baseMean = tab$baseMean
    )
    agree_tab[[nm]] <- d

    smry <- aggregate(baseMean ~ called_by, data = d, FUN = median)
    smry$n <- as.integer(table(d$called_by)[as.character(smry$called_by)])
    smry$contrast <- nm
    message("  ", nm, ":")
    print(smry[, c("contrast", "called_by", "n", "baseMean")], row.names = FALSE)

    # Only worth a test when there is something in both groups to compare.
    disputed <- d$baseMean[d$called_by %in% c("one tool", "two tools")]
    agreed   <- d$baseMean[d$called_by == "all three"]
    if (length(disputed) >= 3 && length(agreed) >= 3) {
        w <- wilcox.test(disputed, agreed)
        message("    disputed vs unanimous, Wilcoxon p = ", signif(w$p.value, 3),
                "  (median ", round(median(disputed)), " vs ",
                round(median(agreed)), ")")
    } else {
        message("    too few genes in one of the groups for a test")
    }
}

if (length(agree_tab)) {
    ad <- do.call(rbind, agree_tab)
    save_tab(ad, "x1_called_by_expression")
    save_fig(
        ggplot(ad, aes(called_by, log2(baseMean + 1), fill = called_by)) +
            geom_boxplot(outlier.size = 0.6) +
            facet_wrap(~ contrast) +
            scale_fill_brewer(palette = "Set2") +
            labs(x = "how many of the three tools called the gene",
                 y = "log2(baseMean + 1)",
                 title = "Disagreement between DE tools against expression level") +
            theme_bw() +
            theme(legend.position = "none",
                  axis.text.x = element_text(angle = 30, hjust = 1)),
        "x1_disagreement_vs_expression", width = 11, height = 5
    )
}

# --- the fold-change scatter ------------------------------------------------
plots <- list()
for (nm in names(CONTRASTS)) {
    tab <- lfc_of(nm)
    plots[[nm]] <- ggplot(tab, aes(DESeq2, edgeR)) +
        geom_point(alpha = 0.3, size = 0.7) +
        geom_abline(slope = 1, intercept = 0, colour = "firebrick", linetype = "dashed") +
        labs(title = nm, x = "log2FC (DESeq2)", y = "log2FC (edgeR)") +
        theme_bw()
}
save_fig(ggarrange(plotlist = plots, nrow = 1), "x1_lfc_scatter",
         width = 14, height = 5)

save(concordance, file = file.path(RDATA_DIR, "x1_concordance.RData"))
banner("x1 done")
