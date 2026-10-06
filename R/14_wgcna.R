#!/usr/bin/env Rscript
#
# 14_wgcna.R — weighted gene co-expression network analysis.
#
# The idea WGCNA is built on: instead of testing twenty thousand genes and
# paying for twenty thousand tests, correlate the genes with each other, cut
# the resulting graph into modules, summarise each module by its eigengene, and
# test those instead. A few dozen tests replace tens of thousands, which buys
# back the power that multiple-testing correction takes away — and lets many
# more sample traits be examined than a per-gene design could support.
#
# Two departures from the course script, both forced by this dataset:
#
#   * The soft threshold is chosen from our own scale-independence curve, not
#     copied. The course materials disagree with themselves about it — the
#     script says 6, the slides say 4, and the lecturer picked 7 live from his
#     own plot. It is a property of the data, so it has to be read off ours.
#
#   * The module-trait heatmap correlates modules against library composition
#     rather than against cell-type fractions. Step 13 established that the
#     deconvolution on chr19 is not interpretable, so using it as a trait would
#     propagate a number we have already shown to be meaningless. Intronic
#     fraction and ribosomal-protein share are measured, not inferred, and step
#     x2 showed they are what actually varies across these libraries.
#
# On sample size, plainly: WGCNA is normally applied to twenty samples or more,
# and correlations estimated from twelve are unstable. The course ran it on
# four. Twelve is better and still too few, so the modules here are worth
# looking at and not worth trusting on their own.
#
suppressPackageStartupMessages({
    library(WGCNA)
    library(DESeq2)   # the saved ddsHTSeq is useless without its methods
})
source(file.path(Sys.getenv("R_DIR", unset = "R"), "00_common.R"))

banner("14 - WGCNA co-expression network")

load(file.path(RDATA_DIR, "11_de.RData"))
load(file.path(RDATA_DIR, "12_annot.RData"))

# WGCNA reads strings as strings, and turns them into factors behind your back
# if this is not set. Its own documentation insists on it.
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = THREADS)

# WGCNA wants samples in rows and genes in columns - the transpose of every
# other matrix in this project.
a_ntd <- t(ntd_rename)
message("input: ", nrow(a_ntd), " samples x ", ncol(a_ntd), " genes")

# ---------------------------------------------------------------------------
# 1. Data quality
# ---------------------------------------------------------------------------
banner("checking for bad genes and samples")

gsg <- goodSamplesGenes(a_ntd, verbose = 0)
message("all OK: ", gsg$allOK)
if (!gsg$allOK) {
    message("dropping ", sum(!gsg$goodGenes), " genes and ",
            sum(!gsg$goodSamples), " samples")
    a_ntd <- a_ntd[gsg$goodSamples, gsg$goodGenes]
    message("now: ", nrow(a_ntd), " x ", ncol(a_ntd))
}

png(file.path(FIG_DIR, "14_sample_clustering.png"),
    width = 1800, height = 1200, res = 200)
par(cex = 0.7, mar = c(0, 4, 2, 0))
plot(hclust(dist(a_ntd), method = "average"),
     main = "Sample clustering (outlier check)", sub = "", xlab = "")
invisible(dev.off())

# ---------------------------------------------------------------------------
# 2. Soft threshold
# ---------------------------------------------------------------------------
banner("choosing the soft threshold")

powers <- c(1:15, seq(16, 20, by = 2))
# pickSoftThreshold prints its own table regardless of verbose; we print a
# tidier one below, so its copy is captured and dropped.
invisible(capture.output(
    sft <- pickSoftThreshold(a_ntd, powerVector = powers, verbose = 0)
))

fit <- data.frame(
    power        = sft$fitIndices[, 1],
    signed_R2    = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    connectivity = sft$fitIndices[, 5]
)
print(fit, digits = 3, row.names = FALSE)
save_tab(fit, "14_soft_threshold")

# The rule the lecturer used: the first power whose signed R^2 clears the
# scale-free cut-off. 0.85 is the middle of the usual 0.8-0.9 range. If nothing
# clears it, say so and take the best available rather than pretending.
R2_CUT <- 0.85
ok <- fit$power[fit$signed_R2 >= R2_CUT]
if (length(ok)) {
    SFT_POWER <- min(ok)
    message("first power reaching signed R2 >= ", R2_CUT, ": ", SFT_POWER)
} else {
    SFT_POWER <- fit$power[which.max(fit$signed_R2)]
    message("NO power reaches signed R2 = ", R2_CUT,
            ". Best is ", round(max(fit$signed_R2), 3),
            " at power ", SFT_POWER, ".")
    message("  The network is then not convincingly scale-free, which is a ",
            "statement about having 12 samples, and the modules below inherit ",
            "that weakness.")
}

png(file.path(FIG_DIR, "14_soft_threshold.png"),
    width = 2000, height = 900, res = 200)
par(mfrow = c(1, 2))
plot(fit$power, fit$signed_R2, type = "n",
     xlab = "soft threshold (power)", ylab = "scale-free topology fit, signed R^2",
     main = "Scale independence")
text(fit$power, fit$signed_R2, labels = fit$power, cex = 0.9, col = "blue")
abline(h = R2_CUT, col = "red")
plot(fit$power, fit$connectivity, type = "n",
     xlab = "soft threshold (power)", ylab = "mean connectivity",
     main = "Mean connectivity")
text(fit$power, fit$connectivity, labels = fit$power, cex = 0.9, col = "blue")
invisible(dev.off())

# ---------------------------------------------------------------------------
# 3. Modules
# ---------------------------------------------------------------------------
banner("building the network")

# WGCNA calls cor() with arguments that only its own version accepts
# (weights.x, weights.y, cosine). Loading DESeq2 pulls in an S4 generic for cor
# that does not have them, and blockwiseModules then dies with
# "unused arguments (weights.x = NULL, weights.y = NULL, cosine = FALSE)".
# Binding WGCNA's version for the duration of the network build is the fix the
# package's own FAQ gives. It is restored immediately afterwards, because the
# trait correlations below want Spearman from stats.
cor <- WGCNA::cor

net <- blockwiseModules(
    a_ntd,
    power              = SFT_POWER,
    networkType        = "unsigned",
    TOMType            = "unsigned",
    minModuleSize      = 30,
    reassignThreshold  = 0,
    mergeCutHeight     = 0.25,
    numericLabels      = TRUE,
    pamRespectsDendro  = FALSE,
    saveTOMs           = FALSE,
    randomSeed         = 999,
    verbose            = 0
)

moduleColors <- labels2colors(net$colors)
mod_sizes <- as.data.frame(sort(table(moduleColors), decreasing = TRUE))
names(mod_sizes) <- c("module", "genes")
message("modules found: ", sum(mod_sizes$module != "grey"),
        " plus the grey bin of unassigned genes")
print(mod_sizes, row.names = FALSE)
save_tab(mod_sizes, "14_module_sizes")

# grey is not a module: it collects the genes that joined nothing. Its size is
# a quality measure of the whole network, so it is reported rather than hidden.
grey_n <- mod_sizes$genes[mod_sizes$module == "grey"]
if (length(grey_n))
    message("unassigned (grey): ", grey_n, " genes, ",
            sprintf("%.0f%%", 100 * grey_n / ncol(a_ntd)), " of the input")

png(file.path(FIG_DIR, "14_module_dendrogram.png"),
    width = 2000, height = 1200, res = 200)
plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],
                    "module colours",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = paste0("Gene clustering dendrogram (power = ", SFT_POWER, ")"))
invisible(dev.off())

save_tab(data.frame(gene = colnames(a_ntd), module = moduleColors),
         "14_gene_modules")

# ---------------------------------------------------------------------------
# 4. Modules against measured sample traits
# ---------------------------------------------------------------------------
banner("module-trait relationships")

MEs <- orderMEs(moduleEigengenes(a_ntd, moduleColors)$eigengenes)

# Network building is done; hand cor() back to stats so the Spearman
# correlations below behave as expected.
cor <- stats::cor

rd <- read_tsv(file.path(PROJECT_DIR, "inputs", "read_distribution.tsv"),
               show_col_types = FALSE) |> as.data.frame()
rownames(rd) <- rd$sample

raw <- counts(ddsHTSeq, normalized = FALSE)
rp_ids <- BIO$ensembl_gene_id[
    grepl("^RP[LS][0-9]", BIO$external_gene_name) &
    BIO$gene_biotype == "protein_coding"
]
rp_fraction <- colSums(raw[strip_version(rownames(raw)) %in% rp_ids, , drop = FALSE]) /
               colSums(raw) * 100

st <- sampleTable[rownames(MEs), ]
traits <- data.frame(
    # Group membership as three indicator columns: with three levels there is
    # no single number whose correlation would mean anything.
    acute        = as.integer(st$CONDITION == "acute"),
    subacute     = as.integer(st$CONDITION == "subacute"),
    control      = as.integer(st$CONDITION == "control"),
    intronic_pct = rd[rownames(MEs), "intronic_pct"],
    exonic_pct   = rd[rownames(MEs), "exonic_pct"],
    rp_fraction  = rp_fraction[rownames(MEs)],
    library_size = colSums(raw)[rownames(MEs)] / 1e6,
    row.names    = rownames(MEs)
)

# Spearman: with twelve samples and indicator traits, assuming the linear
# relationship Pearson tests for would be optimistic.
mtc <- cor(MEs, traits, use = "p", method = "spearman")
mtp <- corPvalueStudent(mtc, nrow(traits))

# Multiple testing is not optional here: this is one test per module per trait.
mtq <- matrix(p.adjust(mtp, method = "BH"), nrow = nrow(mtp),
              dimnames = dimnames(mtp))
message("module-trait tests: ", length(mtp),
        ";  raw p < 0.05: ", sum(mtp < 0.05),
        ";  BH q < 0.05: ", sum(mtq < 0.05))

txt <- paste0(signif(mtc, 2), "\n(", signif(mtq, 1), ")")
dim(txt) <- dim(mtc)

png(file.path(FIG_DIR, "14_module_trait_heatmap.png"),
    width = 1800, height = 1400, res = 200)
par(mar = c(8, 9.5, 3, 3))
labeledHeatmap(Matrix = mtc,
               xLabels = colnames(traits), yLabels = names(MEs),
               ySymbols = names(MEs), colorLabels = FALSE,
               colors = blueWhiteRed(50), textMatrix = txt,
               setStdMargins = FALSE, cex.text = 0.55, zlim = c(-1, 1),
               main = "Module eigengenes against measured traits\neach cell: Spearman rho, then (BH q)")
invisible(dev.off())

mt_long <- data.frame(
    module = rep(rownames(mtc), ncol(mtc)),
    trait  = rep(colnames(mtc), each = nrow(mtc)),
    rho    = as.vector(mtc),
    p      = as.vector(mtp),
    q_BH   = as.vector(mtq)
)
mt_long <- mt_long[order(mt_long$q_BH), ]
print(head(mt_long, 12), row.names = FALSE, digits = 3)
save_tab(mt_long, "14_module_trait")

save(net, moduleColors, MEs, traits, mtc, mtp, mtq, mt_long,
     SFT_POWER, fit, a_ntd,
     file = file.path(RDATA_DIR, "14_wgcna.RData"))
banner("14 done")
