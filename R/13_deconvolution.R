#!/usr/bin/env Rscript
#
# 13_deconvolution.R — immune cell composition with quanTIseq, and whether this
# dataset can support the question at all.
#
# Two things about this analysis are unusual here, one favourable and one not.
#
# Favourable: the tissue is right. quanTIseq's TIL10 signature describes ten
# immune cell types, and these libraries are PBMC — almost pure immune cells.
# The course demonstrated the method on skeletal muscle biopsies, where the
# signature can only speak about the minority of the sample that is immune.
#
# Unfavourable, and decisive: the gene universe is chromosome 19. TIL10 is a
# whole-genome signature of roughly 170 marker genes, and deconvolution works
# by fitting observed expression against all of them at once. Feeding it the
# fraction that happens to live on one chromosome does not give a noisier
# answer to the same question — it gives a confident answer to a different one,
# because the missing markers are silently treated as absent from the mixture.
#
# So this script measures the signature coverage first and reports it as the
# headline number. The deconvolution is still run, because a demonstration that
# the output is not interpretable is more useful than a blank space, but the
# coverage figure is what decides whether any of it may be believed.
#
suppressPackageStartupMessages({
    library(quantiseqr)
    library(SummarizedExperiment)
})
source(file.path(Sys.getenv("R_DIR", unset = "R"), "00_common.R"))

banner("13 - cell type deconvolution")

load(file.path(RDATA_DIR, "11_de.RData"))
load(file.path(RDATA_DIR, "12_annot.RData"))

message("expression matrix: ", nrow(ntd_rename), " genes x ",
        ncol(ntd_rename), " samples, keyed by gene symbol")

# ---------------------------------------------------------------------------
# 1. How much of the TIL10 signature do we actually have?
# ---------------------------------------------------------------------------
banner("TIL10 signature coverage")

# quantiseqr ships the signature as a file in extdata. The exact name has moved
# between versions, so find it rather than hard-coding it.
ext <- system.file("extdata", package = "quantiseqr")
cand <- list.files(ext, pattern = "TIL10", full.names = TRUE, recursive = TRUE)
message("signature files found: ",
        if (length(cand)) paste(basename(cand), collapse = ", ") else "none")

sig_genes <- character(0)
sig_file  <- cand[grepl("signature", basename(cand), ignore.case = TRUE)][1]
if (is.na(sig_file)) sig_file <- cand[1]

if (length(cand) && !is.na(sig_file)) {
    sig <- tryCatch(read.delim(sig_file, check.names = FALSE),
                    error = function(e) NULL)
    if (!is.null(sig) && ncol(sig) >= 2) {
        # First column holds the gene symbol in every layout this file has had.
        sig_genes <- as.character(sig[[1]])
        message("signature file: ", basename(sig_file),
                " with ", length(sig_genes), " marker genes and ",
                ncol(sig) - 1, " cell types")
    }
}

if (!length(sig_genes)) {
    message("WARNING: could not read the TIL10 signature from the installed ",
            "package. Coverage cannot be measured, so the deconvolution below ",
            "must be treated as uninterpretable rather than merely uncertain.")
    coverage <- NA_real_
} else {
    present  <- intersect(sig_genes, rownames(ntd_rename))
    coverage <- 100 * length(present) / length(sig_genes)
    message("markers present in our chr19 matrix: ", length(present),
            " of ", length(sig_genes),
            "  (", sprintf("%.1f%%", coverage), ")")
    if (length(present))
        message("they are: ", paste(sort(present), collapse = ", "))

    save_tab(data.frame(
        signature      = "TIL10",
        markers_total  = length(sig_genes),
        markers_present = length(present),
        coverage_pct   = round(coverage, 2),
        genes_in_matrix = nrow(ntd_rename)
    ), "13_signature_coverage")
}

# ---------------------------------------------------------------------------
# 2. Run it anyway, on both scales
# ---------------------------------------------------------------------------
# The course passes ntd_rename, which is log2(normalised count + 1).
# quanTIseq was designed for linear TPM-like input, so the log scale is a
# departure worth checking rather than inheriting. If the two runs disagree,
# neither can be quoted without saying which one it is.
banner("running quanTIseq")

run_q <- function(mat, label) {
    message("  ", label, " ...")
    tryCatch(
        quantiseqr::run_quantiseq(
            expression_data = mat,
            signature_matrix = "TIL10",
            is_arraydata = FALSE,
            is_tumordata = FALSE,
            scale_mRNA  = TRUE
        ),
        error = function(e) {
            message("    failed: ", conditionMessage(e))
            NULL
        }
    )
}

ti_log <- run_q(ntd_rename, "log2 scale, as in the course script")
# normTransform is log2(x + 1), so this inverts it exactly.
ti_lin <- run_q(2^ntd_rename - 1, "linear scale, as quanTIseq expects")

cell_cols <- function(x) setdiff(colnames(x), "Sample")

if (!is.null(ti_log)) {
    print(ti_log, digits = 3)
    save_tab(ti_log, "13_quantiseq_log_scale")
}
if (!is.null(ti_lin)) save_tab(ti_lin, "13_quantiseq_linear_scale")

# ---------------------------------------------------------------------------
# 3. Do the two scales agree?
# ---------------------------------------------------------------------------
if (!is.null(ti_log) && !is.null(ti_lin)) {
    banner("log scale against linear scale")
    cc <- intersect(cell_cols(ti_log), cell_cols(ti_lin))
    v1 <- as.numeric(as.matrix(ti_log[, cc]))
    v2 <- as.numeric(as.matrix(ti_lin[, cc]))
    st <- spearman_test(v1, v2)
    message("all sample x cell-type estimates: ", st$text)
    message("median absolute difference: ",
            sprintf("%.4f", median(abs(v1 - v2), na.rm = TRUE)))

    # A correlation between two vectors that are mostly zeros is high for a
    # trivial reason. Say how much of the agreement is just matching zeros
    # before anyone reads it as the method being robust to input scale.
    both_zero <- mean(v1 == 0 & v2 == 0, na.rm = TRUE)
    message("share of estimates that are zero in both runs: ",
            sprintf("%.0f%%", 100 * both_zero))
    if (both_zero > 0.5)
        message("  -> most of this agreement is matching zeros, not the method ",
                "being insensitive to the input scale. It is not reassurance.")
    if (!is.na(st$rho) && st$rho < 0.9)
        message("  -> the two scales do not agree. The choice of input scale ",
                "is then part of the result, which is on its own a reason not ",
                "to quote either number.")
}

# ---------------------------------------------------------------------------
# 4. Figure and verdict
# ---------------------------------------------------------------------------
if (!is.null(ti_log)) {
    banner("figure")
    p <- quantiplot(ti_log) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ggtitle("quanTIseq cell fractions",
                subtitle = if (is.na(coverage)) "signature coverage unknown"
                           else sprintf("TIL10 coverage on chr19: %.1f%% of markers - NOT an estimate of composition", coverage))
    save_fig(p, "13_cell_fractions", width = 10, height = 6)

    fr <- ti_log[, cell_cols(ti_log), drop = FALSE]
    other <- if ("Other" %in% colnames(fr)) mean(fr$Other) else NA_real_
    message("mean 'Other' (unassigned) fraction: ",
            if (is.na(other)) "not reported" else sprintf("%.3g%%", 100 * other))

    # A check that turns the biological objection into a number.
    #
    # PBMC are, by definition, the mononuclear fraction of blood: density
    # gradient separation puts neutrophils in the granulocyte pellet, not in
    # the sample. T cells and monocytes should dominate. Anything claiming the
    # opposite is not an imprecise estimate, it is the wrong answer.
    banner("sanity check against what PBMC are")

    MONONUCLEAR <- c("T.cells.CD4", "T.cells.CD8", "Tregs", "B.cells",
                     "Monocytes", "NK.cells", "Dendritic.cells")
    mono_present <- intersect(MONONUCLEAR, colnames(fr))
    mono_share   <- rowSums(fr[, mono_present, drop = FALSE])
    neut_share   <- if ("Neutrophils" %in% colnames(fr)) fr$Neutrophils
                    else rep(NA_real_, nrow(fr))

    check <- data.frame(
        sample            = ti_log$Sample,
        mononuclear_pct   = round(100 * mono_share, 2),
        neutrophil_pct    = round(100 * neut_share, 2),
        unassigned_pct    = if ("Other" %in% colnames(fr)) signif(100 * fr$Other, 3)
                            else NA_real_
    )
    print(check, row.names = FALSE)
    save_tab(check, "13_pbmc_sanity_check")

    message("mean mononuclear fraction: ", sprintf("%.1f%%", 100 * mean(mono_share)))
    message("mean neutrophil fraction:  ", sprintf("%.1f%%", 100 * mean(neut_share)))
    zeroed <- mono_present[colSums(fr[, mono_present, drop = FALSE]) == 0]
    if (length(zeroed))
        message("cell types estimated at exactly zero in every sample: ",
                paste(zeroed, collapse = ", "))
    if (mean(neut_share, na.rm = TRUE) > mean(mono_share))
        message("-> neutrophils outweigh the mononuclear types in a mononuclear ",
                "cell preparation. The estimate is not merely uncertain, it is ",
                "inverted, and it is stated with no unassigned mass at all.")
}

banner("verdict")
if (!is.na(coverage) && coverage < 50) {
    message("TIL10 coverage is ", sprintf("%.1f%%", coverage),
            ". Deconvolution is NOT assessable on this dataset.")
    message("The missing markers are not treated as missing by the method -- ",
            "they enter the fit as genes with no expression, which is a")
    message("statement about the mixture rather than about our reference. ",
            "The numbers above are reported as a demonstration of the method")
    message("and as evidence for this limit; they are not an estimate of the ",
            "cell composition of these samples.")
    message("Lifting this limit needs whole-genome quantification, which the ",
            "chr19 design rules out: the reads were prefiltered against chr19")
    message("exons before alignment, so the counts for the other 22 ",
            "chromosomes do not exist and cannot be recovered from these files.")
} else if (!is.na(coverage)) {
    message("TIL10 coverage is ", sprintf("%.1f%%", coverage),
            " -- high enough to interpret with care.")
}

save(ti_log, ti_lin, coverage, sig_genes,
     file = file.path(RDATA_DIR, "13_deconvolution.RData"))
banner("13 done")
