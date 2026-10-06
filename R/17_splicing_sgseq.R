#!/usr/bin/env Rscript
#
# 17_splicing_sgseq.R — splice event detection and differential splicing.
#
# Transcript-level quantification asks which of a set of known transcript
# models best explains the reads, and answers even when the true transcript is
# not in the set. SGSeq asks a smaller question that short reads can actually
# answer: how many reads cross this exon-exon junction, and how many cover this
# exon. Those are observations rather than reconstructions, which is why the
# course reaches for it after transcript assembly rather than instead of it.
#
# Every heavy stage is cached with saveRDS. analyzeFeatures and analyzeVariants
# on twelve BAMs take a long time, and losing an hour to an error in a plotting
# call further down would be an avoidable waste.
#
# A prediction worth writing down before looking: step x2 established that the
# control libraries carry roughly three times the intronic fraction of the
# patient libraries, so the splice-event analysis should separate the groups
# too, for the same technical reason rather than a new biological one.
#
# Section 6 tests it three ways and they do not all agree, which is why each is
# reported with what it can and cannot show:
#
#   novel features    void. analyzeFeatures ran in annotation-based mode, which
#                     does not predict novel events, so there are 8 of them.
#   events detected   supports it. rho = 0.77 with intronic fraction, p = 0.005,
#                     in the expected direction.
#   junction/exon     prediction failed, and the metric turned out not to
#                     measure splicing efficiency at all. See the comment there.
#
suppressPackageStartupMessages({
    library(SGSeq)
    library(limma)
    library(edgeR)
    library(GenomeInfoDb)
})
source(file.path(Sys.getenv("R_DIR", unset = "R"), "00_common.R"))
source(file.path(Sys.getenv("R_DIR", unset = "R"), "plotFeatures_fixed.R"))

banner("17 - alternative splicing with SGSeq")

load(file.path(RDATA_DIR, "11_de.RData"))

CACHE <- file.path(RESULTS_DIR, "sgseq")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

# SGSeq is extremely chatty: analyzeVariants emits one warning per gene whose
# event count exceeds maxnvariant, and the S4 machinery prints a block of
# [updateObject] lines around every object it touches. Hundreds of lines of it
# bury the actual results. They are counted and summarised rather than simply
# hidden, so nothing disappears without a trace.
quiet_run <- function(expr, warn_pattern = "maxnvariant") {
    n <- 0L
    v <- withCallingHandlers(
        suppressMessages(expr),
        warning = function(w) {
            if (grepl(warn_pattern, conditionMessage(w))) {
                n <<- n + 1L
                invokeRestart("muffleWarning")
            }
        }
    )
    if (n) message("  ", n, " genes exceeded maxnvariant; those events were ",
                   "dropped by SGSeq (warnings summarised, not shown)")
    v
}

cached <- function(name, expr) {
    f <- file.path(CACHE, paste0(name, ".rds"))
    if (file.exists(f)) {
        message("  ", name, ": from cache")
        return(readRDS(f))
    }
    t0 <- Sys.time()
    v <- force(expr)
    saveRDS(v, f)
    message("  ", name, ": computed in ",
            round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
    v
}

# ---------------------------------------------------------------------------
# 1. BAM files
# ---------------------------------------------------------------------------
banner("reading BAM headers")

st <- load_sample_table()
bam_files <- file.path(BAM_DIR, paste0(st$sample, ".bam"))
missing <- !file.exists(bam_files)
if (any(missing))
    stop("missing BAM files: ", paste(st$sample[missing], collapse = ", "))

# The two column names are mandatory and SGSeq will not accept substitutes.
samples_df <- data.frame(sample_name = st$sample, file_bam = bam_files,
                         stringsAsFactors = FALSE)

# cores = 1 on purpose. getBamInfo parallelises across files, and the lecturer's
# own class hit failures when the requested core count did not match what was
# available. This step reads headers only; it is not where the time goes.
si <- cached("si", getBamInfo(samples_df, cores = 1))
print(as.data.frame(si)[, c("sample_name", "paired_end", "read_length",
                            "frag_length", "lib_size")], row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. Annotation
# ---------------------------------------------------------------------------
banner("transcript features")

txf <- cached("txf_chr19", {
    txdb <- importTranscripts(file = ANNOTATION_GTF)
    txdb <- keepStandardChromosomes(txdb, species = "Homo sapiens",
                                    pruning.mode = "coarse")
    txdb <- keepSeqlevels(txdb, CHROM, pruning.mode = "coarse")
    convertToTxFeatures(txdb)
})
message("transcript features: ", length(txf))

# ---------------------------------------------------------------------------
# 3. The slow part
# ---------------------------------------------------------------------------
banner("splice graphs, variants and counts")
message("This is the long stage. On twelve BAMs expect tens of minutes.")

sgfc <- cached("sgfc", quiet_run(analyzeFeatures(si, features = txf, cores = THREADS)))
message("features with counts: ", nrow(sgfc))

sgvc_pred <- cached("sgvc_pred", quiet_run(analyzeVariants(sgfc, cores = THREADS)))
sgv <- rowRanges(sgvc_pred)
message("splice variants predicted: ", length(sgv))

sgvc <- cached("sgvc", quiet_run(getSGVariantCounts(sgv, sample_info = si, cores = THREADS)))

x <- counts(sgvc)
rownames(x) <- paste0("event_", seq_len(nrow(x)))
colnames(x) <- st$sample
message("event count matrix: ", nrow(x), " events x ", ncol(x), " samples")

# ---------------------------------------------------------------------------
# 4. Differential splicing
# ---------------------------------------------------------------------------
banner("differential splicing, limma-voom")

design <- model.matrix(~ 0 + CONDITION, data = st)

dge  <- DGEList(counts = x)
keep <- filterByExpr(dge, design = design)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
dge  <- calcNormFactors(dge)
message("events after filtering: ", nrow(dge), " of ", nrow(x))

png(file.path(FIG_DIR, "17_voom_mean_variance.png"),
    width = 1600, height = 1200, res = 200)
v <- voom(dge, design, plot = TRUE)
invisible(dev.off())

fit <- lmFit(v, design)
splice_res <- list(); splice_counts <- list()
for (nm in names(CONTRASTS)) {
    cv  <- makeContrasts(contrasts = paste(CONTRASTS[[nm]], collapse = " - "),
                         levels = design)
    f2  <- eBayes(contrasts.fit(fit, cv))
    tt  <- topTable(f2, sort.by = "P", n = Inf)
    splice_res[[nm]] <- tt
    n_sig <- sum(tt$adj.P.Val < FDR_CUT & abs(tt$logFC) > LFC_CUT)
    splice_counts[[nm]] <- data.frame(contrast = nm, events_tested = nrow(tt),
                                      significant = n_sig)
    message("  ", nm, ": ", n_sig, " significant events of ", nrow(tt))
    save_tab(tt, paste0("17_splicing_", nm), rowname_col = "event")
}
splice_summary <- do.call(rbind, splice_counts)
rownames(splice_summary) <- NULL
print(splice_summary, row.names = FALSE)
save_tab(splice_summary, "17_splicing_summary")

# ---------------------------------------------------------------------------
# 5. The top event
# ---------------------------------------------------------------------------
banner("visualising the top event")

sgv_counts <- rowRanges(sgvc)
top_contrast <- splice_summary$contrast[which.max(splice_summary$significant)]
tt <- splice_res[[top_contrast]]
message("using the strongest contrast: ", top_contrast)

top_feature <- rownames(tt)[1]
top_index   <- match(top_feature, rownames(x))

if (is.na(top_index)) {
    message("could not locate the top event in the count matrix - skipping plots")
} else {
    top_event_id  <- unique(mcols(sgv_counts)$eventID[top_index])[1]
    top_gene_name <- unique(unlist(mcols(sgv_counts)$geneName[top_index]))[1]
    message("top eventID: ", top_event_id)
    message("top gene:    ", if (is.na(top_gene_name)) "unnamed" else top_gene_name)

    # Each plot is wrapped: a failure in one should not cost the others, and
    # suppressMessages keeps the S4 [updateObject] chatter out of the log.
    draw <- function(file, expr, label, width = 1800, height = 1400) {
        png(file.path(FIG_DIR, file), width = width, height = height, res = 200)
        ok <- tryCatch({ suppressMessages(expr); TRUE },
                       error = function(e) {
                           message(label, " failed: ", conditionMessage(e)); FALSE
                       })
        invisible(dev.off())
        if (ok) message("  figure -> ", file.path(FIG_DIR, file))
    }

    draw("17_top_variant.png",
         plotVariants(sgvc, eventID = top_event_id, cex = 0.9),
         "plotVariants")

    if (!is.na(top_gene_name)) {
        sgfc_gene <- SGSeq:::restrictFeatures(sgfc, geneName = top_gene_name)
        draw("17_splice_graph.png",
             plotSpliceGraph(rowRanges(sgfc_gene), color_novel = "red"),
             "plotSpliceGraph", height = 1000)

        # plotFeatures is broken in current SGSeq; the patched version comes
        # from the course and is sourced at the top of this script.
        draw("17_features_heatmap.png",
             plotFeatures_fixed(sgfc, geneName = top_gene_name,
                                color_novel = "red", include = "both"),
             "plotFeatures_fixed")
    }
}

# ---------------------------------------------------------------------------
# 6. Is the splicing signal the composition confounder again?
# ---------------------------------------------------------------------------
# Novel features are the ones SGSeq found in the reads but not in the
# annotation. Reads inside introns generate exactly those. If the count of
# novel features per sample tracks intronic fraction, the splicing differences
# between groups are the same technical artefact seen a third way.
banner("novel features against library composition")

rd <- read_tsv(file.path(PROJECT_DIR, "inputs", "read_distribution.tsv"),
               show_col_types = FALSE) |> as.data.frame()
rownames(rd) <- rd$sample

feat <- rowRanges(sgfc)

# txName is an accessor on SGFeatures, not a column of mcols(). Reaching for
# mcols(feat)$txName returns NULL, and `NULL %in% ...` is logical(0), which
# silently subsets the count matrix to zero rows and reports "0 of 0" as though
# no novel feature existed. A feature is novel when its txName list is empty:
# SGSeq saw it in the reads and the annotation does not contain it.
tx <- SGSeq::txName(feat)
is_novel <- lengths(tx) == 0
stopifnot(length(is_novel) == nrow(sgfc))

fc <- counts(sgfc)
colnames(fc) <- st$sample

novel_share <- colSums(fc[is_novel, , drop = FALSE]) / colSums(fc) * 100
comp <- data.frame(
    sample       = st$sample,
    group        = st$group,
    novel_pct    = round(novel_share[st$sample], 3),
    intronic_pct = rd[st$sample, "intronic_pct"],
    events_detected = colSums(fc > 0)[st$sample]
)
print(comp, row.names = FALSE)

message("novel features: ", sum(is_novel), " of ", length(is_novel))

# A limitation of this test, stated rather than glossed over: analyzeFeatures
# was called with features = txf, which quantifies the supplied annotation and
# does not predict novel events. De novo prediction happens only when no
# annotation is passed. So the novel-feature count here is near zero by
# construction and the correlation below tests nothing. Running SGSeq de novo
# on twelve BAMs would take another hour; the junction ratio computed after it
# tests the same prediction from data already in hand.
if (sum(is_novel) < 100)
    message("  NOTE: only ", sum(is_novel), " novel features exist because ",
            "analyzeFeatures ran in annotation-based mode. The correlation ",
            "below is void, not negative.")

s1 <- spearman_test(comp$novel_pct, comp$intronic_pct)
s2 <- spearman_test(comp$events_detected, comp$intronic_pct)
message("novel feature share vs intronic fraction: ", s1$text,
        if (sum(is_novel) < 100) "   [void, see note above]" else "")
message("events detected      vs intronic fraction: ", s2$text)

# ---------------------------------------------------------------------------
# 6b. Junction reads against exon reads
# ---------------------------------------------------------------------------
# This started as a replacement for the void novel-feature test, on the
# reasoning that a read spanning a junction must come from a spliced
# transcript, so the junction-to-exon ratio would measure how much of each
# library is properly spliced and should fall as the intronic fraction rises.
#
# That reasoning is wrong, and the result below says so: the correlation is
# positive, not negative. The flaw is in the denominator. SGSeq's exon bins
# come from the splice graph built on the annotation, so a read sitting inside
# an intron enters neither the numerator (it crosses no junction) nor the
# denominator (it lies in no annotated exon) — it is simply absent. Intron
# retention is therefore invisible to this ratio by construction.
#
# The number is kept because it was computed and because the failed prediction
# is part of the record, but it is not evidence for or against the confounder.
# What the positive correlation reflects is not established here. The useful
# measurement in this section is the count of detected events above, which does
# track composition in the expected direction.
banner("junction reads against exon reads")

ftype <- as.character(SGSeq::type(feat))
message("feature types: ",
        paste(names(table(ftype)), table(ftype), sep = "=", collapse = ", "))

junction_reads <- colSums(fc[ftype == "J", , drop = FALSE])
exon_reads     <- colSums(fc[ftype == "E", , drop = FALSE])
comp$junction_exon_ratio <- round(junction_reads[st$sample] / exon_reads[st$sample], 4)

print(comp[, c("sample", "group", "intronic_pct", "junction_exon_ratio")],
      row.names = FALSE)
save_tab(comp, "17_novel_vs_composition")

s3 <- spearman_test(comp$junction_exon_ratio, comp$intronic_pct)
message("junction/exon ratio vs intronic fraction: ", s3$text)

save_fig(
    ggplot(comp, aes(intronic_pct, junction_exon_ratio, colour = group)) +
        geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "grey40",
                    linewidth = 0.5, formula = y ~ x) +
        geom_point(size = 3) +
        ggrepel::geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE) +
        scale_colour_manual(values = GROUP_COLOURS) +
        labs(x = "intronic reads (%)", y = "junction reads / exon reads",
             title = "Spliced-read share against library composition",
             subtitle = paste0("Spearman ", s3$text)) +
        theme_bw(),
    "17_junction_ratio_vs_intronic", width = 8, height = 6
)

message("-> Read this with the caveat in the comment above: annotated exon bins",
        " exclude intronic reads, so this ratio cannot see intron retention.")
if (!is.na(s3$rho) && s3$rho > 0)
    message("   The correlation is positive, the opposite of what the ",
            "splicing-efficiency reading would predict. Its cause is not ",
            "established, and it is not used as evidence either way.")

# No figure is drawn for the novel-feature share. It is void for the reason
# given above, and a plot of a measurement that cannot mean anything is worse
# than no plot: it would end up in showcase/ looking like evidence.

save(splice_res, splice_summary, comp, x,
     file = file.path(RDATA_DIR, "17_splicing.RData"))
banner("17 done")
