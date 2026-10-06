#!/usr/bin/env Rscript
#
# 16_ase_mbased.R — allele-specific expression with MBASED.
#
# The course ran this on one sample. All twelve are run here, because the
# question this project needs answered is not "does this sample show ASE" but
# "does MBASED reach the same verdict the upstream pipeline did".
#
# That verdict was negative and specific. The upstream analysis called 30 059
# sites heterozygous at depth >= 20, and found that only 29% of them carried an
# alternate fraction above 0.4 — the rest were almost certainly homozygous
# sites miscalled as heterozygous. Among the sites that did look like genuine
# heterozygotes, the allelic ratio was 0.5000, i.e. balanced. The imbalance came
# from the genotype calls, and because those genotypes were derived from the
# same reads, the analysis is circular.
#
# MBASED does not check that assumption; it takes the heterozygous calls as
# given and asks whether the two alleles are expressed unequally. So the
# prediction is concrete: MBASED should report a great deal of significant ASE,
# and that signal should track the miscall rate rather than anything
# biological. analyses/x4_ase_crosscheck.R tests exactly that.
#
# Running it on all twelve is what makes the comparison possible: one sample
# would give an anecdote, twelve give a per-sample rate that can be correlated
# with the per-sample genotyping artefact.
#
suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(txdbmaker)
    library(SummarizedExperiment)
    library(MBASED)
    library(BiocParallel)
})
source(file.path(Sys.getenv("R_DIR", unset = "R"), "00_common.R"))

banner("16 - allele-specific expression with MBASED")

NSIM <- as.integer(env_or("MBASED_NSIM", "100000"))
MIN_DEPTH <- 10
message("simulations per gene: ", format(NSIM, big.mark = " "))
message("minimum site depth:   ", MIN_DEPTH)

samples <- load_sample_table()$sample

# ---------------------------------------------------------------------------
# 1. Annotation, cached
# ---------------------------------------------------------------------------
# makeTxDbFromGFF on a 242 MB GTF takes a couple of minutes and the result
# never changes, so it is built once and reused.
banner("gene annotation")

txdb_cache <- file.path(RESULTS_DIR, "chr19_txdb_genes.rds")
if (file.exists(txdb_cache)) {
    genes_gr <- readRDS(txdb_cache)
    message("gene ranges from cache: ", length(genes_gr), " genes")
} else {
    message("building TxDb from ", ANNOTATION_GTF, " (a few minutes)")
    txdb <- suppressWarnings(makeTxDbFromGFF(ANNOTATION_GTF, format = "gtf"))
    genes_gr <- genes(txdb)
    saveRDS(genes_gr, txdb_cache)
    message("gene ranges: ", length(genes_gr), " genes -> cached")
}

# ---------------------------------------------------------------------------
# 2. One sample at a time
# ---------------------------------------------------------------------------
run_one <- function(sid) {
    f <- file.path(ASE_DIR, paste0(sid, ".ase.tsv"))
    if (!file.exists(f)) { message("  ", sid, ": no ASE table, skipped"); return(NULL) }

    ase <- read_tsv(f, show_col_types = FALSE) |>
        filter(totalCount >= MIN_DEPTH,
               !is.na(refAllele), !is.na(altAllele))
    if (!nrow(ase)) { message("  ", sid, ": no site passes the depth filter"); return(NULL) }

    snv_gr <- GRanges(seqnames = ase$contig,
                      ranges   = IRanges(ase$position, width = 1))

    hits <- findOverlaps(snv_gr, genes_gr, ignore.strand = TRUE)
    ase$aseID <- NA_character_
    # A SNV inside two overlapping genes gets the first; findOverlaps returns
    # one row per pair, and assigning by queryHits keeps the last write. Either
    # choice is arbitrary, so it is stated rather than left implicit.
    ase$aseID[queryHits(hits)] <- mcols(genes_gr)$gene_id[subjectHits(hits)]

    keep   <- !is.na(ase$aseID)
    ase    <- ase[keep, ]
    snv_gr <- snv_gr[keep]
    if (!nrow(ase)) { message("  ", sid, ": no site falls inside a gene"); return(NULL) }

    rr <- GRanges(seqnames = seqnames(snv_gr), ranges = ranges(snv_gr),
                  aseID = ase$aseID, allele1 = ase$refAllele, allele2 = ase$altAllele)
    names(rr) <- paste(ase$contig, ase$position, ase$refAllele, ase$altAllele, sep = ":")

    ase_in <- SummarizedExperiment(
        assays = list(
            lociAllele1Counts = matrix(ase$refCount, ncol = 1,
                                       dimnames = list(names(rr), sid)),
            lociAllele2Counts = matrix(ase$altCount, ncol = 1,
                                       dimnames = list(names(rr), sid))
        ),
        rowRanges = rr
    )

    t0 <- Sys.time()
    res <- runMBASED(ase_in, isPhased = FALSE, numSim = NSIM,
                     BPPARAM = SerialParam())
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    out <- data.frame(
        sample = sid,
        geneID = rownames(assays(res)$majorAlleleFrequency),
        MAF    = assays(res)$majorAlleleFrequency[, 1],
        pASE   = assays(res)$pValueASE[, 1],
        pHet   = assays(res)$pValueHeterogeneity[, 1]
    ) |>
        distinct(geneID, .keep_all = TRUE) |>
        mutate(pASE_BH  = p.adjust(pASE, method = "BH"),
               imbalance = abs(MAF - 0.5) * 2,
               sig = pASE_BH < FDR_CUT)

    message("  ", sid, ": ", nrow(ase), " sites in ", nrow(out), " genes, ",
            sum(out$sig), " significant (", sprintf("%.0f%%", 100 * mean(out$sig)),
            "), ", round(secs), " s")
    list(res = out, sites = nrow(ase))
}

banner("running MBASED per sample")
all_out <- list(); site_counts <- integer(0)
for (i in seq_along(samples)) {
    sid <- samples[i]
    r <- run_one(sid)
    if (is.null(r)) next
    all_out[[sid]] <- r$res
    site_counts[sid] <- r$sites
    if (i == 1)
        message("  (", length(samples) - 1, " samples to go at roughly this rate)")
}

if (!length(all_out)) stop("MBASED produced no results for any sample")

ase_all <- bind_rows(all_out)
save_tab(ase_all, "16_mbased_per_gene")

# ---------------------------------------------------------------------------
# 3. Per-sample summary
# ---------------------------------------------------------------------------
banner("per-sample summary")

st <- load_sample_table()
summary_tab <- ase_all |>
    group_by(sample) |>
    summarise(
        genes         = n(),
        significant   = sum(sig),
        pct_significant = round(100 * mean(sig), 1),
        median_MAF    = round(median(MAF), 4),
        median_imbalance = round(median(imbalance), 4),
        .groups = "drop"
    ) |>
    mutate(sites = site_counts[sample],
           group = st$group[match(sample, st$sample)]) |>
    relocate(group, .after = sample) |>
    relocate(sites, .after = group)
print(as.data.frame(summary_tab), row.names = FALSE)
save_tab(summary_tab, "16_mbased_summary")

message("\nacross all samples: ", nrow(ase_all), " gene-sample tests, ",
        sum(ase_all$sig), " significant at BH < ", FDR_CUT,
        " (", sprintf("%.1f%%", 100 * mean(ase_all$sig)), ")")

# A rate this high is the headline. Real ASE affects a few per cent of genes;
# anything approaching a majority is a statement about the input, not about
# allelic regulation. The interpretation is deferred to x4, which has the
# genotype quality numbers to test it against.
if (mean(ase_all$sig) > 0.3)
    message("NOTE: ", sprintf("%.0f%%", 100 * mean(ase_all$sig)),
            " of genes called imbalanced. Published estimates for real ASE sit ",
            "in the low single digits, so this points at the heterozygous ",
            "calls rather than at biology. See analyses/x4_ase_crosscheck.R.")

# The figures are drawn by 16b_ase_figures.R from the tables saved above, so they
# can be redrawn in seconds instead of re-running MBASED (about 30 minutes).

save(ase_all, summary_tab, site_counts, NSIM,
     file = file.path(RDATA_DIR, "16_ase_mbased.RData"))
banner("16 done")
