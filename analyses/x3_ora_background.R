#!/usr/bin/env Rscript
#
# x3_ora_background.R — what the wrong ORA background would have claimed.
#
# Step 12 ran over-representation analysis against a background of the 1292
# genes that actually entered the tests, and said why: the course passes
# referenceSet = "genome", which asks whether a list is enriched relative to all
# ~20 000 human genes. Our list is drawn from chromosome 19 alone, so against
# that background any property of chromosome 19 itself looks like enrichment.
# Chromosome 19 is unusually rich in KRAB zinc-finger genes, so the expected
# symptom is terms about zinc ions, DNA binding and transcription regulation.
#
# This script puts a number on that instead of leaving it as an argument, in
# two parts.
#
#   A. The real DEG lists, run against the genome as the course does it, and
#      compared with the step 12 results: which terms appear only under the
#      genome background, and how many.
#
#   B. A null experiment. Random gene sets are drawn from genes that are NOT
#      differentially expressed in any contrast, so there is nothing to find.
#      Any enrichment reported for them is manufactured by the background.
#      The same draws go through both backgrounds.
#
# Part B is the stronger of the two, because part A cannot say whether a term
# that only the genome background reports is spurious or just missed by the
# narrower test. A random set has no true signal, so the answer is unambiguous.
#
# Both arms use the same WebGestalt FDR threshold, so the only thing that
# differs between them is the reference set.
#
suppressPackageStartupMessages({
    library(WebGestaltR)
})
source(file.path(Sys.getenv("R_DIR", unset = "R"), "00_common.R"))

banner("x3 - ORA background: chromosome 19 against whole genome")

load(file.path(RDATA_DIR, "11_de.RData"))
load(file.path(RDATA_DIR, "12_annot.RData"))

universe <- unique(strip_version(rownames(ntd)))
message("genes in our background: ", length(universe))

# WebGestaltR packages its HTML report with zip, and utils::zip() reads the
# path from R_ZIPCMD, which conda's R leaves empty. See step 12 for the story.
if (!nzchar(Sys.getenv("R_ZIPCMD")) && nzchar(Sys.which("zip")))
    Sys.setenv(R_ZIPCMD = Sys.which("zip"))
have_zip <- nzchar(Sys.getenv("R_ZIPCMD"))

WG_DIR <- file.path(RESULTS_DIR, "WG_x3")
dir.create(WG_DIR, showWarnings = FALSE, recursive = TRUE)

DB <- c("geneontology_Biological_Process", "geneontology_Cellular_Component",
        "geneontology_Molecular_Function", "pathway_KEGG", "pathway_Reactome",
        "network_miRNA_target", "network_Transcription_Factor_target")

quiet <- function(expr) {
    invisible(capture.output(v <- suppressWarnings(suppressMessages(expr))))
    v
}

# One ORA call. An error and an empty result are kept apart, as in step 12: a
# crashed run must never be recorded as "nothing enriched".
ora <- function(interest, db, project, background) {
    args <- list(enrichMethod = "ORA", organism = "hsapiens", enrichDatabase = db,
                 interestGene = interest, interestGeneType = "ensembl_gene_id",
                 outputDirectory = WG_DIR, projectName = project, nThreads = 1)
    if (background == "ours") {
        args$referenceGene     <- universe
        args$referenceGeneType <- "ensembl_gene_id"
    } else {
        args$referenceSet <- "genome"
    }
    # Finished calls are cached by project name. Every one of these is a round
    # trip to a remote server, and a late failure should not cost the earlier
    # queries again. Errors are never cached, so a failed call is retried.
    cache <- file.path(WG_DIR, paste0(project, ".rds"))
    if (file.exists(cache)) return(readRDS(cache))

    go <- function(out) {
        a <- args; a$isOutput <- out
        tryCatch(quiet(do.call(WebGestaltR, a)),
                 error = function(e) structure(list(msg = conditionMessage(e)),
                                               class = "ora_error"))
    }
    res <- go(have_zip)
    if (inherits(res, "ora_error") && grepl("zip", res$msg, ignore.case = TRUE))
        res <- go(FALSE)
    if (inherits(res, "ora_error")) {
        message("    FAILED: ", res$msg)
        return(list(status = "error", table = NULL, n = NA_integer_))
    }
    tab <- if (is.null(res) || !nrow(res)) NULL else as.data.frame(res)
    out <- list(status = "ok", table = tab, n = if (is.null(tab)) 0L else nrow(tab))
    saveRDS(out, cache)
    out
}

empty_tab <- data.frame(geneSet = character(), description = character(),
                        enrichmentRatio = numeric(), FDR = numeric())

# ---------------------------------------------------------------------------
# A. The real DEG lists, genome background against ours
# ---------------------------------------------------------------------------
banner("A. real DEG lists")

s12 <- read_tsv(file.path(TAB_DIR, "12_ora_summary.tsv"), show_col_types = FALSE)
if (any(is.na(s12$enriched)))
    stop("step 12 recorded failed ORA runs; rerun it before comparing against it")

CONTR <- c("acute_vs_control", "subacute_vs_control")
cmp <- list(); only_g <- list(); ours_desc <- character()

for (nm in CONTR) {
    interest <- unique(strip_version(deg_lists[[nm]]$DESeq2))
    message(nm, ": ", length(interest), " DEGs")
    for (db in DB) {
        g <- ora(interest, db, paste0("x3_", nm, "_", db), "genome")
        f <- file.path(TAB_DIR, paste0("12_ora_", nm, "_", db, ".tsv"))
        ours <- if (file.exists(f)) as.data.frame(read_tsv(f, show_col_types = FALSE))
                else empty_tab

        ours_desc <- c(ours_desc, as.character(ours$description))
        n_step12 <- s12$enriched[s12$contrast == nm & s12$database == db]
        if (length(n_step12) && n_step12 != nrow(ours))
            warning("step 12 summary says ", n_step12, " for ", nm, "/", db,
                    " but its table has ", nrow(ours), " rows")

        gt <- if (g$status == "ok" && !is.null(g$table)) g$table else empty_tab
        shared <- intersect(ours$geneSet, gt$geneSet)
        only_genome_ids <- setdiff(gt$geneSet, ours$geneSet)
        message(sprintf("  %-38s ours %3d | genome %3s | shared %3d | only genome %3s",
                        db, nrow(ours),
                        if (is.na(g$n)) "NA" else g$n, length(shared),
                        if (is.na(g$n)) "NA" else length(only_genome_ids)))

        cmp[[length(cmp) + 1]] <- data.frame(
            contrast = nm, database = db, degs = length(interest),
            ours = nrow(ours), genome = g$n, shared = length(shared),
            only_ours = length(setdiff(ours$geneSet, gt$geneSet)),
            only_genome = if (is.na(g$n)) NA_integer_ else length(only_genome_ids)
        )
        if (length(only_genome_ids)) {
            # WebGestalt's result columns differ between database types (the
            # network ones carry no description), so the columns wanted are
            # created when missing instead of assumed, or the rows cannot be
            # stacked afterwards.
            want <- c("geneSet", "description", "enrichmentRatio", "FDR")
            for (w in setdiff(want, names(gt))) gt[[w]] <- NA
            o <- gt[gt$geneSet %in% only_genome_ids, want]
            o$description <- as.character(o$description)
            o$contrast <- nm; o$database <- db
            only_g[[length(only_g) + 1]] <- o
        }
    }
}
cmp <- do.call(rbind, cmp); rownames(cmp) <- NULL
print(cmp, row.names = FALSE)
save_tab(cmp, "x3_background_comparison")

if (length(only_g)) {
    og <- do.call(rbind, only_g)
    og <- og[order(og$FDR), ]
    save_tab(og, "x3_only_genome_terms")

    # Chosen before looking at the list: the vocabulary of a chromosome-19
    # zinc-finger cluster. A fraction well above what the same words account
    # for among the terms found under our own background would support the
    # account; a similar fraction would not.
    pat <- "zinc|KRAB|DNA.binding|nucleic acid|transcription|metal ion|regulation of gene expression"
    is_zf   <- grepl(pat, og$description, ignore.case = TRUE)
    base_zf <- mean(grepl(pat, ours_desc, ignore.case = TRUE))
    message("\nterms reported ONLY under the genome background: ", nrow(og))
    message("  matching the zinc-finger / transcription vocabulary: ",
            sum(is_zf), " (", sprintf("%.0f%%", 100 * mean(is_zf)), ")")
    message("  the same words among the terms found under OUR background: ",
            sprintf("%.0f%%", 100 * base_zf), " (", length(ours_desc), " terms)")
    message("\nstrongest of them:")
    print(head(data.frame(contrast = og$contrast,
                          term = substr(og$description, 1, 52),
                          ratio = round(og$enrichmentRatio, 1),
                          FDR = signif(og$FDR, 2),
                          zf_vocabulary = is_zf), 14), row.names = FALSE)
} else {
    message("\nno term appears only under the genome background")
    og <- NULL
}

# ---------------------------------------------------------------------------
# B. Null experiment: random sets with nothing to find
# ---------------------------------------------------------------------------
banner("B. random gene sets with no signal")

all_degs <- unique(strip_version(unlist(lapply(deg_lists, function(l) l$DESeq2))))
pool <- setdiff(universe, all_degs)
SET_SIZE <- min(150L, floor(length(pool) / 3))
N_DRAWS  <- 8L
message("genes that are DE in no contrast: ", length(pool),
        "  (excluded from nothing else; ", length(all_degs), " DEGs removed)")
message("draws: ", N_DRAWS, " of ", SET_SIZE, " genes each, GO Biological Process")
message("draws share genes (they come from one pool), so they are not ",
        "independent; they are repeated looks at the same null, not eight studies")

set.seed(42)
rand <- list()
for (d in seq_len(N_DRAWS)) {
    pick <- sample(pool, SET_SIZE)
    for (bg in c("ours", "genome")) {
        r <- ora(pick, "geneontology_Biological_Process", paste0("x3_random", d, "_", bg), bg)
        top <- if (!is.null(r$table)) r$table$description[which.min(r$table$FDR)] else NA_character_
        rand[[length(rand) + 1]] <- data.frame(
            draw = d, background = bg, enriched = r$n,
            top_term = if (is.na(top)) "" else substr(top, 1, 48)
        )
    }
    last <- tail(rand, 2)
    message(sprintf("  draw %d: ours %s sets | genome %s sets", d,
                    last[[1]]$enriched, last[[2]]$enriched))
}
rand <- do.call(rbind, rand); rownames(rand) <- NULL
save_tab(rand, "x3_random_sets")

ok <- rand[!is.na(rand$enriched), ]
rsum <- do.call(rbind, lapply(split(ok, ok$background), function(x) data.frame(
    background = x$background[1], draws = nrow(x),
    draws_with_any_enrichment = sum(x$enriched > 0),
    median_sets = median(x$enriched), max_sets = max(x$enriched)
)))
rownames(rsum) <- NULL
message("")
print(rsum, row.names = FALSE)
save_tab(rsum, "x3_random_summary")

gen_top <- rand[rand$background == "genome" & rand$top_term != "", ]
if (nrow(gen_top)) {
    message("\nbest-scoring term in each random draw under the genome background:")
    print(gen_top[, c("draw", "enriched", "top_term")], row.names = FALSE)
}

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
long <- rbind(
    data.frame(cmp[, c("contrast", "database")], background = "chr19 (ours)", sets = cmp$ours),
    data.frame(cmp[, c("contrast", "database")], background = "whole genome", sets = cmp$genome)
)
db_labels <- c(geneontology_Biological_Process = "GO BP", geneontology_Cellular_Component = "GO CC",
               geneontology_Molecular_Function = "GO MF", pathway_KEGG = "KEGG",
               pathway_Reactome = "Reactome", network_miRNA_target = "miRNA targets",
               network_Transcription_Factor_target = "TF targets")
long$database <- factor(db_labels[long$database], levels = unname(db_labels))
save_fig(
    ggplot(long, aes(database, sets, fill = background)) +
        geom_col(position = "dodge") +
        facet_wrap(~ contrast) +
        scale_fill_brewer(palette = "Set1") +
        labs(x = NULL, y = "enriched gene sets (FDR < 0.05)",
             title = "Same DEG lists, two ORA backgrounds") +
        theme_bw() + theme(axis.text.x = element_text(angle = 40, hjust = 1)),
    "x3_background_comparison", width = 12, height = 5
)
save_fig(
    ggplot(ok, aes(factor(draw), enriched, fill = background)) +
        geom_col(position = "dodge") +
        scale_fill_brewer(palette = "Set1") +
        labs(x = "random draw", y = "enriched GO BP sets (FDR < 0.05)",
             title = "Random gene sets with no signal",
             subtitle = "anything reported here is produced by the background, not the genes") +
        theme_bw(),
    "x3_random_sets", width = 8, height = 5
)

# ---------------------------------------------------------------------------
# Verdict
# ---------------------------------------------------------------------------
banner("verdict")
r_ours <- rsum[rsum$background == "ours", ]
r_gen  <- rsum[rsum$background == "genome", ]
if (nrow(r_ours) && nrow(r_gen)) {
    message("Random sets with no signal, ", N_DRAWS, " draws:")
    message("  chr19 background:  ", r_ours$draws_with_any_enrichment,
            " draws reported an enriched term (median ", r_ours$median_sets, " sets)")
    message("  genome background: ", r_gen$draws_with_any_enrichment,
            " draws reported an enriched term (median ", r_gen$median_sets, " sets)")
    if (r_gen$draws_with_any_enrichment > r_ours$draws_with_any_enrichment)
        message("The genome background reports enrichment in gene sets that ",
                "contain no signal, so for chromosome-19 lists its output ",
                "cannot be trusted.")
    else
        message("The two backgrounds do not differ on random sets here, so ",
                "this run does not show the genome background inflating ",
                "enrichment.")
}

save(cmp, og, rand, rsum, file = file.path(RDATA_DIR, "x3_ora_background.RData"))
banner("x3 done")
