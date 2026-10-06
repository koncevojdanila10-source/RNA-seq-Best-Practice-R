#!/usr/bin/env Rscript
#
# 12_annotate_ora.R — gene annotation via biomaRt, and over-representation
# analysis with a reference set that matches the data.
#
# The one deliberate departure from the course script is the ORA background.
# The course passes referenceSet = "genome", which asks whether the DEG list is
# enriched relative to all ~20 000 human genes. Our universe is chromosome 19:
# roughly 3 500 annotated genes, of which only the ones surviving the prefilter
# could ever have been called differentially expressed. The matched question is
# whether the list is enriched relative to the genes that could have been in it,
# so that is the background used here: the genes that actually entered the tests.
#
# That choice was made on principle, and it was not a correction of a visible
# error. analyses/x3_ora_background.R ran both backgrounds and found that, on
# these lists, the choice changes little. The top terms (translation, ribosome,
# rRNA processing) are the same under either; the few terms reported only under
# the genome background are further ribosomal terms just past the threshold, not
# the zinc-finger vocabulary a chromosome-19 artefact would produce; and eight
# random gene sets with no signal produced no enrichment under either background.
# That last test can only rule out a strong inflation (zero of eight leaves a
# per-set rate of up to about 37% inside its 95% interval).
#
suppressPackageStartupMessages({
    library(biomaRt)
    library(WebGestaltR)
})
source(file.path(Sys.getenv("R_DIR", unset = "R"), "00_common.R"))

banner("12 - annotation and over-representation analysis")

load(file.path(RDATA_DIR, "11_de.RData"))

# ---------------------------------------------------------------------------
# 1. biomaRt
# ---------------------------------------------------------------------------
banner("biomaRt")

# Ensembl IDs in a GENCODE GTF carry a version suffix (ENSG00000004776.15).
# biomaRt matches on the unversioned form and returns nothing for the rest,
# silently, so this has to happen before the query rather than after it.
tested_ids   <- rownames(ntd)
tested_clean <- strip_version(tested_ids)

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
BIO <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
    filters    = "ensembl_gene_id",
    values     = tested_clean,
    mart       = ensembl,
    useCache   = FALSE
)
message("queried: ", length(tested_clean), "   returned: ", nrow(BIO))

biotypes <- as.data.frame(sort(table(BIO$gene_biotype), decreasing = TRUE))
names(biotypes) <- c("biotype", "genes")
print(head(biotypes, 12))
save_tab(biotypes, "12_biotypes")

save_fig(
    ggplot(head(biotypes, 12), aes(reorder(biotype, genes), genes)) +
        geom_col(fill = mypal[2]) + coord_flip() +
        labs(x = NULL, y = "genes tested",
             title = "Gene biotypes among the genes that entered the tests") +
        theme_bw(),
    "12_biotypes"
)

# ---------------------------------------------------------------------------
# 2. An expression matrix keyed by gene symbol
# ---------------------------------------------------------------------------
# quantiseqr and WGCNA both want readable names rather than Ensembl IDs, and
# quantiseqr in particular matches its signature by symbol.
banner("expression matrix with gene symbols")

idx <- match(BIO$ensembl_gene_id, tested_clean)
ntd_rename <- ntd[idx, , drop = FALSE]
rownames(ntd_rename) <- BIO$external_gene_name

# Three filters, not the two in the course script. Genes with no symbol come
# back from biomaRt as the empty string, not as NA, so an is.na() test alone
# leaves a row named "" in the matrix — and then a second, and a third, until
# the duplicate filter removes all but one of them and the matrix quietly keeps
# an unnamed gene.
ok <- !is.na(rownames(ntd_rename)) &
       nzchar(rownames(ntd_rename)) &
      !duplicated(rownames(ntd_rename))
message("rows before: ", nrow(ntd_rename), "   after cleaning: ", sum(ok))
ntd_rename <- ntd_rename[ok, , drop = FALSE]

save_tab(as.data.frame(ntd_rename), "12_expression_by_symbol", rowname_col = "symbol")

# ---------------------------------------------------------------------------
# 3. ORA
# ---------------------------------------------------------------------------
banner("over-representation analysis")

# The background: every gene that entered the tests. Anything not in this set
# had no chance of being called significant, so including it in the reference
# would inflate every enrichment.
reference_genes <- unique(strip_version(rownames(ntd)))
message("reference set: ", length(reference_genes), " genes")

# The course selects databases by position in listGeneSet(), which changes
# whenever WebGestalt adds one. Naming them explicitly and checking they exist
# costs one API call and cannot silently analyse the wrong thing.
WANTED_DB <- c(
    "geneontology_Biological_Process",
    "geneontology_Cellular_Component",
    "geneontology_Molecular_Function",
    "pathway_KEGG",
    "pathway_Reactome",
    "network_miRNA_target",
    "network_Transcription_Factor_target"
)

available <- tryCatch(listGeneSet()$name, error = function(e) {
    message("could not reach WebGestalt to list databases: ", conditionMessage(e))
    character(0)
})
DB <- intersect(WANTED_DB, available)
if (!length(DB)) {
    message("no WebGestalt database reachable - skipping ORA. ",
            "Rerun this script when the network is available; ",
            "nothing downstream depends on it.")
} else {
    missing_db <- setdiff(WANTED_DB, DB)
    if (length(missing_db))
        message("not offered by WebGestalt, skipped: ",
                paste(missing_db, collapse = ", "))

    WG_DIR <- file.path(RESULTS_DIR, "WG")
    dir.create(WG_DIR, showWarnings = FALSE, recursive = TRUE)

    # WebGestaltR shells out to `zip` to package its HTML report. That happens
    # *after* the enrichment has been computed, so a zip problem throws an
    # error out of a call that had already done the work we care about.
    #
    # Two separate things have to be true. The binary has to exist, and R has
    # to know where it is: utils::zip() reads the path from R_ZIPCMD and does
    # not fall back to PATH. Conda's R builds ship that variable set to the
    # empty string, which makes utils::zip() fail with "'zip' must be a
    # non-empty character string" even when zip is installed and on PATH.
    if (!nzchar(Sys.getenv("R_ZIPCMD")) && nzchar(Sys.which("zip"))) {
        Sys.setenv(R_ZIPCMD = Sys.which("zip"))
        message("R_ZIPCMD was empty; set to ", Sys.which("zip"))
    }
    have_zip <- nzchar(Sys.getenv("R_ZIPCMD"))
    message("zip binary: ", if (nzchar(Sys.which("zip"))) Sys.which("zip") else "not found",
            "   R_ZIPCMD: ", if (have_zip) Sys.getenv("R_ZIPCMD") else "empty")
    if (!have_zip)
        message("HTML reports will be skipped. The enrichment does not need them.")

    # One call, with a retry that drops the report if the report is what broke.
    # The enrichment result is the deliverable; the HTML is a convenience.
    run_ora <- function(contrast, db, interest, output) {
        WebGestaltR(
            enrichMethod       = "ORA",
            organism           = "hsapiens",
            enrichDatabase     = db,
            interestGene       = interest,
            interestGeneType   = "ensembl_gene_id",
            referenceGene      = reference_genes,
            referenceGeneType  = "ensembl_gene_id",
            isOutput           = output,
            outputDirectory    = WG_DIR,
            projectName        = paste0(contrast, "_", db),
            nThreads           = 1
        )
    }

    ora_summary <- list()
    for (nm in names(CONTRASTS)) {
        interest <- unique(strip_version(deg_lists[[nm]]$DESeq2))
        if (length(interest) < 10) {
            message("  ", nm, ": ", length(interest),
                    " DEGs - too few for ORA, skipped")
            next
        }
        for (db in DB) {
            # An error and an empty result are different outcomes and must not
            # collapse into the same number. Returning NULL on error would make
            # a crashed run indistinguishable from "nothing was enriched", and
            # a fabricated zero is worse than a visible failure.
            res <- tryCatch(
                run_ora(nm, db, interest, have_zip),
                error = function(e)
                    structure(list(msg = conditionMessage(e)), class = "ora_error")
            )

            # If it was the report that failed, the enrichment itself is fine.
            # Redo it without the report rather than throwing the result away.
            if (inherits(res, "ora_error") && grepl("zip", res$msg, ignore.case = TRUE)) {
                message("    report packaging failed - retrying without it")
                res <- tryCatch(
                    run_ora(nm, db, interest, FALSE),
                    error = function(e)
                        structure(list(msg = conditionMessage(e)), class = "ora_error")
                )
            }

            if (inherits(res, "ora_error")) {
                n_hits <- NA_integer_
                message("  ", nm, " / ", db, ": FAILED - ", res$msg)
            } else if (is.null(res) || !nrow(res)) {
                n_hits <- 0L
                message("  ", nm, " / ", db, ": no enriched set at FDR ", FDR_CUT)
            } else {
                n_hits <- nrow(res)
                message("  ", nm, " / ", db, ": ", n_hits, " enriched sets")
                save_tab(res, paste0("12_ora_", nm, "_", db))
            }

            ora_summary[[length(ora_summary) + 1]] <- data.frame(
                contrast = nm, database = db, deg_in = length(interest),
                enriched = n_hits
            )
        }
    }

    if (length(ora_summary)) {
        ora_summary <- do.call(rbind, ora_summary)
        rownames(ora_summary) <- NULL
        print(ora_summary)
        save_tab(ora_summary, "12_ora_summary")

        n_failed <- sum(is.na(ora_summary$enriched))
        if (n_failed)
            message("\nWARNING: ", n_failed, " of ", nrow(ora_summary),
                    " ORA runs failed and are recorded as NA, not as zero. ",
                    "Do not read those rows as a negative result.")
    }
}

# ---------------------------------------------------------------------------
# 4. Persist
# ---------------------------------------------------------------------------
save(BIO, ntd_rename, reference_genes, biotypes,
     file = file.path(RDATA_DIR, "12_annot.RData"))
message("\nenvironment -> ", file.path(RDATA_DIR, "12_annot.RData"))
banner("12 done")
