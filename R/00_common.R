# 00_common.R — paths, sample sheet and helpers shared by every script.
#
# Sourced, never run on its own.
#
# Paths come from the environment when the script is launched through
# run_all.sh (which sources config.sh and exports them), and fall back to the
# same defaults when a script is run by hand in RStudio. config.sh stays the
# single source of truth; these fallbacks only mirror it.

suppressPackageStartupMessages({
    library(tidyverse)
    library(RColorBrewer)
})

# One seed for the whole project. WGCNA, MBASED and every sampling step below
# depend on it, so a rerun reproduces the numbers in SUMMARY.md exactly.
set.seed(42)

env_or <- function(name, default) {
    v <- Sys.getenv(name)
    if (nzchar(v)) v else default
}

HOME_DIR     <- Sys.getenv("HOME")
PROJECT_DIR  <- env_or("PROJECT_DIR",  getwd())
UPSTREAM_DIR <- env_or("UPSTREAM_DIR", file.path(HOME_DIR, "rna_work"))
COUNTS_DIR   <- env_or("COUNTS_DIR",   file.path(UPSTREAM_DIR, "results", "counts"))
BAM_DIR      <- env_or("BAM_DIR",      file.path(UPSTREAM_DIR, "results", "bam"))
ASE_DIR      <- env_or("ASE_DIR",      file.path(UPSTREAM_DIR, "results", "ase"))
REF_DIR      <- env_or("REF_DIR",      file.path(UPSTREAM_DIR, "data", "reference"))
ANNOTATION_GTF <- env_or("ANNOTATION_GTF", file.path(REF_DIR, "chr19.gtf"))

WORK_DIR    <- env_or("WORK_DIR",    file.path(HOME_DIR, "rna_work_r"))
RESULTS_DIR <- env_or("RESULTS_DIR", file.path(WORK_DIR, "results"))
FIG_DIR     <- env_or("FIG_DIR",     file.path(RESULTS_DIR, "figures"))
TAB_DIR     <- env_or("TAB_DIR",     file.path(RESULTS_DIR, "tables"))
RDATA_DIR   <- env_or("RDATA_DIR",   file.path(RESULTS_DIR, "rdata"))
SAMPLE_SHEET <- env_or("SAMPLE_SHEET", file.path(PROJECT_DIR, "samples.tsv"))
THREADS     <- as.integer(env_or("THREADS", "5"))
CHROM       <- env_or("CHROM", "chr19")

for (d in c(RESULTS_DIR, FIG_DIR, TAB_DIR, RDATA_DIR))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)

# --- analysis thresholds ----------------------------------------------------
# The course's own cut-offs, kept here so no script hard-codes its own.
LFC_CUT  <- 1      # |log2 fold change| > 1, i.e. a two-fold difference
FDR_CUT  <- 0.05
MIN_COUNT     <- 10   # prefilter: at least this many counts ...
MIN_FRACTION  <- 0.3  # ... in more than this fraction of samples

mypal <- brewer.pal(8, "Set1")

# One colour per group for every figure that colours by group, so a group looks
# the same wherever it appears. A factor and a character column otherwise get
# different colours from the same palette, because their level order differs.
GROUP_COLOURS <- c(control = "#E41A1C", acute = "#377EB8", subacute = "#4DAF4A")

# --- sample sheet -----------------------------------------------------------
# DESeqDataSetFromHTSeqCount reads column 1 as the sample name and column 2 as
# the file name, so that order is load-bearing rather than cosmetic.
load_sample_table <- function() {
    st <- read_tsv(SAMPLE_SHEET, show_col_types = FALSE) |> as.data.frame()

    st$File <- file.path(paste0(st$sample, ".htseq.counts"))
    present <- file.exists(file.path(COUNTS_DIR, st$File))
    if (!all(present))
        stop("missing count files: ",
             paste(st$sample[!present], collapse = ", "),
             "\nrun scripts/03_counts.sh in the upstream pipeline first")

    # Spaces in a factor level survive into contrast names and break
    # makeContrasts far downstream, so they go now rather than later.
    st$CONDITION <- factor(gsub(" ", "_", st$group))

    # The designs here drop the intercept (~ 0 + CONDITION), so no level is a
    # reference and DESeq2's warning about one is spurious. Putting the control
    # group first costs nothing, silences it, and matches how the contrasts
    # read. Contrasts are named explicitly downstream, so reordering the levels
    # cannot change which comparison is made.
    if ("control" %in% levels(st$CONDITION))
        st$CONDITION <- relevel(st$CONDITION, ref = "control")

    st$PATIENT   <- factor(st$patient)

    st <- st[, c("sample", "File", setdiff(names(st), c("sample", "File")))]
    rownames(st) <- st$sample
    st
}

# --- small helpers ----------------------------------------------------------
save_fig <- function(plot, name, width = 8, height = 6, dpi = 300) {
    path <- file.path(FIG_DIR, paste0(name, ".png"))
    ggsave(path, plot, width = width, height = height, dpi = dpi)
    message("  figure -> ", path)
    invisible(path)
}

save_tab <- function(df, name, rowname_col = NULL) {
    if (!is.null(rowname_col)) df <- tibble::rownames_to_column(df, rowname_col)
    path <- file.path(TAB_DIR, paste0(name, ".tsv"))
    readr::write_tsv(as.data.frame(df), path)
    message("  table  -> ", path)
    invisible(path)
}

# Strip the Ensembl version suffix: ENSG00000004776.15 -> ENSG00000004776.
# biomaRt, WebGestaltR and quantiseqr all reject the versioned form.
strip_version <- function(x) sub("[.].*$", "", x)

# Spearman correlation with a p-value that is never reported as exactly zero.
#
# For small n, cor.test(method = "spearman") computes an exact p-value through
# the AS 89 algorithm, which underflows to 0 for strong correlations. A p-value
# of 0 is not a true statement about any finite sample, and writing one into a
# results table is worse than reporting a slightly conservative bound. When the
# exact routine returns 0, fall back to the asymptotic t approximation and say
# so in the printed line.
spearman_test <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2)
        return(list(rho = NA_real_, p = NA_real_, note = "",
                    text = "not computable (too few points, or no variation)"))

    ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
    p    <- ct$p.value
    note <- ""
    # Only an exact p of zero warrants the fallback. An NA means the test could
    # not be computed at all, and dressing that up as an approximation would
    # report a limitation as if it were a result.
    if (!is.na(p) && p <= 0) {
        ct   <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
        p    <- ct$p.value
        note <- " [exact routine underflowed; asymptotic approximation]"
    }
    rho <- unname(ct$estimate)
    if (is.na(rho) || is.na(p))
        return(list(rho = rho, p = p, note = "",
                    text = "not computable (correlation undefined for this input)"))
    list(rho = rho, p = p, note = note,
         text = paste0("rho = ", round(rho, 3),
                       "  p = ", format.pval(p, digits = 3, eps = 1e-300), note))
}

banner <- function(...) message("\n=== ", ..., " ===")
