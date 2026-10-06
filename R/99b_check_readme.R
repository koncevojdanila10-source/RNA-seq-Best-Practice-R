#!/usr/bin/env Rscript
#
# 99b_check_readme.R — does README.md still say what the tables say?
#
# README.md quotes figures from this analysis, and prose does not update itself.
# A number typed into a sentence goes stale the moment a step is re-run with
# different data, a different seed or a different threshold, and nothing would
# notice. So each quoted figure is recomputed here from the tables the pipeline
# wrote, formatted the way the README writes it, and looked up in the README text.
#
# A mismatch is reported and not repaired. The README is prose written by a
# person, and whether the sentence or the table is the one that needs to change
# is a judgement this script should not make.
#
# Typographic characters (minus sign, en dash, arrow, superscript two) are mapped
# to ASCII on the README side so the lookups stay plain ASCII strings.
#
source(file.path(Sys.getenv("R_DIR", unset = "R"), "00_common.R"))

banner("99b - README consistency check")

readme_path <- file.path(PROJECT_DIR, "README.md")
if (!file.exists(readme_path)) stop("no README.md at ", readme_path)

tab <- function(name) {
    f <- file.path(TAB_DIR, paste0(name, ".tsv"))
    if (!file.exists(f)) return(NULL)
    as.data.frame(read_tsv(f, show_col_types = FALSE))
}

ascii <- function(x) {
    x <- gsub("−", "-",  x, fixed = TRUE)   # minus sign
    x <- gsub("–", "-",  x, fixed = TRUE)   # en dash
    x <- gsub("→", "->", x, fixed = TRUE)   # right arrow
    x <- gsub("²", "2",  x, fixed = TRUE)   # superscript two
    gsub("[[:space:]]+", " ", x)
}
txt <- ascii(paste(readLines(readme_path, warn = FALSE, encoding = "UTF-8"), collapse = " "))

checks <- list()
has <- function(label, needle) {
    checks[[length(checks) + 1]] <<- data.frame(
        check = label, expected = needle, found = grepl(needle, txt, fixed = TRUE))
}
holds <- function(label, cond) {
    checks[[length(checks) + 1]] <<- data.frame(
        check = label, expected = "(condition on the tables)", found = isTRUE(cond))
}
n3   <- function(x) formatC(x, digits = 3, format = "g")
big  <- function(x) formatC(x, big.mark = ",", format = "d")

skipped <- character()
need <- function(name) {
    t <- tab(name)
    if (is.null(t)) skipped <<- c(skipped, name)
    t
}

# --- library composition -----------------------------------------------------
comp <- need("x2_composition")
if (!is.null(comp)) {
    ctrl <- comp$group == "control"
    has("intronic range, patients",
        sprintf("%.1f-%.1f%% in patients", min(comp$intronic_pct[!ctrl]), max(comp$intronic_pct[!ctrl])))
    has("intronic range, controls",
        sprintf("%.1f-%.1f%% in controls", min(comp$intronic_pct[ctrl]), max(comp$intronic_pct[ctrl])))
    has("R2 of intronic fraction on group",
        sprintf("R2 = %.2f", summary(lm(intronic_pct ~ group, data = comp))$r.squared))
    has("rho within the patients",
        sprintf("rho = %.3f", spearman_test(comp$rp_fraction[!ctrl], comp$intronic_pct[!ctrl])$rho))
}

rpe <- need("x2_rp_enrichment")
if (!is.null(rpe)) {
    a <- rpe[rpe$contrast == "acute_vs_control", ]
    s <- rpe[rpe$contrast == "subacute_vs_control", ]
    has("odds ratios", sprintf("odds ratio %.1f and %.1f", a$odds_ratio, s$odds_ratio))
    holds("11 of the 12 ribosomal genes are DEGs in both contrasts",
          a$rp_genes == 12 && a$rp_among_deg == 11 && s$rp_among_deg == 11)
}

cv <- need("x2_deg_with_covariate")
if (!is.null(cv))
    has("DEG count with the covariate",
        sprintf("(%d -> %d and %d -> %d genes)", cv$deg_uncorrected[1], cv$deg_with_covariate[1],
                cv$deg_uncorrected[2], cv$deg_with_covariate[2]))

pcv <- need("11_pca_variance")
if (!is.null(pcv)) has("PC1 variance", sprintf("%d%% of the variance", pcv$percent_variance[1]))

conc <- need("x1_concordance")
if (!is.null(conc)) has("tool agreement", sprintf("Spearman %.2f-1.00", min(conc$spearman_lfc)))

# --- networks ----------------------------------------------------------------
mt <- need("14_module_trait")
if (!is.null(mt))
    has("main module against intronic fraction",
        sprintf("rho = %.3f", mt$rho[mt$module == "MEturquoise" & mt$trait == "intronic_pct"]))

mb <- need("15_module_trait_bionero")
if (!is.null(mb))
    has("BioNERO against intronic fraction",
        sprintf("implementation: -%.3f)", max(abs(mb$rho[mb$trait == "intronic_pct"]))))

ms <- need("14_module_sizes"); pa <- need("15_partition_agreement")
if (!is.null(ms) && !is.null(pa)) {
    has("size of the main module",
        sprintf("%d of the %s genes", ms$genes[ms$module == "turquoise"], big(pa$genes_compared[1])))
    has("ARI at matched settings", sprintf("adjusted Rand index %.2f", pa$adjusted_rand[1]))
    has("ARI at default settings", sprintf("%.2f at each package", pa$adjusted_rand[2]))
}

# --- splicing, deconvolution, ASE --------------------------------------------
sp <- need("17_splicing_summary")
if (!is.null(sp))
    has("differentially spliced events",
        sprintf("%s and %s differentially spliced events", big(sp$significant[1]), big(sp$significant[2])))

nv <- need("17_novel_vs_composition")
if (!is.null(nv))
    has("events detected against intronic fraction",
        sprintf("rho = %.3f", spearman_test(nv$events_detected, nv$intronic_pct)$rho))

cov13 <- need("13_signature_coverage")
if (!is.null(cov13)) {
    has("TIL10 markers present", sprintf("%d of the %d marker genes", cov13$markers_present, cov13$markers_total))
    has("TIL10 coverage", sprintf("(%.2f%%)", cov13$coverage_pct))
}
pb <- need("13_pbmc_sanity_check")
if (!is.null(pb)) has("mean neutrophil share", sprintf("%.1f%% neutrophils", mean(pb$neutrophil_pct)))

ase <- need("16_mbased_summary")
if (!is.null(ase)) {
    has("share of genes called imbalanced",
        sprintf("calls %.0f%% of genes imbalanced", 100 * sum(ase$significant) / sum(ase$genes)))
    holds("sites per sample between 800 and 7,000", min(ase$sites) >= 800 && max(ase$sites) <= 7000)
}
xp <- need("x4_per_sample")
if (!is.null(xp))
    has("ASE rate against doubtful calls",
        sprintf("rho = %.3f)", spearman_test(xp$pct_significant, xp$pct_below_0.4)$rho))
xa <- need("x4_significant_vs_allele_fraction")
if (!is.null(xa))
    has("lowest alternate fraction",
        sprintf("%.3f against %.3f", xa$median_min_alt[xa$sig], xa$median_min_alt[!xa$sig]))

# --- the within-patient contrast ---------------------------------------------
fc <- need("x6_fit_comparison")
if (!is.null(fc)) {
    p <- function(f) n3(fc$min_padj[fc$fit == f])
    has("smallest adjusted p in the three fits",
        sprintf("%s, %s, %s", p("unpaired_8"), p("paired_4"), p("paired_3")))
    r <- function(f) fc$raw_p_lt_0.05_pct[fc$fit == f]
    has("small-p excess with pair P3",
        sprintf("%.0f-%.0f%% below 0.05", r("unpaired_8"), r("paired_4")))
    has("small-p share without pair P3", sprintf("(%.1f%%)", r("paired_3")))
}
pal <- need("18_power_by_alpha")
if (!is.null(pal))
    has("power bracket, edgeR variability",
        sprintf("%.2f to %.2f", min(pal$power_n4), max(pal$power_n4)))
pw <- need("x6_power")
if (!is.null(pw))
    has("power bracket, paired fits",
        sprintf("%.2f to %.2f", min(pw$power_bonferroni), max(pw$power_no_correction)))
pid <- need("x6_pair_identity")
if (!is.null(pid))
    has("genotype concordance of the conflicting pair",
        sprintf("(%.4f)", pid$concordance_nonref[pid$agrees != "OK"]))

# --- report ------------------------------------------------------------------
res <- do.call(rbind, checks); rownames(res) <- NULL
save_tab(res, "99b_readme_check")

message(sprintf("\n%d of %d figures quoted in README.md confirmed against the tables",
                sum(res$found), nrow(res)))
if (length(skipped))
    message("tables not found, so their figures were not checked: ",
            paste(unique(skipped), collapse = ", "))

bad <- res[!res$found, ]
if (nrow(bad)) {
    message("\nNOT FOUND in README.md (the table or the sentence is out of date):")
    print(bad[, c("check", "expected")], row.names = FALSE)
    banner("99b finished with mismatches")
    quit(status = 1)
}
banner("99b done - every checked figure matches")
