#!/usr/bin/env Rscript
#
# 99_summarize_results.R — builds SUMMARY.md and fills showcase/.
#
# Every number in SUMMARY.md is read or recomputed here from the tables the
# pipeline wrote, and nothing is typed in. A figure that a table does not carry
# directly (a correlation, an R-squared, a Jaccard index) is recomputed from the
# table's rows, so it can be traced to its inputs. A section whose input is
# missing is skipped and listed at the end; it is never filled with a guess.
#
# The prose that interprets these numbers lives in README.md. This file holds
# only tables and one-line captions.
#
suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
})
source(file.path(Sys.getenv("R_DIR", unset = "R"), "00_common.R"))

banner("99 - SUMMARY.md and showcase/")
SHOW_DIR <- env_or("SHOWCASE_DIR", file.path(PROJECT_DIR, "showcase"))

missing_inputs <- character()

tab <- function(name) {
    f <- file.path(TAB_DIR, paste0(name, ".tsv"))
    if (!file.exists(f)) { missing_inputs <<- c(missing_inputs, name); return(NULL) }
    as.data.frame(read_tsv(f, show_col_types = FALSE))
}

fmt_num <- function(v) vapply(v, function(x) {
    if (is.na(x)) "NA"
    else if (abs(x - round(x)) < 1e-9 && abs(x) < 1e12)
        format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
    else formatC(x, digits = 4, format = "g")
}, character(1))

md_table <- function(df) {
    if (is.null(df) || !nrow(df)) return("_not available_")
    cells <- lapply(df, function(col) {
        out <- if (is.numeric(col)) fmt_num(col) else as.character(col)
        out[is.na(out)] <- "NA"
        gsub("|", "/", out, fixed = TRUE)
    })
    head <- paste0("| ", paste(names(df), collapse = " | "), " |")
    sep  <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
    rows <- do.call(paste, c(unname(cells), sep = " | "))
    c(head, sep, paste0("| ", rows, " |"))
}

out <- character()
add <- function(...) out <<- c(out, paste0(...))
section <- function(title, caption = NULL) {
    add(""); add("## ", title); add("")
    if (!is.null(caption)) { add(caption); add("") }
}
put <- function(df) { out <<- c(out, md_table(df), "") }
na_skip <- function(...) any(vapply(list(...), is.null, logical(1)))

add("# Results")
add("")
add("Every number below is read or recomputed from the pipeline's own output by")
add("`R/99_summarize_results.R`. Regenerate with `bash summarize_results.sh`.")

# ---------------------------------------------------------------------------
section("Differential expression, three tools",
        "Genes with |log2 FC| > 1 and FDR < 0.05 (`fdr`), or with raw p < 0.05 and no correction (`raw_p`).")
deg <- tab("11_deg_counts")
if (!is.null(deg)) {
    put(deg |> pivot_wider(names_from = tool, values_from = c(fdr, raw_p)))
}

conc <- tab("x1_concordance")
if (!is.null(conc)) {
    add("Agreement between tools on log2 fold change over every shared gene (Spearman), and on the significant sets (Jaccard; NA when both lists are empty).")
    add("")
    put(conc[, c("contrast", "pair", "genes", "spearman_lfc", "jaccard_deg")])
}

# ---------------------------------------------------------------------------
section("Library composition differs between groups",
        "Intronic fraction and the share of reads on ribosomal protein genes, by group.")
comp <- tab("x2_composition")
if (!is.null(comp)) {
    put(comp |> group_by(group) |>
        summarise(samples = n(),
                  intronic_min = min(intronic_pct), intronic_max = max(intronic_pct),
                  rp_share_min = min(rp_fraction),  rp_share_max = max(rp_fraction),
                  .groups = "drop"))

    ctrl <- comp$group == "control"
    sep_complete <- max(comp$intronic_pct[!ctrl]) < min(comp$intronic_pct[ctrl])
    r2 <- summary(lm(intronic_pct ~ group, data = comp))$r.squared
    sp_all <- spearman_test(comp$rp_fraction, comp$intronic_pct)
    sp_pat <- spearman_test(comp$rp_fraction[!ctrl], comp$intronic_pct[!ctrl])
    put(data.frame(
        quantity = c("groups separated completely by intronic fraction",
                     "variance of intronic fraction explained by group (R2)",
                     "ribosomal share vs intronic fraction, all 12 samples",
                     "ribosomal share vs intronic fraction, patients only (n = 8)"),
        value = c(if (sep_complete) "yes" else "no", round(r2, 3),
                  sp_all$text, sp_pat$text)
    ))
}

rpe <- tab("x2_rp_enrichment")
if (!is.null(rpe)) {
    add("Ribosomal protein genes among the DEGs (Fisher's exact test).")
    add("")
    put(rpe)
}
cov <- tab("x2_deg_with_covariate")
if (!is.null(cov)) {
    add("DEG count when intronic fraction is added to the model. The covariate is largely explained by group, so a collapse reflects inseparability, not refutation.")
    add("")
    put(cov)
}

# ---------------------------------------------------------------------------
section("Annotation and enrichment",
        "Gene biotypes among the genes that entered the tests.")
bt <- tab("12_biotypes")
if (!is.null(bt)) put(bt)

ora <- tab("12_ora_summary")
if (!is.null(ora)) {
    add("Enriched gene sets per contrast and database (WebGestalt ORA, FDR < 0.05, background = genes tested). `NA` would mean the run failed.")
    add("")
    put(ora)
}

bgc <- tab("x3_background_comparison")
if (!is.null(bgc)) {
    add("The same DEG lists against the chr19 background and against the whole genome.")
    add("")
    put(bgc)
}
bgr <- tab("x3_random_summary")
if (!is.null(bgr)) {
    add("Random gene sets drawn from genes that are not DE in any contrast, through both backgrounds. Zero draws with enrichment means neither background produced a false one.")
    add("")
    put(bgr)
}

# ---------------------------------------------------------------------------
section("Cell-type deconvolution",
        "Coverage of the quanTIseq TIL10 signature on chromosome 19, and the estimate against what a PBMC sample is.")
cov13 <- tab("13_signature_coverage")
if (!is.null(cov13)) put(cov13)
pb <- tab("13_pbmc_sanity_check")
if (!is.null(pb)) {
    put(pb)
    put(data.frame(
        quantity = c("mean mononuclear fraction", "mean neutrophil fraction"),
        percent  = c(round(mean(pb$mononuclear_pct), 2), round(mean(pb$neutrophil_pct), 2))
    ))
}

# ---------------------------------------------------------------------------
section("Co-expression networks",
        "WGCNA soft threshold: the first power whose signed scale-free R2 reaches 0.85.")
sft <- tab("14_soft_threshold")
if (!is.null(sft)) {
    chosen <- min(sft$power[sft$signed_R2 >= 0.85])
    put(data.frame(quantity = c("chosen power", "signed R2 at that power"),
                   value = c(chosen, round(sft$signed_R2[sft$power == chosen], 3))))
}
ms <- tab("14_module_sizes")
if (!is.null(ms)) { add("WGCNA modules."); add(""); put(ms) }
mt <- tab("14_module_trait")
if (!is.null(mt)) {
    add("Strongest module-trait associations (Spearman, Benjamini-Hochberg q).")
    add("")
    put(head(mt[order(mt$q_BH), ], 8))
}
pa <- tab("15_partition_agreement")
if (!is.null(pa)) {
    add("Agreement between the WGCNA and BioNERO module partitions (adjusted Rand index; 1 = identical, 0 = chance).")
    add("")
    put(pa)
}
mb <- tab("15_module_trait_bionero")
if (!is.null(mb) && !is.null(mt)) {
    best <- function(d) round(max(abs(d$rho[d$trait == "intronic_pct"])), 3)
    put(data.frame(package = c("WGCNA", "BioNERO"),
                   strongest_abs_rho_with_intronic_fraction = c(best(mt), best(mb))))
}
zf <- tab("15_hub_zinc_finger_share")
if (!is.null(zf)) { add("Hub genes from the zinc-finger families."); add(""); put(zf) }

# ---------------------------------------------------------------------------
section("Allele-specific expression",
        "MBASED, all twelve samples, genes with FDR < 0.05.")
ase <- tab("16_mbased_summary")
if (!is.null(ase)) {
    put(ase[, c("sample", "group", "sites", "genes", "significant", "pct_significant")])
    put(data.frame(quantity = "genes called imbalanced, all samples pooled (%)",
                   value = round(100 * sum(ase$significant) / sum(ase$genes), 1)))
}
fl <- tab("16_p_at_floor")
if (!is.null(fl)) {
    add("Share of gene-level p-values below the simulation resolution (1 / numSim); these are reported as p = 0 by MBASED and drawn at the limit.")
    add("")
    put(fl)
}
xp <- tab("x4_per_sample")
if (!is.null(xp)) {
    s1 <- spearman_test(xp$pct_significant, xp$pct_below_0.4)
    put(data.frame(quantity = "per-sample % imbalanced vs % of doubtful heterozygous calls",
                   value = s1$text))
}
xa <- tab("x4_significant_vs_allele_fraction")
if (!is.null(xa)) {
    add("Allele fractions of genes MBASED calls imbalanced against the rest.")
    add("")
    put(xa)
}
xb <- tab("x4_mbased_vs_binomial")
if (!is.null(xb)) {
    xb$MBASED <- as.logical(xb$MBASED); xb$binomial <- as.logical(xb$binomial)
    both   <- sum(xb$Freq[xb$MBASED & xb$binomial])
    either <- sum(xb$Freq[xb$MBASED | xb$binomial])
    put(data.frame(
        quantity = c("gene-sample pairs compared", "called by MBASED", "called by the binomial test",
                     "Jaccard index between the two"),
        value = c(sum(xb$Freq), sum(xb$Freq[xb$MBASED]), sum(xb$Freq[xb$binomial]),
                  round(both / either, 3))
    ))
}

# ---------------------------------------------------------------------------
section("Alternative splicing",
        "Differential splicing events (limma-voom, |logFC| > 1 and FDR < 0.05).")
sp <- tab("17_splicing_summary")
if (!is.null(sp)) put(sp)
nv <- tab("17_novel_vs_composition")
if (!is.null(nv)) {
    s2 <- spearman_test(nv$events_detected, nv$intronic_pct)
    s3 <- spearman_test(nv$junction_exon_ratio, nv$intronic_pct)
    put(data.frame(
        quantity = c("events detected per sample vs intronic fraction",
                     "junction / exon read ratio vs intronic fraction"),
        value = c(s2$text, s3$text)
    ))
}

# ---------------------------------------------------------------------------
section("Statistical power",
        "Within-group coefficient of variation measured from these counts, and from edgeR's model.")
cvm <- tab("18_measured_cv")
if (!is.null(cvm)) put(cvm)
cve <- tab("18_cv_estimators")
if (!is.null(cve)) put(cve)
pal <- tab("18_power_by_alpha")
if (!is.null(pal)) {
    add("Power to detect a twofold change at the median gene with four replicates per group, with and without a multiple-testing correction.")
    add("")
    put(pal)
}
pgv <- tab("18_power_by_gene_variability")
if (!is.null(pgv)) put(pgv)

section("Within-patient contrast, paired design",
        "Acute against subacute in the eight patient samples, with and without modelling the patient pairing, and without patient P3.")
fc <- tab("x6_fit_comparison")
if (!is.null(fc)) put(fc)
pw <- tab("x6_power")
if (!is.null(pw)) put(pw)
pc <- tab("x6_pair_composition")
if (!is.null(pc)) {
    add("Per-pair difference in library composition, acute minus subacute.")
    add("")
    put(pc)
}
pid <- tab("x6_pair_identity")
if (!is.null(pid)) {
    add("Genotype concordance of the pairs the study reports as one person (from the upstream pipeline).")
    add("")
    put(pid)
}

if (length(missing_inputs)) {
    add(""); add("## Sections skipped"); add("")
    add("These inputs were not found, so their sections are absent rather than estimated: ",
        paste0("`", unique(missing_inputs), "`", collapse = ", "), ".")
}

writeLines(out, file.path(PROJECT_DIR, "SUMMARY.md"))
message("wrote ", file.path(PROJECT_DIR, "SUMMARY.md"), " (", length(out), " lines)")
if (length(missing_inputs))
    message("skipped for missing input: ", paste(unique(missing_inputs), collapse = ", "))

# ---------------------------------------------------------------------------
# showcase/: small figures and tables only. Anything large stays out of git.
# ---------------------------------------------------------------------------
banner("showcase")
for (d in c("figures", "tables"))
    dir.create(file.path(SHOW_DIR, d), showWarnings = FALSE, recursive = TRUE)

KEEP_FIGURES <- c(
    "11_pca", "11_normalisation_comparison", "11_venn_acute_vs_control",
    "11_venn_acute_vs_subacute", "12_biotypes", "13_cell_fractions",
    "14_soft_threshold", "14_module_trait_heatmap", "16_ase_volcano", "16_ase_qq",
    "17_junction_ratio_vs_intronic", "17_features_heatmap", "18_power_curves",
    "x1_lfc_scatter", "x1_disagreement_vs_expression", "x2_rp_vs_intronic",
    "x3_background_comparison", "x3_random_sets", "x4_ase_vs_miscall",
    "x4_min_alt_by_significance", "x6_pvalue_histograms"
)
copied <- 0L
for (f in KEEP_FIGURES) {
    src <- file.path(FIG_DIR, paste0(f, ".png"))
    if (file.exists(src)) {
        file.copy(src, file.path(SHOW_DIR, "figures", basename(src)), overwrite = TRUE)
        copied <- copied + 1L
    } else message("  figure not found: ", f)
}

# Tables under 60 KB, minus the per-gene matrices that are large and not a
# summary of anything.
skip_tab <- c("12_expression_by_symbol", "14_gene_modules", "15_gene_modules",
              "15_gene_modules_matched", "16_mbased_per_gene")
all_tabs <- list.files(TAB_DIR, pattern = "\\.tsv$", full.names = TRUE)
tcopied <- 0L
for (f in all_tabs) {
    stem <- sub("\\.tsv$", "", basename(f))
    if (stem %in% skip_tab || file.size(f) > 60 * 1024) next
    file.copy(f, file.path(SHOW_DIR, "tables", basename(f)), overwrite = TRUE)
    tcopied <- tcopied + 1L
}
message("copied ", copied, " figures and ", tcopied, " tables into showcase/")
banner("99 done")
