#!/usr/bin/env Rscript
#
# 16b_ase_figures.R — the MBASED figures, drawn from the tables step 16 saved.
#
# Split out of step 16 for two reasons. MBASED takes about half an hour on
# twelve samples and a plotting mistake should not cost that. And the first
# version of these figures was wrong in a way worth recording.
#
# MBASED estimates p-values by simulation, so the smallest value it can report
# is 1 / numSim. With 10^5 simulations, a gene whose observed imbalance was never
# matched in any simulation comes back with p = 0, which means "below 1e-05", not
# "zero". The first version of these plots replaced the zeros with the smallest
# positive floating-point number, 2e-308, and drew -log10 of that. Most of the
# genes then sat on a line at 308 and the axis said nothing. The floor here is
# the simulation's own resolution, and the figures say how many genes sit on it.
#
suppressPackageStartupMessages({
    library(ggrepel)
})
source(file.path(Sys.getenv("R_DIR", unset = "R"), "00_common.R"))

banner("16b - MBASED figures")

load(file.path(RDATA_DIR, "12_annot.RData"))

ase_all <- as.data.frame(read_tsv(file.path(TAB_DIR, "16_mbased_per_gene.tsv"),
                                  show_col_types = FALSE))
summ    <- as.data.frame(read_tsv(file.path(TAB_DIR, "16_mbased_summary.tsv"),
                                  show_col_types = FALSE))

NSIM  <- as.integer(env_or("MBASED_NSIM", "100000"))
FLOOR <- 1 / NSIM
message("simulations per gene: ", format(NSIM, big.mark = " "),
        "  ->  p-values resolved down to ", format(FLOOR, scientific = TRUE))

# Gene symbols where Ensembl has one, the Ensembl ID otherwise.
sym <- setNames(BIO$external_gene_name, BIO$ensembl_gene_id)
lab <- unname(sym[strip_version(ase_all$geneID)])
ase_all$label <- ifelse(is.na(lab) | !nzchar(lab), ase_all$geneID, lab)

ase_all$p_plot   <- pmax(ase_all$pASE,    FLOOR)
ase_all$fdr_plot <- pmax(ase_all$pASE_BH, FLOOR)
ase_all$at_floor <- ase_all$pASE <= FLOOR

# --- how much of the output sits at the resolution limit -------------------
floor_tab <- do.call(rbind, lapply(split(ase_all, ase_all$sample), function(d)
    data.frame(sample = d$sample[1], genes = nrow(d), at_floor = sum(d$at_floor),
               pct_at_floor = round(100 * mean(d$at_floor), 1))))
rownames(floor_tab) <- NULL
print(floor_tab, row.names = FALSE)
save_tab(floor_tab, "16_p_at_floor")
message("\nacross all samples: ", sprintf("%.1f%%", 100 * mean(ase_all$at_floor)),
        " of gene-level p-values are below the simulation's resolution of ",
        format(FLOOR, scientific = TRUE))

rep_sample <- summ$sample[which.max(summ$genes)]
df <- ase_all[ase_all$sample == rep_sample, ]
message("representative sample (most genes tested): ", rep_sample)

# --- volcano ----------------------------------------------------------------
top_lab <- df[df$sig, ]
top_lab <- top_lab[order(-top_lab$imbalance), ][seq_len(min(10, nrow(top_lab))), ]

save_fig(
    ggplot(df, aes(imbalance, -log10(fdr_plot), colour = sig)) +
        geom_point(alpha = 0.7, size = 1.6) +
        geom_hline(yintercept = -log10(FDR_CUT), linetype = "dashed") +
        geom_hline(yintercept = -log10(FLOOR), linetype = "dotted", colour = "grey40") +
        scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick"),
                            labels = c("no", "yes")) +
        geom_text_repel(data = top_lab, aes(label = label), size = 2.8,
                        max.overlaps = 30, show.legend = FALSE) +
        labs(title = paste0("ASE (MBASED) - ", rep_sample),
             subtitle = paste0(sprintf("%.0f%%", 100 * mean(df$at_floor)),
                               " of genes sit on the dotted line: p cannot be resolved below 1/numSim = ",
                               format(FLOOR, scientific = TRUE)),
             x = "allelic imbalance |MAF - 0.5| x 2",
             y = expression(-log[10]("FDR")),
             colour = paste0("FDR < ", FDR_CUT)) +
        theme_bw(),
    "16_ase_volcano", width = 8, height = 6
)

# --- imbalance across all samples -------------------------------------------
save_fig(
    ggplot(ase_all, aes(imbalance, fill = sig)) +
        geom_histogram(bins = 40, alpha = 0.85) +
        facet_wrap(~ sample, scales = "free_y") +
        scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick"),
                          labels = c("no", "yes")) +
        labs(title = "Allelic imbalance across all samples",
             x = "imbalance (0 = balanced, 1 = monoallelic)",
             y = "genes", fill = paste0("FDR < ", FDR_CUT)) +
        theme_bw(),
    "16_ase_imbalance_hist", width = 12, height = 9
)

# --- QQ plot ------------------------------------------------------------------
# Under the null this follows the diagonal. Here the whole distribution is lifted
# off it and a large share is pinned to the resolution limit, which is a statement
# about the model's assumptions failing on this input.
p  <- sort(df$p_plot)
qq <- data.frame(expected = -log10(ppoints(length(p))), observed = -log10(p))
save_fig(
    ggplot(qq, aes(expected, observed)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        geom_hline(yintercept = -log10(FLOOR), linetype = "dotted", colour = "grey40") +
        geom_point(alpha = 0.6, size = 1) +
        labs(title = paste0("QQ plot of ASE p-values - ", rep_sample),
             subtitle = paste0("dotted line = resolution limit, ",
                               format(FLOOR, scientific = TRUE),
                               "; ", sprintf("%.0f%%", 100 * mean(df$at_floor)),
                               " of genes are on it"),
             x = expression("expected " * -log[10]*"(p)"),
             y = expression("observed " * -log[10]*"(p)")) +
        theme_bw(),
    "16_ase_qq", width = 6, height = 6
)

banner("16b done")
