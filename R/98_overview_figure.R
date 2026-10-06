#!/usr/bin/env Rscript
#
# 98_overview_figure.R — one figure that carries the argument of the whole project.
#
# Four panels, each drawn from a table an earlier step wrote, so every value on
# the figure can be traced to its source and none is typed in:
#
#   A  the libraries differ before any biology      (x2_composition)
#   B  hits appear only against the controls        (11_deg_counts, 17_splicing_summary)
#   C  apparent ASE tracks doubtful genotype calls  (x4_per_sample)
#   D  power is not established                     (18_power_by_alpha)
#
# This script computes nothing new. It is a presentation of results that exist.
#
suppressPackageStartupMessages({
    library(ggpubr)
})
source(file.path(Sys.getenv("R_DIR", unset = "R"), "00_common.R"))

banner("98 - overview figure")

read_table <- function(name) {
    f <- file.path(TAB_DIR, paste0(name, ".tsv"))
    if (!file.exists(f)) stop("missing ", f, " - run the step that writes it first")
    as.data.frame(read_tsv(f, show_col_types = FALSE))
}

comp <- read_table("x2_composition")
deg  <- read_table("11_deg_counts")
spl  <- read_table("17_splicing_summary")
xp   <- read_table("x4_per_sample")
pal  <- read_table("18_power_by_alpha")

BASE <- 11
theme_overview <- function() {
    theme_bw(base_size = BASE) +
        theme(plot.title = element_text(face = "bold", size = BASE + 1),
              plot.subtitle = element_text(colour = "grey30"),
              legend.position = "bottom", legend.title = element_blank(),
              panel.grid.minor = element_blank())
}

# --- A: composition ----------------------------------------------------------
sp_all <- spearman_test(comp$rp_fraction, comp$intronic_pct)
pA <- ggplot(comp, aes(intronic_pct, rp_fraction, colour = group)) +
    geom_point(size = 3.2) +
    scale_colour_manual(values = GROUP_COLOURS) +
    labs(title = "A  The libraries differ before any biology",
         subtitle = sprintf("Spearman rho = %.2f across the 12 samples", sp_all$rho),
         x = "intronic reads (%)", y = "reads on ribosomal protein genes (%)") +
    theme_overview()

# --- B: where the hits are ---------------------------------------------------
ORDER  <- c("acute_vs_control", "subacute_vs_control", "acute_vs_subacute")
LABELS <- c("acute vs\ncontrol", "subacute vs\ncontrol", "acute vs\nsubacute")

b_genes <- data.frame(contrast = deg$contrast[deg$tool == "DESeq2"],
                      n = deg$fdr[deg$tool == "DESeq2"], layer = "genes (DESeq2)")
b_spl   <- data.frame(contrast = spl$contrast, n = spl$significant, layer = "splicing events")
bb <- rbind(b_genes, b_spl)
bb$contrast <- factor(bb$contrast, levels = ORDER, labels = LABELS)

pB <- ggplot(bb, aes(contrast, n, fill = layer)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = formatC(n, big.mark = ",", format = "d")),
              position = position_dodge(width = 0.8), vjust = -0.4, size = 3.3) +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = "B  Hits appear against controls, not within patients",
         subtitle = "|log2 FC| > 1 and FDR < 0.05",
         x = NULL, y = "significant genes or events") +
    theme_overview()

# --- C: ASE against genotype quality -----------------------------------------
sp_ase <- spearman_test(xp$pct_significant, xp$pct_below_0.4)
pC <- ggplot(xp, aes(pct_below_0.4, pct_significant, colour = group)) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "grey40",
                linewidth = 0.5, formula = y ~ x) +
    geom_point(size = 3.2) +
    scale_colour_manual(values = GROUP_COLOURS) +
    labs(title = "C  Apparent ASE tracks doubtful genotype calls",
         subtitle = sprintf("Spearman rho = %.2f across the 12 samples", sp_ase$rho),
         x = "heterozygous calls with alternate fraction < 0.4 (%)",
         y = "genes called imbalanced (%)") +
    theme_overview()

# --- D: the power bracket ----------------------------------------------------
pal$label <- ifelse(pal$alpha >= 0.01, "no correction\n(upper bound)", "Bonferroni\n(lower bound)")
pal$label <- factor(pal$label, levels = c("Bonferroni\n(lower bound)", "no correction\n(upper bound)"))
pD <- ggplot(pal, aes(label, power_n4)) +
    geom_col(fill = "grey55", width = 0.5) +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "firebrick") +
    annotate("text", x = 0.55, y = 0.84, label = "conventional 0.8", hjust = 0,
             colour = "firebrick", size = 3.3) +
    geom_text(aes(label = sprintf("%.2f", power_n4)), vjust = -0.4, size = 3.6) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.04))) +
    labs(title = "D  Power within patients is not established",
         subtitle = "twofold change, median gene, 4 per group, edgeR variability estimate",
         x = NULL, y = "power") +
    theme_overview()

fig <- ggarrange(pA, pB, pC, pD, ncol = 2, nrow = 2)
save_fig(fig, "98_overview", width = 13, height = 10)

banner("98 done")
