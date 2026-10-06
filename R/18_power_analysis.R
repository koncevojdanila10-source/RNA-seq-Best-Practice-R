#!/usr/bin/env Rscript
#
# 18_power_analysis.R — statistical power, parameterised from this dataset.
#
# The question a reviewer asks about any negative RNA-seq result: was there no
# effect, or was the experiment unable to see one? Power is the probability of
# detecting a difference that genuinely exists, and 0.8 is the conventional
# floor — out of 100 real effects, 80 are found.
#
# Four levers move it:
#
#   depth        helps weakly expressed genes, then flattens out
#   effect size  twofold is easy, 20% is disproportionately harder
#   CV           within-group variability; the strongest and least tractable
#   N            replicates per group; the biggest return per unit of effort
#
# The course demonstrated this on invented numbers (cv = 0.4, "typical for
# mice"). The lecturer then said the right thing to do is measure CV from real
# data of the same kind, so that is what happens here: the coefficient of
# variation is computed per gene within each group from our own normalised
# counts, and the curves are drawn at the values this dataset actually has.
#
# That turns the exercise into an answer to a specific question. The
# within-patient contrast — the only one not confounded by library composition
# — returned 2 genes. How many patients would it have taken?
#
suppressPackageStartupMessages({
    library(RNASeqPower)
    library(DESeq2)
})
source(file.path(Sys.getenv("R_DIR", unset = "R"), "00_common.R"))

banner("18 - power analysis on measured parameters")

load(file.path(RDATA_DIR, "11_de.RData"))

# ---------------------------------------------------------------------------
# 1. CV measured from our own counts
# ---------------------------------------------------------------------------
# CV is computed on normalised counts, within group, per gene. Pooling the
# groups would mix the between-group difference into the within-group spread
# and inflate it; that is the quantity the power calculation must not contain.
banner("coefficient of variation, measured per group")

norm_counts <- counts(ddsHTSeq, normalized = TRUE)
groups <- levels(sampleTable$CONDITION)

cv_of <- function(v) {
    m <- mean(v)
    if (!is.finite(m) || m <= 0) return(NA_real_)
    sd(v) / m
}

cv_tab <- list()
for (g in groups) {
    cols <- rownames(sampleTable)[sampleTable$CONDITION == g]
    sub  <- norm_counts[, cols, drop = FALSE]
    # Only genes with usable expression in this group: a CV estimated on a
    # handful of reads is noise about noise.
    keep <- rowMeans(sub) >= MIN_COUNT
    cvs  <- apply(sub[keep, , drop = FALSE], 1, cv_of)
    cvs  <- cvs[is.finite(cvs)]
    cv_tab[[g]] <- data.frame(
        group    = g,
        samples  = length(cols),
        genes    = length(cvs),
        cv_median = round(median(cvs), 3),
        cv_q25   = round(quantile(cvs, 0.25), 3),
        cv_q75   = round(quantile(cvs, 0.75), 3),
        cv_q90   = round(quantile(cvs, 0.90), 3)
    )
}
cv_summary <- do.call(rbind, cv_tab)
rownames(cv_summary) <- NULL
print(cv_summary, row.names = FALSE)
save_tab(cv_summary, "18_measured_cv")

# The patient groups are the ones the uncontaminated contrast compares, so
# their CV is the relevant number for the question being asked.
patient_cv <- median(cv_summary$cv_median[cv_summary$group != "control"])
all_cv     <- median(cv_summary$cv_median)
message("\nmedian within-group CV, patients only: ", round(patient_cv, 3))
message("median within-group CV, all groups:    ", round(all_cv, 3))
message("for reference, the course used 0.4 as a value typical of mice; ",
        "human cohorts are usually quoted at 0.8 and above")

# --- cross-check against edgeR's own dispersion ----------------------------
# The number above is a naive one: the standard deviation of normalised counts
# over the mean, per gene, median over genes. rnapower() wants the biological
# coefficient of variation, and edgeR already estimated it properly in step 11
# by fitting the negative binomial. BCV is the square root of the dispersion.
#
# A CV of 0.24 is low for a human cohort, and the whole conclusion of this
# script rests on it, so it does not get used without a second opinion.
banner("cross-check: edgeR's biological coefficient of variation")

bcv_common  <- sqrt(y$common.dispersion)
bcv_trended <- median(sqrt(y$trended.dispersion), na.rm = TRUE)
bcv_tagwise <- median(sqrt(y$tagwise.dispersion), na.rm = TRUE)

message("naive CV of normalised counts (median gene): ", round(patient_cv, 3))
message("edgeR common BCV:                            ", round(bcv_common, 3))
message("edgeR trended BCV (median gene):             ", round(bcv_trended, 3))
message("edgeR tagwise BCV (median gene):             ", round(bcv_tagwise, 3))

bcv_tab <- data.frame(
    estimator = c("naive CV of normalised counts", "edgeR common",
                  "edgeR trended (median)", "edgeR tagwise (median)"),
    value = round(c(patient_cv, bcv_common, bcv_trended, bcv_tagwise), 4)
)
save_tab(bcv_tab, "18_cv_estimators")

# The naive CV contains technical Poisson noise on top of the biological
# signal, so it should sit at or above the model-based BCV. If it sits well
# below, something is wrong with how it was computed and it must not be the
# number the curves are drawn at.
if (patient_cv < bcv_common * 0.8) {
    message("")
    message("The naive CV is well below edgeR's: it cannot be smaller than the ",
            "biological component it contains, so it is an underestimate.")
    message("Using edgeR's common BCV instead for everything below.")
    patient_cv <- bcv_common
} else {
    message("")
    message("The two agree to within 20%, so the measured CV stands.")
}

# Power is a per-gene quantity and the median gene is the easy case. The genes
# that matter are often the variable ones, so the harder quantiles are carried
# through rather than quietly dropped.
cv_q75 <- median(cv_summary$cv_q75[cv_summary$group != "control"])
cv_q90 <- median(cv_summary$cv_q90[cv_summary$group != "control"])
message("CV at the 75th percentile of genes: ", round(cv_q75, 3))
message("CV at the 90th percentile of genes: ", round(cv_q90, 3))

# Median depth per gene, used as the depth argument. The lower end is taken
# rather than the mean: the mean is dragged upwards by a few very highly
# expressed genes, and power at the mean would describe genes that were never
# the problem.
depth_used <- round(median(rowMeans(norm_counts)[rowMeans(norm_counts) >= MIN_COUNT]))
message("median per-gene depth among expressed genes: ", depth_used, " counts")

# ---------------------------------------------------------------------------
# 2. The course's grid, at our CV
# ---------------------------------------------------------------------------
banner("power curves")

n_range  <- 2:30
fc_range <- c(1.5, 2.0, 3.0)
cv_range <- c(0.4, round(patient_cv, 2), 0.8)
cv_range <- sort(unique(cv_range))

grid <- expand.grid(N = n_range, FoldChange = fc_range, CV = cv_range)
grid$Power <- mapply(function(n, fc, cv) {
    rnapower(depth = depth_used, n = n, cv = cv, effect = fc, alpha = FDR_CUT)
}, grid$N, grid$FoldChange, grid$CV)

grid$FoldChange <- factor(grid$FoldChange)
grid$CV_label <- factor(grid$CV, labels = paste0("CV = ", sort(unique(grid$CV))))

save_fig(
    ggplot(grid, aes(N, Power, colour = FoldChange)) +
        geom_hline(yintercept = 0.8, linetype = "dashed", colour = "red") +
        geom_line(linewidth = 1) +
        geom_point(size = 1.6) +
        facet_wrap(~ CV_label) +
        scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
        scale_colour_brewer(palette = "Set1") +
        labs(title = "RNA-seq power against replicates, effect size and variability",
             subtitle = paste0("depth = ", depth_used,
                               " counts per gene; per-test alpha = ", FDR_CUT,
                               ", no multiple-testing correction (an upper bound)\n",
                               "CV ", round(patient_cv, 2), " is this dataset's edgeR BCV; ",
                               "0.4 is the course's value; 0.8 is typical of human cohorts"),
             x = "biological replicates per group",
             y = "probability of detecting a DEG",
             colour = "fold change") +
        theme_bw(base_size = 12),
    "18_power_curves", width = 13, height = 5
)
save_tab(grid[, c("N", "FoldChange", "CV", "Power")], "18_power_grid")

# ---------------------------------------------------------------------------
# 3. How many replicates would each case have needed?
# ---------------------------------------------------------------------------
banner("replicates required for 80% power")

needed <- function(fc, cv, target = 0.8) {
    for (n in 2:200) {
        p <- rnapower(depth = depth_used, n = n, cv = cv, effect = fc, alpha = FDR_CUT)
        if (is.finite(p) && p >= target) return(n)
    }
    NA_integer_
}

req <- expand.grid(FoldChange = fc_range, CV = cv_range)
req$n_needed <- mapply(needed, req$FoldChange, req$CV)
req$have     <- 4L   # four patients, four controls
req$enough   <- ifelse(is.na(req$n_needed), FALSE, req$n_needed <= req$have)
print(req, row.names = FALSE)
save_tab(req, "18_replicates_required")

# ---------------------------------------------------------------------------
# 4. The actual question
# ---------------------------------------------------------------------------
banner("the within-patient contrast")

n_have   <- sum(sampleTable$CONDITION == "acute")
n_needed <- needed(2.0, patient_cv)
observed <- length(deg_lists$acute_vs_subacute$DESeq2)

message("acute vs subacute is the one contrast free of the composition ",
        "confounder (see analyses/x2_confounder.R).")
message("  replicates available:        ", n_have)
message("  measured CV in those groups: ", round(patient_cv, 3))
message("  depth used:                  ", depth_used, " counts")
message("  DEGs actually found:         ", observed)
power_at_hand <- rnapower(depth = depth_used, n = n_have, cv = patient_cv,
                          effect = 2.0, alpha = FDR_CUT)
message("  power at a twofold effect:   ", round(power_at_hand, 3))
message("  replicates for 80% power:    ",
        if (is.na(n_needed)) "not reachable within 200" else n_needed)

# The median gene is the easy case. The same calculation at the CV of the
# 75th and 90th percentile genes shows how much of the transcriptome the
# headline number actually covers.
power_by_quantile <- data.frame(
    gene_quantile = c("median", "75th percentile", "90th percentile"),
    cv = round(c(patient_cv, cv_q75, cv_q90), 3)
)
power_by_quantile$power_n4 <- round(mapply(function(cv)
    rnapower(depth = depth_used, n = n_have, cv = cv, effect = 2.0, alpha = FDR_CUT),
    power_by_quantile$cv), 3)
power_by_quantile$n_for_80pct <- mapply(needed, 2.0, power_by_quantile$cv)
message("\npower at a twofold effect, by how variable the gene is:")
print(power_by_quantile, row.names = FALSE)
save_tab(power_by_quantile, "18_power_by_gene_variability")

# Sensitivity to multiple testing. rnapower's alpha is the significance level
# of ONE test, and the course passes 0.05 straight in. This dataset tests
# over a thousand genes, so the level actually in force is lower and the
# figure above is an upper bound. Benjamini-Hochberg sits between the two
# rows below: more permissive than Bonferroni, stricter than no correction.
# Pinning it down exactly needs an assumption about how many genes truly
# change, which this dataset cannot supply, so the bracket is reported instead.
n_tests    <- nrow(ddsHTSeq)
alpha_sens <- data.frame(
    assumption = c("no correction: alpha = 0.05 per test (as in the course)",
                   paste0("Bonferroni: alpha = 0.05 / ", n_tests)),
    alpha      = c(FDR_CUT, FDR_CUT / n_tests)
)
alpha_sens$power_n4 <- round(mapply(function(a)
    rnapower(depth = depth_used, n = n_have, cv = patient_cv,
             effect = 2.0, alpha = a), alpha_sens$alpha), 3)
alpha_sens$n_for_80pct <- mapply(function(a) {
    for (n in 2:200)
        if (rnapower(depth = depth_used, n = n, cv = patient_cv,
                     effect = 2.0, alpha = a) >= 0.8) return(n)
    NA_integer_
}, alpha_sens$alpha)
message("\npower at a twofold effect, median gene, by multiple-testing assumption:")
print(alpha_sens, row.names = FALSE)
save_tab(alpha_sens, "18_power_by_alpha")

message("")
# The verdict keys on the whole bracket, not on its generous end. Power with no
# multiple-testing correction is an upper bound, Bonferroni is a lower one, and
# a claim of "informative null" is only earned if the lower end clears 0.8.
power_hi <- alpha_sens$power_n4[1]
power_lo <- alpha_sens$power_n4[2]
message("Power at a twofold change, median gene, 4 replicates per group:")
message("  upper bound (no correction):  ", power_hi)
message("  lower bound (Bonferroni):     ", power_lo)

if (is.na(power_lo) || power_lo < 0.8) {
    message("")
    message("The lower bound is far below 0.8, and the real figure for an FDR ",
            "procedure lies somewhere between the two. The power of this")
    message("experiment to detect a twofold change in a typical gene is ",
            "therefore NOT ESTABLISHED, and it is poor for the more variable")
    message("genes whichever end is taken. The near-empty result for this ",
            "contrast cannot be read as evidence of no change over time, and")
    message("the 2 genes found cannot be read as a finding either. Both are ",
            "what an experiment this size can produce with or without an effect.")
} else {
    message("Even the Bonferroni lower bound clears 0.8, so a near-empty ",
            "result is informative for a twofold change at the median gene.")
}
message("This is power at a twofold change only; a smaller true change is not ",
        "excluded by any of these numbers.")

# One number worth having for anyone designing the follow-up: the depth is
# cheap to change and the cohort is not, so it is worth knowing whether depth
# would have helped at all.
banner("would more sequencing have helped?")
for (d in c(depth_used, depth_used * 2, depth_used * 10, depth_used * 100)) {
    p <- rnapower(depth = d, n = n_have, cv = patient_cv, effect = 2.0, alpha = FDR_CUT)
    message(sprintf("  depth %7.0f counts -> power %.3f", d, p))
}
message("Depth saturates: within-group variability and sample size are what ",
        "bind here, exactly as the lecture's four levers predict.")

save(cv_summary, grid, req, patient_cv, depth_used, power_at_hand,
     file = file.path(RDATA_DIR, "18_power.RData"))
banner("18 done")
