#!/usr/bin/env Rscript
#
# 15_bionero.R — the same co-expression question through BioNERO, and a formal
# comparison with the WGCNA result from step 14.
#
# BioNERO wraps WGCNA rather than replacing it: the same soft-thresholded
# correlation network is underneath. What it adds is a tidier interface, its
# own module merging, and hub-gene extraction. That makes it a genuine
# cross-check on the module structure — two implementations of one idea, run on
# identical input — rather than an independent method.
#
# The lecturer's main practical warning about it is followed here: the network
# type and correlation coefficient have to be chosen so the scale-free fit
# actually converges. On his four samples the signed and signed-hybrid networks
# did not converge and only unsigned did. Rather than assume ours behaves the
# same way, this script fits every combination and reports which ones reach the
# threshold before picking one.
#
# One recommendation of BioNERO's is deliberately not followed: it can remove
# the variance carried by the first principal components as a batch correction.
# On this dataset that would be actively misleading. Step x2 established that
# the leading axis of variation is library composition and that it is 93%
# explained by group, so stripping it would silently delete the group
# difference and leave a clean-looking network that had answered nothing.
#
suppressPackageStartupMessages({
    library(BioNERO)
    library(SummarizedExperiment)
    library(DESeq2)
})
source(file.path(Sys.getenv("R_DIR", unset = "R"), "00_common.R"))

banner("15 - BioNERO co-expression network")

load(file.path(RDATA_DIR, "11_de.RData"))
load(file.path(RDATA_DIR, "12_annot.RData"))
load(file.path(RDATA_DIR, "14_wgcna.RData"))

se_data <- SummarizedExperiment(
    assays  = list(counts = as.matrix(ntd_rename)),
    colData = sampleTable[colnames(ntd_rename), ]
)
message("input: ", nrow(se_data), " genes x ", ncol(se_data), " samples")

# The course filters to the 5000 most variable genes. We have 1109 in total, so
# that filter would be a no-op; keeping every gene also keeps the comparison
# with step 14 exact, which matters more here than speed.
message("no variance filter applied: the whole matrix is smaller than the ",
        "course's 5000-gene cut, and step 14 used all of it")

# ---------------------------------------------------------------------------
# 1. Which network type and correlation actually converge?
# ---------------------------------------------------------------------------
banner("scale-free fit across network types and correlation methods")

COMBOS <- expand.grid(
    net_type   = c("unsigned", "signed", "signed hybrid"),
    cor_method = c("pearson", "spearman"),
    stringsAsFactors = FALSE
)

fits <- list()
for (i in seq_len(nrow(COMBOS))) {
    nt <- COMBOS$net_type[i]; cm <- COMBOS$cor_method[i]
    r <- tryCatch({
        sft <- suppressWarnings(invisible(capture.output(
            out <- SFT_fit(se_data, net_type = nt, cor_method = cm)
        )))
        data.frame(net_type = nt, cor_method = cm,
                   power = out$power,
                   converged = !is.na(out$power))
    }, error = function(e) {
        message("  ", nt, " / ", cm, ": failed - ", conditionMessage(e))
        data.frame(net_type = nt, cor_method = cm,
                   power = NA_integer_, converged = FALSE)
    })
    message("  ", nt, " / ", cm, ": ",
            if (isTRUE(r$converged)) paste0("converged at power ", r$power)
            else "did not converge")
    fits[[i]] <- r
}
sft_table <- do.call(rbind, fits)
print(sft_table, row.names = FALSE)
save_tab(sft_table, "15_sft_by_network_type")

good <- sft_table[sft_table$converged, ]
if (!nrow(good))
    stop("no combination of network type and correlation converged - ",
         "there is nothing to build a network from")

# Two networks get built, because one comparison cannot answer the question on
# its own.
#
#   matched  — BioNERO forced onto step 14's settings: unsigned, Pearson (the
#              WGCNA default, which step 14 did not override), power 5. Any
#              disagreement here is attributable to the implementation.
#   default  — each package left to its own choice: BioNERO on the combination
#              it selected for itself. This is what someone following the two
#              manuals would actually get.
#
# Reporting only the second would blame the tools for a difference in settings;
# reporting only the first would hide how much the defaults matter.
pick <- good[good$net_type == "unsigned" & good$cor_method == "spearman", ]
if (!nrow(pick)) pick <- good[1, ]
NET_TYPE <- pick$net_type[1]; COR_METHOD <- pick$cor_method[1]
SFT_POWER_BN <- pick$power[1]

# What step 14 actually used: blockwiseModules defaults to Pearson and was not
# overridden there, at the power its own scale-free curve selected.
WGCNA_COR   <- "pearson"
WGCNA_POWER <- SFT_POWER

message("BioNERO's own choice:  ", NET_TYPE, " / ", COR_METHOD,
        " at power ", SFT_POWER_BN)
message("step 14 (WGCNA) used:  unsigned / ", WGCNA_COR,
        " at power ", WGCNA_POWER)

# ---------------------------------------------------------------------------
# 2. Build both
# ---------------------------------------------------------------------------
banner("building the BioNERO networks")

build_bn <- function(net_type, cor_method, power, label) {
    message("  ", label, ": ", net_type, " / ", cor_method, " power ", power)
    exp2gcn(
        se_data,
        net_type                 = net_type,
        SFTpower                 = power,
        cor_method               = cor_method,
        module_merging_threshold = 0.8
    )
}

net_bn      <- build_bn(NET_TYPE, COR_METHOD, SFT_POWER_BN, "BioNERO default")
net_matched <- build_bn("unsigned", WGCNA_COR, WGCNA_POWER, "matched to step 14")

modules_of <- function(net) {
    m <- net$genes_and_modules
    names(m) <- c("gene", "module")
    m
}
bn_modules      <- modules_of(net_bn)
matched_modules <- modules_of(net_matched)

bn_sizes <- as.data.frame(sort(table(bn_modules$module), decreasing = TRUE))
names(bn_sizes) <- c("module", "genes")
message("BioNERO modules (own settings): ", nrow(bn_sizes))
print(bn_sizes, row.names = FALSE)
save_tab(bn_sizes, "15_module_sizes")
save_tab(bn_modules, "15_gene_modules")
save_tab(matched_modules, "15_gene_modules_matched")

png(file.path(FIG_DIR, "15_dendrogram.png"), width = 2000, height = 1200, res = 200)
plot_dendro_and_colors(net_bn)
invisible(dev.off())

# ---------------------------------------------------------------------------
# 3. Hub genes
# ---------------------------------------------------------------------------
banner("hub genes")

hubs <- tryCatch(get_hubs_gcn(se_data, net_bn), error = function(e) {
    message("get_hubs_gcn failed: ", conditionMessage(e)); NULL
})
if (!is.null(hubs) && nrow(hubs)) {
    # Column names have moved between BioNERO versions, so find the module
    # column rather than assuming it is called Module.
    mod_col <- grep("module", names(hubs), ignore.case = TRUE, value = TRUE)[1]
    n_mod <- if (!is.na(mod_col)) length(unique(hubs[[mod_col]])) else NA_integer_
    message("hub genes: ", nrow(hubs), " across ", n_mod, " modules",
            "  (columns: ", paste(names(hubs), collapse = ", "), ")")
    print(head(hubs, 15), row.names = FALSE)
    save_tab(hubs, "15_hub_genes")

    # chr19 carries the largest cluster of KRAB zinc-finger genes in the human
    # genome, and they are strongly co-regulated. If the hubs are mostly ZNFs,
    # the network's strongest structure is that cluster — worth quantifying,
    # because it is a property of the chromosome we chose rather than of
    # stroke, and it would not appear this way in a whole-genome analysis.
    gene_col <- grep("gene", names(hubs), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(gene_col)) {
        znf <- grepl("^ZNF|^ZFP|^ZSCAN|^ZKSCAN", hubs[[gene_col]])
        message("hub genes from the zinc-finger families: ", sum(znf), " of ",
                nrow(hubs), sprintf(" (%.0f%%)", 100 * mean(znf)))
        save_tab(data.frame(hub_genes = nrow(hubs), zinc_finger = sum(znf),
                            pct = round(100 * mean(znf), 1)),
                 "15_hub_zinc_finger_share")
    }
} else {
    message("no hub genes returned")
}

# ---------------------------------------------------------------------------
# 4. Do the two implementations agree on the modules?
# ---------------------------------------------------------------------------
# This is the point of running both. Two partitions of the same 1109 genes are
# compared with the adjusted Rand index: 1 means identical grouping, 0 means no
# more agreement than chance, and negative means less than chance. It is
# written out here rather than pulled from a package so the arithmetic is
# visible and adds no dependency.
banner("WGCNA against BioNERO")

adjusted_rand <- function(a, b) {
    tab <- table(a, b)
    n   <- sum(tab)
    choose2 <- function(x) x * (x - 1) / 2
    sum_ij <- sum(choose2(tab))
    sum_a  <- sum(choose2(rowSums(tab)))
    sum_b  <- sum(choose2(colSums(tab)))
    expected <- sum_a * sum_b / choose2(n)
    maximum  <- (sum_a + sum_b) / 2
    (sum_ij - expected) / (maximum - expected)
}

wgcna_modules <- setNames(moduleColors, colnames(a_ntd))
common_genes  <- intersect(names(wgcna_modules), bn_modules$gene)
message("genes in both partitions: ", length(common_genes))
w <- wgcna_modules[common_genes]

compare_to_wgcna <- function(mods, label) {
    b <- setNames(as.character(mods$module), mods$gene)[common_genes]
    ari <- adjusted_rand(w, b)
    message("\n", label)
    message("  BioNERO modules: ", length(unique(b)),
            "   WGCNA modules: ", length(unique(w)))
    message("  adjusted Rand index: ", round(ari, 3))
    print(as.data.frame.matrix(table(WGCNA = w, BioNERO = b)))
    list(ari = ari, b = b, n_mod = length(unique(b)))
}

cmp_matched <- compare_to_wgcna(matched_modules,
    "matched settings (unsigned / pearson / power 5) - implementation only")
cmp_default <- compare_to_wgcna(bn_modules,
    "each package's own settings - what the manuals would give you")

save_tab(tibble::rownames_to_column(
    as.data.frame.matrix(table(WGCNA = w, BioNERO = cmp_default$b)), "wgcna_module"),
    "15_module_crosstab")

agreement <- data.frame(
    comparison      = c("matched settings", "each package's default"),
    bionero_cor     = c(WGCNA_COR, COR_METHOD),
    bionero_power   = c(WGCNA_POWER, SFT_POWER_BN),
    bionero_modules = c(cmp_matched$n_mod, cmp_default$n_mod),
    wgcna_modules   = length(unique(w)),
    genes_compared  = length(common_genes),
    adjusted_rand   = round(c(cmp_matched$ari, cmp_default$ari), 4)
)
print(agreement, row.names = FALSE)
save_tab(agreement, "15_partition_agreement")

ari <- cmp_matched$ari
banner("reading the two numbers")
message("Implementation alone (matched settings): ARI = ", round(cmp_matched$ari, 3))
message("Settings and implementation together:   ARI = ", round(cmp_default$ari, 3))
if (cmp_matched$ari > 0.7) {
    message("-> The two implementations agree once the settings are the same. ",
            "The disagreement at defaults is caused by the defaults, not the code.")
} else if (cmp_matched$ari > 0.3) {
    message("-> Only partial agreement even at identical settings: the same ",
            "broad structure, cut in different places.")
} else {
    message("-> The partitions disagree even at identical settings. With 12 ",
            "samples the correlation matrix is estimated too imprecisely for ",
            "module boundaries to be reproducible, so individual module ",
            "membership here should not be treated as a finding.")
}

# ---------------------------------------------------------------------------
# 5. Does the step 14 finding survive the unstable partition?
# ---------------------------------------------------------------------------
# An unstable partition and an unreliable finding are different things. Step 14
# reported that the largest module's eigengene tracks intronic fraction at
# rho = -0.98. If BioNERO's modules — cut differently, from a different
# correlation — reproduce an association of that size, then the finding is
# robust even though module membership is not. If they do not, both fall.
banner("does the composition association reproduce?")

bn_MEs <- tryCatch({
    me <- net_bn$MEs
    if (is.null(me)) NULL else as.data.frame(me)
}, error = function(e) NULL)

if (is.null(bn_MEs)) {
    # exp2gcn did not hand back eigengenes: compute them from the partition.
    cor <- WGCNA::cor
    bn_MEs <- WGCNA::moduleEigengenes(
        t(as.matrix(ntd_rename)),
        colors = setNames(as.character(bn_modules$module),
                          bn_modules$gene)[rownames(ntd_rename)]
    )$eigengenes
    cor <- stats::cor
    rownames(bn_MEs) <- colnames(ntd_rename)
}

common_samples <- intersect(rownames(bn_MEs), rownames(traits))
if (length(common_samples) < 3) {
    message("could not align BioNERO eigengenes to the sample traits - skipped")
} else {
    bn_cor <- cor(bn_MEs[common_samples, , drop = FALSE],
                  traits[common_samples, , drop = FALSE],
                  use = "p", method = "spearman")
    bn_p   <- WGCNA::corPvalueStudent(bn_cor, length(common_samples))
    bn_q   <- matrix(p.adjust(bn_p, method = "BH"), nrow = nrow(bn_p),
                     dimnames = dimnames(bn_p))

    bn_long <- data.frame(
        module = rep(rownames(bn_cor), ncol(bn_cor)),
        trait  = rep(colnames(bn_cor), each = nrow(bn_cor)),
        rho    = as.vector(bn_cor),
        q_BH   = as.vector(bn_q)
    )
    bn_long <- bn_long[order(-abs(bn_long$rho)), ]
    print(head(bn_long, 8), row.names = FALSE, digits = 3)
    save_tab(bn_long, "15_module_trait_bionero")

    best_bn <- max(abs(bn_cor[, "intronic_pct"]))
    best_wg <- max(abs(mtc[, "intronic_pct"]))
    message("\nstrongest |rho| with intronic fraction:")
    message("  WGCNA:   ", round(best_wg, 3))
    message("  BioNERO: ", round(best_bn, 3))
    if (best_bn > 0.8 && best_wg > 0.8)
        message("-> Both find a module tracking library composition almost ",
                "perfectly. The partition is unstable; this association is not.")
    else
        message("-> The association does not reproduce at comparable strength, ",
                "so it cannot be separated from the instability of the modules.")
}

save(net_bn, net_matched, bn_modules, matched_modules, hubs, sft_table,
     agreement, cmp_matched, cmp_default, ari,
     NET_TYPE, COR_METHOD, SFT_POWER_BN, WGCNA_COR, WGCNA_POWER,
     file = file.path(RDATA_DIR, "15_bionero.RData"))
banner("15 done")
