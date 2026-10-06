#!/usr/bin/env bash
#
# setup_env.sh — one conda environment with everything days 4-6 need.
#
# Installing R and Bioconductor through conda is the part of this project most
# likely to fail on someone else's machine, so the script is defensive rather
# than short:
#
#   1. a CORE tier that every later script needs — if this fails, stop;
#   2. an EXTRA tier of the more specialised packages, installed in one batch,
#      falling back to one-at-a-time so a single unavailable build does not
#      take the whole batch down with it;
#   3. a final check inside R that installs anything still missing with
#      BiocManager (which serves CRAN packages too), and refuses to report
#      success while a required package is absent.
#
# The version manifest at the end is written from R's own package registry, so
# it records what actually got installed rather than what was asked for.
#
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

mkdir -p "${LOG_DIR}" "${PROJECT_DIR}/envs"
exec > >(tee -a "${LOG_DIR}/setup_env.log") 2>&1

log "=== setup_env.sh — environment '${CONDA_ENV}' ==="

# --- solver ----------------------------------------------------------------
# mamba resolves an environment this size in a fraction of the time. Use it
# when it is there, fall back to conda when it is not.
if command -v mamba >/dev/null 2>&1; then
    SOLVER=mamba
else
    SOLVER=conda
    log "NOTE: mamba not found, falling back to conda — expect 20-40 minutes"
fi
log "solver: ${SOLVER}"

CH=(-c conda-forge -c bioconda)

# --- tier 1: core ----------------------------------------------------------
# Differential expression, the transformations used to compare normalisations,
# annotation and the plotting stack. Nothing here is optional.
CORE=(
    r-biocmanager
    bioconductor-deseq2
    bioconductor-edger
    bioconductor-limma
    bioconductor-vsn
    bioconductor-biomart
    bioconductor-summarizedexperiment
    bioconductor-genomicranges
    bioconductor-genomeinfodb
    r-tidyverse
    r-data.table
    r-matrixstats
    r-rcolorbrewer
    r-ggpubr
    r-pheatmap
    r-ggrepel
    r-hexbin
    zip
)

# --- tier 2: the specialised packages --------------------------------------
# One per analysis of days 5-6. Any of these may be missing a build for the
# platform, which is why they are installed separately from the core.
EXTRA=(
    bioconductor-enhancedvolcano
    bioconductor-genomicfeatures
    bioconductor-txdbmaker
    bioconductor-biocparallel
    bioconductor-rtracklayer
    bioconductor-quantiseqr
    bioconductor-mbased
    bioconductor-sgseq
    bioconductor-bionero
    bioconductor-rnaseqpower
    r-wgcna
    r-ggvenndiagram
    r-webgestaltr
)

# Every R package the scripts actually library(), used for the final check.
# Kept separate from the conda names on purpose: the two namespaces do not
# match one-to-one and pretending they do is how silent gaps appear.
REQUIRED_R=(
    DESeq2 edgeR limma vsn biomaRt EnhancedVolcano
    SummarizedExperiment GenomicRanges GenomicFeatures GenomeInfoDb
    txdbmaker BiocParallel rtracklayer
    quantiseqr MBASED SGSeq BioNERO RNASeqPower WGCNA
    tidyverse data.table matrixStats RColorBrewer ggpubr pheatmap
    ggrepel ggVennDiagram WebGestaltR hexbin
)

# --- create ----------------------------------------------------------------
if conda env list | awk '{print $1}' | grep -qx "${CONDA_ENV}"; then
    log "environment '${CONDA_ENV}' already exists — installing into it"
else
    log "creating '${CONDA_ENV}' with the core tier (this is the slow part)"
    "${SOLVER}" create -y -n "${CONDA_ENV}" "${CH[@]}" "${CORE[@]}"
fi

# --- tier 2, batch then one-by-one -----------------------------------------
log "installing the specialised packages"
if ! "${SOLVER}" install -y -n "${CONDA_ENV}" "${CH[@]}" "${EXTRA[@]}"; then
    log "the batch failed — retrying one package at a time to isolate it"
    for pkg in "${EXTRA[@]}"; do
        if "${SOLVER}" install -y -n "${CONDA_ENV}" "${CH[@]}" "${pkg}"; then
            log "  ok      ${pkg}"
        else
            log "  MISSING ${pkg} — will try BiocManager instead"
        fi
    done
fi

# --- final check inside R --------------------------------------------------
# conda reporting success is not the same as library() working. Ask R.
log "checking every required package from inside R"
cat > "${TMPDIR:-/tmp}/check_pkgs.R" <<'REOF'
required <- commandArgs(trailingOnly = TRUE)
missing  <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) {
    message("missing after conda: ", paste(missing, collapse = ", "))
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
    # BiocManager serves CRAN as well, so one call covers both namespaces.
    BiocManager::install(missing, ask = FALSE, update = FALSE)
    missing <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
}

if (length(missing)) {
    message("STILL MISSING: ", paste(missing, collapse = ", "))
    quit(status = 1)
}
message("all ", length(required), " packages load")
REOF

conda run --no-capture-output -n "${CONDA_ENV}" \
    Rscript "${TMPDIR:-/tmp}/check_pkgs.R" "${REQUIRED_R[@]}"

# --- command-line tools the R packages shell out to ------------------------
# WebGestaltR calls `zip` to package its HTML report. It does so only after the
# enrichment has been computed, so a missing zip surfaces as an error thrown by
# a call that had already done the work — which is very easy to misread as
# "nothing was enriched". Check for it here instead.
if conda run -n "${CONDA_ENV}" bash -c 'command -v zip' >/dev/null 2>&1; then
    log "zip: present"
else
    log "WARNING: zip is not in the environment. WebGestaltR will run without"
    log "         writing HTML reports; the enrichment tables are unaffected."
fi

# --- version manifest ------------------------------------------------------
# Written from R's package registry, so it records what is installed rather
# than what was requested.
log "writing envs/versions_installed.txt"
cat > "${TMPDIR:-/tmp}/versions.R" <<'REOF'
pkgs <- commandArgs(trailingOnly = TRUE)
cat("# R and package versions actually installed\n")
cat("#", format(Sys.time(), "%Y-%m-%d"), "\n\n")
cat(R.version.string, "\n")
if (requireNamespace("BiocManager", quietly = TRUE))
    cat("Bioconductor", as.character(BiocManager::version()), "\n")
cat("\n")
for (p in sort(pkgs))
    cat(sprintf("%-22s %s\n", p, as.character(packageVersion(p))))
REOF

conda run --no-capture-output -n "${CONDA_ENV}" \
    Rscript "${TMPDIR:-/tmp}/versions.R" "${REQUIRED_R[@]}" \
    > "${PROJECT_DIR}/envs/versions_installed.txt"

log "=== done ==="
log "environment: ${CONDA_ENV}"
log "manifest:    envs/versions_installed.txt"
