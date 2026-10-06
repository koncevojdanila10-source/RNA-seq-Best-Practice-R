#!/usr/bin/env bash
#
# run_step.sh — run one R script inside the project environment.
#
#   bash run_step.sh 11              # R/11_de_three_tools.R
#   bash run_step.sh x1              # analyses/x1_*.R
#   bash run_step.sh --list          # what is available
#
# Everything config.sh defines is exported before R starts, so the scripts read
# their paths from the environment and config.sh stays the single source of
# truth. Running a script by hand in RStudio still works: R/00_common.R falls
# back to the same defaults when the variables are not set.
#
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

export PROJECT_DIR R_DIR ANALYSES_DIR SHOWCASE_DIR LOG_DIR SAMPLE_SHEET
export UPSTREAM_DIR COUNTS_DIR BAM_DIR ASE_DIR REF_DIR ANNOTATION_GTF
export WORK_DIR RESULTS_DIR FIG_DIR TAB_DIR RDATA_DIR
export CHROM THREADS MBASED_NSIM

mkdir -p "${LOG_DIR}" "${RESULTS_DIR}"

list_scripts() {
    echo "pipeline steps (R/):"
    ls -1 "${R_DIR}"/[0-9]*.R 2>/dev/null | xargs -r -n1 basename | sed 's/^/  /'
    echo "extra analyses (analyses/):"
    ls -1 "${ANALYSES_DIR}"/*.R 2>/dev/null | xargs -r -n1 basename | sed 's/^/  /'
}

[ $# -eq 1 ] || { echo "usage: bash run_step.sh <step|--list>" >&2; list_scripts; exit 1; }
# Written as a full if: under `set -e` a bare `[ ... ] && { ... }` that tests
# false returns non-zero and takes the whole script down with it.
if [ "$1" = "--list" ]; then list_scripts; exit 0; fi

STEP="$1"
# Accept a bare number or prefix: 11, x1, or the full file name.
SCRIPT="$(ls -1 "${R_DIR}/${STEP}"_*.R "${ANALYSES_DIR}/${STEP}"_*.R \
          "${R_DIR}/${STEP}" "${ANALYSES_DIR}/${STEP}" 2>/dev/null | head -1 || true)"
[ -n "${SCRIPT}" ] || { echo "ERROR: no script matches '${STEP}'" >&2; list_scripts; exit 1; }

LOGFILE="${LOG_DIR}/${STEP}.log"
log "running ${SCRIPT}"
log "log -> ${LOGFILE}"

# --no-capture-output so R's progress reaches the terminal as it happens
# instead of arriving in one block at the end.
conda run --no-capture-output -n "${CONDA_ENV}" \
    Rscript "${SCRIPT}" 2>&1 | tee "${LOGFILE}"

log "done: ${SCRIPT}"
