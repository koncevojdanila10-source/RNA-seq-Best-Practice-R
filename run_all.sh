#!/usr/bin/env bash
#
# run_all.sh — every analysis step in order, from the upstream pipeline's output
# to SUMMARY.md and showcase/.
#
#   bash run_all.sh                 everything
#   bash run_all.sh --from 14       resume at a named step
#   bash run_all.sh --help
#
# Steps are numbered in execution order. Each one reads what the earlier ones
# saved under $WORK_DIR/results/rdata, so a failed run is resumed with --from
# rather than started over. The slow stages cache their own results (the SGSeq
# objects, the gene annotation, every WebGestalt query), so re-running a step
# that already finished costs seconds, not hours.
#
# Rough cost on the reference machine, 8 cores:
#   first run of step 17 (SGSeq)       about 70 minutes
#   step 16 (MBASED, 12 samples)       about 30 minutes
#   steps 12 and x3 (WebGestalt)       about 10-20 minutes each
#   everything else                    a few minutes in total
#
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
source config.sh

STEPS=(11 12 13 14 15 16 16b 17 18 x1 x2 x3 x4 x6 98 99 99b)

usage() {
    sed -n '2,22p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    echo "steps: ${STEPS[*]}"
}

FROM=""
while [ $# -gt 0 ]; do
    case "$1" in
        --from)    FROM="${2:?--from needs a step name}"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *)         echo "unknown argument: $1" >&2; usage >&2; exit 1 ;;
    esac
done

if [ -n "${FROM}" ]; then
    ok=0; for s in "${STEPS[@]}"; do [ "$s" = "${FROM}" ] && ok=1; done
    [ "${ok}" = 1 ] || { echo "ERROR: unknown step '${FROM}'. Steps: ${STEPS[*]}" >&2; exit 1; }
fi

# --- the upstream pipeline's output must exist ------------------------------
# This project aligns nothing. Without that output there is nothing to analyse,
# and finding out in step 11 is no better than finding out here.
log "checking the upstream output in ${UPSTREAM_DIR}"
problem=0
n_counts=$(ls "${COUNTS_DIR}"/*.htseq.counts 2>/dev/null | wc -l)
n_bam=$(ls "${BAM_DIR}"/*.bam 2>/dev/null | wc -l)
n_ase=$(ls "${ASE_DIR}"/*.ase.tsv 2>/dev/null | wc -l)
expected=$(awk 'NR>1 && $0 !~ /^#/ {n++} END {print n}' "${SAMPLE_SHEET}")
for pair in "count files:${n_counts}" "BAM files:${n_bam}" "ASE tables:${n_ase}"; do
    what="${pair%%:*}"; have="${pair##*:}"
    if [ "${have}" -ne "${expected}" ]; then
        echo "ERROR: expected ${expected} ${what}, found ${have}" >&2; problem=1
    fi
done
[ -s "${ANNOTATION_GTF}" ] || { echo "ERROR: missing annotation ${ANNOTATION_GTF}" >&2; problem=1; }
if [ "${problem}" = 1 ]; then
    echo "" >&2
    echo "Run the upstream pipeline first:" >&2
    echo "  https://github.com/koncevojdanila10-source/RNA-seq-Best-Practice" >&2
    echo "or point RNA_UPSTREAM_DIR at its work directory." >&2
    exit 1
fi
log "upstream output found: ${n_counts} samples"

# --- environment ------------------------------------------------------------
if conda env list | awk '{print $1}' | grep -qx "${CONDA_ENV}"; then
    log "conda environment '${CONDA_ENV}' exists"
else
    log "creating the '${CONDA_ENV}' environment (this is the slow part)"
    bash setup_env.sh
fi

# --- the steps --------------------------------------------------------------
started=0; [ -z "${FROM}" ] && started=1
t_all=$(date +%s)
for s in "${STEPS[@]}"; do
    [ "$s" = "${FROM}" ] && started=1
    [ "${started}" = 1 ] || continue
    log "=== step ${s} ==="
    t0=$(date +%s)
    bash run_step.sh "${s}"
    log "step ${s} finished in $(( $(date +%s) - t0 )) s"
done

log "all steps finished in $(( ( $(date +%s) - t_all ) / 60 )) min"
log "results: SUMMARY.md and showcase/"
