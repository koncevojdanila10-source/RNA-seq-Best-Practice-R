#!/usr/bin/env bash
#
# config.sh — single source of truth for paths and resources.
# Sourced by setup_env.sh, run_all.sh and every script in R/.
# Nothing here runs analysis; it only defines variables.
#

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_DIR="${PROJECT_DIR}/R"
ANALYSES_DIR="${PROJECT_DIR}/analyses"
SHOWCASE_DIR="${PROJECT_DIR}/showcase"
LOG_DIR="${PROJECT_DIR}/logs"
SAMPLE_SHEET="${PROJECT_DIR}/samples.tsv"

# --- upstream ---------------------------------------------------------------
# This project does not align anything. It consumes the output of the upstream
# pipeline (STAR -> htseq-count -> GATK -> ASEReadCounter):
#
#   https://github.com/koncevojdanila10-source/RNA-seq-Best-Practice
#
# Override if that pipeline's work directory is somewhere else:
#   export RNA_UPSTREAM_DIR=/mnt/d/rna_work && bash run_all.sh
UPSTREAM_DIR="${RNA_UPSTREAM_DIR:-${HOME}/rna_work}"
COUNTS_DIR="${UPSTREAM_DIR}/results/counts"      # <sample>.htseq.counts
BAM_DIR="${UPSTREAM_DIR}/results/bam"            # <sample>.bam + .bai, STAR
ASE_DIR="${UPSTREAM_DIR}/results/ase"            # <sample>.ase.tsv, ase_tested.tsv
REF_DIR="${UPSTREAM_DIR}/data/reference"
ANNOTATION_GTF="${REF_DIR}/chr19.gtf"

# --- our own outputs --------------------------------------------------------
# Kept outside the repository: RData environments and TOM matrices are large
# and every one of them is regenerable from these scripts.
WORK_DIR="${RNA_WORK_DIR:-${HOME}/rna_work_r}"
RESULTS_DIR="${WORK_DIR}/results"
FIG_DIR="${RESULTS_DIR}/figures"
TAB_DIR="${RESULTS_DIR}/tables"
RDATA_DIR="${RESULTS_DIR}/rdata"

# --- scope ------------------------------------------------------------------
CHROM="chr19"
CONDA_ENV="rna_r"

# --- compute budget ---------------------------------------------------------
# Measured on the development machine: 8 cores / 27 GB RAM under WSL2.
THREADS=5              # WGCNA / SGSeq worker threads
MBASED_NSIM=100000    # simulations per gene. The course used 10^6 for a final
                      # run on one sample; with twelve samples and thousands of
                      # genes each that is hours, and BH correction dominates
                      # the smallest p-values anyway. Raise it for a final pass.

log() { echo "[$(date +%H:%M:%S)] $*"; }
