#!/usr/bin/env bash
#
# summarize_results.sh — regenerate SUMMARY.md and showcase/ from the pipeline's
# output tables. Every number in SUMMARY.md is read or recomputed from those
# tables by R/99_summarize_results.R; none is typed in.
#
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
bash run_step.sh 99
