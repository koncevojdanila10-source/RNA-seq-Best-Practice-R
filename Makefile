# Short names for the commands in this repository. Every target only calls a
# shell script; the logic lives there, so nothing here can drift from it.
#
#   make          list the targets
#   make all      every step, from the upstream output to SUMMARY.md
#   make summary  regenerate SUMMARY.md and showcase/ from the saved tables
#   make check    confirm README.md still matches the tables
#
.RECIPEPREFIX = >
.DEFAULT_GOAL := help
.PHONY: help all summary check overview steps

help:
> @sed -n '4,7p' Makefile | sed 's/^# \{0,1\}//'

all:
> bash run_all.sh

summary:
> bash summarize_results.sh

check:
> bash run_step.sh 99b

overview:
> bash run_step.sh 98

steps:
> bash run_step.sh --list
