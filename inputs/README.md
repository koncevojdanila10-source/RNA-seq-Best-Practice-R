# inputs/

Small facts produced by the upstream pipeline that this project depends on but
does not recompute. Everything here is committed so the analyses in
`analyses/` can be reproduced from a clone of this repository alone, without
re-running twelve alignments.

| file | produced by | what it is |
|---|---|---|
| `read_distribution.tsv` | [RNA-seq-Best-Practice](https://github.com/koncevojdanila10-source/RNA-seq-Best-Practice) `scripts/07_qc.sh`, via RSeQC `read_distribution.py` | percentage of each library's reads falling in exonic, intronic and flanking regions |
| `relatedness.tsv` | same repository, `scripts/06_relatedness.sh` | pairwise genotype concordance between samples, used to check that each D1/D2 pair really is one person |

The large inputs — count files, BAMs and the ASEReadCounter tables — are not
here and are not committed. `config.sh` points at the upstream pipeline's work
directory for those.
