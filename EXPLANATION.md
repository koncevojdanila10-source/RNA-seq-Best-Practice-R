# How the scripts work

What each script does, which commands carry the logic, what their arguments mean,
and why each choice was made. Read [`README.md`](README.md) first for what was
found; this file is about how.

## Contents

1. [How the pieces fit together](#1-how-the-pieces-fit-together)
2. [Shell scripts](#2-shell-scripts)
3. [Shared R code: `00_common.R`](#3-shared-r-code-00_commonr)
4. [Step 11 — differential expression](#4-step-11--differential-expression)
5. [Step 12 — annotation and enrichment](#5-step-12--annotation-and-enrichment)
6. [Step 13 — deconvolution](#6-step-13--deconvolution)
7. [Step 14 — WGCNA](#7-step-14--wgcna)
8. [Step 15 — BioNERO](#8-step-15--bionero)
9. [Step 16 — allele-specific expression](#9-step-16--allele-specific-expression)
10. [Step 17 — splicing](#10-step-17--splicing)
11. [Step 18 — power](#11-step-18--power)
12. [The cross-checks, x1 to x6](#12-the-cross-checks-x1-to-x6)
13. [Steps 98, 99 and 99b — figure, summary and README check](#13-steps-98-99-and-99b--figure-summary-and-readme-check)

## 1. How the pieces fit together

```
config.sh  ──export──►  run_step.sh  ──►  Rscript R/NN_*.R  ──►  tables, figures, .RData
 (paths)               (environment)       sources 00_common.R         in ~/rna_work_r
                                                                           │
                                         99_summarize_results.R  ◄────────┘
                                                  │
                                       SUMMARY.md + showcase/
```

Three layers, each with one job:

- **`config.sh`** is the only place paths and resource limits are written down.
- **`run_step.sh`** reads it, exports every variable and starts one R script
  inside the conda environment.
- **The R scripts** take their paths from the environment, through
  `00_common.R`, and hand results to the next step through files, not through a
  shared session. Each step ends with `save(..., file = ".../NN_*.RData")` and
  the next starts with `load(...)`. That is why a step can be re-run alone.

Steps are numbered in execution order. `x1`–`x6` are analyses that test a
finding of the main steps rather than produce one.

## 2. Shell scripts

### `config.sh`

```bash
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
```

`${BASH_SOURCE[0]}` is the path of the file being executed *or sourced*, where
`$0` would give the caller's name when the file is sourced. `dirname` strips the
file name, `cd ... && pwd` turns the result into an absolute path. The project
therefore works wherever it is cloned.

```bash
UPSTREAM_DIR="${RNA_UPSTREAM_DIR:-${HOME}/rna_work}"
```

`${VAR:-default}` uses `VAR` if it is set and non-empty, the default otherwise.
Every path that a user might want to move is written this way, so
`export RNA_UPSTREAM_DIR=...` redirects the whole project without editing a file.

The file only defines variables, so it is safe to source from anywhere.

### `setup_env.sh`

```bash
set -euo pipefail
```

Three options in one: `-e` stops at the first failing command, `-u` makes an
unset variable an error rather than an empty string, `-o pipefail` makes a
pipeline fail if *any* stage fails, not only the last one. Without `pipefail`,
`some_command | tee log` reports success whenever `tee` does.

```bash
exec > >(tee -a "${LOG_DIR}/setup_env.log") 2>&1
```

`exec` with only redirections changes the redirections of the *current* shell
for everything that follows. `>(...)` is process substitution: a command whose
standard input is fed by the redirected output. The effect is that all output,
standard output and standard error (`2>&1`), appears on screen and is appended
(`-a`) to a log.

The installation is split in two tiers and the reason is failure handling:

```bash
mamba create -y -n rna_r -c conda-forge -c bioconda  "${CORE[@]}"
mamba install -y -n rna_r -c conda-forge -c bioconda "${EXTRA[@]}"
```

`-c conda-forge -c bioconda` lists channels in priority order. `"${CORE[@]}"`
expands an array into separate quoted words. The core packages are what every
step needs; the extra ones are the specialised packages that may lack a build
for a platform. If the batch install of the extras fails, the script retries them
one at a time, so one unavailable package does not hide which one it was.

The final check is run *inside R*, because conda reporting success does not mean
`library()` works:

```r
required <- commandArgs(trailingOnly = TRUE)
missing  <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
```

`commandArgs(trailingOnly = TRUE)` returns only the arguments after the script
name. `requireNamespace(pkg, quietly = TRUE)` returns `FALSE` instead of throwing
when a package is absent, so `vapply` yields one logical per package.
Anything still missing is installed with `BiocManager::install`, which serves CRAN
as well as Bioconductor, and the script refuses to report success while a
required package is absent.

The manifest `envs/versions_installed.txt` is written from `packageVersion()`
inside R, so it records what is installed rather than what was requested.

### `run_step.sh`

```bash
SCRIPT="$(ls -1 "${R_DIR}/${STEP}"_*.R "${ANALYSES_DIR}/${STEP}"_*.R \
          "${R_DIR}/${STEP}" "${ANALYSES_DIR}/${STEP}" 2>/dev/null | head -1 || true)"
```

Resolves `11`, `x2` or a full file name to a script. `ls -1` lists one per line,
`2>/dev/null` hides the "no such file" complaints for the patterns that do not
match, `head -1` keeps the first hit. Because `set -e` is active, `|| true`
prevents a harmless non-match from ending the script.

A detail that bites in `set -e` scripts:

```bash
# wrong under set -e: a false test returns non-zero and kills the script
[ "$1" = "--list" ] && { list_scripts; exit 0; }

# right
if [ "$1" = "--list" ]; then list_scripts; exit 0; fi
```

The launch line:

```bash
conda run --no-capture-output -n "${CONDA_ENV}" Rscript "${SCRIPT}" 2>&1 | tee "${LOGFILE}"
```

`conda run` executes a command inside an environment without activating it in the
calling shell. `--no-capture-output` streams output as it happens; by default
conda buffers it and prints everything at the end, which makes a 70-minute step
look hung.

### `run_all.sh`

```bash
STEPS=(11 12 13 14 15 16 16b 17 18 x1 x2 x3 x4 x6 98 99 99b)
```

An array, so the order is data and `--from` can walk it. The loop sets
`started=1` when it reaches the requested step and runs everything after.

Before anything runs, it counts the upstream files:

```bash
n_counts=$(ls "${COUNTS_DIR}"/*.htseq.counts 2>/dev/null | wc -l)
expected=$(awk 'NR>1 && $0 !~ /^#/ {n++} END {print n}' "${SAMPLE_SHEET}")
```

`awk` counts the sample sheet's data rows (`NR>1` skips the header, the regular
expression skips comment lines). The two numbers must agree for counts, BAMs and
allele tables, or the script stops with a message naming the upstream repository.
Failing here costs a second; failing in step 11 or step 17 would cost more.

### `summarize_results.sh`

A two-line wrapper over `run_step.sh 99`. It exists so the command named in the
documentation is the same command a reader would guess.

## 3. Shared R code: `00_common.R`

Sourced by every script, never run on its own.

**Paths** come from the environment with a default:

```r
env_or <- function(name, default) { v <- Sys.getenv(name); if (nzchar(v)) v else default }
```

`Sys.getenv` returns `""` for an unset variable, so `nzchar` is the test for
"set". A script run by hand in RStudio therefore works with the same defaults.

**`set.seed(42)`** is set once, so any step that samples (the random gene sets in
x3, MBASED's simulations, WGCNA) reproduces. Steps that need a specific seed set
it again right before sampling, so they do not depend on what ran earlier.

**`load_sample_table()`** reads `samples.tsv` and builds the table DESeq2 wants:

- column 1 is the sample name and column 2 the count file, an order
  `DESeqDataSetFromHTSeqCount` requires and does not check;
- spaces in group labels become underscores, because a space in a factor level
  surfaces much later as a broken contrast name;
- `CONDITION` is releveled with `control` first. With `~ 0 + CONDITION` no level
  is a reference, so this changes no result, but it silences DESeq2's warning and
  makes control the first colour in every plot.

**`strip_version()`** is `sub("[.].*$", "", x)`: `ENSG00000004776.15` becomes
`ENSG00000004776`. The character class `[.]` matches a literal dot without the
escaping that, written `"\\."` in a shell here-document, was once swallowed on
the way into the file.

**`spearman_test()`** wraps `cor.test(method = "spearman")`. For small samples R
computes an exact p-value that can underflow to exactly `0` for strong
correlations, and a p-value of 0 is not a true statement about any sample. When
that happens the function falls back to `exact = FALSE` (an asymptotic
approximation) and says so in the returned text. It also returns "not
computable" rather than `NA` dressed up as a result when there is no variation.

**`GROUP_COLOURS`** is one named vector of colours. A factor and a character
column are given different colours by the same palette because their level
orders differ, which once made two figures disagree about which group was red.

**`save_fig()` and `save_tab()`** write to `figures/` and `tables/` and print the
path, so a log shows what each step produced.

## 4. Step 11 — differential expression

Three statistical tools are fitted once each and queried for three contrasts,
which is how they are meant to be used: the dispersion is estimated from all
samples and shared.

### DESeq2

```r
dds <- DESeqDataSetFromHTSeqCount(sampleTable = st, directory = COUNTS_DIR,
                                  design = ~ 0 + CONDITION)
```

`~ 0 + CONDITION` removes the intercept. With an intercept R silently makes the
alphabetically first level the reference and every coefficient a difference
against it; with three groups that is a trap. Without it each coefficient is a
group mean and contrasts are stated explicitly.

```r
keep <- rowSums(counts(dds) >= MIN_COUNT) > ncol(dds) * MIN_FRACTION
```

Prefiltering: at least 10 counts in more than 30% of samples (at least 4 of 12).
It does not look at group labels, so it cannot bias the test; it removes genes no
test could resolve and so reduces the multiple-testing burden.

```r
cv <- makeContrasts(contrasts = "CONDITIONacute - CONDITIONcontrol",
                    levels = resultsNames(dds))
res <- results(dds, contrast = as.numeric(cv))
```

`makeContrasts` builds a vector of ±1 aligned to the model's coefficients;
`results(contrast = <numeric vector>)` tests that linear combination. Genes with
`padj = NA` were removed by DESeq2's own filtering and never entered the
correction, so they are dropped rather than counted.

### Transformations

```r
rld <- assay(rlog(dds)); ntd <- assay(normTransform(dds))
meanSdPlot(rld, ranks = FALSE, plot = FALSE)
```

Tests run on raw counts; the *transformed* values are only for plots and
distances. Low-count genes vary enormously in relative terms, which wrecks PCA and
clustering. `normTransform` is `log2(x + 1)` of size-factor-normalised counts;
`rlog` shrinks low-count genes towards the mean and flattens the mean–SD plot,
at the price of a much longer computation.

### edgeR

```r
y <- DGEList(counts(dds))
y <- y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y$samples$lib.size <- colSums(y$counts)
y <- estimateDisp(y, design, tagwise = TRUE)
fit <- glmQLFit(y, design)
lrt <- glmQLFTest(fit, contrast = cv)
```

`filterByExpr` is edgeR's own filter, using the design to decide how many samples
must pass. `keep.lib.sizes = FALSE` plus the explicit `lib.size` line recompute
library sizes after filtering; edgeR does not do this on its own, and TMM would
otherwise normalise against totals that include removed genes. `glmQLFit` and
`glmQLFTest` give the quasi-likelihood F-test, the strictest of the three.

### limma-voom

```r
v   <- voom(y, design)
fit <- lmFit(v, design)
tt  <- topTable(eBayes(contrasts.fit(fit, cv)), sort.by = "none", n = Inf)
```

`voom` models how variance depends on the mean count and returns a *weight* for
every observation; that is what lets a linear model built for continuous
microarray intensities handle counts. `sort.by = "none"` keeps the genes in the
input order so the three tools' tables line up row for row.

### Thresholds and Venn diagrams

A gene is called when `|log2FC| > 1` **and** `FDR < 0.05`: one condition for the
size of the effect, one for its reliability. Each contrast is shown twice, once
with the correction and once with raw p-values, to show what the correction
removes.

## 5. Step 12 — annotation and enrichment

```r
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
BIO <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
             filters = "ensembl_gene_id", values = ids, mart = ensembl)
```

`attributes` are the columns wanted, `filters` the column the query is keyed on,
`values` the list to look up. GENCODE IDs carry a version suffix that biomaRt does
not accept, so `strip_version` runs first, and silently missing genes are checked
by comparing the number queried with the number returned.

**Symbol matrix.** Genes without a symbol come back as the *empty string*, not
`NA`, so three filters are needed: not `NA`, not `""`, not duplicated. Filtering
only on `is.na` leaves a row named `""`.

### WebGestaltR

```r
WebGestaltR(enrichMethod = "ORA", organism = "hsapiens",
            enrichDatabase = db, interestGene = deg_ids, interestGeneType = "ensembl_gene_id",
            referenceGene = background_ids, referenceGeneType = "ensembl_gene_id",
            isOutput = TRUE, outputDirectory = dir, projectName = name)
```

ORA is a Fisher's exact test of whether a gene list overlaps a category more than
a random list from the same *background* would. The background is therefore part
of the question, and here it is the genes that entered the tests, not the whole
genome (see x3 for how much that mattered).

Databases are named (`"geneontology_Biological_Process"`), not chosen by position
in `listGeneSet()`, whose order changes when WebGestalt adds a database.

Two failures worth knowing:

1. WebGestaltR packages its HTML report with `utils::zip()`, which reads the
   path from `R_ZIPCMD`. Conda's R leaves that empty, so the call fails with
   `'zip' must be a non-empty character string` *after* the enrichment has been
   computed. The script sets `R_ZIPCMD` from `Sys.which("zip")`.
2. The first version caught errors with `error = function(e) NULL` and counted
   the result as "0 enriched sets". Ten failed runs were recorded as negative
   findings. An error is now stored as `NA`, kept apart from a genuine empty
   result, and the script warns how many rows are `NA`.

## 6. Step 13 — deconvolution

quanTIseq fits observed expression against a signature of marker genes for ten
immune cell types, by constrained regression. The script first *measures* how
much of the signature exists in the data:

```r
sig <- read.delim(system.file("extdata", "TIL10_signature.txt", package = "quantiseqr"))
present <- intersect(sig[[1]], rownames(ntd_rename))
```

13 of 170 markers are on chromosome 19. The method does not treat a missing
marker as missing; it enters the fit as a gene with no expression, which is a
statement about the mixture. That is why the output is reported with a sanity
check against what the sample is (PBMC contain no neutrophils by definition)
instead of as an estimate.

The call itself:

```r
run_quantiseq(expression_data = mat, signature_matrix = "TIL10",
              is_arraydata = FALSE, is_tumordata = FALSE, scale_mRNA = TRUE)
```

`is_arraydata = FALSE` says the data are RNA-seq, `is_tumordata = FALSE` turns off
tumour-specific signature filtering, `scale_mRNA = TRUE` corrects for the
different amounts of mRNA that cell types contain. It is run on the log scale
(as in the course) and on the linear scale (as the method expects), and the
script reports that most of their agreement is shared zeros.

## 7. Step 14 — WGCNA

The idea: correlate every gene with every other, read the correlation matrix as a
graph, cut it into modules, summarise each module with one number per sample (its
*eigengene*, the first principal component of the module), and test modules
instead of genes. A few dozen tests replace thousands.

```r
a_ntd <- t(ntd_rename)        # WGCNA wants samples in rows, genes in columns
gsg <- goodSamplesGenes(a_ntd)
```

`goodSamplesGenes` finds genes and samples with too many missing or constant
values.

**Soft threshold.**

```r
sft <- pickSoftThreshold(a_ntd, powerVector = c(1:15, seq(16, 20, by = 2)))
signed_R2 <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
```

Correlations are raised to a power so strong ones stay strong and weak ones
collapse towards zero. The power is chosen so the resulting network is
approximately scale-free (a few highly connected genes, many weakly connected).
`fitIndices[, 2]` is the fit R², column 3 the slope; multiplying by `-sign(slope)`
rewards a negative slope, as scale-free requires. The rule is the first power
whose signed R² reaches 0.85: 5 here. The course materials give 6, 4 and 7 in
three places, because it is a property of the data.

```r
cor <- WGCNA::cor
net <- blockwiseModules(a_ntd, power = 5, networkType = "unsigned", TOMType = "unsigned",
                        minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = FALSE, randomSeed = 999)
cor <- stats::cor
```

- `unsigned`: a gene that rises and one that falls together land in one module.
- `TOMType`: topological overlap, which scores two genes by how many neighbours
  they share, not only by their direct correlation.
- `minModuleSize`: smaller groups are merged or left unassigned.
- `mergeCutHeight = 0.25`: modules whose eigengenes correlate above 0.75 merge.
- `reassignThreshold = 0`: no reassignment of genes between modules after cutting.
- `randomSeed`: fixes the tie-breaking for reproducibility.

The `cor <- WGCNA::cor` line is not decoration. DESeq2 loads an S4 generic for
`cor` that lacks the arguments `blockwiseModules` passes (`weights.x`, `cosine`),
and the call dies with "unused arguments". Binding WGCNA's own version during the
build, and handing `stats::cor` back afterwards, is the fix the package documents.

**Module against traits.** The traits are *measured*: group indicators, intronic
fraction, exonic fraction, ribosomal share, library size. Cell-type fractions are
not used, because step 13 showed they are not interpretable here.

```r
mtc <- cor(MEs, traits, use = "p", method = "spearman")
mtp <- corPvalueStudent(mtc, nrow(traits))
mtq <- matrix(p.adjust(mtp, method = "BH"), nrow = nrow(mtp))
```

`use = "p"` is "pairwise complete". Spearman because indicator traits and 12
samples do not justify assuming linearity. The heatmap has one test per module per
trait, so Benjamini–Hochberg is applied across the whole matrix.

## 8. Step 15 — BioNERO

BioNERO wraps WGCNA, so this is two implementations of one idea rather than two
methods. That makes it a check on the module structure.

```r
SFT_fit(se, net_type = "unsigned", cor_method = "spearman")$power
net <- exp2gcn(se, net_type = nt, SFTpower = p, cor_method = cm,
               module_merging_threshold = 0.8)
```

All six combinations of network type (unsigned, signed, signed hybrid) and
correlation (Pearson, Spearman) are tried and reported first: on the lecturer's
four samples only unsigned converged, so convergence is checked rather than
assumed. Here all six converged.

Two BioNERO networks are built: one at its own choice (unsigned, Spearman) and
one forced onto the settings step 14 actually used (unsigned, Pearson, power 5).
Reporting only the first would blame the packages for a difference in settings;
the first version of this script did exactly that and reported an adjusted Rand
index of 0.09 where the matched comparison gives 0.27.

**Adjusted Rand index**, written out in the script:

```r
tab <- table(a, b); n <- sum(tab)
sum_ij <- sum(choose2(tab)); sum_a <- sum(choose2(rowSums(tab))); sum_b <- sum(choose2(colSums(tab)))
expected <- sum_a * sum_b / choose2(n)
(sum_ij - expected) / ((sum_a + sum_b) / 2 - expected)
```

It counts gene pairs placed together in both partitions and corrects for the
agreement expected by chance: 1 is identical, 0 is chance.

BioNERO's offer to remove the variance of the first principal components as a
batch correction is deliberately not used. Here the leading axis is library
composition, 93% explained by group; removing it would delete the group
difference and return a clean-looking network that had answered nothing.

## 9. Step 16 — allele-specific expression

Input per sample is the table GATK's `ASEReadCounter` wrote: per SNP, the reads
supporting the reference and the alternate allele.

```r
ase <- read_tsv(f) |> filter(totalCount >= 10, !is.na(refAllele), !is.na(altAllele))
snv <- GRanges(seqnames = ase$contig, ranges = IRanges(ase$position, width = 1))
hits <- findOverlaps(snv, genes_gr, ignore.strand = TRUE)
```

A SNP is one base, hence `width = 1`. `findOverlaps` assigns each SNP to the gene
it falls in; `ignore.strand = TRUE` because the question is position, not strand.
MBASED works per gene, pooling the SNPs inside it.

```r
ASE_in <- SummarizedExperiment(
    assays = list(lociAllele1Counts = ref_counts, lociAllele2Counts = alt_counts),
    rowRanges = rr)
res <- runMBASED(ASE_in, isPhased = FALSE, numSim = 1e5, BPPARAM = SerialParam())
```

`isPhased = FALSE` because the genotypes are not phased (it is not known which
allele sits on which parental chromosome), which costs information.
`numSim` is the number of simulations behind each p-value. `SerialParam()` runs
one process, slower but stable.

**A p-value of 0 is not zero.** MBASED's p-values are simulated, so the smallest
it can report is `1/numSim`, here 10⁻⁵. Between 22% and 46% of the genes
(by sample) come back at that limit. An earlier version of the plotting code replaced those
zeros with the smallest positive floating-point number (about 2·10⁻³⁰⁸) and
drew `-log10` of it, which put most genes on a line at 308. The figures are now
drawn by `16b_ase_figures.R` from the saved table, with the floor at `1/numSim`,
and say how many genes sit on it. Splitting them out also means a plotting
mistake no longer costs a half-hour MBASED run.

## 10. Step 17 — splicing

SGSeq asks a question short reads can answer: how many reads cross this exon–exon
junction and how many fall on this exon bin. These are counts, not reconstructed
transcripts.

```r
si  <- getBamInfo(samples_df, cores = 1)        # columns sample_name and file_bam, exactly
txf <- convertToTxFeatures(importTranscripts(ANNOTATION_GTF))
sgfc      <- analyzeFeatures(si, features = txf, cores = 5)
sgvc_pred <- analyzeVariants(sgfc, cores = 5)
sgvc      <- getSGVariantCounts(rowRanges(sgvc_pred), sample_info = si, cores = 5)
```

`getBamInfo` is given `cores = 1`: it only reads BAM headers, and it is where
mismatches between the requested and available core count broke the lecturer's own
class. `analyzeFeatures` builds the splice graph and counts reads; `analyzeVariants`
turns the graph into events (skipped exon, alternative donor and so on);
`getSGVariantCounts` counts the reads supporting each event. Each of the three is
cached with `saveRDS`, because together they take over an hour.

The event matrix is then analysed as ordinary expression, with `filterByExpr`,
`voom` and `lmFit`.

`quiet_run()` counts and summarises SGSeq's "events exceed maxnvariant" warnings
(one per affected gene, hundreds in all) instead of letting them bury the results,
and `suppressMessages` removes the S4 `[updateObject]` chatter. They are counted,
not discarded.

Two corrections that are worth knowing:

- `txName` is an *accessor* on an `SGFeatures` object, not a column of
  `mcols()`. `mcols(feat)$txName` is `NULL`, `NULL %in% x` is `logical(0)`, and
  indexing a matrix with it silently returns zero rows, which printed
  "0 of 0 novel features". The fix is `SGSeq::txName(feat)` with a
  `stopifnot` on the length.
- `analyzeFeatures(features = txf)` quantifies the supplied annotation and does
  **not** predict novel events; that needs the annotation left out. The
  novel-feature count (8 of 115,821) is therefore a property of the run mode, so
  the test built on it is void, not negative, and is labelled so in the log.

`plotFeatures` is broken in the installed SGSeq (it matches `1:`/`2:` prefixes
against `J:`/`E:` names and gets all `NA`). The course supplies a patched
function, included as `R/plotFeatures_fixed.R` and marked as the lecturer's.

## 11. Step 18 — power

```r
rnapower(depth = 133, n = 4, cv = 0.317, effect = 2, alpha = 0.05)
```

`depth` is the expected read count of the gene, `n` the replicates per group,
`cv` the within-group coefficient of variation, `effect` the fold change to
detect, `alpha` the significance level of one test.

The course draws curves at `cv = 0.4`, a figure quoted for mice. The lecturer
then said the right value is one measured from data of the same kind, so this
script measures it, per gene within each group, from normalised counts. That
number (0.24) is then checked against edgeR's own estimate, because the naive
standard deviation of four numbers is biased downwards and the whole conclusion
rests on it. edgeR's biological coefficient of variation is the square root of
its dispersion (0.32 for the common dispersion), and the script uses the larger,
more cautious one.

**Alpha.** `rnapower`'s `alpha` is the level of *one* test, and the course passes
0.05. This analysis tests 1,292 genes. Power is therefore reported as a bracket:
`alpha = 0.05` with no correction (an upper bound), and `0.05 / 1292` (Bonferroni,
a lower bound). Benjamini–Hochberg lies between the two, and pinning it down would
need an assumption about how many genes truly change, which these data cannot
supply. The depth loop at the end shows that more sequencing barely helps, so
variability and replicate number are what limit this design.

## 12. The cross-checks, x1 to x6

Each tests something a main step asserted.

**x1 — tool agreement.** Spearman correlation of log2 fold change over *every*
shared gene (not only the significant ones, so the measure does not depend on a
threshold), Jaccard index of the significant sets, and a Wilcoxon test of whether
genes called by only some tools are expressed lower than genes all three agree on.

**x2 — composition.** Four tests that can separate "biology" from "library
preparation": complete separation of groups by intronic fraction; the ribosomal
share against intronic fraction among the patients only, where group cannot
explain it; Fisher's exact test of ribosomal genes among the DEGs; and a refit
with intronic fraction as a covariate, printed with a warning that a collapse
means inseparability, not refutation.

**x3 — ORA background.** The same DEG lists against the chromosome-19 and
whole-genome backgrounds, then random gene sets drawn from genes that are not
DE in any contrast and sent through both. Every WebGestalt call is cached by
project name. The result was that the choice changed little, and that eight random
sets produced no enrichment under either; the script reports it as it came out,
and the header of step 12 was rewritten to match.

**x4 — ASE.** Two hypotheses make opposite predictions: if the imbalance is
biology, it should be reproducible per gene and unrelated to how each library was
made; if the heterozygous calls are wrong, the per-sample rate should follow the
rate of doubtful calls and the flagged genes should sit at extreme allele
fractions. Both predictions of the second held.

**x6 — paired design.** The three-group model ignores that acute and subacute are
the same four people. `~ PATIENT + CONDITION` on the eight patient samples can
subtract each person's baseline. Three fits (unpaired, paired, paired without
pair P3 because its genotypes disagree with the study's claim) are compared on
dispersion, number of hits, power and the p-value histogram. The residual
dispersion of the paired model is then fed to `rnapower`, since in a balanced
paired design the variance of the effect estimate has the same form as for two
independent groups.

## 13. Steps 98, 99 and 99b — figure, summary and README check

**98** draws the four-panel overview figure from tables that earlier steps wrote. It computes nothing new, so every value on the figure can be traced to its source.

**99** (`99_summarize_results.R`) reads the tables the steps wrote and builds `SUMMARY.md`.
Anything a table does not carry directly (an R², a correlation, a Jaccard index,
the soft threshold chosen) is recomputed from the rows of the table, so it traces
back to its inputs. A section whose input is missing is skipped and listed at the
end, never filled with a guess. It also copies the small figures and tables into
`showcase/`, leaving out anything over 60 KB and the per-gene matrices.

**99b** (`99b_check_readme.R`) guards the README. Prose does not update itself, so a
figure typed into a sentence goes stale the moment a step is re-run with a different
seed, threshold or data. Each figure the README quotes is recomputed from the tables,
formatted the way the README writes it, and looked up in the README text. A mismatch
is reported, and the script exits with an error, but nothing is rewritten: the README
is written by a person, and whether the sentence or the table needs to change is a
judgement the script should not make. The typographic characters in the README (minus
sign, en dash, arrow, superscript two) are mapped to ASCII first so the lookups stay
plain strings.
