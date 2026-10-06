# Results

Every number below is read or recomputed from the pipeline's own output by
`R/99_summarize_results.R`. Regenerate with `bash summarize_results.sh`.

**Contents**

- [Differential expression, three tools](#differential-expression-three-tools)
- [Library composition differs between groups](#library-composition-differs-between-groups)
- [Annotation and enrichment](#annotation-and-enrichment)
- [Cell-type deconvolution](#cell-type-deconvolution)
- [Co-expression networks](#co-expression-networks)
- [Allele-specific expression](#allele-specific-expression)
- [Alternative splicing](#alternative-splicing)
- [Statistical power](#statistical-power)
- [Within-patient contrast, paired design](#within-patient-contrast-paired-design)

## Differential expression, three tools

Genes with |log2 FC| > 1 and FDR < 0.05 (`fdr`), or with raw p < 0.05 and no correction (`raw_p`).

| contrast | fdr_DESeq2 | fdr_edgeR | fdr_limma | raw_p_DESeq2 | raw_p_edgeR | raw_p_limma |
|---|---|---|---|---|---|---|
| acute_vs_control | 296 | 291 | 275 | 320 | 317 | 296 |
| subacute_vs_control | 471 | 459 | 448 | 491 | 480 | 479 |
| acute_vs_subacute | 2 | 0 | 0 | 31 | 33 | 28 |

Agreement between tools on log2 fold change over every shared gene (Spearman), and on the significant sets (Jaccard; NA when both lists are empty).

| contrast | pair | genes | spearman_lfc | jaccard_deg |
|---|---|---|---|---|
| acute_vs_control | DESeq2 vs edgeR | 1,284 |     1 | 0.9567 |
| acute_vs_control | DESeq2 vs limma | 1,284 | 0.9896 | 0.8721 |
| acute_vs_control | edgeR vs limma | 1,284 | 0.9894 | 0.8742 |
| subacute_vs_control | DESeq2 vs edgeR | 1,284 |     1 | 0.9456 |
| subacute_vs_control | DESeq2 vs limma | 1,284 | 0.9915 | 0.9027 |
| subacute_vs_control | edgeR vs limma | 1,284 | 0.9916 | 0.9257 |
| acute_vs_subacute | DESeq2 vs edgeR | 1,284 | 0.9998 | 0 |
| acute_vs_subacute | DESeq2 vs limma | 1,284 | 0.9591 | 0 |
| acute_vs_subacute | edgeR vs limma | 1,284 |  0.96 | NA |


## Library composition differs between groups

Intronic fraction and the share of reads on ribosomal protein genes, by group.

| group | samples | intronic_min | intronic_max | rp_share_min | rp_share_max |
|---|---|---|---|---|---|
| acute | 4 | 7.431 |  14.2 | 9.861 | 15.96 |
| control | 4 |  21.2 | 26.99 | 4.713 | 6.662 |
| subacute | 4 | 7.337 |  9.03 | 14.05 | 16.06 |

| quantity | value |
|---|---|
| groups separated completely by intronic fraction | yes |
| variance of intronic fraction explained by group (R2) | 0.93 |
| ribosomal share vs intronic fraction, all 12 samples | rho = -0.902  p = 6e-05 [exact routine underflowed; asymptotic approximation] |
| ribosomal share vs intronic fraction, patients only (n = 8) | rho = -0.881  p = 0.00724 |

Ribosomal protein genes among the DEGs (Fisher's exact test).

| contrast | rp_genes | rp_among_deg | deg_total | odds_ratio | p_value | median_lfc_rp |
|---|---|---|---|---|---|---|
| acute_vs_control | 12 | 11 | 296 | 38.29 | 7.539e-07 | 1.431 |
| subacute_vs_control | 12 | 11 | 471 | 19.57 | 0.0001129 | 1.801 |
| acute_vs_subacute | 12 | 0 | 2 | 0 | 1 | -0.3293 |

DEG count when intronic fraction is added to the model. The covariate is largely explained by group, so a collapse reflects inseparability, not refutation.

| contrast | deg_uncorrected | deg_with_covariate | retained_pct |
|---|---|---|---|
| acute_vs_control | 296 | 4 |   1.4 |
| subacute_vs_control | 471 | 5 |   1.1 |
| acute_vs_subacute | 2 | 1 | 50 |

Share of variance carried by the leading principal components (rlog values, top 500 genes).

| component | percent_variance |
|---|---|
| PC1 | 61 |
| PC2 | 11 |


## Annotation and enrichment

Gene biotypes among the genes that entered the tests.

| biotype | genes |
|---|---|
| protein_coding | 967 |
| lncRNA | 227 |
| processed_pseudogene | 68 |
| TEC | 16 |
| transcribed_processed_pseudogene | 5 |
| transcribed_unprocessed_pseudogene | 3 |
| unprocessed_pseudogene | 3 |
| miRNA | 1 |
| misc_RNA | 1 |
| snRNA | 1 |

Enriched gene sets per contrast and database (WebGestalt ORA, FDR < 0.05, background = genes tested). `NA` would mean the run failed.

| contrast | database | deg_in | enriched |
|---|---|---|---|
| acute_vs_control | geneontology_Biological_Process | 296 | 25 |
| acute_vs_control | geneontology_Cellular_Component | 296 | 7 |
| acute_vs_control | geneontology_Molecular_Function | 296 | 4 |
| acute_vs_control | pathway_KEGG | 296 | 1 |
| acute_vs_control | pathway_Reactome | 296 | 32 |
| acute_vs_control | network_miRNA_target | 296 | 0 |
| acute_vs_control | network_Transcription_Factor_target | 296 | 0 |
| subacute_vs_control | geneontology_Biological_Process | 471 | 12 |
| subacute_vs_control | geneontology_Cellular_Component | 471 | 3 |
| subacute_vs_control | geneontology_Molecular_Function | 471 | 2 |
| subacute_vs_control | pathway_KEGG | 471 | 1 |
| subacute_vs_control | pathway_Reactome | 471 | 28 |
| subacute_vs_control | network_miRNA_target | 471 | 0 |
| subacute_vs_control | network_Transcription_Factor_target | 471 | 0 |

The same DEG lists against the chr19 background and against the whole genome.

| contrast | database | degs | ours | genome | shared | only_ours | only_genome |
|---|---|---|---|---|---|---|---|
| acute_vs_control | geneontology_Biological_Process | 296 | 25 | 14 | 12 | 13 | 2 |
| acute_vs_control | geneontology_Cellular_Component | 296 | 7 | 8 | 5 | 2 | 3 |
| acute_vs_control | geneontology_Molecular_Function | 296 | 4 | 1 | 1 | 3 | 0 |
| acute_vs_control | pathway_KEGG | 296 | 1 | 1 | 1 | 0 | 0 |
| acute_vs_control | pathway_Reactome | 296 | 32 | 34 | 25 | 7 | 9 |
| acute_vs_control | network_miRNA_target | 296 | 0 | 0 | 0 | 0 | 0 |
| acute_vs_control | network_Transcription_Factor_target | 296 | 0 | 0 | 0 | 0 | 0 |
| subacute_vs_control | geneontology_Biological_Process | 471 | 12 | 11 | 11 | 1 | 0 |
| subacute_vs_control | geneontology_Cellular_Component | 471 | 3 | 8 | 3 | 0 | 5 |
| subacute_vs_control | geneontology_Molecular_Function | 471 | 2 | 1 | 1 | 1 | 0 |
| subacute_vs_control | pathway_KEGG | 471 | 1 | 1 | 1 | 0 | 0 |
| subacute_vs_control | pathway_Reactome | 471 | 28 | 28 | 24 | 4 | 4 |
| subacute_vs_control | network_miRNA_target | 471 | 0 | 0 | 0 | 0 | 0 |
| subacute_vs_control | network_Transcription_Factor_target | 471 | 0 | 1 | 0 | 0 | 1 |

Random gene sets drawn from genes that are not DE in any contrast, through both backgrounds. Zero draws with enrichment means neither background produced a false one.

| background | draws | draws_with_any_enrichment | median_sets | max_sets |
|---|---|---|---|---|
| genome | 8 | 0 | 0 | 0 |
| ours | 8 | 0 | 0 | 0 |


## Cell-type deconvolution

Coverage of the quanTIseq TIL10 signature on chromosome 19, and the estimate against what a PBMC sample is.

| signature | markers_total | markers_present | coverage_pct | genes_in_matrix |
|---|---|---|---|---|
| TIL10 | 170 | 13 |  7.65 | 1,109 |

| sample | mononuclear_pct | neutrophil_pct | unassigned_pct |
|---|---|---|---|
| D1_1 | 0 | 100 | 0 |
| D1_2 |  53.5 |  46.5 | 0 |
| D1_3 |  17.8 | 66.68 | 0 |
| D1_4 | 0 | 84.05 | 0 |
| D2_1 | 0 | 100 | 0 |
| D2_2 | 22.42 | 77.58 | 0 |
| D2_3 | 0 | 100 | 0 |
| D2_4 | 0 | 100 | 0 |
| NC_1 | 0 | 100 | 0 |
| NC_3 | 0 | 100 | 0 |
| NC_4 | 0 | 100 | 0 |
| NC_5 |   2.8 | 90.53 | 0 |

| quantity | percent |
|---|---|
| mean mononuclear fraction |  8.04 |
| mean neutrophil fraction | 88.78 |


## Co-expression networks

WGCNA soft threshold: the first power whose signed scale-free R2 reaches 0.85.

| quantity | value |
|---|---|
| chosen power | 5 |
| signed R2 at that power | 0.883 |

WGCNA modules.

| module | genes |
|---|---|
| turquoise | 701 |
| blue | 235 |
| grey | 97 |
| brown | 76 |

Strongest module-trait associations (Spearman, Benjamini-Hochberg q).

| module | trait | rho | p | q_BH |
|---|---|---|---|---|
| MEturquoise | intronic_pct | -0.979 | 3.09e-08 | 4.326e-07 |
| MEturquoise | exonic_pct | 0.979 | 3.09e-08 | 4.326e-07 |
| MEturquoise | library_size | 0.9301 | 1.17e-05 | 0.0001092 |
| MEturquoise | rp_fraction | 0.8671 | 0.0002598 | 0.001819 |
| MEblue | rp_fraction | -0.8601 | 0.0003317 | 0.001857 |
| MEturquoise | control | -0.8193 | 0.001109 | 0.005177 |
| MEblue | intronic_pct | 0.7483 | 0.005124 | 0.01793 |
| MEblue | exonic_pct | -0.7483 | 0.005124 | 0.01793 |

Agreement between the WGCNA and BioNERO module partitions (adjusted Rand index; 1 = identical, 0 = chance).

| comparison | bionero_cor | bionero_power | bionero_modules | wgcna_modules | genes_compared | adjusted_rand |
|---|---|---|---|---|---|---|
| matched settings | pearson | 5 | 9 | 4 | 1,109 | 0.2712 |
| each package's default | spearman | 5 | 6 | 4 | 1,109 | 0.0911 |

| package | strongest_abs_rho_with_intronic_fraction |
|---|---|
| WGCNA | 0.979 |
| BioNERO | 0.937 |

Hub genes from the zinc-finger families.

| hub_genes | zinc_finger | pct |
|---|---|---|
| 61 | 26 |  42.6 |


## Allele-specific expression

MBASED, all twelve samples, genes with FDR < 0.05.

| sample | group | sites | genes | significant | pct_significant |
|---|---|---|---|---|---|
| D1_1 | subacute | 1,823 | 373 | 160 |  42.9 |
| D1_2 | subacute | 1,128 | 274 | 103 |  37.6 |
| D1_3 | subacute | 2,637 | 461 | 218 |  47.3 |
| D1_4 | subacute | 2,370 | 412 | 189 |  45.9 |
| D2_1 | acute | 3,555 | 486 | 256 |  52.7 |
| D2_2 | acute | 817 | 220 | 91 |  41.4 |
| D2_3 | acute | 1,954 | 374 | 161 | 43 |
| D2_4 | acute | 3,165 | 484 | 239 |  49.4 |
| NC_1 | control | 6,965 | 637 | 436 |  68.4 |
| NC_3 | control | 6,612 | 624 | 388 |  62.2 |
| NC_4 | control | 6,325 | 608 | 377 | 62 |
| NC_5 | control | 6,629 | 566 | 362 | 64 |

| quantity | value |
|---|---|
| genes called imbalanced, all samples pooled (%) | 54 |

Share of gene-level p-values below the simulation resolution (1 / numSim); these are reported as p = 0 by MBASED and drawn at the limit.

| sample | genes | at_floor | pct_at_floor |
|---|---|---|---|
| D1_1 | 373 | 84 |  22.5 |
| D1_2 | 274 | 62 |  22.6 |
| D1_3 | 461 | 131 |  28.4 |
| D1_4 | 412 | 108 |  26.2 |
| D2_1 | 486 | 171 |  35.2 |
| D2_2 | 220 | 51 |  23.2 |
| D2_3 | 374 | 96 |  25.7 |
| D2_4 | 484 | 151 |  31.2 |
| NC_1 | 637 | 292 |  45.8 |
| NC_3 | 624 | 262 | 42 |
| NC_4 | 608 | 262 |  43.1 |
| NC_5 | 566 | 249 | 44 |

| quantity | value |
|---|---|
| per-sample % imbalanced vs % of doubtful heterozygous calls | rho = 0.949  p = 2.44e-06 |

Allele fractions of genes MBASED calls imbalanced against the rest.

| sig | genes | median_mean_alt | median_min_alt | pct_with_a_site_below_0.1 |
|---|---|---|---|---|
| FALSE | 2,389 | 0.4827 | 0.4286 |   0.6 |
| TRUE | 2,860 | 0.2772 | 0.0769 |  57.8 |

| quantity | value |
|---|---|
| gene-sample pairs compared | 3,230 |
| called by MBASED | 2,218 |
| called by the binomial test | 2,149 |
| Jaccard index between the two | 0.891 |


## Alternative splicing

Differential splicing events (limma-voom, |logFC| > 1 and FDR < 0.05).

| contrast | events_tested | significant |
|---|---|---|
| acute_vs_control | 4,774 | 1,068 |
| subacute_vs_control | 4,774 | 1,584 |
| acute_vs_subacute | 4,774 | 0 |

| quantity | value |
|---|---|
| events detected per sample vs intronic fraction | rho = 0.769  p = 0.00525 |
| junction / exon read ratio vs intronic fraction | rho = 0.664  p = 0.0222 |


## Statistical power

Within-group coefficient of variation measured from these counts, and from edgeR's model.

| group | samples | genes | cv_median | cv_q25 | cv_q75 | cv_q90 |
|---|---|---|---|---|---|---|
| control | 4 | 1,210 | 0.232 | 0.143 | 0.357 | 0.518 |
| acute | 4 | 1,192 | 0.278 | 0.183 | 0.413 | 0.618 |
| subacute | 4 | 1,173 | 0.203 | 0.125 |  0.33 | 0.528 |

| estimator | value |
|---|---|
| naive CV of normalised counts | 0.2405 |
| edgeR common | 0.3174 |
| edgeR trended (median) | 0.3074 |
| edgeR tagwise (median) | 0.2793 |

Power to detect a twofold change at the median gene with four replicates per group, with and without a multiple-testing correction.

| assumption | alpha | power_n4 | n_for_80pct |
|---|---|---|---|
| no correction: alpha = 0.05 per test (as in the course) |  0.05 | 0.846 | 4 |
| Bonferroni: alpha = 0.05 / 1292 | 3.87e-05 | 0.128 | 12 |

| gene_quantile | cv | power_n4 | n_for_80pct |
|---|---|---|---|
| median | 0.317 | 0.847 | 4 |
| 75th percentile | 0.372 | 0.728 | 5 |
| 90th percentile | 0.573 | 0.394 | 11 |


## Within-patient contrast, paired design

Acute against subacute in the eight patient samples, with and without modelling the patient pairing, and without patient P3.

| fit | samples | genes | median_BCV | padj_lt_0.05 | and_twofold | min_padj | raw_p_lt_0.05_pct |
|---|---|---|---|---|---|---|---|
| unpaired_8 | 8 | 1,256 | 0.276 | 0 | 0 | 0.105 |   9.1 |
| paired_4 | 8 | 1,256 | 0.267 | 0 | 0 | 0.194 |  10.8 |
| paired_3 | 6 | 1,259 |  0.28 | 0 | 0 | 0.902 |   4.5 |

| fit | pairs_or_per_group | BCV | power_no_correction | power_bonferroni |
|---|---|---|---|---|
| unpaired_8 | 4 | 0.276 | 0.923 | 0.235 |
| paired_4 | 4 | 0.267 | 0.938 |  0.27 |
| paired_3 | 3 |  0.28 | 0.826 | 0.113 |

Per-pair difference in library composition, acute minus subacute.

| patient | acute | subacute | d_intronic | d_rp |
|---|---|---|---|---|
| P1 | D2_1 | D1_1 | 5.173 | -4.187 |
| P2 | D2_2 | D1_2 | 0.094 | 1.199 |
| P3 | D2_3 | D1_3 | 2.676 | -2.617 |
| P4 | D2_4 | D1_4 | 1.159 | -3.501 |

Genotype concordance of the pairs the study reports as one person (from the upstream pipeline).

| pair | concordance_nonref | called | agrees |
|---|---|---|---|
| D1_1 / D2_1 | 0.9571 | same person | OK |
| D1_4 / D2_4 | 0.9554 | same person | OK |
| D1_2 / D2_2 | 0.8965 | same person | OK |
| D1_3 / D2_3 | 0.7345 | different people | CONFLICT |

