# showcase/

Small figures and tables that show the pipeline ran and what it produced. They
are copies of files written under `~/rna_work_r/results` and are refreshed by
`bash summarize_results.sh`; none is edited by hand. Everything larger than
60 KB, and every per-gene matrix, stays out of the repository.

## figures/

| file | what it shows |
|---|---|
| `98_overview.png` | the four-panel summary of the project: library composition, where the hits are, ASE against genotype quality, and the power bracket |
| `11_normalisation_comparison.png` | mean against standard deviation under log2, log2 with size factors, and rlog |
| `11_pca.png` | PCA of the twelve libraries on rlog values |
| `11_venn_acute_vs_control.png`, `11_venn_acute_vs_subacute.png` | genes called by DESeq2, edgeR and limma-voom, with and without multiple-testing correction |
| `12_biotypes.png` | gene biotypes among the genes that entered the tests |
| `13_cell_fractions.png` | quanTIseq output on the 13 signature genes present on chr19; shown as evidence of the limit, not as an estimate |
| `14_soft_threshold.png` | WGCNA scale-free fit and mean connectivity against the soft-threshold power |
| `14_module_trait_heatmap.png` | WGCNA module eigengenes against measured library traits |
| `16_ase_volcano.png`, `16_ase_qq.png` | MBASED imbalance against significance and the p-value QQ plot, on one sample, with the simulation's resolution limit marked |
| `17_features_heatmap.png` | read counts on the exons and junctions of the gene at the top of the splicing results |
| `17_junction_ratio_vs_intronic.png` | junction-to-exon read ratio against intronic fraction (see the note in `R/17_splicing_sgseq.R`: this metric cannot see intron retention) |
| `18_power_curves.png` | power against replicates at three values of within-group variability, uncorrected per-test alpha |
| `x1_lfc_scatter.png`, `x1_disagreement_vs_expression.png` | agreement between DESeq2 and edgeR; expression of genes the tools disagree on |
| `x2_rp_vs_intronic.png` | share of reads on ribosomal protein genes against intronic fraction |
| `x3_background_comparison.png`, `x3_random_sets.png` | enriched sets under the two ORA backgrounds; random gene sets through both |
| `x4_ase_vs_miscall.png`, `x4_min_alt_by_significance.png` | apparent ASE against the share of doubtful heterozygous calls; allele fractions of flagged and unflagged genes |
| `x6_pvalue_histograms.png` | p-value calibration of the within-patient contrast in three models |

## tables/

The tables are named `<step>_<content>.tsv`, with the step number or `x` code of
the script that wrote them. [`../SUMMARY.md`](../SUMMARY.md) is built from them
and is the easier place to read them; the files are here so each number in it can
be traced to its source.
