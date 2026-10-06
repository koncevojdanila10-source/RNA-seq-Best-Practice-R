# plotFeatures_fixed.R
#
# PROVENANCE: this function is course material, distributed by the lecturer
# (Blastim RNA-seq course, day 6) as a working replacement for
# SGSeq::plotFeatures. It is reproduced here so that R/17_splicing_sgseq.R
# can be run from a clone of this repository.
#
# The only changes to the lecturer's code are this header and the translation
# of two Russian comments and one Russian error message into English.
#
# It exists because SGSeq::plotFeatures stopped working: plotSpliceGraph now
# returns feature names prefixed "1:/2:/3:/4:" while feature2name still emits
# "J:/E:/D:/A:", so the match() between the two comes back all NA and the
# heatmap has nothing to draw. The patch normalises one naming scheme onto the
# other before matching.
#
# SGSeq itself is Artistic-2.0 licensed; this derivative is included on the
# same terms and is not the work of this repository's author.
#
plotFeatures_fixed <- function(x, geneID = NULL, geneName = NULL, which = NULL,
                               tx_view = FALSE, cex = 1, assay = "FPKM",
                               include = c("junctions", "exons", "both"),
                               transform = function(x) log2(x + 1),
                               Rowv = NULL, distfun = dist, hclustfun = hclust,
                               margin = 0.2, RowSideColors = NULL,
                               square = FALSE, cexRow = 1, cexCol = 1,
                               labRow = colnames(x),
                               col = colorRampPalette(c("black", "gold"))(256),
                               zlim = NULL, heightPanels = c(1, 2), ...) {
  
  include <- match.arg(include)
  
  if (!is(x, "SGFeatureCounts")) {
    stop("x must be an SGFeatureCounts object")
  }
  
  x <- SGSeq:::restrictFeatures(x, geneID = geneID, geneName = geneName, which = which)
  if (nrow(x) == 0) return(invisible(NULL))
  
  features <- rowRanges(x)
  
  if (!is.null(RowSideColors) && !is.list(RowSideColors)) {
    RowSideColors <- list(RowSideColors)
  }
  
  n_sample <- ncol(x)
  include_type <- switch(include,
                         junctions = "J",
                         exons = "E",
                         both = c("E", "J"))
  n_feature <- length(which(type(features) %in% include_type))
  
  pars <- SGSeq:::getLayoutParameters(
    n_sample, n_feature, margin, heightPanels, RowSideColors, square, tx_view
  )
  layout(pars$mat, heights = pars$hei, widths = pars$wid)
  par(cex = cex, mai = c(pars$mai[1], 0, pars$mai[3], 0))
  
  df <- plotSpliceGraph(
    x = features,
    label = "id",
    tx_view = tx_view,
    short_output = FALSE,
    ...
  )
  
  # --- PATCH: map the "1:/2:/3:/4:" prefixes onto "J:/E:/D:/A:"
  normalize_feature_name <- function(z) {
    z <- sub("^1:", "J:", z)
    z <- sub("^2:", "E:", z)
    z <- sub("^3:", "D:", z)
    z <- sub("^4:", "A:", z)
    z
  }
  
  df$name <- normalize_feature_name(df$name)
  feat_names <- normalize_feature_name(SGSeq:::feature2name(features))
  
  i_df <- switch(include,
                 exons = which(df$type == "E"),
                 junctions = which(df$type == "J"),
                 both = seq_len(nrow(df)))
  i_df <- i_df[!duplicated(df$name[i_df])]
  
  i <- match(df$name[i_df], feat_names)
  
  # keep only the features that actually matched
  ok <- !is.na(i)
  i <- i[ok]
  i_df <- i_df[ok]
  
  if (length(i) == 0) {
    stop("no feature matched after the name mapping, nothing to draw in the heatmap")
  }
  
  X <- transform(SummarizedExperiment::assay(x[i, ], assay))
  
  labCol <- df$id[i_df]
  colLabCol <- df$color[i_df]
  
  SGSeq:::plotImage(
    X,
    Rowv = Rowv,
    distfun = distfun,
    hclustfun = hclustfun,
    RowSideColors = RowSideColors,
    cexRow = cexRow,
    cexCol = cexCol,
    labRow = labRow,
    labCol = labCol,
    colLabCol = colLabCol,
    col = col,
    zlim = zlim
  )
  
  frame()
  
  invisible(df[c("id", "name", "type", "featureID")])
}
