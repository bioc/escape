#' Flexible GSEA for Precomputed Gene Lists
#'
#' A convenience front-end to \pkg{fgsea} that lets you point at the
#' \code{avg_log2FC} and \code{p_val_adj} columns coming out of Seurat / DESeq2
#' / edgeR etc. It converts them to a signed -log10(p) ranking, filters on
#' significance / effect size, and then runs fgsea.
#'
#' @param input.data Either:
#'   \itemize{
#'     \item A named numeric vector \strong{already ranked}, or
#'     \item A data.frame/tibble with one row per gene and columns containing
#'       log-fold-change and p-value. If the gene ID is not in
#'       \code{rownames(data)}, supply \code{gene_col}.
#'   }
#' @param gene.sets A named list of character vectors, the result of
#'   \code{\link{getGeneSets}}, or the built-in data object
#'   \code{\link{escape.gene.sets}}.
#' @param gene_col Character or \code{NULL}. Name of the column holding gene
#'   identifiers (ignored when they are row-names). Default is \code{NULL}.
#' @param logFC_col Character. Column name for log-fold-change values. Default
#'   is \code{"avg_log2FC"} (matches Seurat's \code{FindMarkers()}).
#' @param pval_col Character. Column name for p-values (or adjusted p-values).
#'   Default is \code{c("p_val_adj", "p_val")} (first match is used).
#' @param ranking_fun Character. How to build the ranking:
#'   \itemize{
#'     \item \code{"signed_log10_p"} (default): sign(logFC) * -log10(p).
#'     \item \code{"logFC"}: Use log-fold-change values directly.
#'   }
#' @param pval_cutoff Numeric. Filter genes with p-value above this threshold
#'   \strong{before} ranking. Default is \code{1} (no filtering).
#' @param logFC_cutoff Numeric. Filter genes with absolute log-fold-change
#'   below this threshold \strong{before} ranking. Default is \code{0} (no
#'   filtering).
#' @param minSize Integer. Minimum pathway size passed to \pkg{fgsea}. Default
#'   is \code{5}.
#' @param maxSize Integer. Maximum pathway size passed to \pkg{fgsea}. Default
#'   is \code{500}.
#' @param padjust_method Character. Multiple-testing correction method; any
#'   method accepted by \code{\link[stats]{p.adjust}}. Default is \code{"BH"}.
#' @param nproc Integer. Passed to \pkg{fgsea}. Use \code{0} for multithread
#'   if OpenMP is available. Default is \code{0}.
#'
#' @seealso \code{\link[fgsea]{fgsea}}, \code{\link{getGeneSets}},
#'   \code{\link{gseaEnrichment}}
#'
#' @examples
#' pbmc_small <- SeuratObject::pbmc_small
#' 
#' Seurat::Idents(pbmc_small) <- "groups"
#' markers <- Seurat::FindMarkers(pbmc_small, 
#'                                ident.1 = "g1", 
#'                                ident.2 = "g2")
#' 
#' gs <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#'                
#' gsea <- enrichIt(markers, 
#'                  gene.sets = gs)
#' 
#' @return `data.frame` with the usual fgsea columns plus a convenient
#' `leadingEdge` character column collapsed with \";\".
#' @importFrom stats p.adjust
#' @export
enrichIt <- function(input.data,
                     gene.sets,
                     gene_col       = NULL,
                     logFC_col      = "avg_log2FC",
                     pval_col       = c("p_val_adj", "p_val"),
                     ranking_fun    = c("signed_log10_p", "logFC"),
                     pval_cutoff    = 1,
                     logFC_cutoff   = 0,
                     minSize        = 5,
                     maxSize        = 500,
                     padjust_method = "BH",
                     nproc          = 0) {
  
  if (!requireNamespace("fgsea", quietly = TRUE))
    stop("Package 'fgsea' is required.")
  
  ranking_fun <- match.arg(ranking_fun)
  
  ## ------------------------------------------------------------------------
  ## 1. Build/validate the STATISTIC vector
  ## ------------------------------------------------------------------------
  if (is.numeric(input.data) && !is.null(names(input.data))) {
    stats <- sort(input.data[!is.na(input.data)], decreasing = TRUE)
    
  } else if (is.data.frame(input.data)) {
    
    df <- input.data
    
    ## decide which p-value column to use ------------------------
    pval_col <- match.arg(pval_col[ pval_col %in% names(df) ],
                          choices = pval_col)
    
    ## pull gene IDs --------------------------------------------
    if (is.null(gene_col)) {
      if (is.null(rownames(df)))
        stop("Gene IDs must be in row.names or specify 'gene_col'.")
      gene_ids <- rownames(df)
    } else {
      if (!gene_col %in% names(df))
        stop("'gene_col' not found in data.")
      gene_ids <- df[[gene_col]]
    }
    
    ## sanity ----------------------------------------------------
    if (!all(c(logFC_col, pval_col) %in% names(df)))
      stop("Specified 'logFC_col' or 'pval_col' not in data.")
    
    ## filter ----------------------------------------------------
    keep <- !is.na(df[[logFC_col]]) &
      !is.na(df[[pval_col]])  &
      df[[pval_col]] <= pval_cutoff &
      abs(df[[logFC_col]])   >= logFC_cutoff
    df   <- df[keep, ]
    gene_ids <- gene_ids[keep]
    
    if (nrow(df) == 0)
      stop("No genes left after filtering (check cut-offs).")
    
    ## build ranking --------------------------------------------
    stat_vec <- switch(ranking_fun,
                       signed_log10_p = sign(df[[logFC_col]]) * -log10(df[[pval_col]]),
                       logFC          = df[[logFC_col]]
    )
    stats <- setNames(stat_vec, gene_ids)
    stats <- sort(stats, decreasing = TRUE)
    
  } else {
    stop("'data' must be a named numeric vector or a data.frame.")
  }
  
  ## ------------------------------------------------------------------------
  ## 2. Harmonise gene-sets (escape utility) & run fgsea
  ## ------------------------------------------------------------------------
  gene.sets <- .GS.check(gene.sets)
  
  ## Decide scoreType automatically ----------------------------------------
  score_type <- if (all(stats >= 0)) {
    "pos"                     # every value >=0
  } else if (all(stats <= 0)) {
    "neg"                     # every value <=0
  } else {
    "std"                     # mixture of positive and negative
  }
  
  res <- fgsea::fgsea(
    pathways = gene.sets,
    stats    = stats,
    minSize  = minSize,
    maxSize  = maxSize,
    nproc    = nproc,
    scoreType = score_type)
  
  ## tidy --------------------------------------------------------
  res$geneRatio <- vapply(res$leadingEdge, length, integer(1L)) / res$size
  res$leadingEdge <- vapply(res$leadingEdge,
                            paste, collapse = ";", character(1))
  res$padj <- p.adjust(res$pval, method = padjust_method)
  res <- res[order(res$padj, res$pval), ]
  rownames(res) <- NULL
  res
}