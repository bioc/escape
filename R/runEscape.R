#' Calculate Single-Cell Gene-Set Enrichment Scores
#'
#' \code{escape.matrix()} computes per-cell enrichment for arbitrary gene-set
#' collections using one of four scoring back-ends and returns a dense numeric
#' matrix (cells x gene-sets). The expression matrix is processed in
#' user-defined chunks (\code{groups}) so that memory use remains predictable;
#' each chunk is dispatched in parallel via a \pkg{BiocParallel} \code{BPPARAM}
#' backend. Heavy engines (\pkg{GSVA}, \pkg{UCell}, \pkg{AUCell}) are loaded
#' lazily, keeping them in the package's \strong{Suggests} field.
#'
#' @section Supported methods:
#' \describe{
#'   \item{\code{"GSVA"}}{Gene-set variation analysis (Poisson kernel).}
#'   \item{\code{"ssGSEA"}}{Single-sample GSEA.}
#'   \item{\code{"UCell"}}{Rank-based UCell scoring.}
#'   \item{\code{"AUCell"}}{Area-under-the-curve ranking score.}
#' }
#'
#' @param input.data A raw-counts matrix (genes x cells), a
#'   \link[SeuratObject]{Seurat} object, or a
#'   \link[SingleCellExperiment]{SingleCellExperiment}. Gene identifiers must
#'   match those in \code{gene.sets}.
#' @param gene.sets A named list of character vectors, the result of
#'   \code{\link{getGeneSets}}, or the built-in data object
#'   \code{\link{escape.gene.sets}}. List names become column names in the
#'   result.
#' @param method Character. Scoring algorithm (case-insensitive). One of
#'   \code{"GSVA"}, \code{"ssGSEA"}, \code{"UCell"}, or \code{"AUCell"}.
#'   Default is \code{"ssGSEA"}.
#' @param groups Integer. Number of cells per processing chunk. Larger values
#'   reduce overhead but increase memory usage. Default is \code{1000}.
#' @param min.size Integer or \code{NULL}. Minimum number of genes from a set
#'   that must be detected in the expression matrix for that set to be scored.
#'   Default is \code{5}. Use \code{NULL} to disable filtering.
#' @param normalize Logical. If \code{TRUE}, the score matrix is passed to
#'   \code{\link{performNormalization}} (drop-out scaling and optional log
#'   transform). Default is \code{FALSE}.
#' @param make.positive Logical. If \code{TRUE} \emph{and}
#'   \code{normalize = TRUE}, shifts every gene-set column so its global
#'   minimum is zero, facilitating downstream log-ratio analyses. Default is
#'   \code{FALSE}.
#' @param min.expr.cells Numeric. Gene-expression filter threshold. Default is
#'   \code{0} (no gene filtering).
#' @param min.filter.by Character or \code{NULL}. Column name in
#'   \code{meta.data} (Seurat) or \code{colData} (SCE) defining groups within
#'   which the \code{min.expr.cells} rule is applied. Default is \code{NULL}.
#' @param BPPARAM A \pkg{BiocParallel} parameter object describing the
#'   parallel backend. Default is \code{NULL} (serial execution).
#' @param ... Extra arguments passed verbatim to the chosen back-end scoring
#'   function (\code{gsva()}, \code{ScoreSignatures_UCell()}, or
#'   \code{AUCell_calcAUC()}).
#'
#' @return A numeric matrix with one row per cell and one column per gene set,
#'   ordered as in \code{gene.sets}.
#'
#' @author Nick Borcherding, Jared Andrews
#'
#' @seealso \code{\link{runEscape}} to attach scores to a single-cell object;
#'   \code{\link{getGeneSets}} for MSigDB retrieval;
#'   \code{\link{performNormalization}} for the optional normalization workflow.
#'
#' @examples
#' gs <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' 
#' pbmc <- SeuratObject::pbmc_small
#' es   <- escape.matrix(pbmc, 
#'                       gene.sets = gs,
#'                       method = "ssGSEA", 
#'                       groups = 500, 
#'                       min.size = 3)
#'
#' @export
escape.matrix <- function(input.data,
                          gene.sets        = NULL,
                          method           = "ssGSEA",
                          groups           = 1000,
                          min.size         = 5,
                          normalize        = FALSE,
                          make.positive    = FALSE,
                          min.expr.cells   = 0,
                          min.filter.by    = NULL,
                          BPPARAM          = NULL,
                          ...) {
  if(is.null(min.size)) min.size <- 0
  
  # ---- 1) resolve gene-sets & counts ----------------------------------------
  egc  <- .GS.check(gene.sets)
  cnts <- .cntEval(input.data, assay = "RNA", type = "counts")  # dgCMatrix
  
  if (is.null(min.filter.by)) {
    cnts <- .filter_genes(cnts, min.expr.cells)
  } else {
    # get grouping factor from object
    group.vec <- .extract_group_vector(input.data, min.filter.by)
    split.idx <- split(seq_len(ncol(cnts)), group.vec)
    
    cnts <- do.call(cbind, lapply(split.idx, function(cols) {
      sub <- cnts[, cols, drop = FALSE]
      .filter_genes(sub, min.expr.cells)
    }))
  }
  
  # ---- 2) drop undersized gene-sets -----------------------------------------
  keep <- vapply(egc, function(gs) sum(rownames(cnts) %in% gs) >= min.size,
                 logical(1))
  if (!all(keep)) {
    egc <- egc[keep]
    if (!length(egc))
      stop("No gene-sets meet the size threshold (min.size = ", min.size, ")")
  }
  
  # ---- 3) split cells into chunks -------------------------------------------
  chunks <- .split_cols(cnts, groups)
  message("escape.matrix(): processing ", length(chunks), " chunk(s)...")
  
  # ---- 4) compute enrichment in parallel ------------------------------------
  res_list <- .plapply(
    chunks,
    function(mat)
      .compute_enrichment(mat, egc, method, BPPARAM, ...),
    BPPARAM  = BPPARAM
  )
  
  # ---- 5) combine + orient (rows = cells) -----------------------------------
  all_sets <- names(egc)
  res_mat  <- do.call(cbind, lapply(res_list, function(m) {
    m <- as.matrix(m)
    m <- m[match(all_sets, rownames(m)), , drop = FALSE]
    m
  }))
  res_mat <- t(res_mat)
  colnames(res_mat) <- all_sets
  
  # ---- 6) optional dropout scaling ------------------------------------------
  if (normalize) {
    res_mat <- performNormalization(
      input.data      = input.data,
      enrichment.data = res_mat,
      assay           = NULL,
      gene.sets       = gene.sets,
      make.positive   = make.positive,
      groups          = groups
    )
    if (.is_seurat_or_sce(input.data)) {
      res_mat <- .pull.Enrich(res_mat, "escape_normalized")
    }
  }
  
  res_mat
}

#' Calculate Enrichment Scores Using Seurat or SingleCellExperiment Objects
#'
#' \code{runEscape()} is a convenience wrapper around \code{\link{escape.matrix}}
#' that computes enrichment scores and inserts them as a new assay (default
#' \code{"escape"}) in a \pkg{Seurat} or \pkg{SingleCellExperiment} object. All
#' arguments (except \code{new.assay.name}) map directly to their counterparts
#' in \code{escape.matrix()}.
#'
#' @inheritParams escape.matrix
#' @param new.assay.name Character. Name for the assay that will store the
#'   enrichment matrix in the returned object. Default is \code{"escape"}.
#'
#' @return The input single-cell object with an additional assay containing the
#'   enrichment scores (cells x gene-sets). Matrix orientation follows standard
#'   single-cell conventions (gene-sets as rows inside the assay).
#'
#' @author Nick Borcherding, Jared Andrews
#'
#' @seealso \code{\link{escape.matrix}} for the underlying computation;
#'   \code{\link{performNormalization}} to add normalized scores;
#'   \code{\link{heatmapEnrichment}}, \code{\link{ridgeEnrichment}}, and
#'   related plotting helpers for visualization.
#'
#' @examples
#' gs <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' 
#' sce <- SeuratObject::pbmc_small
#' sce <- runEscape(sce, 
#'                  gene.sets = gs, 
#'                  method = "GSVA",
#'                  groups = 1000, 
#'                  min.size = 3,
#'                  new.assay.name = "escape")
#'
#' @export
runEscape <- function(input.data,
                      gene.sets,
                      method = c("ssGSEA", "GSVA", "UCell", "AUCell"),
                      groups = 1e3,
                      min.size = 5,
                      normalize = FALSE,
                      make.positive = FALSE,
                      new.assay.name = "escape",
                      min.expr.cells   = 0,
                      min.filter.by    = NULL,
                      BPPARAM = NULL,
                      ...) {
    method <- match.arg(method)
    .checkSingleObject(input.data)
    esc <- escape.matrix(input.data, gene.sets, method, groups, min.size,
                         normalize, make.positive, min.expr.cells, 
                         min.filter.by, BPPARAM, ...)
    
    input.data <- .adding.Enrich(input.data, esc, new.assay.name)
    return(input.data)
}


.filter_genes <- function(m, min.expr.cells) {
  if (is.null(min.expr.cells) || identical(min.expr.cells, 0))
    return(m)                        # nothing to do
  
  ncells <- ncol(m)
  
  thr <- if (min.expr.cells < 1)
    ceiling(min.expr.cells * ncells)  # proportion → absolute
  else
    as.integer(min.expr.cells)
  
  keep <- Matrix::rowSums(m > 0) >= thr
  m[keep, , drop = FALSE]
}

# helper: pull a column from meta.data / colData no matter the object
#' @importFrom SummarizedExperiment colData
.extract_group_vector <- function(obj, col) {
  if (.is_seurat(obj))
    return(obj[[col, drop = TRUE]])
  if (.is_sce(obj))
    return(colData(obj)[[col]])
  stop("min.filter.by requires a Seurat or SingleCellExperiment object")
}

.filter_genes <- function(m, min.expr.cells) {
  if (is.null(min.expr.cells) || identical(min.expr.cells, 0))
    return(m)                        # nothing to do
  
  ncells <- ncol(m)
  
  thr <- if (min.expr.cells < 1)
    ceiling(min.expr.cells * ncells)  # proportion → absolute
  else
    as.integer(min.expr.cells)
  
  keep <- Matrix::rowSums(m > 0) >= thr
  m[keep, , drop = FALSE]
}
