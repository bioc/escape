#' Visualize the PCA of Enrichment Values
#'
#' This function allows the user to examine the distribution of principal
#' components computed on the enrichment values.
#'
#' @param input.data Single-cell object (Seurat / SCE) \strong{or} the raw list
#'   returned by \code{\link{performPCA}}.
#' @param dimRed Character. Name of the dimensional-reduction slot to pull from
#'   a single-cell object. Ignored when \code{input.data} is the list output.
#' @param x.axis Character. Name of the principal component to display on the
#'   x-axis (e.g., \code{"PC1"}). Default is \code{"PC1"}.
#' @param y.axis Character. Name of the principal component to display on the
#'   y-axis (e.g., \code{"PC2"}). Default is \code{"PC2"}.
#' @param facet.by Character or \code{NULL}. Metadata column used to facet
#'   the plot.
#' @param style Character. Plot style. Options:
#'   \itemize{
#'     \item \code{"point"} (default): Density-aware scatter plot.
#'     \item \code{"hex"}: Hexagonal binning.
#'   }
#' @param add.percent.contribution Logical. If \code{TRUE}, include percent
#'   variance explained in axis labels. Default is \code{TRUE}.
#' @param display.factors Logical. If \code{TRUE}, draw arrows for the top
#'   gene-set loadings. Default is \code{FALSE}.
#' @param number.of.factors Integer. Number of loadings to display when
#'   \code{display.factors = TRUE}. Default is \code{10}.
#' @param palette Character. Color palette name from
#'   \code{\link[grDevices]{hcl.pals}}. Default is \code{"inferno"}.
#'
#' @examples 
#' GS <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' pbmc_small <- SeuratObject::pbmc_small
#' pbmc_small <- runEscape(pbmc_small, 
#'                         gene.sets = GS, 
#'                         min.size = NULL)
#'                         
#' pbmc_small <- performPCA(pbmc_small, 
#'                          assay = "escape")
#'                          
#' pcaEnrichment(pbmc_small,
#'               x.axis = "PC1",
#'               y.axis = "PC2",
#'               dimRed = "escape.PCA")
#'
#' @return A **ggplot2** object.
#' @importFrom utils head
#' @export
pcaEnrichment <- function(input.data,
                          dimRed = NULL,
                          x.axis = "PC1",
                          y.axis = "PC2",
                          facet.by = NULL,
                          style = c("point", "hex"),
                          add.percent.contribution = TRUE,
                          display.factors = FALSE,
                          number.of.factors = 10,
                          palette = "inferno") {
  
  style <- match.arg(style)
  
  # ------------------------------------------------------------------------
  # 1. Extract PCA slots ----------------------------------------------------
  # ------------------------------------------------------------------------
  if (.is_seurat_or_sce(input.data)) {
    pca.values <- .grabDimRed(input.data, dimRed)
  } else if (is.list(input.data) && length(input.data) == 4) {
    pca.values <- input.data
    if (!is.null(facet.by))
      stop("facet.by is only valid with a single-cell object.")
  } else {
    stop("'input.data' must be a Seurat / SCE object or the list from performPCA().")
  }
  
  # Helper to convert PCX to X
  pc_idx <- function(pc) as.integer(sub("PC", "", pc, ignore.case = TRUE))
  x.idx <- pc_idx(x.axis)
  y.idx <- pc_idx(y.axis)
  
  # Axis labels with % variance ------------------------------------------------
  if (add.percent.contribution && "contribution" %in% names(pca.values)) {
    pc.var <- pca.values$contribution
    x.title <- sprintf("%s (%.1f%%)", x.axis, pc.var[x.idx])
    y.title <- sprintf("%s (%.1f%%)", y.axis, pc.var[y.idx])
  } else {
    x.title <- x.axis
    y.title <- y.axis
  }
  
  # ------------------------------------------------------------------------
  # 2. Build plotting data.frame -------------------------------------------
  # ------------------------------------------------------------------------
  plot.df <- as.data.frame(pca.values[[1]])
  
  if (!is.null(facet.by)) {
    meta <- .grabMeta(input.data)
    if (!facet.by %in% colnames(meta))
      stop("'", facet.by, "' not found in the single-cell object metadata.")
    plot.df[[facet.by]] <- meta[[facet.by]]
  }
  
  # ------------------------------------------------------------------------
  # 3. Base ggplot ----------------------------------------------------------
  # ------------------------------------------------------------------------
  aes.map <- ggplot2::aes(x = plot.df[,x.idx], y = plot.df[,y.idx])
  g <- ggplot2::ggplot(plot.df, aes.map)
  
  if (style == "point") {
    if (!requireNamespace("ggpointdensity", quietly = TRUE)) {
      warning("Package `ggpointdensity` not installed, falling back to alpha-blended points")
      g <- g + ggplot2::geom_point(alpha = 0.4, size = 0.6)
    } else {
      g <- g + ggpointdensity::geom_pointdensity() +
        ggplot2::scale_color_gradientn(colors = grDevices::hcl.colors(11, palette)) +
        ggplot2::labs(color = "Density")
    }
  } else {                                    # hex-bin
    if (!requireNamespace("hexbin", quietly = TRUE))
      stop("'hexbin' package required for style = 'hex'.")
    g <- g + ggplot2::stat_binhex() +
      ggplot2::scale_fill_gradientn(colors = grDevices::hcl.colors(11, palette)) +
      ggplot2::labs(fill = "Count")
  }
  
  g <- g + ggplot2::labs(x = x.title, y = y.title) + .themeEscape(grid_lines = "none")
  
  if (!is.null(facet.by))
    g <- g + ggplot2::facet_grid(stats::as.formula(paste(".~", facet.by)))
  
  # ------------------------------------------------------------------------
  # 4. Biplot arrows --------------------------------------------------------
  # ------------------------------------------------------------------------
  if (display.factors) {
    loadings <- as.data.frame(pca.values$rotation)
    sel.score <- (loadings[[x.idx]]^2 + loadings[[y.idx]]^2) / 2
    sel <- head(order(sel.score, decreasing = TRUE), number.of.factors)
    loadings <- loadings[sel, ]
    loadings$names <- rownames(loadings)
    
    # Rescale onto existing plot range (80 % of extents)
    rng.x <- range(plot.df[[x.idx]]) * 0.8
    rng.y <- range(plot.df[[y.idx]]) * 0.8
    rescale <- function(v, to) (v - min(v)) / diff(range(v)) * diff(to) + min(to)
    loadings$xend <- rescale(loadings[[x.idx]], rng.x)
    loadings$yend <- rescale(loadings[[y.idx]], rng.y)
    
    g <- g +
      ggplot2::geom_hline(yintercept = 0, linetype = 2) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) +
      ggplot2::geom_segment(data = loadings,
                            ggplot2::aes(x = 0, y = 0, xend = xend, yend = yend),
                            arrow = ggplot2::arrow(length = grid::unit(0.25, "cm"))) +
      ggplot2::geom_label(data = loadings,
                         ggplot2::aes(x = xend, y = yend, label = names),
                         size = 2)
  }
  
  g
}
