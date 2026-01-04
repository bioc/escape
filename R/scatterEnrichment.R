#' Plot 2D Enrichment Distributions With Density or Hexplots
#'
#' Visualize the relationship between two enrichment scores at single-cell
#' resolution. By default, points are shaded by local 2-D density
#' (\code{color.by = "density"}), but users can instead color by a metadata
#' column (discrete) or by the raw gene-set scores themselves (continuous).
#'
#' @param input.data Output of \code{\link{escape.matrix}} or a single-cell
#'   object previously processed by \code{\link{runEscape}}.
#' @param assay Character. Name of the assay holding enrichment scores when
#'   \code{input.data} is a single-cell object. Ignored otherwise.
#' @param x.axis Character. Gene-set name to plot on the x-axis.
#' @param y.axis Character. Gene-set name to plot on the y-axis.
#' @param facet.by Character or \code{NULL}. Metadata column used to facet
#'   the plot.
#' @param group.by Character. Metadata column used when \code{color.by = "group"}.
#'   Defaults to the Seurat/SCE \code{ident} slot when \code{NULL}.
#' @param color.by Character. Aesthetic mapped to point color. Options:
#'   \itemize{
#'     \item \code{"density"} (default): Shade points by local 2-D density.
#'     \item \code{"group"}: Color by the \code{group.by} metadata column.
#'     \item \code{"x"}: Apply a continuous gradient based on the x-axis values.
#'     \item \code{"y"}: Apply a continuous gradient based on the y-axis values.
#'   }
#' @param style Character. Plot style. Options:
#'   \itemize{
#'     \item \code{"point"} (default): Density-aware scatter plot.
#'     \item \code{"hex"}: Hexagonal binning.
#'   }
#' @param scale Logical. If \code{TRUE}, scores are centered and scaled
#'   (Z-score) prior to plotting. Default is \code{FALSE}.
#' @param bins Integer. Number of hex bins along each axis when
#'   \code{style = "hex"}. Default is \code{40}.
#' @param point.size Numeric. Point size for \code{style = "point"}.
#'   Default is \code{1.2}.
#' @param alpha Numeric. Transparency for points or hexbins.
#'   Default is \code{0.8}.
#' @param add.corr Logical. If \code{TRUE}, add Pearson and Spearman
#'   correlation coefficients to the plot (top-left corner). Default is
#'   \code{FALSE}.
#' @param palette Character. Color palette name from
#'   \code{\link[grDevices]{hcl.pals}}. Default is \code{"inferno"}.
#' 
#' @examples
#' gs <- list(
#'   Bcells = c("MS4A1","CD79B","CD79A","IGH1","IGH2"),
#'   Tcells = c("CD3E","CD3D","CD3G","CD7","CD8A")
#' )
#' pbmc <- SeuratObject::pbmc_small |>
#'   runEscape(gene.sets = gs, min.size = NULL)
#'
#' scatterEnrichment(
#'   pbmc,
#'   assay     = "escape",
#'   x.axis    = "Tcells",
#'   y.axis    = "Bcells",
#'   color.by  = "group",        
#'   group.by  = "groups",
#'   add.corr  = TRUE,
#'   point.size = 1
#' )
#'
#' @return A \pkg{ggplot2} object.
#' @importFrom stats as.formula
#' @export
scatterEnrichment <- function(input.data,
                              assay      = NULL,
                              x.axis,
                              y.axis,
                              facet.by   = NULL,
                              group.by   = NULL,
                              color.by   = c("density", "group", "x", "y"),
                              style      = c("point", "hex"),
                              scale      = FALSE,
                              bins       = 40,
                              point.size = 1.2,
                              alpha      = 0.8,
                              palette    = "inferno",
                              add.corr   = FALSE) {
  
  ## ---- 0  Argument sanity checks -------------------------------------------
  style <- match.arg(style, choices = c("point", "hex"))
  color.by <- match.arg(color.by, choices = c("density", "group", "x", "y"))
  if (is.null(group.by)) group.by <- "ident"
  gene.set <- c(x.axis, y.axis)
  
  ## ---- 1  Assemble long data-frame -----------------------------------------
  enriched <- .prepData(input.data, assay, gene.set, group.by, NULL, facet.by,
                        color.by = NULL)
  
  if (scale) {
    enriched[, gene.set] <- apply(enriched[, gene.set, drop = FALSE], 2, scale)
  }
  
  ## ---- 2  Base ggplot2 object ----------------------------------------------
  aes_base <- ggplot2::aes(x = .data[[x.axis]], y = .data[[y.axis]])
  
  ## ---- 3  Choose colouring strategy ----------------------------------------
  
  if (color.by == "density") {
    aes_combined <- aes_base  # no color aesthetic
  } else if (color.by == "group") {
    aes_combined <- ggplot2::aes(
      x = .data[[x.axis]], 
      y = .data[[y.axis]], 
      color = .data[[group.by]]
    )
  } else {  # "x" or "y"
    sel <- if (color.by == "x") x.axis else y.axis
    aes_combined <- ggplot2::aes(
      x = .data[[x.axis]], 
      y = .data[[y.axis]], 
      color = .data[[sel]]
    )
  }
  
  # Now build the plot
  plt <- ggplot2::ggplot(enriched, aes_combined)
  
  ## ---- 4  Geometry ---------------------------------------------------------
  if (style == "point") {
    if (color.by == "density") {
      plt <- plt +
        ggpointdensity::geom_pointdensity(size = point.size, alpha = alpha) +
        ggplot2::scale_color_gradientn(
          colors = .colorizer(palette, 11),
          name   = "Local density")
    } else {
      geom <- ggplot2::geom_point(size = point.size, alpha = alpha)
      plt  <- plt + geom
    }
  } else {                                # hex-bin
    plt <- plt +
      ggplot2::stat_binhex(bins = bins, alpha = alpha) +
      ggplot2::scale_fill_gradientn(
        colors = .colorizer(palette, 11),
        name   = "Cells / bin")
  }
  
  ## ---- 5  Colour scaling for non-density modes -----------------------------
  if (color.by != "density") {
    sel <- switch(color.by,
                  group = group.by,
                  x     = x.axis,
                  y     = y.axis)
    
    plt <- .colorby(enriched, plt,
                    color.by = sel,
                    palette  = palette,
                    type     = "color")
  }
  
  ## ---- 6  Axes, theme, faceting -------------------------------------------
  plt <- plt +
    ggplot2::labs(x = paste0(x.axis, "\nEnrichment score"),
                  y = paste0(y.axis, "\nEnrichment score")) +
    .themeEscape(grid_lines = "none")
  
  if (!is.null(facet.by)) {
    plt <- plt + ggplot2::facet_grid(as.formula(paste(". ~", facet.by)))
  }
  
  ## ---- 7  Optional correlation overlay -------------------------------------
  if (add.corr) {
    cor_pears   <- stats::cor(enriched[[x.axis]], enriched[[y.axis]],
                              method = "pearson", use = "pairwise.complete.obs")
    cor_spear   <- stats::cor(enriched[[x.axis]], enriched[[y.axis]],
                              method = "spearman", use = "pairwise.complete.obs")
    lbl <- sprintf("Pearson rho = %.2f\nSpearman rho = %.2f", cor_pears, cor_spear)
    plt <- plt +
      ggplot2::annotate("text", x = -Inf, y = Inf, label = lbl,
                        hjust = 0, vjust = 1, size = 3.5,
                        fontface = "italic")
  }
  
  plt
}
