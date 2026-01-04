#' Visualize Enrichment Distributions Using Ridge Plots
#'
#' This function allows the user to examine the distribution of
#' enrichment across groups by generating a ridge plot.
#'
#' @param input.data Output of \code{\link{escape.matrix}} or a single-cell
#'   object previously processed by \code{\link{runEscape}}.
#' @param gene.set.use Character. Name of the gene set to display.
#' @param assay Character. Name of the assay holding enrichment scores when
#'   \code{input.data} is a single-cell object. Ignored otherwise.
#' @param group.by Character. Metadata column plotted on the y-axis. Defaults
#'   to the Seurat/SCE \code{ident} slot when \code{NULL}.
#' @param color.by Character. Aesthetic mapped to fill color. Options:
#'   \itemize{
#'     \item \code{"group"} (default): Uses \code{group.by} for categorical
#'       coloring.
#'     \item \emph{gene-set name}: Use the same value as \code{gene.set.use}
#'       to obtain a numeric gradient.
#'     \item Any other metadata column present in the data.
#'   }
#' @param order.by Character or \code{NULL}. How to arrange the y-axis:
#'   \itemize{
#'     \item \code{"mean"}: Groups ordered by decreasing group mean.
#'     \item \code{"group"}: Natural (alphanumeric) sort of group labels.
#'     \item \code{NULL} (default): Keep original ordering.
#'   }
#' @param facet.by Character or \code{NULL}. Metadata column used to facet
#'   the plot.
#' @param scale Logical. If \code{TRUE}, scores are centered and scaled
#'   (Z-score) prior to plotting. Default is \code{FALSE}.
#' @param add.rug Logical. If \code{TRUE}, draw per-cell tick marks underneath
#'   each ridge. Default is \code{FALSE}.
#' @param palette Character. Color palette name from
#'   \code{\link[grDevices]{hcl.pals}}. Default is \code{"inferno"}.
#'
#' @examples
#' gs <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#'            
#' pbmc <- SeuratObject::pbmc_small |>
#'   runEscape(gene.sets = gs, min.size = NULL)
#'
#' ridgeEnrichment(pbmc, assay = "escape",
#'                 gene.set.use = "Tcells",
#'                 group.by = "groups")
#'
#' @importFrom stats median
#' @return A [ggplot2] object.
#' @export
#'                 
ridgeEnrichment <- function(input.data,
                            gene.set.use,
                            assay      = NULL,
                            group.by   = NULL,
                            color.by   = "group",
                            order.by   = NULL,
                            scale      = FALSE,
                            facet.by   = NULL,
                            add.rug    = FALSE,
                            palette    = "inferno")
{
  ## ---- 0  sanity -------------------------------------------------------
  if (!requireNamespace("ggridges", quietly = TRUE))
    stop("Package 'ggridges' is required for ridge plots; please install it.")
  if (length(gene.set.use) != 1L)
    stop("'gene.set.use' must be length 1.")
  if (is.null(group.by)) group.by <- "ident"
  if (identical(color.by, "group")) color.by <- group.by
  
  ## ---- 1  build long data.frame ---------------------------------------
  df <- .prepData(input.data, assay, gene.set.use, group.by,
                  split.by = NULL, facet.by = facet.by, color.by = color.by)
  
  ## optional scaling (Z-transform per gene-set) -------------------------
  if (scale)
    df[[gene.set.use]] <- as.numeric(scale(df[[gene.set.use]], center = TRUE))
  
  ## optional re-ordering of the y-axis factor ---------------------------
  if (!is.null(order.by))
    df <- .orderFunction(df, order.by, group.by)
  
  ## detect “gradient” mode (numeric color mapped to x) -----------------
  gradient.mode <-
    is.numeric(df[[color.by]]) && identical(color.by, gene.set.use)
  
  if(gradient.mode) {
    fill <- ggplot2::after_stat(df[,color.by])
  } else {
    fill <- df[,color.by]
  }
  
  ## ---- 2  base ggplot --------------------------------------------------
  aes_base <- if (gradient.mode) {
    ggplot2::aes(
      x = .data[[gene.set.use]],
      y = .data[[group.by]],
      fill = after_stat(x)
    )
  } else {
    ggplot2::aes(
      x = .data[[gene.set.use]],
      y = .data[[group.by]],
      fill = .data[[color.by]]
    )
  }
  
  p <- ggplot2::ggplot(df, aes_base)
  
  ## choose ridge geometry + rug -----------------------------------------
  ridge_fun <- if (gradient.mode)
    ggridges::geom_density_ridges_gradient else ggridges::geom_density_ridges
  p <- p + do.call(ridge_fun, c(
    list(
      jittered_points = add.rug,
      point_shape     = '|',
      point_size      = 2.5,
      point_alpha     = 1,
      alpha           = 0.8,
      quantile_lines  = TRUE,
      quantile_fun    = median,
      vline_width     = 0.9
    ),
    if (add.rug) list(
      position = ggridges::position_points_jitter(width = 0.05, height = 0)
    )
  ))
  
  ## ---- 3  scales & labels ---------------------------------------------
  p <- p +
    ylab(group.by) +
    xlab(paste0(gene.set.use, "\nEnrichment Score")) +
    .themeEscape(grid_lines = "none")
  
  p <- .colorby(df, p, color.by, palette, type = "fill") + 
    guides(fill = "none")
  
  ## facetting ------------------------------------------------------------
  if (!is.null(facet.by))
    p <- p + ggplot2::facet_grid(stats::as.formula(paste(". ~", facet.by)))
  
  p
}