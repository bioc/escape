#' Visualize Mean Density Ranking of Genes Across Gene Sets
#'
#' This function allows the user to examine the mean ranking within groups
#' across the gene set. The visualization uses the density function to display
#' the relative position and distribution of rank.
#'
#' @param input.data A \link[SeuratObject]{Seurat} object or a
#'   \link[SingleCellExperiment]{SingleCellExperiment}.
#' @param gene.set.use Character. Name of the gene set to display.
#' @param gene.sets A named list of character vectors, the result of
#'   \code{\link{getGeneSets}}, or the built-in data object
#'   \code{\link{escape.gene.sets}}.
#' @param group.by Character. Metadata column used for grouping. Defaults to
#'   the Seurat/SCE \code{ident} slot when \code{NULL}.
#' @param rug.height Numeric. Vertical spacing of the hit rug as a fraction of
#'   the y-axis. Default is \code{0.02}.
#' @param palette Character. Color palette name from
#'   \code{\link[grDevices]{hcl.pals}}. Default is \code{"inferno"}.
#'
#' @examples
#' gs <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#'            
#' pbmc_small <- SeuratObject::pbmc_small
#'                         
#' densityEnrichment(pbmc_small, 
#'                   gene.set.use = "Tcells",
#'                   gene.sets = gs)
#'
#' @return A `patchwork`/`ggplot2` object.
#' @export
#'
#' @import ggplot2
#' @importFrom stats na.omit
#' @importFrom MatrixGenerics rowMeans2
densityEnrichment <- function(input.data,
                              gene.set.use,
                              gene.sets,
                              group.by = NULL,
                              rug.height  = 0.02,
                              palette  = "inferno") {
  ## -------- 0  Input checks --------------------------------------------------
  .checkSingleObject(input.data)
  if (is.null(group.by)) group.by <- "ident"
  
  gene.sets <- .GS.check(gene.sets)
  if (!gene.set.use %in% names(gene.sets))
    stop("'gene.set.use' not found in 'gene.sets'")
  
  ## -------- 1  Counts & grouping --------------------------------------------
  cnts <- .cntEval(input.data, assay = "RNA", type = "counts") |>
    .filterFeatures()
  
  meta   <- .grabMeta(input.data)
  groups <- na.omit(unique(meta[[group.by]]))
  
  ## -------- 2  Fast rank computation per group ------------------------------
  n.genes <- nrow(cnts)
  weights <- abs(seq(n.genes, 1) - n.genes/2)       # fixed triangular weight
  rank.mat <- matrix(NA_integer_, n.genes, length(groups),
                     dimnames = list(rownames(cnts),
                                     paste0(group.by, ".", groups)))
  
  compute.cdf <- utils::getFromNamespace("compute.gene.cdf", "GSVA")
  
  for (i in seq_along(groups)) {
    cols <- which(meta[[group.by]] == groups[i])
    tmp  <- cnts[, cols, drop = FALSE]
    
    dens <- suppressWarnings(
      compute.cdf(tmp, seq_len(ncol(tmp)), FALSE, FALSE)
    )
    ord  <- apply(dens, 2, order, decreasing = TRUE)          # genes x cells
    scores <- vapply(seq_len(ncol(ord)),
                     function(j) weights[ord[, j]],
                     numeric(n.genes))
    
    mean.score     <- rowMeans2(scores)
    rank.mat[, i]  <- round(rank(mean.score, ties.method = "average") / 2)
  }
  
  ## -------- 3  Long data.frame w/o extra deps -------------------------------
  in.set <- rownames(rank.mat) %in% gene.sets[[gene.set.use]]
  long.df <- data.frame(
    value          = as.vector(rank.mat),
    variable       = rep(colnames(rank.mat), each = n.genes),
    gene.set.query = rep(ifelse(in.set, "yes", NA_character_), times = length(groups)),
    stringsAsFactors = FALSE
  )
  
  ## -------- 4  Plots ---------------------------------------------------------
  cols <- .colorizer(palette, length(groups))

  # Filter to gene set members with valid values
  plot.df <- long.df[long.df$gene.set.query == "yes" & is.finite(long.df$value), ]

  # Alphanumerically sort group labels for consistent ordering
  plot.df$variable <- factor(plot.df$variable,
                             levels = .alphanumericalSort(unique(plot.df$variable)))

  # Density plot panel

  p1 <- ggplot(plot.df,
               aes(x = value, fill = variable)) +
    geom_density(alpha = 0.5, colour = "black", linewidth = 0.4) +
    scale_fill_manual(values = cols, name = "Group") +
    labs(y = "Rank Density",
         title = paste0("Gene Set: ", gene.set.use)) +
    .themeEscape(grid_lines = "Y") +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title   = element_text(hjust = 0.5))

  # Build rug segment data with proper stacking
  seg.df <- plot.df
  seg.df$ord  <- match(seg.df$variable, levels(seg.df$variable))
  seg.df$y    <- -(seg.df$ord - 1) * rug.height
  seg.df$yend <- seg.df$y - rug.height * 0.9

  # Rug plot panel
  p2 <- ggplot(seg.df, aes(x = value, xend = value,
                           y = y, yend = yend,
                           colour = variable)) +
    geom_segment(linewidth = 0.8, alpha = 0.7) +
    scale_colour_manual(values = cols, guide = "none") +
    labs(x = "Mean Rank Order") +
    .themeEscape(grid_lines = "none") +
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          panel.border = element_rect(fill = NA, colour = "grey50", linewidth = 0.5))

  patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(3, 1))
}