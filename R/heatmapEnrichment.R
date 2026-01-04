#' Visualize Enrichment Value Summaries Using Heatmaps
#'
#' This function allows the user to examine a heatmap with the mean enrichment
#' values by group. The heatmap displays gene sets as rows and the grouping
#' variable as columns.
#'
#' @param input.data Output of \code{\link{escape.matrix}} or a single-cell
#'   object previously processed by \code{\link{runEscape}}.
#' @param assay Character. Name of the assay holding enrichment scores when
#'   \code{input.data} is a single-cell object. Ignored otherwise.
#' @param group.by Character. Metadata column plotted on the x-axis. Defaults
#'   to the Seurat/SCE \code{ident} slot when \code{NULL}.
#' @param gene.set.use Character vector or \code{"all"}. Gene-set names to
#'   plot. Use \code{"all"} (default) to show every available gene set.
#' @param cluster.rows Logical. If \code{TRUE}, rows are ordered by Ward-linkage
#'   hierarchical clustering (Euclidean distance). Default is \code{FALSE}.
#' @param cluster.columns Logical. If \code{TRUE}, columns are ordered by
#'   Ward-linkage hierarchical clustering (Euclidean distance). Default is
#'   \code{FALSE}.
#' @param facet.by Character or \code{NULL}. Metadata column used to facet
#'   the plot.
#' @param scale Logical. If \code{TRUE}, Z-transforms each gene-set column
#'   \strong{after} summarization. Default is \code{FALSE}.
#' @param summary.stat Character. Method used to summarize expression within
#'   each group. One of: \code{"mean"} (default), \code{"median"},
#'   \code{"max"}, \code{"sum"}, or \code{"geometric"}.
#' @param palette Character. Color palette name from
#'   \code{\link[grDevices]{hcl.pals}}. Default is \code{"inferno"}.
#'
#' @return A \code{ggplot2} object.
#' @importFrom stats aggregate dist hclust
#' @export
#'
#' @examples
#' gs <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#'            
#' pbmc <- SeuratObject::pbmc_small |>
#'   runEscape(gene.sets = gs, min.size = NULL)
#'   
#' heatmapEnrichment(pbmc, assay = "escape", palette = "viridis")
#' 
heatmapEnrichment <- function(input.data,
                              assay          = NULL,
                              group.by       = NULL,
                              gene.set.use   = "all",
                              cluster.rows   = FALSE,
                              cluster.columns= FALSE,
                              facet.by       = NULL,
                              scale          = FALSE,
                              summary.stat   = "mean",
                              palette        = "inferno")
{
  # ---------- 1. helper to match summary function -------------------------
  summary_fun <- .match_summary_fun(summary.stat)
  
  # ---------- 2. pull / tidy data -----------------------------------------
  if (is.null(group.by)) group.by <- "ident"
  df <- .prepData(input.data, assay, gene.set.use,
                  group.by = group.by,
                  split.by = NULL,
                  facet.by = facet.by, 
                  color.by = NULL)
  
  # Which columns contain gene-set scores?
  if (identical(gene.set.use, "all"))
    gene.set <- setdiff(colnames(df), c(group.by, facet.by))
  else
    gene.set <- gene.set.use
  if (!length(gene.set))
    stop("No gene-set columns found to plot.")
  
  # ---------- 3. summarise with **base aggregate()** ----------------------
  grp_cols   <- c(group.by, facet.by)               # one or two columns
  agg <- aggregate(df[gene.set],
                   by   = df[grp_cols],
                   FUN  = summary_fun,
                   SIMPLIFY = FALSE)
  # aggregate() keeps grouping columns first; ensure correct names
  names(agg)[seq_along(grp_cols)] <- grp_cols
  
  # Optional Z-transform AFTER summary
  if (scale)
    agg[gene.set] <- lapply(agg[gene.set], scale)
  
  # ---------- 4. long format for ggplot (base-R) --------------------------
  long <- data.frame(
    variable = rep(gene.set, each = nrow(agg)),
    value    = unlist(agg[gene.set], use.names = FALSE),
    group    = rep(agg[[group.by]], times = length(gene.set)),
    stringsAsFactors = FALSE
  )
  if (!is.null(facet.by))
    long[[facet.by]] <- rep(agg[[facet.by]], times = length(gene.set))
  
  # ---------- 5. optional clustering --------------------------------------
  if (cluster.rows) {
    ord <- hclust(dist(t(agg[gene.set])), method = "ward.D2")$order
    long$variable <- factor(long$variable, levels = gene.set[ord])
  }
  if (cluster.columns) {
    ord <- hclust(dist(agg[gene.set]), method = "ward.D2")$order
    long$group <- factor(long$group, levels = agg[[group.by]][ord])
  }
  
  # ---------- 6. draw ------------------------------------------------------
  p <- ggplot2::ggplot(long,
                       ggplot2::aes(x = group, y = variable, fill = value)) +
    ggplot2::geom_tile(colour = "black", linewidth = 0.4) +
    ggplot2::scale_fill_gradientn(colours = .colorizer(palette, 11),
                                  name = "Enrichment") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::coord_equal() +
    .themeEscape(grid_lines = "none", legend_position = "bottom") +
    ggplot2::theme(axis.title       = ggplot2::element_blank(),
                   axis.ticks       = ggplot2::element_blank(),
                   legend.direction = "horizontal",
                   panel.border     = ggplot2::element_blank())
  
  if (!is.null(facet.by))
    p <- p + ggplot2::facet_grid(stats::as.formula(paste(". ~", facet.by)))
  
  p
}
