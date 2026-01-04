#' Visualize Enrichment Distributions Using Geyser Plots
#'
#' This function allows the user to examine the distribution of
#' enrichment across groups by generating a geyser plot.
#'
#' @param input.data Output of \code{\link{escape.matrix}} or a single-cell
#'   object previously processed by \code{\link{runEscape}}.
#' @param assay Character. Name of the assay holding enrichment scores when
#'   \code{input.data} is a single-cell object. Ignored otherwise.
#' @param group.by Character. Metadata column plotted on the x-axis. Defaults
#'   to the Seurat/SCE \code{ident} slot when \code{NULL}.
#' @param gene.set.use Character. Name of the gene set to display.
#' @param color.by Character. Aesthetic mapped to point color. Options:
#'   \itemize{
#'     \item \code{"group"} (default): Uses \code{group.by} for categorical
#'       coloring.
#'     \item \emph{gene-set name}: Use the same value as \code{gene.set.use}
#'       to obtain a numeric gradient.
#'     \item Any other metadata column present in the data.
#'   }
#' @param order.by Character or \code{NULL}. How to arrange the x-axis:
#'   \itemize{
#'     \item \code{"mean"}: Groups ordered by decreasing group mean.
#'     \item \code{"group"}: Natural (alphanumeric) sort of group labels.
#'     \item \code{NULL} (default): Keep original ordering.
#'   }
#' @param facet.by Character or \code{NULL}. Metadata column used to facet
#'   the plot.
#' @param summarise.by Character or \code{NULL}. Metadata column used to
#'   summarise data before plotting.
#' @param summary.stat Character. Method used to summarize expression within
#'   each group defined by \code{summarise.by}. One of: \code{"mean"}
#'   (default), \code{"median"}, \code{"max"}, \code{"sum"}, or
#'   \code{"geometric"}.
#' @param scale Logical. If \code{TRUE}, scores are centered and scaled
#'   (Z-score) prior to plotting. Default is \code{FALSE}.
#' @param palette Character. Color palette name from
#'   \code{\link[grDevices]{hcl.pals}}. Default is \code{"inferno"}.
#'
#' @examples
#' gs <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#'
#' pbmc <- SeuratObject::pbmc_small |>
#'   runEscape(gene.sets = gs,
#'             min.size = NULL)
#'
#' geyserEnrichment(pbmc,
#'                  assay = "escape",
#'                  gene.set.use = "Tcells")
#'
#' @import ggplot2
#' @importFrom ggdist stat_pointinterval
#' @importFrom stats as.formula
#' @return A \pkg{ggplot2} object.
#' @export
geyserEnrichment <- function(input.data,
                             assay        = NULL,
                             group.by     = NULL,
                             gene.set.use,
                             color.by     = "group",
                             order.by     = NULL,
                             scale        = FALSE,
                             facet.by     = NULL,
                             summarise.by = NULL,
                             summary.stat = "mean",
                             palette      = "inferno") {
  ## ---- 0) Sanity checks -----------------------------------------------------
  if (missing(gene.set.use) || length(gene.set.use) != 1L)
    stop("Please supply exactly one 'gene.set.use' to plot.")

  if (is.null(group.by))
    group.by <- "ident"

  if (identical(color.by, "group"))
    color.by <- group.by

  if (!is.null(summarise.by) && (identical(summarise.by, group.by) ||
      identical(summarise.by, facet.by)))
    stop("'summarise.by' cannot be the same as 'group.by' or 'facet.by'.
         Please choose a different metadata column.")

  # ---- 1) helper to match summary function -------------------------
  summary_fun <- .match_summary_fun(summary.stat)

  ## ---- 2) Build tidy data.frame -------------------------------------------
  enriched <- .prepData(input.data, assay, gene.set.use, group.by,
                        split.by = summarise.by, facet.by = facet.by, color.by = color.by)

  # Define all grouping variables that must be metadata columns
  grouping_vars <- unique(c(summarise.by, group.by, facet.by))

  # Determine if color.by is a feature
  all_features <- rownames(.cntEval(input.data, assay = assay, type = "data"))

  # Determine if color.by is a feature
  is_feature_color <- !is.null(color.by) &&
    (color.by %in% all_features)

  ## Optionally summarize data with **base aggregate()** ----------------------
  if (!is.null(summarise.by)) {

    # add color.by to summarise_vars if it is a feature, otherwise add to grouping_vars
    summarise_vars <- unique(c(gene.set.use, if (is_feature_color) color.by))
    grouping_vars <- unique(c(grouping_vars, if (!is_feature_color) color.by))

    # Perform aggregation
    enriched <- aggregate(enriched[summarise_vars],
                          by = enriched[grouping_vars],
                          FUN = summary_fun,
                          simplify = TRUE)
  }

  ## Optionally Z-transform ----------------------------------------------------
  if (scale) {
    enriched[[gene.set.use]] <- scale(as.numeric(enriched[[gene.set.use]]))

    # Also scale color.by if it's a feature
    if (is_feature_color) {
      enriched[[color.by]] <- scale(enriched[[color.by]])
    }
  }

  ## Optionally reorder groups -------------------------------------------------
  if (!is.null(order.by))
    enriched <- .orderFunction(enriched, order.by, group.by)

  ## ---- 3) Plot --------------------------------------------------------------
  if (!is.null(color.by))
    plt <- ggplot(enriched, aes(x = .data[[group.by]],
                                y = .data[[gene.set.use]],
                                group = .data[[group.by]],
                                colour = .data[[color.by]]))
  else
    plt <- ggplot(enriched, aes(x = .data[[group.by]],
                                y = .data[[gene.set.use]]),
                                group = .data[[group.by]])

    # Raw points --------------------------------------------------------------
  plt <- plt + geom_jitter(width = 0.25, size = 1.5, alpha = 0.6, na.rm = TRUE) +

    # White base interval + median point -------------------------------------
  stat_pointinterval(interval_size_range = c(2, 3), fatten_point = 1.4,
                     interval_colour = "white", point_colour = "white",
                     position = position_dodge(width = 0.6), show.legend = FALSE) +

    # Black outline for clarity ----------------------------------------------
  stat_pointinterval(interval_size_range = c(1, 2), fatten_point = 1.4,
                     interval_colour = "black", point_colour = "black",
                     position = position_dodge(width = 0.6), show.legend = FALSE) +

    labs(x        = group.by,
         y        = paste0(gene.set.use, "\nEnrichment Score"),
         colour   = color.by) +
    .themeEscape(grid_lines = "Y", legend_position = "bottom") +
    ggplot2::theme(legend.direction = "horizontal")

  ## ---- 4) Colour scale ------------------------------------------------------
  if (!is.null(color.by))
    plt <- .colorby(enriched, plt, color.by, palette, type = "color")

  ## ---- 5) Facetting ---------------------------------------------------------
  if (!is.null(facet.by))
    plt <- plt + facet_grid(as.formula(paste(".~", facet.by)))

  plt
}
