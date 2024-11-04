#' Generate a heatmap to visualize enrichment values
#' 
#' This function allows to the user to examine the heatmap with the mean
#' enrichment values by group. The heatmap will have the gene sets as rows
#' and columns will be the grouping variable.
#'
#' @param input.data Enrichment output from \code{\link{escape.matrix}} or
#' \code{\link{runEscape}}.
#' @param assay Name of the assay to plot if data is a single-cell object.
#' @param group.by Categorical parameter to plot along the x.axis. If input is
#' a single-cell object the default will be cluster.
#' @param gene.set.use Selected gene sets to visualize. If \strong{"all"}, the 
#' heatmap will be generated across all gene sets.
#' @param cluster.rows Use Euclidean distance to order the row values.
#' @param cluster.columns Use Euclidean distance to order the column values.
#' @param facet.by Variable to facet the plot into n distinct graphs.
#' @param scale Visualize raw values \strong{FALSE} or Z-transform 
#' enrichment values \strong{TRUE}.
#' @param summary.stat Use \strong{'median'} or \strong{'mean'} values
#' to display.
#' @param palette Colors to use in visualization - input any 
#' \link[grDevices]{hcl.pals}.
#' 
#' @import ggplot2
#' @importFrom stats dist hclust
#' @importFrom dplyr %>% group_by summarise across
#' @importFrom reshape2 melt
#' 
#' @examples
#' GS <- list(Bcells = c("MS4A1", "CD79B", "CD79A", "IGH1", "IGH2"),
#'            Tcells = c("CD3E", "CD3D", "CD3G", "CD7","CD8A"))
#' pbmc_small <- SeuratObject::pbmc_small
#' pbmc_small <- runEscape(pbmc_small, 
#'                         gene.sets = GS, 
#'                         min.size = NULL)
#'                         
#' heatmapEnrichment(pbmc_small, 
#'                   assay = "escape")
#'
#' @export
#'
#' @return ggplot2 object with heatmap of mean enrichment values

heatmapEnrichment <- function(input.data, 
                              assay = NULL,
                              group.by = NULL, 
                              gene.set.use = "all", 
                              cluster.rows = FALSE,
                              cluster.columns = FALSE,
                              facet.by = NULL,
                              scale = FALSE, 
                              summary.stat = "mean",
                              palette = "inferno") {
  
  options(dplyr.summarise.inform = FALSE)
  if(is.null(group.by)) {
    group.by <- "ident"
  }
  
  enriched <- .prepData(input.data, assay, gene.set.use, group.by, NULL, facet.by) 
  
  if(length(gene.set.use) == 1 && gene.set.use == "all") {
    gene.set <- colnames(enriched)[colnames(enriched) %!in% c(group.by, facet.by)]
  } else {
    gene.set <- gene.set.use
  }
  
  if(summary.stat %!in% c("median", "mean")) {
    message("Please select 'median' or 'mean' for the summary.stat argument. Using mean as a default")
    summary_func <- mean
  }
  # Select the appropriate summary function
  summary_func <- if (summary.stat == "median") median else mean
  
  if(!is.null(facet.by)) {
    enriched.summary <- enriched %>%
                          group_by(.data[[group.by]], .data[[facet.by]]) %>%
                          summarise(across(which(colnames(enriched) %in% gene.set), mean)) %>%
                          as.data.frame()
  } else {
    enriched.summary <- enriched %>%
                          group_by(.data[[group.by]]) %>%
                          summarise(across(which(colnames(enriched) %in% gene.set), mean)) %>%
                          as.data.frame()
  }
  
  if(scale) {
    enriched.summary[,gene.set] <- apply(enriched.summary[,gene.set], 2, scale)
  }
  
  reformated.enriched <-  suppressMessages(melt(enriched.summary))
  
  if(cluster.rows) {
    row.order <- gene.set[hclust(dist(t(enriched.summary[,gene.set]), method = "euclidean"), method = "ward.D2")$order]
    reformated.enriched[,"variable"] <- factor(reformated.enriched[,"variable"], levels = row.order)
  }
  
  if(cluster.columns) {
    column.order <- unique(enriched.summary[,group.by][hclust(dist(enriched.summary[,gene.set], method = "euclidean"), method = "ward.D2")$order])
    reformated.enriched[,group.by] <- factor(reformated.enriched[,group.by], levels = as.vector(column.order))
  }

   plot <- ggplot(reformated.enriched,
                  mapping = aes(x = reformated.enriched[,group.by],
                                 y = variable,
                                 fill = value)) +
                  geom_tile(color = "black", linewidth = 0.5) +
                  scale_y_discrete(expand = c(0, 0)) +
                  scale_x_discrete(expand = c(0, 0)) +
                  labs(fill = "Enrichment Score") +
                  guides(fill = guide_colorbar(title.position = "top",
                                                title.hjust = 0.5)) + 
                  coord_equal() + 
                  scale_fill_gradientn(colors = .colorizer(palette, 11)) +
                  theme_classic() + 
                  theme(axis.title = element_blank(),
                        axis.ticks = element_blank(), 
                        legend.direction = "horizontal", 
                        legend.position = "bottom")
    
    if (!is.null(facet.by)) {
      plot <- plot + 
                facet_grid(as.formula(paste('. ~', facet.by))) 
    }
    return(plot)
  
}
