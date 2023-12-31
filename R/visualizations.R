# Assigning color palette
#' @importFrom grDevices colorRampPalette
assignColor <- function(x, enriched, group) {
    if (length(x) != length(unique(enriched[,group]))) {
        x <- colorRampPalette(x)(length(unique(enriched[,group])))
    } else { x <- x }
    return(x)
}


#' Density plot of the principal components
#'
#' @param PCAout The output of \code{\link{performPCA}}
#' @param PCx The principal component graphed on the x-axis
#' @param PCy The principal component graphed on the y-axis
#' @param colors The color palette for the density plot
#' @param contours Binary classifier to add contours to the density plot
#' @param facet A parameter to separate the graph
#'
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' 
#' @examples 
#' ES2 <- readRDS(url(
#' "https://ncborcherding.github.io/vignettes/escape_enrichment_results.rds"))
#' PCA <- performPCA(enriched = ES2, groups = c("Type", "Cluster"))
#' pcaEnrichment(PCA, PCx = "PC1", PCy = "PC2", contours = TRUE)
#'
#' @export
#'
#' @seealso \code{\link{performPCA}} for generating PCA results.
#' @return ggplot2 object of the results of PCA for the enrichment scores
pcaEnrichment <- function(PCAout, PCx, PCy, 
    colors = c("#0D0887FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF"), 
    contours = TRUE, facet = NULL) 
{
    plot <- ggplot(PCAout, aes(x=PCAout[,PCx], y=PCAout[,PCy])) +
        stat_binhex() +
        scale_fill_gradientn(colours = colorRampPalette(colors)(50)) +
        theme_classic() +
        ylab(PCy) +
        xlab(PCx)

    if (contours == TRUE) {
        plot <- plot + stat_density_2d(color = "black")
    }
    else {
        plot <- plot
    }

    if (!is.null(facet)) {
        plot <- plot + facet_wrap(as.formula(paste('~', facet)))
    }
    plot <- plot +
        geom_hline(yintercept = 0, lty=2) + 
        geom_vline(xintercept = 0, lty=2)

    return(plot)
}

#' Visualize the components of the PCA analysis of the enrichment results
#'
#' Graph the major gene set contributors to the \code{\link{pcaEnrichment}}.
#'
#' @param enriched The output of \code{\link{enrichIt}}.
#' @param gene.sets Names of gene sets to include in the PCA
#' @param PCx The principal component graphed on the x-axis.
#' @param PCy The principal component graphed on the y-axis.
#' @param top.contribution The number of gene sets to graph, organized 
#' by PCA contribution.
#'
#' @importFrom ggplot2 ggplot
#' @importFrom stats prcomp
#' @import dplyr
#' 
#' @examples 
#' ES2 <- readRDS(url(
#' "https://ncborcherding.github.io/vignettes/escape_enrichment_results.rds"))
#' 
#' masterPCAPlot(ES2, PCx = "PC1", PCy = "PC2", gene.sets = colnames(ES2), 
#' top.contribution = 10)
#'
#' @export
#'
#' @seealso \code{\link{enrichIt}} for generating enrichment scores.
#' @return ggplot2 object sumamrizing the PCA for the enrichment scores
masterPCAPlot <- function(enriched, gene.sets, PCx, PCy, top.contribution = 10) {
    input <- select_if(enriched, is.numeric)
    if (!is.null(gene.sets)) {
        input <- input[,colnames(input) %in% gene.sets]
    }
    PCA <- prcomp(input, scale. = TRUE)
    var_explained <- PCA$sdev^2/sum(PCA$sdev^2)
    
    tbl <- data.frame(names = rownames(PCA$rotation), 
            factors.y = PCA$rotation[,PCy]^2/sum(PCA$rotation[,PCy]^2),
            factors.x = PCA$rotation[,PCx]^2/sum(PCA$rotation[,PCx]^2)) 
    names <- tbl %>% top_n(n = 10, (factors.x + factors.y)/2)
    names <- names$names
    df <- as.data.frame(PCA$rotation)
    df <- df[rownames(df) %in% names,]
    df$names <- rownames(df)
    
    plot <- df %>%
        ggplot(aes(x=df[,PCx],y=df[,PCy])) + 
        geom_point() + 
        geom_text(aes_string(label = "names"), 
            size=2, hjust = 0.5, nudge_y = -0.01) + 
        geom_hline(yintercept = 0, lty=2) + 
        geom_vline(xintercept = 0, lty=2) +
        labs(x=paste0(PCx,": ",round(var_explained[1]*100,1),"%"),
            y=paste0(PCy, ": ",round(var_explained[2]*100,1),"%")) +
        theme_classic()
    return(plot)
}

#' Generate a ridge plot to examine enrichment distributions
#' 
#' This function allows to the user to examine the distribution of 
#' enrichment across groups by generating a ridge plot.
#'
#' @param enriched The output of \code{\link{enrichIt}}
#' @param group The parameter to group, displayed on the y-axis.
#' @param gene.set The gene set to graph on the x-axis. 
#' @param scale.bracket This will filter the enrichment scores to remove 
#' extreme outliers. Values entered (1 or 2 numbers) will be the filtering 
#' parameter using z-scores of the selected gene.set. If only 1 value is given, 
#' a seocndary bracket is autommatically selected as the inverse of the number.
#' @param colors The color palette for the ridge plot.
#' @param facet A parameter to separate the graph.
#' @param add.rug Binary classifier to add a rug plot to the x-axis.
#'
#' @import ggplot2
#' @importFrom ggridges geom_density_ridges geom_density_ridges2 position_points_jitter
#' 
#' @examples
#' ES2 <- readRDS(url(
#' "https://ncborcherding.github.io/vignettes/escape_enrichment_results.rds"))
#' ridgeEnrichment(ES2, gene.set = "HALLMARK_DNA_REPAIR", group = "cluster", 
#' facet = "Type", add.rug = TRUE)
#'
#' @export
#'
#' @seealso \code{\link{enrichIt}} for generating enrichment scores.
#' @return ggplot2 object with ridge-based distributions of selected gene.set
ridgeEnrichment <- function(enriched, group = "cluster", gene.set = NULL, 
            scale.bracket = NULL, facet = NULL, add.rug = FALSE,
            colors = c("#0D0887FF", "#7E03A8FF", "#CC4678FF", "#F89441FF", "#F0F921FF")) 
            {
    if (!is.null(scale.bracket)) {
        if (length(scale.bracket) != 1 | length(scale.bracket) != 1) {
            message("Please indicate one or two values for the scale.bracket 
                parameter, such as scale.bracket = c(-2,2)")
        }
        scale.bracket <- order(scale.bracket)
        if(length(scale.bracket) == 1) {
            scale.bracket <- c(scale.bracket, -scale.bracket)
            scale.bracket <- order(scale.bracket)
        } 
        tmp <- enriched
        tmp[,gene.set]<- scale(tmp[,gene.set])
        rows_selected <- rownames(tmp[tmp[,gene.set] >= scale.bracket[1] & 
                            tmp[,gene.set] <= scale.bracket[2],])
        enriched <- enriched[rownames(enriched) %in% rows_selected,]
    }
    colors <- assignColor(colors, enriched, group) 
    plot <- ggplot(enriched, aes(x = enriched[,gene.set], 
                    y = enriched[,group], fill = enriched[,group]))
    
    if (add.rug == TRUE) {
        plot <- plot + geom_density_ridges(
            jittered_points = TRUE,
            position = position_points_jitter(width = 0.05, height = 0),
            point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7) 
        
    } else {
        plot <- plot + 
            geom_density_ridges2(alpha = 0.8) 
    }
    
    plot <- plot + ylab(group) +
        xlab(paste0(gene.set, " (NES)")) +
        labs(fill = group) + 
        scale_fill_manual(values = colors) + 
        theme_classic() +
        guides(fill = "none")
    
    if (!is.null(facet)) {
        plot <- plot + facet_grid(as.formula(paste('. ~', facet))) }
    
    return(plot)
}

#Developing split violin plot
#Code from: https://stackoverflow.com/a/45614547
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
    draw_group = function(self, data, ..., draw_quantiles = NULL) {
        data <- transform(data, xminv = x - violinwidth * (x - xmin), 
            xmaxv = x + violinwidth * (xmax - x))
        grp <- data[1, "group"]
        newdata <- plyr::arrange(transform(data, x = 
            if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
        newdata <- rbind(newdata[1, ], 
            newdata, newdata[nrow(newdata), ], newdata[1, ])
        newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- 
            round(newdata[1, "x"])
        if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
            stopifnot(all(draw_quantiles >= 0), 
                all(draw_quantiles <= 1))
            quantiles <- 
                ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                    aesthetics <- data[rep(1, nrow(quantiles)), 
                    setdiff(names(data), c("x", "y")), drop = FALSE]
                    aesthetics$alpha <- rep(1, nrow(quantiles))
                both <- cbind(quantiles, aesthetics)
                quantile_grob <- GeomPath$draw_panel(both, ...)
                ggplot2:::ggname("geom_split_violin", 
                    grid::grobTree(GeomPolygon$draw_panel(newdata, ...), 
                            quantile_grob))
        } else {
            ggplot2:::ggname("geom_split_violin", 
                GeomPolygon$draw_panel(newdata, ...))}
})

#Defining new geometry
#Code from: https://stackoverflow.com/a/45614547
geom_split_violin <- 
    function(mapping = NULL, data = NULL, 
        stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, 
        trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, 
        inherit.aes = TRUE) {
    layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, 
        inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, 
        draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

#' Generate a split violin plot examine enrichment distributions
#' 
#' This function allows to the user to examine the distribution of 
#' enrichment across groups by generating a split violin plot.
#'
#' @param enriched The output of \code{\link{enrichIt}}
#' @param x.axis Optional parameter for seperation.
#' @param gene.set The gene set to graph on the y-axis. 
#' @param scale.bracket This will filter the enrichment scores to remove 
#' extreme outliers. Values entered (1 or 2 numbers) will be the filtering 
#' parameter using z-scores of the selected gene.set. If only 1 value is given, 
#' a secondary bracket is automatically selected as the inverse of the number.
#' @param split The parameter to split, must be binary.
#' @param colors The color palette for the ridge plot.
#'
#' @import ggplot2
#' 
#' @examples
#' ES2 <- readRDS(url(
#' "https://ncborcherding.github.io/vignettes/escape_enrichment_results.rds"))
#' splitEnrichment(ES2, x.axis = "cluster", split = "Type", 
#' gene.set = "HALLMARK_DNA_REPAIR")
#'
#' @export
#'
#' @seealso \code{\link{enrichIt}} for generating enrichment scores.
#' @return ggplot2 object violin-based distributions of selected gene.set
splitEnrichment <- function(enriched, x.axis = NULL, scale.bracket = NULL,
                            split = NULL, gene.set = NULL, 
                            colors = c("#0D0887FF", "#7E03A8FF", "#CC4678FF", 
                                       "#F89441FF", "#F0F921FF")) {
    
    if (length(unique(enriched[,split])) != 2) {
        message("SplitEnrichment() can only work for binary classification")}
    
    if (!is.null(scale.bracket)) {
        if (length(scale.bracket) != 1 | length(scale.bracket) != 1) {
            message("Please indicate one or two values for the scale.bracket 
                parameter, such as scale.bracket = c(-2,2)")
        }
        scale.bracket <- order(scale.bracket)
        if(length(scale.bracket) == 1) {
            scale.bracket <- c(scale.bracket, -scale.bracket)
            scale.bracket <- order(scale.bracket)
        } 
        tmp <- enriched
        tmp[,gene.set]<- scale(tmp[,gene.set])
        rows_selected <- rownames(tmp[tmp[,gene.set] >= scale.bracket[1] & 
                            tmp[,gene.set] <= scale.bracket[2],])
        enriched <- enriched[rownames(enriched) %in% rows_selected,]
    }
    colors <- assignColor(colors, enriched, split) 
    if (is.null(x.axis)) {
        plot <- ggplot(enriched, aes(x = ".", y = enriched[,gene.set], 
                    fill = enriched[,split])) 
        check = 1
    } else {
        plot <- ggplot(enriched, aes(x = enriched[,x.axis], 
                    y = enriched[,gene.set], 
                    fill = enriched[,split])) + 
            xlab(x.axis) 
        check = NULL}
    plot <- plot + 
        geom_split_violin(alpha=0.8) +
        geom_boxplot(width=0.1, fill = "grey", alpha=0.5, 
            outlier.alpha = 0)  + 
        ylab(paste0(gene.set, " (NES)")) +
        labs(fill = split) + 
        scale_fill_manual(values = colors) + 
        theme_classic() +
        #guides(fill = FALSE)
    if (!is.null(check)) {
        plot <- plot + theme(axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank())}
    return(plot)
}


