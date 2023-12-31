% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualizations.R
\name{ridgeEnrichment}
\alias{ridgeEnrichment}
\title{Generate a ridge plot to examine enrichment distributions}
\usage{
ridgeEnrichment(
  enriched,
  group = "cluster",
  gene.set = NULL,
  scale.bracket = NULL,
  facet = NULL,
  add.rug = FALSE,
  colors = c("#0D0887FF", "#7E03A8FF", "#CC4678FF", "#F89441FF", "#F0F921FF")
)
}
\arguments{
\item{enriched}{The output of \code{\link{enrichIt}}}

\item{group}{The parameter to group, displayed on the y-axis.}

\item{gene.set}{The gene set to graph on the x-axis.}

\item{scale.bracket}{This will filter the enrichment scores to remove 
extreme outliers. Values entered (1 or 2 numbers) will be the filtering 
parameter using z-scores of the selected gene.set. If only 1 value is given, 
a seocndary bracket is autommatically selected as the inverse of the number.}

\item{facet}{A parameter to separate the graph.}

\item{add.rug}{Binary classifier to add a rug plot to the x-axis.}

\item{colors}{The color palette for the ridge plot.}
}
\value{
ggplot2 object with ridge-based distributions of selected gene.set
}
\description{
This function allows to the user to examine the distribution of 
enrichment across groups by generating a ridge plot.
}
\examples{
ES2 <- readRDS(url(
"https://ncborcherding.github.io/vignettes/escape_enrichment_results.rds"))
ridgeEnrichment(ES2, gene.set = "HALLMARK_DNA_REPAIR", group = "cluster", 
facet = "Type", add.rug = TRUE)

}
\seealso{
\code{\link{enrichIt}} for generating enrichment scores.
}
