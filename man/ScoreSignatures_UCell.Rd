% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UCell.R
\name{ScoreSignatures_UCell}
\alias{ScoreSignatures_UCell}
\title{Calculate module enrichment scores from single-cell data}
\usage{
ScoreSignatures_UCell(
  matrix = NULL,
  features,
  precalc.ranks = NULL,
  maxRank = 1500,
  w_neg = 1,
  name = "_UCell",
  assay = "counts",
  chunk.size = 1000,
  ncores = 1,
  ties.method = "average",
  force.gc = FALSE,
  seed = 123
)
}
\arguments{
\item{matrix}{Input matrix,}

\item{features}{A list of signatures}

\item{precalc.ranks}{If you have pre-calculated ranks}

\item{maxRank}{Maximum number of genes to rank per cell; above this rank, a given gene is considered as not expressed.}

\item{w_neg}{Weight on negative genes in signature. e.g. `w_neg=1` weighs equally up- and down-regulated genes, `}

\item{name}{Name suffix appended to signature names}

\item{assay}{The sce object assay where the data is to be found}

\item{chunk.size}{Number of cells to be processed simultaneously (lower size requires slightly more computation but reduces memory demands)}

\item{ncores}{Number of processors to parallelize computation.}

\item{ties.method}{How ranking ties should be resolved (passed on to [data.table::frank])}

\item{force.gc}{Explicitly call garbage collector to reduce memory footprint}

\item{seed}{Integer seed}
}
\value{
Returns input SingleCellExperiment object with UCell scores added to altExp
}
\description{
Calculate module enrichment scores from single-cell data
}
