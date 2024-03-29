% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichment.R
\name{compareGSEA}
\alias{compareGSEA}
\title{Perform GSEA on gene clusters.}
\usage{
compareGSEA(ranked_clusters, ...)
}
\arguments{
\item{ranked_clusters}{a list of ranked logFC vector.}

\item{...}{pass to \code{\link[clusterProfiler]{gseGO}}}
}
\value{
A un-official \code{compareClusterResult} object
}
\description{
Perform GSEA on gene clusters.
}
