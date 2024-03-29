% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichment.R
\name{enrich_add_go_level}
\alias{enrich_add_go_level}
\title{Add GO level data to result of enrichment analysis.}
\usage{
enrich_add_go_level(enr, go_data)
}
\arguments{
\item{enr}{A \code{enrichResult} instance}

\item{go_data}{A data.frame contains four columns:
\enumerate{
  \item First column must be named as \code{Gene_ID}
  \item Second column must be named as \code{GO_ID}
  \item Third column must be named as \code{GO_Name}
  \item Fourth column must be named as \code{GO_Level}
}}
}
\value{
A data.frame
}
\description{
Add GO level data to result of enrichment analysis.
}
