% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_significant_tfs.R
\name{get_significant_tfs}
\alias{get_significant_tfs}
\title{Run DoRothEA analysis to get significant transcription factors}
\usage{
get_significant_tfs(
  seuratobject,
  condition,
  out_path,
  tf_condition_significant,
  pval,
  log2fc
)
}
\arguments{
\item{seuratobject}{Input Seurat Object}

\item{condition}{Experminet condition (e.g. disease, knockout ...)}

\item{out_path}{Output path to save results}

\item{tf_condition_significant}{condition comparison results}

\item{pval}{p-value to filter results}

\item{log2fc}{log fold change value to filter results}
}
\value{
A data frame with transcription factor activity scores per cell type
}
\description{
Description
}
