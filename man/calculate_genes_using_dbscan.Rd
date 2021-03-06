% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Deconvolution.R
\name{calculate_genes_using_dbscan}
\alias{calculate_genes_using_dbscan}
\title{Function that selects the markers genes given a wilcox object by using dbscan algorithm and selecting just the genes that are considered like outliers (Without any cluster)}
\usage{
calculate_genes_using_dbscan(filtered.wilcox, plot.dbscan.results = FALSE)
}
\arguments{
\item{filtered.wilcox}{Wilcox object with the Foldchange analysis base on the comparison between each cluster amoung the others}
}
\value{
List with the marker genes
}
\description{
Function that selects the markers genes given a wilcox object by using dbscan algorithm and selecting just the genes that are considered like outliers (Without any cluster)
}
