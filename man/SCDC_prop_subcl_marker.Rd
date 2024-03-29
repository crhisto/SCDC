% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Deconvolution.R
\name{SCDC_prop_subcl_marker}
\alias{SCDC_prop_subcl_marker}
\title{Tree-guided proportion estimation}
\usage{
SCDC_prop_subcl_marker(
  bulk.eset,
  sc.eset,
  ct.varname,
  fl.varname,
  sample,
  ct.sub = NULL,
  ct.fl.sub,
  iter.max = 3000,
  nu = 1e-04,
  epsilon = 0.001,
  weight.basis = T,
  truep = NULL,
  select.marker = T,
  markers = NULL,
  marker.varname = NULL,
  allgenes.fl = F,
  pseudocount.use = 1,
  LFC.lim = 0.5,
  parallelize = F,
  core_number = NULL,
  fix_number_genes = NULL,
  marker_gene_strategy = "boostrap_outliers",
  iteration.minimun_number_markers = 28,
  iteration.use_maximum = FALSE,
  iteration.maximo_genes = 35,
  iteration.use_final_foldchange = FALSE,
  bootstrap.sample_size = NULL,
  bootstrap.number = NULL,
  additional_genes = NULL,
  ...
)
}
\arguments{
\item{bulk.eset}{ExpressionSet object for bulk samples}

\item{sc.eset}{ExpressionSet object for single cell samples}

\item{ct.varname}{variable name for 'cell types'}

\item{fl.varname}{variable name for first-level 'meta-clusters'}

\item{sample}{variable name for subject/samples}

\item{ct.sub}{a subset of cell types that are selected to construct basis matrix}

\item{ct.fl.sub}{'cell types' for first-level 'meta-clusters'}

\item{iter.max}{the maximum number of iteration in WNNLS}

\item{nu}{a small constant to facilitate the calculation of variance}

\item{epsilon}{a small constant number used for convergence criteria}

\item{weight.basis}{logical, use basis matrix adjusted by MVW, default is T.}

\item{truep}{true cell-type proportions for bulk samples if known}

\item{select.marker}{logical, select marker genes to perform deconvolution in tree-guided steps. Default is T.}

\item{markers}{A set of marker gene that input manually to be used in deconvolution. If NULL, then}

\item{marker.varname}{variable name of cluster groups when selecting marker genes. If NULL, then use ct.varname.}

\item{allgenes.fl}{logical, use all genes in the first-level deconvolution}

\item{pseudocount.use}{a constant number used when selecting marker genes, default is 1.}

\item{LFC.lim}{a threshold of log fold change when selecting genes as input to perform Wilcoxon's test.}

\item{iteration.use_final_foldchange}{TRUE/FALSE. If at the end the cluster has zero genes if this parameter is true, the boostraping is going to be calculated over the foldchange with <0.05, not with zero.}
}
\value{
Estimated proportion, basis matrix, predicted gene expression levels for bulk samples
}
\description{
Proportion estimation function for multi-subject case, and apply tree-guided deconvolution
}
