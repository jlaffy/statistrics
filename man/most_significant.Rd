% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-cluster-significance.R
\name{most_significant}
\alias{most_significant}
\title{Single Most Significant Occurrences of Genes}
\usage{
most_significant(List, by = "small")
}
\arguments{
\item{List}{a list of named numeric vectors. Names: genes; values: their p-values; vectors: clusters derived from statistrics::hcluster.}

\item{by}{'small' or 'large' depending on whether the most significant values are the smallest (eg. p-value) or largest (eg. fold change).}
}
\value{
the inputted List after filtering such that each gene only appears once in the cluster where it is most significant.
}
\description{
Most significant occurrence of each gene.
The data generated allows the parameter n.sig.2 to be generated,
for the number of genes in a cluster with the most significant occurrences.
}
