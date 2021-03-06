% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-cluster-similarity.R
\name{hcsim_cut}
\alias{hcsim_cut}
\title{Cut filters based on jaccard similarity}
\usage{
hcsim_cut(obj, jac.out = NULL, jac.cut = 0.75)
}
\arguments{
\item{obj}{hcsig object. Please refer to the documentation for \code{statistrics::hcsig}.}

\item{jac.out}{the output vector. Usually and by default equal to NULL -- the clusters that pass jaccard similarity will be added to an empty vector.}

\item{jac.cut}{the cutoff for jaccard similarity.}
}
\value{
an hcsig object (the same structure and data types as the input) filtered to include only non-similar clusters according to \code{jac.cut} val.
}
\description{
Cut filters based on jaccard similarity
}
