% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-preprocess.R
\name{center}
\alias{center}
\title{Centering}
\usage{
center(mat, rowWise = TRUE, center.by = "mean")
}
\arguments{
\item{mat}{a matrix of cells-by-genes or genes-by-cells containing gene expression data.}

\item{rowWise}{boolean. Center by rows if TRUE and by columns if FALSE.}
}
\value{
a matrix with centered gene expression data.
}
\description{
Center gene expression data (before computing correlations), such that the average correlation between cells will be 0.
Expression(Gene a, Cell k) - Avg(Expression(Gene a, Cell 1 to n))
mat is a matrix with cell rows and gene columns.
}
