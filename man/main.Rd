% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipe-main.R
\name{main}
\alias{main}
\title{Cell Clustering and Cluster Processing}
\usage{
main(mat, pipeName = "main", cachePath = ".", sep = ":", collapse = "_",
  method.cor = "pearson", method.hc = "average", dissim.dist = 1,
  size.cut = TRUE, min.size = 5, max.size = 0.5, fc.value = 3,
  fc.sort = TRUE, p.value = 10^(-4), p.value.2 = 10^(-5),
  pval.adjust = NULL, pval.sort = FALSE, reorder.by.sig = TRUE,
  n.sig.1 = 50, n.sig.2 = 10, sig.cut.by = "both", jac.cut = 0.75,
  program.cutoff = 50)
}
\arguments{
\item{mat}{a matrix of gene expression data (cells by genes)}

\item{pipeName}{a job ID, a name for the project/pipeline. Defaults to function name.}

\item{cachePath}{passed to \code{cacheCall::cacheCall}. A character string providing path to the Cache directory.}

\item{sep}{passed to \code{cacheCall::cacheCall}. A character to separate arguments and their names. Defaults to ":".}

\item{collapse}{passed to \code{cacheCall::cacheCall}. A character to separate arguments(+names) from one another. Defaults to "__".}

\item{method.cor}{a character string of the distance metric for calculating correlations. Defaults to 'pearson'.}

\item{method.hc}{a character string of the type of linkage used in the hierarchical clustering. Defaults to 'average'.}

\item{dissim.dist}{a numeric value setting the maximum dissimilarity (prior to rescaling between 0 and 1).}

\item{size.cut}{a boolean indicating whether outputted list of all possible clusters (post-hierarchical clustering) should be filtered based on size (too big or too small).}

\item{min.size}{a numeric value (ABSOLUTE) for the minimum cluster size (if size.cut=TRUE). Defaults to 5, such that any cluster with < 5 cells is filtered out.}

\item{max.size}{a numeric value (FRAC. of total) for the maximum cluster size (if size.cut=TRUE). Defaults to 0.5, such that any cluster with more than half of total cell number is filtered out.}

\item{fc.value}{fold change value below which differential gene expression is deemed insignificant.}

\item{fc.sort}{if TRUE, significantly differentially expressed genes are sorted by fold change (highest first). Default is TRUE.}

\item{p.value}{p-value above which differential gene expression is deemed insignificant.}

\item{p.value.2}{p-value above which differential gene expression is deemed insignificant with higher cutoff (sig.2)}

\item{pval.adjust}{NULL or character string. If NULL, do not adjust p-values. If string, adjust p-values using the method specified.}

\item{pval.sort}{if TRUE, significantly differentially expressed genes are sorted by p.value (highest first). pval.sort=TRUE overrides fc.sort=TRUE. Default is FALSE.}

\item{reorder.by.sig}{if TRUE, the list of clusters is reordered by most to least significant.}

\item{n.sig.1}{significance cutoff for the number of significantly differentially expressed genes per cluster. Defaults to 50. Any clusters that do not pass this cutoff OR/AND that of n.sig.2 are filtered out. Also see Details.}

\item{n.sig.2}{significance cutoff for the number of more significantly differentially expressed genes per cluster. Defaults to 10. Any clusters that do not pass this cutoff OR/AND that of n.sig.1 are filtered out. Also see Details.}

\item{sig.cut.by}{Cut based on scores for: 'sig.1', 'sig.2', 'either' or 'both'. Passed to \code{hcsig_cut}.}

\item{jac.cut}{a numeric value indicating the cutoff for jaccard similarity. Of two clusters that are above this cutoff, the one with lower significance will be filtered out. Also see Details.}

\item{program.cutoff}{a numeric value indicating the cutoff for program sizes. Defaults to 50, such that the maximum number of genes that a program can have is 50.}
}
\value{
a list of programs: coherent sets of genes expressed in the data. Every cluster remaining after filtering steps generates a corresponding program.
}
\description{
A pipeline for cluster analysis of single cell RNA sequencing data.
1. Cells are hierarchically clustered and all possible clusters are retrieved.
Cell clusters are then processed as follows:
2. All clusters are independently scored and filtered for significance.
A cluster's significance is calculated according to the number and extent of differentially expressed (DE) genes it contains.
3. Clusters are subsequently filtered out by their Jaccard similarity to one another.
Of 2 too-similar clusters, the one with higher significance scores is kept.
4. The DE genes defining each cell cluster are used to define cell programs.
These programs are sets of coherently expressed genes in the data.
}
\details{
The output of each step in the pipeline is saved as a .rds file in \code{cachePath}.
If the step has been run before, the saved file will be read in in place of step/function execution, saving runtime.
This is enabled by wrapping each function/step in the pipeline with a call to \code{cacheCall:cacheCall()}.
See \code{?cacheCall::cacheCall} for more details.
Note that setting \code{n.sig.1 = 0} and \code{n.sig.2 = 0} is equivalent to skipping \code{hcsig_cut}.
Similarly, setting \code{jac.cut = 0} is equivalent to skipping \code{hcsim_cut}.
}
