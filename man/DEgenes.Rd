% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-cluster-significance.R
\name{DEgenes}
\alias{DEgenes}
\title{Differentially Expressed Genes Between Two Groups}
\usage{
DEgenes(k, mat, fc.value = 3, p.value = 10^(-4), fc.sort = T,
  pval.sort = F, adjust.method = NULL, returning = "pval", ...)
}
\arguments{
\item{k}{a list of character vectors; sets of cell members belonging to each cluster.}

\item{mat}{matrix of vars. vs. obs. matrix is split into two: i) vars that are in k; ii) vars not in k.}

\item{fc.value}{fold change value below which differential gene expression is deemed insignificant.}

\item{p.value}{p-value above which differential gene expression is deemed insignificant.}

\item{fc.sort}{if TRUE, significantly differentially expressed genes are sorted by fold change (highest first). Default is TRUE.}

\item{pval.sort}{if TRUE, significantly differentially expressed genes are sorted by p.value (highest first). pval.sort=TRUE overrides fc.sort=TRUE. Default is FALSE.}

\item{adjust.method}{NULL or character string. If NULL, do not adjust p-values. If string, adjust p-values using the method specified.}

\item{returning}{return one of p-values, fold changes, or both from the output of \code{sig()} call.}

\item{...}{other args to be passed to \code{p_val} or \code{fold_change} through call to \code{sig}.}
}
\value{
numeric vector of p-values named with gene names.
}
\description{
Differentially Expressed Genes Between Two Groups
}
