#' Title
#'
#' @param List List of matrices to apply FUN to.
#' @param cachePath Output path where files will be saved in / loaded from if they already exist.
#' @param FUN function to apply to list of matrices. Defaults to statistrics::main
#' @param args a list of arguments that will define the filename.
#' @param pipeName a job ID, a name for the project/pipeline.
#' @param ... other arguments to be passed to FUN (statistrics::main)
#'
#' @return output of FUN, for each matrix in turn
#' @export
#'
try_apply <- function(List, cachePath=".", FUN="main", args=NULL, pipeName=NULL, ...) {

  if (is.null(pipeName)) {
    pipeName <- FUN
  }

  if (is.null(names(List))) {
    names(List) <- paste("M", 1:length(List), sep="")
  }

  sapply(1:length(List), function(i) {
    pipeName <- paste(names(List)[i], pipeName, sep="-")
    try(do.call(FUN, list(List[[i]], pipeName=pipeName, cachePath=cachePath, ...)))}, simplify=FALSE)

}


#' Cell Clustering and Cluster Processing
#'
#' A pipeline for cluster analysis of single cell RNA sequencing data.
#' 1. Cells are hierarchically clustered and all possible clusters are retrieved.
#' Cell clusters are then processed as follows:
#' 2. All clusters are independently scored and filtered for significance.
#' A cluster's significance is calculated according to the number and extent of differentially expressed (DE) genes it contains.
#' 3. Clusters are subsequently filtered out by their Jaccard similarity to one another.
#' Of 2 too-similar clusters, the one with higher significance scores is kept.
#' 4. The DE genes defining each cell cluster are used to define cell programs.
#' These programs are sets of coherently expressed genes in the data.
#'
#' The output of each step in the pipeline is saved as a .rds file in \code{cachePath}.
#' If the step has been run before, the saved file will be read in in place of step/function execution, saving runtime.
#' This is enabled by wrapping each function/step in the pipeline with a call to \code{cacheCall:cacheCall()}.
#' See \code{?cacheCall::cacheCall} for more details.
#' Note that setting \code{n.sig.1 = 0} and \code{n.sig.2 = 0} is equivalent to skipping \code{hcsig_cut}.
#' Similarly, setting \code{jac.cut = 0} is equivalent to skipping \code{hcsim_cut}.
#'
#'
#' @param mat a matrix of gene expression data (cells by genes)
#' @param pipeName a job ID, a name for the project/pipeline. Defaults to function name.
#' @param cachePath passed to \code{cacheCall::cacheCall}. A character string providing path to the Cache directory.
#' @param sep passed to \code{cacheCall::cacheCall}. A character to separate arguments and their names. Defaults to ":".
#' @param collapse passed to \code{cacheCall::cacheCall}. A character to separate arguments(+names) from one another. Defaults to "__".
#' @param method.cor a character string of the distance metric for calculating correlations. Defaults to 'pearson'.
#' @param method.hc a character string of the type of linkage used in the hierarchical clustering. Defaults to 'average'.
#' @param dissim.dist a numeric value setting the maximum dissimilarity (prior to rescaling between 0 and 1).
#' @param size.cut a boolean indicating whether outputted list of all possible clusters (post-hierarchical clustering) should be filtered based on size (too big or too small).
#' @param min.size a numeric value (ABSOLUTE) for the minimum cluster size (if size.cut=TRUE). Defaults to 5, such that any cluster with < 5 cells is filtered out.
#' @param max.size a numeric value (FRAC. of total) for the maximum cluster size (if size.cut=TRUE). Defaults to 0.5, such that any cluster with more than half of total cell number is filtered out.
#' @param fc.value fold change value below which differential gene expression is deemed insignificant.
#' @param fc.sort if TRUE, significantly differentially expressed genes are sorted by fold change (highest first). Default is TRUE.
#' @param p.value p-value above which differential gene expression is deemed insignificant.
#' @param p.value.2 p-value above which differential gene expression is deemed insignificant with higher cutoff (sig.2)
#' @param pval.adjust NULL or character string. If NULL, do not adjust p-values. If string, adjust p-values using the method specified.
#' @param pval.sort if TRUE, significantly differentially expressed genes are sorted by p.value (highest first). pval.sort=TRUE overrides fc.sort=TRUE. Default is FALSE.
#' @param reorder.by.sig if TRUE, the list of clusters is reordered by most to least significant.
#' @param n.sig.1 significance cutoff for the number of significantly differentially expressed genes per cluster. Defaults to 50. Any clusters that do not pass this cutoff OR/AND that of n.sig.2 are filtered out. Also see Details.
#' @param n.sig.2 significance cutoff for the number of more significantly differentially expressed genes per cluster. Defaults to 10. Any clusters that do not pass this cutoff OR/AND that of n.sig.1 are filtered out. Also see Details.
#' @param sig.cut.by Cut based on scores for: 'sig.1', 'sig.2', 'either' or 'both'. Passed to \code{hcsig_cut}.
#' @param jac.cut a numeric value indicating the cutoff for jaccard similarity. Of two clusters that are above this cutoff, the one with lower significance will be filtered out. Also see Details.
#' @param program.cutoff a numeric value indicating the cutoff for program sizes. Defaults to 50, such that the maximum number of genes that a program can have is 50.
#'
#' @return a list of programs: coherent sets of genes expressed in the data. Every cluster remaining after filtering steps generates a corresponding program.
#' @export
#'
main <- function(mat,
                 pipeName="main",
                 cachePath=".",
                 sep=":",
                 collapse="_",
                 method.cor="pearson",
                 method.hc="average",
                 dissim.dist=1,
                 size.cut=TRUE,
                 min.size=5,
                 max.size=0.50,
                 fc.value=3,
                 fc.sort=TRUE,
                 p.value=10^(-4),
                 p.value.2=10^(-5),
		             pval.adjust=NULL,
                 pval.sort=FALSE,
                 reorder.by.sig=TRUE,
                 n.sig.1=50,
                 n.sig.2=10,
                 sig.cut.by='both',
                 jac.cut=0.75,
                 program.cutoff=50) {

  args <- as.list(environment())[-c(1:5)]

  if (is.null(mat)) stop("<mat> must be assigned a matrix or a list of matrices.")

  if (is.list(mat)) {
    Programs <- try_apply(List=mat,
                          pipeName=pipeName,
                          args=args,
                          cachePath=cachePath,
                          sep=sep,
                          collapse=collapse,
                          method.cor=method.cor,
                          method.hc=method.hc,
                          dissim.dist=dissim.dist,
                          size.cut=size.cut,
                          min.size=min.size,
                          max.size=max.size,
                          fc.value=fc.value,
                          fc.sort=fc.sort,
                          p.value=p.value,
                          p.value.2=p.value.2,
                          pval.adjust=pval.adjust,
                          pval.sort=pval.sort,
                          reorder.by.sig=reorder.by.sig,
                          n.sig.1=n.sig.1,
                          n.sig.2=n.sig.2,
                          jac.cut=jac.cut,
                          program.cutoff=program.cutoff)
  }

  hc <- cacheCall::cacheCall(pipeName=pipeName,
				      fnName='hcluster',
				      args=args,
				      cachePath=cachePath,
              mat=mat,
				      method.cor=method.cor,
				      method.hc=method.hc,
				      dissim.dist=dissim.dist)

  k <- cacheCall::cacheCall(pipeName=pipeName,
				      fnName='hcutree',
				      args=args,
				      cachePath=cachePath,
              hc=hc,
				      h=hc$height,
				      clean=size.cut,
				      min=min.size,
				      max=max.size)

  hcsigObj <- cacheCall::cacheCall(pipeName=pipeName,
					    fnName='hcsig',
					    args=args,
					    cachePath=cachePath,
              k=k,
				      mat=mat,
					    fc.value=fc.value,
					    p.value=p.value,
					    p.value.2=p.value.2,
				      pval.adjust=pval.adjust,
				      fc.sort=fc.sort,
              pval.sort=pval.sort,
					    reorder=reorder.by.sig)

  hcsigCut <- cacheCall::cacheCall(pipeName=pipeName,
					    fnName='hcsig_cut',
					    args=args,
					    cachePath=cachePath,
              obj=hcsigObj,
					    n.sig.1=n.sig.1,
					    n.sig.2=n.sig.2,
              by=sig.cut.by)

  hcsimCut <- cacheCall::cacheCall(pipeName=pipeName,
					    fnName='hcsim_cut',
					    args=args,
					    cachePath=cachePath,
              obj=hcsigCut,
					    jac.cut=jac.cut)

  Programs <- cacheCall::cacheCall(pipeName=pipeName,
					    fnName='program',
					    args=args,
					    cachePath=cachePath,
              List=hcsimCut$sig.1,
					    cutoff=program.cutoff)

  Programs
}

