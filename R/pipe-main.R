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
#' The output of each step in the pipeline is saved as a .rds file in `cachePath`.
#' If the step has been run before, the saved file will be read in in place of step/function execution, saving runtime.
#' This is enabled by wrapping each function/step in the pipeline with a call to `cacheCall:cacheCall()`.
#' See `?cacheCall::cacheCall` for more details.
#'
#'
#' @param mat a matrix of gene expression data (cells by genes)
#' @param cachePath passed to `cacheCall::cacheCall`. A character string providing path to the Cache directory.
#' @param sep passed to `cacheCall::cacheCall`. A character to separate arguments and their names. Defaults to ":".
#' @param collapse passed to `cacheCall::cacheCall`. A character to separate arguments(+names) from one another. Defaults to "__".
#' @param method.cor a character string defining the distance metric for calculating correlations. Defaults to 'pearson'.
#' @param method.hc a character string defining the type of linkage used in the hierarchical clustering. Defaults to 'average'.
#' @param dissim.dist a numeric value setting the maximum dissimilarity (prior to rescaling between 0 and 1).
#' @param size.cut a boolean indicating whether outputted list of all possible clusters (post-hierarchical clustering) should be filtered based on size (too big or too small).
#' @param min.size a numeric value (ABSOLUTE) for the minimum cluster size (if size.cut=TRUE). Defaults to 5, such that any cluster with < 5 cells is filtered out.
#' @param max.size a numeric value (FRAC. of total) for the maximum cluster size (if size.cut=TRUE). Defaults to 0.5, such that any cluster with >50% of total cell number is filtered out.
#' @param fc.value a numeric value for the minimum fold change of a gene that is counted as differentially expressed. Defaults to 3.
#' @param fc.sort if TRUE, significantly differentially expressed genes are sorted by fold change (highest first). Default is TRUE.
#' @param p.value a numeric value for the minimum p.value of a gene that is counted as significantly differentially expressed. Defaults to 10^(-4).
#' @param pval.sort if TRUE, significantly differentially expressed genes are sorted by p.value (highest first). pval.sort=TRUE overrides fc.sort=TRUE. Default is FALSE.
#' @param reorder.by.sig if TRUE, the list of clusters is reordered by most to least significant.
#' @param n.sig.1 significance cutoff for the number of significantly differentially expressed genes per cluster. Defaults to 50. Any clusters that do not pass this cutoff OR/AND that of n.sig.2 are filtered out.
#' @param n.sig.2 significance cutoff for the number of most significantly differentially expressed genes per cluster. Defaults to 10. Any clusters that do not pass this cutoff OR/AND that of n.sig.1 are filtered out.
#' @param jac.cut a numeric value indicating the cutoff for jaccard similarity. Of two clusters that are above this cutoff, the one with lower significance will be filtered out.
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
                 pval.sort=FALSE,
                 reorder.by.sig=TRUE,
                 n.sig.1=50,
                 n.sig.2=10,
                 jac.cut=0.75,
                 program.cutoff=50) {

  if (is.null(mat)) stop("'Mat' must be provided.")

  args <- as.list(environment())[-c(1:4, 20)]

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
								   fc.sort=fc.sort,
                                   pval.sort=pval.sort,
								   reorder=reorder.by.sig)

  hcsigCut <- cacheCall::cacheCall(pipeName=pipeName,
								   fnName='hcsig_cut',
								   args=args,
								   cachePath=cachePath,
                                   obj=hcsigObj,
								   n.sig.1=n.sig.1,
								   n.sig.2=n.sig.2)

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

  Scores <- cacheCall::cacheCall(pipeName=pipeName,
								 fnName='score',
								 args=args,
								 cachePath=cachePath,
								 mat=mat,
								 programs=Programs,
								 many=TRUE,
								 center=FALSE)

  list(Programs=Programs, Scores=Scores)
}
