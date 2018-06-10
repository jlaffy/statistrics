#' PreProcess Expression Matrix
#'
#' If the data is passed in TPM form, the data is first transformed to logTPM.
#' Low quality cells are then removed. They are a filtered according to a cutoff for the number of genes detected.
#' Next, low quality genes are removed. These are filtered according to a cutoff for their avg. expression values across all cells.
#' Finally, the data is centered such that the average expression of each gene across all cells is 0.
#' Note: The standard deviation is left as is -- rather than being normalised to equal 1 --
#' since s.d. values are very skewed due to how sparse the scRNA data is (many 0s).
#'
#' @param mat matrix of vars. by. obs. (in TPM or logTPM).
#' @param pipeName a job ID, a name for the project/pipeline. Defaults to function name.
#' @param cachePath passed to \code{cacheCall::cacheCall}. A character string providing path to the Cache directory.
#' @param complexity.cutoff if a numeric value, apply \code{complexityCut} to matrix with cutoff \code{complexity.cutoff}. Else if FALSE, do not \code{complexityCut}.
#' @param genes.cutoff if a numeric value, apply \code{genesCut} to matrix with cutoff \code{genes.cutoff}. Else if FALSE, do not \code{genesCut}.
#' @param centering if TRUE, apply \code{center} to matrix.
#' @param logTransform if TRUE, apply \code{logTPM} to matrix.
#'
#' @return if all steps are run, return a centered, log-transformed matrix consisting of only high-quality -- user-defined -- cells and genes)
#' @export
#'
preprocess <- function(mat,
                       pipeName="preprocess",
                       cachePath=".",
                       logTransform=TRUE,
                       complexity.cutoff=3000,
                       genes.cutoff=4,
		       centering=TRUE) {

  if (is.null(mat)) stop("'Mat' must be provided.")

  args <- as.list(environment())[-c(1:4)]

  mat <- as.matrix(mat)

  if (isTRUE(logTransform)) {
    mat <- cacheCall::cacheCall(pipeName=pipeName,
                                fnName='logTPM',
                                args=args,
                                cachePath=cachePath,
                                tpm=mat)
  }

  if (complexity.cutoff != F & is.numeric(complexity.cutoff)) {
	mat <- cacheCall::cacheCall(pipeName=pipeName,
                                fnName='complexityCut',
                                args=args,
                                cachePath=cachePath,
                                mat=mat,
                                cutoff=complexity.cutoff)
  }

  if (genes.cutoff != F & is.numeric(genes.cutoff)) {
    mat <- cacheCall::cacheCall(pipeName=pipeName,
                                fnName='genesCut',
                                args=args,
                                cachePath=cachePath,
                                mat=mat,
                                cutoff=genes.cutoff)
  }

  if (isTRUE(centering)) {
    mat <- cacheCall::cacheCall(pipeName=pipeName,
                                fnName='center',
                                args=args,
                                cachePath=cachePath,
                                mat=mat,
                                rowWise=TRUE)
  }

  mat

}
