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
#' @param complexity.cutoff passed to \code{cutoff} parameter in complexityCut.
#' @param genes.cutoff passed to \code{cutoff} parameter in genesCut.
#' @param logTransform if TRUE, apply \code{logTPM} to matrix.
#' @param ... additional arguments to be passed to \code{logTPM}.
#'
#' @return centered, log-transformed matrix consisting of only high-quality -- user-defined -- cells and genes)
#' @export
#'
preprocess.old <- function(mat, complexity.cutoff, genes.cutoff, logTransform=TRUE, ...) {
  if (isTRUE(logTransform)) mat <- logTPM(tpm=as.matrix(mat), ...)
  cut_cells <- complexityCut(mat, cutoff=complexity.cutoff)
  cut_genes <- genesCut(cut_cells, cutoff=genes.cutoff)
  center(cut_genes, cellcols=TRUE)
}

preprocess <- function(mat,
                       logTransform=TRUE,
                       complexity.cutoff=3000,
                       genes.cutoff=4,
                       ...) {

  mat <- as.matrix(mat)

  if (isTRUE(logTransform)) {
    mat <- cacheCall::cacheCall(pipeName=pipeName, fnName='logTPM', args=args, cachePath=cachePath, tpm=mat)
  }
  cut_cells <- cacheCall::cacheCall(pipeName=pipeName, fnName='complexityCut', args=args, cachePath=cachePath, mat=mat, cutoff=complexity.cutoff)

  cut_genes <- cacheCall::cacheCall(pipeName=pipeName, fnName='genesCut', args=args, cachePath=cachePath, mat=cut_cells, cutoff=genes.cutoff)

  centered <- cacheCall::cacheCall(pipeName=pipeName, fnName='center', args=args, cachePath=cachePath, mat=cut_genes, cellcols=TRUE)

  centered
}
