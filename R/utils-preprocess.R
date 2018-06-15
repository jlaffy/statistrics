#' Convert expression data from logTPM to TPM form.
#'
#' TPM aka Transcripts Per Million.
#' If a matrix is provided, \code{TPM} will be applied over each cell and returned as a matrix.
#' Note. \code{TPM(logtpm) == logTPM(tpm)}

#' @param logtpm a numeric vector. Either a matrix or of length 1.
#' @param dividebyten a boolean. Usually TRUE when single cells expression vals are concerned; usually FALSE when bulk (many cell) expression vals are concerned.
#'
#' @return a numeric vector (either a matrix or of length 1) with values converted from logTPM to TPM form.
#' @export
#'
TPM <- function(logtpm, dividebyten=TRUE) {
  if(dividebyten) {
    tpm <- 10*(2^(logtpm)-1)}
  else if(!dividebyten) {
    tpm <- 2^(logtpm)-1}
  return(tpm)
}


#' Convert expression data from TPM to logTPM form.
#'
#' TPM aka Transcripts Per Million.
#' If a matrix is provided, \code{logTPM} will be applied over each cell and returned as a matrix.
#' Note. \code{logTPM(tpm) == TPM(logtpm)}

#' @param tpm a numeric vector. Either a matrix or of length 1.
#' @param dividebyten a boolean. Usually TRUE when single cells expression vals are concerned; usually FALSE when bulk (many cell) expression vals are concerned.
#'
#' @return a numeric vector (either a matrix or of length 1) with values converted from TPM to logTPM form.
#' @export
#'
logTPM <- function(tpm, dividebyten=TRUE) {
  # same as TPM() applies.
  # logTPM(tpm) == TPM(logtpm)

  if(dividebyten) {
    logtpm <- log(tpm/10+1, 2)}
  else if(!dividebyten) {
    logtpm <- log(tpm+1, 2)}
  return(logtpm)
}


#3
#' Cell complexity counter
#'
#' @param mat a matrix of cells-by-genes or genes-by-cells containing gene expression data.
#' @param cellcols boolean. The columns (variables) are cells if TRUE or are genes if FALSE.
#' @param rank boolean. If TRUE, order cells by # of non-zero genes detected
#'
#' @return complexity value per cell (# of genes detected per cell)
#' @export
#'
detected <- function(mat, cellcols=TRUE, rank=TRUE) {
  # Calculates # of genes detected (complexity) per cell.
  # mat is a matrix with GENE ROWS and CELL COLUMNS.

  # if cells are rows and not columns, first transpose the matrix.
  if(!cellcols) {mat <- t(mat)}

  # for each column in the matrix, sums the number of genes with non-zero vals.
  mat2 <- apply(mat, 2, function(x) sum(x!=0))

  if(rank) {return(sort(mat2))}
  else {return(mat2)}

}

#' Filter cells by complexity cutoff
#'
#' Removes all cells whose detected gene number/complexity value is below the cutoff.
#'
#' @param mat a matrix of cells-by-genes containing gene expression data.
#' @param cutoff complexity value below which cells are filtered out.
#' @param justNames boolean. If TRUE, return only the names of cells that passed the filtering.
#'
#' @return cell IDs whose detected gene number >= 3000 and -- if \code{justNames==FALSE} -- their detected values.
#' @export
#'
complexityCut <- function(mat, cutoff=3000, justNames=FALSE) {

  d <- detected(mat, rank=FALSE)

  # [d >= cutoff] is a logical vector of T (>=3000) and F (<3000).
  # It is passed to d and masks/deletes false values in d
  true3000 <- d[d >= cutoff, drop=FALSE]
  goodCells <- names(true3000)

  # return vector of cell IDs with accepted complexity vals.
  if(justNames) {return(goodCells)}

  # return mat with low complexity Cells (columns) removed.
  else {return(mat[,goodCells])}

}

#' Filter genes by cutoff
#'
#' Removes genes whose average expression across cells is below the cutoff.
#'
#' @param mat a matrix of cells-by-genes containing gene expression data.
#' @param cutoff gene expression value below which genes are filtered out.
#'
#' @return matrix without the genes/rows whose aggregate expression across cells was too low.
#' @export
#'
genesCut <- function(mat, cutoff=4) {
  # remove genes with ~no reads across cells

  # remove log
  mat.tpm <- TPM(mat)
  # average expression of each gene across cells
  mat.tpm.avg <- apply(mat.tpm, 1, mean)
  # add 1 and log (dividebyten=F since avg across cells)
  avg.logtpm1 <- logTPM(mat.tpm.avg, dividebyten=F)
  # remove genes with ~no expression
  mat.cut <- avg.logtpm1[avg.logtpm1 >= cutoff]
  # selected gene names
  genes <- names(mat.cut)
  # subset data for selected geenes only
  return(mat[genes,])

}


#' Centering
#'
#' Center gene expression data (before computing correlations), such that the average correlation between cells will be 0.
#' Expression(Gene a, Cell k) - Avg(Expression(Gene a, Cell 1 to n))
#' mat is a matrix with cell rows and gene columns.
#'
#' @param mat a matrix of cells-by-genes or genes-by-cells containing gene expression data.
#' @param rowWise boolean. Center by rows if TRUE and by columns if FALSE.
#'
#' @return a matrix with centered gene expression data.
#' @export
#'
center <- function(mat, rowWise = TRUE) {

  if (isTRUE(rowWise)) {
    mat.centered <- t(scale(t(mat), center = TRUE, scale = FALSE))
  }

  else if (!isTRUE(rowWise)) {
    mat.centered <- scale(mat, center = TRUE, scale = FALSE)
  }

  mat.centered
}

