# ======================
# PREPROCESSING

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
#' @param cellcols boolean. The columns (variables) are cells if TRUE or are genes if FALSE.
#'
#' @return a matric with centered gene expression data.
#' @export
#'
center <- function(mat, cellcols=TRUE) {

  # Genes as columns
  if(cellcols) {
    # save column names since they are lost in apply function below
    colNames <- colnames(mat)
    # centered expression for each gene in each cell relative to column mean
    centered <- apply(mat, 1, scale, center=TRUE, scale=FALSE)
    # revert to Cells as columns
    centered <- t(centered)
  }

  else {
    # save column names since they are lost in apply function below
    colNames <- rownames(mat)
    # centered expression for each gene in each cell relative to column mean
    centered <- apply(mat, 2, scale, center=TRUE, scale=FALSE)
  }

  # add column names back to centered matrix
  colnames(centered) <- colNames
  return(centered)
}


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
preprocess <- function(mat, complexity.cutoff, genes.cutoff, logTransform=TRUE, ...) {
  if (isTRUE(logTransform)) mat <- logTPM(as.matrix(mat), ...)
  cut_cells <- complexityCut(mat, cutoff=complexity.cutoff)
  cut_genes <- genesCut(cut_cells, cutoff=genes.cutoff)
  center(cut_genes, cellcols=TRUE)
}
# ======================

#' Matrix Column/Row Means
#'
#' @param mat matrix of vars. by. obs.
#' @param by integer of 1 or 2. 1 for row means or 2 for column means.
#'
#' @return vector of mean values
#' @export
#'
gene_means <- function(mat, by=1) { apply(mat, by, mean) }

#' Fold Changes between Matrices
#'
#' @param a matrix (vars. by. obs.) for which fold changes are calculated
#' @param b matrix (vars. by. obs.) for which the fold changes of matrix 'a' are calculated relative to
#' @param log logical. if TRUE, calculates fold change in log space.
#' @param log.base base of log. Defaults to 2.
#' @param fc.value NULL or numeric. If NULL, return all FC values. If numeric, return FC values >= fc.value.
#'
#' @return numeric vector containing FC values
#' @export
#'
fold_change <- function(a, b, log=TRUE, log.base=2, fc.value=2){

  a <- as.matrix(a)
  b <- as.matrix(b)

  if (length(dim(a)) == 2) {
    a <- gene_means(a)
    b <- gene_means(b)}

  if (log) {
    logFC <- a - b
    FC <- log.base^logFC}

  else FC <- a/b

  if (is.null(fc.value)) return(FC)

  else if (!is.null(fc.value)) return(FC[FC > fc.value])

}


#' P-values between matrices
#'
#' @param a matrix (vars. by. obs.) for which p-values are calculated.
#' @param b matrix (vars. by. obs.) for which the p-values of matrix 'a' are calculated relative to.
#' @param p.value NULL or numeric. If NULL, return all p-values. If numeric, return p-values <= p.value.
#' @param adjust.method NULL or character string. If NULL, do not adjust p-values. If string, adjust p-values using the method specified.
#'
#' @return numeric vector containing p-values.
#' @export
#'
#' @importFrom stats t.test p.adjust
#'
p_val <- function(a, b, p.value=NULL, adjust.method="bonferroni"){

  p <- sapply(1:nrow(a), function(i) stats::t.test( a[i,], b[i,] )$p.value)

  if (!is.null(adjust.method)) p <- stats::p.adjust(p, method=adjust.method)

  names(p) <- rownames(a)

  if (is.null(p.value)) return(p)

  else if (!is.null(p.value)) return(p[p < p.value])

}


# =============================================
#' Significance according to p-values of fold changes.
#'
#' @param a matrix (vars. by. obs.) for which significance (of each obs.) is calculated
#' @param b matrix (vars. by. obs.) for which significance (of each obs.) in matrix 'a' is calculated relative to.
#' @param p.value arg to \code{p_val}. NULL or numeric. If NULL, return all p-values. If numeric, return p-values <= p.value.
#' @param fc.value arg to \code{fold_change}. NULL or numeric. If NULL, return all FC values. If numeric, return FC values >= fc.value.
#' @param ... other args passed to \code{fold_change} or \code{p_val}.
#'
#' @return one of: i) all p-values for all fold changes (if fc.value == NULL & p.value == NULL) ii/iii) all/signifcant p-values for significant/all fold changes or iv) signifcant p-values for significant fold changes (if fc.value != NULL & p.value != NULL).
#' @export
#'
sig <- function(a, b, p.value, fc.value=2, ...){

  fc <- fold_change(a, b, fc.value=fc.value, ...)

  if (length(fc) == 0) {return()}

  else if (length(fc) > 0){

    a <- a[names(fc), ]
    b <- b[names(fc), ]

    p_val(a, b, p.value=p.value, ...)}

}


#' Find the most significant occurrence of each observation
#'
#' @param x numeric vector of (named) observations
#' @param y numeric vector of (named) observations
#'
#' @return numeric vectors x and y, trimmed, such that each observation that was present in both now appears only in that group where its occurrence was most significant
#' @export
#'
p_2 <- function(x, y){

  sub.x <- x[names(x) %in% names(y)]
  sub.x <- sub.x[sort(names(sub.x))]
  not.sub.x <- x[!(names(x) %in% names(y))]

  sub.y <- y[names(y) %in% names(x)]
  sub.y <- sub.y[sort(names(sub.y))]
  not.sub.y <- y[!(names(y) %in% names(x))]

  subd.x <- sub.x[sub.x <= sub.y]
  subd.y <- sub.y[sub.y <= sub.x]

  sig.x <- c(not.sub.x, subd.x)
  sig.y <- c(not.sub.y, subd.y)

  l <- list(sig.x, sig.y)
  names(l) <- c('x', 'y')

  # return(l)
  return(l$x)

}

# =============================================
#' Basic Similarity Operation (intersection value)
#'
#' @param a character vector
#' @param b character vector
#'
#' @return overlap (intersection) value
#' @export
#'
overlap <- function(a, b) {
  length(intersect(a, b))
}

#' Jaccard Similarity Operation
#'
#' @param a character vector
#' @param b character vector
#'
#' @return jaccard similarity value
#' @export
#'
jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

#' Apply function to every pair between x and a list
#'
#' @param a value or vector
#' @param vect list of values or vectors of the same type as 'a'
#' @param FUN function to be applied on every pair. eg \code{sum(a, vect[[1]])} if a and vect are numeric
#'
#' @return result of function call for each pair
#' @export
#'
pairs_apply <- function(a, vect, FUN=overlap) {
  sapply(vect, function(b) FUN(a=a, b=b))
}

#' Apply function to every pair in a list
#'
#' @param vect a list. All elements must be of the same type.
#' @param FUN function to be applied on every pair. eg sum(a, vect[[1]]) if a and vect are numeric
#'
#' @return result of function call for each pair
#' @export
#'
pairs_apply_all <- function(vect, FUN=overlap) {

  overlaps <- NULL

  for (i in 1:length(vect)) {
    x <- vect[[i]]
    overlap.scores <- list(pairs_apply(a=x,
                                        vect=vect,
                                        FUN=FUN))
    overlaps <- c(overlaps, overlap.scores)
  }

  overlaps
}

#' Jaccard Test
#'
#' @param x character vector
#' @param y character vector
#' @param jac.cut jaccard value cut off.
#'
#' @return logical. TRUE if vectors' similarity is below or equal to cut off. FALSE otherwise.
#' @export
#'
jaccard_test <- function(x, y, jac.cut=0.75) {

  a <- length(intersect(x, y))
  b <- length(union(x,y))
  jac <- a/b

  if(jac <= jac.cut) {return(TRUE)}
  else {return(FALSE)}

}

#' Vector Filtering by Jaccard Similarity
#'
#' One of two vectors with too high jaccard similarity is removed. Priority is given to groups that appear higher in the list of vectors, such that the lower-placed vector is removed.
#'
#' @param x list of named vectors. Vectors are assumed to be ordered from most to least significant, such that of two vectors with a similarity score higher than the cut-off, the higher-placed vector will be kept and the lower-placed vector will be thrown.
#' @param out vectors that passed jaccard test.
#' @param jac.cut jaccard value cut off. arg to \code{jaccard_test}.
#'
#' @return list of i) x: names of vectors that passed jaccard test and ii) out: vectors that passed jaccard test.
#' @export
#'
filter_by_jaccard <- function(x, out=NULL, jac.cut=0.75) {

  if (is.null(names(x))) stop("Vectors must be named.")
  if (length(x) <= 1 && is.null(out)) stop("A minimum of 2 vectors is needed for the jaccard test.")
  else if (length(x) <= 1) return(out)

  # keep first vector
  out <- c(out, x[1])
  # pairs first vector with every other in list
  pairs <- expand.grid(x[1], x[2:length(x)])
  # names of all but the first vector in list
  names <- names(pairs[[2]])
  # good vectors
  boolean <- apply(pairs, 1, function(pair) {
					 jaccard_test(pair[[1]], pair[[2]], jac.cut=jac.cut)})
  names <- names[boolean]

  # proceed with good vectors
  x <- x[names(x) %in% names]

  filter_by_jaccard(x=x, out=out)
}

# =============================================


#' Trim vectors of observations such that only the most significant occurrence of each observation is kept.
#'
#' @param List a list of vectors of observations. E.g. p-values named by the genes they belong to. some genes (observations) will have multiple occurrences in different vectors, and only the vector with the most significant read-out for that gene gets to keep it.
#'
#' @return a list of vectors of observations. Each observation (gene) only appears once.
#' @export
#'
sig_2 <- function(List) {

  remove <- function(l, i){
    l2 <- l
    l2[[i]] <- NULL
    l2}


  sapply(1:length(List), function(i) {

    x <- List[[i]]
    y <- remove(List, i)

    for(i in 1:length(y)) {

      x <- p_2(x, y[[i]]) }

    return(x) })

}


#' Differentially Expressed Genes Between Two Groups
#'
#' @param k character vector of vars, for which a subgroup will have DE genes calculated.
#' @param mat matrix of vars. vs. obs. matrix is split into two: i) vars that are in k; ii) vars not in k.
#' @param fc.value fold change value below which differential gene expression is deemed insignificant.
#' @param p.value p-value above which differential gene expression is deemed insignificant.
#' @param ... other args to be passed to \code{p_val} or \code{fold_change} through call to \code{sig}.
#'
#' @return numeric vector of p-values named with gene names.
#' @export
#'
DEgenes <- function(k, mat, fc.value=2, p.value=10^(-3), ...) {
  a <- mat[, k]
  b <- mat[, !colnames(mat) %in% k]
  out <- sig(a, b, fc.value=2, p.value=10^(-3), ...)
  if (length(out) == 0) stop(paste("No DE genes with p value <=", p.value))
  sort(out)
}


