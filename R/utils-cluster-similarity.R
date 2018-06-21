#' Basic Similarity (intersection value)
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


#' Jaccard Similarity
#'
#' @param x character vector
#' @param y character vector
#' @param jac.cut jaccard value cut off.
#' @param test a boolean indicating whether to test jaccard value against cutoff \code{jac.cut}.
#'
#' @return If no cutoff tested, returns jaccard similarity. If cutoff applied, returns TRUE if jaccard similarity is below or equal to cutoff and FALSE otherwise.
#' @export
#'
jaccard <- function(x, y, jac.cut=0.5, test=TRUE) {

  a <- length(intersect(x, y))
  b <- length(union(x,y))
  jac <- a/b

  if (is.null(jac.cut) | !isTRUE(test)) {

	return(jac)
  }

  if (!is.null(jac.cut) & !is.numeric(jac.cut)) {

	stop("Cutoff arg <jac.cut> is ",
		 class(jac.cut),
		 ". Arg must be numeric.",
		 sep="")
  }

  else if (!is.null(jac.cut)) {
    if(jac <= jac.cut) return(TRUE)
    else return(FALSE)
  }
}


# ==============================================================
# Cluster similarity helper
# ==============================================================

#' Vector Filtering by Jaccard Similarity
#'
#' One of two vectors with too high jaccard similarity is removed. Priority is given to groups that appear higher in the list of vectors, such that the lower-placed vector is removed.
#'
#' @param x list of named vectors. Vectors are assumed to be ordered from most to least significant, such that of two vectors with a similarity score higher than the cut-off, the higher-placed vector will be kept and the lower-placed vector will be thrown.
#' @param out vectors that passed jaccard test.
#' @param jac.cut jaccard value cut off. arg to \code{jaccard}.
#'
#' @return list of i) x: names of vectors that passed jaccard test and ii) out: vectors that passed jaccard test.
#' @export
#'
jaccard_filter <- function(x, out=NULL, jac.cut=0.75) {

  if (is.null(names(x))) {
    stop("Vectors must be named.")
  }
  if (!is.null(out) && length(x) < 1) {
    return(out)
  }
  else if (!is.null(out) && length(x) == 1) {
    out <- c(out, x)
    return(out)
  }
  if (is.null(out) && length(x) <= 1) {
    stop("A minimum of 2 vectors is needed for the jaccard test.")
  }

  # keep first vector
  out <- c(out, x[1])
  # pairs first vector with every other in list
  pairs <- expand.grid(x[1], x[2:length(x)])
  # names of all but the first vector in list
  names <- names(pairs[[2]])
  # good vectors
  boolean <- apply(pairs, 1, function(pair) {
    jaccard(pair[[1]], pair[[2]], jac.cut=jac.cut)})
  names <- names[boolean]

  # proceed with good vectors
  x <- x[names(x) %in% names]

  jaccard_filter(x=x, out=out)
}


# ==============================================================
# CLUSTER SIMILARITY CUT
# ==============================================================

#' Cut filters based on jaccard similarity
#'
#' @param obj hcsig object. Please refer to the documentation for \code{statistrics::hcsig}.
#' @param jac.out the output vector. Usually and by default equal to NULL -- the clusters that pass jaccard similarity will be added to an empty vector.
#' @param jac.cut the cutoff for jaccard similarity.
#'
#' @return an hcsig object (the same structure and data types as the input) filtered to include only non-similar clusters according to \code{jac.cut} val.
#' @export
#'
hcsim_cut <- function(obj, jac.out=NULL, jac.cut=0.75) {

  obj <- input_check(obj)
  k.cut <- jaccard_filter(x=obj$k, out=jac.out, jac.cut=jac.cut)
  kept <- obj$k %in% k.cut

  list(k=k.cut,  sig.1=obj$sig.1[kept], sig.2=obj$sig.2[kept], sig.3=obj$sig.3[kept], fc=obj$fc[kept])
}

# ==============================================================
