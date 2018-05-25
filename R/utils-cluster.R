# =======================
# 2. POST-CLUSTERING
# =======================
# -----------------------
# B. CLEANING CUTS
# -----------------------
#' Remove too big too small clusters
#'
#' @param k chr vector of clusters
#' @param n number of members in cluster
#' @param min cluster size. absolute value
#' @param max cluster size. fraction of total n
#'
#' @return remaining clusters. list of chr vectors
#' @export
#'
hclean <- function(k, n, min=5, max=0.5) {

  # remove duplicate clusters
  k <- k[!duplicated(k)]

  # remove clusters smaller than min val
  k <- k[lapply(k, length) >= min]

  # remove clusters bigger than max val
  k <- k[lapply(k, function(i) length(i)/n) <= max]

  return(k)

}

# ------------------------
# A. RETRIEVING CUTS FROM HC
# ------------------------
#' Get clusters from hc object by tree heights
#'
#' @param hc hierarchical clustering obj.
#' @param h \code{cutree} arg. height at which to cut. hc$height=all heights.
#' @param clean logical. if TRUE, call to \code{hclean}.
#' @param min hclean arg. cluster size. absolute value.
#' @param max hclean arg. cluster size. fraction of total n.
#' @param ... other args passed to \code{cutree}.
#'
#' @return clusters. list of chr vectors.
#' @export
#'
hcutree <- function(hc,
                    h=hc$height,
                    clean=TRUE,
                    min=5,
                    max=0.50,
                    ...) {

  # Cut tree structure. if h=hc$height,
  # all possible sub trees will be generated
  cut <- stats::cutree(hc, h=h, ...)

  # Extract cell names for each cluster
  # Generates nested list
  k <- lapply(1:ncol(cut), function(i) {
    tapply(names(cut[,i]), cut[,i], c)
  })
  # Flatten list. Cell names to top level.
  k <- unlist(k, recursive=FALSE)

  if (isTRUE(clean)) {
    k <- hclean(k, n=ncol(cut) + 1, min=min, max=max)
  }

  names(k) <- 1:length(k)

  k
}

# =======================
# 1. CLUSTERING >> get HC
# =======================
#' Hierarchical clustering
#'
#' @param mat matrix of vars. vs. obs.
#' @param method.cor distance metric for calculating correlations.
#' @param method.hc type of linkage.
#' @param dissim.dist maximum dissimilarity prior to rescaling between 0 and 1.
#' @param ... other args passed to \code{cor} or \code{hclust} calls.
#'
#' @return hc object
#' @export
#'
hcluster <- function(mat,
                     method.cor="pearson",
                     method.hc="average",
                     dissim.dist=1,
                     ...) {

  # Correlation matrix
  cr <- stats::cor(x=mat, method=method.cor, ...)

  # Distance matrix
  d <- stats::as.dist(dissim.dist - cr)

  # Hierarchical cluster analysis
  stats::hclust(d=d, method=method.hc, ...)

}

#' Order variables by hierarchical clustering
#'
#' @param mat matrix of vars. vs. obs.
#'
#' @return chr vector of ordered vars.
#' @export
#'
hcorder <- function(mat) {
  hc <- hcluster(mat)
  hc$labels[hc$order]
}

#' @importFrom reshape2 melt
#' @export
reshape2::melt

#' Calculate (ordered) pairwise correlation matrix
#'
#' @param mat matrix of vars. vs. obs. to be correlated.
#' @param method.cor distance metric for calculating correlations.
#' @param order if TRUE, corr object will be returned with vars ordered by hc.
#' @param melt if TRUE, input matrix will be converted to tidy/long format.
#' @param ... other args passed to \code{cor}.
#'
#' @return matrix of vars. vs. vars correlation scores.
#' @export
#'
hcorr <- function(mat,
                  method.cor="pearson",
                  order=TRUE,
                  melt=FALSE,
                  ...) {
  # Correlation matrix
  corr <- stats::cor(x=mat, method=method.cor, ...)

  # Factor cells by hierarchical clustering
  if (isTRUE(order)) {
    ord <- as.factor(hcorder(mat))
    corr <- corr[ord, ord]
  }

  if (isTRUE(melt)) {
    corr <- reshape2::melt(corr)
  }

  corr
}

#' Get clusters from hc object by predefined cluster number
#'
#' @param hc hierarchical clustering object.
#' @param n number of clusters to return.
#'
#' @return clusters. list of chr vectors.
#' @export
#'
ncutree <- function(hc, n) {
  k <- stats::cutree(hc, k=n)
  k <- as.matrix(k)
  df <- data.frame(k=k[,1], programs=rownames(k))
  s <- split(df, df$k)
  k <- lapply(s, rownames)
  k
}
