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
#' @param greater boolean. If TRUE, return FC values equal to or greater than cutoff, else smaller than or equal to.
#'
#' @return numeric vector containing FC values
#' @export
#'
fold_change <- function(a, b, log=TRUE, log.base=2, fc.value=3, greater=T, return.log = FALSE) {

  a <- as.matrix(a)
  b <- as.matrix(b)

  if (length(dim(a)) == 2) {
    a <- gene_means(a)
    b <- gene_means(b)}

  if (log) {
    logFC <- a - b
    FC <- log.base^logFC
  }

  else FC <- a/b

  if (!is.null(fc.value) & !isTRUE(greater)) {
	return(FC[FC <= fc.value])
  }

  else if (!is.null(fc.value)) {
	return(FC[FC >= fc.value])
  }

  if (isTRUE(return.log)) {
    return(log(FC, log.base))
  }

  FC

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
p_val <- function(a, b, p.value=NULL, adjust.method=NULL){
# p_val <- function(a, b, p.value=NULL, adjust.method="bonferroni"){

  if (nrow(a) < 1) {
	stop("At least one gene/obesrvation is needed for the t.test.")
  }

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
#' @param fc.sort if TRUE, significantly differentially expressed genes are sorted by fold change (highest first). Default is TRUE.
#' @param pval.sort if TRUE, significantly differentially expressed genes are sorted by p.value (highest first). pval.sort=TRUE overrides fc.sort=TRUE. Default is FALSE.
#' @param adjust.method NULL or character string. If NULL, do not adjust p-values. If string, adjust p-values using the method specified.
#'
#' @return one of: i) all p-values for all fold changes (if fc.value == NULL & p.value == NULL) ii/iii) all/signifcant p-values for significant/all fold changes or iv) signifcant p-values for significant fold changes (if fc.value != NULL & p.value != NULL).
#' @export
#'
sig <- function(a, b, p.value, fc.value=3, fc.sort=F, pval.sort=F, adjust.method=NULL, ...){

  fc <- fold_change(a, b, fc.value=fc.value)

  if (length(fc) == 0) {return()}

  else if (length(fc) > 0){
    if (isTRUE(fc.sort)) fc <- sort(fc, decreasing=T)
    a <- a[names(fc), , drop=F]
    b <- b[names(fc), , drop=F]

    pval <- p_val(a, b, p.value=p.value, adjust.method=adjust.method, ...)
    fc <- fc[names(pval)]

    if (isTRUE(pval.sort)) {
      pval.ord <- order(pval, decreasing=F)
      pval <- pval[pval.ord]
      fc <- fc[pval.ord]
    }

    list(fc = fc, pval = pval)
  }

}


#' Differentially Expressed Genes Between Two Groups
#'
#' @param k a list of character vectors; sets of cell members belonging to each cluster.
#' @param mat matrix of vars. vs. obs. matrix is split into two: i) vars that are in k; ii) vars not in k.
#' @param fc.value fold change value below which differential gene expression is deemed insignificant.
#' @param p.value p-value above which differential gene expression is deemed insignificant.
#' @param ... other args to be passed to \code{p_val} or \code{fold_change} through call to \code{sig}.
#' @param fc.sort if TRUE, significantly differentially expressed genes are sorted by fold change (highest first). Default is TRUE.
#' @param pval.sort if TRUE, significantly differentially expressed genes are sorted by p.value (highest first). pval.sort=TRUE overrides fc.sort=TRUE. Default is FALSE.
#' @param adjust.method NULL or character string. If NULL, do not adjust p-values. If string, adjust p-values using the method specified.
#' @param returning return one of p-values, fold changes, or both from the output of \code{sig()} call.
#'
#' @return numeric vector of p-values named with gene names.
#' @export
#'
DEgenes <- function(k, mat, fc.value=3, p.value=10^(-4), fc.sort=T, pval.sort=F, adjust.method=NULL, returning = 'pval', ...) {
  a <- mat[, k, drop=F]
  b <- mat[, !colnames(mat) %in% k, drop=F]
  out <- sig(a, b, fc.value=fc.value, p.value=p.value, fc.sort=fc.sort, pval.sort=pval.sort, adjust.method=adjust.method, ...)

  if (returning %in% c('pval', 'p', 'pvalue')) return(out$pval)
  else if (returning %in% c('fc', 'FC', 'foldchange', 'fold.change')) return(out$fc)
  else if (returning %in% c('both', 'all')) return(out)
  out
}

higher_cutoff <- function(List, cutoff = 10^(-5)) {
  sapply(List, function(L) L[L <= cutoff], USE.NAMES=T, simplify=F)
}

# ==============================================================
# Cluster significance helpers
# ==============================================================

#' Single Most Significant Occurrences of Genes
#'
#' Most significant occurrence of each gene.
#' The data generated allows the parameter n.sig.2 to be generated,
#' for the number of genes in a cluster with the most significant occurrences.
#' @param List a list of named numeric vectors. Names: genes; values: their p-values; vectors: clusters derived from statistrics::hcluster.
#' @param by 'small' or 'large' depending on whether the most significant values are the smallest (eg. p-value) or largest (eg. fold change).
#'
#' @return the inputted List after filtering such that each gene only appears once in the cluster where it is most significant.
#' @export
#'
most_significant <- function(List, by='small') {
  df <- list_to_df(List)

  if (by == 'small') {
    # sort by id and abs(value)
    ordered.df <- df[order(df$name, abs(df$val)), ]
  }
  else if (by == 'large') {
    # sort by id and reverse of abs(value)
    ordered.df <- df[order(df$name, - abs(df$val)), ]
  }

  # take the first row within each id
  return.df <- ordered.df[!duplicated(ordered.df$name), ]
  return.df.2 <- return.df[order(return.df$id), ]

  # add back vectors that were in List but are not in List.out since they had zero genes.
  List.out <- df_to_list(return.df.2, col=return.df.2$id)
  List.zeros <- List[!names(List) %in% names(List.out)]
  List.out <- c(List.out, List.zeros)[names(List)]

  List.out
}


# ==============================================================
# CLUSTER SIGNIFICANCE
# ==============================================================

#' hcsig: hcluster Significance
#'
#' hcluster Significance. For clusters derived from hierarchical clustering (with \code{hclust}), data is retrieved.
#'
#' @param k a list of character vectors; sets of cell members belonging to each cluster.
#' @param mat a matrix of gene expression data (cells by genes)
#' @param fc.value fold change value below which differential gene expression is deemed insignificant.
#' @param p.value p-value above which differential gene expression is deemed insignificant.
#' @param p.value.2 p-value above which differential gene expression is deemed insignificant with higher cutoff (sig.2)
#' @param pval.adjust NULL or character string. If NULL, do not adjust p-values. If string, adjust p-values using the method specified.
#' @param reorder if TRUE, the list of clusters is reordered by most to least significant.
#' @param fc.sort if TRUE, significantly differentially expressed genes are sorted by fold change (highest first). Default is TRUE.
#' @param pval.sort if TRUE, significantly differentially expressed genes are sorted by p.value (highest first). \code{pval.sort=TRUE} overrides \code{fc.sort=TRUE}. Default is FALSE.
#' @param returning return one of p-values, fold changes, or both from call of \code{DEgenes()} to \code{sig()}.
#'
#' @return list of length 4. Each object in the list is also a list. Each list has the same length, which is the length of k arg (the number of clusters). The lists are list$k, same as input; list$sig.1, the significant genes' p-values for each cluster; list$sig.2, list$sig.1 filtered such that only genes with p value higher than p.value.2 are included. list$sig.3, list$sig.1 filtered such that each gene only appears once across the clusters, wherever it had the highest p-value.
#' @export
#'
hcsig <- function(k,
                  mat,
                  fc.value=3,
                  p.value=10^(-4),
                  p.value.2=10^(-5),
                  pval.adjust=NULL,
                  reorder=TRUE,
                  fc.sort=T,
                  pval.sort=F,
                  returning='all') {

  sig.out <- sapply(k, function(kk) DEgenes(k=kk,
                                          mat=mat,
                                          fc.value=fc.value,
                                          p.value=p.value,
                                          adjust.method=pval.adjust,
                                          fc.sort=fc.sort,
                                          pval.sort=pval.sort,
                                          returning=returning),
                  simplify=FALSE,
                  USE.NAMES=TRUE)

  fc <- sapply(sig.out, `[`, 'fc', simplify = T, USE.NAMES = F)
  sig.1 <- sapply(sig.out, `[`, 'pval', simplify = T, USE.NAMES = F)
  sig.2 <- higher_cutoff(sig.1, cutoff = p.value.2)
  sig.3 <- most_significant(sig.1)

    if (isTRUE(reorder)) {
      ord <- cluster_reorder(sig.1, sig.2, sig.3)
      k <- k[ord]
      sig.1 <- sig.1[ord]
      sig.2 <- sig.2[ord]
      sig.3 <- sig.3[ord]
      fc <- fc[ord]
    }

  List <- list(k=k, sig.1=sig.1, sig.2=sig.2, sig.3=sig.3, fc=fc)

  names(List) <- c("k", "sig.1", "sig.2", "sig.3", "fc")

  List

}


# ==============================================================
# CLUSTER SIGNIFICANCE CUT
# ==============================================================

#' Cluster Significance Cut
#'
#' @param obj hcsig object. Please refer to \code{statistrics::hcsig}.
#' @param n.sig.1 significance cutoff for the number of significantly differentially expressed genes per cluster. Defaults to 50. Any clusters that do not pass this cutoff OR/AND that of n.sig.2 are filtered out.
#' @param n.sig.2 significance cutoff for the number of more significantly differentially expressed genes per cluster. Defaults to 10. Any clusters that do not pass this cutoff OR/AND that of n.sig.1 are filtered out.
#' @param by cut based on scores for: 'sig.1', 'sig.2', 'either' or 'both'.
#'
#' @return an hcsig object (the same structure and data types as the input) filtered to include only significant clusters.
#' @export
#'
hcsig_cut <- function(obj, n.sig.1=50, n.sig.2=10, by='both'){

  obj <- input_check(obj)

  sig.1.bool <- lengths(obj$sig.1) >= n.sig.1
  sig.2.bool <- lengths(obj$sig.2) >= n.sig.2

  if (by == 'sig.1') {
    cutter <- sig.1.bool
  }

  else if (by == 'sig.2') {
    cutter <- sig.2.bool
  }

  else if (by=='either') {
    cutter <- sig.1.bool | sig.2.bool
  }

  else if (by=='both') {
    cutter <- sig.1.bool & sig.2.bool
  }

  list(k=obj$k[cutter], sig.1=obj$sig.1[cutter], sig.2=obj$sig.2[cutter], sig.3=obj$sig.3[cutter], fc=obj$fc[cutter])
}
